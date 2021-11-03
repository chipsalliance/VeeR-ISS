// Copyright 2020 Western Digital Corporation or its affiliates.
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cfenv>
#include <array>
#include "CsRegs.hpp"
#include "FpRegs.hpp"
#include "VecRegs.hpp"

using namespace WdRiscv;


template <typename URV>
CsRegs<URV>::CsRegs()
  : regs_(size_t(CsrNumber::MAX_CSR_) + 1)
{
  // Define CSR entries.
  defineMachineRegs();
  defineSupervisorRegs();
  defineUserRegs();
  defineDebugRegs();
  defineVectorRegs();
  defineNonStandardRegs();
}


template <typename URV>
CsRegs<URV>::~CsRegs()
{
  regs_.clear();
  nameToNumber_.clear();
}


template <typename URV>
Csr<URV>*
CsRegs<URV>::defineCsr(const std::string& name, CsrNumber csrn, bool mandatory,
		       bool implemented, URV resetValue, URV writeMask,
		       URV pokeMask, bool isDebug, bool quiet)
{
  size_t ix = size_t(csrn);

  if (ix >= regs_.size())
    return nullptr;

  if (nameToNumber_.count(name))
    {
      if (not quiet)
	std::cerr << "Error: CSR " << name << " already defined\n";
      return nullptr;
    }

  auto& csr = regs_.at(ix);
  if (csr.isDefined())
    {
      if (not quiet)
	std::cerr << "Error: CSR 0x" << std::hex << size_t(csrn) << std::dec
		  << " is already defined as " << csr.getName() << '\n';
      return nullptr;
    }

  PrivilegeMode priv = PrivilegeMode((ix & 0x300) >> 8);
  csr.definePrivilegeMode(priv);

  csr.setDefined(true);

  csr.config(name, csrn, mandatory, implemented, resetValue, writeMask,
	     pokeMask, isDebug);

  nameToNumber_[name] = csrn;
  return &csr;
}


template <typename URV>
Csr<URV>*
CsRegs<URV>::findCsr(const std::string& name)
{
  const auto iter = nameToNumber_.find(name);
  if (iter == nameToNumber_.end())
    return nullptr;

  size_t num = size_t(iter->second);
  if (num >= regs_.size())
    return nullptr;

  return &regs_.at(num);
}


template <typename URV>
Csr<URV>*
CsRegs<URV>::findCsr(CsrNumber number)
{
  size_t ix = size_t(number);
  if (ix >= regs_.size())
    return nullptr;
  return &regs_.at(ix);
}


template <typename URV>
bool
CsRegs<URV>::read(CsrNumber number, PrivilegeMode mode, URV& value) const
{
  auto csr = getImplementedCsr(number);
  if (not csr)
    return false;

  if (mode < csr->privilegeMode())
    return false;

  if (csr->isDebug())
    return false; // Debug-mode register is not accessible by a CSR instruction.

  if (number >= CsrNumber::TDATA1 and number <= CsrNumber::TDATA3)
    return readTdata(number, mode, value);

  if (number == CsrNumber::FFLAGS or number == CsrNumber::FRM)
    {
      auto fcsr = getImplementedCsr(CsrNumber::FCSR);
      if (not fcsr)
        return false;
      value = fcsr->read();
      if (number == CsrNumber::FFLAGS)
        value = value & URV(FpFlags::FcsrMask);
      else
        value = (value & URV(RoundingMode::FcsrMask)) >> URV(RoundingMode::FcsrShift);
      return true;
    }

  // Value of SIP/SIE is that of MIP/MIE modified by SIP/SIE mask
  // and delgation register.
  if (number == CsrNumber::SIP or number == CsrNumber::SIE)
    {
      // Get MIP/MIE
      auto mcsr = getImplementedCsr(CsrNumber(unsigned(number) + 0x200));
      auto deleg = getImplementedCsr(CsrNumber::MIDELEG);
      if (mcsr and deleg)
        value = mcsr->read() & (csr->getReadMask() | deleg->read());
      else
        value = csr->read();
      return true;
    }
          
  value = csr->read();

  if (number >= CsrNumber::PMPADDR0 and number <= CsrNumber::PMPADDR15)
    value = adjustPmpValue(number, value);

  return true;
}


template <typename URV>
void
CsRegs<URV>::enableSupervisorMode(bool flag)
{
  supervisorModeEnabled_ = flag;

  for (auto csrn : { CsrNumber::SSTATUS, CsrNumber::SEDELEG, CsrNumber::SIDELEG,
                      CsrNumber::STVEC, CsrNumber::SIE, CsrNumber::STVEC,
                      CsrNumber::SCOUNTEREN, CsrNumber::SSCRATCH, CsrNumber::SEPC,
                      CsrNumber::SCAUSE, CsrNumber::STVAL, CsrNumber::SIP,
                      CsrNumber::SATP, CsrNumber::MEDELEG, CsrNumber::MIDELEG } )
    {
      auto csr = findCsr(csrn);
      if (not csr)
        {
          std::cerr << "Error: enableSupervisorMode: CSR number 0x"
                    << std::hex << URV(csrn) << " undefined\n";
          assert(0);
        }
      else
        csr->setImplemented(flag);
    }

  if (not flag)
    return;

  typedef InterruptCause IC;

  // In MIP, make writable/pokable bits corresponding to SEIP/STIP/SSIP
  // (supervisor external/timer/software interrupt pending).
  URV extra = URV(1) << unsigned(IC::S_EXTERNAL);
  extra |= URV(1) << unsigned(IC::S_TIMER);
  extra |= URV(1) << unsigned(IC::S_SOFTWARE);

  auto csr = findCsr(CsrNumber::MIP);
  if (csr)
    {
      URV mask = csr->getWriteMask();
      csr->setWriteMask(mask | extra);

      mask = csr->getPokeMask();
      csr->setPokeMask(mask | extra);
    }

  // Same for MIE.
  csr = findCsr(CsrNumber::MIE);
  if (csr)
    {
      URV mask = csr->getWriteMask();
      csr->setWriteMask(mask | extra);

      mask = csr->getPokeMask();
      csr->setPokeMask(mask | extra);
    }
}


template <typename URV>
void
CsRegs<URV>::enableRvf(bool flag)
{
  for (auto csrn : { CsrNumber::FCSR, CsrNumber::FFLAGS, CsrNumber::FRM } )
    {
      auto csr = findCsr(csrn);
      if (not csr)
        {
          std::cerr << "Error: enableRvf: CSR number 0x"
                    << std::hex << URV(csrn) << " undefined\n";
          assert(0);
        }
      else if (not csr->isImplemented())
        csr->setImplemented(flag);
    }
}


template <typename URV>
void
CsRegs<URV>::enableVectorMode(bool flag)
{
  for (auto csrn : { CsrNumber::VSTART, CsrNumber::VXSAT, CsrNumber::VXRM,
		     CsrNumber::VCSR, CsrNumber::VL, CsrNumber::VTYPE,
		     CsrNumber::VLENB } )
    {
      auto csr = findCsr(csrn);
      if (not csr)
        {
          std::cerr << "Error: enableVectorMode: CSR number 0x"
                    << std::hex << URV(csrn) << " undefined\n";
          assert(0);
        }
      else
        csr->setImplemented(flag);
    }
}


template <typename URV>
URV
CsRegs<URV>::legalizeMstatusValue(URV value) const
{
  MstatusFields<URV> fields(value);
  PrivilegeMode mode = PrivilegeMode(fields.bits_.MPP);

  if (fields.bits_.FS == unsigned(FpFs::Dirty) or fields.bits_.XS == unsigned(FpFs::Dirty))
    fields.bits_.SD = 1;
  else
    fields.bits_.SD = 0;

  if (mode == PrivilegeMode::Machine)
    return fields.value_;

  if (mode == PrivilegeMode::Supervisor and not supervisorModeEnabled_)
    mode = PrivilegeMode::User;

  if (mode == PrivilegeMode::Reserved)
    mode = PrivilegeMode::User;

  if (mode == PrivilegeMode::User and not userModeEnabled_)
    mode = PrivilegeMode::Machine;

  fields.bits_.MPP = unsigned(mode);

  return fields.value_;
}


template <typename URV>
bool
CsRegs<URV>::write(CsrNumber number, PrivilegeMode mode, URV value)
{
  Csr<URV>* csr = getImplementedCsr(number);
  if (not csr)
    return false;

  if (mode < csr->privilegeMode() or csr->isReadOnly())
    return false;

  if (csr->isDebug())
    return false; // Debug-mode register is not accessible by a CSR instruction.

  if (isPmpaddrLocked(number))
    {
      recordWrite(number);
      return true;  // Writing a locked PMPADDR register has no effect.
    }

  // fflags and frm are part of fcsr
  if (number == CsrNumber::FFLAGS or number == CsrNumber::FRM or
      number == CsrNumber::FCSR)
    {
      csr->write(value);
      recordWrite(number);
      updateFcsrGroupForWrite(number, value);
      return true;
    }

  // vxsat and vrm are part of vcsr
  if (number == CsrNumber::VXSAT or number == CsrNumber::VXRM or
      number == CsrNumber::VCSR)
    {
      csr->write(value);
      recordWrite(number);
      updateVcsrGroupForWrite(number, value);
      return true;
    }

  if (number >= CsrNumber::TDATA1 and number <= CsrNumber::TDATA3)
    {
      if (not writeTdata(number, mode, value))
	return false;
      recordWrite(number);
      return true;
    }

  // Value of SIP/SIE is that of MIP/MIE modified by SIP/SIE mask
  // and delgation register.
  if (number == CsrNumber::SIP or number == CsrNumber::SIE)
    {
      // Get MIP/MIE
      auto mcsr = getImplementedCsr(CsrNumber(unsigned(number) + 0x200));
      auto deleg = getImplementedCsr(CsrNumber::MIDELEG);
      if (mcsr and deleg)
        {
          URV prevMask = csr->getWriteMask();
          csr->setWriteMask((prevMask | deleg->read()) & mcsr->getWriteMask());
          csr->write(value);
          csr->setWriteMask(prevMask);
        }
      else
        csr->write(value);
      recordWrite(number);
      return true;
    }

  if (number >= CsrNumber::MHPMEVENT3 and number <= CsrNumber::MHPMEVENT31)
    value = legalizeMhpmevent(number, value);
  else if (number >= CsrNumber::PMPCFG0 and number <= CsrNumber::PMPCFG3)
    {
      URV prev = 0;
      peek(number, prev);
      value = legalizePmpcfgValue(prev, value);
    }
  else if (number == CsrNumber::MSTATUS or number == CsrNumber::SSTATUS)
    value = legalizeMstatusValue(value);

  csr->write(value);
  recordWrite(number);

  if (number == CsrNumber::MSTATUS or number == CsrNumber::SSTATUS)
    {
      // Write cannot change SD. Update it with a poke.
      MstatusFields<URV> msf(peekMstatus());
      if (msf.bits_.FS == unsigned(FpFs::Dirty) or msf.bits_.XS == unsigned(FpFs::Dirty))
        msf.bits_.SD = 1;
      else
        msf.bits_.SD = 0;
      csr->poke(msf.value_);
    }

  // Cache interrupt enable.
  if (number == CsrNumber::MSTATUS)
    {
      MstatusFields<URV> fields(csr->read());
      interruptEnable_ = fields.bits_.MIE;
    }

  // Writing MDEAU unlocks mdseac.
  if (number == CsrNumber::MDEAU)
    lockMdseac(false);

  // Writing MEIVT changes the base address in MEIHAP.
  if (number == CsrNumber::MEIVT)
    {
      value = (value >> 10) << 10;  // Clear least sig 10 bits keeping base.
      size_t meihapIx = size_t(CsrNumber::MEIHAP);
      URV meihap = regs_.at(meihapIx).read();
      meihap &= 0x3ff;  // Clear base address bits.
      meihap |= value;  // Copy base address bits from MEIVT.
      regs_.at(meihapIx).poke(meihap);
      recordWrite(CsrNumber::MEIHAP);
    }

  // Writing mcounteren/scounteren changes accessibility of the
  // counters in user/supervisor modes.
  if (number == CsrNumber::MCOUNTEREN or number == CsrNumber::SCOUNTEREN)
    updateCounterPrivilege();

  return true;
}


template <typename URV>
bool
CsRegs<URV>::isWriteable(CsrNumber number, PrivilegeMode mode ) const
{
  const Csr<URV>* csr = getImplementedCsr(number);
  if (not csr)
    return false;

  if (mode < csr->privilegeMode())
    return false;

  if (csr->isReadOnly())
    return false;

  if (csr->isDebug())
    return false;  // Debug-mode register is not accessible by a CSR instruction.

  return true;
}


template <typename URV>
void
CsRegs<URV>::reset()
{
  for (auto& csr : regs_)
    if (csr.isImplemented())
      csr.reset();

  triggers_.reset();
  mPerfRegs_.reset();

  // Cache interrupt enable.
  Csr<URV>* mstatus = getImplementedCsr(CsrNumber::MSTATUS);
  if (mstatus)
    {
      MstatusFields<URV> fields(mstatus->read());
      interruptEnable_ = fields.bits_.MIE;
    }

  mdseacLocked_ = false;
}


template <typename URV>
bool
CsRegs<URV>::configCsr(const std::string& name, bool implemented, URV resetValue,
                       URV mask, URV pokeMask, bool isDebug, bool shared)
{
  auto iter = nameToNumber_.find(name);
  if (iter == nameToNumber_.end())
    return false;

  size_t num = size_t(iter->second);
  if (num >= regs_.size())
    return false;

  return configCsr(CsrNumber(num), implemented, resetValue, mask, pokeMask,
		   isDebug, shared);
}


template <typename URV>
bool
CsRegs<URV>::configCsr(CsrNumber csrNum, bool implemented, URV resetValue,
                       URV mask, URV pokeMask, bool isDebug, bool shared)
{
  if (size_t(csrNum) >= regs_.size())
    {
      std::cerr << "ConfigCsr: CSR number " << size_t(csrNum)
		<< " out of bound\n";
      return false;
    }

  auto& csr = regs_.at(size_t(csrNum));
  if (csr.isMandatory() and not implemented)
    {
      std::cerr << "CSR " << csr.getName() << " is mandatory and is being "
		<< "configured as not-implemented -- configuration ignored.\n";
      return false;
    }

  csr.setImplemented(implemented);
  csr.setInitialValue(resetValue);
  csr.setWriteMask(mask);
  csr.setPokeMask(pokeMask);
  csr.pokeNoMask(resetValue);
  csr.setIsDebug(isDebug);
  csr.setIsShared(shared);

  // Cahche interrupt enable.
  if (csrNum == CsrNumber::MSTATUS)
    {
      MstatusFields<URV> fields(csr.read());
      interruptEnable_ = fields.bits_.MIE;
    }

  return true;
}


template <typename URV>
bool
CsRegs<URV>::configMachineModePerfCounters(unsigned numCounters)
{
  if (numCounters > 29)
    {
      std::cerr << "No more than 29 machine mode performance counters "
		<< "can be defined\n";
      return false;
    }

  unsigned errors = 0;
  bool shared = false;

  for (unsigned i = 0; i < 29; ++i)
    {
      URV resetValue = 0, mask = ~URV(0), pokeMask = ~URV(0);
      if (i >= numCounters)
	mask = pokeMask = 0;

      CsrNumber csrNum = CsrNumber(i + unsigned(CsrNumber::MHPMCOUNTER3));
      bool isDebug = false;
      if (not configCsr(csrNum, true, resetValue, mask, pokeMask, isDebug,
                        shared))
	errors++;

      if (rv32_)
         {
	   csrNum = CsrNumber(i + unsigned(CsrNumber::MHPMCOUNTER3H));
	   if (not configCsr(csrNum, true, resetValue, mask, pokeMask,
                             isDebug, shared))
	     errors++;
	 }

      csrNum = CsrNumber(i + unsigned(CsrNumber::MHPMEVENT3));
      if (not configCsr(csrNum, true, resetValue, mask, pokeMask, isDebug,
                        shared))
	errors++;
    }

  if (errors == 0)
    {
      mPerfRegs_.config(numCounters);
      tiePerfCounters(mPerfRegs_.counters_);
    }

  return errors == 0;
}


template <typename URV>
bool
CsRegs<URV>::configUserModePerfCounters(unsigned numCounters)
{
  if (numCounters > mPerfRegs_.size())
    {
      std::cerr << "User mode number of performance counters (" << numCounters
                << ") cannot exceed that of machine mode ("
                << mPerfRegs_.size() << '\n';
      return false;
    }

  unsigned errors = 0;
  bool shared = false;

  for (unsigned i = 0; i < 29; ++i)
    {
      URV resetValue = 0, mask = ~URV(0), pokeMask = ~URV(0);
      if (i >= numCounters)
	mask = pokeMask = 0;

      CsrNumber csrNum = CsrNumber(i + unsigned(CsrNumber::HPMCOUNTER3));
      bool isDebug = false;
      if (not configCsr(csrNum, true, resetValue, mask, pokeMask, isDebug,
                        shared))
	errors++;

      if (rv32_)
         {
	   csrNum = CsrNumber(i + unsigned(CsrNumber::HPMCOUNTER3H));
	   if (not configCsr(csrNum, true, resetValue, mask, pokeMask,
                             isDebug, shared))
	     errors++;
	 }
    }

  return errors == 0;
}


/// Map a RISCV rounding mode to an fetsetround constant.
static std::array<int, 5> riscvRoundingModeToFe =
  {
   FE_TONEAREST,  // NearsetEven
   FE_TOWARDZERO, // Zero
   FE_DOWNWARD,   // Down
   FE_UPWARD,     // Up
   FE_TONEAREST   // NearestMax
  };


static
inline
int
mapRiscvRoundingModeToFe(RoundingMode mode)
{
  uint32_t ix = uint32_t(mode);
  if (ix < riscvRoundingModeToFe.size())
    return riscvRoundingModeToFe.at(ix);

  // For dynamic mode, it does not matter to what we set the host machine
  // fp mode since it will be changed by the floating point instructions.
  return FE_TONEAREST;
}
  

static
inline
int
setSimulatorRoundingMode(RoundingMode mode)
{
  int previous = std::fegetround();
  int next = mapRiscvRoundingModeToFe(mode);

  if (next != previous)
    std::fesetround(next);

  return previous;
}


template <typename URV>
void
CsRegs<URV>::updateFcsrGroupForWrite(CsrNumber number, URV value)
{
  if (number == CsrNumber::FFLAGS)
    {
      auto fcsr = getImplementedCsr(CsrNumber::FCSR);
      if (fcsr)
	{
          URV mask = URV(FpFlags::FcsrMask);
	  URV fcsrVal = fcsr->read();
          fcsrVal = (fcsrVal & ~mask) | (value & mask);
	  fcsr->write(fcsrVal);
	  // recordWrite(CsrNumber::FCSR);
	}
      return;
    }

  if (number == CsrNumber::FRM)
    {
      auto fcsr = getImplementedCsr(CsrNumber::FCSR);
      if (fcsr)
	{
	  URV fcsrVal = fcsr->read();
          URV mask = URV(RoundingMode::FcsrMask);
          URV shift = URV(RoundingMode::FcsrShift);
          fcsrVal = (fcsrVal & ~mask) | ((value << shift) & mask);
	  fcsr->write(fcsrVal);
	  // recordWrite(CsrNumber::FCSR);
          setSimulatorRoundingMode(RoundingMode((fcsrVal & mask) >> shift));
	}
      return;
    }

  if (number == CsrNumber::FCSR)
    {
      URV newVal = value & URV(FpFlags::FcsrMask);
      auto fflags = getImplementedCsr(CsrNumber::FFLAGS);
      if (fflags and fflags->read() != newVal)
	{
	  fflags->write(newVal);
	  // recordWrite(CsrNumber::FFLAGS);
	}

      newVal = (value & URV(RoundingMode::FcsrMask)) >> URV(RoundingMode::FcsrShift);
      auto frm = getImplementedCsr(CsrNumber::FRM);
      if (frm and frm->read() != newVal)
	{
	  frm->write(newVal);
	  // recordWrite(CsrNumber::FRM);
	}
      setSimulatorRoundingMode(RoundingMode(newVal));
    }
}


template <typename URV>
void
CsRegs<URV>::updateFcsrGroupForPoke(CsrNumber number, URV value)
{
  if (number == CsrNumber::FFLAGS)
    {
      auto fcsr = getImplementedCsr(CsrNumber::FCSR);
      if (fcsr)
	{
          URV mask = URV(FpFlags::FcsrMask);
	  URV fcsrVal = fcsr->read();
          fcsrVal = (fcsrVal & ~mask) | (value & mask);
	  fcsr->poke(fcsrVal);
	}
      return;
    }

  if (number == CsrNumber::FRM)
    {
      auto fcsr = getImplementedCsr(CsrNumber::FCSR);
      if (fcsr)
	{
	  URV fcsrVal = fcsr->read();
          URV mask = URV(RoundingMode::FcsrMask);
          URV shift = URV(RoundingMode::FcsrShift);
          fcsrVal = (fcsrVal & ~mask) | ((value << shift) & mask);
	  fcsr->poke(fcsrVal);
          setSimulatorRoundingMode(RoundingMode((fcsrVal & mask) >> shift));
	}
      return;
    }

  if (number == CsrNumber::FCSR)
    {
      URV newVal = value & URV(FpFlags::FcsrMask);
      auto fflags = getImplementedCsr(CsrNumber::FFLAGS);
      if (fflags and fflags->read() != newVal)
        fflags->poke(newVal);

      newVal = (value & URV(RoundingMode::FcsrMask)) >> URV(RoundingMode::FcsrShift);
      auto frm = getImplementedCsr(CsrNumber::FRM);
      if (frm and frm->read() != newVal)
        frm->poke(newVal);
      setSimulatorRoundingMode(RoundingMode(newVal));
    }
}


template <typename URV>
void
CsRegs<URV>::updateVcsrGroupForWrite(CsrNumber number, URV value)
{
  if (number == CsrNumber::VXSAT)
    {
      auto vcsr = getImplementedCsr(CsrNumber::VCSR);
      if (vcsr)
	{
          URV mask = 1;
	  URV vcsrVal = vcsr->read();
          vcsrVal = (vcsrVal & ~mask) | (value & mask);
	  vcsr->write(vcsrVal);
	  // recordWrite(CsrNumber::VCSR);
	}
      return;
    }

  if (number == CsrNumber::VXRM)
    {
      auto vcsr = getImplementedCsr(CsrNumber::VCSR);
      if (vcsr)
	{
	  URV vcsrVal = vcsr->read();
          URV mask = URV(VecRoundingMode::VcsrMask);
          URV shift = URV(VecRoundingMode::VcsrShift);
          vcsrVal = (vcsrVal & ~mask) | ((value << shift) & mask);
	  vcsr->write(vcsrVal);
	  // recordWrite(CsrNumber::VCSR);
	}
      return;
    }

  if (number == CsrNumber::VCSR)
    {
      URV newVal = value & 1;
      auto vxsat = getImplementedCsr(CsrNumber::VXSAT);
      if (vxsat and vxsat->read() != newVal)
	{
	  vxsat->write(newVal);
	  // recordWrite(CsrNumber::VXSAT);
	}

      newVal = (value & URV(VecRoundingMode::VcsrMask)) >> URV(VecRoundingMode::VcsrShift);
      auto vxrm = getImplementedCsr(CsrNumber::VXRM);
      if (vxrm and vxrm->read() != newVal)
	{
	  vxrm->write(newVal);
	  // recordWrite(CsrNumber::VXRM);
	}
    }
}


template <typename URV>
void
CsRegs<URV>::updateVcsrGroupForPoke(CsrNumber number, URV value)
{
  if (number == CsrNumber::VXSAT)
    {
      auto vcsr = getImplementedCsr(CsrNumber::VCSR);
      if (vcsr)
	{
          URV mask = 1;
	  URV vcsrVal = vcsr->read();
          vcsrVal = (vcsrVal & ~mask) | (value & mask);
	  vcsr->poke(vcsrVal);
	}
      return;
    }

  if (number == CsrNumber::VXRM)
    {
      auto vcsr = getImplementedCsr(CsrNumber::VCSR);
      if (vcsr)
	{
	  URV vcsrVal = vcsr->read();
          URV mask = URV(VecRoundingMode::VcsrMask);
          URV shift = URV(VecRoundingMode::VcsrShift);
          vcsrVal = (vcsrVal & ~mask) | ((value << shift) & mask);
	  vcsr->poke(vcsrVal);
	}
      return;
    }

  if (number == CsrNumber::VCSR)
    {
      URV newVal = value & 1;
      auto vxsat = getImplementedCsr(CsrNumber::VXSAT);
      if (vxsat and vxsat->read() != newVal)
        vxsat->poke(newVal);

      newVal = (value & URV(VecRoundingMode::VcsrMask)) >> URV(VecRoundingMode::VcsrShift);
      auto vxrm = getImplementedCsr(CsrNumber::VXRM);
      if (vxrm and vxrm->read() != newVal)
        vxrm->poke(newVal);
    }
}


template <typename URV>
void
CsRegs<URV>::recordWrite(CsrNumber num)
{
  auto& lwr = lastWrittenRegs_;
  if (std::find(lwr.begin(), lwr.end(), num) == lwr.end())
    lwr.push_back(num);
}


template <typename URV>
void
CsRegs<URV>::defineMachineRegs()
{
  URV rom = 0;        // Read-only mask: no bit writeable.
  URV wam = ~URV(0);  // Write-all mask: all bits writeable.

  bool mand = true;  // Mandatory.
  bool imp = true;   // Implemented.

  using Csrn = CsrNumber;

  // Machine info.
  defineCsr("mvendorid", Csrn::MVENDORID, mand, imp, 0, rom, rom);
  defineCsr("marchid",   Csrn::MARCHID,   mand, imp, 0, rom, rom);
  defineCsr("mimpid",    Csrn::MIMPID,    mand, imp, 0, rom, rom);
  defineCsr("mhartid",   Csrn::MHARTID,   mand, imp, 0, rom, rom);

  // Machine status setup.

  // mstatus
  //           S R        T T T M S M X  F  M  V  S M R S U M R S U
  //           D E        S W V X U P S  S  P  S  P P E P P I E I I
  //             S        R   M R M R       P     P I S I I E S E E
  //                                V               E   E E
  URV mask = 0b0'00000000'1'1'1'1'1'1'11'11'11'11'1'1'0'1'0'1'0'1'0;
  URV val = 0;
  if (not rv32_)
    {
      mask |= uint64_t(0b0000) << 32;  // Mask for SXL and UXL (currently not writable).
      val |= uint64_t(0b1010) << 32;   // Value of SXL and UXL : sxlen=uxlen=64
    }
  URV pokeMask = mask | (URV(1) << (sizeof(URV)*8 - 1));  // Make SD pokable.

  defineCsr("mstatus", Csrn::MSTATUS, mand, imp, val, mask, pokeMask);
  defineCsr("misa", Csrn::MISA, mand,  imp, 0x40001104, rom, rom);

  // Bits corresponding to user-level interrupts are hardwired to zero
  // in medeleg. If N extension is enabled, we will flip those bits
  // (currently N extension is not supported).
  URV userBits = ( (URV(1) << unsigned(InterruptCause::U_SOFTWARE)) |
                   (URV(1) << unsigned(InterruptCause::U_TIMER)) |
                   (URV(1) << unsigned(InterruptCause::U_EXTERNAL)) );
  mask = wam & ~ userBits;
  defineCsr("medeleg", Csrn::MEDELEG, !mand, !imp, 0, mask, mask);

  defineCsr("mideleg", Csrn::MIDELEG, !mand, !imp, 0, wam, wam);

  // Interrupt enable: Least sig 12 bits corresponding to the 12
  // interrupt causes are writable.
  URV mieMask = 0xfff; 
  defineCsr("mie", Csrn::MIE, mand, imp, 0, mieMask, mieMask);

  // Initial value of 0: vectored interrupt. Mask of ~2 to make bit 1
  // non-writable.
  mask = ~URV(2);
  defineCsr("mtvec", Csrn::MTVEC, mand, imp, 0, mask, mask);

  defineCsr("mcounteren", Csrn::MCOUNTEREN, !mand, imp, 0, wam, wam);
  defineCsr("mcountinhibit", Csrn::MCOUNTINHIBIT, !mand, imp, 0, wam, wam);

  // Machine trap handling: mscratch and mepc.
  defineCsr("mscratch", Csrn::MSCRATCH, mand, imp, 0, wam, wam);
  mask = ~URV(1);  // Bit 0 of MEPC is not writable.
  defineCsr("mepc", Csrn::MEPC, mand, imp, 0, mask, mask);

  // All bits of mcause writeable.
  defineCsr("mcause", Csrn::MCAUSE, mand, imp, 0, wam, wam);
  defineCsr("mtval", Csrn::MTVAL, mand, imp, 0, wam, wam);

  // MIP is read-only for CSR instructions but the bits corresponding
  // to defined interrupts are modifiable.
  defineCsr("mip", CsrNumber::MIP, mand, imp, 0, rom, mieMask);

  // Physical memory protection. PMPCFG1 and PMPCFG3 are present only
  // in 32-bit implementations.
  uint64_t cfgMask = 0x9f9f9f9f;
  if (not rv32_)
    cfgMask = 0x9f9f9f9f9f9f9f9fL;
  defineCsr("pmpcfg0",   Csrn::PMPCFG0,   !mand, imp, 0, cfgMask, cfgMask);
  defineCsr("pmpcfg2",   Csrn::PMPCFG2,   !mand, imp, 0, cfgMask, cfgMask);
  if (rv32_)
    {
      defineCsr("pmpcfg1",   Csrn::PMPCFG1,   !mand, imp, 0, cfgMask, cfgMask);
      defineCsr("pmpcfg3",   Csrn::PMPCFG3,   !mand, imp, 0, cfgMask, cfgMask);
    }
  else
    {
      defineCsr("pmpcfg1",   Csrn::PMPCFG1,   !mand, !imp, 0, cfgMask, cfgMask);
      defineCsr("pmpcfg3",   Csrn::PMPCFG3,   !mand, !imp, 0, cfgMask, cfgMask);
    }

  uint64_t pmpMask = 0xffffffff;
  if (not rv32_)
    pmpMask = 0x003f'ffff'ffff'ffffL; // Top 10 bits are zeros

  defineCsr("pmpaddr0",  Csrn::PMPADDR0,  !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr1",  Csrn::PMPADDR1,  !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr2",  Csrn::PMPADDR2,  !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr3",  Csrn::PMPADDR3,  !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr4",  Csrn::PMPADDR4,  !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr5",  Csrn::PMPADDR5,  !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr6",  Csrn::PMPADDR6,  !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr7",  Csrn::PMPADDR7,  !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr8",  Csrn::PMPADDR8,  !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr9",  Csrn::PMPADDR9,  !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr10", Csrn::PMPADDR10, !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr11", Csrn::PMPADDR11, !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr12", Csrn::PMPADDR12, !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr13", Csrn::PMPADDR13, !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr14", Csrn::PMPADDR14, !mand, imp, 0, pmpMask, pmpMask);
  defineCsr("pmpaddr15", Csrn::PMPADDR15, !mand, imp, 0, pmpMask, pmpMask);

  // Machine Counter/Timers.
  defineCsr("mcycle",    Csrn::MCYCLE,    mand, imp, 0, wam, wam);
  defineCsr("minstret",  Csrn::MINSTRET,  mand, imp, 0, wam, wam);
  if (rv32_)
    {
      defineCsr("mcycleh",   Csrn::MCYCLEH,   mand, imp, 0, wam, wam);
      defineCsr("minstreth", Csrn::MINSTRETH, mand, imp, 0, wam, wam);
    }

  // Define mhpmcounter3/mhpmcounter3h to mhpmcounter31/mhpmcounter31h
  // as write-anything/read-zero (user can change that in the config
  // file by setting the number of writeable counters). Same for
  // mhpmevent3/mhpmevent3h to mhpmevent3h/mhpmevent31h.
  for (unsigned i = 3; i <= 31; ++i)
    {
      CsrNumber csrNum = CsrNumber(unsigned(CsrNumber::MHPMCOUNTER3) + i - 3);
      std::string name = "mhpmcounter" + std::to_string(i);
      defineCsr(name, csrNum, mand, imp, 0, rom, rom);

      if (rv32_)
        {
          // High register counterpart of mhpmcounter.
          name += "h";
          csrNum = CsrNumber(unsigned(CsrNumber::MHPMCOUNTER3H) + i - 3);
          bool hmand = rv32_;  // high counters mandatory only in rv32
          defineCsr(name, csrNum, hmand, imp, 0, rom, rom);
        }

      csrNum = CsrNumber(unsigned(CsrNumber::MHPMEVENT3) + i - 3);
      name = "mhpmevent" + std::to_string(i);
      defineCsr(name, csrNum, mand, imp, 0, rom, rom);
    }
}


template <typename URV>
void
CsRegs<URV>::tieSharedCsrsTo(CsRegs<URV>& target)
{
  if (this == &target)
    return;

  assert(regs_.size() == target.regs_.size());
  for (size_t i = 0; i < regs_.size(); ++i)
    {
      CsrNumber csrn = CsrNumber(i);
      auto csr = getImplementedCsr(csrn);
      auto targetCsr = target.getImplementedCsr(csrn);
       if (csr)
        {
          assert(targetCsr);
          if (csr->isShared())
            {
              assert(targetCsr->isShared());
              csr->tie(targetCsr->valuePtr_);
            }
        }
      else
        assert(not targetCsr);
    }
}


template <typename URV>
void
CsRegs<URV>::tiePerfCounters(std::vector<uint64_t>& counters)
{
  // Since the user-mode counters are a shadow of their machine-mode
  // counterparts, we tie them as well regardless of whether or not
  // they are configured.

  if (rv32_)
    {
      // Tie each mhpmcounter CSR value to the least significant 4
      // bytes of the corresponding counters_ entry. Tie each
      // mhpmcounterh CSR value to the most significan 4 bytes of the
      // corresponding counters_ entry.
      for (unsigned num = 3; num <= 31; ++num)
	{
	  unsigned ix = num - 3;
	  if (ix >= counters.size())
	    break;

	  unsigned highIx = ix +  unsigned(CsrNumber::MHPMCOUNTER3H);
	  Csr<URV>& csrHigh = regs_.at(highIx);
	  URV* low = reinterpret_cast<URV*>(&counters.at(ix));
          URV* high = low + 1;
	  csrHigh.tie(high);

	  unsigned lowIx = ix +  unsigned(CsrNumber::MHPMCOUNTER3);
	  Csr<URV>& csrLow = regs_.at(lowIx);
	  csrLow.tie(low);

          // Tie the user-mode performance counter to their
          // machine-mode counterparts.
          highIx = ix +  unsigned(CsrNumber::HPMCOUNTER3H);
          regs_.at(highIx).tie(high);
          lowIx = ix +  unsigned(CsrNumber::HPMCOUNTER3);
          regs_.at(lowIx).tie(low);
	}
    }
  else
    {
      for (unsigned num = 3; num <= 31; ++num)
	{
	  unsigned ix = num - 3;
	  if (ix >= counters.size())
	    break;
	  unsigned csrIx = ix +  unsigned(CsrNumber::MHPMCOUNTER3);
	  Csr<URV>& csr = regs_.at(csrIx);
	  URV* loc = reinterpret_cast<URV*>(&counters.at(ix));
	  csr.tie(loc);

          // Tie user-mode perf register to corresponding machine mode reg.
          csrIx = ix +  unsigned(CsrNumber::HPMCOUNTER3);
          regs_.at(csrIx).tie(loc);
	}
    }
}


template <typename URV>
void
CsRegs<URV>::defineSupervisorRegs()
{
  bool mand = true;   // Mandatory.
  bool imp = true;    // Implemented.
  URV wam = ~ URV(0); // Write-all mask: all bits writeable.

  // Supervisor trap SETUP_CSR.

  using Csrn = CsrNumber;

  // Only bits sie, spie, upie, ube, spp, fs, xs, sum, mxr, uxl (rv64) and sd of
  // sstatus are writeable.  The non-writeable bits read zero.
  uint64_t mask = 0x800de162;
  if (not rv32_)
    mask = 0x80000003000de162L;
  defineCsr("sstatus",    Csrn::SSTATUS,    !mand, !imp, 0, mask, mask);

  auto sstatus = findCsr(Csrn::SSTATUS);
  if (sstatus)
    sstatus->setReadMask(mask);

  // SSTATUS shadows MSTATUS
  auto mstatus = findCsr(Csrn::MSTATUS);
  if (sstatus and mstatus)
    sstatus->tie(mstatus->valuePtr_);

  defineCsr("sedeleg",    Csrn::SEDELEG,    !mand, !imp, 0, wam, wam);
  defineCsr("sideleg",    Csrn::SIDELEG,    !mand, !imp, 0, wam, wam);

  defineCsr("stvec",      Csrn::STVEC,      !mand, !imp, 0, wam, wam);
  defineCsr("scounteren", Csrn::SCOUNTEREN, !mand, !imp, 0, wam, wam);

  // Supervisor Trap Handling 
  defineCsr("sscratch",   Csrn::SSCRATCH,   !mand, !imp, 0, wam, wam);
  mask = ~URV(1);  // Bit 0 of SEPC is not writable.
  defineCsr("sepc",       Csrn::SEPC,       !mand, !imp, 0, mask, mask);
  defineCsr("scause",     Csrn::SCAUSE,     !mand, !imp, 0, wam, wam);
  defineCsr("stval",      Csrn::STVAL,      !mand, !imp, 0, wam, wam);

  // Bits of SIE appear hardwired to zreo unless delegated.
  defineCsr("sie",        Csrn::SIE,        !mand, !imp, 0, wam, wam);
  auto sie = findCsr(Csrn::SIE);
  auto mie = findCsr(Csrn::MIE);
  if (sie and mie)
    sie->tie(mie->valuePtr_);

  // Bits of SIE appear hardwired to zreo unless delegated.
  mask = 0x2;  // Only ssie bit writable (when delegated)
  defineCsr("sip",        Csrn::SIP,        !mand, !imp, 0, mask, mask);

  auto sip = findCsr(Csrn::SIP);
  auto mip = findCsr(Csrn::MIP);
  if (sip and mip)
    sip->tie(mip->valuePtr_); // Sip is a shadow if mip

  // Supervisor Protection and Translation 
  defineCsr("satp",       Csrn::SATP,       !mand, !imp, 0, wam, wam);
}


template <typename URV>
void
CsRegs<URV>::defineUserRegs()
{
  bool mand = true;    // Mandatory.
  bool imp  = true;    // Implemented.
  URV  wam  = ~URV(0); // Write-all mask: all bits writeable.

  using Csrn = CsrNumber;

  // User trap setup.
  URV mask = 0x11; // Only UPIE and UIE bits are writeable.
  defineCsr("ustatus",  Csrn::USTATUS,  !mand, !imp, 0, mask, mask);
  defineCsr("uie",      Csrn::UIE,      !mand, !imp, 0, wam, wam);
  defineCsr("utvec",    Csrn::UTVEC,    !mand, !imp, 0, wam, wam);

  // User Trap Handling
  defineCsr("uscratch", Csrn::USCRATCH, !mand, !imp, 0, wam, wam);
  defineCsr("uepc",     Csrn::UEPC,     !mand, !imp, 0, wam, wam);
  defineCsr("ucause",   Csrn::UCAUSE,   !mand, !imp, 0, wam, wam);
  defineCsr("utval",    Csrn::UTVAL,    !mand, !imp, 0, wam, wam);
  defineCsr("uip",      Csrn::UIP,      !mand, !imp, 0, wam, wam);

  // User Floating-Point CSRs
  defineCsr("fflags",   Csrn::FFLAGS,   !mand, !imp, 0, wam, wam);
  defineCsr("frm",      Csrn::FRM,      !mand, !imp, 0, wam, wam);
  defineCsr("fcsr",     Csrn::FCSR,     !mand, !imp, 0, 0xff, 0xff);

  // User Counter/Timers
  defineCsr("cycle",    Csrn::CYCLE,    !mand, imp,  0, wam, wam);
  defineCsr("time",     Csrn::TIME,     !mand, imp,  0, wam, wam);
  defineCsr("instret",  Csrn::INSTRET,  !mand, imp,  0, wam, wam);
  defineCsr("cycleh",   Csrn::CYCLEH,   !mand, !imp, 0, wam, wam);
  defineCsr("timeh",    Csrn::TIMEH,    !mand, !imp, 0, wam, wam);
  defineCsr("instreth", Csrn::INSTRETH, !mand, !imp, 0, wam, wam);

  // Define hpmcounter3/hpmcounter3h to hpmcounter31/hpmcounter31h
  // as write-anything/read-zero (user can change that in the config
  // file).  Same for mhpmevent3/mhpmevent3h to mhpmevent3h/mhpmevent31h.
  for (unsigned i = 3; i <= 31; ++i)
    {
      CsrNumber csrNum = CsrNumber(unsigned(CsrNumber::HPMCOUNTER3) + i - 3);
      std::string name = "hpmcounter" + std::to_string(i);
      defineCsr(name, csrNum, !mand, !imp, 0, wam, wam);

      // High register counterpart of mhpmcounter.
      name += "h";
      csrNum = CsrNumber(unsigned(CsrNumber::HPMCOUNTER3H) + i - 3);
      defineCsr(name, csrNum, !mand, !imp, 0, wam, wam);
    }
}


template <typename URV>
void
CsRegs<URV>::defineDebugRegs()
{
  bool mand = true; // Mandatory.
  bool imp = true;  // Implemented.
  URV wam = ~URV(0);  // Write-all mask: all bits writeable.

  using Csrn = CsrNumber;

  // Debug/Trace registers.
  defineCsr("tselect", Csrn::TSELECT, !mand, imp,  0, wam, wam);
  defineCsr("tdata1",  Csrn::TDATA1,  !mand, imp,  0, wam, wam);
  defineCsr("tdata2",  Csrn::TDATA2,  !mand, imp,  0, wam, wam);
  defineCsr("tdata3",  Csrn::TDATA3,  !mand, !imp, 0, wam, wam);

  // Define triggers.
  unsigned triggerCount = 4;
  triggers_ = Triggers<URV>(triggerCount);

  Data1Bits<URV> data1Mask(0), data1Val(0);

  // Set the masks of the read-write fields of data1 to all 1.
  data1Mask.mcontrol_.dmode_   = 1;
  data1Mask.mcontrol_.hit_     = 1;
  data1Mask.mcontrol_.select_  = 1;
  data1Mask.mcontrol_.action_  = 1; // Only least sig bit writeable
  data1Mask.mcontrol_.chain_   = 1;
  data1Mask.mcontrol_.match_   = 1; // Only least sig bit of match is writeable.
  data1Mask.mcontrol_.m_       = 1;
  data1Mask.mcontrol_.execute_ = 1;
  data1Mask.mcontrol_.store_   = 1;
  data1Mask.mcontrol_.load_    = 1;

  // Set intitial values of fields of data1.
  data1Val.mcontrol_.type_ = unsigned(TriggerType::AddrData);
  data1Val.mcontrol_.maskMax_ = rv32_ ? 31 : 63;

  // Values, write-masks, and poke-masks of the three components of
  // the triggres.
  URV val1(data1Val.value_), val2(0), val3(0);
  URV wm1(data1Mask.value_), wm2(~URV(0)), wm3(0);
  URV pm1(wm1), pm2(wm2), pm3(wm3);

  triggers_.config(0, val1, val2, val3, wm1, wm2, wm3, pm1, pm2, pm3);
  triggers_.config(1, val1, val2, val3, wm1, wm2, wm3, pm1, pm2, pm3);
  triggers_.config(2, val1, val2, val3, wm1, wm2, wm3, pm1, pm2, pm3);

  Data1Bits<URV> icountMask(0), icountVal(0);

  icountMask.icount_.dmode_  = 1;
  icountMask.icount_.count_  = (~0) & 0x3fff;
  icountMask.icount_.m_      = 1;
  icountMask.icount_.action_ = 0;
  icountMask.icount_.action_ = (~0) & 0x3f;

  icountVal.icount_.type_ = unsigned(TriggerType::InstCount);
  icountVal.icount_.count_ = 0;

  triggers_.config(3, icountVal.value_, 0, 0, icountMask.value_, 0, 0,
		   icountMask.value_, 0, 0);

  hasActiveTrigger_ = triggers_.hasActiveTrigger();
  hasActiveInstTrigger_ = triggers_.hasActiveInstTrigger();

  // Debug mode registers.
  URV dcsrVal = 0x40000003;
  URV dcsrMask = 0x00008e04;
  URV dcsrPokeMask = dcsrMask | 0x1c8; // Cause field modifiable
  bool isDebug = true;
  defineCsr("dcsr", Csrn::DCSR, !mand, imp, dcsrVal, dcsrMask,
	    dcsrPokeMask, isDebug);

  // Least sig bit of dpc is not writeable.
  URV dpcMask = ~URV(1);
  defineCsr("dpc", CsrNumber::DPC, !mand, imp, 0, dpcMask, dpcMask, isDebug);

  defineCsr("dscratch", CsrNumber::DSCRATCH, !mand, !imp, 0, wam, wam,
	    isDebug);
}


template <typename URV>
void
CsRegs<URV>::defineVectorRegs()
{
  bool mand = true;  // Mndatory
  bool imp = true;   // Implemented

  defineCsr("vstart", CsrNumber::VSTART, !mand, !imp, 0, 0, 0);
  defineCsr("vxsat",  CsrNumber::VXSAT,  !mand, !imp, 0, 1, 1);  // 1 bit
  defineCsr("vxrm",   CsrNumber::VXRM,   !mand, !imp, 0, 3, 3);  // 2 bits
  defineCsr("vcsr",   CsrNumber::VCSR,   !mand, !imp, 0, 7, 7);  // 3 bits
  URV pokeMask = ~URV(0);
  defineCsr("vl",     CsrNumber::VL,     !mand, !imp, 0, 0, pokeMask);

  uint64_t mask = 0x800000ff;
  if (not rv32_)
    mask = 0x80000000000000ffL;
  defineCsr("vtype",  CsrNumber::VTYPE,  !mand, !imp, 0, mask, mask);

  defineCsr("vlenb",  CsrNumber::VLENB,  !mand, !imp, 0, 0, 0);
}


template <typename URV>
void
CsRegs<URV>::defineNonStandardRegs()
{
  URV rom = 0;        // Read-only mask: no bit writeable.
  URV wam = ~URV(0);  // Write-all mask: all bits writeable.

  bool mand = true; // Mandatory.
  bool imp = true;  // Implemented.

  using Csrn = CsrNumber;

  // mdseac is read-only to CSR insts but is modifiable with poke.
  defineCsr("mdseac", Csrn::MDSEAC,   !mand, imp, 0, rom, wam);

  // mdeau is write-only, it unlocks mdseac when written, it always
  // reads zero.
  defineCsr("mdeau",  Csrn::MDEAU,    !mand, imp, 0, rom, rom);

  // Least sig 10 bits of interrupt vector table (meivt) are read only.
  URV mask = (~URV(0)) << 10;
  defineCsr("meivt",  Csrn::MEIVT,    !mand, imp, 0, mask, mask);

  // None of the bits are writeable by CSR instructions. All but least
  // sig 2 bis are modifiable.
  defineCsr("meihap", Csrn::MEIHAP,   !mand, imp, 0, rom, ~URV(3));

  defineCsr("mscause",  Csrn::MSCAUSE, !mand, !imp, 0, wam, wam);
}


template <typename URV>
bool
CsRegs<URV>::peek(CsrNumber number, URV& value) const
{
  auto csr = getImplementedCsr(number);
  if (not csr)
    return false;

  if (number >= CsrNumber::TDATA1 and number <= CsrNumber::TDATA3)
    return readTdata(number, PrivilegeMode::Machine, value);

  if (number == CsrNumber::FFLAGS or number == CsrNumber::FRM)
    {
      auto fcsr = getImplementedCsr(CsrNumber::FCSR);
      if (not fcsr)
        return false;
      value = fcsr->read();
      if (number == CsrNumber::FFLAGS)
        value = value & URV(FpFlags::FcsrMask);
      else
        value = (value & URV(RoundingMode::FcsrMask)) >> URV(RoundingMode::FcsrShift);
      return true;
    }

  // Value of SIP/SIE is that of MIP/MIE modified by SIP/SIE mask
  // and delgation register.
  if (number == CsrNumber::SIP or number == CsrNumber::SIE)
    {
      // Get MIP/MIE
      auto mcsr = getImplementedCsr(CsrNumber(unsigned(number) + 0x200));
      auto deleg = getImplementedCsr(CsrNumber::MIDELEG);
      if (mcsr and deleg)
        value = mcsr->read() & (csr->getReadMask() | deleg->read());
      else
        value = csr->read();
      return true;
    }

  value = csr->read();

  if (number >= CsrNumber::PMPADDR0 and number <= CsrNumber::PMPADDR15)
    value = adjustPmpValue(number, value);

  return true;
}
  

template <typename URV>
bool
CsRegs<URV>::poke(CsrNumber number, URV value)
{
  Csr<URV>* csr = getImplementedCsr(number);
  if (not csr)
    return false;

  if (isPmpaddrLocked(number))
    return true;  // Writing a locked PMPADDR register has no effect.

  // fflags and frm are parts of fcsr
  if (number == CsrNumber::FFLAGS or number == CsrNumber::FRM or
      number == CsrNumber::FCSR)
    {
      csr->poke(value);
      updateFcsrGroupForPoke(number, value);
      return true;
    }

  // fflags and frm are parts of fcsr
  if (number == CsrNumber::VXSAT or number == CsrNumber::VXRM or
      number == CsrNumber::VCSR)
    {
      csr->poke(value);
      updateVcsrGroupForPoke(number, value);
      return true;
    }

  if (number >= CsrNumber::TDATA1 and number <= CsrNumber::TDATA3)
    return pokeTdata(number, value);

  // Value of SIE/SIP is anded with mask of SIE/SIP anded with mask of
  // MIE/MIP anded with delgation register. For a bit to be written in
  // SIE/SIP it has to be writeable there and in MIE/MIP and be
  // delegated.
  if (number == CsrNumber::SIE or number == CsrNumber::SIP)
    {
      // Get MIE
      auto mcsr = getImplementedCsr(CsrNumber(unsigned(number) + 0x200));
      auto deleg = getImplementedCsr(CsrNumber::MIDELEG);
      if (mcsr and deleg)
        {
          URV prevMask = csr->getWriteMask();
          csr->setWriteMask(deleg->read() & mcsr->getWriteMask() & csr->getWriteMask());
          csr->write(value);
          csr->setWriteMask(prevMask);
        }
      else
        csr->write(value);
      return true;
    }

  if (number >= CsrNumber::MHPMEVENT3 and number <= CsrNumber::MHPMEVENT31)
    value = legalizeMhpmevent(number, value);
  else if (number >= CsrNumber::PMPCFG0 and number <= CsrNumber::PMPCFG3)
    {
      URV prev = 0;
      peek(number, prev);
      value = legalizePmpcfgValue(prev, value);
    }
  else if (number == CsrNumber::MSTATUS or number == CsrNumber::SSTATUS)
    value = legalizeMstatusValue(value);

  csr->poke(value);

  // Cache interrupt enable.
  if (number == CsrNumber::MSTATUS)
    {
      MstatusFields<URV> fields(csr->read());
      interruptEnable_ = fields.bits_.MIE;
    }

  // Poking MDEAU unlocks mdseac.
  if (number == CsrNumber::MDEAU)
    lockMdseac(false);

  // Poking MEIVT changes the base address in MEIHAP.
  if (number == CsrNumber::MEIVT)
    {
      value = (value >> 10) << 10;  // Clear least sig 10 bits keeping base.
      size_t meihapIx = size_t(CsrNumber::MEIHAP);
      URV meihap = regs_.at(meihapIx).read();
      meihap &= 0x3ff;  // Clear base address bits.
      meihap |= value;  // Copy base address bits from MEIVT.
      regs_.at(meihapIx).poke(meihap);
    }

  // Poking mcounteren/scounteren changes accessibility of the
  // counters in user/supervisor modes.
  if (number == CsrNumber::MCOUNTEREN or number == CsrNumber::SCOUNTEREN)
    updateCounterPrivilege();

  return true;
}


template <typename URV>
bool
CsRegs<URV>::readTdata(CsrNumber number, PrivilegeMode mode, URV& value) const
{
  // Determine currently selected trigger.
  URV trigger = 0;
  if (not read(CsrNumber::TSELECT, mode, trigger))
    return false;

  if (number == CsrNumber::TDATA1)
    return triggers_.readData1(trigger, value);

  if (number == CsrNumber::TDATA2)
    return triggers_.readData2(trigger, value);

  if (number == CsrNumber::TDATA3)
    return triggers_.readData3(trigger, value);

  return false;
}


template <typename URV>
bool
CsRegs<URV>::writeTdata(CsrNumber number, PrivilegeMode mode, URV value)
{
  // Determine currently selected trigger.
  URV trigger = 0;
  if (not read(CsrNumber::TSELECT, mode, trigger))
    return false;

  // The CSR instructions never execute in debug mode.
  bool dMode = false;
  if (number == CsrNumber::TDATA1)
    {
      bool ok = triggers_.writeData1(trigger, dMode, value);
      if (ok) 
	{
	  // TDATA1 modified, update cached values
	  hasActiveTrigger_ = triggers_.hasActiveTrigger();
	  hasActiveInstTrigger_ = triggers_.hasActiveInstTrigger();
	}
      return ok;
    }

  if (number == CsrNumber::TDATA2)
    return triggers_.writeData2(trigger, dMode, value);

  if (number == CsrNumber::TDATA3)
    return triggers_.writeData3(trigger, dMode, value);

  return false;
}


template <typename URV>
bool
CsRegs<URV>::pokeTdata(CsrNumber number, URV value)
{
  // Determine currently selected trigger.
  URV trigger = 0;
  if (not read(CsrNumber::TSELECT, PrivilegeMode::Machine, trigger))
    return false;

  if (number == CsrNumber::TDATA1)
    {
      bool ok = triggers_.pokeData1(trigger, value);
      if (ok) 
	{
	  // TDATA1 modified, update cached values
	  hasActiveTrigger_ = triggers_.hasActiveTrigger();
	  hasActiveInstTrigger_ = triggers_.hasActiveInstTrigger();
	}
      return ok;
    }

  if (number == CsrNumber::TDATA2)
    return triggers_.pokeData2(trigger,value);

  if (number == CsrNumber::TDATA3)
    return triggers_.pokeData3(trigger, value);

  return false;
}


template <typename URV>
unsigned
CsRegs<URV>::getPmpConfigByteFromPmpAddr(CsrNumber csrn) const
{
  if (csrn < CsrNumber::PMPADDR0 or csrn > CsrNumber::PMPADDR15)
    return 0;

  unsigned pmpIx = unsigned(csrn) - unsigned(CsrNumber::PMPADDR0);

  // Determine rank of config register corresponding to pmpIx.
  unsigned cfgOffset = pmpIx / 4; // 0, 1, 2, or 3.

  // Identify byte within config register.
  unsigned byteIx = pmpIx % 4;

  if (not rv32_)
    {
      cfgOffset = (cfgOffset / 2) * 2;  // 0 or 2
      byteIx = pmpIx % 8;
    }

  CsrNumber cfgNum = CsrNumber(unsigned(CsrNumber::PMPCFG0) + cfgOffset);

  URV val = 0;
  if (peek(cfgNum, val))
    return (val >> 8*byteIx) & 0xff;

  return 0;
}


template <typename URV>
URV
CsRegs<URV>::adjustPmpValue(CsrNumber csrn, URV value) const
{
  if (csrn < CsrNumber::PMPADDR0 or csrn > CsrNumber::PMPADDR15)
    return value;   // Not a PMPADDR CSR.

  if (pmpG_ == 0)
    return value;

  unsigned byte = getPmpConfigByteFromPmpAddr(csrn);

  unsigned aField =(byte >> 3) & 3;
  if (aField < 2)
    {
      // A field is OFF or TOR
      if (pmpG_ >= 1)
        value = (value >> pmpG_) << pmpG_; // Clear least sig G bits.
    }
  else
    {
      // A field is NAPOT
      if (pmpG_ >= 2)
        {
          unsigned width = rv32_ ? 32 : 64;
          URV mask = ~URV(0) >> (width - pmpG_ + 1);
          value = value | mask; // Set to 1 least sig G-1 bits
        }
    }

  return value;
}


template <typename URV>
URV
CsRegs<URV>::legalizePmpcfgValue(URV current, URV value) const
{
  URV legal = 0;
  for (unsigned i = 0; i < sizeof(value); ++i)
    {
      uint8_t cb = (current >> (i*8)) & 0xff;  // Current byte.
      uint8_t nb = (value >> (i*8)) & 0xff;    // New byte.

      if (cb >> 7)
        nb = cb; // Field is locked. Use byte from current value.
      else if (pmpG_ != 0)
        {
          // If G is >= 1 then NA4 is not selectable in the A field of
          // the new byte.
          unsigned aField = (nb >> 3) & 3;
          if (aField == 2)
            {
              aField = 3;
              nb = nb | (aField << 3); // Change A field in new byte to 3.
            }
        }

      // w=1 r=0 is not allowed, change to w=0 r=0 
      if ((nb & 3) == 2)
        nb = (nb >> 2) << 2;

      legal = legal | (URV(nb) << i*8);
    }

  return legal;
}  


template <typename URV>
bool
CsRegs<URV>::isPmpaddrLocked(CsrNumber csrn) const
{
  if (csrn < CsrNumber::PMPADDR0 or csrn > CsrNumber::PMPADDR15)
    return false;   // Not a PMPADDR CSR.

  unsigned byte = getPmpConfigByteFromPmpAddr(csrn);
  bool locked = byte & 0x80;
  if (locked)
    return true;

  // If the next PMPADDR is top-of-range and is locked, then the
  // current PMADDR is considered to be locked.
  if (csrn >= CsrNumber::PMPADDR15)
    return false;  // No next PMPADDR register.

  CsrNumber csrn2 = CsrNumber(unsigned(csrn) + 1);
  byte = getPmpConfigByteFromPmpAddr(csrn2);
  locked = byte & 0x80;
  bool tor = ((byte >> 3) & 3) == 1;
  return locked and tor;
}


template <typename URV>
void
CsRegs<URV>::updateCounterPrivilege()
{
  URV mMask = 0;
  if (not peek(CsrNumber::MCOUNTEREN, mMask))
    return;

  URV sMask = 0;
  peek(CsrNumber::SCOUNTEREN, sMask);

  // Bits 0, 1, 2, 3 to 31 of mask correspond to CYCLE, TIME, INSTRET,
  // HPMCOUNTER3 to HPMCOUNTER31
  for (unsigned i = 0; i < 32; ++i)
    {
      bool mFlag = (mMask >> i) & 1;
      PrivilegeMode nextMode = PrivilegeMode::Machine;

      if (mFlag)
        {
          if (supervisorModeEnabled_)
            {
              nextMode = PrivilegeMode::Supervisor;

              bool sFlag = (sMask >> i) & 1;
              if (sFlag and userModeEnabled_)
                nextMode = PrivilegeMode::User;
            }
          else if (userModeEnabled_)
            nextMode = PrivilegeMode::User;
        }

      unsigned num = i + unsigned(CsrNumber::CYCLE);

      CsrNumber csrn = CsrNumber(num);
      auto csr = getImplementedCsr(csrn);
      if (csr)
        csr->setPrivilegeMode(nextMode);

      num = i + unsigned(CsrNumber::CYCLEH);
      csrn = CsrNumber(num);
      csr = getImplementedCsr(csrn);
      if (csr)
        csr->setPrivilegeMode(nextMode);
    }
}



template <typename URV>
URV
CsRegs<URV>::legalizeMhpmevent(CsrNumber number, URV value)
{
  bool enableUser = true;
  bool enableMachine = true;
  URV event = value;

  if (perModeCounterControl_)
    {
      enableUser = ! ((value >> 16) & 1);
      enableMachine = ! ((value >> 19) & 1);
      event = value & URV(0xffff);
    }

  if (hasPerfEventSet_)
    {
      if (not perfEventSet_.count(event))
        event = 0;
    }
  else
    event = std::min(event, maxEventId_);

  if (perModeCounterControl_)
    value = (value & ~URV(0xffff)) | event;
  else
    value = event;

  unsigned counterIx = unsigned(number) - unsigned(CsrNumber::MHPMEVENT3);
  assignEventToCounter(event, counterIx, enableUser, enableMachine);

  return value;
}


template class WdRiscv::CsRegs<uint32_t>;
template class WdRiscv::CsRegs<uint64_t>;
