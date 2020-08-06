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

#include <iomanip>
#include <iostream>
#include <sstream>
#include <cfenv>
#include <cmath>
#include <climits>
#include <map>
#include <mutex>
#include <array>
#include <boost/format.hpp>
#include <emmintrin.h>
#include <sys/time.h>
#include <poll.h>

// On pure 32-bit machines, use boost for 128-bit integer type.
#if __x86_64__
  typedef __int128_t Int128;
  typedef __uint128_t Uint128;
#else
  #include <boost/multiprecision/cpp_int.hpp>
  boost::multiprecision::int128_t Int128;
  boost::multiprecision::uint128_t Uint128;
#endif

#include <cstring>
#include <ctime>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/stat.h>

#include <cassert>
#include <csignal>

#define __STDC_FORMAT_MACROS
#include <cinttypes>
#include <sys/socket.h>
#include <netinet/in.h>

#include "instforms.hpp"
#include "DecodedInst.hpp"
#include "Hart.hpp"
#include "System.hpp"

using namespace WdRiscv;


template <typename TYPE>
static
bool
parseNumber(const std::string& numberStr, TYPE& number)
{
  bool good = not numberStr.empty();

  if (good)
    {
      char* end = nullptr;
      if (sizeof(TYPE) == 4)
	number = strtoul(numberStr.c_str(), &end, 0);
      else if (sizeof(TYPE) == 8)
	number = strtoull(numberStr.c_str(), &end, 0);
      else
	{
	  std::cerr << "parseNumber: Only 32/64-bit RISCV harts supported\n";
	  return false;
	}
      if (end and *end)
	good = false;  // Part of the string are non parseable.
    }
  return good;
}


/// Unsigned-float union: reinterpret bits as uint32_t or float
union Uint32FloatUnion
{
  Uint32FloatUnion(uint32_t u) : u(u)
  { }

  Uint32FloatUnion(float f) : f(f)
  { }

  uint32_t u = 0;
  float f;
};


/// Unsigned-float union: reinterpret bits as uint64_t or double
union Uint64DoubleUnion
{
  Uint64DoubleUnion(uint32_t u) : u(u)
  { }

  Uint64DoubleUnion(double d) : d(d)
  { }

  uint64_t u = 0;
  double d;
};



template <typename URV>
Hart<URV>::Hart(unsigned hartIx, Memory& memory)
  : hartIx_(hartIx), memory_(memory), intRegs_(32),
    fpRegs_(32), syscall_(*this),
    pmpManager_(memory.size(), memory.pageSize()),
    virtMem_(hartIx, memory, memory.pageSize(), 16 /* FIX: TLB size*/)
{
  regionHasLocalMem_.resize(16);
  regionHasLocalDataMem_.resize(16);
  regionHasDccm_.resize(16);
  regionHasMemMappedRegs_.resize(16);
  regionHasLocalInstMem_.resize(16);

  decodeCacheSize_ = 128*1024;  // Must be a power of 2.
  decodeCacheMask_ = decodeCacheSize_ - 1;
  decodeCache_.resize(decodeCacheSize_);

  interruptStat_.resize(size_t(InterruptCause::MAX_CAUSE) + 1);
  exceptionStat_.resize(size_t(ExceptionCause::MAX_CAUSE) + 1);
  for (auto& vec : exceptionStat_)
    vec.resize(size_t(SecondaryCause::MAX_CAUSE) + 1);

  // Tie the retired instruction and cycle counter CSRs to variables
  // held in the hart.
  if constexpr (sizeof(URV) == 4)
    {
      URV* low = reinterpret_cast<URV*> (&retiredInsts_);
      URV* high = low + 1;

      auto& mirLow = csRegs_.regs_.at(size_t(CsrNumber::MINSTRET));
      auto& irLow = csRegs_.regs_.at(size_t(CsrNumber::INSTRET));
      mirLow.tie(low);
      irLow.tie(low);

      auto& mirHigh = csRegs_.regs_.at(size_t(CsrNumber::MINSTRETH));
      auto& irHigh = csRegs_.regs_.at(size_t(CsrNumber::INSTRETH));
      mirHigh.tie(high);
      irHigh.tie(high);

      low = reinterpret_cast<URV*> (&cycleCount_);
      high = low + 1;

      auto& mcycleLow = csRegs_.regs_.at(size_t(CsrNumber::MCYCLE));
      auto& cycleLow = csRegs_.regs_.at(size_t(CsrNumber::CYCLE));
      mcycleLow.tie(low);
      cycleLow.tie(low);

      auto& mcycleHigh = csRegs_.regs_.at(size_t(CsrNumber::MCYCLEH));
      auto& cycleHigh = csRegs_.regs_.at(size_t(CsrNumber::CYCLEH));
      mcycleHigh.tie(high);
      cycleHigh.tie(high);

      // TIME is a read-only shadow of MCYCLE.
      csRegs_.regs_.at(size_t(CsrNumber::TIME)).tie(low);
      csRegs_.regs_.at(size_t(CsrNumber::TIMEH)).tie(high);
    }
  else
    {
      csRegs_.regs_.at(size_t(CsrNumber::MINSTRET)).tie(&retiredInsts_);
      csRegs_.regs_.at(size_t(CsrNumber::MCYCLE)).tie(&cycleCount_);

      // INSTRET and CYCLE are read-only shadows of MINSTRET and MCYCLE.
      csRegs_.regs_.at(size_t(CsrNumber::INSTRET)).tie(&retiredInsts_);
      csRegs_.regs_.at(size_t(CsrNumber::CYCLE)).tie(&cycleCount_);

      // TIME is a read-only shadow of MCYCLE.
      csRegs_.regs_.at(size_t(CsrNumber::TIME)).tie(&cycleCount_);
    }

  // Tie the FCSR register to variable held in the hart.
  csRegs_.regs_.at(size_t(CsrNumber::FCSR)).tie(&fcsrValue_);

  // Add local hart-id to the base-hart-id when MHARTID CSR is
  // constructed.  This will be over-written if MHARTID is configured
  // in the JSON configuration file.
  bool implemented = true, debug = false, shared = false;
  URV base = 0, reset = 0, mask = 0, pokeMask = 0;

  peekCsr(CsrNumber::MHARTID, base, reset, mask, pokeMask);
  URV hartId = base + hartIx;
  csRegs_.configCsr(CsrNumber::MHARTID, implemented, hartId, mask, pokeMask,
                    debug, shared);
}


template <typename URV>
Hart<URV>::~Hart()
{
}


template <typename URV>
void
Hart<URV>::getImplementedCsrs(std::vector<CsrNumber>& vec) const
{
  vec.clear();

  for (unsigned i = 0; i <= unsigned(CsrNumber::MAX_CSR_); ++i)
    {
      CsrNumber csrn = CsrNumber(i);
      if (csRegs_.isImplemented(csrn))
	vec.push_back(csrn);
    }
}


template <typename URV>
bool
Hart<URV>::configureCache(uint64_t size, unsigned lineSize,
                          unsigned setSize)
{
  return memory_.configureCache(size, lineSize, setSize);
}


template <typename URV>
void
Hart<URV>::deleteCache()
{
  memory_.deleteCache();
}


template <typename URV>
void
Hart<URV>::getCacheLineAddresses(std::vector<uint64_t>& addresses)
{
  memory_.getCacheLineAddresses(addresses);
}


template <typename URV>
unsigned
Hart<URV>::countImplementedPmpRegisters() const
{
  unsigned count = 0;

  unsigned num = unsigned(CsrNumber::PMPADDR0);
  for (unsigned ix = 0; ix < 16; ++ix, ++num)
    if (csRegs_.isImplemented(CsrNumber(num)))
      count++;

  if (count and count < 16)
    std::cerr << "Warning: Some but not all PMPADDR CSRs are implemented\n";

  unsigned cfgCount = 0;
  if (mxlen_ == 32)
    {
      num = unsigned(CsrNumber::PMPCFG0);
      for (unsigned ix = 0; ix < 4; ++ix, ++num)
        if (csRegs_.isImplemented(CsrNumber(num)))
          cfgCount++;
      if (count and cfgCount != 4)
        std::cerr << "Warning: Physical memory protection enabled but not all "
                  << "of the config register (PMPCFG) are implemented\n";
    }
  else
    {
      num = unsigned(CsrNumber::PMPCFG0);
      for (unsigned ix = 0; ix < 2; ++ix, num += 2)
        if (csRegs_.isImplemented(CsrNumber(num)))
          cfgCount++;
      if (count and cfgCount != 2)
        std::cerr << "Warning: Physical memory protection enabled but not all "
                  << "of the config register (PMPCFG) are implemented\n";
    }

  return count;
}
      


template <typename URV>
void
Hart<URV>::processExtensions()
{
  // D requires F and is enabled only if F is enabled.
  rva_ = false;
  rvc_ = false;
  rvd_ = false;
  rve_ = false;
  rvf_ = false;
  rvm_ = false;
  rvs_ = false;
  rvu_ = false;

  URV value = 0;
  if (peekCsr(CsrNumber::MISA, value))
    {
      if (value & 1)    // Atomic ('a') option.
	rva_ = true;

      if (value & (URV(1) << ('c' - 'a')))  // Compress option.
	rvc_ = true;

      if (value & (URV(1) << ('f' - 'a')))  // Single precision FP
	{
	  rvf_ = true;

	  bool isDebug = false, shared = true;

	  // Make sure FCSR/FRM/FFLAGS are enabled if F extension is on.
	  if (not csRegs_.isImplemented(CsrNumber::FCSR))
	    csRegs_.configCsr("fcsr", true, 0, 0xff, 0xff, isDebug, shared);
	  if (not csRegs_.isImplemented(CsrNumber::FRM))
	    csRegs_.configCsr("frm", true, 0, 0x7, 0x7, isDebug, shared);
	  if (not csRegs_.isImplemented(CsrNumber::FFLAGS))
	    csRegs_.configCsr("fflags", true, 0, 0x1f, 0x1f, isDebug, shared);
	}

      if (value & (URV(1) << ('d' - 'a')))  // Double precision FP.
	{
	  if (rvf_)
	    rvd_ = true;
	  else
	    std::cerr << "Bit 3 (d) is set in the MISA register but f "
		      << "extension (bit 5) is not enabled -- ignored\n";
	}

      if (value & (URV(1) << ('e' - 'a')))  // Double precision FP.
        {
          rve_ = true;
          intRegs_.regs_.resize(16);
        }

      if (not (value & (URV(1) << ('i' - 'a'))))
	std::cerr << "Bit 8 (i extension) is cleared in the MISA register "
		  << " but extension is mandatory -- assuming bit 8 set\n";

      if (value & (URV(1) << ('m' - 'a')))  // Multiply/divide option.
	rvm_ = true;

      if (value & (URV(1) << ('u' - 'a')))  // User-mode option.
        enableUserMode(true);

      if (value & (URV(1) << ('s' - 'a')))  // Supervisor-mode option.
        enableSupervisorMode(true);

      for (auto ec : { 'b', 'g', 'h', 'j', 'k', 'l', 'n', 'o', 'p',
	    'q', 'r', 't', 'v', 'w', 'x', 'y', 'z' } )
	{
	  unsigned bit = ec - 'a';
	  if (value & (URV(1) << bit))
	    std::cerr << "Bit " << bit << " (" << ec << ") set in the MISA "
		      << "register but extension is not supported "
		      << "-- ignored\n";
	}
    }
}


template <typename URV>
void
Hart<URV>::reset(bool resetMemoryMappedRegs)
{
  privMode_ = PrivilegeMode::Machine;

  intRegs_.reset();
  csRegs_.reset();

  // Suppress resetting memory mapped register on initial resets sent
  // by the test bench. Otherwise, initial resets obliterate memory
  // mapped register data loaded from the ELF file.
  if (resetMemoryMappedRegs)
    memory_.resetMemoryMappedRegisters();
  memory_.invalidateLr(hartIx_);

  clearPendingNmi();
  clearTraceData();

  loadQueue_.clear();

  pc_ = resetPc_;
  currPc_ = resetPc_;

  // Enable extensions if corresponding bits are set in the MISA CSR.
  processExtensions();
  
  perfControl_ = ~uint32_t(0);
  URV value = 0;
  if (peekCsr(CsrNumber::MCOUNTINHIBIT, value))
    perfControl_ = ~value;

  prevPerfControl_ = perfControl_;

  debugMode_ = false;
  debugStepMode_ = false;

  dcsrStepIe_ = false;
  dcsrStep_ = false;

  if (peekCsr(CsrNumber::DCSR, value))
    {
      dcsrStep_ = (value >> 2) & 1;
      dcsrStepIe_ = (value >> 11) & 1;
    }

  updateStackChecker();  // Swerv-specific feature.
  enableWideLdStMode(false);  // Swerv-specific feature.

  hartStarted_ = true;

  // If mhartstart exists then use its bits to decide which hart has
  // started.
  if (hartIx_ != 0)
    {
      auto csr = findCsr("mhartstart");
      if (csr)
        {
          URV value = csr->read();
          hartStarted_ = ((URV(1) << hartIx_) & value) != 0;
        }
    }

  // Enable FP if f/d extension and linux/newlib.
  if ((isRvf() or isRvd()) and (newlib_ or linux_))
    {
      URV val = csRegs_.peekMstatus();
      MstatusFields<URV> fields(val);
      fields.bits_.FS = unsigned(FpFs::Initial);
      csRegs_.write(CsrNumber::MSTATUS, PrivilegeMode::Machine, fields.value_);
    }

  // Update cached values of mstatus.mpp and mstatus.mprv and mstatus.fs.
  updateCachedMstatusFields();

  updateAddressTranslation();

  updateMemoryProtection();
  countImplementedPmpRegisters();

  csRegs_.updateCounterPrivilege();

  alarmLimit_ = alarmInterval_? alarmInterval_ + instCounter_ : ~uint64_t(0);
}


template <typename URV>
void
Hart<URV>::updateCachedMstatusFields()
{
  URV csrVal = csRegs_.peekMstatus();
  MstatusFields<URV> msf(csrVal);
  mstatusMpp_ = PrivilegeMode(msf.bits_.MPP);
  mstatusMprv_ = msf.bits_.MPRV;
  mstatusFs_ = FpFs(msf.bits_.FS);

  virtMem_.setExecReadable(msf.bits_.MXR);
  virtMem_.setSupervisorAccessUser(msf.bits_.SUM);
}


template <typename URV>
bool
Hart<URV>::loadHexFile(const std::string& file)
{
  return memory_.loadHexFile(file);
}


template <typename URV>
bool
Hart<URV>::loadElfFile(const std::string& file, size_t& entryPoint)
{
  unsigned registerWidth = sizeof(URV)*8;

  size_t end = 0;
  if (not memory_.loadElfFile(file, registerWidth, entryPoint, end))
    return false;

  this->pokePc(URV(entryPoint));

  ElfSymbol sym;

  if (not toHostSym_.empty() and this->findElfSymbol(toHostSym_, sym))
    this->setToHostAddress(sym.addr_);

  if (not consoleIoSym_.empty() and this->findElfSymbol(consoleIoSym_, sym))
    this->setConsoleIo(URV(sym.addr_));

  if (this->findElfSymbol("__global_pointer$", sym))
    this->pokeIntReg(RegGp, URV(sym.addr_));

  if (this->findElfSymbol("_end", sym))   // For newlib/linux emulation.
    this->setTargetProgramBreak(URV(sym.addr_));
  else
    this->setTargetProgramBreak(URV(end));

  return true;
}


template <typename URV>
bool
Hart<URV>::peekMemory(size_t address, uint8_t& val) const
{
  if (memory_.read(address, val))
    return true;

  // We may have failed because location is in instruction space.
  return memory_.readInst(address, val);
}
  

template <typename URV>
bool
Hart<URV>::peekMemory(size_t address, uint16_t& val) const
{
  if (memory_.read(address, val))
    return true;

  // We may have failed because location is in instruction space.
  return memory_.readInst(address, val);
}


template <typename URV>
bool
Hart<URV>::peekMemory(size_t address, uint32_t& val) const
{
  if (memory_.read(address, val))
    return true;

  // We may have failed because location is in instruction space.
  return memory_.readInst(address, val);
}


template <typename URV>
bool
Hart<URV>::peekMemory(size_t address, uint64_t& val) const
{
  uint32_t high = 0, low = 0;

  if (memory_.read(address, low) and memory_.read(address + 4, high))
    {
      val = (uint64_t(high) << 32) | low;
      return true;
    }

  // We may have failed because location is in instruction space.
  if (memory_.readInst(address, low) and
      memory_.readInst(address + 4, high))
    {
      val = (uint64_t(high) << 32) | low;
      return true;
    }

  return true;
}


template <typename URV>
bool
Hart<URV>::pokeMemory(size_t addr, uint8_t val)
{
  std::lock_guard<std::mutex> lock(memory_.lrMutex_);

  memory_.invalidateLrs(addr, sizeof(val));

  if (memory_.poke(addr, val))
    {
      invalidateDecodeCache(addr, sizeof(val));
      return true;
    }

  return false;
}


template <typename URV>
bool
Hart<URV>::pokeMemory(size_t addr, uint16_t val)
{
  std::lock_guard<std::mutex> lock(memory_.lrMutex_);

  memory_.invalidateLrs(addr, sizeof(val));

  if (memory_.poke(addr, val))
    {
      invalidateDecodeCache(addr, sizeof(val));
      return true;
    }

  return false;
}


template <typename URV>
bool
Hart<URV>::pokeMemory(size_t addr, uint32_t val)
{
  // We allow poke to bypass masking for memory mapped registers
  // otherwise, there is no way for external driver to clear bits that
  // are read-only to this hart.

  std::lock_guard<std::mutex> lock(memory_.lrMutex_);

  memory_.invalidateLrs(addr, sizeof(val));

  if (memory_.poke(addr, val))
    {
      invalidateDecodeCache(addr, sizeof(val));
      return true;
    }

  return false;
}


template <typename URV>
bool
Hart<URV>::pokeMemory(size_t addr, uint64_t val)
{
  std::lock_guard<std::mutex> lock(memory_.lrMutex_);

  memory_.invalidateLrs(addr, sizeof(val));

  if (memory_.poke(addr, val))
    {
      invalidateDecodeCache(addr, sizeof(val));
      return true;
    }

  return false;
}


template <typename URV>
void
Hart<URV>::setPendingNmi(NmiCause cause)
{
  // First nmi sets the cause. The cause is sticky.
  if (not nmiPending_)
    nmiCause_ = cause;

  nmiPending_ = true;

  // Set the nmi pending bit in the DCSR register.
  URV val = 0;  // DCSR value
  if (peekCsr(CsrNumber::DCSR, val))
    {
      val |= URV(1) << 3;  // nmip bit
      pokeCsr(CsrNumber::DCSR, val);
      recordCsrWrite(CsrNumber::DCSR);
    }
}


template <typename URV>
void
Hart<URV>::clearPendingNmi()
{
  nmiPending_ = false;
  nmiCause_ = NmiCause::UNKNOWN;

  URV val = 0;  // DCSR value
  if (peekCsr(CsrNumber::DCSR, val))
    {
      val &= ~(URV(1) << 3);  // nmip bit
      pokeCsr(CsrNumber::DCSR, val);
      recordCsrWrite(CsrNumber::DCSR);
    }
}


template <typename URV>
void
Hart<URV>::setToHostAddress(size_t address)
{
  toHost_ = URV(address);
  toHostValid_ = true;
}


template <typename URV>
void
Hart<URV>::clearToHostAddress()
{
  toHost_ = 0;
  toHostValid_ = false;
}


template <typename URV>
void
Hart<URV>::putInLoadQueue(unsigned size, size_t addr, unsigned regIx,
			  uint64_t data, bool isWide)
{
  if (not loadQueueEnabled_)
    return;

  if (isAddrInDccm(addr))
    {
      // Blocking load. Invalidate target register in load queue so
      // that it will not be reverted.
      invalidateInLoadQueue(regIx, false);
      return;
    }

  if (loadQueue_.size() >= maxLoadQueueSize_)
    {
      std::cerr << "At #" << instCounter_ << ": Load queue full.\n";
      for (size_t i = 1; i < maxLoadQueueSize_; ++i)
	loadQueue_[i-1] = loadQueue_[i];
      loadQueue_[maxLoadQueueSize_-1] = LoadInfo(size, addr, regIx, data,
						 isWide, instCounter_);
    }
  else
    loadQueue_.push_back(LoadInfo(size, addr, regIx, data, isWide,
				  instCounter_));
}


template <typename URV>
void
Hart<URV>::invalidateInLoadQueue(unsigned regIx, bool isDiv)
{
  if (regIx == lastDivRd_ and not isDiv)
    hasLastDiv_ = false;

  // Replace entry containing target register with x0 so that load exception
  // matching entry will not revert target register.
  for (unsigned i = 0; i < loadQueue_.size(); ++i)
    if (loadQueue_[i].regIx_ == regIx)
      loadQueue_[i].makeInvalid();
}


template <typename URV>
void
Hart<URV>::removeFromLoadQueue(unsigned regIx, bool isDiv)
{
  if (regIx == 0)
    return;

  if (regIx == lastDivRd_ and not isDiv)
    hasLastDiv_ = false;

  // Last (most recent) matching entry is removed. Subsequent entries
  // are invalidated.
  bool last = true;
  size_t removeIx = loadQueue_.size();
  for (size_t i = loadQueue_.size(); i > 0; --i)
    {
      auto& entry = loadQueue_.at(i-1);
      if (not entry.isValid())
	continue;
      if (entry.regIx_ == regIx)
	{
	  if (last)
	    {
	      removeIx = i-1;
	      last = false;
	    }
	  else
	    entry.makeInvalid();
	}
    }

  if (removeIx < loadQueue_.size())
    loadQueue_.erase(loadQueue_.begin() + removeIx);
}


template <typename URV>
inline
void
Hart<URV>::execBeq(const DecodedInst* di)
{
  uint32_t rs1 = di->op0();
  uint32_t rs2 = di->op1();
  if (intRegs_.read(rs1) != intRegs_.read(rs2))
    return;
  pc_ = currPc_ + di->op2As<SRV>();
  pc_ = (pc_ >> 1) << 1;  // Clear least sig bit.
  lastBranchTaken_ = true;
}


template <typename URV>
inline
void
Hart<URV>::execBne(const DecodedInst* di)
{
  if (intRegs_.read(di->op0()) == intRegs_.read(di->op1()))
    return;
  pc_ = currPc_ + di->op2As<SRV>();
  pc_ = (pc_ >> 1) << 1;  // Clear least sig bit.
  lastBranchTaken_ = true;
}


template <typename URV>
inline
void
Hart<URV>::execAddi(const DecodedInst* di)
{
  SRV imm = di->op2As<SRV>();
  SRV v = intRegs_.read(di->op1()) + imm;
  intRegs_.write(di->op0(), v);
}


template <typename URV>
inline
void
Hart<URV>::execAdd(const DecodedInst* di)
{
  URV v = intRegs_.read(di->op1()) + intRegs_.read(di->op2());
  intRegs_.write(di->op0(), v);
}


template <typename URV>
inline
void
Hart<URV>::execAndi(const DecodedInst* di)
{
  SRV imm = di->op2As<SRV>();
  URV v = intRegs_.read(di->op1()) & imm;
  intRegs_.write(di->op0(), v);
}


template <typename URV>
bool
Hart<URV>::isIdempotentRegion(size_t addr) const
{
  unsigned region = unsigned(addr >> (sizeof(URV)*8 - 4));
  URV mracVal = 0;
  if (csRegs_.read(CsrNumber::MRAC, PrivilegeMode::Machine, mracVal))
    {
      unsigned bit = (mracVal >> (region*2 + 1)) & 1;
      return bit == 0  or regionHasLocalMem_.at(region);
    }
  return true;
}


template <typename URV>
bool
Hart<URV>::applyStoreException(URV addr, unsigned& matches)
{
  if (not nmiPending_)
    {
      bool prevLocked = csRegs_.mdseacLocked();
      if (not prevLocked)
        {
          pokeCsr(CsrNumber::MDSEAC, addr); // MDSEAC is read only: Poke it.
          csRegs_.lockMdseac(true);
          setPendingNmi(NmiCause::STORE_EXCEPTION);
        }
    }
  recordCsrWrite(CsrNumber::MDSEAC); // Always record change (per Ajay Nath)

  matches = 1;
  return true;
}


template <typename URV>
bool
Hart<URV>::applyLoadException(URV addr, unsigned tag, unsigned& matches)
{
  // if (not isNmiEnabled())
  //   return false;  // NMI should not have been delivered to this hart.

  if (not nmiPending_)
    {
      bool prevLocked = csRegs_.mdseacLocked();
      if (not prevLocked)
        {
          pokeCsr(CsrNumber::MDSEAC, addr); // MDSEAC is read only: Poke it.
          csRegs_.lockMdseac(true);
          setPendingNmi(NmiCause::LOAD_EXCEPTION);
        }
    }
  recordCsrWrite(CsrNumber::MDSEAC);  // Always record change (per Ajay Nath)

  if (not loadErrorRollback_)
    {
      matches = 1;
      return true;
    }

  // Count matching records. Determine if there is an entry with the
  // same register as the first match but younger.
  bool hasYounger = false;
  unsigned targetReg = 0;  // Register of 1st match.
  matches = 0;
  unsigned iMatches = 0;  // Invalid matching entries.
  for (const LoadInfo& li : loadQueue_)
    {
      if (matches and li.isValid() and targetReg == li.regIx_)
	hasYounger = true;

      if (li.tag_ == tag)
	{
	  if (li.isValid())
	    {
	      targetReg = li.regIx_;
	      matches++;
	    }
	  else
	    iMatches++;
	}
    }

  matches += iMatches;
  if (matches != 1)
    {
      std::cerr << "Error: Load exception addr:0x" << std::hex << addr << std::dec;
      std::cerr << " tag:" << tag;
      if (matches == 0)
	std::cerr << " does not match any entry in the load queue\n";
      else
	std::cerr << " matches " << matches << " entries"
		  << " in the load queue\n";
      return false;
    }

  // Revert register of matching item unless there are younger entries
  // with same resister. Revert with value of older entry with same
  // target register (if multiple such entry, use oldest). Invalidate
  // all older entries with same target. Remove item from queue.
  // Update prev-data of 1st younger item with same target register.
  size_t removeIx = loadQueue_.size();
  for (size_t ix = 0; ix < loadQueue_.size(); ++ix)
    {
      auto& entry = loadQueue_.at(ix);
      if (entry.tag_ == tag)
	{
	  removeIx = ix;
	  if (not entry.isValid())
	    continue;
	}
      else
	continue;

      removeIx = ix;

      URV prev = entry.prevData_;

      // Revert to oldest entry with same target reg. Invalidate older
      // entries with same target reg.
      for (size_t ix2 = removeIx; ix2 > 0; --ix2)
	{
	  auto& entry2 = loadQueue_.at(ix2-1);
	  if (entry2.isValid() and entry2.regIx_ == entry.regIx_)
	    {
	      prev = entry2.prevData_;
	      entry2.makeInvalid();
	    }
	}

      if (not hasYounger)
	{
	  pokeIntReg(entry.regIx_, prev);
	  if (entry.wide_)
            pokeCsr(CsrNumber::MDBHD, entry.prevData_ >> 32);
	}

      // Update prev-data of 1st younger item with same target reg.
      for (size_t ix2 = removeIx + 1; ix2 < loadQueue_.size(); ++ix2)
	{
	  auto& entry2 = loadQueue_.at(ix2);
 	  if (entry2.isValid() and entry2.regIx_ == entry.regIx_)
 	    {
	      entry2.prevData_ = prev;
	      break;
	    }
	}

      break;
    }

  if (removeIx < loadQueue_.size())
    loadQueue_.erase(loadQueue_.begin() + removeIx);

  return true;
}


template <typename URV>
bool
Hart<URV>::applyLoadFinished(URV addr, unsigned tag, unsigned& matches)
{
  if (not loadErrorRollback_)
    {
      matches = 1;
      return true;
    }

  // Count matching records.
  matches = 0;
  size_t matchIx = 0;     // Index of matching entry.
  size_t size = loadQueue_.size();
  for (size_t i = 0; i < size; ++i)
    {
      const LoadInfo& li = loadQueue_.at(i);
      if (li.tag_ == tag)
	{
	  if (not matches)
	    matchIx = i;
	  matches++;
	}
    }

  if (matches == 0)
    {
      std::cerr << "Warning: Load finished addr:0x" << std::hex << addr << std::dec;
      std::cerr << " tag:" << tag << " does not match any entry in the load queue\n";
      return true;
    }

  if (matches > 1)
    {
      std::cerr << "Warning: Load finished at 0x" << std::hex << addr << std::dec;
      std::cerr << " matches multiple intries in the load queue\n";
    }

  LoadInfo& entry = loadQueue_.at(matchIx);

  // Process entries in reverse order (start with oldest)
  // Mark all earlier entries with same target register as invalid.
  // Identify earliest previous value of target register.
  unsigned targetReg = entry.regIx_;
  size_t prevIx = matchIx;
  uint64_t prev = entry.prevData_;  // Previous value of target reg.
  for (size_t j = 0; j < matchIx; ++j)
    {
      LoadInfo& li = loadQueue_.at(j);
      if (not li.isValid())
	continue;
      if (li.regIx_ != targetReg)
	continue;

      li.makeInvalid();
      if (j < prevIx)
	{
	  prevIx = j;
	  prev = li.prevData_;
	}
    }

  // Update prev-data of 1st subsequent entry with same target.
  if (entry.isValid())
    for (size_t j = matchIx + 1; j < size; ++j)
      {
	LoadInfo& li = loadQueue_.at(j);
	if (li.isValid() and li.regIx_ == targetReg)
	  {
	    // Preserve upper 32 bits if wide (64-bit) load.
	    if (li.wide_)
	      prev = ((prev << 32) >> 32) | ((li.prevData_ >> 32) << 32);
	    li.prevData_ = prev;
	    break;
	  }
      }

  // Remove matching entry from queue.
  size_t newSize = 0;
  for (size_t i = 0; i < size; ++i)
    {
      auto& li = loadQueue_.at(i);
      bool remove = i == matchIx;
      if (remove)
	continue;

      if (newSize != i)
	loadQueue_.at(newSize) = li;
      newSize++;
    }
  loadQueue_.resize(newSize);

  return true;
}


static
void
printUnsignedHisto(const char* tag, const std::vector<uint64_t>& histo,
		   FILE* file)
{
  if (histo.size() < 7)
    return;

  if (histo.at(0))
    fprintf(file, "    %s  0          %" PRId64 "\n", tag, histo.at(0));
  if (histo.at(1))
    fprintf(file, "    %s  1          %" PRId64 "\n", tag, histo.at(1));
  if (histo.at(2))
    fprintf(file, "    %s  2          %" PRId64 "\n", tag, histo.at(2));
  if (histo.at(3))
    fprintf(file, "    %s  (2,   16]  %" PRId64 "\n", tag, histo.at(3));
  if (histo.at(4))
    fprintf(file, "    %s  (16,  1k]  %" PRId64 "\n", tag, histo.at(4));
  if (histo.at(5))
    fprintf(file, "    %s  (1k, 64k]  %" PRId64 "\n", tag, histo.at(5));
  if (histo.at(6))
    fprintf(file, "    %s  > 64k      %" PRId64 "\n", tag, histo.at(6));
}


static
void
printSignedHisto(const char* tag, const std::vector<uint64_t>& histo,
		 FILE* file)
{
  if (histo.size() < 13)
    return;

  if (histo.at(0))
    fprintf(file, "    %s <= -64k     %" PRId64 "\n", tag, histo.at(0));
  if (histo.at(1))
    fprintf(file, "    %s (-64k, -1k] %" PRId64 "\n", tag, histo.at(1));
  if (histo.at(2))
    fprintf(file, "    %s (-1k,  -16] %" PRId64 "\n", tag, histo.at(2));
  if (histo.at(3))
    fprintf(file, "    %s (-16,   -3] %" PRId64 "\n", tag, histo.at(3));
  if (histo.at(4))
    fprintf(file, "    %s -2          %" PRId64 "\n", tag, histo.at(4));
  if (histo.at(5))
    fprintf(file, "    %s -1          %" PRId64 "\n", tag, histo.at(5));
  if (histo.at(6))
    fprintf(file, "    %s 0           %" PRId64 "\n", tag, histo.at(6));
  if (histo.at(7))
    fprintf(file, "    %s 1           %" PRId64 "\n", tag, histo.at(7));
  if (histo.at(8))
    fprintf(file, "    %s 2           %" PRId64 "\n", tag, histo.at(8));
  if (histo.at(9))
    fprintf(file, "    %s (2,     16] %" PRId64 "\n", tag, histo.at(9));
  if (histo.at(10))
    fprintf(file, "    %s (16,    1k] %" PRId64 "\n", tag, histo.at(10));
  if (histo.at(11))	              
    fprintf(file, "    %s (1k,   64k] %" PRId64 "\n", tag, histo.at(11));
  if (histo.at(12))	              
    fprintf(file, "    %s > 64k       %" PRId64 "\n", tag, histo.at(12));
}


template <typename URV>
void
Hart<URV>::reportInstructionFrequency(FILE* file) const
{
  struct CompareFreq
  {
    CompareFreq(const std::vector<InstProfile>& profileVec)
      : profileVec(profileVec)
    { }

    bool operator()(size_t a, size_t b) const
    { return profileVec.at(a).freq_ < profileVec.at(b).freq_; }

    const std::vector<InstProfile>& profileVec;
  };

  std::vector<size_t> indices(instProfileVec_.size());
  for (size_t i = 0; i < indices.size(); ++i)
    indices.at(i) = i;
  std::sort(indices.begin(), indices.end(), CompareFreq(instProfileVec_));

  for (size_t i = 0; i < indices.size(); ++i)
    {
      size_t ix = indices.at(i);
      InstId id = InstId(ix);

      const InstEntry& entry = instTable_.getEntry(id);
      const InstProfile& prof = instProfileVec_.at(ix);
      uint64_t freq = prof.freq_;
      if (not freq)
	continue;

      fprintf(file, "%s %" PRId64 "\n", entry.name().c_str(), freq);

      auto regCount = intRegCount();

      uint64_t count = 0;
      for (auto n : prof.rd_) count += n;
      if (count)
	{
	  fprintf(file, "  +rd");
	  for (unsigned i = 0; i < regCount; ++i)
	    if (prof.rd_.at(i))
	      fprintf(file, " %d:%" PRId64, i, prof.rd_.at(i));
	  fprintf(file, "\n");
	}

      uint64_t count1 = 0;
      for (auto n : prof.rs1_) count1 += n;
      if (count1)
	{
	  fprintf(file, "  +rs1");
	  for (unsigned i = 0; i < regCount; ++i)
	    if (prof.rs1_.at(i))
	      fprintf(file, " %d:%" PRId64, i, prof.rs1_.at(i));
	  fprintf(file, "\n");

	  const auto& histo = prof.rs1Histo_;
	  if (entry.isUnsigned())
	    printUnsignedHisto("+hist1", histo, file);
	  else
	    printSignedHisto("+hist1", histo, file);
	}

      uint64_t count2 = 0;
      for (auto n : prof.rs2_) count2 += n;
      if (count2)
	{
	  fprintf(file, "  +rs2");
	  for (unsigned i = 0; i < regCount; ++i)
	    if (prof.rs2_.at(i))
	      fprintf(file, " %d:%" PRId64, i, prof.rs2_.at(i));
	  fprintf(file, "\n");

	  const auto& histo = prof.rs2Histo_;
	  if (entry.isUnsigned())
	    printUnsignedHisto("+hist2", histo, file);
	  else
	    printSignedHisto("+hist2", histo, file);
	}

      if (prof.hasImm_)
	{
	  fprintf(file, "  +imm  min:%d max:%d\n", prof.minImm_, prof.maxImm_);
	  printSignedHisto("+hist ", prof.immHisto_, file);
	}

      if (prof.user_)
        fprintf(file, "  +user %ld\n", prof.user_);
      if (prof.supervisor_)
        fprintf(file, "  +supervisor %ld\n", prof.supervisor_);
      if (prof.machine_)
        fprintf(file, "  +machine %ld\n", prof.machine_);
    }
}


template <typename URV>
void
Hart<URV>::reportTrapStat(FILE* file) const
{
  fprintf(file, "\n");
  fprintf(file, "Interrupts (incuding NMI): %ld\n", interruptCount_);
  for (unsigned i = 0; i < interruptStat_.size(); ++i)
    {
      InterruptCause cause = InterruptCause(i);
      uint64_t count = interruptStat_.at(i);
      if (not count)
        continue;
      switch(cause)
        {
        case InterruptCause::U_SOFTWARE:
          fprintf(file, "  + U_SOFTWARE  : %ld\n", count);
          break;
        case InterruptCause::S_SOFTWARE:
          fprintf(file, "  + S_SOFTWARE  : %ld\n", count);
          break;
        case InterruptCause::M_SOFTWARE:
          fprintf(file, "  + M_SOFTWARE  : %ld\n", count);
          break;
        case InterruptCause::U_TIMER   :
          fprintf(file, "  + U_TIMER     : %ld\n", count);
          break;
        case InterruptCause::S_TIMER   :
          fprintf(file, "  + S_TIMER     : %ld\n", count);
          break;
        case InterruptCause::M_TIMER   :
          fprintf(file, "  + M_TIMER     : %ld\n", count);
          break;
        case InterruptCause::U_EXTERNAL:
          fprintf(file, "  + U_EXTERNAL  : %ld\n", count);
          break;
        case InterruptCause::S_EXTERNAL:
          fprintf(file, "  + S_EXTERNAL  : %ld\n", count);
          break;
        case InterruptCause::M_EXTERNAL:
          fprintf(file, "  + M_EXTERNAL  : %ld\n", count);
          break;
        case InterruptCause::M_INT_TIMER1:
          fprintf(file, "  + M_INT_TIMER1: %ld\n", count);
          break;
        case InterruptCause::M_INT_TIMER0:
          fprintf(file, "  + M_INT_TIMER0: %ld\n", count);
          break;
        case InterruptCause::M_LOCAL :
          fprintf(file, "  + M_LOCAL     : %ld\n", count);
          break;
        default:
          fprintf(file, "  + ????        : %ld\n", count);
        }
    }

  fprintf(file, "\n");
  fprintf(file, "Non maskable interrupts: %ld\n", nmiCount_);

  fprintf(file, "\n");
  fprintf(file, "Exceptions: %ld\n", exceptionCount_);
  for (unsigned i = 0; i < exceptionStat_.size(); ++i)
    {
      ExceptionCause cause = ExceptionCause(i);
      uint64_t count = 0;
      const auto& secCauseVec = exceptionStat_.at(i);
      for (auto n : secCauseVec)
        count += n;

      if (not count)
        continue;

      switch(cause)
        {
        case ExceptionCause::INST_ADDR_MISAL :
          fprintf(file, "  + INST_ADDR_MISAL : %ld\n", count);
          break;
        case ExceptionCause::INST_ACC_FAULT  :
          fprintf(file, "  + INST_ACC_FAULT  : %ld\n", count);
          break;
        case ExceptionCause::ILLEGAL_INST    :
          fprintf(file, "  + ILLEGAL_INST    : %ld\n", count);
          break;
        case ExceptionCause::BREAKP          :
          fprintf(file, "  + BREAKP          : %ld\n", count);
          break;
        case ExceptionCause::LOAD_ADDR_MISAL :
          fprintf(file, "  + LOAD_ADDR_MISAL : %ld\n", count);
          break;
        case ExceptionCause::LOAD_ACC_FAULT  :
          fprintf(file, "  + LOAD_ACC_FAULT  : %ld\n", count);
          break;
        case ExceptionCause::STORE_ADDR_MISAL:
          fprintf(file, "  + STORE_ADDR_MISAL: %ld\n", count);
          break;
        case ExceptionCause::STORE_ACC_FAULT :
          fprintf(file, "  + STORE_ACC_FAULT : %ld\n", count);
          break;
        case ExceptionCause::U_ENV_CALL      :
          fprintf(file, "  + U_ENV_CALL      : %ld\n", count);
          break;
        case ExceptionCause::S_ENV_CALL      :
          fprintf(file, "  + S_ENV_CALL      : %ld\n", count);
          break;
        case ExceptionCause::M_ENV_CALL      :
          fprintf(file, "  + M_ENV_CALL      : %ld\n", count);
          break;
        case ExceptionCause::INST_PAGE_FAULT :
          fprintf(file, "  + INST_PAGE_FAULT : %ld\n", count);
          break;
        case ExceptionCause::LOAD_PAGE_FAULT :
          fprintf(file, "  + LOAD_PAGE_FAULT : %ld\n", count);
          break;
        case ExceptionCause::STORE_PAGE_FAULT:
          fprintf(file, "  + STORE_PAGE_FAULT: %ld\n", count);
          break;
        case ExceptionCause::NONE            :
          fprintf(file, "  + NONE            : %ld\n", count);
          break;
        default:
          fprintf(file, "  + ????            : %ld\n", count);
          break;
        }
      for (unsigned j = 0; j < secCauseVec.size(); ++j)
        {
          uint64_t secCount = secCauseVec.at(j);
          if (secCount)
            fprintf(file, "    + %d: %ld\n", j, secCount);
        }
    }
}


template <typename URV>
void
Hart<URV>::reportPmpStat(FILE* file) const
{
  std::ostringstream oss;
  pmpManager_.printStats(oss);
  fprintf(file, "%s", oss.str().c_str());
}


template <typename URV>
ExceptionCause
Hart<URV>::determineMisalLoadException(URV addr, unsigned accessSize,
                                       SecondaryCause& secCause) const
{
  if (not misalDataOk_)
    {
      secCause = SecondaryCause::NONE;
      return ExceptionCause::LOAD_ADDR_MISAL;
    }

  size_t addr2 = addr + accessSize - 1;

  // Misaligned access to PIC.
  if (isAddrMemMapped(addr))
    {
      secCause = SecondaryCause::LOAD_ACC_PIC;
      return ExceptionCause::LOAD_ACC_FAULT;
    }

  // Crossing 256 MB region boundary.
  if (memory_.getRegionIndex(addr) != memory_.getRegionIndex(addr2))
    {
      secCause = SecondaryCause::LOAD_MISAL_REGION_CROSS;
      return ExceptionCause::LOAD_ADDR_MISAL;
    }

  // Misaligned access to a region with side effect.
  if (not isIdempotentRegion(addr) or not isIdempotentRegion(addr2))
    {
      secCause = SecondaryCause::LOAD_MISAL_IO;
      return ExceptionCause::LOAD_ADDR_MISAL;
    }

  secCause = SecondaryCause::NONE;
  return ExceptionCause::NONE;
}


template <typename URV>
ExceptionCause
Hart<URV>::determineMisalStoreException(URV addr, unsigned accessSize,
                                        SecondaryCause& secCause) const
{
  if (not misalDataOk_)
    {
      secCause = SecondaryCause::NONE;
      return ExceptionCause::STORE_ADDR_MISAL;
    }

  size_t addr2 = addr + accessSize - 1;

  // Misaligned access to PIC.
  if (isAddrMemMapped(addr))
    {
      secCause = SecondaryCause::STORE_ACC_PIC;
      return ExceptionCause::STORE_ACC_FAULT;
    }

  // Misaligned access to a region with side effect.
  if (not isIdempotentRegion(addr) or not isIdempotentRegion(addr2))
    {
      secCause = SecondaryCause::STORE_MISAL_IO;
      return ExceptionCause::STORE_ADDR_MISAL;
    }

  // Crossing 256 MB region boundary.
  if (memory_.getRegionIndex(addr) != memory_.getRegionIndex(addr2))
    {
      secCause = SecondaryCause::STORE_MISAL_REGION_CROSS;
      return ExceptionCause::STORE_ADDR_MISAL;
    }

  secCause = SecondaryCause::NONE;
  return ExceptionCause::NONE;
}


template <typename URV>
void
Hart<URV>::initiateLoadException(ExceptionCause cause, URV addr,
				 SecondaryCause secCause)
{
  initiateException(cause, currPc_, addr, secCause);
}


template <typename URV>
void
Hart<URV>::initiateStoreException(ExceptionCause cause, URV addr,
				  SecondaryCause secCause)
{
  initiateException(cause, currPc_, addr, secCause);
}


template <typename URV>
bool
Hart<URV>::effectiveAndBaseAddrMismatch(URV base, URV addr)
{
  unsigned baseRegion = unsigned(base >> (sizeof(URV)*8 - 4));
  unsigned addrRegion = unsigned(addr >> (sizeof(URV)*8 - 4));
  if (baseRegion == addrRegion)
    return false;

  bool flag1 = regionHasLocalDataMem_.at(baseRegion);
  bool flag2 = regionHasLocalDataMem_.at(addrRegion);
  return flag1 != flag2;
}


template <typename URV>
bool
Hart<URV>::checkStackLoad(URV base, URV addr, unsigned loadSize)
{
  URV low = addr;
  URV high = addr + loadSize - 1;
  URV spVal = intRegs_.read(RegSp);
  bool ok = high <= stackMax_ and low >= spVal;
  ok = ok and base <= stackMax_ and base >= stackMin_;
  return ok;
}


template <typename URV>
bool
Hart<URV>::checkStackStore(URV base, URV addr, unsigned storeSize)
{
  URV low = addr;
  URV high = addr + storeSize - 1;
  bool ok = high <= stackMax_ and low >= stackMin_;
  ok = ok and base <= stackMax_ and base >= stackMin_;
  return ok;
}


template <typename URV>
bool
Hart<URV>::wideLoad(uint32_t rd, URV addr, unsigned ldSize)
{
  auto secCause = SecondaryCause::LOAD_ACC_64BIT;
  auto cause = ExceptionCause::LOAD_ACC_FAULT;

  if ((addr & 7) or ldSize != 4 or not isDataAddressExternal(addr))
    {
      initiateLoadException(cause, addr, secCause);
      return false;
    }

  uint32_t upper = 0, lower = 0;
  if (not memory_.read(addr + 4, upper) or not memory_.read(addr, lower))
    {
      initiateLoadException(cause, addr, secCause);
      return false;
    }

  if (loadQueueEnabled_)
    {
      uint32_t prevLower = peekIntReg(rd);
      URV temp = 0;
      peekCsr(CsrNumber::MDBHD, temp);
      uint64_t prevUpper = temp;
      uint64_t prevWide = (prevUpper << 32) | prevLower;
      putInLoadQueue(8, addr, rd, prevWide, true /*isWide*/);
    }

  intRegs_.write(rd, lower);

  pokeCsr(CsrNumber::MDBHD, upper);
  recordCsrWrite(CsrNumber::MDBHD);  // Bench does not handle this.

  return true;
}


template <typename URV>
ExceptionCause
Hart<URV>::determineLoadException(unsigned rs1, URV base, uint64_t& addr,
				  unsigned ldSize, SecondaryCause& secCause)
{
  secCause = SecondaryCause::NONE;

  // Misaligned load from io section triggers an exception. Crossing
  // dccm to non-dccm causes an exception.
  unsigned alignMask = ldSize - 1;
  bool misal = addr & alignMask;
  misalignedLdSt_ = misal;

  ExceptionCause cause = ExceptionCause::NONE;
  if (misal)
    {
      cause = determineMisalLoadException(addr, ldSize, secCause);
      if (cause == ExceptionCause::LOAD_ADDR_MISAL)
        return cause;  // Misaligned resulting in misaligned-adddress-exception
    }

  // Stack access.
  if (rs1 == RegSp and checkStackAccess_ and
      not checkStackLoad(base, addr, ldSize))
    {
      secCause = SecondaryCause::LOAD_ACC_STACK_CHECK;
      return ExceptionCause::LOAD_ACC_FAULT;
    }

  // Address translation
  if (isRvs())
    {
      PrivilegeMode mode = mstatusMprv_? mstatusMpp_ : privMode_;
      if (mode != PrivilegeMode::Machine)
        {
          uint64_t pa = 0;
          cause = virtMem_.translateForLoad(addr, mode, pa);
          if (cause != ExceptionCause::NONE)
            return cause;
          addr = pa;
        }
    }
  else
    {
      // DCCM unmapped
      if (misal)
        {
          size_t lba = addr + ldSize - 1;  // Last byte address
          if (isAddrInDccm(addr) != isAddrInDccm(lba) or
              isAddrMemMapped(addr) != isAddrMemMapped(lba))
            {
              secCause = SecondaryCause::LOAD_ACC_LOCAL_UNMAPPED;
              return ExceptionCause::LOAD_ACC_FAULT;
            }
        }

      // DCCM unmapped or out of MPU range
      bool isReadable = isAddrReadable(addr);
      if (not isReadable)
        {
          secCause = SecondaryCause::LOAD_ACC_MEM_PROTECTION;
          size_t region = memory_.getRegionIndex(addr);
          if (regionHasLocalDataMem_.at(region))
            {
              if (not isAddrMemMapped(addr))
                {
                  secCause = SecondaryCause::LOAD_ACC_LOCAL_UNMAPPED;
                  return ExceptionCause::LOAD_ACC_FAULT;
                }
            }
          else
            return ExceptionCause::LOAD_ACC_FAULT;
        }

      // 64-bit load
      if (wideLdSt_)
        {
          bool fail = (addr & 7) or ldSize != 4 or ! isDataAddressExternal(addr);
          fail = fail or ! isAddrReadable(addr+4);
          if (fail)
            {
              secCause = SecondaryCause::LOAD_ACC_64BIT;
              return ExceptionCause::LOAD_ACC_FAULT;
            }
        }

      // Region predict (Effective address compatible with base).
      if (eaCompatWithBase_ and effectiveAndBaseAddrMismatch(addr, base))
        {
          secCause = SecondaryCause::LOAD_ACC_REGION_PREDICTION;
          return ExceptionCause::LOAD_ACC_FAULT;
        }
    }

  // PIC access
  if (isAddrMemMapped(addr))
    {
      if (privMode_ != PrivilegeMode::Machine)
        {
          secCause = SecondaryCause::LOAD_ACC_LOCAL_UNMAPPED;
	  return ExceptionCause::LOAD_ACC_FAULT;
        }
      if (misal or ldSize != 4)
	{
	  secCause = SecondaryCause::LOAD_ACC_PIC;
	  return ExceptionCause::LOAD_ACC_FAULT;
	}
    }

  // Misaligned resulting in access fault exception.
  if (misal and cause != ExceptionCause::NONE)
    return cause;

  // Physical memory protection.
  if (pmpEnabled_)
    {
      Pmp pmp = pmpManager_.accessPmp(addr);
      if (not pmp.isRead(privMode_, mstatusMpp_, mstatusMprv_) and
          not isAddrMemMapped(addr))
        {
          secCause = SecondaryCause::LOAD_ACC_PMP;
          return ExceptionCause::LOAD_ACC_FAULT;
        }
    }

  // Fault dictated by test-bench.
  if (forceAccessFail_)
    {
      secCause = forcedCause_;
      return ExceptionCause::LOAD_ACC_FAULT;
    }

  return ExceptionCause::NONE;
}


template <typename URV>
template <typename LOAD_TYPE>
inline
bool
Hart<URV>::fastLoad(uint32_t rd, uint32_t rs1, int32_t imm)
{
  URV base = intRegs_.read(rs1);
  URV addr = base + SRV(imm);

  // Unsigned version of LOAD_TYPE
  typedef typename std::make_unsigned<LOAD_TYPE>::type ULT;

  ULT uval = 0;
  if (memory_.read(addr, uval))
    {
      URV value;
      if constexpr (std::is_same<ULT, LOAD_TYPE>::value)
        value = uval;
      else
        value = SRV(LOAD_TYPE(uval)); // Sign extend.

      intRegs_.write(rd, value);
      return true;  // Success.
    }
  return false;
}


template <typename URV>
template <typename LOAD_TYPE>
inline
bool
Hart<URV>::load(uint32_t rd, uint32_t rs1, int32_t imm)
{
#ifdef FAST_SLOPPY
  return fastLoad<LOAD_TYPE>(rd, rs1, imm);
#else

  URV base = intRegs_.read(rs1);
  uint64_t virtAddr = base + SRV(imm);
  unsigned ldSize = sizeof(LOAD_TYPE);

  ldStAddr_ = virtAddr;   // For reporting ld/st addr in trace-mode.
  ldStAddrValid_ = true;  // For reporting ld/st addr in trace-mode.

  if (loadQueueEnabled_)
    removeFromLoadQueue(rs1, false);

  if (hasActiveTrigger())
    {
      if (ldStAddrTriggerHit(virtAddr, TriggerTiming::Before, true /*isLoad*/,
                             privMode_, isInterruptEnabled()))
	triggerTripped_ = true;
    }

  // Unsigned version of LOAD_TYPE
  typedef typename std::make_unsigned<LOAD_TYPE>::type ULT;

  auto secCause = SecondaryCause::NONE;
  uint64_t addr = virtAddr;
  auto cause = determineLoadException(rs1, base, addr, ldSize, secCause);
  if (cause != ExceptionCause::NONE)
    {
      if (triggerTripped_)
        return false;
      initiateLoadException(cause, virtAddr, secCause);
      return false;
    }

  if (wideLdSt_ and not triggerTripped_)
    return wideLoad(rd, addr, ldSize);

  // Loading from console-io does a standard input read.
  if (conIoValid_ and addr == conIo_ and enableConIn_ and not triggerTripped_)
    {
      SRV val = fgetc(stdin);
      intRegs_.write(rd, val);
      return true;
    }

  ULT uval = 0;
  if (memory_.read(addr, uval))
    {
      URV value;
      if constexpr (std::is_same<ULT, LOAD_TYPE>::value)
        value = uval;  // Loading an unsinged.
      else
        value = SRV(LOAD_TYPE(uval)); // Loading signed: Sign extend.

      // Check for load-data-trigger. Load-data-trigger does not apply
      // to io/region unless address is in local memory. Don't ask.
      if (hasActiveTrigger() and
          (isIdempotentRegion(addr) or
           isAddrMemMapped(addr) or isAddrInDccm(addr)))
        {
          TriggerTiming timing = TriggerTiming::Before;
          bool isLoad = true;
          if (ldStDataTriggerHit(uval, timing, isLoad, privMode_,
                                 isInterruptEnabled()))
            triggerTripped_ = true;
        }

      if (addr >= clintStart_ and addr <= clintLimit_ and addr == 0x200bff8)
        {
          value = instCounter_;
        }

      if (not triggerTripped_)
        {
          // Put entry in load queue with value of rd before this load.
          if (loadQueueEnabled_)
            {
              URV prevRdVal = peekIntReg(rd);
              putInLoadQueue(ldSize, addr, rd, prevRdVal);
            }
          intRegs_.write(rd, value);
          return true;  // Success.
        }
    }

  if (triggerTripped_)
    return false;

  cause = ExceptionCause::LOAD_ACC_FAULT;
  secCause = SecondaryCause::LOAD_ACC_MEM_PROTECTION;
  if (isAddrMemMapped(addr))
    secCause = SecondaryCause::LOAD_ACC_PIC;
  initiateLoadException(cause, virtAddr, secCause);
  return false;
#endif
}


template <typename URV>
inline
void
Hart<URV>::execLw(const DecodedInst* di)
{
  load<int32_t>(di->op0(), di->op1(), di->op2As<int32_t>());
}


template <typename URV>
inline
void
Hart<URV>::execLh(const DecodedInst* di)
{
  load<int16_t>(di->op0(), di->op1(), di->op2As<int32_t>());
}


template <typename URV>
template <typename STORE_TYPE>
inline
bool
Hart<URV>::fastStore(unsigned /*rs1*/, URV /*base*/, URV addr,
                     STORE_TYPE storeVal)
{
  if (memory_.write(hartIx_, addr, storeVal))
    {
      if (toHostValid_ and addr == toHost_ and storeVal != 0)
	{
	  throw CoreException(CoreException::Stop, "write to to-host",
			      toHost_, storeVal);
	}
      return true;
    }

  auto secCause = SecondaryCause::NONE;
  initiateStoreException(ExceptionCause::STORE_ACC_FAULT, addr, secCause);
  return false;
}


template <typename URV>
template <typename STORE_TYPE>
inline
bool
Hart<URV>::store(unsigned rs1, URV base, URV virtAddr, STORE_TYPE storeVal)
{
#ifdef FAST_SLOPPY
  return fastStore(rs1, base, virtAddr, storeVal);
#else

  std::lock_guard<std::mutex> lock(memory_.lrMutex_);

  ldStAddr_ = virtAddr;   // For reporting ld/st addr in trace-mode.
  ldStAddrValid_ = true;  // For reporting ld/st addr in trace-mode.

  // ld/st-address or instruction-address triggers have priority over
  // ld/st access or misaligned exceptions.
  bool hasTrig = hasActiveTrigger();
  TriggerTiming timing = TriggerTiming::Before;
  bool isLd = false;  // Not a load.
  if (hasTrig and ldStAddrTriggerHit(virtAddr, timing, isLd, privMode_,
                                     isInterruptEnabled()))
    triggerTripped_ = true;

  // Determine if a store exception is possible.
  STORE_TYPE maskedVal = storeVal;  // Masked store value.
  auto secCause = SecondaryCause::NONE;
  uint64_t addr = virtAddr;
  ExceptionCause cause = determineStoreException(rs1, base, addr,
						 maskedVal, secCause);

  // Consider store-data  trigger
  if (hasTrig and cause == ExceptionCause::NONE)
    if (ldStDataTriggerHit(maskedVal, timing, isLd, privMode_,
                           isInterruptEnabled()))
      triggerTripped_ = true;
  if (triggerTripped_)
    return false;

  if (cause != ExceptionCause::NONE)
    {
      // For the bench: A precise error does write externl memory.
      if (forceAccessFail_ and memory_.isDataAddressExternal(addr))
        memory_.write(hartIx_, addr, storeVal);
      initiateStoreException(cause, virtAddr, secCause);
      return false;
    }

  unsigned stSize = sizeof(STORE_TYPE);
  if (wideLdSt_)
    return wideStore(addr, storeVal, stSize);

  if (memory_.write(hartIx_, addr, storeVal))
    {
      memory_.invalidateOtherHartLr(hartIx_, addr, stSize);

      invalidateDecodeCache(addr, stSize);

      // If we write to special location, end the simulation.
      if (toHostValid_ and addr == toHost_ and storeVal != 0)
	{
	  throw CoreException(CoreException::Stop, "write to to-host",
			      toHost_, storeVal);
	}

      // If addr is special location, then write to console.
      if (conIoValid_ and addr == conIo_)
        {
          if (consoleOut_)
            {
              fputc(storeVal, consoleOut_);
              if (storeVal == '\n')
                fflush(consoleOut_);
            }
          return true;
	}

      if (addr >= clintStart_ and addr <= clintLimit_)
        processClintWrite(addr, stSize, storeVal);

      return true;
    }

  // Store failed: Take exception. Should not happen but we are paranoid.
  initiateStoreException(ExceptionCause::STORE_ACC_FAULT, virtAddr, secCause);
  return false;
#endif
}


template <typename URV>
void
Hart<URV>::processClintWrite(size_t addr, unsigned stSize, URV storeVal)
{
  if (clintTimerAddrToHart_)
    {
      auto hart = clintTimerAddrToHart_(addr);
      if (hart)
        {
          hart->alarmLimit_ = storeVal;
          URV mipVal = hart->csRegs_.peekMip();
          mipVal = mipVal & ~(URV(1) << URV(InterruptCause::M_TIMER));
          hart->pokeCsr(CsrNumber::MIP, mipVal);
          return;
        }
    }

  if (addr == 0x200bff8)
    {
      std::cerr << "Aie mtime updated\n";
      return;
    }

  if (clintSoftAddrToHart_)
    {
      auto hart = clintSoftAddrToHart_(addr);
      if (not hart)
        return;  // Address is not in the software-interrupt memory mapped locations.

      if (stSize != 4)
        return;  // Must be sw

      if ((storeVal >> 1) != 0)
        return;  // Must write 0 or 1.

      URV mipVal = csRegs_.peekMip();
      if (storeVal)
        mipVal = mipVal | (URV(1) << URV(InterruptCause::M_SOFTWARE));
      else
        mipVal = mipVal & ~(URV(1) << URV(InterruptCause::M_SOFTWARE));
      hart->pokeCsr(CsrNumber::MIP, mipVal);
      recordCsrWrite(CsrNumber::MIP);
    }
}


template <typename URV>
inline
void
Hart<URV>::execSw(const DecodedInst* di)
{
  uint32_t rs1 = di->op1();
  URV base = intRegs_.read(rs1);
  URV addr = base + di->op2As<SRV>();
  uint32_t value = uint32_t(intRegs_.read(di->op0()));

  store<uint32_t>(rs1, base, addr, value);
}


template <typename URV>
bool
Hart<URV>::readInst(size_t address, uint32_t& inst)
{
  inst = 0;

  uint16_t low;  // Low 2 bytes of instruction.
  if (not memory_.readInst(address, low))
    return false;

  inst = low;

  if ((inst & 0x3) == 3)  // Non-compressed instruction.
    {
      uint16_t high;
      if (not memory_.readInst(address + 2, high))
	return false;
      inst |= (uint32_t(high) << 16);
    }

  return true;
}


template <typename URV>
bool
Hart<URV>::defineIccm(size_t addr, size_t size)
{
  bool trim = this->findCsr("mpicbaddr") == nullptr;

  bool ok = memory_.defineIccm(addr, size, trim);
  if (ok and trim)
    {
      size_t region = addr/regionSize();
      regionHasLocalMem_.at(region) = true;
      regionHasLocalInstMem_.at(region) = true;
    }
  return ok;
}
    

template <typename URV>
bool
Hart<URV>::defineDccm(size_t addr, size_t size)
{
  bool trim = this->findCsr("mpicbaddr") == nullptr;

  bool ok = memory_.defineDccm(addr, size, trim);
  if (ok and trim)
    {
      size_t region = addr/regionSize();
      regionHasLocalMem_.at(region) = true;
      regionHasLocalDataMem_.at(region) = true;
      regionHasDccm_.at(region) = true;
    }
  return ok;
}


template <typename URV>
bool
Hart<URV>::defineMemoryMappedRegisterArea(size_t addr, size_t size)
{
  // If mpicbaddr CSR is present, the nothing special is done for 256
  // MB region containing memory-mapped-registers. Otherwise, region
  // is marked non accessible except for memory-mapped-register area.
  bool trim = this->findCsr("mpicbaddr") == nullptr;

  bool ok = memory_.defineMemoryMappedRegisterArea(addr, size, trim);
  if (ok and trim)
    {
      size_t region = addr / memory_.regionSize();
      regionHasLocalMem_.at(region) = true;
      regionHasLocalDataMem_.at(region) = true;
      regionHasMemMappedRegs_.at(region) = true;
    }
  return ok;
}


template <typename URV>
bool
Hart<URV>::defineMemoryMappedRegisterWriteMask(size_t addr, uint32_t mask)
{
  return memory_.defineMemoryMappedRegisterWriteMask(addr, mask);
}


template <typename URV>
bool
Hart<URV>::configMemoryFetch(const std::vector< std::pair<URV,URV> >& windows)
{
  using std::cerr;

  size_t regSize = regionSize(), memSize = memorySize();
  if (windows.empty() or memSize == 0 or regSize == 0)
    return true;

  unsigned errors = 0;

  // Mark all pages in non-iccm regions as non executable.
  for (size_t start = 0; start < memSize; start += regSize)
    {
      size_t end = std::min(start + regSize, memSize);
      size_t region = start / regSize;
      if (not regionHasLocalInstMem_.at(region))
        {
          Pma::Attrib attr = Pma::Attrib(Pma::Exec);
          memory_.pmaMgr_.disable(start, end - 1, attr);
        }
    }

  // Mark pages in configuration windows as executable except when
  // they fall in iccm regions.
  for (auto window : windows)
    {
      if (window.first > window.second)
	{
	  cerr << "Invalid memory range in inst fetch configuration: 0x"
	       << std::hex << window.first << " to 0x" << window.second
	       << '\n' << std::dec;
	  errors++;
          continue;
	}

      // Clip window to memory size.
      size_t addr = window.first, end = window.second;
      addr = std::min(memorySize(), addr);

      // Clip window against regions with iccm. Mark what remains as
      // accessible.
      while (addr < end)
        {
          size_t region = addr / regSize;
          if (regionHasLocalInstMem_.at(region))
            {
              addr += regSize;
              continue;
            }

          Pma::Attrib attr = Pma::Attrib(Pma::Exec);
          size_t addr2 = std::min(addr + regSize, end);
          memory_.pmaMgr_.enable(addr, addr2 - 1, attr);
          if (addr2 == end)
            break;
          addr = addr2;
	}
    }

  return errors == 0;
}


template <typename URV>
bool
Hart<URV>::configMemoryDataAccess(const std::vector< std::pair<URV,URV> >& windows)
{
  using std::cerr;

  size_t regSize = regionSize(), memSize = memorySize();
  if (windows.empty() or memSize == 0 or regSize == 0)
    return true;

  unsigned errors = 0;

  // Mark memory in non-dccm/pic regions as non-read non-write.
  for (size_t start = 0; start < memSize; start += regSize)
    {
      size_t end = std::min(start + regSize, memSize);
      size_t region = start / regSize;
      if (not regionHasLocalDataMem_.at(region))
        {
          Pma::Attrib attr = Pma::Attrib(Pma::Read | Pma::Write);
          memory_.pmaMgr_.disable(start, end - 1, attr);
        }
    }
  

  // Mark pages in configuration windows as accessible except when
  // they fall in dccm/pic regions.
  for (auto window : windows)
    {
      if (window.first > window.second)
	{
	  cerr << "Invalid memory range in data access configuration: 0x"
	       << std::hex << window.first << " to 0x" << window.second
	       << '\n' << std::dec;
	  errors++;
          continue;
	}

      // Clip window to memory size.
      size_t addr = window.first, end = window.second;
      addr = std::min(memorySize(), addr);

      // Clip window against regions with dccm/pic. Mark what remains
      // as accessible.
      while (addr < end)
        {
          size_t region = addr / regSize;
          if (regionHasLocalDataMem_.at(region))
            {
              addr += regSize;
              continue;
            }

          Pma::Attrib attr = Pma::Attrib(Pma::Read | Pma::Write);
          size_t addr2 = std::min(addr + regSize, end);
          memory_.pmaMgr_.enable(addr, addr2 - 1, attr);
          if (addr2 == end)
            break;
          addr = addr2;
	}
    }

  return errors == 0;
}


template <typename URV>
inline
bool
Hart<URV>::fetchInst(URV virtAddr, uint32_t& inst)
{
  uint64_t addr = virtAddr;

  // Inst address translation and memory protection is not affected by MPRV.
  bool instMprv = false;

  if (isRvs() and privMode_ != PrivilegeMode::Machine)
    {
      auto cause = virtMem_.translateForFetch(virtAddr, privMode_, addr);
      if (cause != ExceptionCause::NONE)
        {
          initiateException(cause, virtAddr, virtAddr);
          return false;
        }
    }

  if (virtAddr & 1)
    {
      initiateException(ExceptionCause::INST_ADDR_MISAL, virtAddr, virtAddr);
      return false;
    }

  if (forceFetchFail_)
    {
      forceFetchFail_ = false;
      readInst(addr, inst);
      URV info = pc_ + forceFetchFailOffset_;
      auto cause = ExceptionCause::INST_ACC_FAULT;
      auto secCause = SecondaryCause::INST_PRECISE;
      if (isAddrInIccm(addr))
	secCause = SecondaryCause::INST_DOUBLE_ECC;
      initiateException(cause, pc_, info, secCause);
      return false;
    }

  if ((addr & 3) == 0)   // Word aligned
    {
      if (not memory_.readInst(addr, inst))
        {
          auto secCause = SecondaryCause::INST_MEM_PROTECTION;
          size_t region = memory_.getRegionIndex(addr);
          if (regionHasLocalInstMem_.at(region))
            secCause = SecondaryCause::INST_LOCAL_UNMAPPED;
          initiateException(ExceptionCause::INST_ACC_FAULT, virtAddr, virtAddr,
                            secCause);
          return false;
        }

      if (pmpEnabled_)
        {
          Pmp pmp = pmpManager_.accessPmp(addr);
          if (not pmp.isExec(privMode_, mstatusMpp_, instMprv))
            {
              auto secCause = SecondaryCause::INST_PMP;
              initiateException(ExceptionCause::INST_ACC_FAULT, virtAddr, virtAddr,
                                secCause);
              return false;
            }
        }

      return true;
    }

  uint16_t half;
  if (not memory_.readInst(addr, half))
    {
      auto secCause = SecondaryCause::INST_MEM_PROTECTION;
      size_t region = memory_.getRegionIndex(addr);
      if (regionHasLocalInstMem_.at(region))
	secCause = SecondaryCause::INST_LOCAL_UNMAPPED;
      initiateException(ExceptionCause::INST_ACC_FAULT, virtAddr, virtAddr, secCause);
      return false;
    }

  if (pmpEnabled_)
    {
      Pmp pmp = pmpManager_.accessPmp(addr);
      if (not pmp.isExec(privMode_, mstatusMpp_, instMprv))
        {
          auto secCause = SecondaryCause::INST_PMP;
          initiateException(ExceptionCause::INST_ACC_FAULT, virtAddr, virtAddr,
                            secCause);
          return false;
        }
    }

  inst = half;
  if (isCompressedInst(inst))
    return true;

  if (isRvs() and privMode_ != PrivilegeMode::Machine)
    {
      auto cause = virtMem_.translateForFetch(virtAddr+2, privMode_, addr);
      if (cause != ExceptionCause::NONE)
        {
          initiateException(cause, virtAddr, virtAddr+2);
          return false;
        }
    }
  else
    addr += 2;

  uint16_t upperHalf;
  if (not memory_.readInst(addr, upperHalf))
    {
      // 4-byte instruction: 4-byte fetch failed but 1st 2-byte fetch
      // succeeded. Problem must be in 2nd half of instruction.
      auto secCause = SecondaryCause::INST_MEM_PROTECTION;
      size_t region = memory_.getRegionIndex(addr);
      if (regionHasLocalInstMem_.at(region))
        secCause = SecondaryCause::INST_LOCAL_UNMAPPED;
      initiateException(ExceptionCause::INST_ACC_FAULT, virtAddr, virtAddr + 2,
                        secCause);

      return false;
    }

  if (pmpEnabled_)
    {
      Pmp pmp = pmpManager_.accessPmp(addr);
      if (not pmp.isExec(privMode_, mstatusMpp_, instMprv))
        {
          auto secCause = SecondaryCause::INST_PMP;
          initiateException(ExceptionCause::INST_ACC_FAULT, virtAddr, virtAddr + 2,
                            secCause);
          return false;
        }
    }

  inst = inst | (uint32_t(upperHalf) << 16);
  return true;
}


template <typename URV>
bool
Hart<URV>::fetchInstPostTrigger(URV virtAddr, uint32_t& inst, FILE* traceFile)
{
  uint64_t addr = virtAddr;
  bool pageFault = false;
  if (isRvs() and privMode_ != PrivilegeMode::Machine)
    {
      auto cause = virtMem_.translateForFetch(virtAddr, privMode_, addr);
      pageFault = cause != ExceptionCause::NONE;
    }

  // Fetch will fail if page fault or forced or if address is
  // misaligned or if memory read fails.
  if (not pageFault and not forceFetchFail_ and (addr & 1) == 0)
    {
      if (memory_.readInst(addr, inst))
	return true;  // Read 4 bytes: success.

      uint16_t half;
      if (memory_.readInst(addr, half))
	{
	  if (isCompressedInst(inst))
	    return true; // Read 2 bytes and compressed inst: success.
	}
    }

  // Fetch failed: take pending trigger-exception.
  URV info = virtAddr;
  takeTriggerAction(traceFile, addr, info, instCounter_, true);
  forceFetchFail_ = false;

  return false;
}


template <typename URV>
void
Hart<URV>::illegalInst(const DecodedInst* di)
{
  if (triggerTripped_)
    return;

  // Check if stuck because of lack of illegal instruction exception handler.
  if (instCounter_ == counterAtLastIllegal_ + 1)
    consecutiveIllegalCount_++;
  else
    consecutiveIllegalCount_ = 0;

  if (consecutiveIllegalCount_ > 64)  // FIX: Make a parameter
    {
      throw CoreException(CoreException::Stop,
			  "64 consecutive illegal instructions",
			  0, 0);
    }

  counterAtLastIllegal_ = instCounter_;

  uint32_t inst = di->inst();
  if (isCompressedInst(inst))
    inst = inst & 0xffff;

  initiateException(ExceptionCause::ILLEGAL_INST, currPc_, inst);
}


template <typename URV>
void
Hart<URV>::unimplemented(const DecodedInst* di)
{
  illegalInst(di);
}


// This is a swerv-specific special code that corresponds to special
// hardware that maps the interrupt id (claim id) to a specific
// interrupt handler routine by looking up the routine address in a
// table.
template <typename URV>
void
Hart<URV>::initiateFastInterrupt(InterruptCause cause, URV pcToSave)
{
  // Get the address of the interrupt handler entry from meihap
  // register.
  URV addr = 0;
  if (not csRegs_.read(CsrNumber::MEIHAP, PrivilegeMode::Machine, addr))
    {
      initiateNmi(URV(NmiCause::UNKNOWN), pcToSave);
      return;
    }

  // Check that the entry address is in a DCCM region.
  size_t ix = memory_.getRegionIndex(addr);
  if (not regionHasLocalDataMem_.at(ix))
    {
      initiateNmi(URV(NmiCause::NON_DCCM_ACCESS_ERROR), pcToSave);
      return;
    }

  // If bench has forced an ECC error, honor it.
  if (forceAccessFail_)
    {
      initiateNmi(URV(NmiCause::DOUBLE_BIT_ECC), pcToSave);
      forceAccessFail_ = false;
      return;
    }

  // Fetch the interrupt handler address.
  URV nextPc = 0;
  if (not memory_.read(addr, nextPc))
    {
      initiateNmi(URV(NmiCause::DCCM_ACCESS_ERROR), pcToSave);
      return;
    }

  URV causeVal = URV(cause);
  causeVal |= URV(1) << (mxlen_ - 1);  // Set most sig bit.
  undelegatedInterrupt(causeVal, pcToSave, nextPc);

  if (instFreq_)
    accumulateTrapStats(false /* isNmi*/);

  if (not enableCounters_)
    return;

  PerfRegs& pregs = csRegs_.mPerfRegs_;
  if (cause == InterruptCause::M_EXTERNAL)
    pregs.updateCounters(EventNumber::ExternalInterrupt, prevPerfControl_,
                         lastPriv_);
  else if (cause == InterruptCause::M_TIMER)
    pregs.updateCounters(EventNumber::TimerInterrupt, prevPerfControl_,
                         lastPriv_);
}


// Start an asynchronous exception.
template <typename URV>
void
Hart<URV>::initiateInterrupt(InterruptCause cause, URV pc)
{
  if (fastInterrupts_ and cause == InterruptCause::M_EXTERNAL)
    {
      initiateFastInterrupt(cause, pc);
      return;
    }

  bool interrupt = true;
  URV info = 0;  // This goes into mtval.
  auto secCause = SecondaryCause::NONE;
  initiateTrap(interrupt, URV(cause), pc, info, URV(secCause));

  interruptCount_++;

  if (not enableCounters_)
    return;

  PerfRegs& pregs = csRegs_.mPerfRegs_;
  if (cause == InterruptCause::M_EXTERNAL)
    pregs.updateCounters(EventNumber::ExternalInterrupt, prevPerfControl_,
                         lastPriv_);
  else if (cause == InterruptCause::M_TIMER)
    pregs.updateCounters(EventNumber::TimerInterrupt, prevPerfControl_,
                         lastPriv_);
}


// Start a synchronous exception.
template <typename URV>
void
Hart<URV>::initiateException(ExceptionCause cause, URV pc, URV info,
			     SecondaryCause secCause)
{
  bool interrupt = false;
  exceptionCount_++;
  hasException_ = true;
  initiateTrap(interrupt, URV(cause), pc, info, URV(secCause));

  PerfRegs& pregs = csRegs_.mPerfRegs_;
  if (enableCounters_)
    pregs.updateCounters(EventNumber::Exception, prevPerfControl_,
                         lastPriv_);
}


template <typename URV>
void
Hart<URV>::initiateTrap(bool interrupt, URV cause, URV pcToSave, URV info,
			URV secCause)
{
  forceAccessFail_ = false;

  enableWideLdStMode(false);  // Swerv specific feature.

  memory_.invalidateLr(hartIx_);

  PrivilegeMode origMode = privMode_;

  // Exceptions are taken in machine mode.
  privMode_ = PrivilegeMode::Machine;
  PrivilegeMode nextMode = PrivilegeMode::Machine;

  // But they can be delegated.
  if (isRvs() and origMode != PrivilegeMode::Machine)
    {
      CsrNumber csrn = interrupt? CsrNumber::MIDELEG : CsrNumber::MEDELEG;
      URV delegVal = 0;
      peekCsr(csrn, delegVal);
      if (delegVal & (URV(1) << cause))
        nextMode = PrivilegeMode::Supervisor;
    }

  CsrNumber epcNum = CsrNumber::MEPC;
  CsrNumber causeNum = CsrNumber::MCAUSE;
  CsrNumber scauseNum = CsrNumber::MSCAUSE;
  CsrNumber tvalNum = CsrNumber::MTVAL;
  CsrNumber tvecNum = CsrNumber::MTVEC;

  if (nextMode == PrivilegeMode::Supervisor)
    {
      epcNum = CsrNumber::SEPC;
      causeNum = CsrNumber::SCAUSE;
      tvalNum = CsrNumber::STVAL;
      tvecNum = CsrNumber::STVEC;
    }
  else if (nextMode == PrivilegeMode::User)
    {
      epcNum = CsrNumber::UEPC;
      causeNum = CsrNumber::UCAUSE;
      tvalNum = CsrNumber::UTVAL;
      tvecNum = CsrNumber::UTVEC;
    }

  // Save address of instruction that caused the exception or address
  // of interrupted instruction.
  if (not csRegs_.write(epcNum, privMode_, pcToSave & ~(URV(1))))
    assert(0 and "Failed to write EPC register");

  // Save the exception cause.
  URV causeRegVal = cause;
  if (interrupt)
    causeRegVal |= URV(1) << (mxlen_ - 1);
  if (not csRegs_.write(causeNum, privMode_, causeRegVal))
    assert(0 and "Failed to write CAUSE register");

  // Save secondary exception cause (WD special).
  csRegs_.write(scauseNum, privMode_, secCause);

  // Clear mtval on interrupts. Save synchronous exception info.
  if (not csRegs_.write(tvalNum, privMode_, info))
    assert(0 and "Failed to write TVAL register");

  // Update status register saving xIE in xPIE and previous privilege
  // mode in xPP by getting current value of mstatus ...
  URV status = csRegs_.peekMstatus();

  // ... updating its fields
  MstatusFields<URV> msf(status);

  if (nextMode == PrivilegeMode::Machine)
    {
      msf.bits_.MPP = unsigned(origMode);
      msf.bits_.MPIE = msf.bits_.MIE;
      msf.bits_.MIE = 0;
    }
  else if (nextMode == PrivilegeMode::Supervisor)
    {
      msf.bits_.SPP = unsigned(origMode);
      msf.bits_.SPIE = msf.bits_.SIE;
      msf.bits_.SIE = 0;
    }
  else if (nextMode == PrivilegeMode::User)
    {
      msf.bits_.UPIE = msf.bits_.UIE;
      msf.bits_.UIE = 0;
    }

  // ... and putting it back
  if (not csRegs_.write(CsrNumber::MSTATUS, privMode_, msf.value_))
    assert(0 and "Failed to write MSTATUS register");
  updateCachedMstatusFields();
  
  // Set program counter to trap handler address.
  URV tvec = 0;
  if (not csRegs_.read(tvecNum, privMode_, tvec))
    assert(0 and "Failed to read TVEC register");

  URV base = (tvec >> 2) << 2;  // Clear least sig 2 bits.
  unsigned tvecMode = tvec & 0x3;

  if (tvecMode == 1 and interrupt)
    base = base + 4*cause;

  pc_ = base;

  // Change privilege mode.
  privMode_ = nextMode;

  if (instFreq_)
    accumulateTrapStats(false /*isNmi*/);
}


template <typename URV>
void
Hart<URV>::initiateNmi(URV cause, URV pcToSave)
{
  URV nextPc = nmiPc_;
  undelegatedInterrupt(cause, pcToSave, nextPc);
  nmiCount_++;
  if (instFreq_)
    accumulateTrapStats(true);
}


template <typename URV>
void
Hart<URV>::undelegatedInterrupt(URV cause, URV pcToSave, URV nextPc)
{
  enableWideLdStMode(false);  // Swerv specific feature.

  interruptCount_++;
  memory_.invalidateLr(hartIx_);

  PrivilegeMode origMode = privMode_;

  // NMI is taken in machine mode.
  privMode_ = PrivilegeMode::Machine;

  // Save address of instruction that caused the exception or address
  // of interrupted instruction.
  pcToSave = (pcToSave >> 1) << 1; // Clear least sig bit.
  if (not csRegs_.write(CsrNumber::MEPC, privMode_, pcToSave))
    assert(0 and "Failed to write EPC register");

  // Save the exception cause.
  if (not csRegs_.write(CsrNumber::MCAUSE, privMode_, cause))
    assert(0 and "Failed to write CAUSE register");

  // Save secondary exception cause (WD special).
  csRegs_.write(CsrNumber::MSCAUSE, privMode_, 0);

  // Clear mtval
  if (not csRegs_.write(CsrNumber::MTVAL, privMode_, 0))
    assert(0 and "Failed to write MTVAL register");

  // Update status register saving xIE in xPIE and previous privilege
  // mode in xPP by getting current value of mstatus ...
  URV status = csRegs_.peekMstatus();
  MstatusFields<URV> msf(status);

  // ... updating its fields
  msf.bits_.MPP = unsigned(origMode);
  msf.bits_.MPIE = msf.bits_.MIE;
  msf.bits_.MIE = 0;

  // ... and putting it back
  if (not csRegs_.write(CsrNumber::MSTATUS, privMode_, msf.value_))
    assert(0 and "Failed to write MSTATUS register");
  updateCachedMstatusFields();
  
  // Clear pending nmi bit in dcsr
  URV dcsrVal = 0;
  if (peekCsr(CsrNumber::DCSR, dcsrVal))
    {
      dcsrVal &= ~(URV(1) << 3);
      pokeCsr(CsrNumber::DCSR, dcsrVal);
      recordCsrWrite(CsrNumber::DCSR);
    }

  pc_ = (nextPc >> 1) << 1;  // Clear least sig bit
}


template <typename URV>
bool
Hart<URV>::peekIntReg(unsigned ix, URV& val) const
{ 
  if (ix < intRegs_.size())
    {
      val = intRegs_.read(ix);
      return true;
    }
  return false;
}


template <typename URV>
URV
Hart<URV>::peekIntReg(unsigned ix) const
{ 
  assert(ix < intRegs_.size());
  return intRegs_.read(ix);
}


template <typename URV>
bool
Hart<URV>::peekIntReg(unsigned ix, URV& val, std::string& name) const
{ 
  if (ix < intRegs_.size())
    {
      val = intRegs_.read(ix);
      name = intRegName(ix);
      return true;
    }
  return false;
}


template <typename URV>
bool
Hart<URV>::peekFpReg(unsigned ix, uint64_t& val) const
{ 
  if (not isRvf() and not isRvd())
    return false;

  if (ix < fpRegs_.size())
    {
      val = fpRegs_.readBitsRaw(ix);
      return true;
    }

  return false;
}


template <typename URV>
bool
Hart<URV>::peekUnboxedFpReg(unsigned ix, uint64_t& val) const
{ 
  if (not isRvf() and not isRvd())
    return false;

  if (ix < fpRegs_.size())
    {
      val = fpRegs_.readBitsUnboxed(ix);
      return true;
    }

  return false;
}


template <typename URV>
bool
Hart<URV>::pokeFpReg(unsigned ix, uint64_t val)
{ 
  if (not isRvf() and not isRvd())
    return false;

  if (ix < fpRegs_.size())
    {
      fpRegs_.pokeBits(ix, val);
      return true;
    }

  return false;
}


template <typename URV>
bool
Hart<URV>::pokeIntReg(unsigned ix, URV val)
{ 
  if (ix < intRegs_.size())
    {
      intRegs_.poke(ix, val);
      return true;
    }
  return false;
}


template <typename URV>
bool
Hart<URV>::peekCsr(CsrNumber csrn, URV& val) const
{ 
  return csRegs_.peek(csrn, val);
}


template <typename URV>
bool
Hart<URV>::peekCsr(CsrNumber csrn, URV& val, URV& reset, URV& writeMask,
		   URV& pokeMask) const
{ 
  const Csr<URV>* csr = csRegs_.getImplementedCsr(csrn);
  if (not csr)
    return false;

  if (not peekCsr(csrn, val))
    return false;

  reset = csr->getResetValue();
  writeMask = csr->getWriteMask();
  pokeMask = csr->getPokeMask();
  return true;
}


template <typename URV>
bool
Hart<URV>::peekCsr(CsrNumber csrn, URV& val, std::string& name) const
{ 
  const Csr<URV>* csr = csRegs_.getImplementedCsr(csrn);
  if (not csr)
    return false;

  if (not peekCsr(csrn, val))
    return false;

  name = csr->getName();
  return true;
}


template <typename URV>
bool
Hart<URV>::pokeCsr(CsrNumber csr, URV val)
{ 
  // Direct write to MEIHAP will not affect claimid field. Poking
  // MEIHAP will only affect the claimid field.
  if (csr == CsrNumber::MEIHAP)
    {
      URV claimIdMask = 0x3fc;
      URV prev = 0;
      if (not csRegs_.read(CsrNumber::MEIHAP, PrivilegeMode::Machine, prev))
	return false;
      URV newVal = (prev & ~claimIdMask) | (val & claimIdMask);
      csRegs_.poke(CsrNumber::MEIHAP, newVal);
      return true;
    }

  // Some/all bits of some CSRs are read only to CSR instructions but
  // are modifiable. Use the poke method (instead of write) to make
  // sure modifiable value are changed.
  if (not csRegs_.poke(csr, val))
    return false;

  if (csr == CsrNumber::DCSR)
    {
      dcsrStep_ = (val >> 2) & 1;
      dcsrStepIe_ = (val >> 11) & 1;
    }
  else if (csr >= CsrNumber::MSPCBA and csr <= CsrNumber::MSPCC)
    updateStackChecker();
  else if (csr == CsrNumber::MDBAC)
    enableWideLdStMode(true);
  else if ((csr >= CsrNumber::PMPADDR0 and csr <= CsrNumber::PMPADDR15) or
           (csr >= CsrNumber::PMPCFG0 and csr <= CsrNumber::PMPCFG3))
    updateMemoryProtection();
  else if (csr == CsrNumber::SATP)
    updateAddressTranslation();
  else if (csr == CsrNumber::FCSR or csr == CsrNumber::FRM or csr == CsrNumber::FFLAGS)
    markFsDirty();   // Update FS field of MSTATS if FCSR is written

  // Update cached values of MSTATUS MPP and MPRV.
  if (csr == CsrNumber::MSTATUS or csr == CsrNumber::SSTATUS)
    updateCachedMstatusFields();

  return true;
}


template <typename URV>
URV
Hart<URV>::peekPc() const
{
  return pc_;
}


template <typename URV>
void
Hart<URV>::pokePc(URV address)
{
  pc_ = (address >> 1) << 1; // Clear least sig big
}


template <typename URV>
bool
Hart<URV>::findIntReg(const std::string& name, unsigned& num) const
{
  if (intRegs_.findReg(name, num))
    return true;

  unsigned n = 0;
  if (parseNumber<unsigned>(name, n) and n < intRegs_.size())
    {
      num = n;
      return true;
    }

  return false;
}


template <typename URV>
bool
Hart<URV>::findFpReg(const std::string& name, unsigned& num) const
{
  if (not isRvf())
    return false;   // Floating point extension not enabled.

  if (name.empty())
    return false;

  if (name.at(0) == 'f')
    {
      std::string numStr = name.substr(1);
      unsigned n = 0;
      if (parseNumber<unsigned>(numStr, num) and n < fpRegCount())
	return true;
    }

  unsigned n = 0;
  if (parseNumber<unsigned>(name, n) and n < fpRegCount())
    {
      num = n;
      return true;
    }

  return false;
}



template <typename URV>
Csr<URV>*
Hart<URV>::findCsr(const std::string& name)
{
  Csr<URV>* csr = csRegs_.findCsr(name);

  if (not csr)
    {
      unsigned n = 0;
      if (parseNumber<unsigned>(name, n))
	csr = csRegs_.findCsr(CsrNumber(n));
    }

  return csr;
}


template <typename URV>
bool
Hart<URV>::configCsr(const std::string& name, bool implemented, URV resetValue,
                     URV mask, URV pokeMask, bool debug, bool shared)
{
  return csRegs_.configCsr(name, implemented, resetValue, mask, pokeMask,
			   debug, shared);
}


template <typename URV>
bool
Hart<URV>::defineCsr(const std::string& name, CsrNumber num,
		     bool implemented, URV resetVal, URV mask,
		     URV pokeMask, bool isDebug)
{
  bool mandatory = false, quiet = true;
  auto c = csRegs_.defineCsr(name, num, mandatory, implemented, resetVal,
			     mask, pokeMask, isDebug, quiet);
  return c != nullptr;
}


template <typename URV>
bool
Hart<URV>::configMachineModePerfCounters(unsigned numCounters)
{
  return csRegs_.configMachineModePerfCounters(numCounters);
}


template <typename URV>
bool
Hart<URV>::configUserModePerfCounters(unsigned numCounters)
{
  return csRegs_.configUserModePerfCounters(numCounters);
}


template <typename URV>
bool
Hart<URV>::configMemoryProtectionGrain(uint64_t size)
{
  bool ok = true;

  if (size < 4)
    {
      std::cerr << "Memory protection grain size (" << size << ") is "
                << "smaller than 4. Using 4.\n";
      size = 4;
      ok = false;
    }

  uint64_t log2Size = static_cast<uint64_t>(std::log2(size));
  uint64_t powerOf2 = uint64_t(1) << log2Size;
  if (size != powerOf2)
    {
      std::cerr << "Memory protection grain size (0x" << std::hex
                << size << ") is not a power of 2. Using: 0x"
                << powerOf2 << '\n';
      size = powerOf2;
      ok = false;
    }

  if (log2Size > 33)
    {
      std::cerr << "Memory protection grain size (0x" << std::hex
                << size << ") is larger than 2 to the power 33. "
                << "Using 2 to the power 33.\n";
      size = uint64_t(1) << 33;
      powerOf2 = size;
      log2Size = 33;
      ok = false;
    }

  unsigned pmpG = log2Size - 2;
  csRegs_.setPmpG(pmpG);

  return ok;
}


template <typename URV>
void
formatInstTrace(FILE* out, uint64_t tag, unsigned hartId, URV currPc,
		const char* opcode, char resource, URV addr,
		URV value, const char* assembly);

template <>
void
formatInstTrace<uint32_t>(FILE* out, uint64_t tag, unsigned hartId, uint32_t currPc,
		const char* opcode, char resource, uint32_t addr,
		uint32_t value, const char* assembly)
{
  if (resource == 'r')
    {
      fprintf(out, "#%" PRId64 " %d %08x %8s r %02x         %08x  %s",
              tag, hartId, currPc, opcode, addr, value, assembly);
    }
  else if (resource == 'c')
    {
      if ((addr >> 16) == 0)
        fprintf(out, "#%" PRId64 " %d %08x %8s c %04x       %08x  %s",
                tag, hartId, currPc, opcode, addr, value, assembly);
      else
        fprintf(out, "#%" PRId64 " %d %08x %8s c %08x   %08x  %s",
                tag, hartId, currPc, opcode, addr, value, assembly);
    }
  else
    {
      fprintf(out, "#%" PRId64 " %d %08x %8s %c %08x   %08x  %s", tag, hartId,
              currPc, opcode, resource, addr, value, assembly);
    }
}

template <>
void
formatInstTrace<uint64_t>(FILE* out, uint64_t tag, unsigned hartId, uint64_t currPc,
		const char* opcode, char resource, uint64_t addr,
		uint64_t value, const char* assembly)
{
  fprintf(out, "#%" PRId64 " %d %016" PRIx64 " %8s %c %016" PRIx64 " %016" PRIx64 "  %s",
          tag, hartId, currPc, opcode, resource, addr, value, assembly);
}

template <typename URV>
void
formatFpInstTrace(FILE* out, uint64_t tag, unsigned hartId, URV currPc,
		  const char* opcode, unsigned fpReg,
		  uint64_t fpVal, const char* assembly);

template <>
void
formatFpInstTrace<uint32_t>(FILE* out, uint64_t tag, unsigned hartId, uint32_t currPc,
		  const char* opcode, unsigned fpReg,
		  uint64_t fpVal, const char* assembly)
{
  fprintf(out, "#%" PRId64 " %d %08x %8s f %02x %016" PRIx64 "  %s",
          tag, hartId, currPc, opcode, fpReg, fpVal, assembly);
}

template <>
void
formatFpInstTrace<uint64_t>(FILE* out, uint64_t tag, unsigned hartId, uint64_t currPc,
		  const char* opcode, unsigned fpReg,
		  uint64_t fpVal, const char* assembly)
{
  fprintf(out, "#%" PRId64 " %d %016" PRIx64 " %8s f %016" PRIx64 " %016" PRIx64 "  %s",
          tag, hartId, currPc, opcode, uint64_t(fpReg), fpVal, assembly);
}


static std::mutex printInstTraceMutex;
static std::mutex stderrMutex;

template <typename URV>
void
Hart<URV>::printInstTrace(uint32_t inst, uint64_t tag, std::string& tmp,
			  FILE* out, bool interrupt)
{
  DecodedInst di;
  decode(pc_, inst, di);

  printInstTrace(di, tag, tmp, out, interrupt);
}


template <typename URV>
void
Hart<URV>::printInstTrace(const DecodedInst& di, uint64_t tag, std::string& tmp,
			  FILE* out, bool interrupt)
{
  // Serialize to avoid jumbled output.
  std::lock_guard<std::mutex> guard(printInstTraceMutex);

  disassembleInst(di, tmp);
  if (interrupt)
    tmp += " (interrupted)";

  if (traceLdSt_ and ldStAddrValid_)
    {
      std::ostringstream oss;
      oss << "0x" << std::hex << ldStAddr_;
      tmp += " [" + oss.str() + "]";
    }

  char instBuff[128];
  if (di.instSize() == 4)
    sprintf(instBuff, "%08x", di.inst());
  else
    sprintf(instBuff, "%04x", di.inst() & 0xffff);

  bool pending = false;  // True if a printed line need to be terminated.

  // Process integer register diff.
  int reg = intRegs_.getLastWrittenReg();
  URV value = 0;
  if (reg > 0)
    {
      value = intRegs_.read(reg);
      formatInstTrace<URV>(out, tag, hartIx_, currPc_, instBuff, 'r', reg,
			   value, tmp.c_str());
      pending = true;
    }

  // Process floating point register diff.
  int fpReg = fpRegs_.getLastWrittenReg();
  if (fpReg >= 0)
    {
      uint64_t val = fpRegs_.readBitsRaw(fpReg);
      if (pending) fprintf(out, "  +\n");
      formatFpInstTrace<URV>(out, tag, hartIx_, currPc_, instBuff, fpReg,
			     val, tmp.c_str());
      pending = true;
    }

  // Process CSR diffs.
  std::vector<CsrNumber> csrs;
  std::vector<unsigned> triggers;
  csRegs_.getLastWrittenRegs(csrs, triggers);

  std::vector<bool> tdataChanged(3);

  std::map<URV, URV> csrMap; // Map csr-number to its value.

  for (CsrNumber csr : csrs)
    {
      if (not csRegs_.read(csr, PrivilegeMode::Machine, value))
	continue;

      if (csr >= CsrNumber::TDATA1 and csr <= CsrNumber::TDATA3)
	{
	  size_t ix = size_t(csr) - size_t(CsrNumber::TDATA1);
	  tdataChanged.at(ix) = true;
	  continue; // Debug triggers printed separately below
	}
      csrMap[URV(csr)] = value;
    }

  // Process trigger register diffs.
  for (unsigned trigger : triggers)
    {
      URV data1(0), data2(0), data3(0);
      if (not peekTrigger(trigger, data1, data2, data3))
	continue;
      if (tdataChanged.at(0))
	{
	  URV ecsr = (trigger << 16) | URV(CsrNumber::TDATA1);
	  csrMap[ecsr] = data1;
	}
      if (tdataChanged.at(1))
	{
	  URV ecsr = (trigger << 16) | URV(CsrNumber::TDATA2);
	  csrMap[ecsr] = data2;
	}
      if (tdataChanged.at(2))
	{
	  URV ecsr = (trigger << 16) | URV(CsrNumber::TDATA3);
	  csrMap[ecsr] = data3;
	}
    }

  for (const auto& [key, val] : csrMap)
    {
      if (pending) fprintf(out, "  +\n");
      formatInstTrace<URV>(out, tag, hartIx_, currPc_, instBuff, 'c',
			   key, val, tmp.c_str());
      pending = true;
    }

  // Process memory diff.
  size_t address = 0;
  uint64_t memValue = 0;
  unsigned writeSize = memory_.getLastWriteNewValue(hartIx_, address, memValue);
  if (writeSize > 0)
    {
      if (pending)
	fprintf(out, "  +\n");

      formatInstTrace<URV>(out, tag, hartIx_, currPc_, instBuff, 'm',
			   URV(address), URV(memValue), tmp.c_str());
      pending = true;
    }

  if (pending) 
    fprintf(out, "\n");
  else
    {
      // No diffs: Generate an x0 record.
      formatInstTrace<URV>(out, tag, hartIx_, currPc_, instBuff, 'r', 0, 0,
			  tmp.c_str());
      fprintf(out, "\n");
    }
}


template <typename URV>
void
Hart<URV>::recordDivInst(unsigned rd, URV value)
{
  for (size_t i = loadQueue_.size(); i > 0; --i)
    {
      auto& entry = loadQueue_.at(i-1);
      if (entry.isValid() and entry.regIx_ == rd)
        value = entry.prevData_;
    }

  hasLastDiv_ = true;
  priorDivRdVal_ = value;

  lastDivRd_ = rd;
}


template <typename URV>
bool
Hart<URV>::redirectOutputDescriptor(int fd, const std::string& path)
{
  return syscall_.redirectOutputDescriptor(fd, path);
}


template <typename URV>
bool
Hart<URV>::cancelLastDiv()
{
  if (not hasLastDiv_)
    return false;

  hasLastDiv_ = false;
  return pokeIntReg(lastDivRd_, priorDivRdVal_);
}



template <typename URV>
void
Hart<URV>::undoForTrigger()
{
  unsigned regIx = 0;
  URV value = 0;
  if (intRegs_.getLastWrittenReg(regIx, value))
    pokeIntReg(regIx, value);

  intRegs_.clearLastWrittenReg();

  pc_ = currPc_;
}


void
addToSignedHistogram(std::vector<uint64_t>& histo, int64_t val)
{
  if (histo.size() < 13)
    histo.resize(13);

  if (val < 0)
    {
      if      (val <= -64*1024) histo.at(0)++;
      else if (val <= -1024)    histo.at(1)++;
      else if (val <= -16)      histo.at(2)++;
      else if (val < -2)        histo.at(3)++;
      else if (val == -2)       histo.at(4)++;
      else if (val == -1)       histo.at(5)++;
    }
  else
    {
      if      (val == 0)       histo.at(6)++;
      else if (val == 1)       histo.at(7)++;
      else if (val == 2)       histo.at(8)++;
      else if (val <= 16)      histo.at(9)++;
      else if (val <= 1024)    histo.at(10)++;
      else if (val <= 64*1024) histo.at(11)++;
      else                     histo.at(12)++;
    }
}


void
addToUnsignedHistogram(std::vector<uint64_t>& histo, uint64_t val)
{
  if (histo.size() < 13)
    histo.resize(13);

  if      (val == 0)       histo.at(0)++;
  else if (val == 1)       histo.at(1)++;
  else if (val == 2)       histo.at(2)++;
  else if (val <= 16)      histo.at(3)++;
  else if (val <= 1024)    histo.at(4)++;
  else if (val <= 64*1024) histo.at(5)++;
  else                     histo.at(6)++;
}


/// Return true if given hart is in debug mode and the stop count bit of
/// the DSCR register is set.
template <typename URV>
bool
isDebugModeStopCount(const Hart<URV>& hart)
{
  if (not hart.inDebugMode())
    return false;

  URV dcsrVal = 0;
  if (not hart.peekCsr(CsrNumber::DCSR, dcsrVal))
    return false;

  if ((dcsrVal >> 10) & 1)
    return true;  // stop count bit is set
  return false;
}


template <typename URV>
void
Hart<URV>::updatePerformanceCounters(uint32_t inst, const InstEntry& info,
				     uint32_t op0, uint32_t op1)
{
  InstId id = info.instId();

  if (isDebugModeStopCount(*this))
    return;

  // We do not update the performance counters if an instruction
  // causes an exception unless it is an ebreak or an ecall.
  if (hasException_ and id != InstId::ecall and id != InstId::ebreak and
      id != InstId::c_ebreak)
    return;

  PerfRegs& pregs = csRegs_.mPerfRegs_;
  pregs.updateCounters(EventNumber::InstCommited, prevPerfControl_,
                       lastPriv_);

  if (isCompressedInst(inst))
    pregs.updateCounters(EventNumber::Inst16Commited, prevPerfControl_,
                         lastPriv_);
  else
    pregs.updateCounters(EventNumber::Inst32Commited, prevPerfControl_,
                         lastPriv_);

#if 0
  if ((currPc_ & 3) == 0)
    pregs.updateCounters(EventNumber::InstAligned, prevPerfControl_,
                         lastPriv_);
#endif

  if (info.type() == InstType::Int)
    {
      if (id == InstId::ebreak or id == InstId::c_ebreak)
	pregs.updateCounters(EventNumber::Ebreak, prevPerfControl_,
                             lastPriv_);
      else if (id == InstId::ecall)
	pregs.updateCounters(EventNumber::Ecall, prevPerfControl_,
                             lastPriv_);
      else if (id == InstId::fence)
	pregs.updateCounters(EventNumber::Fence, prevPerfControl_,
                             lastPriv_);
      else if (id == InstId::fencei)
	pregs.updateCounters(EventNumber::Fencei, prevPerfControl_,
                             lastPriv_);
      else if (id == InstId::mret)
	pregs.updateCounters(EventNumber::Mret, prevPerfControl_,
                             lastPriv_);
      else if (id != InstId::illegal)
	pregs.updateCounters(EventNumber::Alu, prevPerfControl_,
                             lastPriv_);
    }
  else if (info.isMultiply())
    {
      pregs.updateCounters(EventNumber::Mult, prevPerfControl_,
                           lastPriv_);
    }
  else if (info.isDivide())
    {
      pregs.updateCounters(EventNumber::Div, prevPerfControl_,
                           lastPriv_);
    }
  else if (info.isLoad())
    {
      pregs.updateCounters(EventNumber::Load, prevPerfControl_,
                           lastPriv_);
      if (misalignedLdSt_)
	pregs.updateCounters(EventNumber::MisalignLoad, prevPerfControl_,
                             lastPriv_);
      if (isDataAddressExternal(ldStAddr_))
	pregs.updateCounters(EventNumber::BusLoad, prevPerfControl_,
                             lastPriv_);
    }
  else if (info.isStore())
    {
      pregs.updateCounters(EventNumber::Store, prevPerfControl_,
                           lastPriv_);
      if (misalignedLdSt_)
	pregs.updateCounters(EventNumber::MisalignStore, prevPerfControl_,
                             lastPriv_);
      size_t addr = 0;
      uint64_t value = 0;
      memory_.getLastWriteOldValue(hartIx_, addr, value);
      if (isDataAddressExternal(addr))
	pregs.updateCounters(EventNumber::BusStore, prevPerfControl_,
                             lastPriv_);
    }
  else if (info.isBitManipulation())
    {
      pregs.updateCounters(EventNumber::Bitmanip, prevPerfControl_,
                           lastPriv_);
    }
  else if (info.isAtomic())
    {
      if (id == InstId::lr_w or id == InstId::lr_d)
	pregs.updateCounters(EventNumber::Lr, prevPerfControl_,
                             lastPriv_);
      else if (id == InstId::sc_w or id == InstId::sc_d)
	pregs.updateCounters(EventNumber::Sc, prevPerfControl_,
                             lastPriv_);
      else
	pregs.updateCounters(EventNumber::Atomic, prevPerfControl_,
                             lastPriv_);
    }
  else if (info.isCsr() and not hasException_)
    {
      if ((id == InstId::csrrw or id == InstId::csrrwi))
	{
	  if (op0 == 0)
	    pregs.updateCounters(EventNumber::CsrWrite, prevPerfControl_,
                                 lastPriv_);
	  else
	    pregs.updateCounters(EventNumber::CsrReadWrite, prevPerfControl_,
                                 lastPriv_);
	}
      else
	{
	  if (op1 == 0)
	    pregs.updateCounters(EventNumber::CsrRead, prevPerfControl_,
                                 lastPriv_);
	  else
	    pregs.updateCounters(EventNumber::CsrReadWrite, prevPerfControl_,
                                 lastPriv_);
	}

      // Counter modified by csr instruction should not count up.
      // Also, counter stops counting after corresponding event reg is
      // written.
      std::vector<CsrNumber> csrs;
      std::vector<unsigned> triggers;
      csRegs_.getLastWrittenRegs(csrs, triggers);
      for (auto& csr : csrs)
        {
          auto csrPtr = csRegs_.getImplementedCsr(csr);
          csrPtr->undoCountUp();

          if (csr >= CsrNumber::MHPMEVENT3 and csr <= CsrNumber::MHPMEVENT31)
            if (not csRegs_.applyPerfEventAssign())
              std::cerr << "Unexpected applyPerfAssign fail\n";
        }
    }
  else if (info.isBranch())
    {
      pregs.updateCounters(EventNumber::Branch, prevPerfControl_,
                           lastPriv_);
      if (lastBranchTaken_)
	pregs.updateCounters(EventNumber::BranchTaken, prevPerfControl_,
                             lastPriv_);
    }
}


template <typename URV>
void
Hart<URV>::accumulateInstructionStats(const DecodedInst& di)
{
  const InstEntry& info = *(di.instEntry());

  if (enableCounters_)
    updatePerformanceCounters(di.inst(), info, di.op0(), di.op1());
  prevPerfControl_ = perfControl_;

  // We do not update the instruction stats if an instruction causes
  // an exception unless it is an ebreak or an ecall.
  InstId id = info.instId();
  if (hasException_ and id != InstId::ecall and id != InstId::ebreak and
      id != InstId::c_ebreak)
    return;

  misalignedLdSt_ = false;
  lastBranchTaken_ = false;

  if (not instFreq_)
    return;

  InstProfile& prof = instProfileVec_.at(size_t(id));

  prof.freq_++;
  if (lastPriv_ == PrivilegeMode::User)
    prof.user_++;
  else if (lastPriv_ == PrivilegeMode::Supervisor)
    prof.supervisor_++;
  else if (lastPriv_ == PrivilegeMode::Machine)
    prof.machine_++;

  bool hasRd = false;

  unsigned rs1 = 0, rs2 = 0;
  bool hasRs1 = false, hasRs2 = false;

  if (info.ithOperandType(0) == OperandType::IntReg)
    {
      hasRd = info.isIthOperandWrite(0);
      if (hasRd)
	prof.rd_.at(di.op0())++;
      else
	{
	  rs1 = di.op0();
	  prof.rs1_.at(rs1)++;
	  hasRs1 = true;
	}
    }

  bool hasImm = false;  // True if instruction has an immediate operand.
  int32_t imm = 0;     // Value of immediate operand.

  if (info.ithOperandType(1) == OperandType::IntReg)
    {
      if (hasRd)
	{
	  rs1 = di.op1();
	  prof.rs1_.at(rs1)++;
	  hasRs1 = true;
	}
      else
	{
	  rs2 = di.op1();
	  prof.rs2_.at(rs2)++;
	  hasRs2 = true;
	}
    }
  else if (info.ithOperandType(1) == OperandType::Imm)
    {
      hasImm = true;
      imm = di.op1();
    }

  if (info.ithOperandType(2) == OperandType::IntReg)
    {
      if (hasRd)
	{
	  rs2 = di.op2();
	  prof.rs2_.at(rs2)++;
	  hasRs2 = true;
	}
      else
	assert(0);
    }
  else if (info.ithOperandType(2) == OperandType::Imm)
    {
      hasImm = true;
      imm = di.op2();
    }

  if (hasImm)
    {
      prof.hasImm_ = true;

      if (prof.freq_ == 1)
	{
	  prof.minImm_ = prof.maxImm_ = imm;
	}
      else
	{
	  prof.minImm_ = std::min(prof.minImm_, imm);
	  prof.maxImm_ = std::max(prof.maxImm_, imm);
	}
      addToSignedHistogram(prof.immHisto_, imm);
    }

  unsigned rd = unsigned(intRegCount() + 1);
  URV rdOrigVal = 0;
  intRegs_.getLastWrittenReg(rd, rdOrigVal);

  if (hasRs1)
    {
      URV val1 = intRegs_.read(rs1);
      if (rs1 == rd)
	val1 = rdOrigVal;
      if (info.isUnsigned())
	addToUnsignedHistogram(prof.rs1Histo_, val1);
      else
	addToSignedHistogram(prof.rs1Histo_, SRV(val1));
    }

  if (hasRs2)
    {
      URV val2 = intRegs_.read(rs2);
      if (rs2 == rd)
	val2 = rdOrigVal;
      if (info.isUnsigned())
	addToUnsignedHistogram(prof.rs2Histo_, val2);
      else
	addToSignedHistogram(prof.rs2Histo_, SRV(val2));
    }
}


template <typename URV>
void
Hart<URV>::accumulateTrapStats(bool isNmi)
{
  URV causeVal = 0;
  peekCsr(CsrNumber::MCAUSE, causeVal);

  // If most sig bit of mcause is 1, we have an interrupt.
  bool isInterrupt = causeVal >> (sizeof(causeVal)*8 - 1);

  if (isNmi)
    ;
  else if (isInterrupt)
    {
      causeVal = (causeVal << 1) >> 1;  // Clear most sig bit.
      if (causeVal < interruptStat_.size())
        interruptStat_.at(causeVal)++;
    }
  else
    {
      URV secVal = 0;
      peekCsr(CsrNumber::MSCAUSE, secVal);
      if (causeVal < exceptionStat_.size())
        {
          auto& secCauseVec = exceptionStat_.at(causeVal);
          if (secVal < secCauseVec.size())
            secCauseVec.at(secVal)++;
        }
    }
}


template <typename URV>
inline
void
Hart<URV>::clearTraceData()
{
  intRegs_.clearLastWrittenReg();
  fpRegs_.clearLastWrittenReg();
  csRegs_.clearLastWrittenRegs();
  memory_.clearLastWriteInfo(hartIx_);
}


template <typename URV>
inline
void
Hart<URV>::setTargetProgramBreak(URV addr)
{
  size_t progBreak = addr;

  size_t pageAddr = memory_.getPageStartAddr(addr);
  if (pageAddr != addr)
    progBreak = pageAddr + memory_.pageSize();

  syscall_.setTargetProgramBreak(progBreak);
}


template <typename URV>
inline
bool
Hart<URV>::setTargetProgramArgs(const std::vector<std::string>& args)
{
  URV sp = 0;

  if (not peekIntReg(RegSp, sp))
    return false;

  // Make sp 16-byte aligned.
  if ((sp & 0xf) != 0)
    sp -= (sp & 0xf);

  // Push the arguments on the stack recording their addresses.
  std::vector<URV> addresses;  // Address of the argv strings.
  for (const auto& arg : args)
    {
      sp -= URV(arg.size() + 1);  // Make room for arg and null char.
      addresses.push_back(sp);

      size_t ix = 0;

      for (uint8_t c : arg)
	if (not memory_.poke(sp + ix++, c))
	  return false;

      if (not memory_.poke(sp + ix++, uint8_t(0))) // Null char.
	return false;
    }

  addresses.push_back(0);  // Null pointer at end of argv.

  // Push on stack null for environment and null for aux vector.
  sp -= sizeof(URV);
  if (not memory_.poke(sp, URV(0)))
    return false;
  sp -= sizeof(URV);
  if (not memory_.poke(sp, URV(0)))
    return false;

  // Push argv entries on the stack.
  sp -= URV(addresses.size() + 1) * sizeof(URV); // Make room for argv & argc
  if ((sp & 0xf) != 0)
    sp -= (sp & 0xf);  // Make sp 16-byte aligned.

  URV ix = 1;  // Index 0 is for argc
  for (const auto addr : addresses)
    {
      if (not memory_.poke(sp + ix++*sizeof(URV), addr))
	return false;
    }

  // Put argc on the stack.
  if (not memory_.poke(sp, URV(args.size())))
    return false;

  if (not pokeIntReg(RegSp, sp))
    return false;

  return true;
}


template <typename URV>
URV
Hart<URV>::lastPc() const
{
  return currPc_;
}


template <typename URV>
int
Hart<URV>::lastIntReg() const
{
  return intRegs_.getLastWrittenReg();
}


template <typename URV>
int
Hart<URV>::lastFpReg() const
{
  return fpRegs_.getLastWrittenReg();
}


template <typename URV>
void
Hart<URV>::lastCsr(std::vector<CsrNumber>& csrs,
		   std::vector<unsigned>& triggers) const
{
  csRegs_.getLastWrittenRegs(csrs, triggers);
}


template <typename URV>
void
Hart<URV>::lastMemory(std::vector<size_t>& addresses,
		      std::vector<uint32_t>& words) const
{
  addresses.clear();
  words.clear();

  size_t address = 0;
  uint64_t value;
  unsigned writeSize = memory_.getLastWriteNewValue(hartIx_, address, value);

  if (not writeSize)
    return;

  addresses.clear();
  words.clear();
  addresses.push_back(address);
  words.push_back(uint32_t(value));

  if (writeSize == 8)
    {
      addresses.push_back(address + 4);
      words.push_back(uint32_t(value >> 32));
    }
}


template <typename URV>
void
handleExceptionForGdb(WdRiscv::Hart<URV>& hart, int fd);


// Return true if debug mode is entered and false otherwise.
template <typename URV>
bool
Hart<URV>::takeTriggerAction(FILE* traceFile, URV pc, URV info,
			     uint64_t& counter, bool beforeTiming)
{
  // Check triggers configuration to determine action: take breakpoint
  // exception or enter debugger.

  bool enteredDebug = false;

  if (csRegs_.hasEnterDebugModeTripped())
    {
      enterDebugMode(DebugModeCause::TRIGGER, pc);
      enteredDebug = true;
    }
  else
    {
      bool doingWide = wideLdSt_;
      auto secCause = SecondaryCause::TRIGGER_HIT;
      initiateException(ExceptionCause::BREAKP, pc, info, secCause);
      if (dcsrStep_)
	{
	  enterDebugMode(DebugModeCause::TRIGGER, pc_);
	  enteredDebug = true;
	  enableWideLdStMode(doingWide);
	}
    }

  if (beforeTiming and traceFile)
    {
      uint32_t inst = 0;
      readInst(currPc_, inst);

      std::string instStr;
      printInstTrace(inst, counter, instStr, traceFile);
    }

  return enteredDebug;
}


template <typename URV>
void
Hart<URV>::copyMemRegionConfig(const Hart<URV>& other)
{
  regionHasLocalMem_ = other.regionHasLocalMem_;
  regionHasLocalDataMem_ = other.regionHasLocalDataMem_;
  regionHasDccm_ = other.regionHasDccm_;
  regionHasMemMappedRegs_ = other.regionHasMemMappedRegs_;
  regionHasLocalInstMem_ = other.regionHasLocalInstMem_;
}


// True if keyboard interrupt (user hit control-c) pending.
static std::atomic<bool> userStop = false;

// Negation of the above. Exists for speed (obsessive compulsive
// engineering).
static std::atomic<bool> noUserStop = true;

void
forceUserStop(int)
{
  userStop = true;
  noUserStop = false;
}


static void
clearUserStop()
{
  userStop = false;
  noUserStop = true;
}


/// Install a signal handler for SIGINT (keyboard) interrupts on
/// construction. Restore to previous handlers on destruction. This
/// allows us to catch a control-c typed by the user in the middle of
/// a long-run and return to the top level interactive command
/// processor.
class SignalHandlers
{
public:

  SignalHandlers()
  {
    clearUserStop();
#ifdef __MINGW64__
  __p_sig_fn_t newKbdAction = forceUserStop;
  prevKbdAction_ = signal(SIGINT, newKbdAction);
#else
  struct sigaction newKbdAction;
  memset(&newKbdAction, 0, sizeof(newKbdAction));
  newKbdAction.sa_handler = forceUserStop;
  sigaction(SIGINT, &newKbdAction, &prevKbdAction_);
#endif
  }

  ~SignalHandlers()
  {
#ifdef __MINGW64__
    signal(SIGINT, prevKbdAction_);
#else
    sigaction(SIGINT, &prevKbdAction_, nullptr);
#endif
  }
  
private:

#ifdef __MINGW64__
  __p_sig_fn_t prevKbdAction_ = nullptr;
#else
  struct sigaction prevKbdAction_;
#endif
};



/// Report the number of retired instruction count and the simulation
/// rate.
static void
reportInstsPerSec(uint64_t instCount, double elapsed, bool userStop)
{
  std::lock_guard<std::mutex> guard(stderrMutex);

  std::cout.flush();

  if (userStop)
    std::cerr << "User stop\n";
  std::cerr << "Retired " << instCount << " instruction"
	    << (instCount > 1? "s" : "") << " in "
	    << (boost::format("%.2fs") % elapsed);
  if (elapsed > 0)
    std::cerr << "  " << size_t(double(instCount)/elapsed) << " inst/s";
  std::cerr << '\n';
}


template <typename URV>
bool
Hart<URV>::logStop(const CoreException& ce, uint64_t counter, FILE* traceFile)
{
  std::lock_guard<std::mutex> guard(stderrMutex);

  bool success = false;
  bool isRetired = false;

  if (ce.type() == CoreException::Stop)
    {
      isRetired = true;
      success = ce.value() == 1; // Anything besides 1 is a fail.
      setTargetProgramFinished(true);
    }
  else if (ce.type() == CoreException::Exit)
    {
      isRetired = true;
      success = ce.value() == 0;
      setTargetProgramFinished(true);
    }

  if (isRetired)
    {
      if (minstretEnabled())
        retiredInsts_++;
      if (traceFile)
	{
	  uint32_t inst = 0;
	  readInst(currPc_, inst);
	  std::string instStr;
	  printInstTrace(inst, counter, instStr, traceFile);
	}
    }

  using std::cerr;

  cerr << std::dec;
  if (ce.type() == CoreException::Stop)
    cerr << (success? "Successful " : "Error: Failed ")
         << "stop: " << ce.what() << ": " << ce.value() << "\n";
  else if (ce.type() == CoreException::Exit)
    cerr << "Target program exited with code " << ce.value() << '\n';
  else
    cerr << "Stopped -- unexpected exception\n";

  return success;
}


static
bool
isInputPending(int fd)
{
  struct pollfd pfds[1];
  pfds[0].fd = fd;
  pfds[0].events = POLLIN;
  pfds[0].revents = 0;
  if (poll(pfds, 1, 0) == 1)
    return pfds[0].revents & POLLIN;
  return false;
}


template <typename URV>
inline
bool
Hart<URV>::fetchInstWithTrigger(URV addr, uint32_t& inst, FILE* file)
{
  // Process pre-execute address trigger and fetch instruction.
  bool hasTrig = hasActiveInstTrigger();
  triggerTripped_ = (hasTrig and
                     instAddrTriggerHit(addr, TriggerTiming::Before,
                                        privMode_, isInterruptEnabled()));
  // Fetch instruction.
  bool fetchOk = true;
  if (triggerTripped_)
    {
      if (not fetchInstPostTrigger(addr, inst, file))
        {
          ++cycleCount_;
          return false;  // Next instruction in trap handler.
        }
    }
  else
    {
      uint32_t ix = (addr >> 1) & decodeCacheMask_;
      DecodedInst* di = &decodeCache_[ix];
      if (not di->isValid() or di->address() != pc_)
        fetchOk = fetchInst(addr, inst);
      else
        inst = di->inst();
    }
  if (not fetchOk)
    {
      ++cycleCount_;
      if (file)
        {
          std::string instStr;
          printInstTrace(inst, instCounter_, instStr, file);
        }

      if (dcsrStep_)
        enterDebugMode(DebugModeCause::STEP, pc_);

      return false;  // Next instruction in trap handler.
    }

  // Process pre-execute opcode trigger.
  if (hasTrig and instOpcodeTriggerHit(inst, TriggerTiming::Before,
                                       privMode_,
                                       isInterruptEnabled()))
    triggerTripped_ = true;

  return true;
}


template <typename URV>
bool
Hart<URV>::untilAddress(size_t address, FILE* traceFile)
{
  std::string instStr;
  instStr.reserve(128);

  // Need csr history when tracing or for triggers
  bool trace = traceFile != nullptr or enableTriggers_;
  clearTraceData();

  uint64_t limit = instCountLim_;
  bool doStats = instFreq_ or enableCounters_;

  // Check for gdb break every 1000000 instructions.
  unsigned gdbCount = 0, gdbLimit = 1000000;

  if (enableGdb_)
    handleExceptionForGdb(*this, gdbInputFd_);

  while (pc_ != address and instCounter_ < limit)
    {
      if (userStop)
        break;

      if (enableGdb_ and ++gdbCount >= gdbLimit)
        {
          gdbCount = 0;
          if (isInputPending(gdbInputFd_))
            {
              handleExceptionForGdb(*this, gdbInputFd_);
              continue;
            }
        }

      try
	{
          uint32_t inst = 0;
	  currPc_ = pc_;

	  ldStAddrValid_ = false;
	  triggerTripped_ = false;
	  hasException_ = false;
          lastPriv_ = privMode_;

	  ++instCounter_;

          if (processExternalInterrupt(traceFile, instStr))
            continue;

          if (not fetchInstWithTrigger(pc_, inst, traceFile))
            {
	      clearTraceData();
              continue;  // Next instruction in trap handler.
            }

	  // Decode unless match in decode cache.
	  uint32_t ix = (pc_ >> 1) & decodeCacheMask_;
	  DecodedInst* di = &decodeCache_[ix];
	  if (not di->isValid() or di->address() != pc_)
	    decode(pc_, inst, *di);

	  bool doingWide = wideLdSt_;

          // Increment pc and execute instruction
	  pc_ += di->instSize();
	  execute(di);

	  ++cycleCount_;

	  if (hasException_)
	    {
              if (doStats)
                accumulateInstructionStats(*di);
	      if (traceFile)
		{
		  printInstTrace(*di, instCounter_, instStr, traceFile);
		  clearTraceData();
		}
	      continue;
	    }

	  if (triggerTripped_)
	    {
	      undoForTrigger();
	      if (takeTriggerAction(traceFile, currPc_, currPc_, instCounter_, true))
		return true;
	      continue;
	    }

	  if (doingWide)
	    enableWideLdStMode(false);

          if (minstretEnabled())
            ++retiredInsts_;

	  if (doStats)
	    accumulateInstructionStats(*di);

	  if (trace)
	    {
	      if (traceFile)
		printInstTrace(*di, instCounter_, instStr, traceFile);
	      clearTraceData();
	    }

	  bool icountHit = (enableTriggers_ and
			    icountTriggerHit(privMode_, isInterruptEnabled()));
	  if (icountHit)
	    if (takeTriggerAction(traceFile, pc_, pc_, instCounter_, false))
	      return true;
          prevPerfControl_ = perfControl_;
	}
      catch (const CoreException& ce)
	{
	  return logStop(ce, instCounter_, traceFile);
	}
    }

  return true;
}


template <typename URV>
bool
Hart<URV>::runUntilAddress(size_t address, FILE* traceFile)
{
  struct timeval t0;
  gettimeofday(&t0, nullptr);

  uint64_t limit = instCountLim_;
  uint64_t counter0 = instCounter_;

  // Setup signal handlers. Restore on destruction.
  SignalHandlers handlers();

  bool success = untilAddress(address, traceFile);
      
  if (instCounter_ == limit)
    std::cerr << "Stopped -- Reached instruction limit\n";
  else if (pc_ == address)
    std::cerr << "Stopped -- Reached end address\n";

  // Simulator stats.
  struct timeval t1;
  gettimeofday(&t1, nullptr);
  double elapsed = (double(t1.tv_sec - t0.tv_sec) +
		    double(t1.tv_usec - t0.tv_usec)*1e-6);

  uint64_t numInsts = instCounter_ - counter0;

  reportInstsPerSec(numInsts, elapsed, userStop);
  return success;
}


template <typename URV>
bool
Hart<URV>::simpleRun()
{
  // For speed: do not record/clear CSR changes.
  enableCsrTrace_ = false;

  bool success = true;

  try
    {
      while (true)
        {
          bool hasLim = (instCountLim_ < ~uint64_t(0));
          if (hasLim)
            simpleRunWithLimit();
          else
            simpleRunNoLimit();

          if (userStop)
            {
              std::cerr << "Stopped -- interrupted\n";
              break;
            }

          if (hasLim)
            std::cerr << "Stopped -- Reached instruction limit\n";
          break;
        }
    }
  catch (const CoreException& ce)
    {
      success = logStop(ce, 0, nullptr);
    }

  enableCsrTrace_ = true;

  return success;
}


template <typename URV>
bool
Hart<URV>::simpleRunWithLimit()
{
  uint64_t limit = instCountLim_;
  while (noUserStop and instCounter_ < limit) 
    {
      currPc_ = pc_;
      ++instCounter_;

      // Fetch/decode unless match in decode cache.
      uint32_t ix = (pc_ >> 1) & decodeCacheMask_;
      DecodedInst* di = &decodeCache_[ix];
      if (not di->isValid() or di->address() != pc_)
        {
          uint32_t inst = 0;
          if (not fetchInst(pc_, inst))
            continue;
          decode(pc_, inst, *di);
        }

      pc_ += di->instSize();
      execute(di);
    }
  return true;
}


template <typename URV>
bool
Hart<URV>::simpleRunNoLimit()
{
  while (noUserStop) 
    {
      currPc_ = pc_;
      ++instCounter_;

      // Fetch/decode unless match in decode cache.
      uint32_t ix = (pc_ >> 1) & decodeCacheMask_;
      DecodedInst* di = &decodeCache_[ix];
      if (not di->isValid() or di->address() != pc_)
        {
          uint32_t inst = 0;
          if (not fetchInst(pc_, inst))
            continue;
          decode(pc_, inst, *di);
        }

      pc_ += di->instSize();
      execute(di);
    }

  return true;
}


template <typename URV>
bool
Hart<URV>::openTcpForGdb()
{
  struct sockaddr_in address;
  int opt = 1;
  int addrlen = sizeof(address);
  address.sin_family = AF_INET;
  address.sin_addr.s_addr = INADDR_ANY;
  address.sin_port = htons( gdbTcpPort_ );
  bool succ = true;
  int gdbFd = socket(AF_INET, SOCK_STREAM, 0);
  succ = (gdbFd > 0) and not setsockopt(gdbFd, SOL_SOCKET,
                                        SO_REUSEADDR | SO_REUSEPORT, &opt,
                                        sizeof(opt));

  succ = succ and bind(gdbFd, (struct sockaddr *)&address,sizeof(address))>=0;
  succ = succ and listen(gdbFd, 3) >= 0;
  if (succ)
    {
      gdbInputFd_ = accept(gdbFd, (struct sockaddr*) &address,
                           (socklen_t*) &addrlen);
      succ = gdbInputFd_ >= 0;
    }
  return succ;

}

/// Run indefinitely.  If the tohost address is defined, then run till
/// a write is attempted to that address.
template <typename URV>
bool
Hart<URV>::run(FILE* file)
{
  if (gdbTcpPort_ >= 0)
    openTcpForGdb();
  else if (enableGdb_)
    gdbInputFd_ = STDIN_FILENO;

  // To run fast, this method does not do much besides
  // straight-forward execution. If any option is turned on, we switch
  // to runUntilAdress which supports all features.
  URV stopAddr = stopAddrValid_? stopAddr_ : ~URV(0); // ~URV(0): No-stop PC.
  bool hasWideLdSt = csRegs_.isImplemented(CsrNumber::MDBAC);
  bool hasClint = clintStart_ < clintLimit_;
  bool complex = (stopAddrValid_ or instFreq_ or enableTriggers_ or enableGdb_
                  or enableCounters_ or alarmInterval_ or file or hasWideLdSt
                  or hasClint or isRvs());
  if (complex)
    return runUntilAddress(stopAddr, file); 

  uint64_t counter0 = instCounter_;

  struct timeval t0;
  gettimeofday(&t0, nullptr);

  // Setup signal handlers. Restore on destruction.
  SignalHandlers handlers();

  bool success = simpleRun();

  // Simulator stats.
  struct timeval t1;
  gettimeofday(&t1, nullptr);
  double elapsed = (double(t1.tv_sec - t0.tv_sec) +
		    double(t1.tv_usec - t0.tv_usec)*1e-6);

  uint64_t numInsts = instCounter_ - counter0;
  reportInstsPerSec(numInsts, elapsed, userStop);
  return success;
}


template <typename URV>
bool
Hart<URV>::isInterruptPossible(InterruptCause& cause)
{
  if (debugMode_ and not debugStepMode_)
    return false;

  URV mip = csRegs_.peekMip();
  URV mie = csRegs_.peekMie();
  if ((mie & mip) == 0)
    return false;  // Nothing enabled that is also pending.

  typedef InterruptCause IC;

  // Check for machine-level interrupts if MIE enabled or if user/supervisor.
  URV mstatus = csRegs_.peekMstatus();
  MstatusFields<URV> fields(mstatus);
  bool globalEnable = fields.bits_.MIE or privMode_ < PrivilegeMode::Machine;
  URV delegVal = 0;
  peekCsr(CsrNumber::MIDELEG, delegVal);
  for (InterruptCause ic : { IC::M_EXTERNAL, IC::M_LOCAL, IC::M_SOFTWARE,
                              IC::M_TIMER, IC::M_INT_TIMER0,
                              IC::M_INT_TIMER1 } )
    {
      URV mask = URV(1) << unsigned(ic);
      bool delegated = (mask & delegVal) != 0;
      bool enabled = globalEnable;
      if (delegated)
        enabled = ((privMode_ == PrivilegeMode::Supervisor and fields.bits_.SIE) or
                   privMode_ < PrivilegeMode::Supervisor);
      if (enabled)
        if (mie & mask & mip)
          {
            cause = ic;
            if (ic == IC::M_TIMER and alarmInterval_ > 0)
              {
                // Reset the timer-interrupt pending bit.
                mip = mip & ~mask;
                pokeCsr(CsrNumber::MIP, mip);
              }
            return true;
          }
    }

  // Supervisor mode interrupts: SIE enabled and supervior mode, or user-mode.
  if (isRvs())
    {
      bool check = ((fields.bits_.SIE and privMode_ == PrivilegeMode::Supervisor)
                    or privMode_ < PrivilegeMode::Supervisor);
      if (check)
        for (auto ic : { IC::S_EXTERNAL, IC::S_SOFTWARE, IC::S_TIMER } )
          {
            URV mask = URV(1) << unsigned(ic);
            if (mie & mask & mip)
              {
                cause = ic;
                return true;
              }
          }
    }

  // User mode interrupts: UIE enabled and user-mode.
  if (isRvu())
    {
      bool check = fields.bits_.UIE and privMode_ == PrivilegeMode::User;
      if (check)
        for (auto ic : { IC::U_EXTERNAL, IC::U_SOFTWARE, IC::U_TIMER } )
          {
            URV mask = URV(1) << unsigned(ic);
            if (mie & mask & mip)
              {
                cause = ic;
                return true;
              }
          }
    }

  return false;
}


template <typename URV>
bool
Hart<URV>::processExternalInterrupt(FILE* traceFile, std::string& instStr)
{
  if (instCounter_ >= alarmLimit_)
    {
      URV mipVal = csRegs_.peekMip();
      mipVal = mipVal | (URV(1) << URV(InterruptCause::M_TIMER));
      csRegs_.poke(CsrNumber::MIP, mipVal);
      alarmLimit_ += alarmInterval_;
    }

  if (debugStepMode_ and not dcsrStepIe_)
    return false;

  // If a non-maskable interrupt was signaled by the test-bench, take it.
  if (nmiPending_)
    {
      initiateNmi(URV(nmiCause_), pc_);
      nmiPending_ = false;
      nmiCause_ = NmiCause::UNKNOWN;
      uint32_t inst = 0; // Load interrupted inst.
      readInst(currPc_, inst);
      if (traceFile)  // Trace interrupted instruction.
	printInstTrace(inst, instCounter_, instStr, traceFile, true);
      return true;
    }

  // If interrupts enabled and one is pending, take it.
  InterruptCause cause;
  if (isInterruptPossible(cause))
    {
      // Attach changes to interrupted instruction.
      initiateInterrupt(cause, pc_);
      uint32_t inst = 0; // Load interrupted inst.
      readInst(currPc_, inst);
      if (traceFile)  // Trace interrupted instruction.
	printInstTrace(inst, instCounter_, instStr, traceFile, true);
      ++cycleCount_;
      return true;
    }
  return false;
}


template <typename URV>
void
Hart<URV>::invalidateDecodeCache(URV addr, unsigned storeSize)
{
  // Consider putting this in a callback associated with memory
  // write/poke. This way it can be applied only to pages marked
  // execute.

  // We want to check the location before the address just in case it
  // contains a 4-byte instruction that overlaps what was written.
  storeSize += 3;
  addr -= 3;

  for (unsigned i = 0; i < storeSize; i += 2)
    {
      URV instAddr = (addr + i) >> 1;
      uint32_t cacheIx = instAddr & decodeCacheMask_;
      auto& entry = decodeCache_[cacheIx];
      if ((entry.address() >> 1) == instAddr)
	entry.invalidate();
    }
}


template <typename URV>
void
Hart<URV>::invalidateDecodeCache()
{
  for (auto& entry : decodeCache_)
    entry.invalidate();
}


template <typename URV>
void
Hart<URV>::singleStep(FILE* traceFile)
{
  std::string instStr;

  // Single step is mostly used for follow-me mode where we want to
  // know the changes after the execution of each instruction.
  bool doStats = instFreq_ or enableCounters_;

  try
    {
      uint32_t inst = 0;
      currPc_ = pc_;

      ldStAddrValid_ = false;
      triggerTripped_ = false;
      hasException_ = false;
      ebreakInstDebug_ = false;
      lastPriv_ = privMode_;

      ++instCounter_;

      if (processExternalInterrupt(traceFile, instStr))
	return;  // Next instruction in interrupt handler.

      if (not fetchInstWithTrigger(pc_, inst, traceFile))
        return;

      DecodedInst di;
      decode(pc_, inst, di);

      bool doingWide = wideLdSt_;

      // Increment pc and execute instruction
      pc_ += di.instSize();
      execute(&di);

      ++cycleCount_;

      // A ld/st must be seen within 2 steps of a forced access fault.
      if (forceAccessFail_ and (instCounter_ > forceAccessFailMark_ + 1))
	{
	  std::cerr << "Spurious exception command from test-bench.\n";
	  forceAccessFail_ = false;
	}

      if (hasException_)
	{
	  if (doStats)
	    accumulateInstructionStats(di);
	  if (traceFile)
	    printInstTrace(di, instCounter_, instStr, traceFile);
	  if (dcsrStep_ and not ebreakInstDebug_)
	    enterDebugMode(DebugModeCause::STEP, pc_);
	  return;
	}

      if (triggerTripped_)
	{
	  undoForTrigger();
	  takeTriggerAction(traceFile, currPc_, currPc_, instCounter_, true);
	  return;
	}

      if (doingWide)
	enableWideLdStMode(false);

      if (minstretEnabled() and not ebreakInstDebug_)
        ++retiredInsts_;

      if (doStats)
	accumulateInstructionStats(di);

      if (traceFile)
	printInstTrace(inst, instCounter_, instStr, traceFile);

      // If a register is used as a source by an instruction then any
      // pending load with same register as target is removed from the
      // load queue (because in such a case the hardware will stall
      // till load is completed). Source operands of load instructions
      // are handled in the load and loadRserve methods.
      const InstEntry* entry = di.instEntry();
      if (not entry->isLoad())
	{
	  if (entry->isIthOperandIntRegSource(0))
	    removeFromLoadQueue(di.op0(), entry->isDivide());
	  if (entry->isIthOperandIntRegSource(1))
	    removeFromLoadQueue(di.op1(), entry->isDivide());
	  if (entry->isIthOperandIntRegSource(2))
	    removeFromLoadQueue(di.op2(), entry->isDivide());

	  // If a register is written by a non-load instruction, then
	  // its entry is invalidated in the load queue.
	  int regIx = intRegs_.getLastWrittenReg();
	  if (regIx > 0)
	    invalidateInLoadQueue(regIx, entry->isDivide());
	}

      bool icountHit = (enableTriggers_ and 
			icountTriggerHit(privMode_, isInterruptEnabled()));
      if (icountHit)
	{
	  takeTriggerAction(traceFile, pc_, pc_, instCounter_, false);
	  return;
	}

      // If step bit set in dcsr then enter debug mode unless already there.
      if (dcsrStep_ and not ebreakInstDebug_)
	enterDebugMode(DebugModeCause::STEP, pc_);

      prevPerfControl_ = perfControl_;
    }
  catch (const CoreException& ce)
    {
      // If step bit set in dcsr then enter debug mode unless already there.
      // This is for the benefit of the test bench.
      if (dcsrStep_ and not ebreakInstDebug_)
	enterDebugMode(DebugModeCause::STEP, pc_);

      logStop(ce, instCounter_, traceFile);
    }
}


template <typename URV>
void
Hart<URV>::postDataAccessFault(URV offset, SecondaryCause secCause)
{
  forceAccessFail_ = true;
  forcedCause_ = secCause;
  forceAccessFailOffset_ = offset;
  forceAccessFailMark_ = instCounter_;
}


template <typename URV>
bool
Hart<URV>::whatIfSingleStep(uint32_t inst, ChangeRecord& record)
{
  uint64_t prevExceptionCount = exceptionCount_;
  URV prevPc = pc_;

  clearTraceData();
  triggerTripped_ = false;

  // Note: triggers not yet supported.

  DecodedInst di;
  decode(pc_, inst, di);

  // Execute instruction
  pc_ += di.instSize();
  execute(&di);

  bool result = exceptionCount_ == prevExceptionCount;

  // If step bit set in dcsr then enter debug mode unless already there.
  if (dcsrStep_ and not ebreakInstDebug_)
    enterDebugMode(DebugModeCause::STEP, pc_);

  // Collect changes. Undo each collected change.
  exceptionCount_ = prevExceptionCount;

  collectAndUndoWhatIfChanges(prevPc, record);
  
  return result;
}


template <typename URV>
bool
Hart<URV>::whatIfSingleStep(URV whatIfPc, uint32_t inst, ChangeRecord& record)
{
  URV prevPc = pc_;
  pc_ = whatIfPc;

  // Note: triggers not yet supported.
  triggerTripped_ = false;

  // Fetch instruction. We don't care about what we fetch. Just checking
  // if there is a fetch exception.
  uint32_t dummyInst = 0;
  bool fetchOk = fetchInst(pc_, dummyInst);

  if (not fetchOk)
    {
      collectAndUndoWhatIfChanges(prevPc, record);
      return false;
    }

  bool res = whatIfSingleStep(inst, record);

  pc_ = prevPc;
  return res;
}


template <typename URV>
bool
Hart<URV>::whatIfSingStep(const DecodedInst& di, ChangeRecord& record)
{
  clearTraceData();
  uint64_t prevExceptionCount = exceptionCount_;
  URV prevPc  = pc_, prevCurrPc = currPc_;

  currPc_ = pc_ = di.address();

  // Note: triggers not yet supported.
  triggerTripped_ = false;

  // Save current value of operands.
  uint64_t prevRegValues[4];
  for (unsigned i = 0; i < 4; ++i)
    {
      URV prev = 0;
      prevRegValues[i] = 0;
      uint32_t operand = di.ithOperand(i);

      switch (di.ithOperandType(i))
	{
	case OperandType::None:
	case OperandType::Imm:
	  break;
	case OperandType::IntReg:
	  peekIntReg(operand, prev);
	  prevRegValues[i] = prev;
	  break;
	case OperandType::FpReg:
	  peekFpReg(operand, prevRegValues[i]);
	  break;
	case OperandType::CsReg:
	  peekCsr(CsrNumber(operand), prev);
	  prevRegValues[i] = prev;
	  break;
	}
    }

  // Temporarily set value of operands to what-if values.
  for (unsigned i = 0; i < 4; ++i)
    {
      uint32_t operand = di.ithOperand(i);

      switch (di.ithOperandType(i))
	{
	case OperandType::None:
	case OperandType::Imm:
	  break;
	case OperandType::IntReg:
	  pokeIntReg(operand, di.ithOperandValue(i));
	  break;
	case OperandType::FpReg:
	  pokeFpReg(operand, di.ithOperandValue(i));
	  break;
	case OperandType::CsReg:
	  pokeCsr(CsrNumber(operand), di.ithOperandValue(i));
	  break;
	}
    }

  // Execute instruction.
  pc_ += di.instSize();
  if(di.instEntry()->instId() != InstId::illegal)
	  execute(&di);
  bool result = exceptionCount_ == prevExceptionCount;

  // Collect changes. Undo each collected change.
  exceptionCount_ = prevExceptionCount;
  collectAndUndoWhatIfChanges(prevPc, record);

  // Restore temporarily modified registers.
  for (unsigned i = 0; i < 4; ++i)
    {
      uint32_t operand = di.ithOperand(i);

      switch (di.ithOperandType(i))
	{
	case OperandType::None:
	case OperandType::Imm:
	  break;
	case OperandType::IntReg:
	  pokeIntReg(operand, prevRegValues[i]);
	  break;
	case OperandType::FpReg:
	  pokeFpReg(operand, prevRegValues[i]);
	  break;
	case OperandType::CsReg:
	  pokeCsr(CsrNumber(operand), prevRegValues[i]);
	  break;
	}
    }

  pc_ = prevPc;
  currPc_ = prevCurrPc;

  return result;
}


template <typename URV>
void
Hart<URV>::collectAndUndoWhatIfChanges(URV prevPc, ChangeRecord& record)
{
  record.clear();

  record.newPc = pc_;
  pc_ = prevPc;

  unsigned regIx = 0;
  URV oldValue = 0;
  if (intRegs_.getLastWrittenReg(regIx, oldValue))
    {
      URV newValue = 0;
      peekIntReg(regIx, newValue);
      pokeIntReg(regIx, oldValue);

      record.hasIntReg = true;
      record.intRegIx = regIx;
      record.intRegValue = newValue;
    }

  uint64_t oldFpValue = 0;
  if (fpRegs_.getLastWrittenReg(regIx, oldFpValue))
    {
      uint64_t newFpValue = 0;
      peekFpReg(regIx, newFpValue);
      pokeFpReg(regIx, oldFpValue);

      record.hasFpReg = true;
      record.fpRegIx = regIx;
      record.fpRegValue = newFpValue;
    }

  record.memSize = memory_.getLastWriteNewValue(hartIx_, record.memAddr,
                                                record.memValue);

  size_t addr = 0;
  uint64_t value = 0;
  size_t byteCount = memory_.getLastWriteOldValue(hartIx_, addr, value);
  for (size_t i = 0; i < byteCount; ++i)
    {
      uint8_t byte = value & 0xff;
      memory_.poke(addr, byte);
      addr++;
      value = value >> 8;
    }

  std::vector<CsrNumber> csrNums;
  std::vector<unsigned> triggerNums;
  csRegs_.getLastWrittenRegs(csrNums, triggerNums);

  for (auto csrn : csrNums)
    {
      Csr<URV>* csr = csRegs_.getImplementedCsr(csrn);
      if (not csr)
	continue;

      URV newVal = csr->read();
      URV oldVal = csr->prevValue();
      csr->write(oldVal);

      record.csrIx.push_back(csrn);
      record.csrValue.push_back(newVal);
    }

  clearTraceData();
}


template <typename URV>
void
Hart<URV>::setFcsrFlags(FpFlags flag)
{
  auto prev = getFpFlags();
  auto val = prev | unsigned(flag);
  if (val != prev)
    {
      setFpFlags(val);
      recordCsrWrite(CsrNumber::FCSR);
    }
}


template <typename URV>
inline
void
Hart<URV>::execLui(const DecodedInst* di)
{
  intRegs_.write(di->op0(), SRV(int32_t(di->op1())));
}


template <typename URV>
void
Hart<URV>::execute(const DecodedInst* di)
{
#pragma GCC diagnostic ignored "-Wpedantic"

  // Setup an array of labels to index it by the decoded instruction
  // id to jump to execute function of that instruction. A table of
  // methods would be a lot cleaner (no goto) but it would cost 20%
  // more in execution time.
  static void* labels[] =
    {
     &&illegal,
     &&lui,
     &&auipc,
     &&jal,
     &&jalr,
     &&beq,
     &&bne,
     &&blt,
     &&bge,
     &&bltu,
     &&bgeu,
     &&lb,
     &&lh,
     &&lw,
     &&lbu,
     &&lhu,
     &&sb,
     &&sh,
     &&sw,
     &&addi,
     &&slti,
     &&sltiu,
     &&xori,
     &&ori,
     &&andi,
     &&slli,
     &&srli,
     &&srai,
     &&add,
     &&sub,
     &&sll,
     &&slt,
     &&sltu,
     &&xor_,
     &&srl,
     &&sra,
     &&or_,
     &&and_,
     &&fence,
     &&fencei,
     &&ecall,
     &&ebreak,
     &&csrrw,
     &&csrrs,
     &&csrrc,
     &&csrrwi,
     &&csrrsi,
     &&csrrci,
     &&lwu,
     &&ld,
     &&sd,
     &&addiw,
     &&slliw,
     &&srliw,
     &&sraiw,
     &&addw,
     &&subw,
     &&sllw,
     &&srlw,
     &&sraw,
     &&mul,
     &&mulh,
     &&mulhsu,
     &&mulhu,
     &&div,
     &&divu,
     &&rem,
     &&remu,
     &&mulw,
     &&divw,
     &&divuw,
     &&remw,
     &&remuw,
     &&lr_w,
     &&sc_w,
     &&amoswap_w,
     &&amoadd_w,
     &&amoxor_w,
     &&amoand_w,
     &&amoor_w,
     &&amomin_w,
     &&amomax_w,
     &&amominu_w,
     &&amomaxu_w,
     &&lr_d,
     &&sc_d,
     &&amoswap_d,
     &&amoadd_d,
     &&amoxor_d,
     &&amoand_d,
     &&amoor_d,
     &&amomin_d,
     &&amomax_d,
     &&amominu_d,
     &&amomaxu_d,
     &&flw,
     &&fsw,
     &&fmadd_s,
     &&fmsub_s,
     &&fnmsub_s,
     &&fnmadd_s,
     &&fadd_s,
     &&fsub_s,
     &&fmul_s,
     &&fdiv_s,
     &&fsqrt_s,
     &&fsgnj_s,
     &&fsgnjn_s,
     &&fsgnjx_s,
     &&fmin_s,
     &&fmax_s,
     &&fcvt_w_s,
     &&fcvt_wu_s,
     &&fmv_x_w,
     &&feq_s,
     &&flt_s,
     &&fle_s,
     &&fclass_s,
     &&fcvt_s_w,
     &&fcvt_s_wu,
     &&fmv_w_x,
     &&fcvt_l_s,
     &&fcvt_lu_s,
     &&fcvt_s_l,
     &&fcvt_s_lu,
     &&fld,
     &&fsd,
     &&fmadd_d,
     &&fmsub_d,
     &&fnmsub_d,
     &&fnmadd_d,
     &&fadd_d,
     &&fsub_d,
     &&fmul_d,
     &&fdiv_d,
     &&fsqrt_d,
     &&fsgnj_d,
     &&fsgnjn_d,
     &&fsgnjx_d,
     &&fmin_d,
     &&fmax_d,
     &&fcvt_s_d,
     &&fcvt_d_s,
     &&feq_d,
     &&flt_d,
     &&fle_d,
     &&fclass_d,
     &&fcvt_w_d,
     &&fcvt_wu_d,
     &&fcvt_d_w,
     &&fcvt_d_wu,
     &&fcvt_l_d,
     &&fcvt_lu_d,
     &&fmv_x_d,
     &&fcvt_d_l,
     &&fcvt_d_lu,
     &&fmv_d_x,
     &&mret,
     &&uret,
     &&sret,
     &&wfi,
     &&sfence_vma,
     &&c_addi4spn,
     &&c_fld,
     &&c_lq,
     &&c_lw,
     &&c_flw,
     &&c_ld,
     &&c_fsd,
     &&c_sq,
     &&c_sw,
     &&c_fsw,
     &&c_sd,
     &&c_addi,
     &&c_jal,
     &&c_li,
     &&c_addi16sp,
     &&c_lui,
     &&c_srli,
     &&c_srli64,
     &&c_srai,
     &&c_srai64,
     &&c_andi,
     &&c_sub,
     &&c_xor,
     &&c_or,
     &&c_and,
     &&c_subw,
     &&c_addw,
     &&c_j,
     &&c_beqz,
     &&c_bnez,
     &&c_slli,
     &&c_slli64,
     &&c_fldsp,
     &&c_lwsp,
     &&c_flwsp,
     &&c_ldsp,
     &&c_jr,
     &&c_mv,
     &&c_ebreak,
     &&c_jalr,
     &&c_add,
     &&c_fsdsp,
     &&c_swsp,
     &&c_fswsp,
     &&c_addiw,
     &&c_sdsp,
     &&clz,
     &&ctz,
     &&pcnt,
     &&andn,
     &&orn,
     &&xnor,
     &&slo,
     &&sro,
     &&sloi,
     &&sroi,
     &&min,
     &&max,
     &&minu,
     &&maxu,
     &&rol,
     &&ror,
     &&rori,
     &&rev8,
     &&rev,
     &&pack,
     &&addwu,
     &&subwu,
     &&addiwu,
     &&sext_b,
     &&sext_h,
     &&addu_w,
     &&subu_w,
     &&slliu_w,
     &&packh,
     &&packu,
     &&packw,
     &&packuw,
     &&grev,
     &&grevi,
     &&gorc,
     &&gorci,
     &&shfl,
     &&shfli,
     &&unshfl,
     &&unshfli,
     &&sbset,
     &&sbclr,
     &&sbinv,
     &&sbext,
     &&sbseti,
     &&sbclri,
     &&sbinvi,
     &&sbexti,
     &&bdep,
     &&bext,
     &&bfp,
     &&clmul,
     &&clmulh,
     &&clmulr,
     &&sh1add,
     &&sh2add,
     &&sh3add,
     &&sh1addu_w,
     &&sh2addu_w,
     &&sh3addu_w,

     // zbr
     &&crc32_b,
     &&crc32_h,
     &&crc32_w,
     &&crc32_d,
     &&crc32c_b,
     &&crc32c_h,
     &&crc32c_w,
     &&crc32c_d,

     // zbm
     &&bmator,
     &&bmatxor,
     &&bmatflip,

     // zbt
     &&cmov,
     &&cmix,
     &&fsl,
     &&fsr,
     &&fsri
    };

  const InstEntry* entry = di->instEntry();
  size_t id = size_t(entry->instId());
  assert(id < sizeof(labels));
  goto *labels[id];

 illegal:
  illegalInst(di);
  return;

 lui:
  execLui(di);
  return;

 auipc:
  execAuipc(di);
  return;

 jal:
  execJal(di);
  return;

 jalr:
  execJalr(di);
  return;

 beq:
  execBeq(di);
  return;

 bne:
  execBne(di);
  return;

 blt:
  execBlt(di);
  return;

 bge:
  execBge(di);
  return;

 bltu:
  execBltu(di);
  return;

 bgeu:
  execBgeu(di);
  return;

 lb:
  execLb(di);
  return;

 lh:
  execLh(di);
  return;

 lw:
  execLw(di);
  return;

 lbu:
  execLbu(di);
  return;

 lhu:
  execLhu(di);
  return;

 sb:
  execSb(di);
  return;

 sh:
  execSh(di);
  return;

 sw:
  execSw(di);
  return;

 addi:
  execAddi(di);
  return;

 slti:
  execSlti(di);
  return;

 sltiu:
  execSltiu(di);
  return;

 xori:
  execXori(di);
  return;

 ori:
  execOri(di);
  return;

 andi:
  execAndi(di);
  return;

 slli:
  execSlli(di);
  return;

 srli:
  execSrli(di);
  return;

 srai:
  execSrai(di);
  return;

 add:
  execAdd(di);
  return;

 sub:
  execSub(di);
  return;

 sll:
  execSll(di);
  return;

 slt:
  execSlt(di);
  return;

 sltu:
  execSltu(di);
  return;

 xor_:
  execXor(di);
  return;

 srl:
  execSrl(di);
  return;

 sra:
  execSra(di);
  return;

 or_:
  execOr(di);
  return;

 and_:
  execAnd(di);
  return;

 fence:
  execFence(di);
  return;

 fencei:
  execFencei(di);
  return;

 ecall:
  execEcall(di);
  return;

 ebreak:
  execEbreak(di);
  return;

 csrrw:
  execCsrrw(di);
  return;

 csrrs:
  execCsrrs(di);
  return;

 csrrc:
  execCsrrc(di);
  return;

 csrrwi:
  execCsrrwi(di);
  return;

 csrrsi:
  execCsrrsi(di);
  return;

 csrrci:
  execCsrrci(di);
  return;

 lwu:
  execLwu(di);
  return;

 ld:
  execLd(di);
  return;

 sd:
  execSd(di);
  return;

 addiw:
  execAddiw(di);
  return;

 slliw:
  execSlliw(di);
  return;

 srliw:
  execSrliw(di);
  return;

 sraiw:
  execSraiw(di);
  return;

 addw:
  execAddw(di);
  return;

 subw:
  execSubw(di);
  return;

 sllw:
  execSllw(di);
  return;

 srlw:
  execSrlw(di);
  return;

 sraw:
  execSraw(di);
  return;

 mul:
  execMul(di);
  return;

 mulh:
  execMulh(di);
  return;

 mulhsu:
  execMulhsu(di);
  return;

 mulhu:
  execMulhu(di);
  return;

 div:
  execDiv(di);
  return;

 divu:
  execDivu(di);
  return;

 rem:
  execRem(di);
  return;

 remu:
  execRemu(di);
  return;

 mulw:
  execMulw(di);
  return;

 divw:
  execDivw(di);
  return;

 divuw:
  execDivuw(di);
  return;

 remw:
  execRemw(di);
  return;

 remuw:
  execRemuw(di);
  return;

 lr_w:
  execLr_w(di);
  return;

 sc_w:
  execSc_w(di);
  return;

 amoswap_w:
  execAmoswap_w(di);
  return;

 amoadd_w:
  execAmoadd_w(di);
  return;

 amoxor_w:
  execAmoxor_w(di);
  return;

 amoand_w:
  execAmoand_w(di);
  return;

 amoor_w:
  execAmoor_w(di);
  return;

 amomin_w:
  execAmomin_w(di);
  return;

 amomax_w:
  execAmomax_w(di);
  return;

 amominu_w:
  execAmominu_w(di);
  return;

 amomaxu_w:
  execAmomaxu_w(di);
  return;

 lr_d:
  execLr_d(di);
  return;

 sc_d:
  execSc_d(di);
  return;

 amoswap_d:
  execAmoswap_d(di);
  return;

 amoadd_d:
  execAmoadd_d(di);
  return;

 amoxor_d:
  execAmoxor_d(di);
  return;

 amoand_d:
  execAmoand_d(di);
  return;

 amoor_d:
  execAmoor_d(di);
  return;

 amomin_d:
  execAmomin_d(di);
  return;

 amomax_d:
  execAmomax_d(di);
  return;

 amominu_d:
  execAmominu_d(di);
  return;

 amomaxu_d:
  execAmomaxu_d(di);
  return;

 flw:
  execFlw(di);
  return;

 fsw:
  execFsw(di);
  return;

 fmadd_s:
  execFmadd_s(di);
  return;

 fmsub_s:
  execFmsub_s(di);
  return;

 fnmsub_s:
  execFnmsub_s(di);
  return;

 fnmadd_s:
  execFnmadd_s(di);
  return;

 fadd_s:
  execFadd_s(di);
  return;

 fsub_s:
  execFsub_s(di);
  return;

 fmul_s:
  execFmul_s(di);
  return;

 fdiv_s:
  execFdiv_s(di);
  return;

 fsqrt_s:
  execFsqrt_s(di);
  return;

 fsgnj_s:
  execFsgnj_s(di);
  return;

 fsgnjn_s:
  execFsgnjn_s(di);
  return;

 fsgnjx_s:
  execFsgnjx_s(di);
  return;

 fmin_s:
  execFmin_s(di);
  return;

 fmax_s:
  execFmax_s(di);
  return;

 fcvt_w_s:
  execFcvt_w_s(di);
  return;

 fcvt_wu_s:
  execFcvt_wu_s(di);
  return;

 fmv_x_w:
  execFmv_x_w(di);
  return;

 feq_s:
  execFeq_s(di);
  return;

 flt_s:
  execFlt_s(di);
  return;

 fle_s:
  execFle_s(di);
  return;

 fclass_s:
  execFclass_s(di);
  return;

 fcvt_s_w:
  execFcvt_s_w(di);
  return;

 fcvt_s_wu:
  execFcvt_s_wu(di);
  return;

 fmv_w_x:
  execFmv_w_x(di);
  return;

 fcvt_l_s:
  execFcvt_l_s(di);
  return;

 fcvt_lu_s:
  execFcvt_lu_s(di);
  return;

 fcvt_s_l:
  execFcvt_s_l(di);
  return;

 fcvt_s_lu:
  execFcvt_s_lu(di);
  return;

 fld:
  execFld(di);
  return;

 fsd:
  execFsd(di);
  return;

 fmadd_d:
  execFmadd_d(di);
  return;

 fmsub_d:
  execFmsub_d(di);
  return;

 fnmsub_d:
  execFnmsub_d(di);
  return;

 fnmadd_d:
  execFnmadd_d(di);
  return;

 fadd_d:
  execFadd_d(di);
  return;

 fsub_d:
  execFsub_d(di);
  return;

 fmul_d:
  execFmul_d(di);
  return;

 fdiv_d:
  execFdiv_d(di);
  return;

 fsqrt_d:
  execFsqrt_d(di);
  return;

 fsgnj_d:
  execFsgnj_d(di);
  return;

 fsgnjn_d:
  execFsgnjn_d(di);
  return;

 fsgnjx_d:
  execFsgnjx_d(di);
  return;

 fmin_d:
  execFmin_d(di);
  return;

 fmax_d:
  execFmax_d(di);
  return;

 fcvt_s_d:
  execFcvt_s_d(di);
  return;

 fcvt_d_s:
  execFcvt_d_s(di);
  return;

 feq_d:
  execFeq_d(di);
  return;

 flt_d:
  execFlt_d(di);
  return;

 fle_d:
  execFle_d(di);
  return;

 fclass_d:
  execFclass_d(di);
  return;

 fcvt_w_d:
  execFcvt_w_d(di);
  return;

 fcvt_wu_d:
  execFcvt_wu_d(di);
  return;

 fcvt_d_w:
  execFcvt_d_w(di);
  return;

 fcvt_d_wu:
  execFcvt_d_wu(di);
  return;

 fcvt_l_d:
  execFcvt_l_d(di);
  return;

 fcvt_lu_d:
  execFcvt_lu_d(di);
  return;

 fmv_x_d:
  execFmv_x_d(di);
  return;

 fcvt_d_l:
  execFcvt_d_l(di);
  return;

 fcvt_d_lu:
  execFcvt_d_lu(di);
  return;

 fmv_d_x:
  execFmv_d_x(di);
  return;

 mret:
  execMret(di);
  return;

 uret:
  execUret(di);
  return;

 sret:
  execSret(di);
  return;

 wfi:
  execWfi(di);
  return;

 sfence_vma:
  execSfence_vma(di);
  return;

 c_addi4spn:
  execAddi(di);
  return;

 c_fld:
  execFld(di);
  return;

 c_lq:
  illegalInst(di);
  return;

 c_lw:
  execLw(di);
  return;

 c_flw:
  execFlw(di);
  return;

 c_ld:
  execLd(di);
  return;

 c_fsd:
  execFsd(di);
  return;

 c_sq:
  illegalInst(di);
  return;

 c_sw:
  execSw(di);
  return;

 c_fsw:
  execFsw(di);
  return;

 c_sd:
  execSd(di);
  return;

 c_addi:
  execAddi(di);
  return;

 c_jal:
  execJal(di);
  return;

 c_li:
  // execAddi(di);
  intRegs_.write(di->op0(), di->op2As<SRV>());
  return;

 c_addi16sp:
  execAddi(di);
  return;

 c_lui:
  execLui(di);
  return;

 c_srli:
  execSrli(di);
  return;

 c_srli64:
  execSrli(di);
  return;

 c_srai:
  execSrai(di);
  return;

 c_srai64:
  execSrai(di);
  return;

 c_andi:
  execAndi(di);
  return;

 c_sub:
  execSub(di);
  return;

 c_xor:
  execXor(di);
  return;

 c_or:
  execOr(di);
  return;

 c_and:
  execAnd(di);
  return;

 c_subw:
  execSubw(di);
  return;

 c_addw:
  execAddw(di);
  return;

 c_j:
  execJal(di);
  return;

 c_beqz:
  execBeq(di);
  return;

 c_bnez:
  execBne(di);
  return;

 c_slli:
  execSlli(di);
  return;

 c_slli64:
  execSlli(di);
  return;

 c_fldsp:
  execFld(di);
  return;

 c_lwsp:
  execLw(di);
  return;

 c_flwsp:
  execFlw(di);
  return;

 c_ldsp:
  execLd(di);
  return;

 c_jr:
  execJalr(di);
  return;

 c_mv:
  intRegs_.write(di->op0(), intRegs_.read(di->op2()));
  //execAdd(di);
  return;

 c_ebreak:
  execEbreak(di);
  return;

 c_jalr:
  execJalr(di);
  return;

 c_add:
  execAdd(di);
  return;

 c_fsdsp:
  execFsd(di);
  return;

 c_swsp:
  execSw(di);
  return;

 c_fswsp:
  execFsw(di);
  return;

 c_addiw:
  execAddiw(di);
  return;

 c_sdsp:
  execSd(di);
  return;

 clz:
  execClz(di);
  return;

 ctz:
  execCtz(di);
  return;

 pcnt:
  execPcnt(di);
  return;

 andn:
  execAndn(di);
  return;

 orn:
  execOrn(di);
  return;

 xnor:
  execXnor(di);
  return;

 slo:
  execSlo(di);
  return;

 sro:
  execSro(di);
  return;

 sloi:
  execSloi(di);
  return;

 sroi:
  execSroi(di);
  return;

 min:
  execMin(di);
  return;

 max:
  execMax(di);
  return;

 minu:
  execMinu(di);
  return;

 maxu:
  execMaxu(di);
  return;

 rol:
  execRol(di);
  return;

 ror:
  execRor(di);
  return;

 rori:
  execRori(di);
  return;

 rev8:
  execRev8(di);
  return;

 rev:
  execRev(di);
  return;

 pack:
  execPack(di);
  return;

 addwu:
  execAddwu(di);
  return;

 subwu:
  execSubwu(di);
  return;

 addiwu:
  execAddiwu(di);
  return;

 sext_b:
  execSext_b(di);
  return;

 sext_h:
  execSext_h(di);
  return;

 addu_w:
  execAddu_w(di);
  return;

 subu_w:
  execSubu_w(di);
  return;

 slliu_w:
  execSlliu_w(di);
  return;

 packh:
  execPackh(di);
  return;

 packu:
  execPacku(di);
  return;

 packw:
  execPackw(di);
  return;

 packuw:
  execPackuw(di);
  return;

 grev:
  execGrev(di);
  return;

 grevi:
  execGrevi(di);
  return;

 gorc:
  execGorc(di);
  return;

 gorci:
  execGorci(di);
  return;

 shfl:
  execShfl(di);
  return;

 shfli:
  execShfli(di);
  return;

 unshfl:
  execUnshfl(di);
  return;

 unshfli:
  execUnshfli(di);
  return;

 sbset:
  execSbset(di);
  return;

 sbclr:
  execSbclr(di);
  return;

 sbinv:
  execSbinv(di);
  return;

 sbext:
  execSbext(di);
  return;

 sbseti:
  execSbseti(di);
  return;

 sbclri:
  execSbclri(di);
  return;

 sbinvi:
  execSbinvi(di);
  return;

 sbexti:
  execSbexti(di);
  return;

 bdep:
  execBdep(di);
  return;

 bext:
  execBext(di);
  return;

 bfp:
  execBfp(di);
  return;

 clmul:
  execClmul(di);
  return;

 clmulh:
  execClmulh(di);
  return;

 clmulr:
  execClmulr(di);
  return;

 sh1add:
  execSh1add(di);
  return;

 sh2add:
  execSh2add(di);
  return;

 sh3add:
  execSh3add(di);
  return;

 sh1addu_w:
  execSh1addu_w(di);
  return;

 sh2addu_w:
  execSh2addu_w(di);
  return;

 sh3addu_w:
  execSh3addu_w(di);
  return;

 crc32_b:
  execCrc32_b(di);
  return;

 crc32_h:
  execCrc32_h(di);
  return;

 crc32_w:
  execCrc32_w(di);
  return;

 crc32_d:
  execCrc32_d(di);
  return;

 crc32c_b:
  execCrc32c_b(di);
  return;

 crc32c_h:
  execCrc32c_h(di);
  return;

 crc32c_w:
  execCrc32c_w(di);
  return;

 crc32c_d:
  execCrc32c_d(di);
  return;

 bmator:
  execBmator(di);
  return;

 bmatxor:
  execBmatxor(di);
  return;

 bmatflip:
  execBmatflip(di);
  return;

 cmov:
  execCmov(di);
  return;

 cmix:
  execCmix(di);
  return;

 fsl:
  execFsl(di);
  return;

 fsr:
  execFsr(di);
  return;

 fsri:
  execFsri(di);
  return;
}


template <typename URV>
void
Hart<URV>::enableInstructionFrequency(bool b)
{
  instFreq_ = b;
  if (b)
    {
      instProfileVec_.resize(size_t(InstId::maxId) + 1);

      auto regCount = intRegCount();
      for (auto& inst : instProfileVec_)
	{
	  inst.rd_.resize(regCount);
	  inst.rs1_.resize(regCount);
	  inst.rs2_.resize(regCount);
	  inst.rs1Histo_.resize(13);  // FIX: avoid magic 13
	  inst.rs2Histo_.resize(13);  // FIX: avoid magic 13
	  inst.immHisto_.resize(13);  // FIX: avoid magic 13
	}
    }
}


template <typename URV>
void
Hart<URV>::enterDebugMode(DebugModeCause cause, URV pc)
{
  // Entering debug modes loses LR reservation.
  memory_.invalidateLr(hartIx_);

  if (debugMode_)
    {
      if (debugStepMode_)
	debugStepMode_ = false;
      else
	std::cerr << "Error: Entering debug-halt while in debug-halt\n";
    }
  else
    {
      debugMode_ = true;
      if (debugStepMode_)
	std::cerr << "Error: Entering debug-halt with debug-step true\n";
      debugStepMode_ = false;
    }

  URV value = 0;
  if (peekCsr(CsrNumber::DCSR, value))
    {
      value &= ~(URV(7) << 6);        // Clear cause field (starts at bit 6).
      value |= URV(cause) << 6;       // Set cause field
      value = (value >> 2) << 2;      // Clear privilege mode bits.
      value |= URV(privMode_) & 0x3;  // Set privelge mode bits.

      if (nmiPending_)
	value |= URV(1) << 3;    // Set nmip bit.
      csRegs_.poke(CsrNumber::DCSR, value);

      csRegs_.poke(CsrNumber::DPC, pc);
    }
}


template <typename URV>
void
Hart<URV>::enterDebugMode(URV pc)
{
  if (forceAccessFail_)
    {
      std::cerr << "Entering debug mode with a pending forced exception from"
		<< " test-bench. Exception cleared.\n";
      forceAccessFail_ = false;
    }

  // This method is used by the test-bench to make the simulator
  // follow it into debug-halt or debug-stop mode. Do nothing if the
  // simulator got into debug mode on its own.
  if (debugMode_)
    return;   // Already in debug mode.

  if (debugStepMode_)
    std::cerr << "Error: Enter-debug command finds hart in debug-step mode.\n";

  debugStepMode_ = false;
  debugMode_ = false;

  enterDebugMode(DebugModeCause::DEBUGGER, pc);
}


template <typename URV>
void
Hart<URV>::exitDebugMode()
{
  if (not debugMode_)
    {
      std::cerr << "Error: Bench sent exit debug while not in debug mode.\n";
      return;
    }

  peekCsr(CsrNumber::DPC, pc_);

  // If in debug-step go to debug-halt. If in debug-halt go to normal
  // or debug-step based on step-bit in DCSR.
  if (debugStepMode_)
    debugStepMode_ = false;
  else
    {
      if (dcsrStep_)
	debugStepMode_ = true;
      else
        {
          debugMode_ = false;
        }
    }

  // If pending nmi bit is set in dcsr, set pending nmi in the hart
  // object.
  URV dcsrVal = 0;
  if (not peekCsr(CsrNumber::DCSR, dcsrVal))
    std::cerr << "Error: Failed to read DCSR in exit debug.\n";

  if ((dcsrVal >> 3) & 1)
    setPendingNmi(nmiCause_);
}


template <typename URV>
void
Hart<URV>::execBlt(const DecodedInst* di)
{
  SRV v1 = intRegs_.read(di->op0()),  v2 = intRegs_.read(di->op1());
  if (v1 < v2)
    {
      pc_ = currPc_ + di->op2As<SRV>();
      pc_ = (pc_ >> 1) << 1;  // Clear least sig bit.
      lastBranchTaken_ = true;
    }
}


template <typename URV>
void
Hart<URV>::execBltu(const DecodedInst* di)
{
  URV v1 = intRegs_.read(di->op0()),  v2 = intRegs_.read(di->op1());
  if (v1 < v2)
    {
      pc_ = currPc_ + di->op2As<SRV>();
      pc_ = (pc_ >> 1) << 1;  // Clear least sig bit.
      lastBranchTaken_ = true;
    }
}


template <typename URV>
void
Hart<URV>::execBge(const DecodedInst* di)
{
  SRV v1 = intRegs_.read(di->op0()),  v2 = intRegs_.read(di->op1());
  if (v1 >= v2)
    {
      pc_ = currPc_ + di->op2As<SRV>();
      pc_ = (pc_ >> 1) << 1;  // Clear least sig bit.
      lastBranchTaken_ = true;
    }
}


template <typename URV>
void
Hart<URV>::execBgeu(const DecodedInst* di)
{
  URV v1 = intRegs_.read(di->op0()),  v2 = intRegs_.read(di->op1());
  if (v1 >= v2)
    {
      pc_ = currPc_ + di->op2As<SRV>();
      pc_ = (pc_ >> 1) << 1;  // Clear least sig bit.
      lastBranchTaken_ = true;
    }
}


template <typename URV>
void
Hart<URV>::execJalr(const DecodedInst* di)
{
  URV temp = pc_;  // pc has the address of the instruction after jalr
  pc_ = intRegs_.read(di->op1()) + di->op2As<SRV>();
  pc_ = (pc_ >> 1) << 1;  // Clear least sig bit.
  intRegs_.write(di->op0(), temp);
  lastBranchTaken_ = true;
}


template <typename URV>
void
Hart<URV>::execJal(const DecodedInst* di)
{
  intRegs_.write(di->op0(), pc_);
  pc_ = currPc_ + SRV(int32_t(di->op1()));
  pc_ = (pc_ >> 1) << 1;  // Clear least sig bit.
  lastBranchTaken_ = true;
}


template <typename URV>
void
Hart<URV>::execAuipc(const DecodedInst* di)
{
  intRegs_.write(di->op0(), currPc_ + SRV(int32_t(di->op1())));
}


template <typename URV>
inline
bool
Hart<URV>::checkShiftImmediate(const DecodedInst* di, URV imm)
{
  bool bad = isRv64()? imm > 63 : imm > 31;

  if (bad)
    {
      illegalInst(di);
      return false;
    }
  return true;
}


template <typename URV>
void
Hart<URV>::execSlli(const DecodedInst* di)
{
  URV amount = di->op2();
  if (not checkShiftImmediate(di, amount))
    return;

  URV v = intRegs_.read(di->op1()) << amount;
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execSlti(const DecodedInst* di)
{
  SRV imm = di->op2As<SRV>();
  URV v = SRV(intRegs_.read(di->op1())) < imm ? 1 : 0;
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execSltiu(const DecodedInst* di)
{
  URV imm = di->op2As<SRV>();   // We sign extend then use as unsigned.
  URV v = intRegs_.read(di->op1()) < imm ? 1 : 0;
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execXori(const DecodedInst* di)
{
  URV v = intRegs_.read(di->op1()) ^ di->op2As<SRV>();
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execSrli(const DecodedInst* di)
{
  URV amount(di->op2());
  if (not checkShiftImmediate(di, amount))
    return;

  URV v = intRegs_.read(di->op1()) >> amount;
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execSrai(const DecodedInst* di)
{
  uint32_t amount(di->op2());
  if (not checkShiftImmediate(di, amount))
    return;

  URV v = SRV(intRegs_.read(di->op1())) >> amount;
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execOri(const DecodedInst* di)
{
  URV v = intRegs_.read(di->op1()) | di->op2As<SRV>();
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execSub(const DecodedInst* di)
{
  URV v = intRegs_.read(di->op1()) - intRegs_.read(di->op2());
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execSll(const DecodedInst* di)
{
  URV mask = intRegs_.shiftMask();
  URV v = intRegs_.read(di->op1()) << (intRegs_.read(di->op2()) & mask);
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execSlt(const DecodedInst* di)
{
  SRV v1 = intRegs_.read(di->op1());
  SRV v2 = intRegs_.read(di->op2());
  URV v = v1 < v2 ? 1 : 0;
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execSltu(const DecodedInst* di)
{
  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV v = v1 < v2 ? 1 : 0;
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execXor(const DecodedInst* di)
{
  URV v = intRegs_.read(di->op1()) ^ intRegs_.read(di->op2());
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execSrl(const DecodedInst* di)
{
  URV mask = intRegs_.shiftMask();
  URV v = intRegs_.read(di->op1()) >> (intRegs_.read(di->op2()) & mask);
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execSra(const DecodedInst* di)
{
  URV mask = intRegs_.shiftMask();
  URV v = SRV(intRegs_.read(di->op1())) >> (intRegs_.read(di->op2()) & mask);
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execOr(const DecodedInst* di)
{
  URV v = intRegs_.read(di->op1()) | intRegs_.read(di->op2());
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execAnd(const DecodedInst* di)
{
  URV v = intRegs_.read(di->op1()) & intRegs_.read(di->op2());
  intRegs_.write(di->op0(), v);
}


template <typename URV>
void
Hart<URV>::execFence(const DecodedInst*)
{
  loadQueue_.clear();
}


template <typename URV>
void
Hart<URV>::execFencei(const DecodedInst*)
{
  // invalidateDecodeCache();  // No need for this. We invalidate on each write.
}


template <typename URV>
ExceptionCause
Hart<URV>::validateAmoAddr(uint32_t rs1, uint64_t& addr, unsigned accessSize,
                           SecondaryCause& secCause)
{
  URV mask = URV(accessSize) - 1;

  auto cause = ExceptionCause::NONE;
  if (accessSize == 4)
    {
      uint32_t dummy = 0;
      cause = determineStoreException(rs1, addr, addr, dummy, secCause);
    }
  else
    {
      uint64_t dummy = 0;
      cause = determineStoreException(rs1, addr, addr, dummy, secCause);
    }

  // Address must be word aligned for word access and double-word
  // aligned for double-word access.
  bool fail = (addr & mask) != 0;

  // Check if invalid outside DCCM.
  if (amoInDccmOnly_ and not isAddrInDccm(addr))
    fail = true;

  if (fail)
    {
      // AMO secondary cause has priority over ECC.
      if (cause == ExceptionCause::NONE or
          secCause == SecondaryCause::STORE_ACC_DOUBLE_ECC)
        {
          // Per spec cause is store-access-fault.
          cause = ExceptionCause::STORE_ACC_FAULT;
          secCause = SecondaryCause::STORE_ACC_AMO;
        }
    }

  return cause;
}


template <typename URV>
bool
Hart<URV>::amoLoad32(uint32_t rs1, URV& value)
{
  enableWideLdStMode(false);

  URV virtAddr = intRegs_.read(rs1);

  ldStAddr_ = virtAddr;   // For reporting load addr in trace-mode.
  ldStAddrValid_ = true;  // For reporting load addr in trace-mode.

  if (loadQueueEnabled_)
    removeFromLoadQueue(rs1, false);

  if (hasActiveTrigger())
    {
      if (ldStAddrTriggerHit(virtAddr, TriggerTiming::Before, false /*isLoad*/,
			     privMode_, isInterruptEnabled()))
	triggerTripped_ = true;
    }

  unsigned ldSize = 4;

  auto secCause = SecondaryCause::STORE_ACC_AMO;
  
  uint64_t addr = virtAddr;
  auto cause = validateAmoAddr(rs1, addr, ldSize, secCause);

  if (cause != ExceptionCause::NONE)
    {
      if (not triggerTripped_)
        initiateLoadException(cause, virtAddr, secCause);
      return false;
    }

  uint32_t uval = 0;
  if (memory_.read(addr, uval))
    {
      value = SRV(int32_t(uval)); // Sign extend.
      return true;  // Success.
    }

  cause = ExceptionCause::STORE_ACC_FAULT;
  initiateLoadException(cause, virtAddr, secCause);
  return false;
}


template <typename URV>
bool
Hart<URV>::amoLoad64(uint32_t rs1, URV& value)
{
  URV virtAddr = intRegs_.read(rs1);

  ldStAddr_ = virtAddr;   // For reporting load addr in trace-mode.
  ldStAddrValid_ = true;  // For reporting load addr in trace-mode.

  if (loadQueueEnabled_)
    removeFromLoadQueue(rs1, false);

  if (hasActiveTrigger())
    {
      if (ldStAddrTriggerHit(virtAddr, TriggerTiming::Before, false /*isLoad*/,
			     privMode_, isInterruptEnabled()))
	triggerTripped_ = true;
    }

  unsigned ldSize = 8;

  auto secCause = SecondaryCause::STORE_ACC_AMO;
  uint64_t addr = virtAddr;
  auto cause = validateAmoAddr(rs1, addr, ldSize, secCause);

  if (cause != ExceptionCause::NONE)
    {
      if (not triggerTripped_)
        initiateLoadException(cause, virtAddr, secCause);
      return false;
    }

  uint64_t uval = 0;
  if (memory_.read(addr, uval))
    {
      value = SRV(int64_t(uval)); // Sign extend.
      return true;  // Success.
    }

  cause = ExceptionCause::STORE_ACC_FAULT;
  initiateLoadException(cause, virtAddr, secCause);
  return false;
}


template <typename URV>
void
Hart<URV>::execEcall(const DecodedInst*)
{
  if (triggerTripped_)
    return;

  if (newlib_ or linux_)
    {
      URV a0 = syscall_.emulate();
      intRegs_.write(RegA0, a0);
      return;
    }

  auto secCause = SecondaryCause::NONE;

  if (privMode_ == PrivilegeMode::Machine)
    initiateException(ExceptionCause::M_ENV_CALL, currPc_, 0, secCause);
  else if (privMode_ == PrivilegeMode::Supervisor)
    initiateException(ExceptionCause::S_ENV_CALL, currPc_, 0, secCause);
  else if (privMode_ == PrivilegeMode::User)
    initiateException(ExceptionCause::U_ENV_CALL, currPc_, 0, secCause);
  else
    assert(0 and "Invalid privilege mode in execEcall");
}


template <typename URV>
void
Hart<URV>::execEbreak(const DecodedInst*)
{
  if (triggerTripped_)
    return;

  if (enableGdb_)
    {
      pc_ = currPc_;
      handleExceptionForGdb(*this, gdbInputFd_);
      return;
    }

  // If in machine/supervisor/user mode and DCSR bit ebreakm/s/u is
  // set, then enter debug mode.
  URV dcsrVal = 0;
  if (peekCsr(CsrNumber::DCSR, dcsrVal))
    {
      bool ebm = (dcsrVal >> 15) & 1;
      bool ebs = (dcsrVal >> 13) & 1;
      bool ebu = (dcsrVal >> 12) & 1;

      bool debug = ( (ebm and privMode_ == PrivilegeMode::Machine) or
                     (ebs and privMode_ == PrivilegeMode::Supervisor) or
                     (ebu and privMode_ == PrivilegeMode::User) );

      if (debug)
        {
          // The documentation (RISCV external debug support) does
          // not say whether or not we set EPC and MTVAL.
          enterDebugMode(DebugModeCause::EBREAK, currPc_);
          ebreakInstDebug_ = true;
          recordCsrWrite(CsrNumber::DCSR);
          return;
        }
    }

  URV savedPc = currPc_;  // Goes into MEPC.
  URV trapInfo = currPc_;  // Goes into MTVAL.

  auto cause = ExceptionCause::BREAKP;
  auto secCause = SecondaryCause::BREAKP;
  initiateException(cause, savedPc, trapInfo, secCause);
}


template <typename URV>
void
Hart<URV>::execSfence_vma(const DecodedInst* di)
{
  if (not isRvs())
    {
      illegalInst(di);
      return;
    }

  URV status = csRegs_.peekMstatus();

  MstatusFields<URV> fields(status);
  if (fields.bits_.TVM)
    {
      illegalInst(di);
      return;
    }

  // Invalidate whole TLB. This is overkill. TBD FIX: Improve.
  virtMem_.tlb_.invalidate();

  // std::cerr << "sfence.vma " << di->op1() << ' ' << di->op2() << '\n';
  if (di->op1() == 0)
    invalidateDecodeCache();
  else
    {
      uint64_t va = intRegs_.read(di->op1());
      uint64_t pageStart = virtMem_.pageStartAddress(va);
      uint64_t last = pageStart + virtMem_.pageSize();
      for (uint64_t addr = pageStart; addr < last; addr += 4)
        invalidateDecodeCache(addr, 4);
    }
}


template <typename URV>
void
Hart<URV>::execMret(const DecodedInst* di)
{
  if (privMode_ < PrivilegeMode::Machine)
    {
      illegalInst(di);
      return;
    }

  if (triggerTripped_)
    return;

  memory_.invalidateLr(hartIx_); // Clear LR reservation (if any).

  // Restore privilege mode and interrupt enable by getting
  // current value of MSTATUS, ...
  URV value = csRegs_.peekMstatus();

  // ... updating/unpacking its fields,
  MstatusFields<URV> fields(value);
  PrivilegeMode savedMode = PrivilegeMode(fields.bits_.MPP);
  fields.bits_.MIE = fields.bits_.MPIE;
  fields.bits_.MPP = unsigned(PrivilegeMode::User);
  fields.bits_.MPIE = 1;
  if (savedMode != PrivilegeMode::Machine)
    fields.bits_.MPRV = 0;

  // ... and putting it back
  if (not csRegs_.write(CsrNumber::MSTATUS, privMode_, fields.value_))
    assert(0 and "Failed to write MSTATUS register\n");
  updateCachedMstatusFields();

  // Restore program counter from MEPC.
  URV epc;
  if (not csRegs_.read(CsrNumber::MEPC, privMode_, epc))
    illegalInst(di);
  pc_ = (epc >> 1) << 1;  // Restore pc clearing least sig bit.
      
  // Update privilege mode.
  privMode_ = savedMode;
}


template <typename URV>
void
Hart<URV>::execSret(const DecodedInst* di)
{
  if (not isRvs())
    {
      illegalInst(di);
      return;
    }

  if (privMode_ < PrivilegeMode::Supervisor)
    {
      illegalInst(di);
      return;
    }

  if (triggerTripped_)
    return;

  // Restore privilege mode and interrupt enable by getting
  // current value of MSTATUS, ...
  URV value = 0;
  if (not csRegs_.read(CsrNumber::SSTATUS, privMode_, value))
    {
      illegalInst(di);
      return;
    }

  memory_.invalidateLr(hartIx_); // Clear LR reservation (if any).

  // ... updating/unpacking its fields,
  MstatusFields<URV> fields(value);
  PrivilegeMode savedMode = fields.bits_.SPP? PrivilegeMode::Supervisor :
    PrivilegeMode::User;
  fields.bits_.SIE = fields.bits_.SPIE;
  fields.bits_.SPP = 0;
  fields.bits_.SPIE = 1;
  if (savedMode != PrivilegeMode::Machine)
    fields.bits_.MPRV = 0;

  // ... and putting it back
  if (not csRegs_.write(CsrNumber::SSTATUS, privMode_, fields.value_))
    {
      illegalInst(di);
      return;
    }
  updateCachedMstatusFields();

  // Restore program counter from SEPC.
  URV epc;
  if (not csRegs_.read(CsrNumber::SEPC, privMode_, epc))
    {
      illegalInst(di);
      return;
    }
  pc_ = (epc >> 1) << 1;  // Restore pc clearing least sig bit.

  // Update privilege mode.
  privMode_ = savedMode;
}


template <typename URV>
void
Hart<URV>::execUret(const DecodedInst* di)
{
  if (not isRvu())
    {
      illegalInst(di);
      return;
    }

  if (privMode_ != PrivilegeMode::User)
    {
      illegalInst(di);
      return;
    }

  if (triggerTripped_)
    return;

  // Restore privilege mode and interrupt enable by getting
  // current value of MSTATUS, ...
  URV value = 0;
  if (not csRegs_.read(CsrNumber::USTATUS, privMode_, value))
    {
      illegalInst(di);
      return;
    }

  // ... updating/unpacking its fields,
  MstatusFields<URV> fields(value);
  fields.bits_.UIE = fields.bits_.UPIE;
  fields.bits_.UPIE = 1;

  // ... and putting it back
  if (not csRegs_.write(CsrNumber::USTATUS, privMode_, fields.value_))
    {
      illegalInst(di);
      return;
    }
  updateCachedMstatusFields();

  // Restore program counter from UEPC.
  URV epc;
  if (not csRegs_.read(CsrNumber::UEPC, privMode_, epc))
    {
      illegalInst(di);
      return;
    }
  pc_ = (epc >> 1) << 1;  // Restore pc clearing least sig bit.
}


template <typename URV>
void
Hart<URV>::execWfi(const DecodedInst*)
{
  return;   // Currently implemented as a no-op.
}


template <typename URV>
bool
Hart<URV>::doCsrRead(const DecodedInst* di, CsrNumber csr, URV& value)
{
  if (csRegs_.read(csr, privMode_, value))
    return true;

  illegalInst(di);
  return false;
}


template <typename URV>
void
Hart<URV>::updateStackChecker()
{  
  peekCsr(CsrNumber::MSPCBA, stackMax_);
  peekCsr(CsrNumber::MSPCTA, stackMin_);

  URV val = 0;
  checkStackAccess_ = peekCsr(CsrNumber::MSPCC, val) and val != 0;
}


template <typename URV>
void
Hart<URV>::doCsrWrite(const DecodedInst* di, CsrNumber csr, URV csrVal,
                      unsigned intReg, URV intRegVal)
{
  if (not csRegs_.isWriteable(csr, privMode_))
    {
      illegalInst(di);
      return;
    }

  // Make auto-increment happen before write for minstret and cycle.
  if (csr == CsrNumber::MINSTRET or csr == CsrNumber::MINSTRETH)
    if (minstretEnabled())
      retiredInsts_++;
  if (csr == CsrNumber::MCYCLE or csr == CsrNumber::MCYCLEH)
    cycleCount_++;

  // Update CSR and integer register.
  csRegs_.write(csr, privMode_, csrVal);
  intRegs_.write(intReg, intRegVal);

  if (csr == CsrNumber::DCSR)
    {
      dcsrStep_ = (csrVal >> 2) & 1;
      dcsrStepIe_ = (csrVal >> 11) & 1;
    }
  else if (csr >= CsrNumber::MSPCBA and csr <= CsrNumber::MSPCC)
    updateStackChecker();
  else if (csr == CsrNumber::MDBAC)
    enableWideLdStMode(true);
  else if ((csr >= CsrNumber::PMPADDR0 and csr <= CsrNumber::PMPADDR15) or
           (csr >= CsrNumber::PMPCFG0 and csr <= CsrNumber::PMPCFG3))
    updateMemoryProtection();
  else if (csr == CsrNumber::SATP)
    updateAddressTranslation();
  else if (csr == CsrNumber::FCSR or csr == CsrNumber::FRM or csr == CsrNumber::FFLAGS)
    markFsDirty(); // Update FS field of MSTATS if FCSR is written

  // Update cached values of MSTATUS MPP and MPRV.
  if (csr == CsrNumber::MSTATUS or csr == CsrNumber::SSTATUS)
    updateCachedMstatusFields();

  // Csr was written. If it was minstret, compensate for
  // auto-increment that will be done by run, runUntilAddress or
  // singleStep method.
  if (csr == CsrNumber::MINSTRET or csr == CsrNumber::MINSTRETH)
    if (minstretEnabled())
      retiredInsts_--;

  // Same for mcycle.
  if (csr == CsrNumber::MCYCLE or csr == CsrNumber::MCYCLEH)
    cycleCount_--;
}


// Set control and status register csr (op2) to value of register rs1
// (op1) and save its original value in register rd (op0).
template <typename URV>
void
Hart<URV>::execCsrrw(const DecodedInst* di)
{
  if (triggerTripped_)
    return;

  CsrNumber csr = CsrNumber(di->op2());

  URV prev = 0;
  if (not doCsrRead(di, csr, prev))
    return;

  URV next = intRegs_.read(di->op1());

  doCsrWrite(di, csr, next, di->op0(), prev);
}


template <typename URV>
void
Hart<URV>::execCsrrs(const DecodedInst* di)
{
  if (triggerTripped_)
    return;

  CsrNumber csr = CsrNumber(di->op2());

  URV prev = 0;
  if (not doCsrRead(di, csr, prev))
    return;

  URV next = prev | intRegs_.read(di->op1());
  if (di->op1() == 0)
    {
      intRegs_.write(di->op0(), prev);
      return;
    }

  doCsrWrite(di, csr, next, di->op0(), prev);
}


template <typename URV>
void
Hart<URV>::execCsrrc(const DecodedInst* di)
{
  if (triggerTripped_)
    return;

  CsrNumber csr = CsrNumber(di->op2());

  URV prev = 0;
  if (not doCsrRead(di, csr, prev))
    return;

  URV next = prev & (~ intRegs_.read(di->op1()));
  if (di->op1() == 0)
    {
      intRegs_.write(di->op0(), prev);
      return;
    }

  doCsrWrite(di, csr, next, di->op0(), prev);
}


template <typename URV>
void
Hart<URV>::execCsrrwi(const DecodedInst* di)
{
  if (triggerTripped_)
    return;

  CsrNumber csr = CsrNumber(di->op2());

  URV prev = 0;
  if (di->op0() != 0)
    if (not doCsrRead(di, csr, prev))
      return;

  doCsrWrite(di, csr, di->op1(), di->op0(), prev);
}


template <typename URV>
void
Hart<URV>::execCsrrsi(const DecodedInst* di)
{
  if (triggerTripped_)
    return;

  CsrNumber csr = CsrNumber(di->op2());

  URV prev = 0;
  if (not doCsrRead(di, csr, prev))
    return;

  uint32_t imm = di->op1();

  URV next = prev | imm;
  if (imm == 0)
    {
      intRegs_.write(di->op0(), prev);
      return;
    }

  doCsrWrite(di, csr, next, di->op0(), prev);
}


template <typename URV>
void
Hart<URV>::execCsrrci(const DecodedInst* di)
{
  if (triggerTripped_)
    return;

  CsrNumber csr = CsrNumber(di->op2());

  URV prev = 0;
  if (not doCsrRead(di, csr, prev))
    return;

  uint32_t imm = di->op1();

  URV next = prev & (~ imm);
  if (imm == 0)
    {
      intRegs_.write(di->op0(), prev);
      return;
    }

  doCsrWrite(di, csr, next, di->op0(), prev);
}


template <typename URV>
void
Hart<URV>::execLb(const DecodedInst* di)
{
  load<int8_t>(di->op0(), di->op1(), di->op2As<int32_t>());
}


template <typename URV>
void
Hart<URV>::execLbu(const DecodedInst* di)
{
  load<uint8_t>(di->op0(), di->op1(), di->op2As<int32_t>());
}


template <typename URV>
void
Hart<URV>::execLhu(const DecodedInst* di)
{
  load<uint16_t>(di->op0(), di->op1(), di->op2As<int32_t>());
}


template <typename URV>
bool
Hart<URV>::wideStore(URV addr, URV storeVal, unsigned storeSize)
{
  if ((addr & 7) or storeSize != 4 or not isDataAddressExternal(addr))
    {
      auto cause = ExceptionCause::STORE_ACC_FAULT;
      auto secCause = SecondaryCause::STORE_ACC_64BIT;
      initiateStoreException(cause, addr, secCause);
      return false;
    }

  uint32_t lower = storeVal;

  URV temp = 0;
  peekCsr(CsrNumber::MDBHD, temp);
  uint32_t upper = temp;

  if (not memory_.write(hartIx_, addr + 4, upper) or
      not memory_.write(hartIx_, addr, lower))
    {
      // FIX: Clear last written data if 1st 4 bytes written.
      auto cause = ExceptionCause::STORE_ACC_FAULT;
      auto secCause = SecondaryCause::STORE_ACC_64BIT;
      initiateStoreException(cause, addr, secCause);
      return false;
    }

  return true;
}


template <typename URV>
template <typename STORE_TYPE>
ExceptionCause
Hart<URV>::determineStoreException(unsigned rs1, URV base, uint64_t& addr,
				   STORE_TYPE& storeVal,
				   SecondaryCause& secCause)
{
  unsigned stSize = sizeof(STORE_TYPE);

  // Misaligned store to io section causes an exception. Crossing
  // dccm to non-dccm causes an exception.
  constexpr unsigned alignMask = sizeof(STORE_TYPE) - 1;
  bool misal = addr & alignMask;
  misalignedLdSt_ = misal;

  ExceptionCause cause = ExceptionCause::NONE;
  if (misal)
    {
      cause = determineMisalStoreException(addr, stSize, secCause);
      if (cause == ExceptionCause::STORE_ADDR_MISAL)
        return cause;
    }

  // Stack access.
  if (rs1 == RegSp and checkStackAccess_ and
      not checkStackStore(base, addr, stSize))
    {
      secCause = SecondaryCause::STORE_ACC_STACK_CHECK;
      return ExceptionCause::STORE_ACC_FAULT;
    }

  bool writeOk = false;

  // Address translation
  if (isRvs())
    {
      PrivilegeMode mode = mstatusMprv_? mstatusMpp_ : privMode_;
      if (mode != PrivilegeMode::Machine)
        {
          uint64_t pa = 0;
          cause = virtMem_.translateForStore(addr, mode, pa);
          if (cause != ExceptionCause::NONE)
            return cause;
          addr = pa;
        }
      writeOk = memory_.checkWrite(addr, storeVal);
    }
  else
    {
      // DCCM unmapped or out of MPU windows. Invalid PIC access handled later.
      writeOk = memory_.checkWrite(addr, storeVal);
      if (not writeOk and not isAddrMemMapped(addr))
        {
          secCause = SecondaryCause::STORE_ACC_MEM_PROTECTION;
          size_t region = memory_.getRegionIndex(addr);
          if (regionHasLocalDataMem_.at(region))
            secCause = SecondaryCause::STORE_ACC_LOCAL_UNMAPPED;
          return ExceptionCause::STORE_ACC_FAULT;
        }

      // 64-bit store
      if (wideLdSt_)
        {
          bool fail = (addr & 7) or stSize != 4 or ! isDataAddressExternal(addr);
          uint64_t val = 0;
          fail = fail or ! memory_.checkWrite(addr, val);
          if (fail)
            {
              secCause = SecondaryCause::STORE_ACC_64BIT;
              return ExceptionCause::STORE_ACC_FAULT;
            }
        }

      // Region predict (Effective address compatible with base).
      if (eaCompatWithBase_ and effectiveAndBaseAddrMismatch(addr, base))
        {
          secCause = SecondaryCause::STORE_ACC_REGION_PREDICTION;
          return ExceptionCause::STORE_ACC_FAULT;
        }
    }

  // PIC access
  if (isAddrMemMapped(addr))
    {
      if (privMode_ != PrivilegeMode::Machine)
        {
          secCause = SecondaryCause::STORE_ACC_LOCAL_UNMAPPED;
	  return ExceptionCause::STORE_ACC_FAULT;
        }
      if (not writeOk)
        {
          secCause = SecondaryCause::STORE_ACC_PIC;
          return ExceptionCause::STORE_ACC_FAULT;
        }
    }

  if (not writeOk)
    {
      secCause = SecondaryCause::STORE_ACC_MEM_PROTECTION;
      return ExceptionCause::STORE_ACC_FAULT;
    }

  if (misal and cause != ExceptionCause::NONE)
    return cause;

  // Physical memory protection.
  if (pmpEnabled_)
    {
      Pmp pmp = pmpManager_.accessPmp(addr);
      if (not pmp.isWrite(privMode_, mstatusMpp_, mstatusMprv_) and
          not isAddrMemMapped(addr))
        {
          secCause = SecondaryCause::STORE_ACC_PMP;
          return ExceptionCause::STORE_ACC_FAULT;
        }
    }

  // Fault dictated by test-bench
  if (forceAccessFail_)
    {
      secCause = forcedCause_;
      return ExceptionCause::STORE_ACC_FAULT;
    }

  return ExceptionCause::NONE;
}


template <typename URV>
void
Hart<URV>::execSb(const DecodedInst* di)
{
  uint32_t rs1 = di->op1();
  URV base = intRegs_.read(rs1);
  URV addr = base + di->op2As<SRV>();
  uint8_t value = uint8_t(intRegs_.read(di->op0()));
  store<uint8_t>(rs1, base, addr, value);
}


template <typename URV>
void
Hart<URV>::execSh(const DecodedInst* di)
{
  uint32_t rs1 = di->op1();
  URV base = intRegs_.read(rs1);
  URV addr = base + di->op2As<SRV>();
  uint16_t value = uint16_t(intRegs_.read(di->op0()));
  store<uint16_t>(rs1, base, addr, value);
}


namespace WdRiscv
{

  template<>
  void
  Hart<uint32_t>::execMul(const DecodedInst* di)
  {
    int32_t a = intRegs_.read(di->op1());
    int32_t b = intRegs_.read(di->op2());

    int32_t c = a * b;
    intRegs_.write(di->op0(), c);
  }


  template<>
  void
  Hart<uint32_t>::execMulh(const DecodedInst* di)
  {
    int64_t a = int32_t(intRegs_.read(di->op1()));  // sign extend.
    int64_t b = int32_t(intRegs_.read(di->op2()));
    int64_t c = a * b;
    int32_t high = static_cast<int32_t>(c >> 32);

    intRegs_.write(di->op0(), high);
  }


  template <>
  void
  Hart<uint32_t>::execMulhsu(const DecodedInst* di)
  {
    int64_t a = int32_t(intRegs_.read(di->op1()));
    uint64_t b = intRegs_.read(di->op2());
    int64_t c = a * b;
    int32_t high = static_cast<int32_t>(c >> 32);

    intRegs_.write(di->op0(), high);
  }


  template <>
  void
  Hart<uint32_t>::execMulhu(const DecodedInst* di)
  {
    uint64_t a = intRegs_.read(di->op1());
    uint64_t b = intRegs_.read(di->op2());
    uint64_t c = a * b;
    uint32_t high = static_cast<uint32_t>(c >> 32);

    intRegs_.write(di->op0(), high);
  }


  template<>
  void
  Hart<uint64_t>::execMul(const DecodedInst* di)
  {
    Int128 a = int64_t(intRegs_.read(di->op1()));  // sign extend to 64-bit
    Int128 b = int64_t(intRegs_.read(di->op2()));

    int64_t c = static_cast<int64_t>(a * b);
    intRegs_.write(di->op0(), c);
  }


  template<>
  void
  Hart<uint64_t>::execMulh(const DecodedInst* di)
  {
    Int128 a = int64_t(intRegs_.read(di->op1()));  // sign extend.
    Int128 b = int64_t(intRegs_.read(di->op2()));
    Int128 c = a * b;
    int64_t high = static_cast<int64_t>(c >> 64);

    intRegs_.write(di->op0(), high);
  }


  template <>
  void
  Hart<uint64_t>::execMulhsu(const DecodedInst* di)
  {

    Int128 a = int64_t(intRegs_.read(di->op1()));
    Int128 b = intRegs_.read(di->op2());
    Int128 c = a * b;
    int64_t high = static_cast<int64_t>(c >> 64);

    intRegs_.write(di->op0(), high);
  }


  template <>
  void
  Hart<uint64_t>::execMulhu(const DecodedInst* di)
  {
    Uint128 a = intRegs_.read(di->op1());
    Uint128 b = intRegs_.read(di->op2());
    Uint128 c = a * b;
    uint64_t high = static_cast<uint64_t>(c >> 64);

    intRegs_.write(di->op0(), high);
  }

}


template <typename URV>
void
Hart<URV>::execDiv(const DecodedInst* di)
{
  SRV a = intRegs_.read(di->op1());
  SRV b = intRegs_.read(di->op2());
  SRV c = -1;   // Divide by zero result
  if (b != 0)
    {
      SRV minInt = SRV(1) << (mxlen_ - 1);
      if (a == minInt and b == -1)
	c = a;
      else
	c = a / b;  // Per spec: User-Level ISA, Version 2.3, Section 6.2
    }

  recordDivInst(di->op0(), peekIntReg(di->op0()));

  intRegs_.write(di->op0(), c);
}


template <typename URV>
void
Hart<URV>::execDivu(const DecodedInst* di)
{
  URV a = intRegs_.read(di->op1());
  URV b = intRegs_.read(di->op2());
  URV c = ~ URV(0);  // Divide by zero result.
  if (b != 0)
    c = a / b;

  recordDivInst(di->op0(), peekIntReg(di->op0()));

  intRegs_.write(di->op0(), c);
}


// Remainder instruction.
template <typename URV>
void
Hart<URV>::execRem(const DecodedInst* di)
{
  SRV a = intRegs_.read(di->op1());
  SRV b = intRegs_.read(di->op2());
  SRV c = a;  // Divide by zero remainder.
  if (b != 0)
    {
      SRV minInt = SRV(1) << (mxlen_ - 1);
      if (a == minInt and b == -1)
	c = 0;   // Per spec: User-Level ISA, Version 2.3, Section 6.2
      else
	c = a % b;
    }

  recordDivInst(di->op0(), peekIntReg(di->op0()));

  intRegs_.write(di->op0(), c);
}


// Unsigned remainder instruction.
template <typename URV>
void
Hart<URV>::execRemu(const DecodedInst* di)
{
  URV a = intRegs_.read(di->op1());
  URV b = intRegs_.read(di->op2());
  URV c = a;  // Divide by zero remainder.
  if (b != 0)
    c = a % b;

  recordDivInst(di->op0(), peekIntReg(di->op0()));

  intRegs_.write(di->op0(), c);
}


template <typename URV>
void
Hart<URV>::execLwu(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }
  load<uint32_t>(di->op0(), di->op1(), di->op2As<int32_t>());
}


template <>
void
Hart<uint32_t>::execLd(const DecodedInst* di)
{
  illegalInst(di);
}


template <>
void
Hart<uint64_t>::execLd(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }
  load<uint64_t>(di->op0(), di->op1(), di->op2As<int32_t>());
}


template <typename URV>
void
Hart<URV>::execSd(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  unsigned rs1 = di->op1();

  URV base = intRegs_.read(rs1);
  URV addr = base + di->op2As<SRV>();
  URV value = intRegs_.read(di->op0());

  store<uint64_t>(rs1, base, addr, value);
}


template <typename URV>
void
Hart<URV>::execSlliw(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint32_t amount(di->op2());

  if (amount > 0x1f)
    {
      illegalInst(di);   // Bits 5 and 6 of immeidate must be zero.
      return;
    }

  int32_t word = int32_t(intRegs_.read(di->op1()));
  word <<= amount;

  SRV value = word; // Sign extend to 64-bit.
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSrliw(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint32_t amount(di->op2());

  if (amount > 0x1f)
    {
      illegalInst(di);   // Bits 5 and 6 of immediate must be zero.
      return;
    }

  uint32_t word = uint32_t(intRegs_.read(di->op1()));
  word >>= amount;

  SRV value = int32_t(word); // Sign extend to 64-bit.
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSraiw(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint32_t amount(di->op2());

  if (amount > 0x1f)
    {
      illegalInst(di);   // Bits 5 and 6 of immeddiate must be zero.
      return;
    }

  int32_t word = int32_t(intRegs_.read(di->op1()));
  word >>= amount;

  SRV value = word; // Sign extend to 64-bit.
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execAddiw(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  int32_t word = int32_t(intRegs_.read(di->op1()));
  word += di->op2As<int32_t>();
  SRV value = word;  // sign extend to 64-bits
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execAddw(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  int32_t word = int32_t(intRegs_.read(di->op1()) + intRegs_.read(di->op2()));
  SRV value = word;  // sign extend to 64-bits
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSubw(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  int32_t word = int32_t(intRegs_.read(di->op1()) - intRegs_.read(di->op2()));
  SRV value = word;  // sign extend to 64-bits
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSllw(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint32_t shift = intRegs_.read(di->op2()) & 0x1f;
  int32_t word = int32_t(intRegs_.read(di->op1()) << shift);
  SRV value = word;  // sign extend to 64-bits
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSrlw(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint32_t word = uint32_t(intRegs_.read(di->op1()));
  uint32_t shift = uint32_t(intRegs_.read(di->op2()) & 0x1f);
  word >>= shift;
  SRV value = int32_t(word);  // sign extend to 64-bits
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSraw(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  int32_t word = int32_t(intRegs_.read(di->op1()));
  uint32_t shift = uint32_t(intRegs_.read(di->op2()) & 0x1f);
  word >>= shift;
  SRV value = word;  // sign extend to 64-bits
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execMulw(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  int32_t word1 = int32_t(intRegs_.read(di->op1()));
  int32_t word2 = int32_t(intRegs_.read(di->op2()));
  int32_t word = int32_t(word1 * word2);
  SRV value = word;  // sign extend to 64-bits
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execDivw(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  int32_t word1 = int32_t(intRegs_.read(di->op1()));
  int32_t word2 = int32_t(intRegs_.read(di->op2()));

  int32_t word = -1;  // Divide by zero result
  if (word2 != 0)
    {
      int32_t minInt = int32_t(1) << 31;
      if (word1 == minInt and word2 == -1)
	word = word1;
      else
	word = word1 / word2;
    }

  recordDivInst(di->op0(), peekIntReg(di->op0()));

  SRV value = word;  // sign extend to 64-bits
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execDivuw(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint32_t word1 = uint32_t(intRegs_.read(di->op1()));
  uint32_t word2 = uint32_t(intRegs_.read(di->op2()));

  uint32_t word = ~uint32_t(0);  // Divide by zero result.
  if (word2 != 0)
    word = word1 / word2;

  recordDivInst(di->op0(), peekIntReg(di->op0()));

  URV value = SRV(int32_t(word));  // Sign extend to 64-bits
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execRemw(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  int32_t word1 = int32_t(intRegs_.read(di->op1()));
  int32_t word2 = int32_t(intRegs_.read(di->op2()));

  int32_t word = word1;  // Divide by zero remainder
  if (word2 != 0)
    {
      int32_t minInt = int32_t(1) << 31;
      if (word1 == minInt and word2 == -1)
	word = 0;   // Per spec: User-Level ISA, Version 2.3, Section 6.2
      else
	word = word1 % word2;
    }

  recordDivInst(di->op0(), peekIntReg(di->op0()));

  SRV value = word;  // sign extend to 64-bits
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execRemuw(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint32_t word1 = uint32_t(intRegs_.read(di->op1()));
  uint32_t word2 = uint32_t(intRegs_.read(di->op2()));

  uint32_t word = word1;  // Divide by zero remainder
  if (word2 != 0)
    word = word1 % word2;

  recordDivInst(di->op0(), peekIntReg(di->op0()));

  URV value = SRV(int32_t(word));  // Sign extend to 64-bits
  intRegs_.write(di->op0(), value);
}


template <typename URV>
inline
RoundingMode
Hart<URV>::effectiveRoundingMode(RoundingMode instMode)
{
  if (instMode != RoundingMode::Dynamic)
    return instMode;

  return getFpRoundingMode();
}


#ifdef FAST_SLOPPY

template <typename URV>
inline
void
Hart<URV>::updateAccruedFpBits(bool /* skipUnderflow */)
{
}


template <typename URV>
inline
void
Hart<URV>::markFsDirty()
{
}

#else

template <typename URV>
inline
void
Hart<URV>::updateAccruedFpBits(bool skipUnderflow)
{
  URV val = getFpFlags();
  URV prev = val;

  int flags = fetestexcept(FE_ALL_EXCEPT);
      
  if (flags)
    {
      if (flags & FE_INEXACT)
        val |= URV(FpFlags::Inexact);

      if (flags & FE_UNDERFLOW)
        if (not skipUnderflow)
          val |= URV(FpFlags::Underflow);

      if (flags & FE_OVERFLOW)
        val |= URV(FpFlags::Overflow);
      
      if (flags & FE_DIVBYZERO)
        val |= URV(FpFlags::DivByZero);

      if (flags & FE_INVALID)
        val |= URV(FpFlags::Invalid);

      if (val != prev)
        {
          setFpFlags(val);
          recordCsrWrite(CsrNumber::FCSR);
        }
    }
}


template <typename URV>
inline
void
Hart<URV>::markFsDirty()
{
  if (mstatusFs_ == FpFs::Dirty)
    return;

  URV val = csRegs_.peekMstatus();
  MstatusFields<URV> fields(val);
  fields.bits_.FS = unsigned(FpFs::Dirty);

  csRegs_.poke(CsrNumber::MSTATUS, fields.value_);

  updateCachedMstatusFields();
}

#endif


/// Map a RISCV rounding mode to an fetsetround constant.
static std::array<int, 5> riscvRoungingModeToFe =
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
  return riscvRoungingModeToFe.at(ix);
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


/// Clear the floating point flags in the machine running this
/// simulator. Do nothing in the simuated RISCV machine.
static
inline
void
clearSimulatorFpFlags()
{
  uint32_t val = _mm_getcsr();
  val &= ~uint32_t(0x3f);
  _mm_setcsr(val);
  // std::feclearexcept(FE_ALL_EXCEPT);
}


template <typename URV>
inline
bool
Hart<URV>::checkRoundingModeSp(const DecodedInst* di)
{
  if (not isRvf() or not isFpEnabled())
    {
      illegalInst(di);
      return false;
    }

  RoundingMode riscvMode = effectiveRoundingMode(di->roundingMode());
  if (riscvMode >= RoundingMode::Invalid1)
    {
      illegalInst(di);
      return false;
    }

  clearSimulatorFpFlags();
  setSimulatorRoundingMode(riscvMode);
  return true;
}


template <typename URV>
inline
bool
Hart<URV>::checkRoundingModeDp(const DecodedInst* di)
{
  if (not isRvd() or not isFpEnabled())
    {
      illegalInst(di);
      return false;
    }

  RoundingMode riscvMode = effectiveRoundingMode(di->roundingMode());
  if (riscvMode >= RoundingMode::Invalid1)
    {
      illegalInst(di);
      return false;
    }

  clearSimulatorFpFlags();
  setSimulatorRoundingMode(riscvMode);
  return true;
}


template <typename URV>
void
Hart<URV>::execFlw(const DecodedInst* di)
{
  if (not isRvf() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  uint32_t rd = di->op0(), rs1 = di->op1();
  SRV imm = di->op2As<SRV>();

  URV base = intRegs_.read(rs1);
  URV virtAddr = base + imm;

  ldStAddr_ = virtAddr;   // For reporting load addr in trace-mode.
  ldStAddrValid_ = true;  // For reporting load addr in trace-mode.

  auto secCause = SecondaryCause::NONE;
  auto cause = ExceptionCause::NONE;

#ifndef FAST_SLOPPY
  if (hasActiveTrigger())
    {
      if (ldStAddrTriggerHit(virtAddr, TriggerTiming::Before, true /*isLoad*/,
			     privMode_, isInterruptEnabled()))
	triggerTripped_ = true;
      if (triggerTripped_)
	return;
    }

  unsigned ldSize = 8;

  uint64_t addr = virtAddr;
  cause = determineLoadException(rs1, base, addr, ldSize, secCause);
  if (cause != ExceptionCause::NONE)
    {
      initiateLoadException(cause, virtAddr, secCause);
      return;
    }
#endif

  uint32_t word = 0;
  if (memory_.read(addr, word))
    {
      Uint32FloatUnion ufu(word);
      fpRegs_.writeSingle(rd, ufu.f);
      markFsDirty();
      return;
    }

  cause = ExceptionCause::LOAD_ACC_FAULT;
  secCause = SecondaryCause::LOAD_ACC_MEM_PROTECTION;
  if (isAddrMemMapped(addr))
    secCause = SecondaryCause::LOAD_ACC_PIC;

  initiateLoadException(cause, virtAddr, secCause);
}


template <typename URV>
void
Hart<URV>::execFsw(const DecodedInst* di)
{
  if (not isRvf() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  uint32_t rs1 = di->op1(), rs2 = di->op0();
  SRV imm = di->op2As<SRV>();

  URV base = intRegs_.read(rs1);
  URV addr = base + imm;

  // This operation does not check for proper NAN boxing. We read raw bits.
  uint64_t val = fpRegs_.readBitsRaw(rs2);

  store<uint32_t>(rs1, base, addr, uint32_t(val));
}


/// Use fused mutiply-add to perform x*y + z.
/// Set invalid to true if x and y are zero and infinity or
/// vice versa since RISCV consider that as an invalid operation.
static
float
fusedMultiplyAdd(float x, float y, float z, bool& invalid)
{
#ifdef __FP_FAST_FMA
  float res = x*y + z;
#else
  float res = std::fma(x, y, z);
#endif
  invalid = (std::isinf(x) and y == 0) or (x == 0 and std::isinf(y));
  return res;
}


/// Use fused mutiply-add to perform x*y + z.
static
double
fusedMultiplyAdd(double x, double y, double z, bool& invalid)
{
#ifdef __FP_FAST_FMA
  double res = x*y + z;
#else
  double res = std::fma(x, y, z);
#endif
  invalid = (std::isinf(x) and y == 0) or (x == 0 and std::isinf(y));
  return res;
}


/// Return the exponent bits of the given floating point value.
unsigned
spExponentBits(float sp)
{
  Uint32FloatUnion uf(sp);
  return (uf.u << 1) >> 24;
}


/// Return the exponent bits of the given double precision value.
unsigned
dpExponentBits(double dp)
{
  Uint64DoubleUnion ud(dp);
  return (ud.u << 1) >> 53;
}



template <typename URV>
void
Hart<URV>::execFmadd_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float f3 = fpRegs_.readSingle(di->op3());

  bool invalid = false;
  float res = fusedMultiplyAdd(f1, f2, f3, invalid);
  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  bool noUnderflow = spExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);
  if (invalid)
    setFcsrFlags(FpFlags::Invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmsub_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float f3 = -fpRegs_.readSingle(di->op3());

  bool invalid = false;
  float res = fusedMultiplyAdd(f1, f2, f3, invalid);
  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  bool noUnderflow = spExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);
  if (invalid)
    setFcsrFlags(FpFlags::Invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFnmsub_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = -fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float f3 = fpRegs_.readSingle(di->op3());

  bool invalid = false;
  float res = fusedMultiplyAdd(f1, f2, f3, invalid);
  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  bool noUnderflow = spExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);
  if (invalid)
    setFcsrFlags(FpFlags::Invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFnmadd_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  // we want -(f[op1] * f[op2]) - f[op3]

  float f1 = -fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float f3 = -fpRegs_.readSingle(di->op3());

  bool invalid = false;
  float res = fusedMultiplyAdd(f1, f2, f3, invalid);
  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  bool noUnderflow = spExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);
  if (invalid)
    setFcsrFlags(FpFlags::Invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFadd_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float res = f1 + f2;
  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  bool noUnderflow = spExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsub_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float res = f1 - f2;
  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  bool noUnderflow = spExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);
}


template <typename URV>
void
Hart<URV>::execFmul_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float res = f1 * f2;
  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  bool noUnderflow = spExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFdiv_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float res = f1 / f2;
  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  bool noUnderflow = spExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsqrt_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  float res = std::sqrt(f1);
  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  bool noUnderflow = spExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);
}


template <typename URV>
void
Hart<URV>::execFsgnj_s(const DecodedInst* di)
{
  if (not isRvf() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float res = std::copysignf(f1, f2);  // Magnitude of rs1 and sign of rs2
  fpRegs_.writeSingle(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsgnjn_s(const DecodedInst* di)
{
  if (not isRvf() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float res = std::copysignf(f1, f2);  // Magnitude of rs1 and sign of rs2
  res = -res;  // Magnitude of rs1 and negative the sign of rs2
  fpRegs_.writeSingle(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsgnjx_s(const DecodedInst* di)
{
  if (not isRvf() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());

  int sign1 = (std::signbit(f1) == 0) ? 0 : 1;
  int sign2 = (std::signbit(f2) == 0) ? 0 : 1;
  int sign = sign1 ^ sign2;

  float x = sign? -1 : 1;

  float res = std::copysignf(f1, x);  // Magnitude of rs1 and sign of x
  fpRegs_.writeSingle(di->op0(), res);

  markFsDirty();
}


/// Return true if given float is a signaling not-a-number.
static
bool
issnan(float f)
{
  if (std::isnan(f))
    {
      Uint32FloatUnion ufu(f);

      // Most sig bit of significant must be zero.
      return ((ufu.u >> 22) & 1) == 0;
    }
  return false;
}


/// Return true if given double is a signaling not-a-number.
static
bool
issnan(double d)
{
  if (std::isnan(d))
    {
      Uint64DoubleUnion udu(d);

      // Most sig bit of significant must be zero.
      return ((udu.u >> 51) & 1) == 0;
    }
  return false;
}


template <typename URV>
void
Hart<URV>::execFmin_s(const DecodedInst* di)
{
  if (not isRvf() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  float in1 = fpRegs_.readSingle(di->op1());
  float in2 = fpRegs_.readSingle(di->op2());
  float res = 0;

  bool isNan1 = std::isnan(in1), isNan2 = std::isnan(in2);
  if (isNan1 and isNan2)
    res = std::numeric_limits<float>::quiet_NaN();
  else if (isNan1)
    res = in2;
  else if (isNan2)
    res = in1;
  else
    res = std::fminf(in1, in2);

  if (issnan(in1) or issnan(in2))
    setFcsrFlags(FpFlags::Invalid);
  else if (std::signbit(in1) != std::signbit(in2) and in1 == in2)
    res = std::copysign(res, -1.0F);  // Make sure min(-0, +0) is -0.

  fpRegs_.writeSingle(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmax_s(const DecodedInst* di)
{
  if (not isRvf() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  float in1 = fpRegs_.readSingle(di->op1());
  float in2 = fpRegs_.readSingle(di->op2());
  float res = 0;

  bool isNan1 = std::isnan(in1), isNan2 = std::isnan(in2);
  if (isNan1 and isNan2)
    res = std::numeric_limits<float>::quiet_NaN();
  else if (isNan1)
    res = in2;
  else if (isNan2)
    res = in1;
  else
    res = std::fmaxf(in1, in2);

  if (issnan(in1) or issnan(in2))
    setFcsrFlags(FpFlags::Invalid);
  else if (std::signbit(in1) != std::signbit(in2) and in1 == in2)
    res = std::copysign(res, 1.0F);  // Make sure max(-0, +0) is +0.

  fpRegs_.writeSingle(di->op0(), res);

  markFsDirty();
}


/// Return sign bit (0 or 1) of given float. Works for all float
/// values including NANs and INFINITYs.
static
unsigned
signOf(float f)
{
  Uint32FloatUnion ufu(f);
  return ufu.u >> 31;
}


/// Return sign bit (0 or 1) of given double. Works for all double
/// values including NANs and INFINITYs.
static
unsigned
signOf(double d)
{
  Uint64DoubleUnion udu(d);
  return udu.u >> 63;
}


template <typename URV>
void
Hart<URV>::execFcvt_w_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  SRV result = 0;
  bool valid = false;

  int32_t minInt = int32_t(1) << 31;
  int32_t maxInt = (~uint32_t(0)) >> 1;

  unsigned signBit = signOf(f1);
  if (std::isinf(f1))
    result = signBit ? minInt : maxInt;
  else if (std::isnan(f1))
    result = maxInt;
  else
    {
      float near = std::nearbyint(f1);
      if (near >= float(maxInt))
	result = maxInt;
      else if (near < float(minInt))
	result = SRV(minInt);
      else
	{
	  valid = true;
          result = int32_t(std::lrintf(f1));
	}
    }

  intRegs_.write(di->op0(), result);

  updateAccruedFpBits(true);
  if (not valid)
    setFcsrFlags(FpFlags::Invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_wu_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  SRV result = 0;
  bool valid = false;
  bool exact = true;

  uint32_t maxUint32 = ~uint32_t(0);
  if (std::isnan(f1))
    {
      result = ~URV(0);
    }
  else
    {
      double near = std::nearbyint(f1);
      if (near > double(maxUint32))
        {
          result = ~URV(0);
        }
      else if (near == 0)
        {
          result = 0;
          valid = true;
          exact = near == f1;
        }
      else if (near < 0)
        {
          result = 0;
        }
      else
        {
          result = SRV(int32_t(std::lrint(f1)));
          valid = true;
          exact = near == f1;
        }
    }

  intRegs_.write(di->op0(), result);

  if (not valid)
    setFcsrFlags(FpFlags::Invalid);
  if (not exact)
    setFcsrFlags(FpFlags::Inexact);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmv_x_w(const DecodedInst* di)
{
  if (not isRvf() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  // This operation does not check for proper NAN boxing. We read raw bits.
  uint64_t v1 = fpRegs_.readBitsRaw(di->op1());
  int32_t s1 = v1;  // Keep lower 32 bits

  SRV value = SRV(s1); // Sign extend.

  intRegs_.write(di->op0(), value);
}

 
template <typename URV>
void
Hart<URV>::execFeq_s(const DecodedInst* di)
{
  if (not isRvf() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());

  URV res = 0;

  if (std::isnan(f1) or std::isnan(f2))
    {
      if (issnan(f1) or issnan(f2))
	setFcsrFlags(FpFlags::Invalid);
    }
  else
    res = (f1 == f2)? 1 : 0;

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFlt_s(const DecodedInst* di)
{
  if (not isRvf() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());

  URV res = 0;

  if (std::isnan(f1) or std::isnan(f2))
    setFcsrFlags(FpFlags::Invalid);
  else
    res = (f1 < f2)? 1 : 0;
    
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFle_s(const DecodedInst* di)
{
  if (not isRvf() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());

  URV res = 0;

  if (std::isnan(f1) or std::isnan(f2))
    setFcsrFlags(FpFlags::Invalid);
  else
    res = (f1 <= f2)? 1 : 0;

  intRegs_.write(di->op0(), res);
}


bool
mostSignificantFractionBit(float x)
{
  Uint32FloatUnion ufu(x);
  return (ufu.u >> 22) & 1;
}


bool
mostSignificantFractionBit(double x)
{
  Uint64DoubleUnion udu(x);
  return (udu.u >> 51) & 1;
}



template <typename URV>
void
Hart<URV>::execFclass_s(const DecodedInst* di)
{
  if (not isRvf() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  float f1 = fpRegs_.readSingle(di->op1());
  URV result = 0;

  bool pos = not std::signbit(f1);
  int type = std::fpclassify(f1);

  if (type == FP_INFINITE)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosInfinity);
      else
	result |= URV(FpClassifyMasks::NegInfinity);
    }
  else if (type == FP_NORMAL)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosNormal);
      else
	result |= URV(FpClassifyMasks::NegNormal);
    }
  else if (type == FP_SUBNORMAL)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosSubnormal);
      else
	result |= URV(FpClassifyMasks::NegSubnormal);
    }
  else if (type == FP_ZERO)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosZero);
      else
	result |= URV(FpClassifyMasks::NegZero);
    }
  else if (type == FP_NAN)
    {
      bool quiet = mostSignificantFractionBit(f1);
      if (quiet)
	result |= URV(FpClassifyMasks::QuietNan);
      else
	result |= URV(FpClassifyMasks::SignalingNan);
    }

  intRegs_.write(di->op0(), result);
}


template <typename URV>
void
Hart<URV>::execFcvt_s_w(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  int32_t i1 = intRegs_.read(di->op1());
  float res = float(i1);
  fpRegs_.writeSingle(di->op0(), res);

  bool noUnderflow = spExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_s_wu(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  uint32_t u1 = intRegs_.read(di->op1());
  float res = float(u1);
  fpRegs_.writeSingle(di->op0(), res);

  bool noUnderflow = spExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmv_w_x(const DecodedInst* di)
{
  if (not isRvf() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  uint32_t u1 = intRegs_.read(di->op1());

  Uint32FloatUnion ufu(u1);
  fpRegs_.writeSingle(di->op0(), ufu.f);

  markFsDirty();
}


template <>
void
Hart<uint32_t>::execFcvt_l_s(const DecodedInst* di)
{
  illegalInst(di);  // fcvt.l.s is not an RV32 instruction.
}


template <>
void
Hart<uint64_t>::execFcvt_l_s(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  SRV result = 0;
  bool valid = false;

  int64_t maxInt = (~uint64_t(0)) >> 1;
  int64_t minInt = int64_t(1) << 63;

  unsigned signBit = signOf(f1);
  if (std::isinf(f1))
    {
      if (signBit)
	result = minInt;
      else
	result = maxInt;
    }
  else if (std::isnan(f1))
    result = maxInt;
  else
    {
      double near = std::nearbyint(double(f1));
      if (near >= double(maxInt))
	result = maxInt;
      else if (near < double(minInt))
	result = minInt;
      else
	{
	  valid = true;
          result = std::lrint(f1);
	}
    }

  intRegs_.write(di->op0(), result);

  updateAccruedFpBits(true);
  if (not valid)
    setFcsrFlags(FpFlags::Invalid);

  markFsDirty();
}


template <>
void
Hart<uint32_t>::execFcvt_lu_s(const DecodedInst* di)
{
  illegalInst(di);  // RV32 does not have fcvt.lu.s
}


template <>
void
Hart<uint64_t>::execFcvt_lu_s(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  uint64_t result = 0;
  bool valid = false;
  bool exact = true;

  uint64_t maxUint = ~uint64_t(0);

  unsigned signBit = signOf(f1);
  if (std::isinf(f1))
    {
      if (signBit)
	result = 0;
      else
	result = maxUint;
    }
  else if (std::isnan(f1))
    result = maxUint;
  else
    {
      double near = std::nearbyint(double(f1));
      if (near == 0)
        {
          result = 0;
          valid = true;
          exact = near == f1;
        }
      else if (near < 0)
        {
          result = 0;
        }
      else
        {
          // Using "near > maxUint" will not work beacuse of rounding.
          if (near >= 2*double(uint64_t(1)<<63))
            result = maxUint;
          else
            {
              // std::lprint will produce an overflow if most sig bit
              // of result is 1 (it thinks there's an overflow).  We
              // compensate with the divide multiply by 2.
              if (f1 < (uint64_t(1) << 63))
                result = std::llrint(f1);
              else
                {
                  result = std::llrint(f1/2);
                  result *= 2;
                }
              valid = true;
              exact = near == f1;
            }
        }
    }

  intRegs_.write(di->op0(), result);

  if (not valid)
    setFcsrFlags(FpFlags::Invalid);
  if (not exact)
    setFcsrFlags(FpFlags::Inexact);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_s_l(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeSp(di))
    return;

  SRV i1 = intRegs_.read(di->op1());
  float res = float(i1);
  fpRegs_.writeSingle(di->op0(), res);

  bool noUnderflow = spExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_s_lu(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeSp(di))
    return;

  URV i1 = intRegs_.read(di->op1());
  float res = float(i1);
  fpRegs_.writeSingle(di->op0(), res);

  bool noUnderflow = spExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFld(const DecodedInst* di)
{
  if (not isRvd() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  unsigned rs1 = di->op1();
  SRV imm = di->op2As<SRV>();

  URV base = intRegs_.read(rs1);
  URV virtAddr = base + imm;

  ldStAddr_ = virtAddr;   // For reporting load addr in trace-mode.
  ldStAddrValid_ = true;  // For reporting load addr in trace-mode.

  auto secCause = SecondaryCause::NONE;
  auto cause = ExceptionCause::NONE;

#ifndef FAST_SLOPPY
  if (hasActiveTrigger())
    {
      if (ldStAddrTriggerHit(virtAddr, TriggerTiming::Before, true /*isLoad*/,
			     privMode_, isInterruptEnabled()))
	triggerTripped_ = true;
      if (triggerTripped_)
	return;
    }

  unsigned ldSize = 8;

  uint64_t addr = virtAddr;
  cause = determineLoadException(rs1, base, addr, ldSize, secCause);
  if (cause != ExceptionCause::NONE)
    {
      initiateLoadException(cause, virtAddr, secCause);
      return;
    }
#endif

  union UDU  // Unsigned double union: reinterpret bits as unsigned or double
  {
    uint64_t u;
    double d;
  };

  uint64_t val64 = 0;
  if (memory_.read(addr, val64))
    {
      UDU udu;
      udu.u = val64;
      fpRegs_.write(di->op0(), udu.d);

      markFsDirty();
      return;
    }

  cause = ExceptionCause::LOAD_ACC_FAULT;
  secCause = SecondaryCause::LOAD_ACC_MEM_PROTECTION;
  if (isAddrMemMapped(addr))
    secCause = SecondaryCause::LOAD_ACC_PIC;

  initiateLoadException(cause, virtAddr, secCause);
}


template <typename URV>
void
Hart<URV>::execFsd(const DecodedInst* di)
{
  if (not isRvd() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  uint32_t rs1 = di->op1();
  uint32_t rs2 = di->op0();

  URV base = intRegs_.read(rs1);
  URV addr = base + di->op2As<SRV>();
  double val = fpRegs_.read(rs2);

  union UDU  // Unsigned double union: reinterpret bits as unsigned or double
  {
    uint64_t u;
    double d;
  };

  UDU udu;
  udu.d = val;

  store<uint64_t>(rs1, base, addr, udu.u);
}


template <typename URV>
void
Hart<URV>::execFmadd_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double f1 = fpRegs_.read(di->op1());
  double f2 = fpRegs_.read(di->op2());
  double f3 = fpRegs_.read(di->op3());

  bool invalid = false;
  double res = fusedMultiplyAdd(f1, f2, f3, invalid);
  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  bool noUnderflow = dpExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);
  if (invalid)
    setFcsrFlags(FpFlags::Invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmsub_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double f1 = fpRegs_.read(di->op1());
  double f2 = fpRegs_.read(di->op2());
  double f3 = -fpRegs_.read(di->op3());

  bool invalid = false;
  double res = fusedMultiplyAdd(f1, f2, f3, invalid);

  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  bool noUnderflow = dpExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);
  if (invalid)
    setFcsrFlags(FpFlags::Invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFnmsub_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double f1 = -fpRegs_.read(di->op1());
  double f2 = fpRegs_.read(di->op2());
  double f3 = fpRegs_.read(di->op3());

  bool invalid = false;
  double res = fusedMultiplyAdd(f1, f2, f3, invalid);
  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  bool noUnderflow = dpExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);
  if (invalid)
    setFcsrFlags(FpFlags::Invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFnmadd_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  // we want -(f[op1] * f[op2]) - f[op3]

  double f1 = -fpRegs_.read(di->op1());
  double f2 = fpRegs_.read(di->op2());
  double f3 = -fpRegs_.read(di->op3());

  bool invalid = false;
  double res = fusedMultiplyAdd(f1, f2, f3, invalid);
  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  bool noUnderflow = dpExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);
  if (invalid)
    setFcsrFlags(FpFlags::Invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFadd_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());
  double res = d1 + d2;
  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  bool noUnderflow = dpExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsub_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());
  double res = d1 - d2;
  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  bool noUnderflow = dpExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmul_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());
  double res = d1 * d2;
  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  bool noUnderflow = dpExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFdiv_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());
  double res = d1 / d2;
  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  bool noUnderflow = dpExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsgnj_d(const DecodedInst* di)
{
  if (not isRvd() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());
  double res = copysign(d1, d2);  // Magnitude of rs1 and sign of rs2
  fpRegs_.write(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsgnjn_d(const DecodedInst* di)
{
  if (not isRvd() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());
  double res = copysign(d1, d2);  // Magnitude of rs1 and sign of rs2
  res = -res;  // Magnitude of rs1 and negative the sign of rs2
  fpRegs_.write(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsgnjx_d(const DecodedInst* di)
{
  if (not isRvd() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());

  int sign1 = (std::signbit(d1) == 0) ? 0 : 1;
  int sign2 = (std::signbit(d2) == 0) ? 0 : 1;
  int sign = sign1 ^ sign2;

  double x = sign? -1 : 1;

  double res = copysign(d1, x);  // Magnitude of rs1 and sign of x
  fpRegs_.write(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmin_d(const DecodedInst* di)
{
  if (not isRvd() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  double in1 = fpRegs_.read(di->op1());
  double in2 = fpRegs_.read(di->op2());
  double res = 0;

  bool isNan1 = std::isnan(in1), isNan2 = std::isnan(in2);
  if (isNan1 and isNan2)
    res = std::numeric_limits<double>::quiet_NaN();
  else if (isNan1)
    res = in2;
  else if (isNan2)
    res = in1;
  else
    res = fmin(in1, in2);

  if (issnan(in1) or issnan(in2))
    setFcsrFlags(FpFlags::Invalid);
  else if (std::signbit(in1) != std::signbit(in2) and in1 == in2)
    res = std::copysign(res, -1.0);  // Make sure min(-0, +0) is -0.

  fpRegs_.write(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmax_d(const DecodedInst* di)
{
  if (not isRvd() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  double in1 = fpRegs_.read(di->op1());
  double in2 = fpRegs_.read(di->op2());
  double res = 0;

  bool isNan1 = std::isnan(in1), isNan2 = std::isnan(in2);
  if (isNan1 and isNan2)
    res = std::numeric_limits<double>::quiet_NaN();
  else if (isNan1)
    res = in2;
  else if (isNan2)
    res = in1;
  else
    res = std::fmax(in1, in2);

  if (issnan(in1) or issnan(in2))
    setFcsrFlags(FpFlags::Invalid);
  else if (std::signbit(in1) != std::signbit(in2) and in1 == in2)
    res = std::copysign(res, 1.0);  // Make sure max(-0, +0) is +0.

  fpRegs_.write(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_d_s(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  double res = f1;
  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  bool noUnderflow = dpExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_s_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());
  float res = float(d1);
  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  bool noUnderflow = spExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsqrt_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());
  double res = std::sqrt(d1);
  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  bool noUnderflow = dpExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFle_d(const DecodedInst* di)
{
  if (not isRvd() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());

  URV res = 0;

  if (std::isnan(d1) or std::isnan(d2))
    setFcsrFlags(FpFlags::Invalid);
  else
    res = (d1 <= d2)? 1 : 0;

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFlt_d(const DecodedInst* di)
{
  if (not isRvd() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());

  URV res = 0;

  if (std::isnan(d1) or std::isnan(d2))
    setFcsrFlags(FpFlags::Invalid);
  else
    res = (d1 < d2)? 1 : 0;

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFeq_d(const DecodedInst* di)
{
  if (not isRvd() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());

  URV res = 0;

  if (std::isnan(d1) or std::isnan(d2))
    {
      if (issnan(d1) or issnan(d2))
	setFcsrFlags(FpFlags::Invalid);
    }
  else
    res = (d1 == d2)? 1 : 0;

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFcvt_w_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());
  SRV result = 0;
  bool valid = false;

  int32_t minInt = int32_t(1) << 31;
  int32_t maxInt = (~uint32_t(0)) >> 1;

  unsigned signBit = signOf(d1);
  if (std::isinf(d1))
    result = signBit? minInt : maxInt;
  else if (std::isnan(d1))
    result = maxInt;
  else
    {
      double near = std::nearbyint(d1);
      if (near > double(maxInt))
	result = maxInt;
      else if (near < double(minInt))
	result = minInt;
      else
	{
	  valid = true;
	  result = SRV(int32_t(std::lrint(d1)));
	}
    }

  intRegs_.write(di->op0(), result);

  updateAccruedFpBits(true);
  if (not valid)
    setFcsrFlags(FpFlags::Invalid);

  markFsDirty();
}



template <typename URV>
void
Hart<URV>::execFcvt_wu_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());
  SRV result = 0;
  bool valid = false;
  bool exact = true;

  uint32_t maxUint32 = ~uint32_t(0);
  if (std::isnan(d1))
    {
      result = ~URV(0);
    }
  else
    {
      double near = std::nearbyint(d1);
      if (near > double(maxUint32))
        {
          result = ~URV(0);
        }
      else if (near == 0)
        {
          result = 0;
          valid = true;
          exact = near == d1;
        }
      else if (near < 0)
        {
          result = 0;
        }
      else
        {
          result = SRV(int32_t(std::lrint(d1)));
          valid = true;
          exact = near == d1;
        }
    }

  intRegs_.write(di->op0(), result);

  if (not valid)
    setFcsrFlags(FpFlags::Invalid);
  if (not exact)
    setFcsrFlags(FpFlags::Inexact);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_d_w(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  int32_t i1 = intRegs_.read(di->op1());
  double res = i1;
  fpRegs_.write(di->op0(), res);

  bool noUnderflow = dpExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_d_wu(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  uint32_t i1 = intRegs_.read(di->op1());
  double res = i1;
  fpRegs_.write(di->op0(), res);

  bool noUnderflow = dpExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFclass_d(const DecodedInst* di)
{
  if (not isRvd() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  double d1 = fpRegs_.read(di->op1());
  URV result = 0;

  bool pos = not std::signbit(d1);
  int type = std::fpclassify(d1);

  if (type == FP_INFINITE)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosInfinity);
      else
	result |= URV(FpClassifyMasks::NegInfinity);
    }
  else if (type == FP_NORMAL)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosNormal);
      else
	result |= URV(FpClassifyMasks::NegNormal);
    }
  else if (type == FP_SUBNORMAL)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosSubnormal);
      else
	result |= URV(FpClassifyMasks::NegSubnormal);
    }
  else if (type == FP_ZERO)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosZero);
      else
	result |= URV(FpClassifyMasks::NegZero);
    }
  else if(type == FP_NAN)
    {
      bool quiet = mostSignificantFractionBit(d1);
      if (quiet)
	result |= URV(FpClassifyMasks::QuietNan);
      else
	result |= URV(FpClassifyMasks::SignalingNan);
    }

  intRegs_.write(di->op0(), result);
}


template <>
void
Hart<uint32_t>::execFcvt_l_d(const DecodedInst* di)
{
  illegalInst(di);  // fcvt.l.d not available in RV32
}


template <>
void
Hart<uint64_t>::execFcvt_l_d(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeDp(di))
    return;

  double f1 = fpRegs_.read(di->op1());
  SRV result = 0;
  bool valid = false;

  int64_t maxInt = (~uint64_t(0)) >> 1;
  int64_t minInt = int64_t(1) << 63;

  unsigned signBit = signOf(f1);
  if (std::isinf(f1))
    {
      if (signBit)
	result = minInt;
      else
	result = maxInt;
    }
  else if (std::isnan(f1))
    result = maxInt;
  else
    {
      double near = std::nearbyint(f1);

      // Note "near > double(maxInt)" will not work because of
      // rounding.
      if (near >= double(uint64_t(1) << 63))
	result = maxInt;
      else if (near < double(minInt))
	result = minInt;
      else
	{
	  valid = true;
	  result = std::lrint(f1);
	}
    }

  intRegs_.write(di->op0(), result);

  updateAccruedFpBits(true);
  if (not valid)
    setFcsrFlags(FpFlags::Invalid);

  markFsDirty();
}


template <>
void
Hart<uint32_t>::execFcvt_lu_d(const DecodedInst* di)
{
  illegalInst(di);  /// fcvt.lu.d is not available in RV32.
}


template <>
void
Hart<uint64_t>::execFcvt_lu_d(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeDp(di))
    return;

  double f1 = fpRegs_.read(di->op1());
  uint64_t result = 0;
  bool valid = false;
  bool exact = true;

  uint64_t maxUint = ~uint64_t(0);

  unsigned signBit = signOf(f1);
  if (std::isinf(f1))
    {
      if (signBit)
	result = 0;
      else
	result = maxUint;
    }
  else if (std::isnan(f1))
    result = maxUint;
  else
    {
      double near = std::nearbyint(f1);
      if (near == 0)
        {
          result = 0;
          valid = true;
          exact = near == f1;
        }
      else if (near < 0)
        {
          result = 0;
        }
      else
        {
          // Using "near > maxUint" will not work beacuse of rounding.
          if (near >= 2*double(uint64_t(1)<<63))
            result = maxUint;
          else
            {
              // std::llrint will produce an overflow if most sig bit
              // of result is 1 (it thinks there's an overflow).  We
              // compensate with the divide multiply by 2.
              if (f1 < (uint64_t(1) << 63))
                result = std::llrint(f1);
              else
                {
                  result = std::llrint(f1/2);
                  result *= 2;
                }
              valid = true;
              exact = near == f1;
            }
        }
    }

  intRegs_.write(di->op0(), result);

  if (not valid)
    setFcsrFlags(FpFlags::Invalid);
  if (not exact)
    setFcsrFlags(FpFlags::Inexact);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_d_l(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeDp(di))
    return;

  SRV i1 = intRegs_.read(di->op1());
  double res = double(i1);
  fpRegs_.write(di->op0(), res);

  bool noUnderflow = dpExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_d_lu(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeDp(di))
    return;

  URV i1 = intRegs_.read(di->op1());
  double res = double(i1);
  fpRegs_.write(di->op0(), res);

  bool noUnderflow = dpExponentBits(res) != 0;
  updateAccruedFpBits(noUnderflow);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmv_d_x(const DecodedInst* di)
{
  if (not isRv64() or not isRvd() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  uint64_t u1 = intRegs_.read(di->op1());

  union UDU  // Unsigned double union: reinterpret bits as unsigned or double
  {
    uint64_t u;
    double d;
  };

  UDU udu;
  udu.u = u1;

  fpRegs_.write(di->op0(), udu.d);

  markFsDirty();
}


// In 32-bit harts, fmv_x_d is an illegal instruction.
template <>
void
Hart<uint32_t>::execFmv_x_d(const DecodedInst* di)
{
  illegalInst(di);
}


template <>
void
Hart<uint64_t>::execFmv_x_d(const DecodedInst* di)
{
  if (not isRv64() or not isRvd() or not isFpEnabled())
    {
      illegalInst(di);
      return;
    }

  uint64_t v1 = fpRegs_.readBitsRaw(di->op1());
  intRegs_.write(di->op0(), v1);
}


template <typename URV>
template <typename LOAD_TYPE>
bool
Hart<URV>::loadReserve(uint32_t rd, uint32_t rs1, uint64_t& physAddr)
{
  enableWideLdStMode(false);

  URV virtAddr = intRegs_.read(rs1);

  ldStAddr_ = virtAddr;   // For reporting load addr in trace-mode.
  ldStAddrValid_ = true;  // For reporting load addr in trace-mode.

  if (loadQueueEnabled_)
    removeFromLoadQueue(rs1, false);

  if (hasActiveTrigger())
    {
      typedef TriggerTiming Timing;

      bool isLd = true;
      if (ldStAddrTriggerHit(virtAddr, Timing::Before, isLd,
                             privMode_, isInterruptEnabled()))
	triggerTripped_ = true;
      if (triggerTripped_)
	return false;
    }

  // Unsigned version of LOAD_TYPE
  typedef typename std::make_unsigned<LOAD_TYPE>::type ULT;

  auto secCause = SecondaryCause::NONE;
  unsigned ldSize = sizeof(LOAD_TYPE);
  uint64_t addr = virtAddr;
  auto cause = determineLoadException(rs1, virtAddr, addr, ldSize, secCause);
  if (cause != ExceptionCause::NONE)
    {
      if (cause == ExceptionCause::LOAD_ADDR_MISAL and
	  misalAtomicCauseAccessFault_)
        {
          cause = ExceptionCause::LOAD_ACC_FAULT;
          secCause = SecondaryCause::LOAD_ACC_AMO;
        }
    }

  // Address outside DCCM causes an exception (this is swerv specific).
  bool fail = amoInDccmOnly_ and not isAddrInDccm(addr);

  // Access must be naturally aligned.
  if ((addr & (ldSize - 1)) != 0)
    fail = true;

  if (fail)
    {
      // AMO secondary cause has priority over ECC.
      if (cause == ExceptionCause::NONE or
          secCause == SecondaryCause::LOAD_ACC_DOUBLE_ECC)
        {
          // Per spec cause is store-access-fault.
          cause = ExceptionCause::LOAD_ACC_FAULT;
          secCause = SecondaryCause::LOAD_ACC_AMO;
        }
    }

  if (cause != ExceptionCause::NONE)
    {
      initiateLoadException(cause, virtAddr, secCause);
      return false;
    }

  ULT uval = 0;
  if (not memory_.read(addr, uval))
    {  // Should never happen.
      initiateLoadException(cause, virtAddr, secCause);
      return false;
    }      

  URV value = uval;
  if (not std::is_same<ULT, LOAD_TYPE>::value)
    value = SRV(LOAD_TYPE(uval)); // Sign extend.

  // Put entry in load queue with value of rd before this load.
  if (loadQueueEnabled_)
    putInLoadQueue(ldSize, addr, rd, peekIntReg(rd));

  intRegs_.write(rd, value);

  physAddr = addr;
  return true;
}


template <typename URV>
void
Hart<URV>::execLr_w(const DecodedInst* di)
{
  std::lock_guard<std::mutex> lock(memory_.lrMutex_);
  uint64_t physAddr = 0;
  if (not loadReserve<int32_t>(di->op0(), di->op1(), physAddr))
    return;

  memory_.makeLr(hartIx_, physAddr, 4 /*size*/);
}


/// STORE_TYPE is either uint32_t or uint64_t.
template <typename URV>
template <typename STORE_TYPE>
bool
Hart<URV>::storeConditional(unsigned rs1, URV virtAddr, STORE_TYPE storeVal)
{
  enableWideLdStMode(false);

  ldStAddr_ = virtAddr;   // For reporting ld/st addr in trace-mode.
  ldStAddrValid_ = true;  // For reporting ld/st addr in trace-mode.

  // ld/st-address or instruction-address triggers have priority over
  // ld/st access or misaligned exceptions.
  bool hasTrig = hasActiveTrigger();
  TriggerTiming timing = TriggerTiming::Before;
  bool isLoad = false;
  if (hasTrig)
    if (ldStAddrTriggerHit(virtAddr, timing, isLoad, privMode_,
                           isInterruptEnabled()))
      triggerTripped_ = true;

  // Misaligned store causes an exception.
  constexpr unsigned alignMask = sizeof(STORE_TYPE) - 1;
  bool misal = virtAddr & alignMask;
  misalignedLdSt_ = misal;

  auto secCause = SecondaryCause::NONE;
  uint64_t addr = virtAddr;
  auto cause = determineStoreException(rs1, virtAddr, addr, storeVal, secCause);

  bool fail = misal or (amoInDccmOnly_ and not isAddrInDccm(addr));
  if (fail)
    {
      // AMO secondary cause has priority over ECC.
      if (cause == ExceptionCause::NONE or
          secCause == SecondaryCause::STORE_ACC_DOUBLE_ECC)
        {
          // Per spec cause is store-access-fault.
          cause = ExceptionCause::STORE_ACC_FAULT;
          secCause = SecondaryCause::STORE_ACC_AMO;
        }
    }

  // If no exception: consider store-data  trigger
  if (cause == ExceptionCause::NONE and hasTrig)
    if (ldStDataTriggerHit(storeVal, timing, isLoad, privMode_,
                           isInterruptEnabled()))
      triggerTripped_ = true;
  if (triggerTripped_)
    return false;

  if (cause != ExceptionCause::NONE)
    {
      initiateStoreException(cause, virtAddr, secCause);
      return false;
    }

  if (not memory_.hasLr(hartIx_, addr))
    return false;

  if (memory_.write(hartIx_, addr, storeVal))
    {
      invalidateDecodeCache(addr, sizeof(STORE_TYPE));

      // If we write to special location, end the simulation.
      if (toHostValid_ and addr == toHost_ and storeVal != 0)
        throw CoreException(CoreException::Stop, "write to to-host",
                            toHost_, storeVal);
      return true;
    }

  // Should never happen.
  secCause = SecondaryCause::STORE_ACC_AMO;
  initiateStoreException(ExceptionCause::STORE_ACC_FAULT, virtAddr, secCause);
  return false;
}


template <typename URV>
void
Hart<URV>::execSc_w(const DecodedInst* di)
{
  std::lock_guard<std::mutex> lock(memory_.lrMutex_);

  uint32_t rs1 = di->op1();
  URV value = intRegs_.read(di->op2());
  URV addr = intRegs_.read(rs1);

  uint64_t prevCount = exceptionCount_;

  bool ok = storeConditional(rs1, addr, uint32_t(value));
  memory_.invalidateLr(hartIx_);

  if (ok)
    {
      memory_.invalidateOtherHartLr(hartIx_, addr, 4);
      intRegs_.write(di->op0(), 0); // success
      return;
    }

  // If exception or trigger tripped then rd is not modified.
  if (triggerTripped_ or exceptionCount_ != prevCount)
    return;

  intRegs_.write(di->op0(), 1);  // fail
}


template <typename URV>
void
Hart<URV>::execAmoadd_w(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad32(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(di->op2());
      URV result = rs2Val + rdVal;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoswap_w(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad32(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(di->op2());
      URV result = rs2Val;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoxor_w(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad32(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(di->op2());
      URV result = rs2Val ^ rdVal;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoor_w(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad32(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(di->op2());
      URV result = rs2Val | rdVal;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoand_w(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad32(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(di->op2());
      URV result = rs2Val & rdVal;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmomin_w(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad32(rs1, loadedValue);

  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(di->op2());
      URV result = (SRV(rs2Val) < SRV(rdVal))? rs2Val : rdVal;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmominu_w(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad32(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(di->op2());

      uint32_t w1 = uint32_t(rs2Val), w2 = uint32_t(rdVal);
      uint32_t result = (w1 < w2)? w1 : w2;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmomax_w(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad32(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(di->op2());
      URV result = (SRV(rs2Val) > SRV(rdVal))? rs2Val : rdVal;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmomaxu_w(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad32(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(di->op2());

      uint32_t w1 = uint32_t(rs2Val), w2 = uint32_t(rdVal);

      URV result = (w1 > w2)? w1 : w2;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execLr_d(const DecodedInst* di)
{
  std::lock_guard<std::mutex> lock(memory_.lrMutex_);

  uint64_t physAddr = 0;
  if (not loadReserve<int64_t>(di->op0(), di->op1(), physAddr))
    return;

  memory_.makeLr(hartIx_, physAddr, 8 /*size*/);
}


template <typename URV>
void
Hart<URV>::execSc_d(const DecodedInst* di)
{
  std::lock_guard<std::mutex> lock(memory_.lrMutex_);

  uint32_t rs1 = di->op1();
  URV value = intRegs_.read(di->op2());
  URV addr = intRegs_.read(rs1);

  uint64_t prevCount = exceptionCount_;

  bool ok = storeConditional(rs1, addr, uint64_t(value));
  memory_.invalidateLr(hartIx_);

  if (ok)
    {
      memory_.invalidateOtherHartLr(hartIx_, addr, 8);
      intRegs_.write(di->op0(), 0); // success
      return;
    }

  // If exception or trigger tripped then rd is not modified.
  if (triggerTripped_ or exceptionCount_ != prevCount)
    return;

  intRegs_.write(di->op0(), 1);  // fail
}


template <typename URV>
void
Hart<URV>::execAmoadd_d(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad64(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(di->op2());
      URV result = rs2Val + rdVal;

      bool storeOk = store<uint32_t>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoswap_d(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad64(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(di->op2());
      URV result = rs2Val;

      bool storeOk = store<URV>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoxor_d(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad64(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(di->op2());
      URV result = rs2Val ^ rdVal;

      bool storeOk = store<URV>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoor_d(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad64(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(di->op2());
      URV result = rs2Val | rdVal;

      bool storeOk = store<URV>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoand_d(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad64(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(di->op2());
      URV result = rs2Val & rdVal;

      bool storeOk = store<URV>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmomin_d(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad64(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(di->op2());
      URV result = (SRV(rs2Val) < SRV(rdVal))? rs2Val : rdVal;

      bool storeOk = store<URV>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmominu_d(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad64(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(di->op2());
      URV result = (rs2Val < rdVal)? rs2Val : rdVal;

      bool storeOk = store<URV>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmomax_d(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad64(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(di->op2());
      URV result = (SRV(rs2Val) > SRV(rdVal))? rs2Val : rdVal;

      bool storeOk = store<URV>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmomaxu_d(const DecodedInst* di)
{
  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rs1 = di->op1();
  bool loadOk = amoLoad64(rs1, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(di->op2());
      URV result = (rs2Val > rdVal)? rs2Val : rdVal;

      bool storeOk = store<URV>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(di->op0(), rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execClz(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());

  if (v1 == 0)
    v1 = 8*sizeof(URV);
  else
    {
      if constexpr (sizeof(URV) == 4)
        v1 = __builtin_clz(v1);
      else
        v1 = __builtin_clzl(v1);
    }

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execCtz(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());

  if constexpr (sizeof(URV) == 4)
    v1 = __builtin_ctz(v1);
  else
    v1 = __builtin_ctzl(v1);

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execPcnt(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV res = __builtin_popcount(v1);
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execAndn(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV res = v1 & ~v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execOrn(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV res = v1 | ~v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execXnor(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV res = v1 ^ ~v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execSlo(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV mask = intRegs_.shiftMask();
  URV shift = intRegs_.read(di->op2()) & mask;

  URV v1 = intRegs_.read(di->op1());
  URV res = ~((~v1) << shift);
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execSro(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV mask = intRegs_.shiftMask();
  URV shift = intRegs_.read(di->op2()) & mask;

  URV v1 = intRegs_.read(di->op1());
  URV res = ~((~v1) >> shift);
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execSloi(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV imm = di->op2();
  if (not checkShiftImmediate(di, imm))
    return;

  URV v1 = intRegs_.read(di->op1());
  URV res = ~((~v1) << imm);
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execSroi(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  uint32_t imm = di->op2();
  if (not checkShiftImmediate(di, imm))
    return;

  URV v1 = intRegs_.read(di->op1());
  URV res = ~((~v1) >> imm);
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execMin(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  SRV v1 = intRegs_.read(di->op1());
  SRV v2 = intRegs_.read(di->op2());
  SRV res = v1 < v2? v1 : v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execMax(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  SRV v1 = intRegs_.read(di->op1());
  SRV v2 = intRegs_.read(di->op2());
  SRV res = v1 > v2? v1 : v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execMinu(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV res = v1 < v2? v1 : v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execMaxu(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV res = v1 > v2? v1 : v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execRol(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV mask = intRegs_.shiftMask();
  URV rot = intRegs_.read(di->op2()) & mask;  // Rotate amount

  URV v1 = intRegs_.read(di->op1());
  URV res = (v1 << rot) | (v1 >> (mxlen_ - rot));

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execRor(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV mask = intRegs_.shiftMask();
  URV rot = intRegs_.read(di->op2()) & mask;  // Rotate amount

  URV v1 = intRegs_.read(di->op1());
  URV res = (v1 >> rot) | (v1 << (mxlen_ - rot));

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execRori(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV rot = di->op2();
  if (not checkShiftImmediate(di, rot))
    return;

  URV v1 = intRegs_.read(di->op1());
  URV res = (v1 >> rot) | (v1 << (mxlen_ - rot));

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execRev8(const DecodedInst* di)
{
  // Byte swap.

  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());

  if constexpr (sizeof(URV) == 4)
    v1 = __builtin_bswap32(v1);
  else
    v1 = __builtin_bswap64(v1);

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execRev(const DecodedInst* di)
{
  // Bit reverse.

  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());

  if constexpr (sizeof(URV) == 4)
    {
      v1 = ((v1 & 0xaaaaaaaa) >> 1) | ((v1 & 0x55555555) << 1);
      v1 = ((v1 & 0xcccccccc) >> 2) | ((v1 & 0x33333333) << 2);
      v1 = ((v1 & 0xf0f0f0f0) >> 4) | ((v1 & 0x0f0f0f0f) << 4);
      v1 = __builtin_bswap32(v1);
    }
  else
    {
      v1 = ((v1 & 0xaaaaaaaaaaaaaaaa) >> 1) | ((v1 & 0x5555555555555555) << 1);
      v1 = ((v1 & 0xcccccccccccccccc) >> 2) | ((v1 & 0x3333333333333333) << 2);
      v1 = ((v1 & 0xf0f0f0f0f0f0f0f0) >> 4) | ((v1 & 0x0f0f0f0f0f0f0f0f) << 4);
      v1 = __builtin_bswap64(v1);
    }

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execPack(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  unsigned halfXlen = mxlen_ >> 1;
  URV lower = (intRegs_.read(di->op1()) << halfXlen) >> halfXlen;
  URV upper = intRegs_.read(di->op2()) << halfXlen;
  URV res = upper | lower;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execAddwu(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV value = uint32_t(intRegs_.read(di->op1()) + intRegs_.read(di->op2()));
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSubwu(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV value = uint32_t(intRegs_.read(di->op1()) - intRegs_.read(di->op2()));
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execAddiwu(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV value = uint32_t(intRegs_.read(di->op1()) + di->op2As<int32_t>());
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSext_b(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  int8_t byte = intRegs_.read(di->op1());
  SRV value = byte;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSext_h(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  int16_t half = intRegs_.read(di->op1());
  SRV value = half;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execAddu_w(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV value =   intRegs_.read(di->op1()) + uint32_t(intRegs_.read(di->op2()));
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSubu_w(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV value =   intRegs_.read(di->op1()) - uint32_t(intRegs_.read(di->op2()));
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSlliu_w(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  uint32_t amount(di->op2());

  if (amount > 0x1f)
    {
      illegalInst(di);   // Bits 5 and 6 of immediate must be zero.
      return;
    }

  uint32_t word = int32_t(intRegs_.read(di->op1()));
  word <<= amount;

  URV value = word;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execPackh(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV lower = intRegs_.read(di->op1()) & 0xff;
  URV upper = (intRegs_.read(di->op2()) & 0xff) << 8;
  URV value = lower | upper;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execPacku(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  unsigned halfXlen = mxlen_ >> 1;

  URV lower = intRegs_.read(di->op1()) >> halfXlen;
  URV upper = (intRegs_.read(di->op2()) >> halfXlen) << halfXlen;
  URV value = lower | upper;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execPackw(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  URV lower = intRegs_.read(di->op1()) & 0xffff;
  URV upper = (intRegs_.read(di->op2()) & 0xffff) << 16;
  URV value = lower | upper;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execPackuw(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  if (not isRv64())
    {
      illegalInst(di);
      return;
    }


  unsigned halfXlen = mxlen_ >> 1;

  URV lower = (intRegs_.read(di->op1()) >> halfXlen);
  URV upper = (intRegs_.read(di->op2()) >> halfXlen);

  lower = lower & 0xffff;
  upper = (upper & 0xffff) << 16;
  URV value = lower | upper;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execGrev(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  if constexpr (sizeof(URV) == 4)
    {
      unsigned shamt = v2 & 31;
      if (shamt & 1)
        v1 = ((v1 & 0x55555555) << 1)  | ((v1 & 0xaaaaaaaa) >> 1);
      if (shamt & 2)
        v1 = ((v1 & 0x33333333) << 2)  | ((v1 & 0xcccccccc) >> 2);
      if (shamt & 4)
        v1 = ((v1 & 0x0f0f0f0f) << 4)  | ((v1 & 0xf0f0f0f0) >> 4);
      if (shamt & 8)
        v1 = ((v1 & 0x00ff00ff) << 8)  | ((v1 & 0xff00ff00) >> 8);
      if (shamt & 16)
        v1 = ((v1 & 0x0000ffff) << 16) | ((v1 & 0xffff0000) >> 16);
    }
  else
    {
      int shamt = v2 & 63;
      if (shamt & 1)
        v1 = ((v1 & 0x5555555555555555ll) << 1)  | ((v1 & 0xaaaaaaaaaaaaaaaall) >> 1);
      if (shamt & 2)
        v1 = ((v1 & 0x3333333333333333ll) << 2)  | ((v1 & 0xccccccccccccccccll) >> 2);
      if (shamt & 4)
        v1 = ((v1 & 0x0f0f0f0f0f0f0f0fll) << 4)  | ((v1 & 0xf0f0f0f0f0f0f0f0ll) >> 4);
      if (shamt & 8)
        v1 = ((v1 & 0x00ff00ff00ff00ffll) << 8)  | ((v1 & 0xff00ff00ff00ff00ll) >> 8);
      if (shamt & 16)
        v1 = ((v1 & 0x0000ffff0000ffffll) << 16) | ((v1 & 0xffff0000ffff0000ll) >> 16);
      if (shamt & 32)
        v1 = ((v1 & 0x00000000ffffffffll) << 32) | ((v1 & 0xffffffff00000000ll) >> 32);
    }

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execGrevi(const DecodedInst* di)
{
  URV shamt = di->op2();

  bool zbb = false;  // True if variant is also a zbb instruction.
  if (isRv64())
    zbb = shamt == 0x38 or shamt == 0x3f;  // rev8 and rev are also in zbb
  else
    zbb = shamt == 0x18 or shamt == 0x1f;  // rev8 and rev are also in zbb

  bool illegal = not isRvzbp();
  if (zbb)
    illegal = not isRvzbb() and not isRvzbp();

  if (illegal)
    {
      illegalInst(di);
      return;
    }

  if (not checkShiftImmediate(di, shamt))
    return;

  URV v1 = intRegs_.read(di->op1());

  if constexpr (sizeof(URV) == 4)
    {
      if (shamt & 1)
        v1 = ((v1 & 0x55555555) << 1)  | ((v1 & 0xaaaaaaaa) >> 1);
      if (shamt & 2)
        v1 = ((v1 & 0x33333333) << 2)  | ((v1 & 0xcccccccc) >> 2);
      if (shamt & 4)
        v1 = ((v1 & 0x0f0f0f0f) << 4)  | ((v1 & 0xf0f0f0f0) >> 4);
      if (shamt & 8)
        v1 = ((v1 & 0x00ff00ff) << 8)  | ((v1 & 0xff00ff00) >> 8);
      if (shamt & 16)
        v1 = ((v1 & 0x0000ffff) << 16) | ((v1 & 0xffff0000) >> 16);
    }
  else
    {
      if (shamt & 1)
        v1 = ((v1 & 0x5555555555555555ll) << 1)  | ((v1 & 0xaaaaaaaaaaaaaaaall) >> 1);
      if (shamt & 2)
        v1 = ((v1 & 0x3333333333333333ll) << 2)  | ((v1 & 0xccccccccccccccccll) >> 2);
      if (shamt & 4)
        v1 = ((v1 & 0x0f0f0f0f0f0f0f0fll) << 4)  | ((v1 & 0xf0f0f0f0f0f0f0f0ll) >> 4);
      if (shamt & 8)
        v1 = ((v1 & 0x00ff00ff00ff00ffll) << 8)  | ((v1 & 0xff00ff00ff00ff00ll) >> 8);
      if (shamt & 16)
        v1 = ((v1 & 0x0000ffff0000ffffll) << 16) | ((v1 & 0xffff0000ffff0000ll) >> 16);
      if (shamt & 32)
        v1 = ((v1 & 0x00000000ffffffffll) << 32) | ((v1 & 0xffffffff00000000ll) >> 32);
    }

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execGorci(const DecodedInst* di)
{
  URV shamt = di->op2();

  bool zbb = false;  // True if variant is also a zbb instruction.
  if (isRv64())
    zbb = shamt == 0x7 or shamt == 0x30;  // orc.b and orc16 are also in zbb
  else
    zbb = shamt == 0x7 or shamt == 0x10;  // orc.b and orc16 are also in zbb

  bool illegal = not isRvzbp();
  if (zbb)
    illegal = not isRvzbb() and not isRvzbp();

  if (illegal)
    {
      illegalInst(di);
      return;
    }

  if (not checkShiftImmediate(di, shamt))
    return;

  URV v1 = intRegs_.read(di->op1());

  if constexpr (sizeof(URV) == 4)
    {
      if (shamt & 1)
        v1 |= ((v1 & 0xaaaaaaaa) >>  1) | ((v1 & 0x55555555) <<  1);
      if (shamt & 2)
        v1 |= ((v1 & 0xcccccccc) >>  2) | ((v1 & 0x33333333) <<  2);
      if (shamt & 4)
        v1 |= ((v1 & 0xf0f0f0f0) >>  4) | ((v1 & 0x0f0f0f0f) <<  4);
      if (shamt & 8)
        v1 |= ((v1 & 0xff00ff00) >>  8) | ((v1 & 0x00ff00ff) <<  8);
      if (shamt & 16)
        v1 |= ((v1 & 0xffff0000) >> 16) | ((v1 & 0x0000ffff) << 16);
    }
  else
    {
      if (shamt & 1)
        v1 |= ((v1 & 0xaaaaaaaaaaaaaaaa) >>  1) | ((v1 & 0x5555555555555555) <<  1);
      if (shamt & 2)
        v1 |= ((v1 & 0xcccccccccccccccc) >>  2) | ((v1 & 0x3333333333333333) <<  2);
      if (shamt & 4)
        v1 |= ((v1 & 0xf0f0f0f0f0f0f0f0) >>  4) | ((v1 & 0x0f0f0f0f0f0f0f0f) <<  4);
      if (shamt & 8)
        v1 |= ((v1 & 0xff00ff00ff00ff00) >>  8) | ((v1 & 0x00ff00ff00ff00ff) <<  8);
      if (shamt & 16)
        v1 |= ((v1 & 0xffff0000ffff0000) >> 16) | ((v1 & 0x0000ffff0000ffff) << 16);
      if (shamt & 32)
        v1 |= ((v1 & 0xffffffff00000000) >> 32) | ((v1 & 0x00000000ffffffff) << 32);
    }

  intRegs_.write(di->op0(), v1);
}


static
uint32_t
shuffleStage32(uint32_t src, uint32_t maskL, uint32_t maskR, unsigned n)
{
  uint32_t x = src & ~(maskL | maskR);
  x |= ((src << n) & maskL) | ((src >> n) & maskR);
  return x;
}


static
uint32_t
shuffle32(uint32_t x, unsigned shamt)
{
  if (shamt & 8)
    x = shuffleStage32(x, 0x00ff0000, 0x0000ff00, 8);
  if (shamt & 4)
    x = shuffleStage32(x, 0x0f000f00, 0x00f000f0, 4);
  if (shamt & 2)
    x = shuffleStage32(x, 0x30303030, 0x0c0c0c0c, 2);
  if (shamt & 1)
    x = shuffleStage32(x, 0x44444444, 0x22222222, 1);

  return x;
}


static
uint32_t
unshuffle32(uint32_t x, unsigned shamt)
{
  if (shamt & 1)
    x = shuffleStage32(x, 0x44444444, 0x22222222, 1);
  if (shamt & 2)
    x = shuffleStage32(x, 0x30303030, 0x0c0c0c0c, 2);
  if (shamt & 4)
    x = shuffleStage32(x, 0x0f000f00, 0x00f000f0, 4);
  if (shamt & 8)
    x = shuffleStage32(x, 0x00ff0000, 0x0000ff00, 8);

  return x;
}


static
uint64_t
shuffleStage64(uint64_t src, uint64_t maskL, uint64_t maskR, unsigned n)
{
  uint64_t x = src & ~(maskL | maskR);
  x |= ((src << n) & maskL) | ((src >> n) & maskR);
  return x;
}


static
uint64_t
shuffle64(uint64_t x, unsigned shamt)
{
  if (shamt & 16)
    x = shuffleStage64(x, 0x0000ffff00000000LL, 0x00000000ffff0000LL, 16);
  if (shamt & 8)
    x = shuffleStage64(x, 0x00ff000000ff0000LL, 0x0000ff000000ff00LL, 8);
  if (shamt & 4)
    x = shuffleStage64(x, 0x0f000f000f000f00LL, 0x00f000f000f000f0LL, 4);
  if (shamt & 2)
    x = shuffleStage64(x, 0x3030303030303030LL, 0x0c0c0c0c0c0c0c0cLL, 2);
  if (shamt & 1)
    x = shuffleStage64(x, 0x4444444444444444LL, 0x2222222222222222LL, 1);

  return x;
}


static
uint64_t
unshuffle64(uint64_t x, unsigned shamt)
{
  if (shamt & 1)
    x = shuffleStage64(x, 0x4444444444444444LL, 0x2222222222222222LL, 1);
  if (shamt & 2)
    x = shuffleStage64(x, 0x3030303030303030LL, 0x0c0c0c0c0c0c0c0cLL, 2);
  if (shamt & 4)
    x = shuffleStage64(x, 0x0f000f000f000f00LL, 0x00f000f000f000f0LL, 4);
  if (shamt & 8)
    x = shuffleStage64(x, 0x00ff000000ff0000LL, 0x0000ff000000ff00LL, 8);
  if (shamt & 16)
    x = shuffleStage64(x, 0x0000ffff00000000LL, 0x00000000ffff0000LL, 16);

  return x;
}


template <typename URV>
void
Hart<URV>::execSbset(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV mask = intRegs_.shiftMask();
  unsigned bitIx = intRegs_.read(di->op2()) & mask;

  URV value = intRegs_.read(di->op1()) | (URV(1) << bitIx);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSbclr(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV mask = intRegs_.shiftMask();
  unsigned bitIx = intRegs_.read(di->op2()) & mask;

  URV value = intRegs_.read(di->op1()) & ~(URV(1) << bitIx);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSbinv(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV mask = intRegs_.shiftMask();
  unsigned bitIx = intRegs_.read(di->op2()) & mask;

  URV value = intRegs_.read(di->op1()) ^ (URV(1) << bitIx);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSbext(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV mask = intRegs_.shiftMask();
  unsigned bitIx = intRegs_.read(di->op2()) & mask;

  URV value = (intRegs_.read(di->op1()) >> bitIx) & 1;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSbseti(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV bitIx = di->op2();
  if (not checkShiftImmediate(di, bitIx))
    return;

  URV value = intRegs_.read(di->op1()) | (URV(1) << bitIx);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSbclri(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV bitIx = di->op2();
  if (not checkShiftImmediate(di, bitIx))
    return;

  URV value = intRegs_.read(di->op1()) & ~(URV(1) << bitIx);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSbinvi(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV bitIx = di->op2();
  if (not checkShiftImmediate(di, bitIx))
    return;

  URV value = intRegs_.read(di->op1()) ^ (URV(1) << bitIx);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSbexti(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV bitIx = di->op2();
  if (not checkShiftImmediate(di, bitIx))
    return;

  URV value = (intRegs_.read(di->op1()) >> bitIx) & 1;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execBext(const DecodedInst* di)
{
  if (not isRvzbe())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV res = 0;
  for (unsigned i = 0, j = 0; i < mxlen_; ++i)
    if ((v2 >> i) & 1)
      {
        if ((v1 >> i) & 1)
          res |= URV(1) << j;
        ++j;
      }

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execBdep(const DecodedInst* di)
{
  if (not isRvzbe())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV res = 0;
  for (unsigned i = 0, j = 0; i < mxlen_; ++i)
    if ((v2 >> i) & 1)
      {
        if ((v1 >> j) & 1)
          res |= URV(1) << i;
        j++;
      }

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execBfp(const DecodedInst* di)
{
  if (not isRvzbf())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  unsigned off = (v2 >> 16) & intRegs_.shiftMask();
  unsigned len = (v2 >> 24) & 0xf;
  if (len == 0)
    len = 16;

  URV mask = (URV(1) << len) - 1;
  mask = (mask << off) | (mask >> (mxlen_ - off));
  URV data = (v2 << off) | (v2 >> (mxlen_ - off));

  URV res = (data & mask) | (v1 & ~mask);
  intRegs_.write(di->op0(), res);
}
    

template <typename URV>
void
Hart<URV>::execGorc(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  uint32_t shamt = intRegs_.read(di->op2()) & 0x1f;

  if constexpr (sizeof(URV) == 4)
    {
      if (shamt & 1)
        v1 |= ((v1 & 0xaaaaaaaa) >>  1) | ((v1 & 0x55555555) <<  1);
      if (shamt & 2)
        v1 |= ((v1 & 0xcccccccc) >>  2) | ((v1 & 0x33333333) <<  2);
      if (shamt & 4)
        v1 |= ((v1 & 0xf0f0f0f0) >>  4) | ((v1 & 0x0f0f0f0f) <<  4);
      if (shamt & 8)
      v1 |= ((v1 & 0xff00ff00) >>  8) | ((v1 & 0x00ff00ff) <<  8);
      if (shamt & 16)
        v1 |= ((v1 & 0xffff0000) >> 16) | ((v1 & 0x0000ffff) << 16);
    }
  else
    {
      if (shamt & 1)
        v1 |= ((v1 & 0xaaaaaaaaaaaaaaaa) >>  1) | ((v1 & 0x5555555555555555) <<  1);
      if (shamt & 2)
        v1 |= ((v1 & 0xcccccccccccccccc) >>  2) | ((v1 & 0x3333333333333333) <<  2);
      if (shamt & 4)
        v1 |= ((v1 & 0xf0f0f0f0f0f0f0f0) >>  4) | ((v1 & 0x0f0f0f0f0f0f0f0f) <<  4);
      if (shamt & 8)
        v1 |= ((v1 & 0xff00ff00ff00ff00) >>  8) | ((v1 & 0x00ff00ff00ff00ff) <<  8);
      if (shamt & 16)
        v1 |= ((v1 & 0xffff0000ffff0000) >> 16) | ((v1 & 0x0000ffff0000ffff) << 16);
      if (shamt & 32)
        v1 |= ((v1 & 0xffffffff00000000) >> 32) | ((v1 & 0x00000000ffffffff) << 32);
    }

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execClmul(const DecodedInst* di)
{
  if (not isRvzbc())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV x = 0;
  for (unsigned i = 0; i < mxlen_; ++i)
    if ((v2 >> i) & 1)
      x ^= v1 << i;

  intRegs_.write(di->op0(), x);
}


template <typename URV>
void
Hart<URV>::execClmulh(const DecodedInst* di)
{
  if (not isRvzbc())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV x = 0;
  for (unsigned i = 1; i < mxlen_; ++i)
    if ((v2 >> i) & 1)
      x ^= v1 >> (mxlen_ - i);

  intRegs_.write(di->op0(), x);
}


template <typename URV>
void
Hart<URV>::execShfl(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV val = 0;

  if constexpr (sizeof(URV) == 4)
    {
      unsigned shamt = v2 & 15;
      val = shuffle32(v1, shamt);
    }
  else
    {
      unsigned shamt = v2 & 31;
      val = shuffle64(v1, shamt);
    }

  intRegs_.write(di->op0(), val);
}


template <typename URV>
void
Hart<URV>::execShfli(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV amt = di->op2();
  URV val = 0;

  if constexpr (sizeof(URV) == 4)
    {
      if (amt > 15)
        {
          illegalInst(di);
          return;
        }
      val = shuffle32(v1, amt);
    }
  else
    {
      if (amt > 31)
        {
          illegalInst(di);
          return;
        }
      val = shuffle64(v1, amt);
    }

  intRegs_.write(di->op0(), val);
}


template <typename URV>
void
Hart<URV>::execUnshfl(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV val = 0;

  if constexpr (sizeof(URV) == 4)
    {
      unsigned shamt = v2 & 15;
      val = unshuffle32(v1, shamt);
    }
  else
    {
      unsigned shamt = v2 & 31;
      val = unshuffle64(v1, shamt);
    }

  intRegs_.write(di->op0(), val);
}


template <typename URV>
void
Hart<URV>::execUnshfli(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV amt = di->op2();
  URV val = 0;

  if constexpr (sizeof(URV) == 4)
    {
      if (amt > 15)
        {
          illegalInst(di);
          return;
        }
      val = unshuffle32(v1, amt);
    }
  else
    {
      if (amt > 31)
        {
          illegalInst(di);
          return;
        }
      val = unshuffle64(v1, amt);
    }

  intRegs_.write(di->op0(), val);
}


template <typename URV>
void
Hart<URV>::execClmulr(const DecodedInst* di)
{
  if (not isRvzbc())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV x = 0;
  for (unsigned i = 0; i < mxlen_; ++i)
    if ((v2 >> i) & 1)
      x ^= v1 >> (mxlen_ - i - 1);

  intRegs_.write(di->op0(), x);
}


template <typename URV>
void
Hart<URV>::execSh1add(const DecodedInst* di)
{
  if (not isRvzba())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV res = (v1 << 1) + v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execSh2add(const DecodedInst* di)
{
  if (not isRvzba())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV res = (v1 << 2) + v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execSh3add(const DecodedInst* di)
{
  if (not isRvzba())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV res = (v1 << 3) + v2;
  intRegs_.write(di->op0(), res);
}


static
Pmp::Mode
getModeFromPmpconfigByte(uint8_t byte)
{
  unsigned m = 0;

  if (byte & 1) m = Pmp::Read  | m;
  if (byte & 2) m = Pmp::Write | m;
  if (byte & 4) m = Pmp::Exec  | m;

  return Pmp::Mode(m);
}


template <typename URV>
void
Hart<URV>::execSh1addu_w(const DecodedInst* di)
{
  if (not isRv64() or not isRvzba())
    {
      illegalInst(di);
      return;
    }

  URV v1 = uint32_t(intRegs_.read(di->op1()));
  URV v2 = intRegs_.read(di->op2());

  URV res = (v1 << 1) + v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execSh2addu_w(const DecodedInst* di)
{
  if (not isRv64() or not isRvzba())
    {
      illegalInst(di);
      return;
    }

  URV v1 = uint32_t(intRegs_.read(di->op1()));
  URV v2 = intRegs_.read(di->op2());

  URV res = (v1 << 2) + v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execSh3addu_w(const DecodedInst* di)
{
  if (not isRv64() or not isRvzba())
    {
      illegalInst(di);
      return;
    }

  URV v1 = uint32_t(intRegs_.read(di->op1()));
  URV v2 = intRegs_.read(di->op2());

  URV res = (v1 << 3) + v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
static
URV
crc32(URV x, unsigned nbits)
{
  for (unsigned i = 0; i < nbits; ++i)
    x = (x >> 1) ^ (0xedb88320 & ~((x & 1) - 1));
  return x;
}


template <typename URV>
static
URV
crc32c(URV x, unsigned nbits)
{
  for (unsigned i = 0; i < nbits; ++i)
    x = (x >> 1) ^ (0x82F63B78  & ~((x & 1) - 1));
  return x;
}



template <typename URV>
void
Hart<URV>::execCrc32_b(const DecodedInst* di)
{
  if (not isRvzbr())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32(intRegs_.read(di->op1()), 8);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execCrc32_h(const DecodedInst* di)
{
  if (not isRvzbr())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32(intRegs_.read(di->op1()), 16);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execCrc32_w(const DecodedInst* di)
{
  if (not isRvzbr())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32(intRegs_.read(di->op1()), 32);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execCrc32_d(const DecodedInst* di)
{
  if (not isRvzbr() or not isRv64())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32(intRegs_.read(di->op1()), 64);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execCrc32c_b(const DecodedInst* di)
{
  if (not isRvzbr())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32c(intRegs_.read(di->op1()), 8);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execCrc32c_h(const DecodedInst* di)
{
  if (not isRvzbr())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32c(intRegs_.read(di->op1()), 16);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execCrc32c_w(const DecodedInst* di)
{
  if (not isRvzbr())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32c(intRegs_.read(di->op1()), 32);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execCrc32c_d(const DecodedInst* di)
{
  if (not isRvzbr() or not isRv64())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32c(intRegs_.read(di->op1()), 64);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execBmator(const DecodedInst* di)
{
  if (not isRvzbm() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint8_t u[8]; // rows of rs1
  uint8_t v[8]; // cols of rs2

  uint64_t rs1 = intRegs_.read(di->op1());
  uint64_t rs2 = intRegs_.read(di->op2());

  uint64_t rs2t = rs2;
  rs2t = shuffle64(rs2t, 31);
  rs2t = shuffle64(rs2t, 31);
  rs2t = shuffle64(rs2t, 31);

  for (int i = 0; i < 8; i++)
    {
      u[i] = rs1 >> (i*8);
      v[i] = rs2t >> (i*8);
    }

  uint64_t x = 0;
  for (int i = 0; i < 64; i++)
    {
      if ((u[i / 8] & v[i % 8]) != 0)
        x |= 1LL << i;
    }

  intRegs_.write(di->op0(), x);
}


template <typename URV>
void
Hart<URV>::updateMemoryProtection()
{
  pmpManager_.reset();

  const unsigned count = 16;
  unsigned impCount = 0;  // Count of implemented PMP registers

  // Process the pmp entries in reverse order (since they are supposed to
  // be checked in first to last priority). Apply memory protection to
  // the range defined by each entry allowing lower numbered entries to
  // over-ride higher numberd ones.
  unsigned pmpG = csRegs_.getPmpG();
  unsigned num = unsigned(CsrNumber::PMPADDR15);
  for (unsigned ix = 0; ix < count; ++ix, --num)
    {
      unsigned pmpIx = count - ix - 1;
      CsrNumber csrn = CsrNumber(num);
      if (not csRegs_.isImplemented(csrn))
        continue;
      impCount++;

      unsigned config = csRegs_.getPmpConfigByteFromPmpAddr(csrn);
      Pmp::Type type = Pmp::Type((config >> 3) & 3);
      bool lock = config & 0x80;

      Pmp::Mode mode = getModeFromPmpconfigByte(config);
      if (type == Pmp::Type::Off)
        continue;   // Entry is off.

      URV pmpVal = 0;
      if (not peekCsr(csrn, pmpVal))
        continue;  // Unimplemented PMPADDR reg. Should not happen.

      if (type == Pmp::Type::Tor)    // Top of range
        {
          uint64_t low = 0;
          if (pmpIx > 0)
            {
              URV prevVal = 0;
              CsrNumber lowerCsrn = CsrNumber(num - 1);
              peekCsr(lowerCsrn, prevVal);
              low = prevVal;
              low = (low >> pmpG) << pmpG;  // Clear least sig G bits.
              low = low << 2;
            }
              
          uint64_t high = pmpVal;
          high = (high >> pmpG) << pmpG;
          high = high << 2;
          if (low < high)
            pmpManager_.setMode(low, high - 1, type, mode, pmpIx, lock);
          continue;
        }

      uint64_t sizeM1 = 3;     // Size minus 1
      uint64_t napot = pmpVal;  // Naturally aligned power of 2.
      if (type == Pmp::Type::Napot)  // Naturally algined power of 2.
        {
          unsigned rzi = 0;  // Righmost-zero-bit index in pmpval.
          if (pmpVal == URV(-1))
            {
              // Handle special case where pmpVal is set to maximum value
              napot = 0;
              rzi = mxlen_;
            }
          else
            {
              rzi = __builtin_ctzl(~pmpVal); // rightmost-zero-bit ix.
              napot = (napot >> rzi) << rzi; // Clear bits below rightmost zero bit.
            }

          // Avoid overflow when computing 2 to the power 64 or
          // higher. This is incorrect but should work in practice
          // where the physical address space is 64-bit wide or less.
          if (rzi + 3 >= 64)
            sizeM1 = -1L;
          else
            sizeM1 = (uint64_t(1) << (rzi + 3)) - 1;
        }
      else
        assert(type == Pmp::Type::Na4);

      uint64_t low = napot;
      low = (low >> pmpG) << pmpG;
      low = low << 2;
      uint64_t high = low + sizeM1;
      pmpManager_.setMode(low, high, type, mode, pmpIx, lock);
    }

  pmpEnabled_ = impCount > 0;
}


template <typename URV>
void
Hart<URV>::updateAddressTranslation()
{
  if (not isRvs())
    return;

  URV value = 0;
  if (not peekCsr(CsrNumber::SATP, value))
    return;

  URV mode = 0, asid = 0, ppn = 0;
  if constexpr (sizeof(URV) == 4)
    {
      mode = value >> 31;
      asid = (value >> 22) & 0x1ff;
      ppn = value & 0x3fffff;  // Least sig 22 bits
    }
  else
    {
      mode = value >> 60;
      if ((mode >= 1 and mode <= 7) or mode >= 12)
        mode = 0;  // 1-7 and 12-15 are reserved in version 1.12 of sepc.
      asid = (value >> 44) & 0xffff;
      ppn = value & 0xfffffffffffll;  // Least sig 44 bits
    }

  virtMem_.setMode(VirtMem::Mode(mode));
  virtMem_.setAddressSpace(asid);
  virtMem_.setPageTableRootPage(ppn);
}


template <typename URV>
void
Hart<URV>::execBmatxor(const DecodedInst* di)
{
  if (not isRvzbm() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint8_t u[8]; // rows of rs1
  uint8_t v[8]; // cols of rs2

  uint64_t rs1 = intRegs_.read(di->op1());
  uint64_t rs2 = intRegs_.read(di->op2());

  uint64_t rs2t = rs2;
  rs2t = shuffle64(rs2t, 31);
  rs2t = shuffle64(rs2t, 31);
  rs2t = shuffle64(rs2t, 31);

  for (int i = 0; i < 8; i++)
    {
      u[i] = rs1 >> (i*8);
      v[i] = rs2t >> (i*8);
    }

  uint64_t x = 0;
  for (int i = 0; i < 64; i++)
    {
      if (__builtin_popcount(u[i / 8] & v[i % 8]) & 1)
        x ^= 1LL << i;
    }

  intRegs_.write(di->op0(), x);
}


template <typename URV>
void
Hart<URV>::execBmatflip(const DecodedInst* di)
{
  if (not isRvzbm() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint64_t rs1 = intRegs_.read(di->op1());
  rs1 = shuffle64(rs1, 31);
  rs1 = shuffle64(rs1, 31);
  rs1 = shuffle64(rs1, 31);

  intRegs_.write(di->op0(), rs1);
}


template <typename URV>
void
Hart<URV>::execCmov(const DecodedInst* di)
{
  if (not isRvzbt())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV v3 = intRegs_.read(di->op3());

  URV res = v2 ? v1 : v3;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execCmix(const DecodedInst* di)
{
  if (not isRvzbt())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV v3 = intRegs_.read(di->op3());

  URV res = (v1 & v2) | (v3 & ~v2);
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFsl(const DecodedInst* di)
{
  if (not isRvzbt())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV v3 = intRegs_.read(di->op3());

  unsigned shamt = v2 & (2*mxlen_ - 1);

  URV aa = v1, bb = v3;

  if (shamt >= mxlen_)
    {
      shamt -= mxlen_;
      aa = v3;
      bb = v1;
    }

  URV res = shamt ? (aa << shamt) | (bb >> (mxlen_ - shamt)) : aa;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFsr(const DecodedInst* di)
{
  if (not isRvzbt())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV v3 = intRegs_.read(di->op3());

  unsigned shamt = v2 & (2*mxlen_ - 1);

  URV aa = v1, bb = v3;

  if (shamt >= mxlen_)
    {
      shamt -= mxlen_;
      aa = v3;
      bb = v1;
    }

  URV res = shamt ? (aa >> shamt) | (bb << (mxlen_ - shamt)) : aa;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFsri(const DecodedInst* di)
{
  if (not isRvzbt())
    {
      illegalInst(di);
      return;
    }

  URV aa = intRegs_.read(di->op1());
  URV bb = intRegs_.read(di->op2());
  URV imm = intRegs_.read(di->op3());

  unsigned shamt = imm & (2*mxlen_ - 1);

  if (shamt >= mxlen_)
    {
      shamt -= mxlen_;
      std::swap(aa, bb);
    }

  URV res = shamt ? (aa >> shamt) | (bb << (mxlen_ - shamt)) : aa;
  intRegs_.write(di->op0(), res);
}


template class WdRiscv::Hart<uint32_t>;
template class WdRiscv::Hart<uint64_t>;
