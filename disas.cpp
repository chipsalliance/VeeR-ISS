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
#include <cstring>
#include "Hart.hpp"
#include "instforms.hpp"


using namespace WdRiscv;


static
std::string
roundingModeString(RoundingMode mode)
{
  switch (mode)
    {
    case RoundingMode::NearestEven: return "rne";
    case RoundingMode::Zero:        return "rtz";
    case RoundingMode::Down:        return "rdn";
    case RoundingMode::Up:          return "rup";
    case RoundingMode::NearestMax:  return "rmm";
    case RoundingMode::Invalid1:    return "inv1";
    case RoundingMode::Invalid2:    return "inv2";
    case RoundingMode::Dynamic:     return "dyn";
    default:                        return "inv";
    }
  return "inv";
}


/// Helper to disassemble method. Print on the given stream given
/// instruction which is of the form:  inst reg1, imm(reg2)
template <typename URV>
static
void
printLdSt(const Hart<URV>& hart, std::ostream& stream, const DecodedInst& di)
{
  unsigned rd = di.op0(), rs1 = di.op1();
  int32_t imm = di.op2As<int32_t>();

  const InstEntry* entry = di.instEntry();

  stream << std::left << std::setw(8) << entry->name() << ' ';

  const char* sign = imm < 0? "-" : "";
  if (imm < 0)
    imm = -imm;

  // Keep least sig 12 bits.
  imm = imm & 0xfff;

  stream << (entry->isFp() ? hart.fpRegName(rd) : hart.intRegName(rd))
	 << ", " << sign << "0x" << std::hex << imm
	 << "(" << hart.intRegName(rs1) << ")" << std::dec;
}


/// Helper to disassemble method. Print on the given stream the
/// disassembly of the given instruction.
template <typename URV>
static
void
printInst(const Hart<URV>& hart, std::ostream& out, const DecodedInst& di)
{
  const InstEntry* entry = di.instEntry();
  if (entry->isLoad() or entry->isStore())
    {
      printLdSt(hart, out, di);
      return;
    }

  out << std::left << std::setw(9) << entry->name();
  unsigned opCount = di.operandCount();

  const char* sep = "";

  for (unsigned i = 0; i < opCount; ++i)
    {
      out << sep; sep = ", ";

      auto type = di.ithOperandType(i);
      switch (type)
	{
	case OperandType::IntReg:
	  out << hart.intRegName(di.ithOperand(i));
	  break;
	case OperandType::FpReg:
	  out << hart.fpRegName(di.ithOperand(i));
	  break;
	case OperandType::VecReg:
	  out << "v" << di.ithOperand(i);
	  break;
	case OperandType::Imm:
	  out << di.ithOperandAsInt(i);
	  break;
	default:
	  out << "??";
	  break;
	}
    }

  if (entry->hasRoundingMode())
    out << sep << roundingModeString(di.roundingMode());
}


/// Helper to disassemble method. Print on the given stream given
/// instruction which is of the form: inst rd, rs2, rs1, rs3
template <typename URV>
static
void
printRdRs2Rs1Rs3(const Hart<URV>& hart, std::ostream& stream, const char* inst,
                 const DecodedInst& di)
{
  unsigned rd = di.op0(), rs1 = di.op1(), rs2 = di.op2(), rs3 = di.op3();

  // Print instruction in a 9 character field.
  stream << std::left << std::setw(9) << inst;

  stream << hart.intRegName(rd) << ", " << hart.intRegName(rs2)
         << ", " << hart.intRegName(rs1) << ", " << hart.intRegName(rs3);
}


/// Helper to disassemble method. Print on the given stream given
/// instruction which is of the form: inst rd, rs1, rs3, rs2
template <typename URV>
static
void
printRdRs1Rs3Rs2(const Hart<URV>& hart, std::ostream& stream, const char* inst,
                 const DecodedInst& di)
{
  unsigned rd = di.op0(), rs1 = di.op1(), rs2 = di.op2(), rs3 = di.op3();

  // Print instruction in a 9 character field.
  stream << std::left << std::setw(9) << inst;

  stream << hart.intRegName(rd) << ", " << hart.intRegName(rs1)
         << ", " << hart.intRegName(rs3) << ", " << hart.intRegName(rs2);
}


/// Helper to disassemble method. Print on the given stream given
/// instruction which is of the form: inst rd, rs1, rs3, immed
template <typename URV>
static
void
printRdRs1Rs3Imm(const Hart<URV>& hart, std::ostream& stream, const char* inst,
                   const DecodedInst& di)
{
  unsigned rd = di.op0(), rs1 = di.op1(), rs3 = di.op2();
  unsigned imm = di.op3();

  // Print instruction in a 9 character field.
  stream << std::left << std::setw(9) << inst;

  stream << hart.intRegName(rd) << ", " << hart.intRegName(rs1)
         << ", " << hart.intRegName(rs3) << ", 0x" << std::hex << imm
         << std::dec;
}


/// Helper to disassemble method. Print on the given stream given
/// instruction which is of the form: csrinst rd, csrn, rs1
template <typename URV>
static
void
printCsr(Hart<URV>& hart, std::ostream& stream,	 const DecodedInst& di)
{
  unsigned rd = di.op0(), csrn = di.op2();

  stream << std::left << std::setw(9) << di.instEntry()->name();
  stream << hart.intRegName(rd) << ", ";

  auto csr = hart.findCsr(CsrNumber(csrn));
  if (csr)
    stream << csr->getName();
  else
    stream << "illegal";

  if (di.ithOperandType(1) == OperandType::Imm)
    stream << ", 0x" << std::hex << di.op1() << std::dec;
  else
    stream << ", " << hart.intRegName(di.op1());
}


/// Helper to disassemble method. Print on the given stream given
/// instruction which is of the form: inst reg, imm where inst is a
/// compressed instruction.
template <typename URV>
static
void
printRegImm(const Hart<URV>& hart, std::ostream& stream, const char* inst,
	    unsigned rs1, int32_t imm)
{
  // Print instruction in a 8 character field.
  stream << std::left << std::setw(8) << inst << ' ';

  stream << hart.intRegName(rs1) << ", ";

  if (imm < 0)
    stream << "-0x" << std::hex << (-imm) << std::dec;
  else
    stream << "0x" << std::hex << imm << std::dec;
}


/// Helper to disassemble method. Print on the given stream given 3
/// operand branch instruction which is of the form: inst reg, reg,
/// imm where imm is a 12 bit constant.
template <typename URV>
static
void
printBranch3(const Hart<URV>& hart, std::ostream& stream,
	     const DecodedInst& di)
{
  unsigned rs1 = di.op0(), rs2 = di.op1();

  stream << std::left << std::setw(8) << di.instEntry()->name() << ' ';
  stream << hart.intRegName(rs1) << ", " << hart.intRegName(rs2) << ", . ";

  char sign = '+';
  int32_t imm = di.op2As<int32_t>();
  if (imm < 0)
    {
      sign = '-';
      imm = -imm;
    }
      
  stream << sign << " 0x" << std::hex << imm << std::dec;
}


/// Helper to disassemble method. Print on the given stream given
/// 2 operand  branch instruction which is of the form: inst reg, imm.
template <typename URV>
static
void
printBranch2(const Hart<URV>& hart, std::ostream& stream, const DecodedInst& di)
{
  unsigned rs1 = di.op0();
  int32_t imm = di.op2As<int32_t>();

  stream << std::left << std::setw(8) << di.instEntry()->name() << ' ';
  stream << hart.intRegName(rs1) << ", . ";

  char sign = '+';
  if (imm < 0)
    {
      sign = '-';
      imm = -imm;
    }
  stream << sign << " 0x" << std::hex << imm << std::dec;
}


/// Helper to disassemble method.
template <typename URV>
static
void
printAmo(const Hart<URV>& hart, std::ostream& stream, const DecodedInst& di)
{
  unsigned rd = di.op0(), rs1 = di.op1(), rs2 = di.op2();
  bool aq = di.isAtomicAcquire(), rl = di.isAtomicRelease();

  stream << di.instEntry()->name();

  if (aq)
    stream << ".aq";

  if (rl)
    stream << ".rl";

  stream << ' ' << hart.intRegName(rd) << ", " << hart.intRegName(rs2) << ", ("
	 << hart.intRegName(rs1) << ")";
}


/// Helper to disassemble method.
template <typename URV>
static
void
printLr(const Hart<URV>& hart, std::ostream& stream, const char* inst,
	const DecodedInst& di)
{
  unsigned rd = di.op0(), rs1 = di.op1();
  bool aq = di.isAtomicAcquire(), rl = di.isAtomicRelease();

  stream << inst;

  if (aq)
    stream << ".aq";

  if (rl)
    stream << ".rl";

  stream << ' ' << hart.intRegName(rd) << ", (" << hart.intRegName(rs1) << ")";
}


/// Helper to disassemble method.
template <typename URV>
static
void
printSc(const Hart<URV>& hart, std::ostream& stream, const char* inst,
	const DecodedInst& di)
{
  unsigned rd = di.op0(), rs1 = di.op1(), rs2 = di.op2();
  bool aq = di.isAtomicAcquire(), rl = di.isAtomicRelease();

  stream << inst;

  if (aq)
    stream << ".aq";

  if (rl)
    stream << ".rl";

  stream << ' ' << hart.intRegName(rd) << ", " << hart.intRegName(rs2)
	 << ", (" << hart.intRegName(rs1) << ")";
}


static
std::string
insertFieldCountInName(const std::string& name, unsigned count, unsigned n)
{
  std::string res = name.substr(0, n) + std::to_string(count) + name.substr(n);
  return res;
}


template <typename URV>
static
void
printVecInst(Hart<URV>& hart, std::ostream& out, const DecodedInst& di)
{
  uint32_t opcode7 = di.inst() & 0x7f;  // Least sig 7 bits
  InstId id = di.instEntry()->instId();

  if (opcode7 == 0x7 or opcode7 == 0x27)
    {  // Vector load store
      std::string name = di.instEntry()->name();
      if (id >= InstId::vlre8_v and id <= InstId::vlre1024_v)
	name = insertFieldCountInName(name, di.vecFieldCount(), 2);
      else if ((id >= InstId::vlsege8_v and id <= InstId::vssege1024_v) or
	       (id >= InstId::vlsege8ff_v and id <= InstId::vlsege1024ff_v))
	name = insertFieldCountInName(name, di.vecFieldCount(), 5);
      else if (id >= InstId::vlssege8_v and id <= InstId::vsssege1024_v)
	name = insertFieldCountInName(name, di.vecFieldCount(), 6);
      else if (id >= InstId::vluxsegei8_v and id <= InstId::vsoxsegei1024_v)
	name = insertFieldCountInName(name, di.vecFieldCount(), 7);
      out << name << " v" << di.op0();
      out << ", ("  << hart.intRegName(di.op1()) << ")";
      if (di.operandCount() == 3)
	{
	  if (di.ithOperandType(2) == OperandType::IntReg)
	    out << ", " << hart.intRegName(di.ithOperand(2));
	  else
	    out << ", v" << di.op2();
	}
      if (di.isMasked())
	out << ", v0.t";
      return;
    }

  if (id == InstId::vsetvli or id == InstId::vsetivli)
    {
      out << di.instEntry()->name() << ' ' << hart.intRegName(di.op0()) << ", ";
      if (id == InstId::vsetivli)
	out << di.op1();
      else
	out << hart.intRegName(di.op1());
      out << ", ";
      std::string mm = ((di.op2() >> 7) & 1) ? "ma" : "mu";
      std::string tt = ((di.op2() >> 6) & 1) ? "ta" : "tu";
      auto gm = VecRegs::to_string(GroupMultiplier(di.op2() & 7));
      auto ew = VecRegs::to_string(ElementWidth((di.op2() >> 3) & 7));
      out << ew << ',' << gm << ',' << tt << ',' << mm;
      return;
    }

  if (id == InstId::vsetvl)
    {
      out << "vsetvl " << hart.intRegName(di.op0()) << ", "
	  << hart.intRegName(di.op1()) << ", " << hart.intRegName(di.op2());
      return;
    }

  std::string name = di.instEntry()->name();
  if (id >= InstId::vmadc_vvm and id <= InstId::vmsbc_vxm and not di.isMasked())
    name = name.substr(0, name.size() - 1);
  out << name;

  const char* sep = " ";

  for (unsigned i = 0; i < di.operandCount(); ++i)
    {
      out << sep; sep = ", ";

      auto type = di.ithOperandType(i);
      switch (type)
	{
	case OperandType::IntReg:
	  out << hart.intRegName(di.ithOperand(i));
	  break;
	case OperandType::FpReg:
	  out << hart.fpRegName(di.ithOperand(i));
	  break;
	case OperandType::VecReg:
	  out << "v" << di.ithOperand(i);
	  break;
	case OperandType::Imm:
	  out << di.ithOperandAsInt(i);
	  break;
	default:
	  out << "??";
	  break;
	}
    }

  if (di.isMasked())
    {
      if ((id >= InstId::vadc_vvm and id <= InstId::vmsbc_vxm) or
	  (id >= InstId::vmerge_vvm and id <= InstId::vmerge_vim) or
	  (id == InstId::vfmerge_vfm))
	out << sep << "v0";
      else
	out << sep << "v0.t";
    }
}
	  

template <typename URV>
void
Hart<URV>::disassembleInst(uint32_t inst, std::ostream& stream)
{
  DecodedInst di;
  uint64_t physPc = pc_;
  decode(pc_, physPc, inst, di);
  disassembleInst(di, stream);
}


template <typename URV>
void
Hart<URV>::disassembleInst(uint32_t inst, std::string& str)
{
  str.clear();

  std::ostringstream oss;
  disassembleInst(inst, oss);
  str = oss.str();
}


template <typename URV>
void
Hart<URV>::disassembleInst(const DecodedInst& di, std::ostream& out)
{
  InstId id = di.instEntry()->instId();
  switch(id)
    {
    case InstId::illegal:
      out << "illegal";
      break;

    case InstId::lui:
      printRegImm(*this, out, "lui", di.op0(), di.op1As<int32_t>() >> 12);
      break;

    case InstId::auipc:
      out << "auipc    " << intRegName(di.op0())
	  << ", 0x" << std::hex << ((di.op1() >> 12) & 0xfffff) << std::dec;
      break;

    case InstId::jal:
      {
	if (di.op0() == 0)
	  out << "j        ";
	else
	  out << "jal      " << intRegName(di.op0()) << ", ";
	char sign = '+';
	int32_t imm = di.op1As<int32_t>();
	if (imm < 0) { sign = '-'; imm = -imm; }
	out << ". " << sign << " 0x" << std::hex << (imm & 0xfffff) << std::dec;
      }
      break;

    case InstId::jalr:
      printLdSt(*this, out, di);
      break;

    case InstId::beq:
    case InstId::bne:
    case InstId::blt:
    case InstId::bge:
    case InstId::bltu:
    case InstId::bgeu:
      printBranch3(*this, out, di);
      break;

    case InstId::csrrw:
    case InstId::csrrs:
    case InstId::csrrc:
    case InstId::csrrwi:
    case InstId::csrrsi:
    case InstId::csrrci:
      printCsr(*this, out, di);
      break;

    case InstId::lr_w:
      printLr(*this, out, "lr.w", di);
      break;

    case InstId::sc_w:
      printSc(*this, out, "sc.w", di);
      break;

    case InstId::lr_d:
      printLr(*this, out, "lr.d", di);
      break;

    case InstId::sc_d:
      printSc(*this, out, "sc.d", di);
      break;

    case InstId::c_addi4spn:
      printRegImm(*this, out, "c.addi4spn ", di.op0(), di.op2As<int32_t>() >> 2);
      break;

    case InstId::c_lq:
      out << "illegal";
      break;

    case InstId::c_sq:
      out << "illegal";
      break;

    case InstId::c_addi:
      if (di.op0() == 0)
	out << "c.nop";
      else
	printRegImm(*this, out, "c.addi", di.op0(), di.op2As<int32_t>());
      break;

    case InstId::c_jal:
      {
	out << "c.jal    . ";
	int32_t imm = di.op1As<int32_t>();
	char sign = '+';
	if (imm < 0) { sign = '-'; imm = -imm; }
	out << sign << " 0x" << std::hex << imm << std::dec;
      }
      break;

    case InstId::c_li:
      printRegImm(*this, out, "c.li", di.op0(), di.op2As<int32_t>());
      break;

    case InstId::c_addi16sp:
      {
	int32_t imm = di.op2As<int32_t>();
	out << "c.addi16sp ";
	if (imm < 0) { out << "-"; imm = -imm; }
	out << "0x" << std::hex << (imm >> 4) << std::dec;
      }
      break;

    case InstId::c_lui:
      printRegImm(*this, out, "c.lui", di.op0(), di.op1() >> 12);
      break;

    case InstId::c_srli:
      printRegImm(*this, out, "c.srli", di.op0(), di.op2As<int32_t>());
      break;

    case InstId::c_srli64:
      printRegImm(*this, out, "c.srli64", di.op0(), di.op2As<int32_t>());
      break;

    case InstId::c_srai:
      printRegImm(*this, out, "c.srai", di.op0(), di.op2As<int32_t>());
      break;

    case InstId::c_srai64:
      printRegImm(*this, out, "c.srai64", di.op0(), di.op2As<int32_t>());
      break;

    case InstId::c_andi:
      printRegImm(*this, out, "c.andi", di.op0(), di.op2As<int32_t>());
      break;

    case InstId::c_sub:
      out << "c.sub    " << intRegName(di.op0()) << ", " << intRegName(di.op2());
      break;

    case InstId::c_xor:
      out << "c.xor    " << intRegName(di.op0()) << ", " << intRegName(di.op2());
      break;

    case InstId::c_or:
      out << "c.or     " << intRegName(di.op0()) << ", " << intRegName(di.op2());
      break;

    case InstId::c_and:
      out << "c.and    " << intRegName(di.op0()) << ", " << intRegName(di.op2());
      break;

    case InstId::c_subw:
      out << "c.subw   " << intRegName(di.op0()) << ", " << intRegName(di.op2());
      break;

    case InstId::c_addw:
      out << "c.addw   " << intRegName(di.op0()) << ", " << intRegName(di.op2());
      break;

    case InstId::c_j:
      {
	out << "c.j      . ";
	int32_t imm = di.op1As<int32_t>();
	char sign = '+';
	if (imm < 0) { sign = '-'; imm = -imm; }
	out << sign << " 0x" << std::hex << imm << std::dec;
      }
      break;

    case InstId::c_beqz:
    case InstId::c_bnez:
      printBranch2(*this, out, di);
      break;

    case InstId::c_slli:
      out << "c.slli   " << intRegName(di.op0()) << ", " << di.op2();
      break;

    case InstId::c_slli64:
      out << "c.slli64 " << intRegName(di.op0()) << ", " << di.op2();
      break;

    case InstId::c_fldsp:
      out << "c.fldsp   " << fpRegName(di.op0()) << ", 0x" << std::hex
	  << di.op2As<int32_t>() << std::dec;
      break;

    case InstId::c_lwsp:
      out << "c.lwsp   " << intRegName(di.op0()) << ", 0x" << std::hex
	  << di.op2As<int32_t>() << std::dec;
      break;

    case InstId::c_flwsp:
      out << "c.flwsp   " << fpRegName(di.op0()) << ", 0x" << std::hex
	  << di.op2As<int32_t>() << std::dec;
      break;

    case InstId::c_ldsp:
      out << "c.ldsp   " << intRegName(di.op0()) << ", 0x" << std::hex
	  << di.op2As<int32_t>() << std::dec;
      break;

    case InstId::c_jr:
      out << "c.jr     " << intRegName(di.op1());
      break;

    case InstId::c_mv:
      out << "c.mv     " << intRegName(di.op0()) << ", " << intRegName(di.op2());
      break;

    case InstId::c_ebreak:
      out << "c.ebreak";
      break;

    case InstId::c_jalr:
      out << "c.jalr   " << intRegName(di.op1());
      break;

    case InstId::c_add:
      out << "c.add    " << intRegName(di.op0()) << ", " << intRegName(di.op2());
      break;

    case InstId::c_fsdsp:
      out << "c.fsdsp   " << fpRegName(di.op0())
	  << ", 0x" << std::hex << di.op2As<int32_t>() << std::dec;
      break;

    case InstId::c_swsp:
      out << "c.swsp   " << intRegName(di.op0()) << ", 0x"
	  << std::hex << di.op2As<int32_t>() << std::dec;
      break;

    case InstId::c_fswsp:
      out << "c.swsp   " << fpRegName(di.op0()) << ", 0x"
	  << std::hex << di.op2As<int32_t>() << std::dec;
      break;

    case InstId::c_addiw:
      printRegImm(*this, out, "c.addiw", di.op0(), di.op2As<int32_t>());
      break;

    case InstId::c_sdsp:
      out << "c.sdsp   " << intRegName(di.op0()) << ", 0x"
	  << std::hex << di.op2As<int32_t>() << std::dec;
      break;

    case InstId::cmov:
      printRdRs2Rs1Rs3(*this, out, "cmov", di);
      break;

    case InstId::cmix:
      printRdRs2Rs1Rs3(*this, out, "cmix", di);
      break;

    case InstId::fsl:
      printRdRs1Rs3Rs2(*this, out, "fsl", di);
      break;

    case InstId::fsr:
      printRdRs1Rs3Rs2(*this, out, "fsr", di);
      break;

    case InstId::fsri:
      printRdRs1Rs3Imm(*this, out, "fsri", di);
      break;

    case InstId::fslw:
      printRdRs1Rs3Rs2(*this, out, "fslw", di);
      break;

    case InstId::fsrw:
      printRdRs1Rs3Rs2(*this, out, "fsrw", di);
      break;

    case InstId::fsriw:
      printRdRs1Rs3Imm(*this, out, "fsriw", di);
      break;

    case InstId::load64:
      printLdSt(*this, out, di);
      break;

    case InstId::store64:
      printLdSt(*this, out, di);
      break;

    default:
      if (di.instEntry()->isAtomic())
	printAmo(*this, out, di);
      else if (di.instEntry()->isVector())
	printVecInst(*this, out, di);
      else
	printInst(*this, out, di);
    }
}


template <typename URV>
void
Hart<URV>::disassembleInst(const DecodedInst& di, std::string& str)
{
  str.clear();

  std::ostringstream oss;
  disassembleInst(di, oss);
  str = oss.str();
}


template class WdRiscv::Hart<uint32_t>;
template class WdRiscv::Hart<uint64_t>;
