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
#include <fstream>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include "Interactive.hpp"
#include "linenoise.hpp"

using namespace WdRiscv;


/// Return format string suitable for printing an integer of type URV
/// in hexadecimal form.
template <typename URV>
static
const char*
getHexForm()
{
  if (sizeof(URV) == 4)
    return "0x%08x";
  if (sizeof(URV) == 8)
    return "0x%016x";
  if (sizeof(URV) == 16)
    return "0x%032x";
  return "0x%x";
}


/// Convert the command line string numberStr to a number using
/// strotull and a base of zero (prefixes 0 and 0x are
/// honored). Return true on success and false on failure (string does
/// not represent a number). TYPE is an integer type (e.g
/// uint32_t). Option is the command line option associated with the
/// string and is used for diagnostic messages.
template <typename TYPE>
static
bool
parseCmdLineNumber(const std::string& option,
		   const std::string& numberStr,
		   TYPE& number)
{
  bool good = not numberStr.empty();

  if (good)
    {
      char* end = nullptr;
      uint64_t value = strtoull(numberStr.c_str(), &end, 0);
      number = static_cast<TYPE>(value);
      if (number != value)
	{
	  std::cerr << "parseCmdLineNumber: Number too large: " << numberStr
		    << '\n';
	  return false;
	}
      if (end and *end)
	good = false;  // Part of the string are non parseable.
    }

  if (not good)
    std::cerr << "Invalid command line " << option << " value: " << numberStr
	      << '\n';
  return good;
}


static
bool
parseCmdLineBool(const std::string& option, const std::string& str, bool& val)
{
  bool good = true;

  if (str == "0" or str == "false")
    val = false;
  else if (str == "1" or str == "true")
    val = true;
  else
    good = false;

  if (not good)
    std::cerr << "Invalid command line " << option << " value: " << str
	      << '\n';

  return good;
}


template <typename URV>
Interactive<URV>::Interactive(System<URV>& system)
  : system_(system)
{
}


template <typename URV>
bool
Interactive<URV>::untilCommand(Hart<URV>& hart, const std::string& line,
			       const std::vector<std::string>& tokens,
			       FILE* traceFile)
{
  if (tokens.size() != 2)
    {
      std::cerr << "Invalid until command: " << line << '\n';
      std::cerr << "Expecting: until address\n";
      return false;
    }

  size_t addr = 0;
  if (not parseCmdLineNumber("address", tokens.at(1), addr))
    return false;

  if (addr >= hart.memorySize())
    std::cerr << "Warning: Address outside memory range.\n";

  return hart.untilAddress(addr, traceFile);
}


template <typename URV>
bool
Interactive<URV>::stepCommand(Hart<URV>& hart, const std::string& /*line*/,
			      const std::vector<std::string>& tokens,
			      FILE* traceFile)
{
  if (not hart.isStarted())
    {
      // WD special.
      std::cerr << "Cannot step a non-started hart: Consider writing "
                << "the mhartstart CSR\n";
      return false;
    }

  if (tokens.size() == 1)
    {
      hart.singleStep(traceFile);
      hart.clearTraceData();
      return true;
    }

  uint64_t count = 0;
  if (not parseCmdLineNumber("instruction-count", tokens.at(1), count))
    return false;

  if (count == 0)
    return true;

  for (uint64_t i = 0; i < count; ++i)
    {
      hart.singleStep(traceFile);
      hart.clearTraceData();
    }

  return true;
}


template <typename URV>
void
Interactive<URV>::peekAllFpRegs(Hart<URV>& hart, std::ostream& out)
{
  for (unsigned i = 0; i < hart.fpRegCount(); ++i)
    {
      uint64_t val = 0;
      if (hart.peekFpReg(i, val))
	{
	  out << "f" << i << ": "
              << (boost::format("0x%016x") % val) << '\n';
	}
    }
}


template <typename URV>
void
Interactive<URV>::peekAllIntRegs(Hart<URV>& hart, std::ostream& out)
{
  bool abiNames = hart.abiNames();
  auto hexForm = getHexForm<URV>(); // Format string for printing a hex val

  for (unsigned i = 0; i < hart.intRegCount(); ++i)
    {
      std::string name;
      URV val = 0;
      if (hart.peekIntReg(i, val, name))
	{
	  std::string tag = name;
	  if (abiNames)
	    tag += "(" + std::to_string(i) + ")";
	  tag += ":";

          out << (boost::format("%-9s") % tag)
              << (boost::format(hexForm) % val) << '\n';
	}
    }
}


template <typename URV>
extern void
unpackMacoValue(URV value, URV mask, bool rv32, uint64_t& start, uint64_t& end,
                bool& idempotent, bool& cacheable);


template <typename URV>
void
Interactive<URV>::peekAllCsrs(Hart<URV>& hart, std::ostream& out)
{
  auto hexForm = getHexForm<URV>(); // Format string for printing a hex val

  out << (boost::format("%-23s") % "csr");
  if (sizeof(URV) == 4)
    out << (boost::format("%-10s %-10s %-10s %-10s\n") % "value" %
            "reset" % "mask" % "pokemask");
  else
    out << (boost::format("%-18s %-18s %-18s %-10s\n") % "value" %
            "reset" % "mask" % "pokemask");

  for (size_t i = 0; i <= size_t(CsrNumber::MAX_CSR_); ++i)
    {
      CsrNumber csr = CsrNumber(i);
      std::string name;
      URV val = 0;
      if (hart.peekCsr(csr, val, name))
	{
	  std::ostringstream oss;
	  oss << name << "(0x" << std::hex << i << "):"  << std::dec;

	  out << (boost::format("%-23s") % oss.str())
              << (boost::format(hexForm) % val);

	  URV reset = 0, writeMask = 0, pokeMask = 0;
	  if (hart.peekCsr(csr, val, reset, writeMask, pokeMask))
	    {
	      out << ' ' << (boost::format(hexForm) % reset);
	      out << ' ' << (boost::format(hexForm) % writeMask);
	      out << ' ' << (boost::format(hexForm) % pokeMask);
	    }
	  out << '\n';
	}
    }

  out << '\n';

  PrivilegeMode pm = hart.privilegeMode();
  out << "Privilege mode: ";
  switch(pm)
    {
    case PrivilegeMode::User:       out << "user\n";       break;
    case PrivilegeMode::Supervisor: out << "supervisor\n"; break;
    case PrivilegeMode::Reserved:   out << "reserved\n";   break;
    case PrivilegeMode::Machine:    out << "machine\n";    break;
    }

  out << '\n';

  out << "pmpaddr  type mode locked low                high\n";

  uint64_t low = 0, high = 0;
  Pmp::Type type = Pmp::Type::Off;
  Pmp::Mode mode = Pmp::Mode::None;
  bool locked = false;

  for (unsigned ix = 0; ix < 16; ++ix)
    {
      if (not hart.unpackMemoryProtection(ix, type, mode, locked, low, high))
        continue;

      std::string typeStr = Pmp::toString(type);
      std::string modeStr = Pmp::toString(mode);
      const char* lockStr = locked? "y" : "n";
      out << 
        (boost::format("%7d %5s %4s %6s 0x%016x 0x%016x") % ix % typeStr %
         modeStr % lockStr % low % high) << '\n';
    }

  bool headerPrinted = false;
  bool rv32 = sizeof(URV) == 4;

  for (unsigned ix = 0; ix < 16; ++ix)
    {
      std::string name = std::string("maco") + std::to_string(ix);
      auto maco = hart.findCsr(name);
      URV value = 0;
      if (maco and hart.peekCsr(maco->getNumber(), value))
        {
          if (not headerPrinted)
            {
              out << "\n";
              out << "maco io cacheable low                high\n";
              headerPrinted = true;
            }
          bool idempotent = false, cacheable = false;
          uint64_t low = 0, high = 0;
          unpackMacoValue(value, maco->getWriteMask(), rv32, low, high,
                          idempotent, cacheable);
          std::string ioStr = idempotent? "n" : "y";
          std::string cacheStr = cacheable? "y" : "n";
          out << 
            (boost::format("%4d %2s %9s 0x%016x 0x%016x") % ix % ioStr %
             cacheStr % low % high) << '\n';
        }
    }
}


template <typename URV>
void
Interactive<URV>::peekAllTriggers(Hart<URV>& hart, std::ostream& out)
{
  auto hexForm = getHexForm<URV>(); // Format string for printing a hex val

  out << (boost::format("%-12s") % "trigger");
  if (sizeof(URV) == 4)
    out << (boost::format("%-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n") %
            "value1" % "value2" % "value3" %
            "mask1" % "mask2" % "mask3" %
            "poke-mask1" % "poke-mask2"  % "poke-mask3");
  else
    out << (boost::format("%-18s %-18s %-18s %-18s %-18s %-18s %-18s %-18s %-18s\n") %
            "value1" % "value2" % "value3" %
            "mask1" % "mask2" % "mask3" %
            "poke-mask1" % "poke-mask2"  % "poke-mask3");


  // value/reset/write-mask/poke-mask
  URV tselVal = 0, tselReset, tselWm = 0, tselPm = 0;

  if (hart.peekCsr(CsrNumber::TSELECT, tselVal, tselReset, tselWm, tselPm))
    {
      URV maxTrigger = tselWm;
      for (URV trigger = 0; trigger <= maxTrigger; ++trigger)
	{
	  uint64_t v1(0), v2(0), v3(0), wm1(0), wm2(0), wm3(0);
	  uint64_t pm1(0), pm2(0), pm3(0);

	  if (hart.peekTrigger(trigger, v1, v2, v3, wm1, wm2, wm3,
			       pm1, pm2, pm3))
	    {
	      std::string name = "trigger" + std::to_string(trigger) + ":";
	      out << (boost::format("%-11s") % name);
	      out << ' ' << (boost::format(hexForm) % v1);
	      out << ' ' << (boost::format(hexForm) % v2);
	      out << ' ' << (boost::format(hexForm) % v3);
	      out << ' ' << (boost::format(hexForm) % wm1);
	      out << ' ' << (boost::format(hexForm) % wm2);
	      out << ' ' << (boost::format(hexForm) % wm3);
	      out << ' ' << (boost::format(hexForm) % pm1);
	      out << ' ' << (boost::format(hexForm) % pm2);
	      out << ' ' << (boost::format(hexForm) % pm3);
	      out << '\n';
	    }
	  else
	    break;
	}
    }
}


template <typename URV>
static
bool
peekMemory(Hart<URV>& hart, uint64_t addr0, uint64_t addr1, std::ostream& out)
{
  auto hexForm = getHexForm<URV>(); // Format string for printing a hex val

  uint32_t word = 0;
  bool usePma = false;

  for (uint64_t addr = addr0; addr <= addr1; addr += 4)
    {
      if (not hart.peekMemory(addr, word, usePma))
        {
          std::cerr << "Peek memory address out of bounds: 0x"
                    << std::hex << addr << std::dec << '\n';
          return false;
        }
      out << (boost::format(hexForm) % addr) << ": ";
      out << (boost::format("0x%08x") % word) << std::endl;
    }

  return true;
}


template <typename URV>
bool
Interactive<URV>::peekCommand(Hart<URV>& hart, const std::string& line,
			      const std::vector<std::string>& tokens,
                              std::ostream& out)
{
  if (tokens.size() < 2)
    {
      std::cerr << "Invalid peek command: " << line << '\n';
      std::cerr << "Expecting: peek <item> <addr>  or  peek pc  or  peek all\n";
      std::cerr << "  Item is one of r, f, c, t , pc, or m for integer, floating point,\n";
      std::cerr << "  CSR, trigger register, program counter, or memory location respectively\n";

      std::cerr << "  example:  peek r x3\n";
      std::cerr << "  example:  peek c mtval\n";
      std::cerr << "  example:  peek m 0x4096\n";
      std::cerr << "  example:  peek t 0\n";
      std::cerr << "  example:  peek pc\n";
      return false;
    }

  auto hexForm = getHexForm<URV>(); // Format string for printing a hex val
  URV val = 0;

  const std::string& resource = tokens.at(1);

  if (resource == "all")
    {
      out << "pc: " << (boost::format(hexForm) % hart.peekPc()) << '\n';
      out << "\n";

      peekAllIntRegs(hart, out);
      out << "\n";

      peekAllCsrs(hart, out);
      out << "\n";

      peekAllTriggers(hart, out);
      return true;
    }

  if (resource == "pc")
    {
      URV pc = hart.peekPc();
      out << (boost::format(hexForm) % pc) << std::endl;
      return true;
    }

  if (tokens.size() < 3)
    {
      std::cerr << "Invalid peek command: " << line << '\n';
      std::cerr << "Expecting: peek <resource> <address>\n";
      return false;
    }

  const std::string& addrStr = tokens.at(2);

  if (resource == "m")
    {
      uint64_t addr0 = 0;
      if (not parseCmdLineNumber("memory-address", addrStr, addr0))
	return false;

      uint64_t addr1 = addr0;
      if (tokens.size() >= 4)
	if (not parseCmdLineNumber("memory-address", tokens.at(3), addr1))
	  return false;

      if (tokens.size() >= 5)
        {
          std::ofstream out(tokens.at(4));
          if (not out)
            {
              std::cerr << "Failed to open " << tokens.at(4) << " for write of peek command output\n";
              return false;
            }
          return peekMemory(hart, addr0, addr1, out);
        }

      return peekMemory(hart, addr0, addr1, out);
    }

  if (resource == "r")
    {
      if (addrStr == "all")
	{
	  peekAllIntRegs(hart, out);
	  return true;
	}

      unsigned intReg = 0;
      if (not hart.findIntReg(addrStr, intReg))
	{
	  std::cerr << "No such integer register: " << addrStr << '\n';
	  return false;
	}
      if (hart.peekIntReg(intReg, val))
	{
	  out << (boost::format(hexForm) % val) << std::endl;
	  return true;
	}
      std::cerr << "Failed to read integer register: " << addrStr << '\n';
      return false;
    }

  if (resource == "f")
    {
      if (not hart.isRvf())
	{
	  std::cerr << "Floating point extension is no enabled\n";
	  return false;
	}

      if (addrStr == "all")
	{
	  peekAllFpRegs(hart, out);
	  return true;
	}

      unsigned fpReg = 0;
      if (not hart.findFpReg(addrStr, fpReg))
	{
	  std::cerr << "No such integer register: " << addrStr << '\n';
	  return false;
	}
      uint64_t fpVal = 0;
      if (hart.peekFpReg(fpReg, fpVal))
	{
	  out << (boost::format("0x%016x") % fpVal) << std::endl;
	  return true;
	}
      std::cerr << "Failed to read fp register: " << addrStr << '\n';
      return false;
    }

  if (resource == "c")
    {
      if (addrStr == "all")
	{
	  peekAllCsrs(hart, out);
	  return true;
	}

      auto csr = hart.findCsr(addrStr);
      if (not csr)
	{
	  std::cerr << "No such CSR: " << addrStr << '\n';
	  return false;
	}
      if (hart.peekCsr(csr->getNumber(), val))
	{
	  out << (boost::format(hexForm) % val) << std::endl;
	  return true;
	}
      std::cerr << "Failed to read CSR: " << addrStr << '\n';
      return false;
    }

  if (resource == "t")
    {
      if (addrStr == "all")
	{
	  peekAllTriggers(hart, out);
	  return true;
	}

      URV trigger = 0;
      if (not parseCmdLineNumber("trigger-number", addrStr, trigger))
	return false;
      uint64_t v1(0), v2(0), v3(0);
      if (hart.peekTrigger(trigger, v1, v2, v3))
	{
	  out << (boost::format(hexForm) % v1) << ' '
              << (boost::format(hexForm) % v2) << ' '
              << (boost::format(hexForm) % v3) << std::endl;
	  return true;
	}
      std::cerr << "Trigger number out of bounds: " << addrStr << '\n';
      return false;
    }

  std::cerr << "No such resource: " << resource
	    << " -- expecting r, m, c, t, or pc\n";
  return false;
}


template <typename URV>
bool
Interactive<URV>::pokeCommand(Hart<URV>& hart, const std::string& line,
			      const std::vector<std::string>& tokens)
{
  if (tokens.size() < 3)
    {
      std::cerr << "Invalid poke command: " << line << '\n';
      std::cerr << "  Expecting: poke pc <value>\n";
      std::cerr << "    or       poke <resource> <address> <value>\n";
      std::cerr << "    or       poke t <number> <value1> <value2> <value3>\n";
      std::cerr << "  where <resource> is one of r, f, c, t, pc or m\n";
      return false;
    }

  const std::string& resource = tokens.at(1);
  uint64_t value = 0;

  if (resource == "pc")
    {
      if (not parseCmdLineNumber("pc", tokens.at(2), value))
	return false;
      hart.pokePc(value);
      return true;
    }

  size_t count = tokens.size();
  if ((resource == "t" and count != 6) or (resource != "t" and count != 4))
    {
      std::cerr << "Invalid poke command: " << line << '\n';
      std::cerr << "  Expecting: poke <resource> <address> <value>\n";
      std::cerr << "    or       poke t <number> <value1> <value2> <value3>\n";
      std::cerr << "  where <resource> is one of r, f, c, t, pc, or m\n";
      return false;
    }

  const std::string& addrStr = tokens.at(2);
  const std::string& valueStr = tokens.at(3);

  if (not parseCmdLineNumber("poke", valueStr, value))
    return false;

  if (resource == "r")
    {
      unsigned intReg = 0;
      if (hart.findIntReg(addrStr, intReg))
	{
	  if (hart.pokeIntReg(intReg, value))
	    return true;
	  std::cerr << "Failed to write integer register " << addrStr << '\n';
	  return false;
	}

      std::cerr << "No such integer register " << addrStr << '\n';
      return false;
    }

  if (resource == "f")
    {
      unsigned fpReg = 0;
      if (hart.findFpReg(addrStr, fpReg))
	{
	  if (hart.pokeFpReg(fpReg, value))
	    return true;
	  std::cerr << "Failed to write FP register " << addrStr << '\n';
	  return false;
	}

      std::cerr << "No such FP register " << addrStr << '\n';
      return false;
    }

  if (resource == "c")
    {
      auto csr = hart.findCsr(addrStr);
      if (csr)
	{
	  if (hart.pokeCsr(csr->getNumber(), value))
	    return true;
	  std::cerr << "Failed to write CSR " << addrStr << '\n';
	  return false;
	}

      std::cerr << "No such CSR " << addrStr << '\n';
      return false;
    }

  if (resource == "t")
    {
      URV trigger = 0, v1 = 0, v2 = 0, v3 = 0;
      if (not parseCmdLineNumber("trigger", addrStr, trigger))
	return false;
      if (not parseCmdLineNumber("value1", tokens.at(3), v1))
	return false;
      if (not parseCmdLineNumber("value2", tokens.at(4), v2))
	return false;
      if (not parseCmdLineNumber("value3", tokens.at(5), v3))
	return false;
      if (hart.pokeTrigger(trigger, v1, v2, v3))
	return true;
      std::cerr << "Trigger out of bounds: " << addrStr << '\n';
      return false;
    }

  if (resource == "m")
    {
      size_t addr = 0;
      if (not parseCmdLineNumber("address", addrStr, addr))
	return false;
      bool usePma = false; // Ignore physicla memory attributes.
      uint32_t word = value; // Memory peek/poke in words.
      if (hart.pokeMemory(addr, word, usePma))
	return true;
      std::cerr << "Address out of bounds: " << addrStr << '\n';
      return false;
    }

  std::cerr << "No such resource: " << resource <<
    " -- expecting r, c, m or pc\n";
  return false;
}


template <typename URV>
bool
Interactive<URV>::disassCommand(Hart<URV>& hart, const std::string& line,
				const std::vector<std::string>& tokens)
{
  if (tokens.size() >= 2 and tokens.at(1) == "opcode")
    {
      for (size_t i = 2; i < tokens.size(); ++i)
	{
	  uint32_t code = 0;
	  if (not parseCmdLineNumber("opcode", tokens[i], code))
	    return false;
	  std::string str;
	  hart.disassembleInst(code, str);
	  std::cout << "  " << tokens[i] << ":  " << str << '\n';
	}
      return true;
    }

  auto hexForm = getHexForm<URV>(); // Format string for printing a hex val

  if (tokens.size() == 3 and (tokens.at(1) == "func" or
			      tokens.at(1) == "function"))
    {
      std::string item = tokens.at(2);
      std::string name;  // item (if a symbol) or function name containing item
      ElfSymbol symbol;
      if (hart.findElfSymbol(item, symbol))
	name = item;
      else
	{
	  // See if address falls in a function, then disassemble function.
	  URV addr = 0;
	  if (not parseCmdLineNumber("address", item, addr))
	    return false;

	  // Find function containing address.
	  hart.findElfFunction(addr, name, symbol);
	}

      if (name.empty())
	{
	  std::cerr << "Not a function or an address withing a function: " << item
		    << '\n';
	  return false;
	}

      std::cout << "disassemble function " << name << ":\n";

      size_t start = symbol.addr_, end = symbol.addr_ + symbol.size_;
      for (size_t addr = start; addr < end; )
	{
	  uint32_t inst = 0;
          bool usePma = false;
	  if (not hart.peekMemory(addr, inst, usePma))
	    {
	      std::cerr << "Address out of bounds: 0x" << std::hex << addr
			<< '\n' << std::dec;
	      return false;
	    }

	  unsigned instSize = instructionSize(inst);
	  if (instSize == 2)
	    inst = (inst << 16) >> 16; // Clear top 16 bits.

	  std::string str;
	  hart.disassembleInst(inst, str);
	  std::cout << "  " << (boost::format(hexForm) % addr) << ' '
		    << (boost::format(hexForm) % inst) << ' ' << str << '\n';

	  addr += instSize;
	}
      return true;
    }

  if (tokens.size() != 3)
    {
      std::cerr << "Invalid disass command: " << line << '\n';
      std::cerr << "Expecting: disass opcode <number> ...\n";
      std::cerr << "       or: disass function <name>\n";
      std::cerr << "       or: disass function <addr>\n";
      std::cerr << "       or: disass <addr1> <addr2>\n";
      return false;
    }

  URV addr1, addr2;
  if (not parseCmdLineNumber("address", tokens[1], addr1))
    return false;

  if (not parseCmdLineNumber("address", tokens[2], addr2))
    return false;

  for (URV addr = addr1; addr <= addr2; )
    {
      uint32_t inst = 0;
      bool usePma = false;
      if (not hart.peekMemory(addr, inst, usePma))
	{
	  std::cerr << "Address out of bounds: 0x" << std::hex << addr
		    << '\n' << std::dec;
	  return false;
	}

      unsigned instSize = instructionSize(inst);
      if (instSize == 2)
	inst = (inst << 16) >> 16; // Clear top 16 bits.

      std::string str;
      hart.disassembleInst(inst, str);
      std::cout << (boost::format(hexForm) % addr) << ' '
		<< (boost::format(hexForm) % inst) << ' '
		<< str << '\n';

      addr += instSize;
    }

  return true;
}


template <typename URV>
bool
Interactive<URV>::elfCommand(Hart<URV>& hart, const std::string& line,
			     const std::vector<std::string>& tokens)
{
  if (tokens.size() != 2)
    {
      std::cerr << "Invalid elf command: " << line << '\n';
      std::cerr << "Expecting: elf <file-name>\n";
      return false;
    }

  std::string filePath = tokens.at(1);

  size_t entryPoint = 0;

  if (not hart.loadElfFile(filePath, entryPoint))
    return false;

  hart.pokePc(URV(entryPoint));

  return true;
}


template <typename URV>
bool
Interactive<URV>::hexCommand(Hart<URV>& hart, const std::string& line,
			     const std::vector<std::string>& tokens)
{
  if (tokens.size() != 2)
    {
      std::cerr << "Invalid hex command: " << line << '\n';
      std::cerr << "Expecting: hex <file-name>\n";
      return false;
    }

  std::string fileName = tokens.at(1);

  if (not hart.loadHexFile(fileName))
    return false;

  return true;
}


template <typename URV>
bool
Interactive<URV>::resetCommand(Hart<URV>& hart, const std::string& /*line*/,
			       const std::vector<std::string>& tokens)
{
  if (tokens.size() == 1)
    {
      hart.reset(resetMemoryMappedRegs_);
      return true;
    }

  if (tokens.size() == 2)
    {
      URV resetPc = 0;
      if (not parseCmdLineNumber("reset-pc", tokens[1], resetPc))
	return false;

      hart.defineResetPc(resetPc);
      hart.reset(resetMemoryMappedRegs_);
      return true;
    }

  std::cerr << "Invalid reset command (extra arguments)\n";
  return false;
}


template <typename URV>
bool
Interactive<URV>::replayFileCommand(const std::string& line,
				    const std::vector<std::string>& tokens,
				    std::ifstream& stream)
{
  if (tokens.size() != 2)
    {
      std::cerr << "Invalid replay_file command: " << line << '\n';
      std::cerr << "Expecting: replay_file <file-name>\n";
      return false;
    }

  std::string fileName = tokens.at(1);

  stream.close();
  stream.open(fileName.c_str());
  if (not stream.good())
    {
      std::cerr << "Failed to open replay-file '" << fileName << "'\n";
      return false;
    }

  return true;
}


template <typename URV>
bool
Interactive<URV>::exceptionCommand(Hart<URV>& hart, const std::string& line,
				   const std::vector<std::string>& tokens)
{
  using std::cerr;

  bool bad = false;
  URV addr = 0;

  if (tokens.size() < 2)
    bad = true;
  else
    {
      const std::string& tag = tokens.at(1);
      if (tag == "inst")
	{
	  if (tokens.size() == 2)
	    hart.postInstAccessFault(0);
	  else if (tokens.size() == 3)
	    {
	      if (parseCmdLineNumber("exception inst offset", tokens.at(2), addr))
		hart.postInstAccessFault(addr);
	      else
		bad = true;
	    }
	  else
	    bad = true;
	}

      else if (tag == "data")
	{
	  if (tokens.size() == 2)
	    hart.postDataAccessFault(0, SecondaryCause::LOAD_ACC_DOUBLE_ECC);
	  else if (tokens.size() == 3)
	    {
	      if (parseCmdLineNumber("exception data offset", tokens.at(2), addr))
		hart.postDataAccessFault(addr, SecondaryCause::LOAD_ACC_DOUBLE_ECC);
	      else
		bad = true;
	    }
	  else
	    bad = true;
	}

      else if (tag == "store")
	{
	  bad = tokens.size() != 3;
	  if (not bad)
	    {
	      bad = not parseCmdLineNumber("exception store address",
					   tokens.at(2), addr);
	      if (not bad)
                {
                  unsigned matches = 0;
                  if (not hart.applyStoreException(addr, matches))
                    {
                      cerr << "Invalid exception store command: " << line << '\n';
                      if (matches == 0)
                        cerr << "  No pending store or invalid address\n";
                      else
                        cerr << "  Multiple matching addresses (unsupported)\n";
                      return false;
                    }
                }
            }
        }

      else if (tag == "load")
	{
	  bad = tokens.size() != 4;
	  if (not bad)
	    {
	      bad = not parseCmdLineNumber("exception load address",
					   tokens.at(2), addr);
	      unsigned tag = 0;
	      bad = bad or not parseCmdLineNumber("exception load tag",
						  tokens.at(3), tag);
	      if (not bad)
                {
                  unsigned matches = 0;
                  if (not hart.applyLoadException(addr, tag, matches))
                    {
                      cerr << "Invalid exception load command: " << line << '\n';
                      if (matches == 0)
                        cerr << "  No pending load or invalid address/tag\n";
                      else
                        cerr << "  Multiple matching tags\n";
                      return false;
                    }
                }
	    }
	}

      else if (tag == "precise_load")
        {
	  if (tokens.size() == 2)
	    hart.postDataAccessFault(0, SecondaryCause::LOAD_ACC_PRECISE);
	  else if (tokens.size() == 3)
	    {
	      if (parseCmdLineNumber("exception precise_load offset", tokens.at(2), addr))
		hart.postDataAccessFault(addr, SecondaryCause::LOAD_ACC_PRECISE);
	      else
		bad = true;
	    }
	  else
	    bad = true;
        }

      else if (tag == "precise_store")
        {
	  if (tokens.size() == 2)
            hart.postDataAccessFault(0, SecondaryCause::STORE_ACC_PRECISE);
	  else if (tokens.size() == 3)
	    {
	      if (parseCmdLineNumber("exception precise_store offset", tokens.at(2), addr))
		hart.postDataAccessFault(addr, SecondaryCause::STORE_ACC_PRECISE);
	      else
		bad = true;
	    }
	  else
	    bad = true;
        }

      else if (tag == "nmi")
	{
	  bad = tokens.size() != 3;
	  if (not bad)
	    {
	      bad = not parseCmdLineNumber("nmi", tokens.at(2), addr);
              if (hart.isNmiEnabled())
                hart.setPendingNmi(NmiCause(addr));
	    }
	}

      else if (tag == "memory_data")
	{
	  if (parseCmdLineNumber("memory_data", tokens.at(2), addr))
	    {
	      return true;
	    }
	}

      else if (tag == "memory_inst")
	{
	  if (parseCmdLineNumber("memory_inst", tokens.at(2), addr))
	    {
	      return true;
	    }
	}

      else
	bad = true;
    }

  if (bad)
    {
      std::cerr << "Invalid exception command: " << line << '\n';
      std::cerr << "  Expecting: exception inst [<offset>]\n";
      std::cerr << "   or:       exception data [<offset>]\n";
      std::cerr << "   or:       exception load <address> <tag>\n";
      std::cerr << "   or:       exception store <address>\n";
      std::cerr << "   or:       exception nmi <cause>\n";
      return false;
    }

  return true;
}


template <typename URV>
bool
Interactive<URV>::loadFinishedCommand(Hart<URV>& hart, const std::string& line,
				      const std::vector<std::string>& tokens)
{
  if (tokens.size() != 3)
    {
      std::cerr << "Invalid load_finished command: " << line << '\n';
      std::cerr << "  Expecting: load_finished address tag\n";
      return false;
    }

  URV addr = 0;
  if (not parseCmdLineNumber("address", tokens.at(1), addr))
    return false;

  unsigned tag = 0;
  if (not parseCmdLineNumber("tag", tokens.at(2), tag))
      return false;

  unsigned matchCount = 0;
  hart.applyLoadFinished(addr, tag, matchCount);

  return true;
}


template <typename URV>
bool
Interactive<URV>::dumpMemoryCommand(const std::string& line,
                                    const std::vector<std::string>& tokens)
{
  if (tokens.size() != 2)
    {
      std::cerr << "Invalid dump_memory command: " << line << '\n';
      std::cerr << "  Expecting: dump_memory path\n";
      return false;
    }

  return system_.writeAccessedMemory(tokens[1]);
}


/// If tokens contain a string of the form hart=<id> then remove that
/// token from tokens and set hartId to <id> returning true. Return
/// false if no hart=<id> token is found or if there is an error (<id>
/// is not an integer value) in which case error is set to true.
static
bool
getCommandHartId(std::vector<std::string>& tokens, unsigned& hartId,
		 bool& error)
{
  error = false;
  if (tokens.empty())
    return false;

  bool hasHart = false;

  // Remaining tokens after removal of hart=<id> tokens.
  std::vector<std::string> rest;

  for (const auto& token : tokens)
    {
      if (boost::starts_with(token, "hart="))
	{
	  std::string value = token.substr(strlen("hart="));
	  try
	    {
	      hartId = boost::lexical_cast<unsigned>(value);
	      hasHart = true;
	    }
	  catch(...)
	    {
	      std::cerr << "Bad hart id: " << value << '\n';
	      error = true;
	      return false;
	    }
	}
      else
	rest.push_back(token);
    }

  tokens = rest;
  return hasHart;
}


/// Interactive "help" command.
static
void
printInteractiveHelp()
{
  using std::cout;
  cout << "The argument hart=<id> may be used with any command.\n";
  cout << "help [<command>]\n";
  cout << "  Print help for given command or for all commands if no command given.\n\n";
  cout << "run\n";
  cout << "  Run till interrupted.\n\n";
  cout << "until <address>\n";
  cout << "  Run until address or interrupted.\n\n";
  cout << "step [<n>]\n";
  cout << "  Execute n instructions (1 if n is missing).\n\n";
  cout << "peek <res> <addr>\n";
  cout << "  Print value of resource res (one of r, f, c, m) and address addr.\n";
  cout << "  For memory (m) up to 2 addresses may be provided to define a range\n";
  cout << "  of memory locations to be printed; also, an optional filename after\n";
  cout << "  the two addresses writes the command output to that file.\n";
  cout << "  examples: peek r x1   peek c mtval   peek m 0x4096\n";
  cout << "            peek m 0x10 0x40 out\n\n";
  cout << "peek pc\n";
  cout << "  Print value of the program counter.\n\n";
  cout << "peek all\n";
  cout << "  Print value of all non-memory resources\n\n";
  cout << "poke res addr value\n";
  cout << "  Set value of resource res (one of r, c or m) and address addr\n";
  cout << "  Examples: poke r x1 0xff  poke c 0x4096 0xabcd\n\n";
  cout << "disass opcode <code> <code> ...\n";
  cout << "  Disassemble opcodes. Example: disass opcode 0x3b 0x8082\n\n";
  cout << "disass function <name>\n";
  cout << "  Disassemble function with given name. Example: disas func main\n\n";
  cout << "disass <addr1> <addr2>>\n";
  cout << "  Disassemble memory locations between addr1 and addr2.\n\n";
  cout << "elf file\n";
  cout << "  Load elf file into simulated memory.\n\n";
  cout << "hex file\n";
  cout << "  Load hex file into simulated memory.\n\n";
  cout << "replay_file file\n";
  cout << "  Open command file for replay.\n\n";
  cout << "replay n\n";
  cout << "  Execute the next n commands in the replay file or all the\n";
  cout << "  remaining commands if n is missing.\n\n";
  cout << "replay step n\n";
  cout << "  Execute consecutive commands from the replay file until n\n";
  cout << "  step commands are executed or the file is exhausted\n\n";
  cout << "reset [<reset_pc>]\n";
  cout << "  Reset hart.  If reset_pc is given, then change the reset program\n";
  cout << "  counter to the given reset_pc before resetting the hart.\n\n";
  cout << "symbols\n";
  cout << "  List all the symbols in the loaded ELF file(s).\n\n";
  cout << "pagetable\n";
  cout << "  Print the entries of the address tanslation table.\n\n";
  cout << "exception inst [<offset>]\n";
  cout << "  Take an instruction access fault on the subsequent step command. Given\n";
  cout << "  offset (defaults to zero) is added to the instruction PC to form the address\n";
  cout << "  responsible for the fault (that address is placed in the mtval CSR).\n\n";
  cout << "exception data [<offset>]\n";
  cout << "  Take a data access fault on the subsequent load/store instruction executed\n";
  cout << "  by a step command. The offset value is currently not used.\n\n";
  cout << "quit\n";
  cout << "  Terminate the simulator\n\n";
}


template <typename URV>
void
Interactive<URV>::helpCommand(const std::vector<std::string>& tokens)
{
  using std::cout;

  if (tokens.size() <= 1)
    {
      printInteractiveHelp();
      return;
    }

  auto& tag = tokens.at(1);
  if (tag == "help")
    {
      cout << "help [<command>]\n"
	   << "  Print information about interactive commands. If a command\n"
	   << "  argument is given, print info about that command.\n";
      return;
    }

  if (tag == "run")
    {
      cout << "run\n"
	   << "  Run the target program until it exits (in newlib emulation mode),\n"
	   << "  it writes into the \"tohost\" location, or the user interrupts\n"
	   << "  it by pressing control-c on the keyboard.\n";
      return;
    }

  if (tag == "until")
    {
      cout << "until <address>\n"
	   << "  Same as run but the target program will also stop when the\n"
	   << "  instruction at the given address is reached (but before it is\n"
	   << "  executed).\n";
      return;
    }

  if (tag == "step")
    {
      cout << "step [<n>]\n"
	   << "  Execute a single instruction. If an integer argument <n> is\n"
	   << "  given, then execute up to n instructions or until a stop\n"
	   << "  condition (see run command) is encountered\n";
      return;
    }

  if (tag == "peek")
    {
      cout << "peek <res> <addr>\n"
	   << "peek pc\n"
	   << "  Show contents of given resource having given address. Possible\n"
	   << "  resources are r, f, c, or m for integer, floating-point,\n"
	   << "  control-and-status register or for memory respectively.\n"
	   << "  Addr stands for a register number, register name or memory\n"
	   << "  address. If resource is memory (m), then an additional address\n"
	   << "  may be provided to define a range of memory locations to be\n"
	   << "  display and an optional filename after 2nd address may be\n"
           << "  provided to write memory contents to a file.  Examples\n"
	   << "    peek pc\n"
	   << "    peek r t0\n"
	   << "    peek r x12\n"
	   << "    peek c mtval\n"
	   << "    peek m 0x80000000\n"
	   << "    peek m 0x80000000 0x80000010\n"
      	   << "    peek m 0x80000000 0x80000010 out\n";
      return;
    }

  if (tag == "poke")
    {
      cout << "poke <res> <addr> <value>\n"
	   << "poke pc <value>\n"
	   << "  Set the contents of given resource having given address to the\n"
	   << "  given value. Possible resources are r, f, c, or m for integer,\n"
	   << "  floating-point, control-and-status register or for memory\n"
	   << "  respectively. Addr stands for a register number, register name\n"
	   << "  or memory address.  Examples:\n"
	   << "    poke r t0 0\n"
	   << "    poke r x12 0x44\n"
	   << "    poke c mtval 0xff\n"
	   << "    poke m 0x80000000 0xabdcffff\n";
      return;
    }

  if (tag == "disas" or tag == "disass")
    {
      cout << "disass opcode <op0> <op1> ...\n"
	   << "disass func <address>\n"
	   << "disass <addr1> <addr2>\n"
	   << "  The first form will disassemble the given opcodes.\n"
	   << "  The second form will disassemble the instructions of the\n"
	   << "  function containing the given address.\n"
	   << "  The third form will disassemble the memory contents between\n"
	   << "  addresses addr1 and addr2 inclusive.\n";
      return;
    }

  if (tag == "elf")
    {
      cout << "elf <file> ...\n"
	   << "  Load into memory the contents of the given ELF file.\n"
	   << "  Set the program counter to the value of the ELF file entry point.\n"
	   << "  If the file contains the symbol \"tohost\" then subsequent writes\n"
	   << "  to the corresponding address will stop the simulation.\n";
      return;
    }

  if (tag == "replay_file")
    {
      cout << "replay_file <file> ...\n"
	   << "  Define the input replay file to serve as input for the replay\n"
	   << "  command. The user would typically load the commands of a session\n"
	   << "  and replays them in a subsequent session.\n";
      return;
    }

  if (tag == "replay")
    {
      cout << "replay [step] [<n>]\n"
	   << "  Without any arguments, replay all remaining commands in the\n"
	   << "  replay file (defined by the replay_file command).\n"
	   << "  With the keyword step, key-in on step commands in the replay\n"
	   << "  file. With an integer number n, replay n commands (or n step\n"
	   << "  commands if step keyword is present).\n";
      return;
    }

  if (tag == "reset")
    {
      cout << "reset [<reset_pc>]\n"
	   << "  Reset simulated processor. If reset_pc is given, then change\n"
	   << "  reset program counter to the given reset_pc before resetting\n"
	   << "  the processor.\n";
      return;
    }

  if (tag == "reset")
    {
      cout << "quit\n"
	   << "  Terminate the simulator.\n";
      return;
    }

  std::cerr << "No such command: " << tag	<< '\n';
}


/// Command line interpreter: Execute a command line.
template <typename URV>
bool
Interactive<URV>::executeLine(unsigned& currentHartId,
			      const std::string& inLine, FILE* traceFile,
			      FILE* commandLog,
			      std::ifstream& replayStream, bool& done)
{
  // Remove comments (anything starting with #).
  std::string line = inLine;
  auto sharpIx = line.find_first_of('#');
  if (sharpIx != std::string::npos)
    line = line.substr(0, sharpIx);

  // Remove leading/trailing white space
  boost::algorithm::trim_if(line, boost::is_any_of(" \t"));

  if (line.empty())
    return true;

  // Break line into tokens.
  std::vector<std::string> tokens;
  boost::split(tokens, line, boost::is_any_of(" \t"),
	       boost::token_compress_on);
  if (tokens.empty())
    return true;

  std::string outLine;   // Line to print on command log.

  // Recover hart id (if any) removing hart=<id> token from tokens.
  unsigned hartId = 0;
  bool error = false;
  bool hasHart = getCommandHartId(tokens, hartId, error);
  if (error)
    return false;

  if (hasHart)
    {
      outLine = line;
      currentHartId = hartId;
    }
  else
    {
      hartId = currentHartId;
      outLine = std::string("hart=") + std::to_string(hartId) + " " + line;
    }

  const std::string& command = tokens.front();
  if (command == "q" or command == "quit")
    {
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      done = true;
      return true;
    }

  auto hartPtr = system_.findHartByHartId(hartId);
  if (not hartPtr)
    {
      std::cerr << "Hart id out of bounds: " << hartId << '\n';
      return false;
    }

  Hart<URV>& hart = *hartPtr;

  if (not hart.isStarted())
    {
      if (command != "peek" and command != "poke" and command != "reset")
        {
          std::cerr << "Error: Command " << command << " received for a "
                    << "non-started hart.\n";
          return false;
        }
    }

  // After the first step/run/until command, a reset command will reset
  // the memory mapped registers.
  if (command == "step" or command == "run" or command == "until")
    resetMemoryMappedRegs_ = true;

  if (command == "run")
    {
      bool success = hart.run(traceFile);
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      return success;
    }

  if (command == "u" or command == "until")
    {
      if (not untilCommand(hart, line, tokens, traceFile))
	return false;
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      return true;
    }

  if (command == "s" or command == "step")
    {
      if (hart.inDebugMode() and not hart.inDebugStepMode())
	{
	  std::cerr << "Error: Single step while in debug-halt mode\n";
	  return false;
	}
      if (not stepCommand(hart, line, tokens, traceFile))
	return false;
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      return true;
    }

  if (command == "peek")
    {
      if (not peekCommand(hart, line, tokens, std::cout))
	return false;
       if (commandLog)
	 fprintf(commandLog, "%s\n", outLine.c_str());
       return true;
    }

  if (command == "poke")
    {
      if (not pokeCommand(hart, line, tokens))
	return false;
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      return true;
    }

  if (command == "d" or command == "disas" or command == "disass")
    {
      if (not disassCommand(hart, line, tokens))
	return false;
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      return true;
    }

  if (command == "elf")
    {
      if (not elfCommand(hart, line, tokens))
	return false;
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      return true;
    }

  if (command == "hex")
    {
      if (not hexCommand(hart, line, tokens))
	return false;
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      return true;
    }

  if (command == "reset")
    {
      if (not resetCommand(hart, line, tokens))
	return false;
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      return true;
    }

  if (command == "exception")
    {
      if (not exceptionCommand(hart, line, tokens))
	return false;
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      return true;
    }

  if (command == "enter_debug")
    {
      bool force = false;
      if (tokens.size() == 2)
        if (not parseCmdLineBool("force", tokens[1], force))
          return false;
      hart.enterDebugMode(hart.peekPc(), force);
      if (commandLog)
	fprintf(commandLog, "%s %s\n", outLine.c_str(), force? "true" : "false");
      return true;
    }

  if (command == "exit_debug")
    {
      hart.exitDebugMode();
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      return true;
    }

  if (command == "load_finished")
    {
      if (not loadFinishedCommand(hart, line, tokens))
	return false;
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      return true;
    }

  if (command == "cancel_div")
    {
      if (not hart.cancelLastDiv())
        std::cerr << "Warning: Unexpected cancel_div\n";
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      return true;
    }

  if (command == "cancel_lr")
    {
      hart.cancelLr();
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      return true;
    }

  if (command == "replay_file")
    {
      if (not replayFileCommand(line, tokens, replayStream))
	return false;
      return true;
    }

  if (command == "replay")
    {
      if (not replayStream.is_open())
	{
	  std::cerr << "No replay file defined. Use the replay_file to define one\n";
	  return false;
	}
      bool replayDone = false;
      if (not replayCommand(currentHartId, line, tokens, traceFile, commandLog,
			    replayStream, replayDone))
	return false;
      return true;
    }

  if (command == "symbols")
    {
      hart.printElfSymbols(std::cout);
      return true;
    }

  if (command == "pagetable")
    {
      hart.printPageTable(std::cout);
      return true;
    }

  if (command == "dump_memory")
    {
      if (not dumpMemoryCommand(line, tokens))
        return false;
      if (commandLog)
	fprintf(commandLog, "%s\n", outLine.c_str());
      return true;
    }

  if (command == "h" or command == "?" or command == "help")
    {
      helpCommand(tokens);
      return true;
    }

  std::cerr << "No such command: " << line << '\n';
  return false;
}


/// Interactive "replay" command.
template <typename URV>
bool
Interactive<URV>::replayCommand(unsigned& currentHartId,
				const std::string& line,
				const std::vector<std::string>& tokens,
				FILE* traceFile, FILE* commandLog,
				std::ifstream& replayStream, bool& done)
{
  std::string replayLine;
  uint64_t maxCount = ~uint64_t(0);  // Unlimited

  if (tokens.size() <= 2)    // Either replay or replay n.
    {
      if (tokens.size() == 2)
	if (not parseCmdLineNumber("command-count", tokens.at(1), maxCount))
	  return false;

      uint64_t count = 0;
      while (count < maxCount  and  not done  and
	     std::getline(replayStream, replayLine))
	{
	  if (not executeLine(currentHartId, replayLine, traceFile,
			      commandLog, replayStream, done))
	    return false;
	  count++;
	}
      return true;
    }

  if (tokens.size() == 3)
    {
      if (tokens.at(1) != "step")
	{
	  std::cerr << "Invalid command: " << line << '\n';
	  std::cerr << "Expecting: replay <step> <count>\n";
	  return false;
	}

      if (not parseCmdLineNumber("step-count", tokens.at(2), maxCount))
	return false;
      
      uint64_t count = 0;
      while (count < maxCount  and  not done   and
	     std::getline(replayStream, replayLine))
	{
	  if (not executeLine(currentHartId, replayLine, traceFile,
			      commandLog, replayStream, done))
	    return false;

	  std::vector<std::string> tokens;
	  boost::split(tokens, replayLine, boost::is_any_of(" \t"),
		       boost::token_compress_on);
	  if (tokens.size() > 0 and tokens.at(0) == "step")
	    count++;
	  else if (tokens.size() > 1 and tokens.at(1) == "step")
	    count++;
	}

      return true;
    }

  std::cerr << "Invalid command: " << line << '\n';
  std::cerr << "Expecting: replay, replay <count>, or replay step <count>\n";
  return false;    
}


template <typename URV>
bool
Interactive<URV>::interact(FILE* traceFile, FILE* commandLog)
{
  linenoise::SetHistoryMaxLen(1024);

  uint64_t errors = 0;
  unsigned currentHartId = 0;
  std::string replayFile;
  std::ifstream replayStream;

  const char* prompt = isatty(0) ? "whisper> " : "";

  auto hartPtr = system_.ithHart(0);
  if (hartPtr)
    {
      URV value = 0;
      if (hartPtr->peekCsr(CsrNumber::MHARTID, value))
        currentHartId = value;
    }

  bool done = false;
  while (not done)
    {
      errno = 0;
      std::string line = linenoise::Readline(prompt);

      if (line.empty())
	{
	  if (std::cin.eof())
	    return true;
	  continue;
	}

      linenoise::AddHistory(line.c_str());

      if (not executeLine(currentHartId, line, traceFile, commandLog,
			  replayStream, done))
	errors++;
    }

  return errors == 0;
}


template class WdRiscv::Interactive<uint32_t>;
template class WdRiscv::Interactive<uint64_t>;
