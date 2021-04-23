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

#include <cstdio>
#include <csignal>
#include <iostream>
#include <sstream>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <sys/socket.h>
#include <netinet/in.h>
#include "Hart.hpp"


#ifdef __MINGW64__
#define SIGTRAP 5
#endif

static
int
putDebugChar(char c, int fd)
{
  return write(fd, &c, sizeof(uint8_t));
}


static
uint8_t
getDebugChar(int fd)
{
  uint8_t res;
  if (read(fd, &res, sizeof(res)) == sizeof(res))
    return res;
  return uint8_t(-1); // TODO throw exception?
}


static
int
hexCharToInt(unsigned char c)
{
  if (c >= 'a' and c <= 'f')
    return 10 + c - 'a';
  if (c >= 'A' and c <= 'F')
    return 10 + c - 'A';
  if (c >= '0' and c <= '9')
    return c - '0';
  return -1;
}


static
bool
hexCharToInt(char c, unsigned& value)
{
  if (c >= 'a' and c <= 'f')
    {
      value = 10 + c - 'a';
      return true;
    }

  if (c >= 'A' and c <= 'F')
    {
      value = 10 + c - 'A';
      return true;
    }

  if (c >= '0' and c <= '9')
    {
      value = c - '0';
      return true;
    }

  return false;
}


/// Convert an 8-bit value to 2 headecimal characters.
static
inline
void
byteToHexChars(uint8_t& byte, uint8_t& highDigit, uint8_t& lowDigit)
{
  static char hexDigits[] = "0123456789abcdef";
  highDigit = hexDigits[(byte >> 4)];
  lowDigit = hexDigits[byte & 0xf];
}
  

static
bool
getStringComponents(const std::string& str, char delim,
		    std::string& comp1, std::string& comp2)
{
  auto sepIx = str.find(delim);
  if (sepIx == std::string::npos)
    return false;

  comp1 = str.substr(0, sepIx);
  comp2 = str.substr(sepIx + 1);
  return true;
}


static
bool
getStringComponents(const std::string& str, char delim1, char delim2,
		    std::string& comp1, std::string& comp2, std::string& comp3)
{
  auto delim1Ix = str.find(delim1);
  if (delim1Ix == std::string::npos)
    return false;

  comp1 = str.substr(0, delim1Ix);

  auto delim2Ix = str.find(delim2, delim1Ix + 1);
  if (delim1Ix == std::string::npos)
    return false;

  comp2 = str.substr(delim1Ix + 1, delim2Ix - delim1Ix - 1);
  comp3 = str.substr(delim2Ix + 1);
  return true;
}


// Receive a packet from gdb. Request a retransmit from gdb if packet
// checksum is incorrect. Return succesfully received packet.
static
void
receivePacketFromGdb(int fd, std::string& packet)
{
  unsigned char ch = ' '; // Anything besides $ will do.

  char buffer[1024];
  uint8_t sum = 0;  // checksum

  ssize_t count = read(fd, buffer, sizeof(buffer));
  if (count < 0)
    return;

  // Packet starts with a '$'
  ssize_t ix = 0;
  while (ix < count and buffer[ix] != '$')
    ++ix;
      
  ++ix; // Skip '$'

  while (ix < count and buffer[ix] != '#')
    {
      char cc = buffer[ix++];
      sum += cc;
      packet.push_back(cc);
    }

  ++ix;  // Skip '#'

  if (ix + 1 < count)
    {
      ch = buffer[ix++];
      uint8_t pacSum = static_cast<uint8_t>(hexCharToInt(ch) << 4); // Packet checksum
      ch = buffer[ix++];
      pacSum = static_cast<uint8_t>(pacSum + hexCharToInt(ch));

      if (sum != pacSum)
        {
          std::cerr << "Bad checksum form gdb: "
                    << (boost::format("%02x v %02x") % unsigned(sum) %
                        unsigned(pacSum))
                    << '\n';
          putDebugChar('-', fd); // Signal failed reception.
        }
      else
        putDebugChar('+', fd);  // Signal successul reception.
    }
}


// Send given data string as a gdb remote packet. Resend until a
// positive ack is received.
//
// Format of packet:  $<data>#<checksum>
//
// TODO: quote special characters.
static void
sendPacketToGdb(const std::string& data, int fd)
{
  const char hexDigit[] = "0123456789abcdef";

  unsigned char checksum = 0;
  for (unsigned char c : data)
    checksum = static_cast<uint8_t>(checksum + c);

  std::string packet;
  packet.reserve(8 + data.size());

  packet.push_back('$');
  packet += data;
  packet.push_back('#');
  packet.push_back(hexDigit[checksum >> 4]);
  packet.push_back(hexDigit[checksum & 0xf]);

  while (true)
    {
      write(fd, packet.data(), packet.size());
      char c = getDebugChar(fd);
      if (c == '+')
	return;
    }
}


/// Return hexadecimal representation of given integer register value.
template <typename T>
std::string
littleEndianIntToHex(T val)
{
  std::ostringstream oss;

  for (size_t i = 0; i < sizeof(T); ++i)
    {
      unsigned byte = val & 0xff;
      val = val >> 8;
      oss << (boost::format("%02x") % byte);
    }

  return oss.str();
}


/// Convert given little-endian hexadecimal string to an integer value
/// of type T. Return true on success and false if the string does not
/// contain a hexadecimal string or if it is too large for type T.
template <typename T>
bool
littleEndianHexToInt(const std::string& str, T& value)
{
  size_t byteCount = 0;

  value = 0;

  for (size_t i = 0; i < str.size(); i += 2)
    {
      char c = str.at(i);
      unsigned byte =  0;
      if (not hexCharToInt(c, byte))
	return false;

      if (i + 1 < str.size())
	{
	  byte <<= 4;

	  char c2 = str.at(i+1);
	  unsigned x = 0;
	  if (not hexCharToInt(c2, x))
	    return false;
	  byte |= x;
	}
      value |= byte << (byteCount*8);
      byteCount++;
    }

  return byteCount <= sizeof(T);
}


/// Convert given hexadecimal string to an integer value of type
/// T. Return true on success and false if the string does not contain
/// a hexadecimal string or if it is too large for type T.
template <typename T>
bool
hexToInt(const std::string& str, T& value)
{
  value = 0;

  if (str.empty())
    return false;

  const char* data = str.c_str();
  char* end = nullptr;
  uint64_t v = strtoull(data, &end, 16);
  if (*end)
    return false;

  value = static_cast<T>(v);
  if (v != value)
    return false; // Overflow

  return true;
}


// GBD uses indices 0-31 for integer registers, 32 for pc, 33-64 for
// floating-point registers, 65 and higher for CSRs.
static const unsigned fpRegOffset = 33, pcOffset = 32, csrOffset = 65;


template <typename URV>
void
handlePeekRegisterForGdb(WdRiscv::Hart<URV>& hart, unsigned regNum,
			 std::ostream& stream)
{
  URV value = 0; bool ok = true, fp = false;
  if (regNum < pcOffset)
    ok = hart.peekIntReg(regNum, value);
  else if (regNum == pcOffset)
    value = hart.peekPc();
  else if (regNum >= fpRegOffset and regNum < csrOffset)
    {
      fp = true;
      if (hart.isRvf() or hart.isRvd())
	{
	  unsigned fpReg = regNum - fpRegOffset;
	  uint64_t val64 = 0;
	  hart.peekUnboxedFpReg(fpReg, val64);
          stream << littleEndianIntToHex(val64);
	}
      else
	stream << littleEndianIntToHex(0L); // Old versions of gdb cannot handle error "E03";
    }
  else
    {
      URV csr = regNum - csrOffset;
      // Undocumented: gdb uses 0x1000 to obtain privilege mode.
      if (csr == 0x1000)
        value = URV(hart.privilegeMode());
      else
        hart.peekCsr(WdRiscv::CsrNumber(csr), value);
    }

  if (ok)
    {
      if (not fp)
	stream << littleEndianIntToHex(value);
    }
  else
    stream << "E04";
}


template <typename URV>
void
getGdbTargetXml(WdRiscv::Hart<URV>& hart, std::string& xml)
{
  xml.clear();

  xml += R"zzz(<?xml version="1.0"?>
<!DOCTYPE target SYSTEM "gdb-target.dtd">
<target version="1.0">
  <feature name="org.gnu.gdb.riscv.cpu">
)zzz";

  // Integer registers.
  std::string width = std::to_string(8*sizeof(URV));

  for (unsigned ix = 0; ix < hart.intRegCount(); ++ix)
    {
      std::string name = hart.intRegName(ix);
      std::string num = std::to_string(ix);
      xml += "    <reg name=\"" + name + "\" bitsize=\"" + width + "\" regnum=\"";
      xml += num + "\" save-restore=\"yes\" type=\"int\" group=\"general\"/>\n";
    }

  std::string pcNum = std::to_string(pcOffset);
  xml += "    <reg name=\"pc\" bitsize=\"" + width + "\" regnum=\"";
  xml += pcNum + "\" save-restore=\"yes\" type=\"int\" group=\"general\"/>\n";

  xml += "  </feature>\n";

  // FP registers.
  if (hart.isRvf())
    {
      xml += "  <feature name=\"org.gnu.gdb.riscv.fpu\">\n";

      for (unsigned ix = 0; ix < hart.fpRegCount(); ++ix)
        {
          std::string name = hart.fpRegName(ix);
          std::string num = std::to_string(ix + fpRegOffset);
          xml += "    <reg name=\"" + name + "\" bitsize=\"64\" regnum=\"";
          xml += num + "\" save-restore=\"yes\" type=\"float\" group=\"fp\"/>\n";
        }

      xml += "  </feature>\n";
    }

  // CSRs
  std::vector<WdRiscv::CsrNumber> csrNumbers;
  hart.getImplementedCsrs(csrNumbers);

  xml += "  <feature name=\"org.gnu.gdb.riscv.csr\">\n";
  for (auto csrn : csrNumbers)
    {
      URV val = 0;
      std::string name;
      if (not hart.peekCsr(csrn, val, name))
        continue;

      std::string num = std::to_string(unsigned(csrn) + csrOffset);

      xml += "    <reg name=\"" + name + "\" bitsize=\"" + width + "\" regnum=\"";
      xml += num + "\" save-restore=\"no\" type=\"int\" group=\"csr\"/>\n";
    }
  xml += "  </feature>\n";

  // Virtual CSRs (privilege mode pseudo-csr)
  // CSR 0x1000 is used by gdb to obtain privilege mode.
  std::string num = std::to_string(0x1000 + csrOffset);
  xml += "  <feature name=\"org.gnu.gdb.riscv.virtual\">\n";
  xml += "    <reg name=\"priv\" bitsize=\"" + width + "\" regnum=\"";
  xml += num + "\" save-restore=\"no\" type=\"int\" group=\"general\"/>\n";
  xml += "  </feature>\n";

  xml += "</target>";
}


// XML describing this RISCV processor to gdb.
static std::string targetXml;
static std::mutex xmlMutex;

template <typename URV>
void
processXferQuery(const std::string& packet, WdRiscv::Hart<URV>& hart,
                 std::ostream& reply)
{
  if (not boost::starts_with(packet, "qXfer:features:read:target.xml:"))
    {
      reply << "";
      return;
    }

  // Request is of the form: qXfer:features:read:target.xml:<offset>,<length>
  // Recover offset and length substring.
  auto ix = packet.find_last_of(':');
  auto ol = packet.substr(ix + 1);  // offset and length
  unsigned offset = 0, length = 0;
  if (sscanf(ol.c_str(), "%x,%x", &offset, &length) != 2)
    {
      std::cerr << "Malformed qXfer request: " << packet << '\n';
      reply << "";
      return;
    }

  {
    std::lock_guard<std::mutex> lock(xmlMutex);

    // Fill xml string on first call to this function.
    if (targetXml.empty())
      getGdbTargetXml(hart, targetXml);
  }

  auto part = targetXml.substr(offset, length);
  if (offset + length < targetXml.size())
    reply << 'm' << part;
  else
    reply << 'l' << part;
}


// Called after a stop (caused by an exception). Notify GDB of the
// stop.  Return the signal number corresponding to the exception.
template <typename URV>
unsigned
notifyGdbAfterStop(WdRiscv::Hart<URV>& hart, int fd)
{
  // Construct a reply of the form T xx n1:r1;n2:r2;... where xx is
  // the trap cause and ni is a resource (e.g. register number) and ri
  // is the resource data (e.g. content of register).
  std::ostringstream reply;

  unsigned signalNum = SIGTRAP;
  URV cause = 0;
  if (hart.peekCsr(WdRiscv::CsrNumber::MCAUSE, cause))
    {
      if (cause == URV(WdRiscv::ExceptionCause::BREAKP))
	signalNum = SIGTRAP;
      // FIX:  implement other causes.
    }

  reply << "T" << (boost::format("%02x") % signalNum);

  URV spVal = 0;
  unsigned spNum = WdRiscv::RegSp;
  hart.peekIntReg(spNum, spVal);
  reply << (boost::format("%02x") % spNum) << ':'
	<< littleEndianIntToHex(spVal) << ';';
  sendPacketToGdb(reply.str(), fd);

  return signalNum;
}


template <typename URV>
void
handleExceptionForGdb(WdRiscv::Hart<URV>& hart, int fd)
{
  // The trap handler is expected to set the PC to point to the instruction
  // after the one with the exception if necessary/possible.
  unsigned signalNum = notifyGdbAfterStop(hart, fd);

  bool gotQuit = false;

  std::string buffer;
  buffer.reserve(128);
  std::ostringstream reply(buffer);

  std::string packet;
  packet.reserve(128);

  while (1)
    {
      reply.str("");
      reply.clear();

      packet.clear();
      receivePacketFromGdb(fd, packet);
      if (packet.empty())
	continue;

      switch (packet.at(0))
	{
	case '?':
	  // Return signal number
	  reply << "S" << (boost::format("%02x") % signalNum);
	  break;

	case 'g':  // return the value of the CPU registers
	  {
	    for (unsigned i = 0; i < hart.intRegCount(); ++i)
	      {
		URV val = 0;
		hart.peekIntReg(i, val);
		reply << littleEndianIntToHex(val);
	      }
            reply << littleEndianIntToHex(hart.peekPc());
	  }
	  break;

	case 'G':  // set the value of the CPU registers - return OK
	  {
	    std::string status = "OK";
	    size_t len = packet.length();
	    if (len + 1 < hart.intRegCount() * sizeof(URV) * 2)
	      status = "E01";
	    else
	      {
		const char* ptr = packet.c_str() + 1;
		for (unsigned i = 0; i < hart.intRegCount(); ++i)
		  {
		    std::string buffer;
		    for (unsigned i = 0; i < 2*sizeof(URV); ++i)
		      buffer += *ptr++;
		    URV val = 0;
		    if (littleEndianHexToInt(buffer, val))
		      hart.pokeIntReg(i, val);
		    else
		      {
			status  = "E01";
			break;
		      }
		  }
	      }

	    reply << status;
	  }
	  break;

	case 'H':   // Hc<thread> or Hg<thread>
	  {
	    if (packet.length() < 2)
	      reply << "E01";
	    else
	      {
		unsigned threadId = 0;
		if (packet[1] != 'c' or packet[1] != 'g')
		  reply << "E01";
		else if (not hexToInt(packet.substr(2), threadId))
		  reply << "E01";
		else if (threadId != 0)
		  reply << "E01";  // Multi-thread not supported yet.
		else
		  reply << "OK";
	      }
	  }
	  break;

	case 'm': // mAA..AA,LLLL  Read LLLL bytes at address AA..AA
	  {
	    std::string addrStr, lenStr;
	    if (not getStringComponents(packet.substr(1), ',', addrStr, lenStr))
	      reply << "E01";
	    else
	      {
		URV addr = 0, len = 0;
		if (not hexToInt(addrStr, addr) or not hexToInt(lenStr, len))
		  reply << "E02";
		else
		  {
		    for (URV ix = 0; ix < len; ++ix)
		      {
			uint8_t byte = 0, high = 0, low = 0;
                        bool usePma = false; // Ignore physical memory attributes.
			hart.peekMemory(addr++, byte, usePma);
                        byteToHexChars(byte, high, low);
                        reply.write(reinterpret_cast<char*> (&high), 1);
                        reply.write((char*) &low, 1);
		      }
		  }
	      }
	  }
	  break;

	case 'M': // MAA..AA,LLLL: Write LLLL bytes at address AA.AA return OK
	  {
	    std::string addrStr, lenStr, data;
	    if (not getStringComponents(packet.substr(1), ',', ':', addrStr,
					lenStr, data))
	      reply << "E01";
	    else
	      {
		URV addr = 0, len = 0;
		if (not hexToInt(addrStr, addr) or not hexToInt(lenStr, len))
		  reply << "E02";
		else
		  {
		    if (data.size() < len*2)
		      reply << "E03";
		    else
		      {
			for (URV ix = 0; ix < len; ++ix)
			  {
			    int bb = hexCharToInt(data.at(2*ix));
			    bb = (bb << 4) | hexCharToInt(data.at(2*ix+1));
                            uint8_t val = bb;
                            bool usePma = false; // Ignore physical memory attributes.
			    hart.pokeMemory(addr, val, usePma);
                            addr++;
			  }
			reply << "OK";
		      }
		  }
	      }
	  }
	  break;

	case 'c':  // cAA..AA    Continue at address AA..AA(optional)
	  {
            if (not hart.isStarted())
              {
                std::cerr << "Hart not yet started (need to write mhartstart)\n";
                reply << "E01";
                break;
              }

	    if (packet.size() == 1)
	      return;

	    URV newPc = 0;
	    if (hexToInt(packet.substr(1), newPc))
	      {
		hart.pokePc(newPc);
		return;
	      }

	    reply << "E01";
	  }
	  break;

	case 'p':  // pn    Read value of register n
	  {
	    std::string regNumStr = packet.substr(1);
	    unsigned regNum = 0;

	    if (not hexToInt(regNumStr, regNum))
	      reply << "E01";
	    else
	      handlePeekRegisterForGdb(hart, regNum, reply);
	  }
	  break;

	case 'P':   // Pn=v   Set register n to v
	  {
	    auto eqIx = packet.find('=');
	    if (eqIx == std::string::npos or eqIx == 1)
	      reply << "E01";
	    else
	      {
		std::string regNumStr = packet.substr(1, eqIx - 1);
		std::string valueStr = packet.substr(eqIx + 1);
		unsigned regNum = 0;
		URV value = 0;
		if (not hexToInt(regNumStr, regNum))
		  reply << "E02";
		else if (not littleEndianHexToInt(valueStr, value))
		  reply << "E03";
		else
		  {
                    bool ok = true;
		    if (regNum < pcOffset)
		      ok = hart.pokeIntReg(regNum, value);
		    else if (regNum == pcOffset)
		      hart.pokePc(value);
		    else if (regNum >= fpRegOffset and regNum < csrOffset)
		      ok = hart.pokeFpReg(regNum - fpRegOffset, value);
                    else if (regNum >= csrOffset)
                      ok = hart.pokeCsr(WdRiscv::CsrNumber(regNum - csrOffset), value);
		    reply << (ok? "OK" : "E04");
		  }
	      }
	  }
	  break;

	case 's':
          if (not hart.isStarted())
            {
              std::cerr << "Hart not yet started (need to write mhartstart)\n";
              reply << "E01";
              break;
            }

	  hart.singleStep(nullptr);
	  notifyGdbAfterStop(hart, fd);
	  continue;
	  break;

	case 'k':  // kill
	  reply << "OK";
	  gotQuit = true;
	  break;

	case 'q':
	  if (packet == "qC")
	    reply << "QC 0";
	  else if (packet == "qAttached")
	    reply << "0";
	  else if (packet == "qOffsets")
	    reply << "Text=0;Data=0;Bss=0";
	  else if (packet == "qSymbol::")
	    reply << "OK";
	  else if (packet == "qfThreadInfo")
	    reply << "m0";
	  else if (packet == "qsThreadInfo")
	    reply << "l";
	  else if (packet == "qTStatus")
	    reply << "T0;tnotrun:0";
          else if (boost::starts_with(packet, "qSupported"))
            reply << "qXfer:features:read+";
          else if (boost::starts_with(packet, "qXfer"))
            processXferQuery(packet, hart, reply);
	  else
	    {
	      std::cerr << "Unhandled gdb request: " << packet << '\n';
	      reply << ""; // Unsupported: Empty response.
	    }
	  break;

	case 'v':
	  if (packet == "vMustReplyEmpty")
	    reply << "";
	  else if (packet.find("vKill;") == 0)
	    {
	      reply << "OK";
	      gotQuit = true;
	    }
	  else if (packet == "vCont?")
	    {
	      reply << "";
	    }
	  else
	    {
	      std::cerr << "Unhandled gdb request: " << packet << '\n';
	      reply << ""; // Unsupported: Empty response.
	    }
	  break;

	default:
	  std::cerr << "Unhandled gdb request: " << packet << '\n';
	  reply << "";   // Unsupported comand: Empty response.
	}

      // Reply to the request
      sendPacketToGdb(reply.str(), fd);

      if (gotQuit)
	exit(0);
    }
}


template void handleExceptionForGdb<uint32_t>(WdRiscv::Hart<uint32_t>&, int);
template void handleExceptionForGdb<uint64_t>(WdRiscv::Hart<uint64_t>&, int);
