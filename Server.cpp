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
#include <sstream>
#include <fstream>
#include <map>
#include <algorithm>
#include <boost/format.hpp>
#include <cstring>
#ifdef __MINGW64__
#include <winsock2.h>
#else
#include <arpa/inet.h>
#include <sys/socket.h>
#endif

#define __STDC_FORMAT_MACROS
#include <cinttypes>

#include "WhisperMessage.h"
#include "Hart.hpp"
#include "Server.hpp"
#include "Interactive.hpp"


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


/// Unpack socket message (received in server mode) into the given
/// WhisperMessage object.
static
void
deserializeMessage(const char buffer[], size_t bufferLen,
		   WhisperMessage& msg)

{
  assert (bufferLen >= sizeof(msg));

  const char* p = buffer;
  uint32_t x = ntohl(* reinterpret_cast<const uint32_t*> (p));
  msg.hart = x;
  p += sizeof(x);

  x = ntohl(* reinterpret_cast<const uint32_t*> (p));
  msg.type = x;
  p += sizeof(x);

  x = ntohl(* reinterpret_cast<const uint32_t*> (p));
  msg.resource = x;
  p += sizeof(x);

  x = ntohl(* reinterpret_cast<const uint32_t*> (p));
  msg.flags = x;
  p += sizeof(x);

  uint32_t part = ntohl(* reinterpret_cast<const uint32_t*> (p));
  msg.rank = uint64_t(part) << 32;
  p += sizeof(part);

  part = ntohl(* reinterpret_cast<const uint32_t*> (p));
  msg.rank |= part;
  p += sizeof(part);

  part = ntohl(* reinterpret_cast<const uint32_t*> (p));
  msg.address = uint64_t(part) << 32;
  p += sizeof(part);

  part = ntohl(* reinterpret_cast<const uint32_t*> (p));
  msg.address |= part;
  p += sizeof(part);

  part = ntohl(* reinterpret_cast<const uint32_t*> (p));
  msg.value = uint64_t(part) << 32;
  p += sizeof(part);

  part = ntohl(* reinterpret_cast<const uint32_t*> (p));
  msg.value |= part;
  p += sizeof(part);

  memcpy(msg.buffer, p, sizeof(msg.buffer));
  p += sizeof(msg.buffer);

  memcpy(msg.tag, p, sizeof(msg.tag));
  p += sizeof(msg.tag);

  assert(size_t(p - buffer) <= bufferLen);
}


/// Serialize the given WhisperMessage into the given buffer in
/// preparation for socket send. Return the number of bytes written
/// into buffer.
static
size_t
serializeMessage(const WhisperMessage& msg, char buffer[],
		 size_t bufferLen)
{
  assert (bufferLen >= sizeof(msg));

  char* p = buffer;
  uint32_t x = htonl(msg.hart);
  memcpy(p, &x, sizeof(x));
  p += sizeof(x);

  x = htonl(msg.type);
  memcpy(p, &x, sizeof(x));
  p += sizeof(x);

  x = htonl(msg.resource);
  memcpy(p, &x, sizeof(x));
  p += sizeof(x);

  x = htonl(msg.flags);
  memcpy(p, &x, sizeof(x));
  p += sizeof(x);

  uint32_t part = static_cast<uint32_t>(msg.rank >> 32);
  x = htonl(part);
  memcpy(p, &x, sizeof(x));
  p += sizeof(x);

  part = (msg.rank) & 0xffffffff;
  x = htonl(part);
  memcpy(p, &x, sizeof(x));
  p += sizeof(x);

  part = static_cast<uint32_t>(msg.address >> 32);
  x = htonl(part);
  memcpy(p, &x, sizeof(x));
  p += sizeof(x);

  part = (msg.address) & 0xffffffff;
  x = htonl(part);
  memcpy(p, &x, sizeof(x));
  p += sizeof(x);

  part = static_cast<uint32_t>(msg.value >> 32);
  x = htonl(part);
  memcpy(p, &x, sizeof(x));
  p += sizeof(x);

  part = msg.value & 0xffffffff;
  x = htonl(part);
  memcpy(p, &x, sizeof(x));
  p += sizeof(x);

  memcpy(p, msg.buffer, sizeof(msg.buffer));
  p += sizeof(msg.buffer);

  memcpy(p, msg.tag, sizeof(msg.tag));
  p += sizeof(msg.tag);

  size_t len = p - buffer;
  assert(len <= bufferLen);
  assert(len <= sizeof(msg));
  for (size_t i = len; i < sizeof(msg); ++i)
    buffer[i] = 0;

  return sizeof(msg);
}


static bool
receiveMessage(int soc, WhisperMessage& msg)
{
  char buffer[sizeof(msg)];
  char* p = buffer;

  size_t remain = sizeof(msg);

  while (remain > 0)
    {
      ssize_t l = recv(soc, p, remain, 0);
      if (l < 0)
	{
	  if (errno == EINTR)
	    continue;
	  std::cerr << "Failed to receive socket message\n";
	  return false;
	}
      if (l == 0)
	{
	  msg.type = Quit;
	  return true;
	}
      remain -= l;
      p += l;
    }

  deserializeMessage(buffer, sizeof(buffer), msg);

  return true;
}


static bool
sendMessage(int soc, WhisperMessage& msg)
{
  char buffer[sizeof(msg)];

  serializeMessage(msg, buffer, sizeof(buffer));

  // Send command.
  ssize_t remain = sizeof(msg);
  char* p = buffer;
  while (remain > 0)
    {
      ssize_t l = send(soc, p, remain , 0);
      if (l < 0)
	{
	  if (errno == EINTR)
	    continue;
	  std::cerr << "Failed to send socket command\n";
	  return false;
	}
      remain -= l;
      p += l;
    }

  return true;
}


template <typename URV>
Server<URV>::Server(System<URV>& system)
  : system_(system)
{
}
  

template <typename URV>
bool
Server<URV>::pokeCommand(const WhisperMessage& req, WhisperMessage& reply)
{
  reply = req;

  if (not checkHartId(req, reply))
    return false;

  uint32_t hartId = req.hart;
  auto hartPtr = system_.findHartByHartId(hartId);
  if (not hartPtr)
    return false;
  auto& hart = *hartPtr;

  switch (req.resource)
    {
    case 'r':
      {
	unsigned reg = static_cast<unsigned>(req.address);
	URV val = static_cast<URV>(req.value);
	if (reg == req.address)
	  if (hart.pokeIntReg(reg, val))
	    return true;
      }
      break;

    case 'f':
      {
	unsigned reg = static_cast<unsigned>(req.address);
	uint64_t val = req.value;
	if (reg == req.address)
	  if (hart.pokeFpReg(reg, val))
	    return true;
      }
      break;

    case 'c':
      {
	URV val = static_cast<URV>(req.value);
	if (hart.pokeCsr(CsrNumber(req.address), val))
	  return true;
      }
      break;

    case 'm':
      {
        bool usePma = false; // Ignore phsical memory attributes.

        if (sizeof(URV) == 4)
          {
            // Poke a word in 32-bit harts.
            if (hart.pokeMemory(req.address, uint32_t(req.value), usePma))
              return true;
          }
        else if (hart.pokeMemory(req.address, req.value, usePma))
          return true;
      }
      break;

    case 'p':
      {
	URV val = static_cast<URV>(req.value);
	hart.pokePc(val);
	return true;
      }
      break;
    }

  reply.type = Invalid;
  return false;
}


template <typename URV>
bool
Server<URV>::peekCommand(const WhisperMessage& req, WhisperMessage& reply)
{
  reply = req;

  if (not checkHartId(req, reply))
    return false;

  uint32_t hartId = req.hart;
  auto hartPtr = system_.findHartByHartId(hartId);
  if (not hartPtr)
    return false;
  auto& hart = *hartPtr;

  URV value;

  switch (req.resource)
    {
    case 'r':
      {
	unsigned reg = static_cast<unsigned>(req.address);
	if (reg == req.address)
	  if (hart.peekIntReg(reg, value))
	    {
	      reply.value = value;
	      return true;
	    }
      }
      break;
    case 'f':
      {
	unsigned reg = static_cast<unsigned>(req.address);
	uint64_t fpVal = 0;
	if (reg == req.address)
	  if (hart.peekFpReg(reg, fpVal))
	    {
	      reply.value = fpVal;
	      return true;
	    }
      }
      break;
    case 'c':
      if (hart.peekCsr(CsrNumber(req.address), value))
	{
	  reply.value = value;
	  return true;
	}
      break;
    case 'm':
      if (hart.peekMemory(req.address, value, false /*usePma*/))
	{
	  reply.value = value;
	  return true;
	}
      break;
    case 'p':
      reply.value = hart.peekPc();
      return true;
    }

  reply.type = Invalid;
  return true;
}


template <typename URV>
void
Server<URV>::disassembleAnnotateInst(Hart<URV>& hart,
                                     uint32_t inst, bool interrupted,
				     bool hasPreTrigger, bool hasPostTrigger,
				     std::string& text)
{
  hart.disassembleInst(inst, text);
  uint32_t op0 = 0, op1 = 0, op2 = 0, op3 = 0;
  const InstEntry& entry = hart.decode(inst, op0, op1, op2, op3);
  if (entry.isBranch())
    {
      if (hart.lastPc() + instructionSize(inst) != hart.peekPc())
       text += " (T)";
      else
       text += " (NT)";
    }

  if (not interrupted)
    if (entry.isLoad() or entry.isStore() or entry.isAtomic())
      {
        URV addr = hart.lastLdStAddress();
        std::ostringstream oss;
        oss << " [0x" << std::hex << addr << "]" << std::dec;
        text += oss.str();
        bool cacheable = hart.isAddrCacheable(addr);
        bool io = not hart.isAddrCacheable(addr);
        if (cacheable or io)
          {
            text += " (";
            if (cacheable)
              text += "C";
            if (io)
              text += "S";
            text += ")";
          }
      }

  if (interrupted)
    text += " (interrupted)";
  else if (hasPreTrigger)
    text += " (pre-trigger)";
  else if (hasPostTrigger)
    text += " (post-trigger)";
}


/// Collect the double-words overlapping the areas of memory modified
/// by the most recent emulated system call. It is assumed that
/// memory protection boundaries are double-word aligned.
template <typename URV>
static void
collectSyscallMemChanges(Hart<URV>& hart,
                         const std::vector<std::pair<uint64_t, uint64_t>>& scVec,
                         std::vector<WhisperMessage>& changes,
                         uint64_t slamAddr)
{
  auto hexForm = getHexForm<URV>(); // Format string for printing a hex val

  for (auto al : scVec)
    {
      uint64_t addr = al.first;
      uint64_t len = al.second;

      uint64_t offset = addr & 7;  // Offset from double word address.
      if (offset)
        {
          addr -= offset;   // Make double word aligned.
          len += offset;    // Keep last byte the same.
        }

      for (uint64_t ix = 0; ix < len; ix += 8, addr += 8)
        {
          uint64_t val = 0;
          if (not hart.peekMemory(addr, val, true))
            {
              std::cerr << "Peek-memory fail at 0x%x"
                        << (boost::format(hexForm) % addr)
                        << " in collectSyscallMemChanges\n";
              break;
            }

          if (not slamAddr)
            continue;

          bool ok = hart.pokeMemory(slamAddr, addr, true);
          if (ok)
            {
              changes.push_back(WhisperMessage{0, Change, 'm', slamAddr, addr});

              slamAddr += 8;
              ok = hart.pokeMemory(slamAddr, val, true);
              if (ok)
                {
                  changes.push_back(WhisperMessage{0, Change, 'm', slamAddr, val});
                  slamAddr += 8;
                }
            }

          if (not ok)
            {
              std::cerr << "Poke-memory fail at 0x%x"
                        << (boost::format(hexForm) % slamAddr)
                        << " in collectSyscallMemChanges\n";
              break;
            }
        }
    }

  // Put a zero to mark end of syscall changes in slam area.
  if (slamAddr)
    {
      if (hart.pokeMemory(slamAddr, uint64_t(0), true))
        changes.push_back(WhisperMessage{0, Change, 'm', slamAddr, 0});
      if (hart.pokeMemory(slamAddr + 8, uint64_t(0), true))
        changes.push_back(WhisperMessage{0, Change, 'm', slamAddr+8, 0});
    }
}


template <typename URV>
void
Server<URV>::processStepCahnges(Hart<URV>& hart,
				uint32_t inst,
				std::vector<WhisperMessage>& pendingChanges,
				bool interrupted, bool hasPre, bool hasPost,
				WhisperMessage& reply)
{
  // Get executed instruction.
  URV pc = hart.lastPc();

  // Add pc and instruction to reply.
  reply.type = ChangeCount;
  reply.address = pc;
  reply.resource = inst;

  // Add disassembly of instruction to reply.
  std::string text;
  disassembleAnnotateInst(hart, inst, interrupted, hasPre, hasPost, text);

  strncpy(reply.buffer, text.c_str(), sizeof(reply.buffer) - 1);
  reply.buffer[sizeof(reply.buffer) -1] = 0;

  // Order of changes: rfcvm (int reg, fp reg, csr, vec reg, memory, csr)

  // Collect integer register change caused by execution of instruction.
  pendingChanges.clear();
  int regIx = hart.lastIntReg();
  if (regIx > 0)
    {
      URV value = 0;
      if (hart.peekIntReg(regIx, value))
	{
	  WhisperMessage msg;
	  msg.type = Change;
	  msg.resource = 'r';
	  msg.address = regIx;
	  msg.value = value;
	  pendingChanges.push_back(msg);
	}
    }

  // Collect floating point register change.
  int fpRegIx = hart.lastFpReg();
  if (fpRegIx >= 0)
    {
      uint64_t val = 0;
      if (hart.peekFpReg(fpRegIx, val))
	{
	  WhisperMessage msg;
	  msg.type = Change;
	  msg.resource = 'f';
	  msg.address = fpRegIx;
	  msg.value = val;
	  pendingChanges.push_back(msg);
	}
    }

  // Collect vector register change (format not yet defined).

  // Collect CSR and trigger changes.
  std::vector<CsrNumber> csrs;
  std::vector<unsigned> triggers;
  hart.lastCsr(csrs, triggers);

  // Map to keep CSRs in order and to drop duplicate entries.
  std::map<URV,URV> csrMap;

  // Components of the triggers that changed (if any).
  std::vector<bool> tdataChanged(3);

  // Collect changed CSRs and their values. Collect components of
  // changed trigger.
  for (CsrNumber csr : csrs)
    {
      URV value;
      if (hart.peekCsr(csr, value))
	{
	  if (csr >= CsrNumber::TDATA1 and csr <= CsrNumber::TDATA3)
	    {
	      size_t ix = size_t(csr) - size_t(CsrNumber::TDATA1);
	      tdataChanged.at(ix) = true;
	    }
	  else
	    csrMap[URV(csr)] = value;
	}
    }

  // Collect changes associated with trigger register.
  for (unsigned trigger : triggers)
    {
      uint64_t data1(0), data2(0), data3(0);
      if (not hart.peekTrigger(trigger, data1, data2, data3))
	continue;

      // Components of trigger that changed.
      bool t1 = false, t2 = false, t3 = false;
      hart.getTriggerChange(trigger, t1, t2, t3);

      if (t1)
	{
	  URV addr = (trigger << 16) | unsigned(CsrNumber::TDATA1);
	  csrMap[addr] = data1;
	}
      if (t2)
	{
	  URV addr = (trigger << 16) | unsigned(CsrNumber::TDATA2);
	  csrMap[addr] = data2;
	}
      if (t3)
	{
	  URV addr = (trigger << 16) | unsigned(CsrNumber::TDATA3);
	  csrMap[addr] = data3;
	}
    }

  for (const auto& [key, val] : csrMap)
    {
      WhisperMessage msg(0, Change, 'c', key, val);
      pendingChanges.push_back(msg);
    }

  // Collect memory change.
  uint64_t memAddr = 0, memVal = 0;
  unsigned size = hart.lastMemory(memAddr, memVal);
  if (size)
    {
      WhisperMessage msg(0, Change, 'm', memAddr, memVal);
      msg.flags = size;
      pendingChanges.push_back(msg);
    }

  // Collect emulated system call changes.
  uint64_t slamAddr = hart.syscallSlam();
  if (slamAddr)
    {
      std::vector<std::pair<uint64_t, uint64_t>> scVec;
      hart.lastSyscallChanges(scVec);
      if (scVec.size())
        collectSyscallMemChanges(hart, scVec, pendingChanges, slamAddr);
    }

  // Add count of changes to reply.
  reply.value = pendingChanges.size();

  // The changes will be retrieved one at a time from the back of the
  // pendigChanges vector: Put the vector in reverse order. Changes
  // are retrieved using a Change request (see interactUsingSocket).
  std::reverse(pendingChanges.begin(), pendingChanges.end());
}


template <typename URV>
bool
Server<URV>::checkHartId(const WhisperMessage& req, WhisperMessage& reply)
{
  uint32_t hartId = req.hart;
  auto hartPtr = system_.findHartByHartId(hartId);
  if (not hartPtr)
    {
      std::cerr << "Error: Hart ID (" << std::dec << hartId
                << ") out of bounds\n";
      reply.type = Invalid;
      return false;
    }
  return true;
}


template <typename URV>
bool
Server<URV>::checkHart(const WhisperMessage& req, const std::string& command,
                       WhisperMessage& reply)
{
  if (not checkHartId(req, reply))
    return false;

  uint32_t hartId = req.hart;
  auto hartPtr = system_.findHartByHartId(hartId);
  if (not hartPtr)
    return false;

  auto& hart = *hartPtr;
  if (not hart.isStarted())
    {
      std::cerr << "Error: Command " << command
                << " received for a non-started hart\n";
      reply.type = Invalid;
      return false;
    }

  return true;
}


// Server mode step command.
template <typename URV>
bool
Server<URV>::stepCommand(const WhisperMessage& req, 
			 std::vector<WhisperMessage>& pendingChanges,
			 WhisperMessage& reply,
			 FILE* traceFile)
{
  reply = req;

  // Hart id must be valid. Hart must be started.
  if (not checkHart(req, "step", reply))
    return false;

  uint32_t hartId = req.hart;
  auto hartPtr = system_.findHartByHartId(hartId);
  if (not hartPtr)
    return false;
  auto& hart = *hartPtr;

  // Step is not allowed in debug mode unless we are in debug_step as
  // well.
  if (hart.inDebugMode() and not hart.inDebugStepMode())
    {
      std::cerr << "Error: Single step while in debug-halt mode\n";
      reply.type = Invalid;
      return false;
    }

  // Get instruction before execution (in case code is self-modifying).
  uint32_t inst = 0;
  hart.readInst(hart.peekPc(), inst);

  // Get privilege mode.
  int privMode = int(hart.privilegeMode());

  // Execute instruction. Determine if an interrupt was taken or if a
  // trigger got tripped.
  uint64_t interruptCount = hart.getInterruptCount();

  hart.singleStep(traceFile);

  bool interrupted = hart.getInterruptCount() != interruptCount;

  unsigned preCount = 0, postCount = 0;
  hart.countTrippedTriggers(preCount, postCount);

  bool hasPre = preCount > 0;
  bool hasPost = postCount > 0;

  processStepCahnges(hart, inst, pendingChanges, interrupted, hasPre,
		     hasPost, reply);

  // Send privilege mode in reply.flags
  reply.flags = privMode;

#if 0
  // Send floating point flags in bits 16 to 19 of reply.flags.
  unsigned fpFlags = 0;
  reply.flags |= (fpFlags << 16);
#endif

  hart.clearTraceData();
  return true;
}


// Server mode exception command.
template <typename URV>
bool
Server<URV>::exceptionCommand(const WhisperMessage& req, 
			      WhisperMessage& reply,
			      std::string& text)
{
  reply = req;

  if (not checkHart(req, "exception", reply))
    return false;

  uint32_t hartId = req.hart;
  auto hartPtr = system_.findHartByHartId(hartId);
  if (not hartPtr)
    return false;
  auto& hart = *hartPtr;

  std::ostringstream oss;

  bool ok = true;
  URV addr = static_cast<URV>(req.address);
  if (addr != req.address)
    std::cerr << "Error: Address too large (" << std::hex << req.address
	      << ") in exception command.\n" << std::dec;

  WhisperExceptionType expType = WhisperExceptionType(req.value);
  switch (expType)
    {
    case InstAccessFault:
      hart.postInstAccessFault(addr);
      oss << "exception inst " << addr;
      break;

    case DataAccessFault:
      hart.postDataAccessFault(addr, SecondaryCause::LOAD_ACC_DOUBLE_ECC);
      oss << "exception data " << addr;
      break;

    case ImpreciseStoreFault:
      {
        unsigned count = 0;
        ok = hart.applyStoreException(addr, count);
        reply.value = count;
      }
      oss << "exception store 0x" << std::hex << addr << std::dec;
      break;

    case ImpreciseLoadFault:
      {
	unsigned tag = reply.flags;  // Tempoary.
        unsigned count = 0;
        ok = hart.applyLoadException(addr, tag, count);
	reply.value = count;
        oss << "exception load 0x" << std::hex << addr << std::dec << ' '
            << tag;
      }
      break;

    case NonMaskableInterrupt:
      if (hart.isNmiEnabled())
        hart.setPendingNmi(NmiCause(addr));
      oss << "exception nmi 0x" << std::hex << addr << std::dec;
      break;

    case PreciseLoadFault:
      oss << "exception precise_load 0x" << std::hex << addr << std::dec;
      hart.postDataAccessFault(addr, SecondaryCause::LOAD_ACC_PRECISE);
      ok = false;
      break;

    case PreciseStoreFault:
      hart.postDataAccessFault(addr, SecondaryCause::STORE_ACC_PRECISE);
      oss << "exception precise_store 0x" << std::hex << addr << std::dec;
      ok = false;
      break;

    default:
      oss << "exception ? 0x" << std::hex << addr << std::dec;
      ok = false;
      break;
    }

  if (not ok)
    reply.type = Invalid;

  text = oss.str();
  return ok;
}


/// Dump all registers contents in tracefile.
template <typename URV>
static void
serverPrintFinalRegisterState(std::shared_ptr<Hart<URV>> hartPtr)
{
  std::ofstream out("issfinal.log");
  if (not out)
    return;
  Interactive<URV>::peekAllIntRegs(*hartPtr, out);
  out << "\n";
  Interactive<URV>::peekAllFpRegs(*hartPtr, out);
  out << "\n";
  Interactive<URV>::peekAllTriggers(*hartPtr, out);
  out << "\n";
  Interactive<URV>::peekAllCsrs(*hartPtr, out);
}


// Server mode loop: Receive command and send reply till a quit
// command is received. Return true on successful termination (quit
// received). Return false otherwise.
template <typename URV>
bool
Server<URV>::interact(int soc, FILE* traceFile, FILE* commandLog)
{
  std::vector<WhisperMessage> pendingChanges;

  auto hexForm = getHexForm<URV>(); // Format string for printing a hex val

  // Initial resets do not reset memory mapped registers.
  bool resetMemoryMappedReg = false;

  while (true)
    {
      WhisperMessage msg;
      if (not receiveMessage(soc, msg))
	return false;
      WhisperMessage reply = msg;

      if (checkHartId(msg, reply))
	{
          std::string timeStamp = std::to_string(msg.rank);

          uint32_t hartId = msg.hart;
          auto hartPtr = system_.findHartByHartId(hartId);
          assert(hartPtr);
          auto& hart = *hartPtr;

	  if (msg.type == Step or msg.type == Until)
	    resetMemoryMappedReg = true;

	  switch (msg.type)
	    {
	    case Quit:
	      if (commandLog)
		fprintf(commandLog, "hart=%d quit\n", hartId);
              serverPrintFinalRegisterState(hartPtr);
	      return true;

	    case Poke:
	      pokeCommand(msg, reply);
	      if (commandLog)
		{
		  if (msg.resource == 'p')
                    fprintf(commandLog, "hart=%d poke pc %s # ts=%s tag=%s\n", hartId,
			    (boost::format(hexForm) % msg.value).str().c_str(),
			    timeStamp.c_str(), msg.tag);
		  else
		    fprintf(commandLog, "hart=%d poke %c %s %s # ts=%s tag=%s\n", hartId,
			    msg.resource,
			    (boost::format(hexForm) % msg.address).str().c_str(),
			    (boost::format(hexForm) % msg.value).str().c_str(),
			    timeStamp.c_str(), msg.tag);
		}
	      break;

	    case Peek:
	      peekCommand(msg, reply);
	      if (commandLog)
                {
                  if (msg.resource == 'p')
                    fprintf(commandLog, "hart=%d peek pc # ts=%s tag=%s\n",
                            hartId, timeStamp.c_str(), msg.tag);
                  else
                    fprintf(commandLog, "hart=%d peek %c %s # ts=%s tag=%s\n",
                            hartId,
                            msg.resource,
                            (boost::format(hexForm) % msg.address).str().c_str(),
                            timeStamp.c_str(), msg.tag);
                }
	      break;

	    case Step:
	      stepCommand(msg, pendingChanges, reply, traceFile);
	      if (commandLog)
		fprintf(commandLog, "hart=%d step #%" PRId64 " # ts=%s\n",
                        hartId, hart.getInstructionCount(), timeStamp.c_str());
	      break;

	    case ChangeCount:
	      reply.type = ChangeCount;
	      reply.value = pendingChanges.size();
	      reply.address = hart.lastPc();
	      {
		uint32_t inst = 0;
		hart.readInst(hart.lastPc(), inst);
		reply.resource = inst;
		std::string text;
		hart.disassembleInst(inst, text);
		uint32_t op0 = 0, op1 = 0, op2 = 0, op3 = 0;
		const InstEntry& entry = hart.decode(inst, op0, op1, op2, op3);
		if (entry.isBranch())
		  {
		    if (hart.lastPc() + instructionSize(inst) != hart.peekPc())
		      text += " (T)";
		    else
		      text += " (NT)";
		  }
		strncpy(reply.buffer, text.c_str(), sizeof(reply.buffer) - 1);
		reply.buffer[sizeof(reply.buffer) -1] = 0;
	      }
	      break;

	    case Change:
	      if (pendingChanges.empty())
		reply.type = Invalid;
	      else
		{
		  reply = pendingChanges.back();
		  pendingChanges.pop_back();
		}
	      break;

	    case Reset:
	      {
		URV addr = static_cast<URV>(msg.address);
		if (addr != msg.address)
		  std::cerr << "Error: Address too large (" << std::hex
			    << msg.address << ") in reset command.\n" << std::dec;
		pendingChanges.clear();
		if (msg.value != 0)
		  hart.defineResetPc(addr);
		hart.reset(resetMemoryMappedReg);
		if (commandLog)
		  {
		    if (msg.value != 0)
		      fprintf(commandLog, "hart=%d reset %s # ts=%s\n", hartId,
			      (boost::format(hexForm) % addr).str().c_str(),
			      timeStamp.c_str());
		    else
		      fprintf(commandLog, "hart=%d reset # ts=%s\n", hartId,
			      timeStamp.c_str());
		  }
	      }
	      break;

	    case Exception:
	      {
                std::string text;
		exceptionCommand(msg, reply, text);
		if (commandLog)
		  fprintf(commandLog, "hart=%d %s # ts=%s\n", hartId,
			  text.c_str(), timeStamp.c_str());
	      }
	      break;

	    case EnterDebug:
              {
                bool force = msg.flags;
                if (checkHart(msg, "enter_debug", reply))
                  hart.enterDebugMode(hart.peekPc(), force);
                if (commandLog)
                  fprintf(commandLog, "hart=%d enter_debug %s # ts=%s\n", hartId,
                          force? "true" : "false", timeStamp.c_str());
              }
	      break;

	    case ExitDebug:
              if (checkHart(msg, "exit_debug", reply))
                hart.exitDebugMode();
              if (commandLog)
                fprintf(commandLog, "hart=%d exit_debug # ts=%s\n", hartId,
                        timeStamp.c_str());
	      break;

	    case LoadFinished:
              {
                URV addr = static_cast<URV>(msg.address);
                unsigned tag = msg.flags;
                if (checkHart(msg, "load_finished", reply))
                  {
                    if (addr != msg.address)
                      std::cerr << "Error: Address too large (" << std::hex
                                << msg.address << ") in load finished command.\n"
                                << std::dec;
                    unsigned matchCount = 0;
                    hart.applyLoadFinished(addr, tag, matchCount);
                    reply.value = matchCount;
                  }
                if (commandLog)
                  fprintf(commandLog, "hart=%d load_finished 0x%0*" PRIx64 " %d # ts=%s\n",
                          hartId,
                          ( (sizeof(URV) == 4) ? 8 : 16 ), uint64_t(addr),
                          tag, timeStamp.c_str());
              }
              break;

            case CancelDiv:
              if (checkHart(msg, "cancel_div", reply))
                if (not hart.cancelLastDiv())
                  reply.type = Invalid;
              if (commandLog)
                fprintf(commandLog, "hart=%d cancel_div # ts=%s\n", hartId,
                        timeStamp.c_str());
              break;

            case CancelLr:
              if (checkHart(msg, "cancel_lr", reply))
                hart.cancelLr();
              if (commandLog)
                fprintf(commandLog, "hart=%d cancel_lr # ts=%s\n", hartId,
                        timeStamp.c_str());
              break;

            case DumpMemory:
              if (not system_.writeAccessedMemory(msg.buffer))
                reply.type = Invalid;
              if (commandLog)
                fprintf(commandLog, "hart=%d dump_memory %s # ts=%s\n",
                        hartId, msg.buffer, timeStamp.c_str());
              break;

	    default:
              std::cerr << "Unknown command\n";
	      reply.type = Invalid;
	    }
	}

      if (not sendMessage(soc, reply))
	return false;
    }

  return false;
}


template class WdRiscv::Server<uint32_t>;
template class WdRiscv::Server<uint64_t>;
