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

#include <climits>
#include <mutex>
#include <atomic>
#include <cassert>

#include "instforms.hpp"
#include "DecodedInst.hpp"
#include "Hart.hpp"
#include "System.hpp"

using namespace WdRiscv;


template <typename URV>
ExceptionCause
Hart<URV>::validateAmoAddr(uint32_t rs1, uint64_t& addr, unsigned accessSize,
                           SecondaryCause& secCause)
{
  URV mask = URV(accessSize) - 1;

  auto cause = ExceptionCause::NONE;
  bool forcedFail = false;
  if (accessSize == 4)
    {
      uint32_t storeVal = 0;
      cause = determineStoreException(rs1, addr, addr, storeVal, secCause, forcedFail);
    }
  else
    {
      uint64_t storeVal = 0;
      cause = determineStoreException(rs1, addr, addr, storeVal, secCause, forcedFail);
    }

  if (cause == ExceptionCause::STORE_ADDR_MISAL and
      misalAtomicCauseAccessFault_)
    {
      cause = ExceptionCause::STORE_ACC_FAULT;
      secCause = SecondaryCause::STORE_ACC_AMO;
    }

  // Check if invalid unless cacheable.
  if (amoInCacheableOnly_ and not isAddrCacheable(addr))
    if (cause == ExceptionCause::NONE)
      {
        cause = ExceptionCause::STORE_ACC_FAULT;
        secCause = SecondaryCause::STORE_ACC_AMO_UNCACHED;
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
Hart<URV>::amoLoad32(uint32_t rd, uint32_t rs1, uint32_t rs2, URV& value)
{
  URV virtAddr = intRegs_.read(rs1);

  ldStAddr_ = virtAddr;   // For reporting load addr in trace-mode.
  ldStPhysAddr_ = ldStAddr_;
  ldStAddrValid_ = true;  // For reporting load addr in trace-mode.

  if (loadQueueEnabled_)
    {
      removeFromLoadQueue(rs1, false);
      removeFromLoadQueue(rs2, false);
    }

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
  ldStPhysAddr_ = addr;

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

      URV prevRdVal = peekIntReg(rd);
      putInLoadQueue(4 /*ldSize*/, addr, rd, prevRdVal);

      return true;  // Success.
    }

  cause = ExceptionCause::STORE_ACC_FAULT;
  initiateLoadException(cause, virtAddr, secCause);
  return false;
}


template <typename URV>
bool
Hart<URV>::amoLoad64(uint32_t rd, uint32_t rs1, uint32_t rs2, URV& value)
{
  URV virtAddr = intRegs_.read(rs1);

  ldStAddr_ = virtAddr;   // For reporting load addr in trace-mode.
  ldStPhysAddr_ = ldStAddr_;
  ldStAddrValid_ = true;  // For reporting load addr in trace-mode.

  if (loadQueueEnabled_)
    {
      removeFromLoadQueue(rs1, false);
      removeFromLoadQueue(rs2, false);
    }

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
  ldStPhysAddr_ = addr;

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

      URV prevRdVal = peekIntReg(rd);
      putInLoadQueue(8 /*ldSize*/, addr, rd, prevRdVal);

      return true;  // Success.
    }

  cause = ExceptionCause::STORE_ACC_FAULT;
  initiateLoadException(cause, virtAddr, secCause);
  return false;
}


template <typename URV>
template <typename LOAD_TYPE>
bool
Hart<URV>::loadReserve(uint32_t rd, uint32_t rs1, uint64_t& physAddr)
{
  URV virtAddr = intRegs_.read(rs1);

  ldStAddr_ = virtAddr;   // For reporting load addr in trace-mode.
  ldStPhysAddr_ = ldStAddr_;
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
  ldStPhysAddr_ = addr;

  // Check if invalid unless cacheable.
  if (amoInCacheableOnly_ and not isAddrCacheable(addr))
    if (cause == ExceptionCause::NONE)
      {
        cause = ExceptionCause::LOAD_ACC_FAULT;
        secCause = SecondaryCause::LOAD_ACC_AMO_UNCACHED;
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
  if (not isRva())
    {
      illegalInst(di);
      return;
    }

  std::lock_guard<std::mutex> lock(memory_.lrMutex_);

  lrCount_++;
  uint64_t physAddr = 0;
  if (not loadReserve<int32_t>(di->op0(), di->op1(), physAddr))
    return;
  memory_.makeLr(hartIx_, physAddr, 4 /*size*/);
  lrSuccess_++;
}


/// STORE_TYPE is either uint32_t or uint64_t.
template <typename URV>
template <typename STORE_TYPE>
bool
Hart<URV>::storeConditional(uint32_t rs1, URV virtAddr, STORE_TYPE storeVal)
{
  ldStAddr_ = virtAddr;   // For reporting ld/st addr in trace-mode.
  ldStPhysAddr_ = ldStAddr_;
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
  bool forcedFail = false;
  auto cause = determineStoreException(rs1, virtAddr, addr, storeVal, secCause,
                                       forcedFail);

  if (cause == ExceptionCause::STORE_ADDR_MISAL and
      misalAtomicCauseAccessFault_)
    {
      cause = ExceptionCause::STORE_ACC_FAULT;
      secCause = SecondaryCause::STORE_ACC_AMO;
    }
  ldStPhysAddr_ = addr;

  // Check if invalid unless cacheable.
  if (amoInCacheableOnly_ and not isAddrCacheable(addr))
    if (cause == ExceptionCause::NONE)
      {
        cause = ExceptionCause::STORE_ACC_FAULT;
        secCause = SecondaryCause::STORE_ACC_AMO_UNCACHED;
      }

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

  if (not memory_.hasLr(hartIx_, addr, sizeof(storeVal)))
    return false;

  if (memory_.write(hartIx_, addr, storeVal))
    {
      invalidateDecodeCache(virtAddr, sizeof(STORE_TYPE));

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
  if (not isRva())
    {
      illegalInst(di);
      return;
    }

  std::lock_guard<std::mutex> lock(memory_.lrMutex_);

  uint32_t rd = di->op0(), rs1 = di->op1();
  URV value = intRegs_.read(di->op2());
  URV addr = intRegs_.read(rs1);
  scCount_++;

  if (loadQueueEnabled_)
    removeFromLoadQueue(rs1, false);

  uint64_t prevCount = exceptionCount_;

  bool ok = storeConditional(rs1, addr, uint32_t(value));
  cancelLr(); // Clear LR reservation (if any).

  if (ok)
    {
      if (loadQueueEnabled_)
        putInLoadQueue(4 /* size*/, addr, rd, peekIntReg(rd));

      memory_.invalidateOtherHartLr(hartIx_, addr, 4);
      intRegs_.write(rd, 0); // success
      scSuccess_++;

      return;
    }

  // If exception or trigger tripped then rd is not modified.
  if (triggerTripped_ or exceptionCount_ != prevCount)
    return;

  if (loadQueueEnabled_)
    putInLoadQueue(4 /* size*/, addr, rd, peekIntReg(rd));

  intRegs_.write(di->op0(), 1);  // fail
}


template <typename URV>
void
Hart<URV>::execAmoadd_w(const DecodedInst* di)
{
  if (not isRva())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad32(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(rs2);
      URV result = rs2Val + rdVal;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoswap_w(const DecodedInst* di)
{
  if (not isRva())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad32(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(rs2);
      URV result = rs2Val;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoxor_w(const DecodedInst* di)
{
  if (not isRva())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad32(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(rs2);
      URV result = rs2Val ^ rdVal;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoor_w(const DecodedInst* di)
{
  if (not isRva())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad32(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(rs2);
      URV result = rs2Val | rdVal;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoand_w(const DecodedInst* di)
{
  if (not isRva())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad32(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(rs2);
      URV result = rs2Val & rdVal;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmomin_w(const DecodedInst* di)
{
  if (not isRva())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad32(rd, rs1, rs2, loadedValue);

  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);
      URV rs2Val = intRegs_.read(rs2);

      int32_t w1 = int32_t(rs2Val);
      int32_t w2 = int32_t(loadedValue);

      int32_t result = (w1 < w2)? w1 : w2;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, loadedValue);
    }
}


template <typename URV>
void
Hart<URV>::execAmominu_w(const DecodedInst* di)
{
  if (not isRva())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad32(rd, rs1, rs2, loadedValue);

  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rs2Val = intRegs_.read(rs2);

      uint32_t w1 = uint32_t(rs2Val);
      uint32_t w2 = uint32_t(loadedValue);

      uint32_t result = (w1 < w2)? w1 : w2;

      bool storeOk = store<uint32_t>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, loadedValue);
    }
}


template <typename URV>
void
Hart<URV>::execAmomax_w(const DecodedInst* di)
{
  if (not isRva())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad32(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rs2Val = intRegs_.read(rs2);

      int32_t w1 = int32_t(rs2Val);
      int32_t w2 = int32_t(loadedValue);

      int32_t result = w1 > w2 ? w1 : w2;

      bool storeOk = store<uint32_t>(rs1, addr, addr, uint32_t(result));

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, loadedValue);
    }
}


template <typename URV>
void
Hart<URV>::execAmomaxu_w(const DecodedInst* di)
{
  if (not isRva())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad32(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      // Sign extend least significant word of register value.
      SRV rdVal = SRV(int32_t(loadedValue));

      URV rs2Val = intRegs_.read(rs2);

      uint32_t w1 = uint32_t(rs2Val);
      uint32_t w2 = uint32_t(rdVal);

      URV result = (w1 > w2)? w1 : w2;

      bool storeOk = store<uint32_t>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execLr_d(const DecodedInst* di)
{
  if (not isRva() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  std::lock_guard<std::mutex> lock(memory_.lrMutex_);

  lrCount_++;
  uint64_t physAddr = 0;
  if (not loadReserve<int64_t>(di->op0(), di->op1(), physAddr))
    return;
  memory_.makeLr(hartIx_, physAddr, 8 /*size*/);
  lrSuccess_++;
}


template <typename URV>
void
Hart<URV>::execSc_d(const DecodedInst* di)
{
  if (not isRva() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  std::lock_guard<std::mutex> lock(memory_.lrMutex_);

  uint32_t rd = di->op0(), rs1 = di->op1();
  URV value = intRegs_.read(di->op2());
  URV addr = intRegs_.read(rs1);
  scCount_++;

  if (loadQueueEnabled_)
    removeFromLoadQueue(rs1, false);

  uint64_t prevCount = exceptionCount_;

  bool ok = storeConditional(rs1, addr, uint64_t(value));
  cancelLr(); // Clear LR reservation (if any).

  if (ok)
    {
      if (loadQueueEnabled_)
        putInLoadQueue(8 /* size*/, addr, rd, peekIntReg(rd));

      memory_.invalidateOtherHartLr(hartIx_, addr, 8);
      intRegs_.write(rd, 0); // success
      scSuccess_++;

      return;
    }

  // If exception or trigger tripped then rd is not modified.
  if (triggerTripped_ or exceptionCount_ != prevCount)
    return;

  if (loadQueueEnabled_)
    putInLoadQueue(8 /* size*/, addr, rd, peekIntReg(rd));

  intRegs_.write(di->op0(), 1);  // fail
}


template <typename URV>
void
Hart<URV>::execAmoadd_d(const DecodedInst* di)
{
  if (not isRva() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad64(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(rs2);
      URV result = rs2Val + rdVal;

      bool storeOk = store<uint64_t>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoswap_d(const DecodedInst* di)
{
  if (not isRva() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad64(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(rs2);
      URV result = rs2Val;

      bool storeOk = store<uint64_t>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoxor_d(const DecodedInst* di)
{
  if (not isRva() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad64(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(rs2);
      URV result = rs2Val ^ rdVal;

      bool storeOk = store<uint64_t>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoor_d(const DecodedInst* di)
{
  if (not isRva() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad64(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(rs2);
      URV result = rs2Val | rdVal;

      bool storeOk = store<uint64_t>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmoand_d(const DecodedInst* di)
{
  if (not isRva() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad64(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(rs2);
      URV result = rs2Val & rdVal;

      bool storeOk = store<uint64_t>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmomin_d(const DecodedInst* di)
{
  if (not isRva() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad64(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(rs2);
      URV result = (SRV(rs2Val) < SRV(rdVal))? rs2Val : rdVal;

      bool storeOk = store<uint64_t>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmominu_d(const DecodedInst* di)
{
  if (not isRva() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad64(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(rs2);
      URV result = (rs2Val < rdVal)? rs2Val : rdVal;

      bool storeOk = store<uint64_t>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmomax_d(const DecodedInst* di)
{
  if (not isRva() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad64(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(rs2);
      URV result = (SRV(rs2Val) > SRV(rdVal))? rs2Val : rdVal;

      bool storeOk = store<uint64_t>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template <typename URV>
void
Hart<URV>::execAmomaxu_d(const DecodedInst* di)
{
  if (not isRva() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  // Lock mutex to serialize AMO instructions. Unlock automatically on
  // exit from this scope.
  std::lock_guard<std::mutex> lock(memory_.amoMutex_);

  URV loadedValue = 0;
  uint32_t rd = di->op0(), rs1 = di->op1(), rs2 = di->op2();
  bool loadOk = amoLoad64(rd, rs1, rs2, loadedValue);
  if (loadOk)
    {
      URV addr = intRegs_.read(rs1);

      URV rdVal = loadedValue;
      URV rs2Val = intRegs_.read(rs2);
      URV result = (rs2Val > rdVal)? rs2Val : rdVal;

      bool storeOk = store<uint64_t>(rs1, addr, addr, result);

      if (storeOk and not triggerTripped_)
	intRegs_.write(rd, rdVal);
    }
}


template class WdRiscv::Hart<uint32_t>;
template class WdRiscv::Hart<uint64_t>;
