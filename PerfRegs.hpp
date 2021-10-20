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

#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <unordered_map>
#include <type_traits>
#include <cassert>
#include "trapEnums.hpp"


namespace WdRiscv
{

  /// Symbolic names for performance events.
  enum class EventNumber
    {
      None,
      ClockActive,       // 1:  Cycles clock active
      ICacheHits,        // 2:  Instruction cache hits
      ICacheMisses,      // 3:  Instruction cache misses

      InstCommited,      // 4:  Instructions committed
      Inst16Commited,    // 5:  16-bit instructions committed
      Inst32Commited,    // 6:  32-bit instructions committed
      InstAligned,       // 7   4-byte aligned instructions

      InstDecode,        // 8:  Instructions decoded

      Mult,              // 9:  Multiply instructions committed
      Div,               // 10: Divide  instructions committed
      Load,              // 11: Loads committed
      Store,             // 12: stores committed
      MisalignLoad,      // 13: misaligned loads
      MisalignStore,     // 14: misaligned stores
      Alu,               // 15: alu instructions committed
      CsrRead,           // 16: Csr read instructions committed
      CsrReadWrite,      // 17: Csr read/write instructions committed
      CsrWrite,          // 18: Csr write instructions committed
      Ebreak,            // 19: Ebreak instructions committed
      Ecall,             // 20: Ecall instructions committed
      Fence,             // 21: Fence instructions committed
      Fencei,            // 22: Fence.i instructions committed
      Mret,              // 23: Mret instructions committed
      Branch,            // 24: Branch instructions committed

      BranchMiss,        // 25: Mis-predicted branches

      BranchTaken,       // 26: Taken branches

      BranchUnpredict,   // 27: Unpredictable branches
      FetchStall,        // 28: Fetcher stall cycles
      AlignStall,        // 29: Aligner stall cycles
      DecodeStall,       // 30: Decoder stall cycles
      PostSyncStall,     // 31: Post sync stall cycles
      PreSynchStall,     // 32: Pre sync stall cycles
      PipeFrozen,        // 33: Cycles pipeline is frozen
      StoreStall,        // 34: LSU store stalls cycles
      DmaDccmStall,      // 35: DMA DCCM stall cycles
      DmaIccmStall,      // 36: DMA ICCM stall cycles

      Exception,         // 37: Exception count

      TimerInterrupt,    // 38: Timer interrupts

      ExternalInterrupt, // 39: External interrupts

      TluFlush,          // 40: TLU flushes (flush lower) 
      TluFlushError,     // 41: Branch error flushes
      BusFetch,          // 42: Fetch bus transactions
      BusTransactions,   // 43: Load/store bus transactions
      BusMisalign,       // 44: Misaligned load/store bus transactions
      IbusError,         // 45: I-bus errors
      DbusError,         // 46: D-bus errors
      IbusBusy,          // 47: Cycles stalled due to Ibus busy 
      DbusBusy,          // 48: Cycles stalled due to Dbus busy 
      InetrruptDisabled, // 49: Cycles interrupts disabled 
      InterrutpStall,    // 50: Cycles interrupts stalled while disabled
      Atomic,            // 51: Atomic (amo) instruction
      Lr,                // 52: Load-reserve instruction
      Sc,                // 53: Store-conditional instruction

      Bitmanip,          // 54: Bit-manipulation
      BusLoad,           // 55: Bus load instructions committed
      BusStore,          // 56: Bus store instructions committed

      MultDiv,           // 57: M-extension instruction (Multiply/divide)
      FpHalf,            // 58: Half precision instruction
      FpSingle,          // 59: Single precision instruction
      FpDouble,          // 60: Double precision instruction
      Vector,            // 61: Vector instruction
      _End               // 62: Non-event serving as count of events
    };


  template <typename URV>
  class CsRegs;

  template <typename URV>
  class Hart;


  /// Model a set of consecutive performance counters. Theses
  /// correspond to a set of consecutive performance counter CSR.
  class PerfRegs
  {
  public:

    friend class Hart<uint32_t>;
    friend class Hart<uint64_t>;
    friend class CsRegs<uint32_t>;
    friend class CsRegs<uint64_t>;

    /// Define numCounters counters. These correspond to mhp
    PerfRegs(unsigned numCounters = 0);

    /// Configure numCounters counters initialized to zero.  This
    /// should not be used if some CSR registers are tied to the
    /// counters in here.
    void config(unsigned numCounters, unsigned maxEventId);

    /// Update (count-up) all the performance counters currently
    /// associated with the given event and enabled for the given
    /// mode.
    bool updateCounters(EventNumber event, uint32_t perfControl,
                        PrivilegeMode mode)
    {
      size_t eventIx = size_t(event);
      if (eventIx >= countersOfEvent_.size())
	return false;
      const auto& counterIndices = countersOfEvent_.at(eventIx);
      bool user = (mode == PrivilegeMode::User);
      bool machine = (mode == PrivilegeMode::Machine);
      for (auto counterIx : counterIndices)
	{
          // Performance counters handeled in here are MHPMCOUNTER3 to
          // MHPMCOUNTER31 and they are indexed 0 to 29.
          if ((perfControl >> (3+counterIx)) & 1)
            {
              bool enable = ((user and enableUser_.at(counterIx)) or
                             (machine and enableMachine_.at(counterIx)));
              if (enable)
                counters_.at(counterIx)++;
            }
	}
      return true;
    }

    /// Associate given event number with given counter.  Subsequent
    /// calls to updatePerofrmanceCounters(en) will cause given
    /// counter to count up by 1 in user mode if enableUser is true
    /// and in machine mode if enableMachine is true. Return true on
    /// success. Return false if counter number is out of bounds.
    bool assignEventToCounter(EventNumber event, unsigned counter,
                              bool enableUser, bool enableMachine)
    {
      pendingEvent_ = event;
      pendingCounter_ = counter;
      pendingUser_ = enableUser;
      pendingMachine_ = enableMachine;
      hasPending_ = true;
      return true;
    }

    /// Return the number of perormance counters.
    size_t size() const
    { return counters_.size(); }

  protected:

    bool applyPerfEventAssign();

    /// Reset all assosiations among events and counters.
    void reset();

  private:

    // Map counter index to event currently associated with counter.
    std::vector<EventNumber> eventOfCounter_;

    // Map counter index to enable flag in user mode.
    std::vector<bool> enableUser_;

    // Map counter index to enable flag in machine mode.
    std::vector<bool> enableMachine_;

    // Map an event number to a vector containing the indices of the
    // counters currently associated with that event.
    std::vector< std::vector<unsigned> > countersOfEvent_;

    std::vector<uint64_t> counters_;

    // Pending event assignment to counter.
    EventNumber pendingEvent_ = EventNumber::None;
    unsigned pendingCounter_ = 0;
    bool pendingUser_ = false;
    bool pendingMachine_ = false;
    bool hasPending_ = false;
  };
}
