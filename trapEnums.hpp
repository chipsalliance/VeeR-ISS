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

namespace WdRiscv
{

  /// Privilige mode.
  enum class PrivilegeMode : uint32_t
    {
      User = 0,
      Supervisor = 1,
      Reserved = 2,
      Machine = 3
    };


  /// RISCV interrupt cause.
  enum class InterruptCause : uint32_t
    {
      U_SOFTWARE   = 0,  // User mode software interrupt
      S_SOFTWARE   = 1,  // Supervisor mode software interrupt
      M_SOFTWARE   = 3,  // Machine mode software interrupt
      U_TIMER      = 4,  // User mode timer interrupt
      S_TIMER      = 5,  // Supervisor
      M_TIMER      = 7,  // Machine
      U_EXTERNAL   = 8,  // User mode external interrupt
      S_EXTERNAL   = 9,  // Supervisor
      M_EXTERNAL   = 11, // Machine
      M_INT_TIMER1 = 28, // Internal timer 1 (WD extension) bit position.
      M_INT_TIMER0 = 29, // Internal timer 0 (WD extension) bit position.
      M_LOCAL      = 30, // Correctable error local interrupt (WD extension)
      MAX_CAUSE    = M_LOCAL
    };


  /// RISCV exception cause.
  enum class ExceptionCause : uint32_t
    {
      INST_ADDR_MISAL   = 0,  // Instruction address misaligned
      INST_ACC_FAULT    = 1,  // Instruction access fault
      ILLEGAL_INST      = 2,  // Illegal instruction
      BREAKP            = 3,  // Breakpoint
      LOAD_ADDR_MISAL   = 4,  // Load address misaligned
      LOAD_ACC_FAULT    = 5,  // Load access fault
      STORE_ADDR_MISAL  = 6,  // Store address misaligned
      STORE_ACC_FAULT   = 7,  // Store access fault.
      U_ENV_CALL        = 8,  // Environment call from user mode
      S_ENV_CALL        = 9,  // Environment call from supervisor mode
      M_ENV_CALL        = 11, // Environment call from machine mode
      INST_PAGE_FAULT   = 12, // Instruction page fault
      LOAD_PAGE_FAULT   = 13, // Load page fault
      STORE_PAGE_FAULT  = 15, // Store page fault
      NONE              = 16,
      MAX_CAUSE         = NONE
    };


  /// Non-maskable interrupt cause.
  enum class NmiCause : uint32_t
    {
      UNKNOWN               = 0,
      STORE_EXCEPTION       = 0xf0000000,
      LOAD_EXCEPTION        = 0xf0000001,
      DOUBLE_BIT_ECC        = 0xf0001000,
      DCCM_ACCESS_ERROR     = 0xf0001001,
      NON_DCCM_ACCESS_ERROR = 0xf0001002
    };


  /// Secondary exception cause values (WD special).
  enum class SecondaryCause : uint32_t
    {
      NONE = 0,

      // Cause = INST_ACC_FAULT
      INST_DOUBLE_ECC = 1,
      INST_LOCAL_UNMAPPED = 2,
      INST_MEM_PROTECTION = 3,
      INST_PMP = 8,
      INST_PRECISE = 9,  // precise bus error
      INST_OUT_OF_BOUNDS = 0xb,

      // Cause = BREAKP
      TRIGGER_HIT = 1,
      BREAKP = 2,

      // Cause = LOAD_ADDR_MISAL
      LOAD_MISAL_IO = 1,
      LOAD_MISAL_REGION_CROSS = 2,

      // Cause = LOAD_ACC_FAULT
      LOAD_ACC_DOUBLE_ECC = 1,
      LOAD_ACC_LOCAL_UNMAPPED = 2,
      LOAD_ACC_STACK_CHECK = 0xa,
      LOAD_ACC_MEM_PROTECTION = 3,
      LOAD_ACC_64BIT = 4,
      LOAD_ACC_REGION_PREDICTION = 5,
      LOAD_ACC_PIC = 6,
      LOAD_ACC_AMO = 7,
      LOAD_ACC_PMP = 8,
      LOAD_ACC_PRECISE = 9,  // precise bus error
      LOAD_ACC_OUT_OF_BOUNDS = 0xb,
      LOAD_ACC_AMO_UNCACHED = 0xc,

      // Cause = STORE_ADDR_MISAL
      STORE_MISAL_IO = 1,
      STORE_MISAL_REGION_CROSS = 0x2,

      // Cause = STORE_ACC_FAULT
      STORE_ACC_DOUBLE_ECC = 1,
      STORE_ACC_LOCAL_UNMAPPED = 2,
      STORE_ACC_MEM_PROTECTION = 3,
      STORE_ACC_64BIT = 4,
      STORE_ACC_REGION_PREDICTION = 5,
      STORE_ACC_PIC = 6,
      STORE_ACC_AMO = 7,
      STORE_ACC_PMP = 8,
      STORE_ACC_PRECISE = 9,
      STORE_ACC_STACK_CHECK = 0xa,
      STORE_ACC_OUT_OF_BOUNDS = 0xb,
      STORE_ACC_AMO_UNCACHED = 0xc,

      MAX_CAUSE = INST_OUT_OF_BOUNDS
    };


  /// Reason for entering debug mode (value stored in cause field
  /// of dcsr)
  enum class DebugModeCause
    {
      EBREAK = 1, TRIGGER = 2, DEBUGGER = 3, STEP = 4
    };

}
