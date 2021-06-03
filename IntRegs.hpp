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
#include <string>
#include <unordered_map>
#include <type_traits>
#include <cassert>

namespace WdRiscv
{

    /// Symbolic names of the integer registers.
    enum IntRegNumber
      {
	RegX0 = 0,
	RegX1 = 1,
	RegX2 = 2,
	RegX3 = 3,
	RegX4 = 4,
	RegX5 = 5,
	RegX6 = 6,
	RegX7 = 7,
	RegX8 = 8,
	RegX9 = 9,
	RegX10 = 10,
	RegX11 = 11,
	RegX12 = 12,
	RegX13 = 13,
	RegX14 = 14,
	RegX15 = 15,
	RegX16 = 16,
	RegX17 = 17,
	RegX18 = 18,
	RegX19 = 19,
	RegX20 = 20,
	RegX21 = 21,
	RegX22 = 22,
	RegX23 = 23,
	RegX24 = 24,
	RegX25 = 25,
	RegX26 = 26,
	RegX27 = 27,
	RegX28 = 28,
	RegX29 = 29,
	RegX30 = 30,
	RegX31 = 31,
	RegZero = RegX0,
	RegRa = RegX1, // return address
	RegSp = RegX2, // stack pointer
	RegGp = RegX3, // global pointer
	RegTp = RegX4, // thread pointer
	RegFp = RegX8, // frame pointer
	RegS0 = RegX8, // Callee saved registers
	RegS1 = RegX9, 
	RegA0 = RegX10, // Call arguments (caller save)
	RegA1 = RegX11,
	RegA2 = RegX12,
	RegA3 = RegX13,
	RegA4 = RegX14,
	RegA5 = RegX15,
	RegA6 = RegX16,
	RegA7 = RegX17,
	RegS2 = RegX18,  // Callee saved registers.
	RegS3 = RegX19,
	RegS4 = RegX20,
	RegS5 = RegX21,
	RegS6 = RegX22,
	RegS7 = RegX23,
	RegS8 = RegX24,
	RegS9 = RegX25,
	RegS10 = RegX26,
	RegS11 = RegX27,
	RegT0 = RegX5, // temporary
	RegT1 = RegX6,
	RegT2 = RegX7,
	RegT3 = RegX28,
	RegT4 = RegX29,
	RegT5 = RegX30,
	RegT6 = RegX31
      };


  template <typename URV>
  class Hart;

  /// Model a RISCV integer register file.
  /// URV (unsigned register value) is the register value type. For
  /// 32-bit registers, URV must be uint32_t. For 64-bit integers,
  /// it must be uint64_t.
  template <typename URV>
  class IntRegs
  {
  public:

    friend class Hart<URV>;  // To access reset and last-written-reg methods.

    /// Constructor: Define a register file with the given number of
    /// registers. Each register is of type URV. All registers initialized
    /// to zero.
    IntRegs(unsigned registerCount);

    /// Destructor.
    ~IntRegs()
    { regs_.clear(); }
    
    /// Return value of ith register. Register zero always yields zero.
    URV read(unsigned i) const
    { return regs_[i]; }

    /// Set value of ith register to the given value. Setting register
    /// zero has no effect.
    void write(unsigned i, URV value)
    {
      originalValue_ = regs_[i];
      if (i != 0)
	regs_[i] = value;
      lastWrittenReg_ = i;
    }

    /// Similar to write but does not record a change.
    void poke(unsigned i, URV value)
    {
      if (i != 0)
	regs_.at(i) = value;
    }

    /// Return the count of registers in this register file.
    size_t size() const
    { return regs_.size(); }

    /// Set ix to the number of the register corresponding to the
    /// given name returning true on success and false if no such
    /// register.  For example, if name is "x2" then ix will be set to
    /// 2. If name is "tp" then ix will be set to 4.
    bool findReg(const std::string& name, unsigned& ix) const;

    /// Return the name of the given register.
    std::string regName(unsigned i, bool abiNames = false) const
    {
      if (abiNames)
	{
	  if (i < numberToAbiName_.size())
	    return numberToAbiName_[i];
	  return std::string("x?");
	}
      if (i < numberToName_.size())
	return numberToName_[i];
      return std::string("x?");
    }

  protected:

    /// Clear all regisers.
    void reset()
    {
      clearLastWrittenReg();
      for (auto& reg : regs_)
	reg = 0;
    }

    /// Clear the number denoting the last written register.
    void clearLastWrittenReg()
    { lastWrittenReg_ = -1; }

    /// Return the number of the last written register or -1 if no register has
    /// been written since the last clearLastWrittenReg.
    int getLastWrittenReg() const
    { return lastWrittenReg_; }

    /// Set regIx and regValue to the index and previous value (before
    /// write) of the last written register returning true on success
    /// and false if no integer was written by the last executed
    /// instruction (in which case regIx and regVal are left
    /// unmodified).
    bool getLastWrittenReg(unsigned& regIx, URV& regValue) const
    {
      if (lastWrittenReg_ < 0) return false;
      regIx = lastWrittenReg_;
      regValue = originalValue_;
      return true;
    }

  private:

    std::vector<URV> regs_;
    int lastWrittenReg_ = -1;  // Register accessed in most recent write.
    URV originalValue_ = 0;    // Original value of last written reg.
    std::unordered_map<std::string, IntRegNumber> nameToNumber_;
    std::vector<std::string> numberToAbiName_;
    std::vector<std::string> numberToName_;
  };
}
