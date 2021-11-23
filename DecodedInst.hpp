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

#include <vector>
#include <string>
#include <unordered_map>
#include "InstEntry.hpp"
#include "FpRegs.hpp"
#include "InstId.hpp"


namespace WdRiscv
{

  /// Model a decoded instruction: instruction address, opcode, and
  /// operand fields. All instructions are assumed to have the form
  ///   inst op0, op1, op2, op3
  /// where op0 to op3 are optional. For example, in "add x2, x1, x0"
  /// op0 is x2, op1 is x1 and op2 is x0.
  ///
  /// Load instructions of the form "load rd, offset(rs1)" get mapped
  /// to "load rd, rs1, offset" assigning rd to op0 and offset to op2.
  ///
  /// Store instructions of the form "store rs2, offset(rs1)" get mapped
  /// to "store rs2, rs1, offset" assigning rs2 to op0 and offset to op2.
  ///   
  class DecodedInst
  {
  public:
    
    /// Default contructor: Define an invalid object.
    DecodedInst()
      : addr_(0), physAddr_(0), inst_(0), size_(0), entry_(nullptr), op0_(0),
	op1_(0), op2_(0), op3_(0), valid_(false), masked_(false), vecFields_(0)
    { values_[0] = values_[1] = values_[2] = values_[3] = 0; }

    /// Constructor.
    DecodedInst(uint64_t addr, uint32_t inst, const InstEntry* entry,
		uint32_t op0, uint32_t op1, uint32_t op2, uint32_t op3)
      : addr_(addr), physAddr_(0), inst_(inst), size_(instructionSize(inst)),
	entry_(entry), op0_(op0), op1_(op1), op2_(op2), op3_(op3),
	valid_(entry != nullptr), masked_(false), vecFields_(0)
    { values_[0] = values_[1] = values_[2] = values_[3] = 0; }

    /// Return instruction size in bytes.
    uint32_t instSize() const
    { return size_; }

    /// Return the virtual address of the instruction.
    uint64_t address() const
    { return addr_; }

    /// Return the physical address of the instruction.
    uint64_t physAddress() const
    { return physAddr_; }

    /// Return instruction code.
    uint32_t inst() const
    { return inst_; }

    /// Return the 1st operand (zero if instruction has no operands).
    /// First operand is typically the destination register.
    uint32_t op0() const
    { return op0_; }

    /// Return 2nd operand (zero if instruction has no 2nd operand).
    /// Second operand is typically source register rs1.
    uint32_t op1() const
    { return op1_; }

    /// Return 2nd operand as a signed integer. This is useful
    /// for instructions where the 2nd operand is a signed immediate
    /// value.
    template <typename SI>
    SI op1As() const
    { return int32_t(op1_); }

    /// Return 3rd operand (zero if instruction has no 3rd operand).
    /// Third operand is typically source register rs2 or immediate
    /// value.
    uint32_t op2() const
    { return op2_; }

    /// Return 3rd operand as a signed 32-bit integer. This is useful
    /// for instructions where the 3rd operand is a signed immediate
    /// value.
    template <typename SI>
    SI op2As() const
    { return int32_t(op2_); }

    /// Return 4th operand (zero if instruction has no 4th operand).
    /// Fourth operand is typically source register rs3 for
    /// multiply-add like floating point instructions.
    uint32_t op3() const
    { return op3_; }

    /// Return the operand count associated with this
    /// instruction. Immediate values are counted as operands. For
    /// example, in "addi x3, x4, 10", there are 3 operands: 3, 4, and
    /// 10 with types IntReg, IntReg and Imm respectively.
    unsigned operandCount() const
    { return isValid()? entry_->operandCount() : 0; }

    /// Return the ith operands or zero if i is out of bounds. For exmaple, if
    /// decode insruction is "addi x3, x4, 10" then the 0th operand would be 3
    /// and the second operands would be 10.
    uint32_t ithOperand(unsigned i) const;

    /// Return the ith operands or zero if i is out of bounds. For exmaple, if
    /// decode insruction is "addi x3, x4, 10" then the 0th operand would be 3
    /// and the second operands would be 10.
    int32_t ithOperandAsInt(unsigned i) const;

    /// Return the type of the ith operand or None if i is out of
    /// bounds. Object must be valid.
    OperandType ithOperandType(unsigned i) const
    { return isValid()? entry_->ithOperandType(i) : OperandType::None; }

    /// Return true if this object is valid.
    bool isValid() const
    { return valid_; }

    /// Make invalid.
    void invalidate()
    { valid_ = false; }

    /// Return associated instruction table information.
    const InstEntry* instEntry() const
    { return entry_; }

    /// Relevant for floating point instructions with rounding mode.
    RoundingMode roundingMode() const
    { return RoundingMode((inst_ >> 12) & 7); }

    /// Relevant to atomic instructions: Return true if acquire bit is
    /// set.
    bool isAtomicAcquire() const
    { return (inst_ >> 26) & 1; }

    /// Relevant to atomic instructions: Return true if release bit is
    /// set.
    bool isAtomicRelease() const
    { return (inst_ >> 25) & 1; }

    /// Associate a value with each operand by fetching
    /// registers. After this method, the value of an immediate
    /// operand x is x. The value of register operand y is the value
    /// currently stored in register x. The value of a non-existing
    /// operand is zero. Note that the association is only in this
    /// object and that no register value is changed by this method.
    template <typename URV>
    void fetchOperands(const Hart<URV>& hart);

    /// Associated a value with the ith operand. This has no effect if
    /// i is out of bounds or if the ith operand is an immediate. Note
    /// that the association is only in this object and that no
    /// register value is changed by this method.
    void setIthOperandValue(unsigned i, uint64_t value);

    /// Return value associated with ith operand.
    uint64_t ithOperandValue(unsigned i) const
    { return i < 4? values_[i] : 0; }

    /// Return true if this is a vector instruction with masking enabled, 
    bool isMasked() const
    { return masked_; }

    /// Return number of fields in vector ld/st instruction. Return zero
    /// if this is not a vector ld/st.
    unsigned vecFieldCount() const
    { return vecFields_; }

  protected:

    friend class Hart<uint32_t>;
    friend class Hart<uint64_t>;

    void setAddr(uint64_t addr)
    { addr_ = addr; }

    void setInst(uint32_t inst)
    { inst_ = inst; size_ = instructionSize(inst); }

    void setEntry(const InstEntry* e)
    { entry_ = e; if (not e) valid_ = false; }

    void setOp0(uint32_t op0)
    { op0_ = op0; }

    void setOp1(uint32_t op1)
    { op1_ = op1; }

    void setOp2(uint32_t op2)
    { op2_ = op2; }

    void setOp3(uint32_t op3)
    { op3_ = op3; }

    void setMasked(bool flag)
    { masked_ = flag; }

    void setVecFieldCount(uint32_t count)
    { vecFields_ = count; }

    void reset(uint64_t addr, uint64_t physAddr, uint32_t inst,
	       const InstEntry* entry,
	       uint32_t op0, uint32_t op1, uint32_t op2, uint32_t op3)
    {
      addr_ = addr;
      physAddr_ = physAddr;
      inst_ = inst;
      entry_ = entry;
      op0_ = op0; op1_ = op1; op2_ = op2; op3_ = op3;
      size_ = instructionSize(inst);
      valid_ = entry != nullptr;
    }

  private:

    uint64_t addr_;
    uint64_t physAddr_;
    uint32_t inst_;
    uint32_t size_;
    const InstEntry* entry_;
    uint32_t op0_;    // 1st operand (typically a register number)
    uint32_t op1_;    // 2nd operand (register number or immediate value)
    uint32_t op2_;    // 3rd operand (register number or immediate value)
    uint32_t op3_;    // 4th operand (typically a register number)

    uint64_t values_[4];  // Values of operands.
    bool valid_;
    bool masked_;     // For vector instructions.
    uint8_t vecFields_;   // For vector ld/st instructions.
  };


  /// Return 2nd operand as a signed 64-bit integer. This is useful
  /// for instructions where the 2nd operand is a signed immediate
  /// value.
  template <>
  inline
  int64_t
  DecodedInst::op1As<int64_t>() const
  { return int32_t(op1_); }


  /// Return 3rd operand as a signed 64-bit integer. This is useful
  /// for instructions where the 3rd operand is a signed immediate
  /// value.
  template <>
  inline
  int64_t
  DecodedInst::op2As<int64_t>() const
  { return int32_t(op2_); }

}
