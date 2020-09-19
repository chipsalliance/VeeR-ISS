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
#include <cassert>

namespace WdRiscv
{

  enum class VecGroupMultiplier : uint32_t
    {
     One      = 0,
     Two      = 1,
     Four     = 2,
     Eight    = 3,
     Reserved = 4,
     Eighth   = 5,
     Quarter  = 6,
     Half     = 7
    };


  /// Selected element width.
  enum class ElementWidth : uint32_t
    {
     Byte         = 0,
     HalfWord     = 1,
     Word         = 2,
     DoubleWord   = 3,
     QuadWord     = 4,
     OctWord      = 5,
     HalfKbits    = 6,
     Kbits        = 7
    };


  /// Symbolic names of the integer registers.
  enum VecRegNumber
    {
     RegV0 = 0,
     RegV1 = 1,
     RegV2 = 2,
     RegV3 = 3,
     RegV4 = 4,
     RegV5 = 5,
     RegV6 = 6,
     RegV7 = 7,
     RegV8 = 8,
     RegV9 = 9,
     RegV10 = 10,
     RegV11 = 11,
     RegV12 = 12,
     RegV13 = 13,
     RegV14 = 14,
     RegV15 = 15,
     RegV16 = 16,
     RegV17 = 17,
     RegV18 = 18,
     RegV19 = 19,
     RegV20 = 20,
     RegV21 = 21,
     RegV22 = 22,
     RegV23 = 23,
     RegV24 = 24,
     RegV25 = 25,
     RegV26 = 26,
     RegV27 = 27,
     RegV28 = 28,
     RegV29 = 29,
     RegV30 = 30,
     RegV31 = 31
    };


  template <typename URV>
  class Hart;

  /// Model a RISCV vector register file.
  class VecRegs
  {
  public:

    friend class Hart<uint32_t>;
    friend class Hart<uint64_t>;

    /// Constructor: Define an empty vector regidter file which may be
    /// reconfigured later using the config method.
    VecRegs();

    /// Destructor.
    ~VecRegs();
    
    /// Return vector register count.
    unsigned registerCount() const
    { return regCount_; }

    /// Return the number of bytes per register. This is independent
    /// of group multiplier.
    unsigned bytesPerRegister() const
    { return bytesPerReg_; }

    /// Set value to that of the element with given index within the
    /// vector register of the given number returning true on sucess
    /// and false if the combination of element index, vector number
    /// and group multipier is invalid.
    template<typename T>
    bool read(unsigned regNum, unsigned elementIx, unsigned groupMultiplier,
              T& value) const
    {
      if ((elementIx + 1) * sizeof(T) > bytesPerReg_*groupMultiplier)
        return false;
      if (regNum*bytesPerReg_ + (elementIx + 1)*sizeof(T) > bytesPerReg_*regCount_)
        return false;
      const T* data = reinterpret_cast<const T*>(data_ + regNum*bytesPerReg_);
      value = data[elementIx];
      return true;
    }

    template<typename T>
    bool write(unsigned regNum, unsigned elementIx, unsigned groupMultiplier,
               const T& value)
    {
      if ((elementIx + 1) * sizeof(T) > bytesPerReg_*groupMultiplier)
        return false;
      if (regNum*bytesPerReg_ + (elementIx + 1)*sizeof(T) > bytesPerReg_*regCount_)
        return false;
      T* data = reinterpret_cast<T*>(data_ + regNum*bytesPerReg_);
      data[elementIx] = value;
      return true;
    }

    /// Return the count of registers in this register file.
    size_t size() const
    { return regCount_; }

    /// Return the number of bits in a register in this register file.
    uint32_t bitsPerReg() const
    { return 8*bytesPerReg_; }

    /// Convert the given symbolic element width to a byte count.
    static uint32_t elementWidthInBytes(ElementWidth sew)
    { return uint32_t(1) << uint32_t(sew); }

    /// Convert the given symbolic group multiplier to a number scaled by
    /// eight (e.g. One is converted to 8, and Eighth to 1). Return 0 if
    /// given symbolic multiplier is not valid.
    static uint32_t groupMultiplierX8(VecGroupMultiplier vm)
    {
      if (vm < VecGroupMultiplier::Reserved)
        return 8*(uint32_t(1) << uint32_t(vm));
      if (vm > VecGroupMultiplier::Reserved)
        return 8 >> (8 - unsigned(vm));
      return 0;
    }

  protected:

    /// It is convenient to contruct an empty regiter file (bytesPerReg = 0)
    /// and configure it later. Old configuration is lost. Register of
    /// newly configured file are initlaized to zero.
    void config(unsigned bytesPerRegm, unsigned maxBytesPerElem);

    void reset();

    unsigned groupMultiplier() const
    { return group_; }

    unsigned startIndex() const
    { return start_; }

    unsigned elemCount() const
    { return elems_; }

    ElementWidth elemWidth() const
    { return sew_; }

    // Return true if current vtype configuration is legal. This is a cached
    // value of VTYPE.VILL.
    bool legalConfig() const
    { return not vill_; }

  private:

    unsigned regCount_ = 0;
    unsigned bytesPerReg_ = 0;
    unsigned bytesPerElem_ = 0;
    unsigned bytesInRegFile_ = 0;
    uint8_t* data_ = nullptr;

    unsigned group_ = 1;                    // Group multiplier -- cached VTYPE.VLMUL
    unsigned start_ = 0;                    // Cached VSTART
    unsigned elems_ = 0;                    // Cached VL
    ElementWidth sew_ = ElementWidth::Byte; // Cached VTYPE.SEW
    bool vill_ = false;                     // Cached VTYPE.VILL
  };
}
