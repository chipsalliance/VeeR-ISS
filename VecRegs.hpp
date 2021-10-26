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

  enum class GroupMultiplier : uint32_t
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
     Byte      = 0,
     Half      = 1,
     Word      = 2,
     Word2     = 3,
     Word4     = 4,
     Word8     = 5,
     Word16    = 6,
     Word32    = 7
    };


  enum class VecRoundingMode
    {
     NearestUp   = 0,
     NearestEven = 1,
     Down        = 2,
     Odd         = 3,
     VcsrMask    = 6,  // Mask of rounding mode bits in VCSR
     VcsrShift   = 1   // Index of least-sig rounding mode bit in VCSR
    };


  enum class VecEnums : uint32_t
    {
     GroupLimit = 8, // One past largest VecGroupMultiplier value
     WidthLimit = 8  // One past largest ElementWidth value
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
    uint32_t registerCount() const
    { return regCount_; }

    /// Return the number of bytes per register. This is independent
    /// of group multiplier.
    uint32_t bytesPerRegister() const
    { return bytesPerReg_; }

    /// Set value to that of the element with given index within the
    /// vector register of the given number returning true on sucess
    /// and false if the combination of element index, vector number
    /// and group multipier (presecaled by 8) is invalid. We pre-scale
    /// the group multiplier to avoid passing a fraction.
    template<typename T>
    bool read(uint32_t regNum, uint32_t elemIx, uint32_t groupX8,
              T& value) const
    {
      if (elemIx*sizeof(T) > ((bytesPerReg_*groupX8) >> 3) - sizeof(T))
        return false;
      if (regNum*bytesPerReg_ + elemIx*sizeof(T) > bytesPerReg_*regCount_ - sizeof(T))
        return false;
      const T* data = reinterpret_cast<const T*>(data_ + regNum*bytesPerReg_);
      value = data[elemIx];
      return true;
    }

    /// Set the element with given index within the vector register of
    /// the given number to the given value returning true on success
    /// and false if the combination of element index, vector number
    /// and group multipier (presecaled by 8) is invalid. We pre-scale
    /// the group multiplier to avoid passing a fraction.
    template<typename T>
    bool write(uint32_t regNum, uint32_t elemIx, uint32_t groupX8,
               const T& value)
    {
      if ((elemIx + 1) * sizeof(T) > ((bytesPerReg_*groupX8) >> 3))
        return false;
      if (regNum*bytesPerReg_ + (elemIx + 1)*sizeof(T) > bytesPerReg_*regCount_)
        return false;
      T* data = reinterpret_cast<T*>(data_ + regNum*bytesPerReg_);
      data[elemIx] = value;
      lastWrittenReg_ = regNum;
      lastGroupX8_ = groupX8;
      return true;
    }

    /// Read offset for a load-indexed/sore-indexed instruction.
    bool readIndex(uint32_t regNum, uint32_t elemIx, ElementWidth eew,
                   uint32_t groupX8, uint64_t& offset) const
    {
      switch(eew)
        {
        case ElementWidth::Byte:
          {
            uint8_t temp = 0;
            if (not read(regNum, elemIx, groupX8, temp))
              return false;
            offset = temp;
            return true;
          }
        case ElementWidth::Half:
          {
            uint16_t temp = 0;
            if (not read(regNum, elemIx, groupX8, temp))
              return false;
            offset = temp;
            return true;
          }
        case ElementWidth::Word:
          {
            uint32_t temp = 0;
            if (not read(regNum, elemIx, groupX8, temp))
              return false;
            offset = temp;
            return true;
          }
        case ElementWidth::Word2:
          {
            uint64_t temp = 0;
            if (not read(regNum, elemIx, groupX8, temp))
              return false;
            offset = temp;
            return true;
          }
        default:
          return false;
        }
      return false;
    }

    /// Return the count of registers in this register file.
    size_t size() const
    { return regCount_; }

    /// Return the number of bits in a register in this register file.
    uint32_t bitsPerRegister() const
    { return 8*bytesPerReg_; }

    /// Return the currently configured element width.
    ElementWidth elemWidth() const
    { return sew_; }

    /// Return the current configured group multiplier.
    GroupMultiplier groupMultiplier() const
    { return group_; }

    /// Return the currently configured element width in bits (for example
    /// if SEW is Byte, then this reurns 8).
    uint32_t elementWidthInBits() const
    { return sewInBits_; }

    /// Return the element width in bit given the symbolic element width.
    /// Rturn 0 if symbolic value is out of bounds.
    static uint32_t elementWidthInBits(ElementWidth ew)
    { return ew > ElementWidth::Word32 ? 0 : uint32_t(8) << uint32_t(ew); }

    /// Return the currently configured group multiplier as a unsigned
    /// integer scaled by 8.  For example if group multiplier is One,
    /// this reurns 8. If group multiplier is Eigth, this returns 1.
    /// We pre-scale by 8 to avoid division when the multiplier is a
    /// fraction.
    uint32_t groupMultiplierX8() const
    { return groupX8_; }

    /// Return true if double the given element width (eew=2*sew) is
    /// legal with the given group multiplier (prescaled by 8).
    bool isDoubleWideLegal(ElementWidth sew, uint32_t groupX8) const
    {
      uint32_t wideGroup = groupX8 * 2;
      GroupMultiplier emul = GroupMultiplier::One;
      if (not groupNumberX8ToSymbol(wideGroup, emul))
        return false;

      ElementWidth eew = sew;
      if (not doubleSew(sew, eew))
        return false;

      return legalConfig(eew, emul);
    }

    /// Set symbol to the symbolic value of the given numeric group
    /// multiplier (premultiplied by 8). Return true on success and
    /// false if groupX8 is out of bounds.
    static inline
    bool groupNumberX8ToSymbol(uint32_t groupX8, GroupMultiplier& symbol)
    {
      if (groupX8 == 1)  { symbol = GroupMultiplier::Eighth;   return true; }
      if (groupX8 == 2)  { symbol = GroupMultiplier::Quarter;  return true; }
      if (groupX8 == 4)  { symbol = GroupMultiplier::Half;     return true; }
      if (groupX8 == 8)  { symbol = GroupMultiplier::One;      return true; }
      if (groupX8 == 16) { symbol = GroupMultiplier::Two;      return true; }
      if (groupX8 == 32) { symbol = GroupMultiplier::Four;     return true; }
      if (groupX8 == 64) { symbol = GroupMultiplier::Eight;    return true; }
      return false;
    }

    /// Set dsew to the double of the given sew returning true on
    /// succes and false if given sew cannot be doubled.  false if
    /// groupX8 is out of bounds.
    static
    bool doubleSew(ElementWidth sew, ElementWidth& dsew)
    {
      typedef ElementWidth EW;
      if (sew == EW::Byte   ) { dsew = EW:: Half;   return true; }
      if (sew == EW::Half   ) { dsew = EW:: Word;   return true; }
      if (sew == EW::Word   ) { dsew = EW:: Word2;  return true; }
      if (sew == EW::Word2  ) { dsew = EW:: Word4;  return true; }
      if (sew == EW::Word4  ) { dsew = EW:: Word8;  return true; }
      if (sew == EW::Word8  ) { dsew = EW:: Word16; return true; }
      if (sew == EW::Word16 ) { dsew = EW:: Word32; return true; }
      return false;
    }

    /// Convert the given symbolic element width to a byte count.
    static uint32_t elementWidthInBytes(ElementWidth sew)
    { return uint32_t(1) << uint32_t(sew); }

    /// Convert the given symbolic group multiplier to a number scaled by
    /// eight (e.g. One is converted to 8, and Eigth to 1). Return 0 if
    /// given symbolic multiplier is not valid.
    static uint32_t groupMultiplierX8(GroupMultiplier gm)
    {
      if (gm < GroupMultiplier::Reserved)
        return 8*(uint32_t(1) << uint32_t(gm));
      if (gm > GroupMultiplier::Reserved and gm <= GroupMultiplier::Half)
        return 8 >> (8 - uint32_t(gm));
      return 0;
    }

    static std::string to_string(GroupMultiplier group)
    {
      static std::vector<std::string> vec =
        {"m1", "m2", "m4", "m8", "m?", "mf8", "mf4", "mf2"};
      return size_t(group) < vec.size()? vec.at(size_t(group)) : "m?";
    }

    static std::string to_string(ElementWidth ew)
    {
      static std::vector<std::string> vec =
        {"e8", "e16", "e32", "e64", "e128", "e256", "e512", "e1024"};
      return size_t(ew) < vec.size()? vec.at(size_t(ew)) : "e?";
    }


  protected:

    /// Clear load/address and store data used for logging/tracing./
    void clearTraceData()
    {
      ldStAddr_.clear();
      stData_.clear();
      clearLastWrittenReg();
      opsEmul_.assign(opsEmul_.size(), 1);
    }

    /// Clear the number denoting the last written register.
    void clearLastWrittenReg()
    { lastWrittenReg_ = -1; }
    
    /// Return the number of the last written vector regsiter or -1 if no
    /// no register has been written since the last clearLastWrittenReg.
    int getLastWrittenReg(uint32_t& groupX8) const
    {
      if (lastWrittenReg_ >= 0)
        {
          groupX8 = lastGroupX8_;
          return lastWrittenReg_;
        }
      return -1;
    }

    /// For instructions that do not use the write method, mark the
    /// last written register and the effective element widht.
    void touchReg(uint32_t reg, uint32_t groupX8)
    { lastWrittenReg_ = reg; lastGroupX8_ = groupX8; }

    /// Same as above for mask registers
    void touchMask(uint32_t reg)
    { touchReg(reg, 8); }  // Grouping of of 1

    /// Return true if element of given index is active with respect
    /// to the given mask vector register. Element is active if the
    /// corresponding mask bit is 1.
    bool isActive(uint32_t maskReg, uint32_t ix) const
    {
      if (maskReg >= regCount_)
        return false;

      uint32_t byteIx = ix >> 3;
      uint32_t bitIx = ix & 7;  // bit in byte
      if (byteIx >= bytesPerReg_)
        return false;

      const uint8_t* data = data_ + maskReg*bytesPerReg_;
      return (data[byteIx] >> bitIx) & 1;
    }

    /// Set the ith bit of the given mask regiser to the given value.
    /// Return true on success and false on failure (register or
    /// element index out of bound.
    bool writeMaskRegister(uint32_t maskReg, uint32_t i, bool value)
    {
      if (maskReg >= regCount_)
        return false;

      uint32_t byteIx = i >> 3;
      uint32_t bitIx = i & 7;  // bit in byte
      if (byteIx >= bytesPerReg_)
        return false;

      uint8_t* data = data_ + maskReg*bytesPerReg_;
      uint8_t mask = uint8_t(1) << bitIx;
      if (value)
        data[byteIx] |= mask;
      else
        data[byteIx] &= ~mask;
      lastWrittenReg_ = maskReg;
      lastGroupX8_ = 8;
      return true;
    }

    /// Return the pointers to the 1st byte of the memory area
    /// associated with the given vector. Return nullptr if
    /// vector index is out of bounds.
    uint8_t* getVecData(uint32_t vecIx)
    {
      if (vecIx >= regCount_)
        return nullptr;
      return data_ + vecIx*bytesPerReg_;
    }

    const uint8_t* getVecData(uint32_t vecIx) const
    {
      if (vecIx >= regCount_)
        return nullptr;
      return data_ + vecIx*bytesPerReg_;
    }

    /// It is convenient to contruct an empty regiter file (bytesPerReg = 0)
    /// and configure it later. Old configuration is lost. Register of
    /// newly configured file are initlaized to zero.
    void config(uint32_t bytesPerReg, uint32_t minBytesPerElem,
		uint32_t maxBytesPerElem);

    void reset();

    uint32_t startIndex() const
    { return start_; }

    void setStartIndex(uint32_t start)
    { start_ = start; }

    /// Return currently configure element count (cached valye of VL).
    uint32_t elemCount() const
    { return elems_; }

    /// Set currently configure element count (cached valye of VL).
    void elemCount(uint32_t n)
    { elems_ = n; }

    /// Set the currently configured element width.
    void elemWidth(ElementWidth ew)
    { sew_ = ew; }

    /// Set the currently configured group multiplier.
    void groupMultiplier(GroupMultiplier gm)
    { group_ = gm; groupX8_ = groupMultiplierX8(gm); }

    /// Return true if current vtype configuration is legal. This is a cached
    /// value of VTYPE.VILL.
    bool legalConfig() const
    { return not vill_; }

    /// Return true if the given element width and grouping
    /// combination is legal.
    bool legalConfig(ElementWidth ew, GroupMultiplier mul) const
    {
      if (size_t(ew) >= legalConfigs_.size()) return false;
      const auto& groupFlags = legalConfigs_.at(size_t(ew));
      if (size_t(mul) >= groupFlags.size()) return false;
      return groupFlags.at(size_t(mul));
    }

    /// Update cached vtype fields. This is called when Vsetvli is
    /// executed.
    void updateConfig(ElementWidth sew, GroupMultiplier gm,
                      bool maskAgn, bool tailAgn, bool illegal)
    {
      sew_ = sew;
      group_ = gm;
      maskAgn_ = maskAgn;
      tailAgn_ = tailAgn;
      vill_ = illegal;

      groupX8_ = groupMultiplierX8(gm);
      sewInBits_ = elementWidthInBits(sew);
    }

  private:

    /// Map an vector group multiplier to a flag indicating whether given
    /// group is supported.
    typedef std::vector<bool> GroupFlags;

    /// Map an element width to a vector of flags indicating supported groups.
    typedef std::vector<GroupFlags> GroupsForWidth;

    uint32_t regCount_ = 0;
    uint32_t bytesPerReg_ = 0;
    uint32_t minBytesPerElem_ = 0;
    uint32_t maxBytesPerElem_ = 0;
    uint32_t bytesInRegFile_ = 0;
    uint8_t* data_ = nullptr;

    uint32_t start_ = 0;                           // Cached VSTART
    uint32_t elems_ = 0;                           // Cached VL
    ElementWidth sew_ = ElementWidth::Byte;        // Cached VTYPE.SEW
    GroupMultiplier group_ = GroupMultiplier::One; // Cached VTYPE.VLMUL
    bool maskAgn_ = false;                         // Cached VTYPE.ma
    bool tailAgn_ = false;                         // Cached VTYPE.ta
    bool vill_ = false;                            // Cached VTYPE.VILL

    uint32_t groupX8_ = 8;    // Group multipler as a number scaled by 8.
    uint32_t sewInBits_ = 8;  // SEW expressed in bits (Byte corresponds to 8).

    GroupsForWidth legalConfigs_;

    int lastWrittenReg_ = -1;
    uint32_t lastGroupX8_ = 8;   // 8 times last grouping factor

    // Following used for logging/tracing. Cleared before each instruction.
    // Collected by a vector load/store instruction.
    std::vector<uint64_t> ldStAddr_;  // Addresses of vector load/store instruction
    std::vector<uint64_t> stData_;    // Data of vector store instruction
    std::vector<unsigned> opsEmul_;   // Effecive grouping of vector operands.
  };
}
