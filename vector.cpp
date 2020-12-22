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
#include <cfenv>
#include <cmath>
#include <climits>
#include <cassert>
#include <boost/multiprecision/cpp_int.hpp>
#include "wideint.hpp"
#include "instforms.hpp"
#include "DecodedInst.hpp"
#include "Hart.hpp"


// make_unsigned/make_signed do work on our types -- compensate.
namespace std
{
  template <>
  struct
  make_unsigned<WdRiscv::Int128>
  {
    typedef WdRiscv::Uint128 type;
  };

  template <>
  struct
  make_unsigned<WdRiscv::Int256>
  {
    typedef WdRiscv::Uint256 type;
  };

  template <>
  struct
  make_unsigned<WdRiscv::Int512>
  {
    typedef WdRiscv::Uint512 type;
  };

  template <>
  struct
  make_unsigned<WdRiscv::Int1024>
  {
    typedef WdRiscv::Uint1024 type;
  };

  template <>
  struct
  make_unsigned<WdRiscv::Uint128>
  {
    typedef WdRiscv::Uint128 type;
  };

  template <>
  struct
  make_unsigned<WdRiscv::Uint256>
  {
    typedef WdRiscv::Uint256 type;
  };

  template <>
  struct
  make_unsigned<WdRiscv::Uint512>
  {
    typedef WdRiscv::Uint512 type;
  };

  template <>
  struct
  make_unsigned<WdRiscv::Uint1024>
  {
    typedef WdRiscv::Uint1024 type;
  };

  template <>
  struct
  make_signed<WdRiscv::Uint128>
  {
    typedef WdRiscv::Int128 type;
  };

  template <>
  struct
  make_signed<WdRiscv::Uint256>
  {
    typedef WdRiscv::Int256 type;
  };

  template <>
  struct
  make_signed<WdRiscv::Uint512>
  {
    typedef WdRiscv::Int512 type;
  };

  template <>
  struct
  make_signed<WdRiscv::Uint1024>
  {
    typedef WdRiscv::Int1024 type;
  };

  template <>
  struct
  make_signed<WdRiscv::Int128>
  {
    typedef WdRiscv::Int128 type;
  };

  template <>
  struct
  make_signed<WdRiscv::Int256>
  {
    typedef WdRiscv::Int256 type;
  };

  template <>
  struct
  make_signed<WdRiscv::Int512>
  {
    typedef WdRiscv::Int512 type;
  };

  template <>
  struct
  make_signed<WdRiscv::Int1024>
  {
    typedef WdRiscv::Int1024 type;
  };
}


namespace WdRiscv
{
  /// Return the width in bits of the given integer type T. This is
  /// usually 8*sizeof(T) but is different for wide types implemented
  /// using boost multiprecision types.
  template <typename T>
  unsigned
  integerWidth()
  {
    if constexpr (std::is_same<Int128, T>::value)   return 128;
    if constexpr (std::is_same<Int256, T>::value)   return 256;
    if constexpr (std::is_same<Int512, T>::value)   return 512;
    if constexpr (std::is_same<Int1024, T>::value)  return 1024;
    if constexpr (std::is_same<Uint128, T>::value)  return 128;
    if constexpr (std::is_same<Uint256, T>::value)  return 256;
    if constexpr (std::is_same<Uint512, T>::value)  return 512;
    if constexpr (std::is_same<Uint1024, T>::value) return 1024;

    return 8*sizeof(T);
  }


  /// Return the integral type that is twice as wide as the given
  /// type. For example:
  ///    makeDoubleWide<uint16_t>::type
  /// yields the type
  ///    uint32_t.
  template <typename T>
  struct makeDoubleWide
  {
  };

  template <> struct makeDoubleWide<uint8_t>    { typedef uint16_t type; };
  template <> struct makeDoubleWide<uint16_t>   { typedef uint32_t type; };
  template <> struct makeDoubleWide<uint32_t>   { typedef uint64_t type; };
  template <> struct makeDoubleWide<uint64_t>   { typedef Uint128  type; };
  template <> struct makeDoubleWide<Uint128>    { typedef Uint256  type; };
  template <> struct makeDoubleWide<Uint256>    { typedef Uint512  type; };
  template <> struct makeDoubleWide<Uint512>    { typedef Uint1024 type; };

  template <> struct makeDoubleWide<int8_t>     { typedef int16_t type; };
  template <> struct makeDoubleWide<int16_t>    { typedef int32_t type; };
  template <> struct makeDoubleWide<int32_t>    { typedef int64_t type; };
  template <> struct makeDoubleWide<int64_t>    { typedef Int128  type; };
  template <> struct makeDoubleWide<Int128>     { typedef Int256  type; };
  template <> struct makeDoubleWide<Int256>     { typedef Int512  type; };
  template <> struct makeDoubleWide<Int512>     { typedef Int1024 type; };


  template <typename T>
  T
  minVal()
  {
    typedef typename std::make_unsigned<T>::type UT;
    if constexpr (std::is_same<T, UT>::value)
      return T(0);
    else
      {
        unsigned amount = sizeof(T)*8 - 1;
        return T(1) << amount;
      }
  }


  template <typename T>
  T
  maxVal()
  {
    typedef typename std::make_unsigned<T>::type UT;
    if constexpr (std::is_same<T, UT>::value)
      return ~T(0);
    else
      return (~UT(0) >> 1);
  }


  /// Set result to the upper half of a*b computed in double width
  /// intermediate.
  template <typename T>
  void mulh(const T& a, const T& b, T& result)
  {
    typedef typename makeDoubleWide<T>::type T2; // Double wide type

    unsigned tbits = integerWidth<T> (); // Number of bits in T

    T2 temp = a;
    temp *= b;
    temp >>= tbits;
    result = T(temp);
  }


  /// Specialized mulh for 1024-bit unsigned operands.
  template <>
  void mulh(const Uint1024& a, const Uint1024& b, Uint1024& result)
  {
    // Unpack b into 512-bit pieces
    Uint512 a0 = Uint512(a), a1 = Uint512(a >> 512);
    Uint512 b0 = Uint512(b), b1 = Uint512(b >> 512);

    // Multiply the 4 pieces accumulating results. Maintain upper 1024
    // bits of resuts.
    Uint1024 temp = a0;
    temp *= b0;
    result = temp;

    result >>= 512;
    temp = a1;
    temp *= b0;
    result += temp;

    result >>= 512;
    temp = a0;
    temp *= b1;
    result += temp;
    
    result >>= 512;
    temp = a1;
    temp *= b1;
    result += temp;
  }


  /// Specialized mulh for 1024-bit signed operands.
  template <>
  void mulh(const Int1024& a, const Int1024& b, Int1024& result)
  {
    // Unpack b into 512-bit pieces
    Int512 a0 = Int512(a), a1 = Int512(a >> 512);
    Int512 b0 = Int512(b), b1 = Int512(b >> 512);

    // Multiply the 4 pieces accumulating results. Maintain upper 1024
    // bits of resuts.
    Int1024 temp = a0;
    temp *= b0;
    result = temp;

    result >>= 512;
    temp = a1;
    temp *= b0;
    result += temp;

    result >>= 512;
    temp = a0;
    temp *= b1;
    result += temp;
    
    result >>= 512;
    temp = a1;
    temp *= b1;
    result += temp;
  }


  /// Set result to the upper half of a*b computed in double width
  /// intermediate. TS is a signed integer type (e.g. int8_t).
  /// TU is the corresponding unsigned integer type (e.g. uint8_t).
  template <typename TS, typename TU>
  void mulhsu(const TS& a, const TU& b, TS& result)
  {
    typedef typename makeDoubleWide<TS>::type TS2; // Double wide signed type

    unsigned bits = integerWidth<TS> (); // Number of bits in TS and TU

    TS2 temp = a;
    temp *= TS(b);
    temp >>= bits;
    result = TS(temp);
  }


  /// Specialized mulhsu for 1024-bit unsigned operands.
  template <>
  void mulhsu(const Int1024& a, const Uint1024& b, Int1024& result)
  {
    if (a > Int1024(0))
      {
        Uint1024 ua(a);
        Uint1024 temp = 0;
        mulh(ua, b, temp);
        result = temp;
        return;
      }

    Int1024 smallest = Int1024(1) << 1023;  // Smallest 1024-bit integer.
    if (a == smallest)
      {
        // Result consists of bits 1024 to 2047 of -(power(2,1023)*b)
        // computed in unlimited precision.
        Uint1024 temp = b;
        temp >>= 1;
        result = temp;
        result = -result;
        return;
      }

    Uint1024 nega(-a), temp(0);
    mulh(nega, b, temp);
    result = temp;
    result = -result;
  }

  /// Set result to the product of a and b where a is signed and
  /// is an unsigned and where a and b have the same width.
  template <typename TS, typename TU>
  void mulsu(const TS& a, const TU& b, TS& result)
  {
    TU aa = TU(a);
    aa *= b;
    result = TS(a);
    if (a < 0)
      result = - result;
  }
}


using namespace WdRiscv;


template <typename URV>
bool
Hart<URV>::checkMaskableInst(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return false;
    }

  if (di->isMasked() and di->op0() == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return false;
    }

  return true;
}


template <typename URV>
void
Hart<URV>::vsetvl(unsigned rd, unsigned rs1, URV vtypeVal)
{
  bool ma = (vtypeVal >> 7) & 1;  // Mask agnostic
  bool ta = (vtypeVal >> 6) & 1;  // Tail agnostic
  GroupMultiplier gm = GroupMultiplier(vtypeVal & 7);
  ElementWidth ew = ElementWidth((vtypeVal >> 3) & 7);

  bool vill = (vtypeVal >> (8*sizeof(URV) - 1)) & 1;
  vill = vill or not vecRegs_.legalConfig(ew, gm);

  // Determine vl
  URV elems = 0;

  if (gm == GroupMultiplier::Reserved)
    vill = true;
  else
    {
      uint32_t gm8 = vecRegs_.groupMultiplierX8(gm);
      unsigned bitsPerElem = vecRegs_.elementWidthInBits(ew);
      unsigned vlmax = (gm8*vecRegs_.bitsPerRegister()/bitsPerElem) / 8;
      if (vlmax == 0)
        vill = true;
      else
        {
          if (rd != 0 and rs1 == 0)
            elems = vlmax;
          else if (rd == 0 and rs1 == 0)
            peekCsr(CsrNumber::VL, elems);  // Keep current value of VL.
          else  // strip mining
            {
              URV avl = intRegs_.read(rs1);  // Applical vl
              if (avl <= vlmax)
                elems = avl;
              else if (avl >= 2*vlmax)
                elems = vlmax;
              else
                elems = (avl + 1) / 2;
            }
        }
    }

  if (vill)
    {
      ma = false; ta = false; gm = GroupMultiplier(0); ew = ElementWidth(0);
      elems = 0;
    }

  if (vill or (rd != 0 or rs1 != 0))
    csRegs_.write(CsrNumber::VL, PrivilegeMode::Machine,elems);  // Update VL.
  vecRegs_.elemCount(elems);  // Update cached value of VL.

  intRegs_.write(rd, elems);

  // Pack vtype values and update vtype
  URV vtype = 0;
  vtype |= URV(gm) | (URV(ew) << 3) | (URV(ta) << 6) | (URV(ma) << 6);
  vtype |= (URV(vill) << (8*sizeof(URV) - 1));
  csRegs_.write(CsrNumber::VTYPE, PrivilegeMode::Machine, vtype);

  // Update cached vtype fields in vecRegs_.
  vecRegs_.updateConfig(ew, gm, ma, ta, vill);
}


template <typename URV>
void
Hart<URV>::execVsetvli(const DecodedInst* di)
{
  if (not isVecLegal())
    {
      illegalInst(di);
      return;
    }

  unsigned rd = di->op0();
  unsigned rs1 = di->op1();
  unsigned imm = di->op2();
  
  URV vtypeVal = imm;
  vsetvl(rd, rs1, vtypeVal);
}


template <typename URV>
void
Hart<URV>::execVsetvl(const DecodedInst* di)
{
  if (not isVecLegal())
    {
      illegalInst(di);
      return;
    }
  
  unsigned rd = di->op0();
  unsigned rs1 = di->op1();

  URV vtypeVal = intRegs_.read(di->op2());
  vsetvl(rd, rs1, vtypeVal);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vadd_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 + e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVadd_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vadd_vv<int8_t> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half:   vadd_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word:   vadd_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2:  vadd_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4:  vadd_vv<Int128> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8:  vadd_vv<Int256> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vadd_vv<Int512> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vadd_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vadd_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 + e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVadd_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vadd_vx<int8_t> (vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Half:   vadd_vx<int16_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word:   vadd_vx<int32_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word2:  vadd_vx<int64_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word4:  vadd_vx<Int128> (vd, vs1, Int128(e2),  group, start, elems, masked); break;
    case EW::Word8:  vadd_vx<Int256> (vd, vs1, Int256(e2),  group, start, elems, masked); break;
    case EW::Word16: vadd_vx<Int512> (vd, vs1, Int512(e2),  group, start, elems, masked); break;
    case EW::Word32: vadd_vx<Int1024>(vd, vs1, Int1024(e2), group, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVadd_vi(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(), vs1 = di->op1();
  int32_t imm = di->op2As<int32_t>();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vadd_vx<int8_t> (vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Half:   vadd_vx<int16_t>(vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Word:   vadd_vx<int32_t>(vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Word2:  vadd_vx<int64_t>(vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Word4:  vadd_vx<Int128> (vd, vs1, Int128(imm),  group, start, elems, masked); break;
    case EW::Word8:  vadd_vx<Int256> (vd, vs1, Int256(imm),  group, start, elems, masked); break;
    case EW::Word16: vadd_vx<Int512> (vd, vs1, Int512(imm),  group, start, elems, masked); break;
    case EW::Word32: vadd_vx<Int1024>(vd, vs1, Int1024(imm), group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vsub_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 - e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVsub_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vsub_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vsub_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vsub_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vsub_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vsub_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vsub_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vsub_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vsub_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vsub_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = SRV(intRegs_.read(rs2)), dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 - e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVsub_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vsub_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vsub_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vsub_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vsub_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vsub_vx<Int128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vsub_vx<Int256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vsub_vx<Int512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vsub_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vrsub_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = SRV(intRegs_.read(rs2)), dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e2 - e1;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVrsub_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(), vs1 = di->op1(), rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vrsub_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vrsub_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vrsub_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vrsub_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vrsub_vx<Int128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vrsub_vx<Int256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vrsub_vx<Int512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vrsub_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vrsub_vi(unsigned vd, unsigned vs1, int32_t imm, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = imm, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e2 - e1;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVrsub_vi(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  int32_t imm = di->op2As<int32_t>();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vrsub_vi<int8_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Half: vrsub_vi<int16_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word: vrsub_vi<int32_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word2: vrsub_vi<int64_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word4: vrsub_vi<Int128>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word8: vrsub_vi<Int256>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word16: vrsub_vi<Int512>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word32: vrsub_vi<Int1024>(vd, vs1, imm, group, start, elems, masked); break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vwadd_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type DWT; // Double wide type
  unsigned errors = 0, wideGroup = group*2;

  ELEM_TYPE e1 = 0, e2 = 0;
  DWT dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = DWT(e1);
          dest += DWT(e2);
          if (not vecRegs_.write(vd, ix, wideGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwaddu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwadd_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vwadd_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vwadd_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vwadd_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vwadd_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vwadd_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vwadd_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVwadd_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwadd_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vwadd_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vwadd_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vwadd_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vwadd_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vwadd_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vwadd_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vwadd_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type DWT; // Double wide type
  unsigned errors = 0, doubleGroup = group*2;

  ELEM_TYPE e1 = 0;
  DWT dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = DWT(e1);
          dest += DWT(e2);
          if (not vecRegs_.write(vd, ix, doubleGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwaddu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  URV e2 = intRegs_.read(di->op2());

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwadd_vx<uint8_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Half: vwadd_vx<uint16_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word: vwadd_vx<uint32_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word2: vwadd_vx<uint64_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word4: vwadd_vx<Uint128>(vd, vs1, Uint128(e2), group, start, elems, masked); break;
    case EW::Word8: vwadd_vx<Uint256>(vd, vs1, Uint256(e2), group, start, elems, masked); break;
    case EW::Word16: vwadd_vx<Uint512>(vd, vs1, Uint512(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVwadd_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  SRV e2 = SRV(intRegs_.read(di->op2()));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwadd_vx<int8_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Half: vwadd_vx<int16_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word: vwadd_vx<int32_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word2: vwadd_vx<int64_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word4: vwadd_vx<Int128>(vd, vs1, Int128(e2), group, start, elems, masked); break;
    case EW::Word8: vwadd_vx<Int256>(vd, vs1, Int256(e2), group, start, elems, masked); break;
    case EW::Word16: vwadd_vx<Int512>(vd, vs1, Int512(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVwsubu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  URV e2 = intRegs_.read(di->op2()); // FIX: Spec says sign extened. We differ.

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwadd_vx<uint8_t>(vd, vs1, -uint8_t(e2), group, start, elems, masked); break;
    case EW::Half: vwadd_vx<uint16_t>(vd, vs1, -uint16_t(e2), group, start, elems, masked); break;
    case EW::Word: vwadd_vx<uint32_t>(vd, vs1, -uint32_t(e2), group, start, elems, masked); break;
    case EW::Word2: vwadd_vx<uint64_t>(vd, vs1, -uint64_t(e2), group, start, elems, masked); break;
    case EW::Word4: vwadd_vx<Uint128>(vd, vs1, -Uint128(e2), group, start, elems, masked); break;
    case EW::Word8: vwadd_vx<Uint256>(vd, vs1, Uint256(0)-Uint256(e2), group, start, elems, masked); break;
    case EW::Word16: vwadd_vx<Uint512>(vd, vs1, Uint512(0)-Uint512(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVwsub_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  SRV e2 = SRV(intRegs_.read(di->op2()));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwadd_vx<int8_t>(vd, vs1, -e2, group, start, elems, masked); break;
    case EW::Half: vwadd_vx<int16_t>(vd, vs1, -e2, group, start, elems, masked); break;
    case EW::Word: vwadd_vx<int32_t>(vd, vs1, -e2, group, start, elems, masked); break;
    case EW::Word2: vwadd_vx<int64_t>(vd, vs1, -e2, group, start, elems, masked); break;
    case EW::Word4: vwadd_vx<Int128>(vd, vs1, -Int128(e2), group, start, elems, masked); break;
    case EW::Word8: vwadd_vx<Int256>(vd, vs1, -Int256(e2), group, start, elems, masked); break;
    case EW::Word16: vwadd_vx<Int512>(vd, vs1, -Int512(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template<typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vwsub_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type DWT; // Double wide type
  unsigned errors = 0, doubleGroup = group*2;

  ELEM_TYPE e1 = 0, e2 = 0;
  DWT dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = DWT(e1);
          dest -= DWT(e2);
          if (not vecRegs_.write(vd, ix, doubleGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}

                        
template <typename URV>
void
Hart<URV>::execVwsubu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwsub_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vwsub_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vwsub_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vwsub_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vwsub_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vwsub_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vwsub_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVwsub_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwsub_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vwsub_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vwsub_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vwsub_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vwsub_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vwsub_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vwsub_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vwadd_wv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type DWT; // Double wide type
  unsigned errors = 0, doubleGroup = group*2;

  ELEM_TYPE e2 = 0;
  DWT e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, doubleGroup, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1;
          dest += DWT(e2);
          if (not vecRegs_.write(vd, ix, doubleGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwaddu_wv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwadd_wv<uint8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vwadd_wv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vwadd_wv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vwadd_wv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vwadd_wv<Uint128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vwadd_wv<Uint256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vwadd_wv<Uint512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVwadd_wv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwadd_wv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vwadd_wv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vwadd_wv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vwadd_wv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vwadd_wv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vwadd_wv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vwadd_wv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVwaddu_wx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  URV e2 = intRegs_.read(di->op2());

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vadd_vx<uint16_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Half: vadd_vx<uint32_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word: vadd_vx<uint64_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word2: vadd_vx<Uint128>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word4: vadd_vx<Uint256>(vd, vs1, Uint256(e2), group, start, elems, masked); break;
    case EW::Word8: vadd_vx<Uint512>(vd, vs1, Uint512(e2), group, start, elems, masked); break;
    case EW::Word16: vadd_vx<Uint1024>(vd, vs1, Uint1024(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVwadd_wx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  SRV e2 = SRV(intRegs_.read(di->op2()));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vadd_vx<int16_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Half: vadd_vx<int32_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word: vadd_vx<int64_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word2: vadd_vx<Int128>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word4: vadd_vx<Int256>(vd, vs1, Int256(e2), group, start, elems, masked); break;
    case EW::Word8: vadd_vx<Int512>(vd, vs1, Int512(e2), group, start, elems, masked); break;
    case EW::Word16: vadd_vx<Int1024>(vd, vs1, Int1024(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVwsubu_wx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  URV e2 = intRegs_.read(di->op2()); // FIX: Spec says sign extened. We differ.

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vadd_vx<uint16_t>(vd, vs1, -e2, group, start, elems, masked); break;
    case EW::Half: vadd_vx<uint32_t>(vd, vs1, -e2, group, start, elems, masked); break;
    case EW::Word: vadd_vx<uint64_t>(vd, vs1, -e2, group, start, elems, masked); break;
    case EW::Word2: vadd_vx<Uint128>(vd, vs1, -e2, group, start, elems, masked); break;
    case EW::Word4: vadd_vx<Uint256>(vd, vs1, Uint256(0)-Uint256(e2), group, start, elems, masked); break;
    case EW::Word8: vadd_vx<Uint512>(vd, vs1, Uint512(0)-Uint512(e2), group, start, elems, masked); break;
    case EW::Word16: vadd_vx<Uint1024>(vd, vs1, Uint1024(0)-Uint1024(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVwsub_wx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  SRV e2 = SRV(intRegs_.read(di->op2()));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vadd_vx<int16_t>(vd, vs1, -e2, group, start, elems, masked); break;
    case EW::Half: vadd_vx<int32_t>(vd, vs1, -e2, group, start, elems, masked); break;
    case EW::Word: vadd_vx<int64_t>(vd, vs1, -e2, group, start, elems, masked); break;
    case EW::Word2: vadd_vx<Int128>(vd, vs1, -e2, group, start, elems, masked); break;
    case EW::Word4: vadd_vx<Int256>(vd, vs1, -Int256(e2), group, start, elems, masked); break;
    case EW::Word8: vadd_vx<Int512>(vd, vs1, -Int512(e2), group, start, elems, masked); break;
    case EW::Word16: vadd_vx<Int1024>(vd, vs1, -Int1024(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vwsub_wv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type DWT; // Double wide type
  unsigned errors = 0, doubleGroup = group*2;

  ELEM_TYPE e2 = 0;
  DWT e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, doubleGroup, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1;
          dest -= DWT(e2);
          if (not vecRegs_.write(vd, ix, doubleGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwsubu_wv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwsub_wv<uint8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vwsub_wv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vwsub_wv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vwsub_wv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vwsub_wv<Uint128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vwsub_wv<Uint256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vwsub_wv<Uint512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVwsub_wv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  if (not vecRegs_.isDoubleWideLegal(sew, group))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwsub_wv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vwsub_wv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vwsub_wv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vwsub_wv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vwsub_wv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vwsub_wv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vwsub_wv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vminu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 < e2 ? e1 : e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVminu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vminu_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vminu_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vminu_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vminu_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vminu_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vminu_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vminu_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vminu_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vminu_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  // Spec (sep 24, 2020) says this should be sign extended. We hope
  // they come to their senses.
  ELEM_TYPE e2 = intRegs_.read(rs2);
  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 < e2 ? e1 : e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVminu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vminu_vx<uint8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vminu_vx<uint16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vminu_vx<uint32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vminu_vx<uint64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vminu_vx<Uint128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vminu_vx<Uint256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vminu_vx<Uint512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vminu_vx<Uint1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmin_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 < e2 ? e1 : e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmin_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmin_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vmin_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vmin_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vmin_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vmin_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vmin_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vmin_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vmin_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmin_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = SRV(intRegs_.read(rs2)), dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 < e2 ? e1 : e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmin_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(), rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmin_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vmin_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vmin_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vmin_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vmin_vx<Int128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vmin_vx<Int256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vmin_vx<Int512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vmin_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmaxu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 > e2 ? e1 : e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmaxu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmaxu_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vmaxu_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vmaxu_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vmaxu_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vmaxu_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vmaxu_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vmaxu_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vmaxu_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmaxu_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

   // Spec (sep 24, 2020) says this should be sign extended. We hope
   // they come to their senses.
  ELEM_TYPE e2 = intRegs_.read(rs2);
  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 > e2 ? e1 : e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmaxu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmaxu_vx<uint8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vmaxu_vx<uint16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vmaxu_vx<uint32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vmaxu_vx<uint64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vmaxu_vx<Uint128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vmaxu_vx<Uint256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vmaxu_vx<Uint512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vmaxu_vx<Uint1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmax_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 > e2 ? e1 : e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmax_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmax_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vmax_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vmax_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vmax_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vmax_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vmax_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vmax_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vmax_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmax_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = SRV(intRegs_.read(rs2)), dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 > e2 ? e1 : e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmax_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmax_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vmax_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vmax_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vmax_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vmax_vx<Int128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vmax_vx<Int256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vmax_vx<Int512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vmax_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vand_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 & e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVand_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vand_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vand_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vand_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vand_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vand_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vand_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vand_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vand_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vand_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  // Spec says sign extend scalar register. We comply. Looks foolish.
  ELEM_TYPE e1 = 0, e2 = SRV(intRegs_.read(rs2)), dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 & e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVand_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vand_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vand_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vand_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vand_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vand_vx<Int128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vand_vx<Int256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vand_vx<Int512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vand_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vand_vi(unsigned vd, unsigned vs1, int32_t imm, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  // Spec says sign extend immediate. We comply. Looks foolish.
  ELEM_TYPE e1 = 0, e2 = imm, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 & e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVand_vi(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(), vs1 = di->op1();
  int32_t imm = di->op2As<int32_t>();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vand_vi<int8_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Half: vand_vi<int16_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word: vand_vi<int32_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word2: vand_vi<int64_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word4: vand_vi<Int128>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word8: vand_vi<Int256>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word16: vand_vi<Int512>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word32: vand_vi<Int1024>(vd, vs1, imm, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vor_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                  unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 | e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVor_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vor_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vor_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vor_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vor_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vor_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vor_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vor_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vor_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vor_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                  unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  // Spec says sign extend scalar register. We comply. Looks foolish.
  ELEM_TYPE e1 = 0, e2 = SRV(intRegs_.read(rs2)), dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 | e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVor_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(), rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vor_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vor_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vor_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vor_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vor_vx<Int128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vor_vx<Int256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vor_vx<Int512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vor_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vor_vi(unsigned vd, unsigned vs1, int32_t imm, unsigned group,
                  unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  ELEM_TYPE val2 = imm;
  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 | val2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVor_vi(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  int32_t imm = di->op2As<int32_t>();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vor_vi<int8_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Half: vor_vi<int16_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word: vor_vi<int32_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word2: vor_vi<int64_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word4: vor_vi<Int128>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word8: vor_vi<Int256>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word16: vor_vi<Int512>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word32: vor_vi<Int1024>(vd, vs1, imm, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vxor_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 ^ e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVxor_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vxor_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vxor_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vxor_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vxor_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vxor_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vxor_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vxor_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vxor_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vxor_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  // Spec says sign extend scalar register. We comply. Looks foolish.
  ELEM_TYPE e1 = 0, e2 = SRV(intRegs_.read(rs2)), dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 ^ e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVxor_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(), rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vxor_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vxor_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vxor_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vxor_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vxor_vx<Int128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vxor_vx<Int256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vxor_vx<Int512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vxor_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vxor_vi(unsigned vd, unsigned vs1, int32_t imm, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  ELEM_TYPE val2 = imm;
  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 ^ val2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVxor_vi(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  int32_t imm = di->op2As<int32_t>();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vxor_vi<int8_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Half: vxor_vi<int16_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word: vxor_vi<int32_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word2: vxor_vi<int64_t>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word4: vxor_vi<Int128>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word8: vxor_vi<Int256>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word16: vxor_vi<Int512>(vd, vs1, imm, group, start, elems, masked); break;
    case EW::Word32: vxor_vi<Int1024>(vd, vs1, imm, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vrgather_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                       unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.read(vs2, ix, group, e2))
        {
          unsigned vs1Ix = unsigned(e2);
          dest = 0;
          if (vecRegs_.read(vs1, vs1Ix, group, e1))
            dest = e1;

          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVrgather_vv(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vrgather_vv<uint8_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Half: vrgather_vv<uint16_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word: vrgather_vv<uint32_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word2: vrgather_vv<uint64_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word4: vrgather_vv<Uint128>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word8: vrgather_vv<Uint256>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word16: vrgather_vv<Uint512>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word32: vrgather_vv<Uint1024>(vd, vs1, vs2, group, start, elems); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vrgather_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                       unsigned start, unsigned elems)
{
  unsigned errors = 0;

  unsigned bytesPerElem = sizeof(ELEM_TYPE);
  unsigned vlmax = group*vecRegs_.bitsPerRegister()/bytesPerElem;

  URV rv2 = intRegs_.read(rs2);
  URV vs1Ix = rv2 < vlmax ? rv2 : vlmax;

  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      dest = 0;
      if (vecRegs_.read(vs1, vs1Ix, group, e1))
        dest = e1;

      if (not vecRegs_.write(vd, ix, group, dest))
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVrgather_vx(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vrgather_vx<uint8_t>(vd, vs1, rs2, group, start, elems); break;
    case EW::Half: vrgather_vx<uint16_t>(vd, vs1, rs2, group, start, elems); break;
    case EW::Word: vrgather_vx<uint32_t>(vd, vs1, rs2, group, start, elems); break;
    case EW::Word2: vrgather_vx<uint64_t>(vd, vs1, rs2, group, start, elems); break;
    case EW::Word4: vrgather_vx<Uint128>(vd, vs1, rs2, group, start, elems); break;
    case EW::Word8: vrgather_vx<Uint256>(vd, vs1, rs2, group, start, elems); break;
    case EW::Word16: vrgather_vx<Uint512>(vd, vs1, rs2, group, start, elems); break;
    case EW::Word32: vrgather_vx<Uint1024>(vd, vs1, rs2, group, start, elems); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vrgather_vi(unsigned vd, unsigned vs1, uint32_t imm, unsigned group,
                       unsigned start, unsigned elems)
{
  uint32_t vs1Ix = imm;
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      dest = 0;
      if (vecRegs_.read(vs1, vs1Ix, group, e1))
        dest = e1;

      if (not vecRegs_.write(vd, ix, group, dest))
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVrgather_vi(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1();
  uint32_t imm = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vrgather_vi<uint8_t>(vd, vs1, imm, group, start, elems); break;
    case EW::Half: vrgather_vi<uint16_t>(vd, vs1, imm, group, start, elems); break;
    case EW::Word: vrgather_vi<uint32_t>(vd, vs1, imm, group, start, elems); break;
    case EW::Word2: vrgather_vi<uint64_t>(vd, vs1, imm, group, start, elems); break;
    case EW::Word4: vrgather_vi<Uint128>(vd, vs1, imm, group, start, elems); break;
    case EW::Word8: vrgather_vi<Uint256>(vd, vs1, imm, group, start, elems); break;
    case EW::Word16: vrgather_vi<Uint512>(vd, vs1, imm, group, start, elems); break;
    case EW::Word32: vrgather_vi<Uint1024>(vd, vs1, imm, group, start, elems); break;
    }
}




template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vrgatherei16_vv(unsigned vd, unsigned vs1, unsigned vs2,
                           unsigned group, unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;
  uint16_t e2 = 0;
  unsigned e2Group = (16*group)/(8*sizeof(ELEM_TYPE));

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.read(vs2, ix, e2Group, e2))
        {
          unsigned vs1Ix = e2;
          dest = 0;
          if (vecRegs_.read(vs1, vs1Ix, group, e1))
            dest = e1;

          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVrgatherei16_vv(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vrgatherei16_vv<uint8_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Half: vrgatherei16_vv<uint16_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word: vrgatherei16_vv<uint32_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word2: vrgatherei16_vv<uint64_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word4: vrgatherei16_vv<Uint128>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word8: vrgatherei16_vv<Uint256>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word16: vrgatherei16_vv<Uint512>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word32: vrgatherei16_vv<Uint1024>(vd, vs1, vs2, group, start, elems); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vcompress_vm(unsigned vd, unsigned vs1, unsigned vs2,
                        unsigned group, unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;
  unsigned destIx = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.isActive(vs2, ix))
        {
          if (vecRegs_.read(vs1, ix, group, e1))
            {
              dest = e1;
              if (not vecRegs_.write(vd, destIx++, group, dest))
                errors++;
            }
          else
            errors++;
        }
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVcompress_vm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig() or di->isMasked())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vcompress_vm<uint8_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Half: vcompress_vm<uint16_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word: vcompress_vm<uint32_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word2: vcompress_vm<uint64_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word4: vcompress_vm<Uint128>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word8: vcompress_vm<Uint256>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word16: vcompress_vm<Uint512>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word32: vcompress_vm<Uint1024>(vd, vs1, vs2, group, start, elems); break;
    }
}



template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vredsum_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                      unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e2 = 0;
  unsigned scalarElemIx = 0, scalarElemGroupX8 = 8;

  if (not vecRegs_.read(vs2, scalarElemIx, scalarElemGroupX8, e2))
    errors++;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result += e1;
    else
      errors++;

  if (not vecRegs_.write(vd, scalarElemIx, scalarElemGroupX8, result))
    errors++;

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVredsum_vs(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vredsum_vs<int8_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Half: vredsum_vs<int16_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word: vredsum_vs<int32_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word2: vredsum_vs<int64_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word4: vredsum_vs<Int128>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word8: vredsum_vs<Int256>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word16: vredsum_vs<Int512>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word32: vredsum_vs<Int1024>(vd, vs1, vs2, group, start, elems); break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vredand_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                      unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e2 = 0;
  unsigned scalarElemIx = 0, scalarElemGroupX8 = 8;

  if (not vecRegs_.read(vs2, scalarElemIx, scalarElemGroupX8, e2))
    errors++;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result = result & e1;
    else
      errors++;

  if (not vecRegs_.write(vd, scalarElemIx, scalarElemGroupX8, result))
    errors++;

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVredand_vs(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vredand_vs<uint8_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Half: vredand_vs<uint16_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word: vredand_vs<uint32_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word2: vredand_vs<uint64_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word4: vredand_vs<Uint128>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word8: vredand_vs<Uint256>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word16: vredand_vs<Uint512>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word32: vredand_vs<Uint1024>(vd, vs1, vs2, group, start, elems); break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vredor_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                     unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e2 = 0;
  unsigned scalarElemIx = 0, scalarElemGroupX8 = 8;

  if (not vecRegs_.read(vs2, scalarElemIx, scalarElemGroupX8, e2))
    errors++;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result = result | e1;
    else
      errors++;

  if (not vecRegs_.write(vd, scalarElemIx, scalarElemGroupX8, result))
    errors++;

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVredor_vs(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vredor_vs<uint8_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Half: vredor_vs<uint16_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word: vredor_vs<uint32_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word2: vredor_vs<uint64_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word4: vredor_vs<Uint128>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word8: vredor_vs<Uint256>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word16: vredor_vs<Uint512>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word32: vredor_vs<Uint1024>(vd, vs1, vs2, group, start, elems); break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vredxor_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                      unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e2 = 0;
  unsigned scalarElemIx = 0, scalarElemGroupX8 = 8;

  if (not vecRegs_.read(vs2, scalarElemIx, scalarElemGroupX8, e2))
    errors++;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result = result ^ e1;
    else
      errors++;

  if (not vecRegs_.write(vd, scalarElemIx, scalarElemGroupX8, result))
    errors++;

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVredxor_vs(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vredxor_vs<uint8_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Half: vredxor_vs<uint16_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word: vredxor_vs<uint32_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word2: vredxor_vs<uint64_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word4: vredxor_vs<Uint128>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word8: vredxor_vs<Uint256>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word16: vredxor_vs<Uint512>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word32: vredxor_vs<Uint1024>(vd, vs1, vs2, group, start, elems); break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vredminu_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                       unsigned start, unsigned elems)
{
  ELEM_TYPE e2 = 0;
  unsigned errors = 0, scalarElemIx = 0, scalarElemGroupX8 = 8;

  if (not vecRegs_.read(vs2, scalarElemIx, scalarElemGroupX8, e2))
    errors++;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result = result < e1 ? result : e1;
    else
      errors++;

  if (not vecRegs_.write(vd, scalarElemIx, scalarElemGroupX8, result))
    errors++;

  assert(errors == 0);

}


template <typename URV>
void
Hart<URV>::execVredminu_vs(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vredminu_vs<uint8_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Half: vredminu_vs<uint16_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word: vredminu_vs<uint32_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word2: vredminu_vs<uint64_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word4: vredminu_vs<Uint128>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word8: vredminu_vs<Uint256>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word16: vredminu_vs<Uint512>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word32: vredminu_vs<Uint1024>(vd, vs1, vs2, group, start, elems); break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vredmin_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                      unsigned start, unsigned elems)
{
  ELEM_TYPE e2 = 0;
  unsigned errors = 0, scalarElemIx = 0, scalarElemGroupX8 = 8;

  if (not vecRegs_.read(vs2, scalarElemIx, scalarElemGroupX8, e2))
    errors++;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result = result < e1 ? result : e1;
    else
      errors++;

  if (not vecRegs_.write(vd, scalarElemIx, scalarElemGroupX8, result))
    errors++;

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVredmin_vs(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vredmin_vs<int8_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Half: vredmin_vs<int16_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word: vredmin_vs<int32_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word2: vredmin_vs<int64_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word4: vredmin_vs<Int128>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word8: vredmin_vs<Int256>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word16: vredmin_vs<Int512>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word32: vredmin_vs<Int1024>(vd, vs1, vs2, group, start, elems); break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vredmaxu_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                       unsigned start, unsigned elems)
{
  ELEM_TYPE e2 = 0;
  unsigned errors = 0, scalarElemIx = 0, scalarElemGroupX8 = 8;

  if (not vecRegs_.read(vs2, scalarElemIx, scalarElemGroupX8, e2))
    errors++;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result = result > e1 ? result : e1;
    else
      errors++;

  if (not vecRegs_.write(vd, scalarElemIx, scalarElemGroupX8, result))
    errors++;

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVredmaxu_vs(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vredmaxu_vs<uint8_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Half: vredmaxu_vs<uint16_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word: vredmaxu_vs<uint32_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word2: vredmaxu_vs<uint64_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word4: vredmaxu_vs<Uint128>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word8: vredmaxu_vs<Uint256>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word16: vredmaxu_vs<Uint512>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word32: vredmaxu_vs<Uint1024>(vd, vs1, vs2, group, start, elems); break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vredmax_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                      unsigned start, unsigned elems)
{
  ELEM_TYPE e2 = 0;
  unsigned errors = 0, scalarElemIx = 0, scalarElemGroupX8 = 8;

  if (not vecRegs_.read(vs2, scalarElemIx, scalarElemGroupX8, e2))
    errors++;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result = result > e1 ? result : e1;
    else
      errors++;

  if (not vecRegs_.write(vd, scalarElemIx, scalarElemGroupX8, result))
    errors++;

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVredmax_vs(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(), start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vredmax_vs<int8_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Half: vredmax_vs<int16_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word: vredmax_vs<int32_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word2: vredmax_vs<int64_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word4: vredmax_vs<Int128>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word8: vredmax_vs<Int256>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word16: vredmax_vs<Int512>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word32: vredmax_vs<Int1024>(vd, vs1, vs2, group, start, elems); break;
    }
}


template <typename URV>
void
Hart<URV>::execVmand_mm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  if (elems > vecRegs_.bytesPerRegister())
    assert(0);

  uint8_t* vdData = vecRegs_.getVecData(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecData(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecData(di->op2());
  if (not vs1Data or not vs2Data or not vdData)
    assert(0);

  // If bits indices are byte aligned process bytes
  if ((start & 7) == 0 and (elems & 7) == 0)
    {
      start = start >> 3;
      elems = elems >> 3;
      for (unsigned i = start; i < elems; ++i)
        vdData[i] = vs1Data[i] & vs2Data[i];
    }
  else     // Bit indices are not byte aligned.
    for (unsigned i = start; i < elems; ++i)
      {
        unsigned byteIx = i >> 3;
        unsigned bitIx = i & 7; // Bit index in byte
        uint8_t mask = 1 << bitIx;
        vdData[byteIx] = ( (vdData[byteIx] & ~mask) |
                           (vs1Data[byteIx] & vs2Data[byteIx] & mask) );
      }

  vecRegs_.setLastWrittenReg(di->op0(), elems-1, 1);
}


template <typename URV>
void
Hart<URV>::execVmnand_mm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  if (elems > vecRegs_.bytesPerRegister())
    assert(0);

  uint8_t* vdData = vecRegs_.getVecData(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecData(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecData(di->op2());
  if (not vs1Data or not vs2Data or not vdData)
    assert(0);

  // If bits indices are byte aligned process bytes
  if ((start & 7) == 0 and (elems & 7) == 0)
    {
      start = start >> 3;
      elems = elems >> 3;
      for (unsigned i = start; i < elems; ++i)
        vdData[i] = ~ (vs1Data[i] & vs2Data[i]);
    }
  else    // Bit indices are not byte aligned.
    for (unsigned i = start; i < elems; ++i)
      {
        unsigned byteIx = i >> 3;
        unsigned bitIx = i & 7; // Bit index in byte
        uint8_t mask = 1 << bitIx;
        vdData[byteIx] = ( (vdData[byteIx] & ~mask) |
                           (~(vs1Data[byteIx] & vs2Data[byteIx]) & mask) );
      }

  vecRegs_.setLastWrittenReg(di->op0(), elems-1, 1);
}


template <typename URV>
void
Hart<URV>::execVmandnot_mm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  if (elems > vecRegs_.bytesPerRegister())
    assert(0);

  uint8_t* vdData = vecRegs_.getVecData(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecData(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecData(di->op2());
  if (not vs1Data or not vs2Data or not vdData)
    assert(0);

  // If bits indices are byte aligned process bytes
  if ((start & 7) == 0 and (elems & 7) == 0)
    {
      start = start >> 3;
      elems = elems >> 3;
      for (unsigned i = start; i < elems; ++i)
        vdData[i] = vs1Data[i] & ~vs2Data[i];
    }
  else    // Bit indices are not byte aligned.
    for (unsigned i = start; i < elems; ++i)
      {
        unsigned byteIx = i >> 3;
        unsigned bitIx = i & 7; // Bit index in byte
        uint8_t mask = 1 << bitIx;
        vdData[byteIx] = ( (vdData[byteIx] & ~mask) |
                           ((vs1Data[byteIx] & ~vs2Data[byteIx]) & mask) );
      }

  vecRegs_.setLastWrittenReg(di->op0(), elems-1, 1);
}


template <typename URV>
void
Hart<URV>::execVmxor_mm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  if (elems > vecRegs_.bytesPerRegister())
    assert(0);

  uint8_t* vdData = vecRegs_.getVecData(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecData(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecData(di->op2());
  if (not vs1Data or not vs2Data or not vdData)
    assert(0);

  // If bits indices are byte aligned process bytes
  if ((start & 7) == 0 and (elems & 7) == 0)
    {
      start = start >> 3;
      elems = elems >> 3;
      for (unsigned i = start; i < elems; ++i)
        vdData[i] = vs1Data[i] ^ vs2Data[i];
      return;
    }
  else    // Bit indices are not byte aligned.
    for (unsigned i = start; i < elems; ++i)
      {
        unsigned byteIx = i >> 3;
        unsigned bitIx = i & 7; // Bit index in byte
        uint8_t mask = 1 << bitIx;
        vdData[byteIx] = ( (vdData[byteIx] & ~mask) |
                           ((vs1Data[byteIx] ^ vs2Data[byteIx]) & mask) );
      }

  vecRegs_.setLastWrittenReg(di->op0(), elems-1, 1);
}


template <typename URV>
void
Hart<URV>::execVmor_mm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  if (elems > vecRegs_.bytesPerRegister())
    assert(0);

  uint8_t* vdData = vecRegs_.getVecData(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecData(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecData(di->op2());
  if (not vs1Data or not vs2Data or not vdData)
    assert(0);

  // If bits indices are byte aligned process bytes
  if ((start & 7) == 0 and (elems & 7) == 0)
    {
      start = start >> 3;
      elems = elems >> 3;
      for (unsigned i = start; i < elems; ++i)
        vdData[i] = vs1Data[i] | vs2Data[i];
      return;
    }
  else    // Bit indices are not byte aligned.
    for (unsigned i = start; i < elems; ++i)
      {
        unsigned byteIx = i >> 3;
        unsigned bitIx = i & 7; // Bit index in byte
        uint8_t mask = 1 << bitIx;
        vdData[byteIx] = ( (vdData[byteIx] & ~mask) |
                           ((vs1Data[byteIx] | vs2Data[byteIx]) & mask) );
      }

  vecRegs_.setLastWrittenReg(di->op0(), elems-1, 1);
}


template <typename URV>
void
Hart<URV>::execVmnor_mm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  if (elems > vecRegs_.bytesPerRegister())
    assert(0);

  uint8_t* vdData = vecRegs_.getVecData(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecData(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecData(di->op2());
  if (not vs1Data or not vs2Data or not vdData)
    assert(0);

  // If bits indices are byte aligned process bytes
  if ((start & 7) == 0 and (elems & 7) == 0)
    {
      start = start >> 3;
      elems = elems >> 3;
      for (unsigned i = start; i < elems; ++i)
        vdData[i] = ~(vs1Data[i] | vs2Data[i]);
    }
  else    // Bit indices are not byte aligned.
    for (unsigned i = start; i < elems; ++i)
      {
        unsigned byteIx = i >> 3;
        unsigned bitIx = i & 7; // Bit index in byte
        uint8_t mask = 1 << bitIx;
        vdData[byteIx] = ( (vdData[byteIx] & ~mask) |
                           ( ~(vs1Data[byteIx] | vs2Data[byteIx]) & mask) );
      }

  vecRegs_.setLastWrittenReg(di->op0(), elems-1, 1);
}


template <typename URV>
void
Hart<URV>::execVmornot_mm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  if (elems > vecRegs_.bytesPerRegister())
    assert(0);

  uint8_t* vdData = vecRegs_.getVecData(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecData(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecData(di->op2());
  if (not vs1Data or not vs2Data or not vdData)
    assert(0);

  // If bits indices are byte aligned process bytes
  if ((start & 7) == 0 and (elems & 7) == 0)
    {
      start = start >> 3;
      elems = elems >> 3;
      for (unsigned i = start; i < elems; ++i)
        vdData[i] = vs1Data[i] | ~vs2Data[i];
    }
  else     // Bit indices are not byte aligned.
    for (unsigned i = start; i < elems; ++i)
      {
        unsigned byteIx = i >> 3;
        unsigned bitIx = i & 7; // Bit index in byte
        uint8_t mask = 1 << bitIx;
        vdData[byteIx] = ( (vdData[byteIx] & ~mask) |
                           ((vs1Data[byteIx] | ~vs2Data[byteIx]) & mask) );
      }

  vecRegs_.setLastWrittenReg(di->op0(), elems-1, 1);
}


template <typename URV>
void
Hart<URV>::execVmxnor_mm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  unsigned start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  if (elems > vecRegs_.bytesPerRegister())
    assert(0);

  uint8_t* vdData = vecRegs_.getVecData(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecData(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecData(di->op2());
  if (not vs1Data or not vs2Data or not vdData)
    assert(0);

  // If bits indices are byte aligned process bytes
  if ((start & 7) == 0 and (elems & 7) == 0)
    {
      start = start >> 3;
      elems = elems >> 3;
      for (unsigned i = start; i < elems; ++i)
        vdData[i] = vs1Data[i] ^ ~vs2Data[i];
    }
  else     // Bit indices are not byte aligned.
    for (unsigned i = start; i < elems; ++i)
      {
        unsigned byteIx = i >> 3;
        unsigned bitIx = i & 7; // Bit index in byte
        uint8_t mask = 1 << bitIx;
        vdData[byteIx] = ( (vdData[byteIx] & ~mask) |
                           ((vs1Data[byteIx] ^ ~vs2Data[byteIx]) & mask) );
      }

  vecRegs_.setLastWrittenReg(di->op0(), elems-1, 1);
}


template <typename URV>
void
Hart<URV>::execVpopc_m(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  uint32_t start = vecRegs_.startIndex();
  if (start > 0)
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned rd = di->op0(),  vs1 = di->op1(),  elems = vecRegs_.elemCount();;

  uint32_t count = 0;
  for (uint32_t ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;
      if (vecRegs_.isActive(vs1, ix))
        count++;
    }

  intRegs_.write(rd, count);
}


template <typename URV>
void
Hart<URV>::execVfirst_m(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  uint32_t start = vecRegs_.startIndex();
  if (start > 0)
    {
      illegalInst(di);
      return;
    }

  unsigned rd = di->op0(),  vs1 = di->op1(),  elems = vecRegs_.elemCount();

  SRV first = -1;

  for (uint32_t ix = start; ix < elems; ++ix)
    if (vecRegs_.isActive(vs1, ix))
      {
        first = ix;
        break;
      }

  intRegs_.write(rd, first);
}


template <typename URV>
void
Hart<URV>::execVmsbf_m(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  uint32_t start = vecRegs_.startIndex();
  if (start > 0)
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  elems = vecRegs_.elemCount();;

  if (vd == vs1 or (masked and vd == 0))
    {
      illegalInst(di);
      return;
    }

  uint8_t* vdData = vecRegs_.getVecData(vd);
  uint8_t* vs1Data = vecRegs_.getVecData(vs1);

  bool found = false;  // true if set bit is found in vs1

  for (uint32_t ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;
      if (vecRegs_.isActive(vs1, ix))
        {
          unsigned byteIx = ix >> 3;
          unsigned bitIx = ix & 7; // Bit index in byte
          uint8_t mask = 1 << bitIx;
          vdData[ix] = vdData[ix] & ~mask;
          found = found or (vs1Data[byteIx] & mask);
          if (not found)
            vdData[byteIx] |= mask;
        }
    }

  vecRegs_.setLastWrittenReg(vd, elems-1, 1);
}


template <typename URV>
void
Hart<URV>::execVmsif_m(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  uint32_t start = vecRegs_.startIndex();
  if (start > 0)
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  elems = vecRegs_.elemCount();;

  if (vd == vs1 or (masked and vd == 0))
    {
      illegalInst(di);
      return;
    }

  uint8_t* vdData = vecRegs_.getVecData(vd);
  uint8_t* vs1Data = vecRegs_.getVecData(vs1);

  bool found = false;  // true if set bit is found in vs1

  for (uint32_t ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;
      if (vecRegs_.isActive(vs1, ix))
        {
          unsigned byteIx = ix >> 3;
          unsigned bitIx = ix & 7; // Bit index in byte
          uint8_t mask = 1 << bitIx;
          vdData[ix] = vdData[ix] & ~mask;
          if (not found)
            vdData[byteIx] |= mask;
          found = vs1Data[byteIx] & mask;
        }
    }

  vecRegs_.setLastWrittenReg(vd, elems-1, 1);
}


template <typename URV>
void
Hart<URV>::execVmsof_m(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  uint32_t start = vecRegs_.startIndex();
  if (start > 0)
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  elems = vecRegs_.elemCount();;

  if (vd == vs1 or (masked and vd == 0))
    {
      illegalInst(di);
      return;
    }

  uint8_t* vdData = vecRegs_.getVecData(vd);
  uint8_t* vs1Data = vecRegs_.getVecData(vs1);

  bool found = false;  // true if set bit is found in vs1

  for (uint32_t ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;
      if (vecRegs_.isActive(vs1, ix))
        {
          unsigned byteIx = ix >> 3;
          unsigned bitIx = ix & 7; // Bit index in byte
          uint8_t mask = 1 << bitIx;
          vdData[ix] = vdData[ix] & ~mask;
          if (not found)
            {
              found = vs1Data[byteIx] & mask;
              if (found)
                vdData[byteIx] |= mask;
            }
        }
    }

  vecRegs_.setLastWrittenReg(vd, elems-1, 1);
}


template <typename URV>
void
Hart<URV>::execViota_m(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  // Spec does not explicitly state this.   FIX double check.
  uint32_t start = vecRegs_.startIndex();
  if (start > 0)
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  elems = vecRegs_.elemCount();;

  if (masked and vd == 0)
    {
      illegalInst(di);
      return;
    }

  unsigned sum = 0;

  for (uint32_t ix = start; ix < elems; ++ix)
    {
      bool sourceSet = vecRegs_.isActive(vs1, ix);

      if (masked and not vecRegs_.isActive(0, ix))
        {
          if (sourceSet)
            sum++;
          continue;
        }

      switch (sew)
        {
        case ElementWidth::Byte: vecRegs_.write(vd, ix, group, int8_t(sum)); break;
        case ElementWidth::Half: vecRegs_.write(vd, ix, group, int16_t(sum)); break;
        case ElementWidth::Word: vecRegs_.write(vd, ix, group, int32_t(sum)); break;
        case ElementWidth::Word2: vecRegs_.write(vd, ix, group, int64_t(sum)); break;
        case ElementWidth::Word4: vecRegs_.write(vd, ix, group, Int128(sum)); break;
        case ElementWidth::Word8: vecRegs_.write(vd, ix, group, Int256(sum)); break;
        case ElementWidth::Word16: vecRegs_.write(vd, ix, group, Int512(sum)); break;
        case ElementWidth::Word32: vecRegs_.write(vd, ix, group, Int1024(sum)); break;
        }

      if (sourceSet)
        sum++;
    }
}


template <typename URV>
void
Hart<URV>::execVid_v(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  // Spec does not mention vstart > 0. Got a clarification saying it is ok not
  // to ake an exception in that case.
  uint32_t start = vecRegs_.startIndex();

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  elems = vecRegs_.elemCount();;

  if (masked and vd == 0)
    {
      illegalInst(di);
      return;
    }

  for (uint32_t ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      switch (sew)
        {
        case ElementWidth::Byte: vecRegs_.write(vd, ix, group, int8_t(ix)); break;
        case ElementWidth::Half: vecRegs_.write(vd, ix, group, int16_t(ix)); break;
        case ElementWidth::Word: vecRegs_.write(vd, ix, group, int32_t(ix)); break;
        case ElementWidth::Word2: vecRegs_.write(vd, ix, group, int64_t(ix)); break;
        case ElementWidth::Word4: vecRegs_.write(vd, ix, group, Int128(ix)); break;
        case ElementWidth::Word8: vecRegs_.write(vd, ix, group, Int256(ix)); break;
        case ElementWidth::Word16: vecRegs_.write(vd, ix, group, Int512(ix)); break;
        case ElementWidth::Word32: vecRegs_.write(vd, ix, group, Int1024(ix)); break;
        }
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vslideup(unsigned vd, unsigned vs1, URV amount, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (ix < amount)
        continue;

      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      unsigned from = ix - amount;

      if (vecRegs_.read(vs1, from, group, e1))
        {
          dest = e1;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVslideup_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(), rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned dist = vd > vs1 ? vd - vs1 : vs1 - vd;
  if (dist > group)
    {
      illegalInst(di);  // Source/dest vecs cannot overlap
      return;
    }

  URV amount = intRegs_.read(rs2);

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vslideup<uint8_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Half: vslideup<uint16_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word: vslideup<uint32_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word2: vslideup<uint64_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word4: vslideup<Uint128>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word8: vslideup<Uint256>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word16: vslideup<Uint512>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word32: vslideup<Uint1024>(vd, vs1, amount, group, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVslideup_vi(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(), imm = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned dist = vd > vs1 ? vd - vs1 : vs1 - vd;
  if (dist > group)
    {
      illegalInst(di);  // Source/dest vecs cannot overlap
      return;
    }

  URV amount = imm;

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vslideup<uint8_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Half: vslideup<uint16_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word: vslideup<uint32_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word2: vslideup<uint64_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word4: vslideup<Uint128>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word8: vslideup<Uint256>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word16: vslideup<Uint512>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word32: vslideup<Uint1024>(vd, vs1, amount, group, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVslide1up_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(), rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned dist = vd > vs1 ? vd - vs1 : vs1 - vd;
  if (dist > group)
    {
      illegalInst(di);  // Source/dest vecs cannot overlap
      return;
    }

  URV amount = 1;

  // Sign extend scalar value
  SRV replacement = SRV(intRegs_.read(rs2));

  switch (sew)
    {
    case ElementWidth::Byte:
      vslideup<uint8_t>(vd, vs1, amount, group, start, elems, masked);
      vecRegs_.write(vd, 0, group, int8_t(replacement));
      break;

    case ElementWidth::Half:
      vslideup<uint16_t>(vd, vs1, amount, group, start, elems, masked);
      vecRegs_.write(vd, 0, group, int16_t(replacement));
      break;

    case ElementWidth::Word:
      vslideup<uint32_t>(vd, vs1, amount, group, start, elems, masked);
      vecRegs_.write(vd, 0, group, int32_t(replacement));
      break;

    case ElementWidth::Word2:
      vslideup<uint64_t>(vd, vs1, amount, group, start, elems, masked);
      vecRegs_.write(vd, 0, group, int64_t(replacement));
      break;

    case ElementWidth::Word4:
      vslideup<Uint128>(vd, vs1, amount, group, start, elems, masked);
      vecRegs_.write(vd, 0, group, Int128(replacement));
      break;

    case ElementWidth::Word8:
      vslideup<Uint256>(vd, vs1, amount, group, start, elems, masked);
      vecRegs_.write(vd, 0, group, Int256(replacement));
      break;

    case ElementWidth::Word16:
      vslideup<Uint512>(vd, vs1, amount, group, start, elems, masked);
      vecRegs_.write(vd, 0, group, Int512(replacement));
      break;

    case ElementWidth::Word32:
      vslideup<Uint1024>(vd, vs1, amount, group, start, elems, masked);
      vecRegs_.write(vd, 0, group, Int1024(replacement));
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vslidedown(unsigned vd, unsigned vs1, URV amount, unsigned group,
                      unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      unsigned from = ix + amount;
      e1 = 0;
      vecRegs_.read(vs1, from, group, e1);
      dest = e1;
      if (not vecRegs_.write(vd, ix, group, dest))
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVslidedown_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(), rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned dist = vd > vs1 ? vd - vs1 : vs1 - vd;
  if (dist > group)
    {
      illegalInst(di);  // Source/dest vecs cannot overlap
      return;
    }

  URV amount = intRegs_.read(rs2);

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vslidedown<uint8_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Half: vslidedown<uint16_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word: vslidedown<uint32_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word2: vslidedown<uint64_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word4: vslidedown<Uint128>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word8: vslidedown<Uint256>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word16: vslidedown<Uint512>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word32: vslidedown<Uint1024>(vd, vs1, amount, group, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVslidedown_vi(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(), imm = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned dist = vd > vs1 ? vd - vs1 : vs1 - vd;
  if (dist > group)
    {
      illegalInst(di);  // Source/dest vecs cannot overlap
      return;
    }

  URV amount = imm;

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vslidedown<uint8_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Half: vslidedown<uint16_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word: vslidedown<uint32_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word2: vslidedown<uint64_t>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word4: vslidedown<Uint128>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word8: vslidedown<Uint256>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word16: vslidedown<Uint512>(vd, vs1, amount, group, start, elems, masked); break;
    case EW::Word32: vslidedown<Uint1024>(vd, vs1, amount, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmul_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 * e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmul_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmul_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vmul_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vmul_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vmul_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vmul_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vmul_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vmul_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vmul_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmul_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = SRV(intRegs_.read(rs2)), dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 * e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmul_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(), vs1 = di->op1(), rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmul_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vmul_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vmul_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vmul_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vmul_vx<Int128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vmul_vx<Int256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vmul_vx<Int512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vmul_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmulh_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          mulh<ELEM_TYPE>(e1, e2, dest);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmulh_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmulh_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vmulh_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vmulh_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vmulh_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vmulh_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vmulh_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vmulh_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vmulh_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmulh_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = SRV(intRegs_.read(rs2)), dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          mulh(e1, e2, dest);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmulh_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmulh_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vmulh_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vmulh_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vmulh_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vmulh_vx<Int128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vmulh_vx<Int256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vmulh_vx<Int512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vmulh_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVmulhu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmulh_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vmulh_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vmulh_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vmulh_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vmulh_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vmulh_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vmulh_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vmulh_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmulhu_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = intRegs_.read(rs2), dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          mulh(e1, e2, dest);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmulhu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmulhu_vx<uint8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vmulhu_vx<uint16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vmulhu_vx<uint32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vmulhu_vx<uint64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vmulhu_vx<Uint128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vmulhu_vx<Uint256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vmulhu_vx<Uint512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vmulhu_vx<Uint1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmulhsu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                      unsigned start, unsigned elems, bool masked)
{
  typedef typename std::make_unsigned<ELEM_TYPE>::type  U_ELEM_TYPE;

  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;
  U_ELEM_TYPE e2 = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          mulhsu(e1, e2, dest);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmulhsu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmulhsu_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vmulhsu_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vmulhsu_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vmulhsu_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vmulhsu_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vmulhsu_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vmulhsu_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vmulhsu_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmulhsu_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                      unsigned start, unsigned elems, bool masked)
{
  typedef typename std::make_unsigned<ELEM_TYPE>::type  U_ELEM_TYPE;

  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;
  U_ELEM_TYPE e2 = intRegs_.read(rs2);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          mulhsu(e1, e2, dest);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmulhsu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmulhsu_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vmulhsu_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vmulhsu_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vmulhsu_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vmulhsu_vx<Int128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vmulhsu_vx<Int256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vmulhsu_vx<Int512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vmulhsu_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vwmulu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type  ELEM_TYPE_X2;

  ELEM_TYPE e1 = 0, e2 = 0;
  ELEM_TYPE_X2 dest = 0;

  unsigned errors = 0, wideGroup = group * 2;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = ELEM_TYPE_X2(e1);
          dest *= e2;
          if (not vecRegs_.write(vd, ix, wideGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwmulu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  // Double wide legal. Destination register multiple of emul.
  if (not vecRegs_.isDoubleWideLegal(sew, group) or ((vd*8) % (group*2)) != 0)
    {
      illegalInst(di);
      return;
    }

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwmulu_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vwmulu_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vwmulu_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vwmulu_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vwmulu_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vwmulu_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vwmulu_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vwmulu_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type  ELEM_TYPE_X2;

  ELEM_TYPE e1 = 0;
  ELEM_TYPE_X2 dest = 0;
  unsigned errors = 0, wideGroup = group * 2;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = ELEM_TYPE_X2(e1);
          dest *= e2;
          if (not vecRegs_.write(vd, ix, wideGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwmulu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  // Double wide legal. Destination register multiple of emul.
  if (not vecRegs_.isDoubleWideLegal(sew, group) or ((vd*8) % (group*2)) != 0)
    {
      illegalInst(di);
      return;
    }

  SRV e2 = SRV(intRegs_.read(rs2));  // Spec says sign extend. Bogus.

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwmulu_vx<uint8_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Half: vwmulu_vx<uint16_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word: vwmulu_vx<uint32_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word2: vwmulu_vx<uint64_t>(vd, vs1, int64_t(e2), group, start, elems, masked); break;
    case EW::Word4: vwmulu_vx<Uint128>(vd, vs1, Int128(e2), group, start, elems, masked); break;
    case EW::Word8: vwmulu_vx<Uint256>(vd, vs1, Int256(e2), group, start, elems, masked); break;
    case EW::Word16: vwmulu_vx<Uint512>(vd, vs1, Int512(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vwmul_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type  ELEM_TYPE_X2;

  ELEM_TYPE e1 = 0, e2 = 0;
  ELEM_TYPE_X2 dest = 0;
  unsigned errors = 0,  wideGroup = group * 2;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = ELEM_TYPE_X2(e1);
          dest *= ELEM_TYPE_X2(e2);
          if (not vecRegs_.write(vd, ix, wideGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwmul_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  // Double wide legal. Destination register multiple of emul.
  if (not vecRegs_.isDoubleWideLegal(sew, group) or ((vd*8) % (group*2)) != 0)
    {
      illegalInst(di);
      return;
    }

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwmul_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vwmul_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vwmul_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vwmul_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vwmul_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vwmul_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vwmul_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vwmul_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type  ELEM_TYPE_X2;

  ELEM_TYPE e1 = 0;
  ELEM_TYPE_X2 dest = 0;
  ELEM_TYPE_X2 e2Wide(e2);
  unsigned errors = 0, wideGroup = group * 2;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = ELEM_TYPE_X2(e1);
          dest *= e2Wide;
          if (not vecRegs_.write(vd, ix, wideGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwmul_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  // Double wide legal. Destination register multiple of emul.
  if (not vecRegs_.isDoubleWideLegal(sew, group) or ((vd*8) % (group*2)) != 0)
    {
      illegalInst(di);
      return;
    }

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwmul_vx<int8_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Half: vwmul_vx<int16_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word: vwmul_vx<int32_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word2: vwmul_vx<int64_t>(vd, vs1, int64_t(e2), group, start, elems, masked); break;
    case EW::Word4: vwmul_vx<Int128>(vd, vs1, Int128(e2), group, start, elems, masked); break;
    case EW::Word8: vwmul_vx<Int256>(vd, vs1, Int256(e2), group, start, elems, masked); break;
    case EW::Word16: vwmul_vx<Int512>(vd, vs1, Int512(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vwmulsu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type  ELEM_TYPE_X2;
  typedef typename std::make_unsigned<ELEM_TYPE>::type ELEM_TYPE_U;

  ELEM_TYPE e1 = 0;
  ELEM_TYPE_U e2u = 0;
  ELEM_TYPE_X2 dest = 0;
  unsigned errors = 0,  wideGroup = group * 2;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2u))
        {
          dest = ELEM_TYPE_X2(e1);
          ELEM_TYPE_X2 tmp2(e2u);
          dest *= tmp2;
          if (not vecRegs_.write(vd, ix, wideGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwmulsu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  // Double wide legal. Destination register multiple of emul.
  if (not vecRegs_.isDoubleWideLegal(sew, group) or ((vd*8) % (group*2)) != 0)
    {
      illegalInst(di);
      return;
    }

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwmulsu_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vwmulsu_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vwmulsu_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vwmulsu_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vwmulsu_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vwmulsu_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vwmulsu_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vwmulsu_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type  ELEM_TYPE_X2;
  typedef typename std::make_unsigned<ELEM_TYPE>::type ELEM_TYPE_U;

  ELEM_TYPE e1 = 0;
  ELEM_TYPE_X2 dest = 0;
  ELEM_TYPE_U e2u = ELEM_TYPE_U(e2);
  ELEM_TYPE_X2 e2Wide(e2u);

  unsigned errors = 0,  wideGroup = group * 2;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = ELEM_TYPE_X2(e1);
          dest *= e2Wide;
          if (not vecRegs_.write(vd, ix, wideGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwmulsu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  // Double wide legal. Destination register multiple of emul.
  if (not vecRegs_.isDoubleWideLegal(sew, group) or ((vd*8) % (group*2)) != 0)
    {
      illegalInst(di);
      return;
    }

  SRV e2 = SRV(intRegs_.read(rs2));   // Spec says sign extend. Bogus.

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwmulsu_vx<int8_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Half: vwmulsu_vx<int16_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word: vwmulsu_vx<int32_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word2: vwmulsu_vx<int64_t>(vd, vs1, int64_t(e2), group, start, elems, masked); break;
    case EW::Word4: vwmulsu_vx<Int128>(vd, vs1, Int128(e2), group, start, elems, masked); break;
    case EW::Word8: vwmulsu_vx<Int256>(vd, vs1, Int256(e2), group, start, elems, masked); break;
    case EW::Word16: vwmulsu_vx<Int512>(vd, vs1, Int512(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vwmacc_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type DWT; // Double wide type
  unsigned errors = 0, wideGroup = group*2;

  ELEM_TYPE e1 = 0, e2 = 0;
  DWT dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2)
          and vecRegs_.read(vd, ix, wideGroup, dest))
        {
          dest += DWT(e2) * DWT(e2);
          if (not vecRegs_.write(vd, ix, wideGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwmaccu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  // Double wide legal. Destination register multiple of emul.
  if (not vecRegs_.isDoubleWideLegal(sew, group) or ((vd*8) % (group*2)) != 0)
    {
      illegalInst(di);
      return;
    }

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwmacc_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vwmacc_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vwmacc_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vwmacc_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vwmacc_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vwmacc_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vwmacc_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vwmaccu_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                      unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type DWT; // Double wide type
  typedef typename std::make_signed<DWT>::type SDWT; // Signed double wide type
  unsigned errors = 0, doubleGroup = group*2;

  ELEM_TYPE e1 = 0;
  DWT dest = 0;
  SDWT de2 = SDWT(e2);  // sign extend (spec is foolish)

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs1, ix, doubleGroup, dest))
        {
          dest += DWT(e1) * DWT(de2);
          if (not vecRegs_.write(vd, ix, doubleGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwmaccu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  // Double wide legal. Destination register multiple of emul.
  if (not vecRegs_.isDoubleWideLegal(sew, group) or ((vd*8) % (group*2)) != 0)
    {
      illegalInst(di);
      return;
    }

  SRV e2 = SRV(intRegs_.read(rs2));  // Spec says sign extend. Bogus.

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwmaccu_vx<uint8_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Half: vwmaccu_vx<uint16_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word: vwmaccu_vx<uint32_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word2: vwmaccu_vx<uint64_t>(vd, vs1, int64_t(e2), group, start, elems, masked); break;
    case EW::Word4: vwmaccu_vx<Uint128>(vd, vs1, Int128(e2), group, start, elems, masked); break;
    case EW::Word8: vwmaccu_vx<Uint256>(vd, vs1, Int256(e2), group, start, elems, masked); break;
    case EW::Word16: vwmaccu_vx<Uint512>(vd, vs1, Int512(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVwmacc_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  // Double wide legal. Destination register multiple of emul.
  if (not vecRegs_.isDoubleWideLegal(sew, group) or ((vd*8) % (group*2)) != 0)
    {
      illegalInst(di);
      return;
    }

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwmacc_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vwmacc_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vwmacc_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vwmacc_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vwmacc_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vwmacc_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vwmacc_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vwmacc_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type DWT; // Double wide type
  unsigned errors = 0, doubleGroup = group*2;

  ELEM_TYPE e1 = 0;
  DWT dest = 0;
  DWT de2 = DWT(e2);  // sign extend

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs1, ix, doubleGroup, dest))
        {
          dest += DWT(e1) * de2;
          if (not vecRegs_.write(vd, ix, doubleGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwmacc_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  // Double wide legal. Destination register multiple of emul.
  if (not vecRegs_.isDoubleWideLegal(sew, group) or ((vd*8) % (group*2)) != 0)
    {
      illegalInst(di);
      return;
    }

  SRV e2 = SRV(intRegs_.read(rs2));  // Sign extend.

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwmacc_vx<int8_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Half: vwmacc_vx<int16_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word: vwmacc_vx<int32_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word2: vwmacc_vx<int64_t>(vd, vs1, int64_t(e2), group, start, elems, masked); break;
    case EW::Word4: vwmacc_vx<Int128>(vd, vs1, Int128(e2), group, start, elems, masked); break;
    case EW::Word8: vwmacc_vx<Int256>(vd, vs1, Int256(e2), group, start, elems, masked); break;
    case EW::Word16: vwmacc_vx<Int512>(vd, vs1, Int512(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vwmaccsu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                       unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type DWT; // Double wide type
  typedef typename std::make_unsigned<DWT>::type DWTU; // Double wide type unsigned

  unsigned errors = 0, wideGroup = group*2;

  ELEM_TYPE e1 = 0, e2 = 0;
  DWT dest = 0, temp = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2)
          and vecRegs_.read(vd, ix, wideGroup, dest))
        {
          DWTU de2 = DWT(e2);
          mulsu(DWT(e1), de2, temp);
          dest += temp;
          if (not vecRegs_.write(vd, ix, wideGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwmaccsu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  ElementWidth sew = vecRegs_.elemWidth();

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  // Double wide legal. Destination register multiple of emul.
  if (not vecRegs_.isDoubleWideLegal(sew, group) or ((vd*8) % (group*2)) != 0)
    {
      illegalInst(di);
      return;
    }

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwmaccsu_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vwmaccsu_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vwmaccsu_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vwmaccsu_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vwmaccsu_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vwmaccsu_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vwmaccsu_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vwmaccsu_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                       unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type DWT; // Double wide type
  typedef typename std::make_unsigned<DWT>::type DWTU; // Double wide type unsigned

  unsigned errors = 0, wideGroup = group*2;

  ELEM_TYPE e1 = 0;
  DWT de2 = DWT(e2);  // Sign extend.  Spec is bogus.
  DWTU de2u = DWTU(de2); // Then make unsigned,
  DWT dest = 0, temp = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vd, ix, wideGroup, dest))
        {
          mulsu(DWT(e1), de2u, temp);
          dest += temp;
          if (not vecRegs_.write(vd, ix, wideGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwmaccsu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  // Double wide legal. Destination register multiple of emul.
  if (not vecRegs_.isDoubleWideLegal(sew, group) or ((vd*8) % (group*2)) != 0)
    {
      illegalInst(di);
      return;
    }

  SRV e2 = SRV(intRegs_.read(rs2));  // Spec says sign extend. Bogus.

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwmaccsu_vx<int8_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Half: vwmaccsu_vx<int16_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word: vwmaccsu_vx<int32_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word2: vwmaccsu_vx<int64_t>(vd, vs1, int64_t(e2), group, start, elems, masked); break;
    case EW::Word4: vwmaccsu_vx<Int128>(vd, vs1, Int128(e2), group, start, elems, masked); break;
    case EW::Word8: vwmaccsu_vx<Int256>(vd, vs1, Int256(e2), group, start, elems, masked); break;
    case EW::Word16: vwmaccsu_vx<Int512>(vd, vs1, Int512(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vwmaccus_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                       unsigned start, unsigned elems, bool masked)
{
  typedef typename makeDoubleWide<ELEM_TYPE>::type DWT; // Double wide type
  typedef typename std::make_unsigned<ELEM_TYPE>::type ELEM_TYPEU;
  typedef typename std::make_unsigned<DWT>::type DWTU; // Double wide type unsigned

  unsigned errors = 0, wideGroup = group*2;

  ELEM_TYPEU e1 = 0;
  DWT de2 = DWT(e2);  // Sign extend.
  DWT dest = 0, temp = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vd, ix, wideGroup, dest)
          and vecRegs_.read(vd, ix, wideGroup, dest))
        {
          mulsu(de2, DWTU(e1), temp);
          dest += temp;
          if (not vecRegs_.write(vd, ix, wideGroup, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVwmaccus_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  // Double wide legal. Destination register multiple of emul.
  if (not vecRegs_.isDoubleWideLegal(sew, group) or ((vd*8) % (group*2)) != 0)
    {
      illegalInst(di);
      return;
    }

  SRV e2 = SRV(intRegs_.read(rs2));  // Sign extend.

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vwmaccus_vx<int8_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Half: vwmaccus_vx<int16_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word: vwmaccus_vx<int32_t>(vd, vs1, e2, group, start, elems, masked); break;
    case EW::Word2: vwmaccus_vx<int64_t>(vd, vs1, int64_t(e2), group, start, elems, masked); break;
    case EW::Word4: vwmaccus_vx<Int128>(vd, vs1, Int128(e2), group, start, elems, masked); break;
    case EW::Word8: vwmaccus_vx<Int256>(vd, vs1, Int256(e2), group, start, elems, masked); break;
    case EW::Word16: vwmaccus_vx<Int512>(vd, vs1, Int512(e2), group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vdivu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = ~ ELEM_TYPE(0); // divide by zero result
          if (e2 != ELEM_TYPE(0))
            dest = e1 / e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVdivu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vdivu_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vdivu_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vdivu_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vdivu_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vdivu_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vdivu_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vdivu_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vdivu_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vdivu_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  // Spec (sep 24, 2020) says scalar register value should be sign
  // extended. We hope they come to their senses.
  ELEM_TYPE e1 = 0, e2 = intRegs_.read(rs2), dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = ~ ELEM_TYPE(0); // divide by zero result
          if (e2 != ELEM_TYPE(0))
            dest = e1 / e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVdivu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vdivu_vv<uint8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vdivu_vv<uint16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vdivu_vv<uint32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vdivu_vv<uint64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vdivu_vv<Uint128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vdivu_vv<Uint256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vdivu_vv<Uint512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vdivu_vv<Uint1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vdiv_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  int elemBits = integerWidth<ELEM_TYPE> ();
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;
  ELEM_TYPE minInt = ELEM_TYPE(1) << (elemBits - 1);
  ELEM_TYPE negOne = ELEM_TYPE(-1);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = negOne; // Divide by zero result
          if (e2 != ELEM_TYPE(0))
            {
              if (e1 == minInt and e2 == negOne)
                dest = e1;
              else
                dest = e1 / e2;
            }
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVdiv_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vdiv_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vdiv_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vdiv_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vdiv_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vdiv_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vdiv_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vdiv_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vdiv_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vdiv_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  int elemBits = integerWidth<ELEM_TYPE> ();
  ELEM_TYPE e1 = 0, e2 = SRV(intRegs_.read(rs2)), dest = 0;
  ELEM_TYPE minInt = ELEM_TYPE(1) << (elemBits - 1);
  ELEM_TYPE negOne = ELEM_TYPE(-1);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = negOne; // Divide by zero result
          if (e2 != ELEM_TYPE(0))
            {
              if (e1 == minInt and e2 == negOne)
                dest = e1;
              else
                dest = e1 / e2;
            }
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVdiv_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vdiv_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vdiv_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vdiv_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vdiv_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vdiv_vx<Int128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vdiv_vx<Int256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vdiv_vx<Int512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vdiv_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vremu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1; // divide by zero result
          if (e2 != ELEM_TYPE(0))
            dest = e1 % e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVremu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vremu_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vremu_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vremu_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vremu_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vremu_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vremu_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vremu_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vremu_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vremu_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  // Spec (sep 24, 2020) says scalar register value should be sign
  // extended. We hope they come to their senses.
  ELEM_TYPE e1 = 0, e2 = intRegs_.read(rs2), dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1; // divide by zero result
          if (e2 != ELEM_TYPE(0))
            dest = e1 % e2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVremu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vremu_vv<uint8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vremu_vv<uint16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vremu_vv<uint32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vremu_vv<uint64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vremu_vv<Uint128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vremu_vv<Uint256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vremu_vv<Uint512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vremu_vv<Uint1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vrem_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  int elemBits = integerWidth<ELEM_TYPE> ();
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;
  ELEM_TYPE minInt = ELEM_TYPE(1) << (elemBits - 1);
  ELEM_TYPE negOne = ELEM_TYPE(-1);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1; // Divide by zero remainder
          if (e2 != ELEM_TYPE(0))
            {
              if (e1 == minInt and e2 == negOne)
                dest = 0;
              else
                dest = e1 % e2;
            }
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVrem_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vrem_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half: vrem_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word: vrem_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2: vrem_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4: vrem_vv<Int128>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8: vrem_vv<Int256>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vrem_vv<Int512>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vrem_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vrem_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0,  elemBits = integerWidth<ELEM_TYPE> ();
  ELEM_TYPE e1 = 0, e2 = SRV(intRegs_.read(rs2)), dest = 0;
  ELEM_TYPE minInt = ELEM_TYPE(1) << (elemBits - 1);
  ELEM_TYPE negOne = ELEM_TYPE(-1);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1; // Divide by zero remainder
          if (e2 != ELEM_TYPE(0))
            {
              if (e1 == minInt and e2 == negOne)
                dest = 0;
              else
                dest = e1 % e2;
            }
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVrem_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vremu_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Half: vremu_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word: vremu_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word2: vremu_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word4: vremu_vx<Int128>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word8: vremu_vx<Int256>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word16: vremu_vx<Int512>(vd, vs1, rs2, group, start, elems, masked); break;
    case EW::Word32: vremu_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE, typename FROM_TYPE>
void
Hart<URV>::vsext(unsigned vd, unsigned vs1, unsigned group, unsigned fromGroup,
                 unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  FROM_TYPE e1 = 0;
  ELEM_TYPE dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, fromGroup, e1))
        {
          dest = e1;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVsext_vf2(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;


  unsigned group = vecRegs_.groupMultiplierX8();
  unsigned fromGroup = group/2;
  if (fromGroup == 0)
    {
      illegalInst(di);
      return;
    }
  GroupMultiplier emul = GroupMultiplier::One;
  if (not vecRegs_.groupNumberX8ToSymbol(fromGroup, emul))
    {
      illegalInst(di);
      return;
    }

  typedef ElementWidth EW;

  EW sew = vecRegs_.elemWidth();
  EW eew = sew;  // Effective elem width of source.
  switch (sew)
    {
    case EW::Byte: illegalInst(di); return;
    case EW::Half: eew = EW::Byte; break;
    case EW::Word: eew = EW::Half; break;
    case EW::Word2: eew = EW::Word; break;
    case EW::Word4: eew = EW::Word2; break;
    case EW::Word8: eew = EW::Word4; break;
    case EW::Word16: eew = EW::Word8; break;
    case EW::Word32: eew = EW::Word16; break;
    }
  if (not vecRegs_.legalConfig(eew, emul))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  switch (sew)
    {
    case EW::Byte: illegalInst(di); break;
    case EW::Half: vsext<int16_t,int8_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word: vsext<int32_t, int16_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word2: vsext<int64_t, int32_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word4: vsext<Int128, int64_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word8: vsext<Int256, Int128>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word16: vsext<Int512, Int256>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word32: vsext<Int1024, Int512>(vd, vs1, group, fromGroup, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVsext_vf4(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  unsigned fromGroup = group/4;
  if (fromGroup == 0)
    {
      illegalInst(di);
      return;
    }
  GroupMultiplier emul = GroupMultiplier::One;
  if (not vecRegs_.groupNumberX8ToSymbol(fromGroup, emul))
    {
      illegalInst(di);
      return;
    }

  typedef ElementWidth EW;

  EW sew = vecRegs_.elemWidth();
  EW eew = sew;  // Effective elem width of source.
  switch (sew)
    {
    case EW::Byte: illegalInst(di); return;
    case EW::Half: illegalInst(di); return;
    case EW::Word: eew = EW::Byte; break;
    case EW::Word2: eew = EW::Half; break;
    case EW::Word4: eew = EW::Word; break;
    case EW::Word8: eew = EW::Word2; break;
    case EW::Word16: eew = EW::Word4; break;
    case EW::Word32: eew = ElementWidth::Word8; break;
    }
  if (not vecRegs_.legalConfig(eew, emul))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  switch (sew)
    {
    case EW::Byte: illegalInst(di); break;
    case EW::Half: illegalInst(di); break;
    case EW::Word: vsext<int32_t, int8_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word2: vsext<int64_t, int16_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word4: vsext<Int128, int32_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word8: vsext<Int256, int64_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word16: vsext<Int512, Int128>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word32: vsext<Int1024, Int256>(vd, vs1, group, fromGroup, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVsext_vf8(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  unsigned fromGroup = group/8;
  if (fromGroup == 0)
    {
      illegalInst(di);
      return;
    }
  GroupMultiplier emul = GroupMultiplier::One;
  if (not vecRegs_.groupNumberX8ToSymbol(fromGroup, emul))
    {
      illegalInst(di);
      return;
    }

  typedef ElementWidth EW;

  EW sew = vecRegs_.elemWidth();
  EW eew = sew;  // Effective elem width of source.
  switch (sew)
    {
    case EW::Byte: illegalInst(di); return;
    case EW::Half: illegalInst(di); return;
    case EW::Word: illegalInst(di); return;
    case EW::Word2: eew = EW::Byte; break;
    case EW::Word4: eew = EW::Half; break;
    case EW::Word8: eew = EW::Word; break;
    case EW::Word16: eew = EW::Word2; break;
    case EW::Word32: eew = EW::Word4; break;
    }
  if (not vecRegs_.legalConfig(eew, emul))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  switch (sew)
    {
    case EW::Byte: illegalInst(di); return;
    case EW::Half: illegalInst(di); return;
    case EW::Word: illegalInst(di); return;
    case EW::Word2: vsext<int64_t, int8_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word4: vsext<Int128, int16_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word8: vsext<Int256, int32_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word16: vsext<Int512, int64_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word32: vsext<Int1024, Int128>(vd, vs1, group, fromGroup, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE, typename FROM_TYPE>
void
Hart<URV>::vzext(unsigned vd, unsigned vs1, unsigned group, unsigned fromGroup,
                 unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  FROM_TYPE e1 = 0;
  ELEM_TYPE dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, fromGroup, e1))
        {
          dest = e1;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVzext_vf2(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  unsigned fromGroup = group/2;
  if (fromGroup == 0)
    {
      illegalInst(di);
      return;
    }
  GroupMultiplier emul = GroupMultiplier::One;
  if (not vecRegs_.groupNumberX8ToSymbol(fromGroup, emul))
    {
      illegalInst(di);
      return;
    }

  typedef ElementWidth EW;

  EW sew = vecRegs_.elemWidth();
  EW eew = sew;  // Effective elem width of source.
  switch (sew)
    {
    case EW::Byte: illegalInst(di); return;
    case EW::Half: eew = EW::Byte; break;
    case EW::Word: eew = EW::Half; break;
    case EW::Word2: eew = EW::Word; break;
    case EW::Word4: eew = EW::Word2; break;
    case EW::Word8: eew = EW::Word4; break;
    case EW::Word16: eew = EW::Word8; break;
    case EW::Word32: eew = EW::Word16; break;
    }
  if (not vecRegs_.legalConfig(eew, emul))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  switch (sew)
    {
    case EW::Byte: illegalInst(di); break;
    case EW::Half: vzext<uint16_t,uint8_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word: vzext<uint32_t, uint16_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word2: vzext<uint64_t, uint32_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word4: vzext<Uint128, uint64_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word8: vzext<Uint256, Uint128>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word16: vzext<Uint512, Uint256>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word32: vzext<Uint1024, Uint512>(vd, vs1, group, fromGroup, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVzext_vf4(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  unsigned fromGroup = group/4;
  if (fromGroup == 0)
    {
      illegalInst(di);
      return;
    }
  GroupMultiplier emul = GroupMultiplier::One;
  if (not vecRegs_.groupNumberX8ToSymbol(fromGroup, emul))
    {
      illegalInst(di);
      return;
    }

  typedef ElementWidth EW;

  EW sew = vecRegs_.elemWidth();
  EW eew = sew;  // Effective elem width of source.
  switch (sew)
    {
    case EW::Byte: illegalInst(di); return;
    case EW::Half: illegalInst(di); return;
    case EW::Word: eew = EW::Byte; break;
    case EW::Word2: eew = EW::Half; break;
    case EW::Word4: eew = EW::Word; break;
    case EW::Word8: eew = EW::Word2; break;
    case EW::Word16: eew = EW::Word4; break;
    case EW::Word32: eew = EW::Word8; break;
    }
  if (not vecRegs_.legalConfig(eew, emul))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  switch (sew)
    {
    case EW::Byte: illegalInst(di); break;
    case EW::Half: illegalInst(di); break;
    case EW::Word: vzext<uint32_t, uint8_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word2: vzext<uint64_t, uint16_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word4: vzext<Uint128, uint32_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word8: vzext<Uint256, uint64_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word16: vzext<Uint512, Uint128>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word32: vzext<Uint1024, Uint256>(vd, vs1, group, fromGroup, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVzext_vf8(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  unsigned group = vecRegs_.groupMultiplierX8();
  unsigned fromGroup = group/8;
  if (fromGroup == 0)
    {
      illegalInst(di);
      return;
    }
  GroupMultiplier emul = GroupMultiplier::One;
  if (not vecRegs_.groupNumberX8ToSymbol(fromGroup, emul))
    {
      illegalInst(di);
      return;
    }

  typedef ElementWidth EW;

  EW sew = vecRegs_.elemWidth();
  EW eew = sew;  // Effective elem width of source.
  switch (sew)
    {
    case EW::Byte: illegalInst(di); return;
    case EW::Half: illegalInst(di); return;
    case EW::Word: illegalInst(di); return;
    case EW::Word2: eew = EW::Byte; break;
    case EW::Word4: eew = EW::Half; break;
    case EW::Word8: eew = EW::Word; break;
    case EW::Word16: eew = EW::Word2; break;
    case EW::Word32: eew = EW::Word4; break;
    }
  if (not vecRegs_.legalConfig(eew, emul))
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1();
  unsigned elems = vecRegs_.elemCount(), start = vecRegs_.startIndex();

  switch (sew)
    {
    case EW::Byte: illegalInst(di); break;
    case EW::Half: illegalInst(di); break;
    case EW::Word: illegalInst(di); break;
    case EW::Word2: vzext<uint64_t, uint8_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word4: vzext<Uint128, uint16_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word8: vzext<Uint256, uint32_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word16: vzext<Uint512, uint64_t>(vd, vs1, group, fromGroup, start, elems, masked); break;
    case EW::Word32: vzext<Uint1024, Uint128>(vd, vs1, group, fromGroup, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vadc_vvm(unsigned vd, unsigned vs1, unsigned vs2, unsigned vcin,
                    unsigned group, unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 + e2;
          if (vecRegs_.isActive(vcin, ix))
            dest += ELEM_TYPE(1);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vadc_vxm(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned vcin,
                    unsigned group, unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 + e2;
          if (vecRegs_.isActive(vcin, ix))
            dest += ELEM_TYPE(1);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vsbc_vvm(unsigned vd, unsigned vs1, unsigned vs2, unsigned vbin,
                    unsigned group, unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 - e2;
          if (vecRegs_.isActive(vbin, ix))
            dest -= ELEM_TYPE(1);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vsbc_vxm(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned vbin,
                    unsigned group, unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 - e2;
          if (vecRegs_.isActive(vbin, ix))
            dest -= ELEM_TYPE(1);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmadc_vvm(unsigned vcout, unsigned vs1, unsigned vs2, bool carry, unsigned vcin,
                     unsigned group, unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 + e2;
          if (carry and vecRegs_.isActive(vcin, ix))
            dest += ELEM_TYPE(1);

          bool cout = dest < e1;
          if (not vecRegs_.setMaskRegister(vcout, ix, cout))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmadc_vxm(unsigned vcout, unsigned vs1, ELEM_TYPE e2, bool carry, unsigned vcin,
                     unsigned group, unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 + e2;
          if (carry and vecRegs_.isActive(vcin, ix))
            dest += ELEM_TYPE(1);

          bool cout = dest < e1;
          if (not vecRegs_.setMaskRegister(vcout, ix, cout))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmsbc_vvm(unsigned vbout, unsigned vs1, unsigned vs2, bool borrow, unsigned vbin,
                     unsigned group, unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 - e2;
          if (borrow and vecRegs_.isActive(vbin, ix))
            dest -= ELEM_TYPE(1);

          bool bout = e1 < e2;
          if (not vecRegs_.setMaskRegister(vbout, ix, bout))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmsbc_vxm(unsigned vbout, unsigned vs1, ELEM_TYPE e2, bool borrow, unsigned vbin,
                     unsigned group, unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 - e2;
          if (borrow and vecRegs_.isActive(vbin, ix))
            dest -= ELEM_TYPE(1);

          bool bout = e1 < e2;
          if (not vecRegs_.setMaskRegister(vbout, ix, bout))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVadc_vvm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2(),  vcin = 0;
  if (vd == vcin or not masked)   // cannot overlap vcin, unmasked verion reserved
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vadc_vvm<uint8_t>(vd, vs1, vs2, vcin, group, start, elems); break;
    case EW::Half: vadc_vvm<uint16_t>(vd, vs1, vs2, vcin, group, start, elems); break;
    case EW::Word: vadc_vvm<uint32_t>(vd, vs1, vs2, vcin, group, start, elems); break;
    case EW::Word2: vadc_vvm<uint64_t>(vd, vs1, vs2, vcin, group, start, elems); break;
    case EW::Word4: vadc_vvm<Uint128>(vd, vs1, vs2, vcin, group, start, elems); break;
    case EW::Word8: vadc_vvm<Uint256>(vd, vs1, vs2, vcin, group, start, elems); break;
    case EW::Word16: vadc_vvm<Uint512>(vd, vs1, vs2, vcin, group, start, elems); break;
    case EW::Word32: vadc_vvm<Uint1024>(vd, vs1, vs2, vcin, group, start, elems); break;
    }
}


template <typename URV>
void
Hart<URV>::execVadc_vxm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vcin = 0;
  if (vd == vcin or not masked)   // cannot overlap vcin, unmasked verion reserved
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(di->op2()));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vadc_vxm<uint8_t>(vd, vs1, e2, vcin, group, start, elems); break;
    case EW::Half: vadc_vxm<uint16_t>(vd, vs1, e2, vcin, group, start, elems); break;
    case EW::Word: vadc_vxm<uint32_t>(vd, vs1, int32_t(e2), vcin, group, start, elems); break;
    case EW::Word2: vadc_vxm<uint64_t>(vd, vs1, int64_t(e2), vcin, group, start, elems); break;
    case EW::Word4: vadc_vxm<Uint128>(vd, vs1, Int128(e2), vcin, group, start, elems); break;
    case EW::Word8: vadc_vxm<Uint256>(vd, vs1, Int256(e2), vcin, group, start, elems); break;
    case EW::Word16: vadc_vxm<Uint512>(vd, vs1, Int512(e2), vcin, group, start, elems); break;
    case EW::Word32: vadc_vxm<Uint1024>(vd, vs1, Int1024(e2), vcin, group, start, elems); break;
    }
}


template <typename URV>
void
Hart<URV>::execVadc_vim(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vcin = 0;
  if (vd == vcin or not masked)   // cannot overlap vcin, unmasked verion reserved
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = di->op2As<int32_t>();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vadc_vxm<uint8_t>(vd, vs1, e2, vcin, group, start, elems); break;
    case EW::Half: vadc_vxm<uint16_t>(vd, vs1, e2, vcin, group, start, elems); break;
    case EW::Word: vadc_vxm<uint32_t>(vd, vs1, int32_t(e2), vcin, group, start, elems); break;
    case EW::Word2: vadc_vxm<uint64_t>(vd, vs1, int64_t(e2), vcin, group, start, elems); break;
    case EW::Word4: vadc_vxm<Uint128>(vd, vs1, Int128(e2), vcin, group, start, elems); break;
    case EW::Word8: vadc_vxm<Uint256>(vd, vs1, Int256(e2), vcin, group, start, elems); break;
    case EW::Word16: vadc_vxm<Uint512>(vd, vs1, Int512(e2), vcin, group, start, elems); break;
    case EW::Word32: vadc_vxm<Uint1024>(vd, vs1, Int1024(e2), vcin, group, start, elems); break;
    }
}


template <typename URV>
void
Hart<URV>::execVsbc_vvm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2(),  vbin = 0;
  if (vd == vbin or not masked)   // cannot overlap borrow-in, unmasked verion reserved
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vsbc_vvm<uint8_t>(vd, vs1, vs2, vbin, group, start, elems); break;
    case EW::Half: vsbc_vvm<uint16_t>(vd, vs1, vs2, vbin, group, start, elems); break;
    case EW::Word: vsbc_vvm<uint32_t>(vd, vs1, vs2, vbin, group, start, elems); break;
    case EW::Word2: vsbc_vvm<uint64_t>(vd, vs1, vs2, vbin, group, start, elems); break;
    case EW::Word4: vsbc_vvm<Uint128>(vd, vs1, vs2, vbin, group, start, elems); break;
    case EW::Word8: vsbc_vvm<Uint256>(vd, vs1, vs2, vbin, group, start, elems); break;
    case EW::Word16: vsbc_vvm<Uint512>(vd, vs1, vs2, vbin, group, start, elems); break;
    case EW::Word32: vsbc_vvm<Uint1024>(vd, vs1, vs2, vbin, group, start, elems); break;
    }
}


template <typename URV>
void
Hart<URV>::execVsbc_vxm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vbin = 0;
  if (vd == vbin or not masked)   // cannot overlap borrow-in, unmasked verion reserved
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(di->op2()));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vsbc_vxm<uint8_t>(vd, vs1, e2, vbin, group, start, elems); break;
    case EW::Half: vsbc_vxm<uint16_t>(vd, vs1, e2, vbin, group, start, elems); break;
    case EW::Word: vsbc_vxm<uint32_t>(vd, vs1, int32_t(e2), vbin, group, start, elems); break;
    case EW::Word2: vsbc_vxm<uint64_t>(vd, vs1, int64_t(e2), vbin, group, start, elems); break;
    case EW::Word4: vsbc_vxm<Uint128>(vd, vs1, Int128(e2), vbin, group, start, elems); break;
    case EW::Word8: vsbc_vxm<Uint256>(vd, vs1, Int256(e2), vbin, group, start, elems); break;
    case EW::Word16: vsbc_vxm<Uint512>(vd, vs1, Int512(e2), vbin, group, start, elems); break;
    case EW::Word32: vsbc_vxm<Uint1024>(vd, vs1, Int1024(e2), vbin, group, start, elems); break;
    }
}


template <typename URV>
void
Hart<URV>::execVmadc_vvm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool carry = di->isMasked();
  unsigned vcout = di->op0(),  vs1 = di->op1(),  vs2 = di->op2(),  vcin = 0;
  if (vcout == vcin)   // cannot overlap vcin
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmadc_vvm<uint8_t>(vcout, vs1, vs2, carry, vcin, group, start, elems); break;
    case EW::Half: vmadc_vvm<uint16_t>(vcout, vs1, vs2, carry, vcin, group, start, elems); break;
    case EW::Word: vmadc_vvm<uint32_t>(vcout, vs1, vs2, carry, vcin, group, start, elems); break;
    case EW::Word2: vmadc_vvm<uint64_t>(vcout, vs1, vs2, carry, vcin, group, start, elems); break;
    case EW::Word4: vmadc_vvm<Uint128>(vcout, vs1, vs2, carry, vcin, group, start, elems); break;
    case EW::Word8: vmadc_vvm<Uint256>(vcout, vs1, vs2, carry, vcin, group, start, elems); break;
    case EW::Word16: vmadc_vvm<Uint512>(vcout, vs1, vs2, carry, vcin, group, start, elems); break;
    case EW::Word32: vmadc_vvm<Uint1024>(vcout, vs1, vs2, carry, vcin, group, start, elems); break;
    }
}


template <typename URV>
void
Hart<URV>::execVmadc_vxm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool carry = di->isMasked();
  unsigned vcout = di->op0(),  vs1 = di->op1(),  vcin = 0;
  if (vcout == vcin)   // cannot overlap vcin
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(di->op2()));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmadc_vxm<uint8_t>(vcout, vs1, e2, carry, vcin, group, start, elems); break;
    case EW::Half: vmadc_vxm<uint16_t>(vcout, vs1, e2, carry, vcin, group, start, elems); break;
    case EW::Word: vmadc_vxm<uint32_t>(vcout, vs1, int32_t(e2), carry, vcin, group, start, elems); break;
    case EW::Word2: vmadc_vxm<uint64_t>(vcout, vs1, int64_t(e2), carry, vcin, group, start, elems); break;
    case EW::Word4: vmadc_vxm<Uint128>(vcout, vs1, Int128(e2), carry, vcin, group, start, elems); break;
    case EW::Word8: vmadc_vxm<Uint256>(vcout, vs1, Int256(e2), carry, vcin, group, start, elems); break;
    case EW::Word16: vmadc_vxm<Uint512>(vcout, vs1, Int512(e2), carry, vcin, group, start, elems); break;
    case EW::Word32: vmadc_vxm<Uint1024>(vcout, vs1, Int1024(e2), carry, vcin, group, start, elems); break;
    }
}


template <typename URV>
void
Hart<URV>::execVmadc_vim(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool carry = di->isMasked();
  unsigned vcout = di->op0(),  vs1 = di->op1(),  vcin = 0;
  if (vcout == vcin)   // cannot overlap vcin
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = di->op2As<int32_t>();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmadc_vxm<uint8_t>(vcout, vs1, e2, carry, vcin, group, start, elems); break;
    case EW::Half: vmadc_vxm<uint16_t>(vcout, vs1, e2, carry, vcin, group, start, elems); break;
    case EW::Word: vmadc_vxm<uint32_t>(vcout, vs1, int32_t(e2), carry, vcin, group, start, elems); break;
    case EW::Word2: vmadc_vxm<uint64_t>(vcout, vs1, int64_t(e2), carry, vcin, group, start, elems); break;
    case EW::Word4: vmadc_vxm<Uint128>(vcout, vs1, Int128(e2), carry, vcin, group, start, elems); break;
    case EW::Word8: vmadc_vxm<Uint256>(vcout, vs1, Int256(e2), carry, vcin, group, start, elems); break;
    case EW::Word16: vmadc_vxm<Uint512>(vcout, vs1, Int512(e2), carry, vcin, group, start, elems); break;
    case EW::Word32: vmadc_vxm<Uint1024>(vcout, vs1, Int1024(e2), carry, vcin, group, start, elems); break;
    }
}


template <typename URV>
void
Hart<URV>::execVmsbc_vvm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool borrow = di->isMasked();
  unsigned vbout = di->op0(),  vs1 = di->op1(),  vs2 = di->op2(),  vbin = 0;
  if (vbout == vbin)   // cannot overlap borrow-in
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmsbc_vvm<uint8_t>(vbout, vs1, vs2, borrow, vbin, group, start, elems); break;
    case EW::Half: vmsbc_vvm<uint16_t>(vbout, vs1, vs2, borrow, vbin, group, start, elems); break;
    case EW::Word: vmsbc_vvm<uint32_t>(vbout, vs1, vs2, borrow, vbin, group, start, elems); break;
    case EW::Word2: vmsbc_vvm<uint64_t>(vbout, vs1, vs2, borrow, vbin, group, start, elems); break;
    case EW::Word4: vmsbc_vvm<Uint128>(vbout, vs1, vs2, borrow, vbin, group, start, elems); break;
    case EW::Word8: vmsbc_vvm<Uint256>(vbout, vs1, vs2, borrow, vbin, group, start, elems); break;
    case EW::Word16: vmsbc_vvm<Uint512>(vbout, vs1, vs2, borrow, vbin, group, start, elems); break;
    case EW::Word32: vmsbc_vvm<Uint1024>(vbout, vs1, vs2, borrow, vbin, group, start, elems); break;
    }
}


template <typename URV>
void
Hart<URV>::execVmsbc_vxm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool borrow = di->isMasked();
  unsigned vbout = di->op0(),  vs1 = di->op1(),  vbin = 0;
  if (vbout == vbin)   // cannot overlap borrow-in
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(di->op2()));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmsbc_vxm<uint8_t>(vbout, vs1, e2, borrow, vbin, group, start, elems); break;
    case EW::Half: vmsbc_vxm<uint16_t>(vbout, vs1, e2, borrow, vbin, group, start, elems); break;
    case EW::Word: vmsbc_vxm<uint32_t>(vbout, vs1, int32_t(e2), borrow, vbin, group, start, elems); break;
    case EW::Word2: vmsbc_vxm<uint64_t>(vbout, vs1, int64_t(e2), borrow, vbin, group, start, elems); break;
    case EW::Word4: vmsbc_vxm<Uint128>(vbout, vs1, Int128(e2), borrow, vbin, group, start, elems); break;
    case EW::Word8: vmsbc_vxm<Uint256>(vbout, vs1, Int256(e2), borrow, vbin, group, start, elems); break;
    case EW::Word16: vmsbc_vxm<Uint512>(vbout, vs1, Int512(e2), borrow, vbin, group, start, elems); break;
    case EW::Word32: vmsbc_vxm<Uint1024>(vbout, vs1, Int1024(e2), borrow, vbin, group, start, elems); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmerge_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                     unsigned start, unsigned elems)
{
  unsigned errors = 0;

  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = vecRegs_.isActive(0, ix) ? e2 : e1;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmerge_vv(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (not masked)
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmerge_vv<int8_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Half: vmerge_vv<int16_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word: vmerge_vv<int32_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word2: vmerge_vv<int64_t>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word4: vmerge_vv<Int128>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word8: vmerge_vv<Int256>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word16: vmerge_vv<Int512>(vd, vs1, vs2, group, start, elems); break;
    case EW::Word32: vmerge_vv<Int1024>(vd, vs1, vs2, group, start, elems); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmerge_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                     unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = vecRegs_.isActive(0, ix) ? e2 : e1;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmerge_vx(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  if (not masked)
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmerge_vx<int8_t>(vd, vs1, e2, group, start, elems); break;
    case EW::Half: vmerge_vx<int16_t>(vd, vs1, e2, group, start, elems); break;
    case EW::Word: vmerge_vx<int32_t>(vd, vs1, e2, group, start, elems); break;
    case EW::Word2: vmerge_vx<int64_t>(vd, vs1, e2, group, start, elems); break;
    case EW::Word4: vmerge_vx<Int128>(vd, vs1, Int128(e2), group, start, elems); break;
    case EW::Word8: vmerge_vx<Int256>(vd, vs1, Int256(e2), group, start, elems); break;
    case EW::Word16: vmerge_vx<Int512>(vd, vs1, Int512(e2), group, start, elems); break;
    case EW::Word32: vmerge_vx<Int1024>(vd, vs1, Int1024(e2), group, start, elems); break;
    }
}


template <typename URV>
void
Hart<URV>::execVmerge_vi(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  int32_t imm = di->op2As<int32_t>();
  if (not masked)
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmerge_vx<int8_t>(vd, vs1, imm, group, start, elems); break;
    case EW::Half: vmerge_vx<int16_t>(vd, vs1, imm, group, start, elems); break;
    case EW::Word: vmerge_vx<int32_t>(vd, vs1, imm, group, start, elems); break;
    case EW::Word2: vmerge_vx<int64_t>(vd, vs1, imm, group, start, elems); break;
    case EW::Word4: vmerge_vx<Int128>(vd, vs1, Int128(imm), group, start, elems); break;
    case EW::Word8: vmerge_vx<Int256>(vd, vs1, Int256(imm), group, start, elems); break;
    case EW::Word16: vmerge_vx<Int512>(vd, vs1, Int512(imm), group, start, elems); break;
    case EW::Word32: vmerge_vx<Int1024>(vd, vs1, Int1024(imm), group, start, elems); break;
    }
}


template <typename URV>
void
Hart<URV>::execVmv_x_s(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig() or di->isMasked())
    {
      illegalInst(di);
      return;
    }

  unsigned rd = di->op0(), vs1 = di->op1(), groupX8 = 8;

  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      {
        int8_t val = 0;
        vecRegs_.read(vs1, 0, groupX8, val);
        intRegs_.write(rd, SRV(val));
      }
      break;

    case ElementWidth::Half:
      {
        int16_t val = 0;
        vecRegs_.read(vs1, 0, groupX8, val);
        intRegs_.write(rd, SRV(val));
      }
      break;

    case ElementWidth::Word:
      {
        uint32_t val = 0;
        vecRegs_.read(vs1, 0, groupX8, val);
        intRegs_.write(rd, SRV(val));
      }
      break;

    case ElementWidth::Word2:
      {
        int64_t val = 0;
        vecRegs_.read(vs1, 0, groupX8, val);
        intRegs_.write(rd, SRV(val));
      }
      break;

    case ElementWidth::Word4:
      {
        Int128 val = 0;
        vecRegs_.read(vs1, 0, groupX8, val);
        intRegs_.write(rd, SRV(val));
      }
      break;

    case ElementWidth::Word8:
      {
        Int256 val = 0;
        vecRegs_.read(vs1, 0, groupX8, val);
        intRegs_.write(rd, SRV(val));
      }
      break;

    case ElementWidth::Word16:
      {
        Int512 val = 0;
        vecRegs_.read(vs1, 0, groupX8, val);
        intRegs_.write(rd, SRV(val));
      }
      break;

    case ElementWidth::Word32:
      {
        Int1024 val = 0;
        vecRegs_.read(vs1, 0, groupX8, val);
        intRegs_.write(rd, SRV(val));
      }
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVmv_s_x(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig() or di->isMasked())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(), rs1 = di->op1(), groupX8 = 8;

  ElementWidth sew = vecRegs_.elemWidth();

  SRV val = intRegs_.read(rs1);

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vecRegs_.write(vd, 0, groupX8, int8_t(val)); break;
    case EW::Half: vecRegs_.write(vd, 0, groupX8, int16_t(val)); break;
    case EW::Word: vecRegs_.write(vd, 0, groupX8, int32_t(val)); break;
    case EW::Word2: vecRegs_.write(vd, 0, groupX8, int64_t(val)); break;
    case EW::Word4: vecRegs_.write(vd, 0, groupX8, Int128(val)); break;
    case EW::Word8: vecRegs_.write(vd, 0, groupX8, Int256(val)); break;
    case EW::Word16: vecRegs_.write(vd, 0, groupX8, Int512(val)); break;
    case EW::Word32: vecRegs_.write(vd, 0, groupX8, Int1024(val)); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmv_v_v(unsigned vd, unsigned vs1, unsigned group,
                   unsigned start, unsigned elems)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmv_v_v(const DecodedInst* di)
{
  if (di->isMasked() or (not isVecLegal()) or (not vecRegs_.legalConfig()))
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmv_v_v<int8_t>(vd, vs1, group, start, elems); break;
    case EW::Half: vmv_v_v<int16_t>(vd, vs1, group, start, elems); break;
    case EW::Word: vmv_v_v<int32_t>(vd, vs1, group, start, elems); break;
    case EW::Word2: vmv_v_v<int64_t>(vd, vs1, group, start, elems); break;
    case EW::Word4: vmv_v_v<Int128>(vd, vs1, group, start, elems); break;
    case EW::Word8: vmv_v_v<Int256>(vd, vs1, group, start, elems); break;
    case EW::Word16: vmv_v_v<Int512>(vd, vs1, group, start, elems); break;
    case EW::Word32: vmv_v_v<Int1024>(vd, vs1, group, start, elems); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vmv_v_x(unsigned vd, ELEM_TYPE e1, unsigned group,
                   unsigned start, unsigned elems)
{
  unsigned errors = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    if (not vecRegs_.write(vd, ix, group, e1))
      errors++;

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVmv_v_x(const DecodedInst* di)
{
  if (di->isMasked() or (not isVecLegal()) or not (vecRegs_.legalConfig()))
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0();
  unsigned rs1 = di->op1();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  int e1 = SRV(intRegs_.read(rs1));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmv_v_x<int8_t>(vd, e1, group, start, elems); break;
    case EW::Half: vmv_v_x<int16_t>(vd, e1, group, start, elems); break;
    case EW::Word: vmv_v_x<int32_t>(vd, e1, group, start, elems); break;
    case EW::Word2: vmv_v_x<int64_t>(vd, e1, group, start, elems); break;
    case EW::Word4: vmv_v_x<Int128>(vd, Int128(e1), group, start, elems); break;
    case EW::Word8: vmv_v_x<Int256>(vd, Int256(e1), group, start, elems); break;
    case EW::Word16: vmv_v_x<Int512>(vd, Int512(e1), group, start, elems); break;
    case EW::Word32: vmv_v_x<Int1024>(vd, Int1024(e1), group, start, elems); break;
    }
}


template <typename URV>
void
Hart<URV>::execVmv_v_i(const DecodedInst* di)
{
  if (di->isMasked() or (not isVecLegal()) or (not vecRegs_.legalConfig()))
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  int e1 = 0;

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte: vmv_v_x<int8_t>(vd, e1, group, start, elems); break;
    case EW::Half: vmv_v_x<int16_t>(vd, e1, group, start, elems); break;
    case EW::Word: vmv_v_x<int32_t>(vd, e1, group, start, elems); break;
    case EW::Word2: vmv_v_x<int64_t>(vd, e1, group, start, elems); break;
    case EW::Word4: vmv_v_x<Int128>(vd, Int128(e1), group, start, elems); break;
    case EW::Word8: vmv_v_x<Int256>(vd, Int256(e1), group, start, elems); break;
    case EW::Word16: vmv_v_x<Int512>(vd, Int512(e1), group, start, elems); break;
    case EW::Word32: vmv_v_x<Int1024>(vd, Int1024(e1), group, start, elems); break;
    }
}


template <typename URV>
void
Hart<URV>::execVmv1r_v(const DecodedInst* di)
{
  if (di->isMasked() or (not isVecLegal()) or (not vecRegs_.legalConfig()))
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(), vs1 = di->op1();
  if (vd == vs1)
    return;

  unsigned bytes = vecRegs_.bytesPerRegister();

  uint8_t* dest = vecRegs_.getVecData(vd);
  uint8_t* source = vecRegs_.getVecData(vs1);
  assert(dest);
  assert(source);

  memcpy(dest, source, bytes);

  vecRegs_.setLastWrittenReg(vd, bytes-1, 8);
}


template <typename URV>
void
Hart<URV>::execVmv2r_v(const DecodedInst* di)
{
  if (di->isMasked() or (not isVecLegal()) or (not vecRegs_.legalConfig()))
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(), vs1 = di->op1();
  if ((vd & 1) != 0 or (vs1 & 1) != 0)
    {
      illegalInst(di);   // Vec indices must be even
      return;
    }

  if (vd == vs1)
    return;

  unsigned bytes = vecRegs_.bytesPerRegister() * 2;

  uint8_t* dest = vecRegs_.getVecData(vd);
  uint8_t* source = vecRegs_.getVecData(vs1);
  assert(dest);
  assert(source);

  memcpy(dest, source, bytes);

  vecRegs_.setLastWrittenReg(vd, bytes-1, 8);
}


template <typename URV>
void
Hart<URV>::execVmv4r_v(const DecodedInst* di)
{
  if (di->isMasked() or (not isVecLegal()) or (not vecRegs_.legalConfig()))
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(), vs1 = di->op1();
  if ((vd & 3) != 0 or (vs1 & 3) != 0)
    {
      illegalInst(di);   // Vec indices must be multiples of 4
      return;
    }

  if (vd == vs1)
    return;

  unsigned bytes = vecRegs_.bytesPerRegister() * 4;

  uint8_t* dest = vecRegs_.getVecData(vd);
  uint8_t* source = vecRegs_.getVecData(vs1);
  assert(dest);
  assert(source);

  memcpy(dest, source, bytes);

  vecRegs_.setLastWrittenReg(vd, bytes-1, 8);
}


template <typename URV>
void
Hart<URV>::execVmv8r_v(const DecodedInst* di)
{
  if (di->isMasked() or (not isVecLegal()) or (not vecRegs_.legalConfig()))
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(), vs1 = di->op1();
  if ((vd & 7) != 0 or (vs1 & 7) != 0)
    {
      illegalInst(di);   // Vec indices must be multiples of 8
      return;
    }

  if (vd == vs1)
    return;

  unsigned bytes = vecRegs_.bytesPerRegister() * 8;

  uint8_t* dest = vecRegs_.getVecData(vd);
  uint8_t* source = vecRegs_.getVecData(vs1);
  assert(dest);
  assert(source);

  memcpy(dest, source, bytes);

  vecRegs_.setLastWrittenReg(vd, bytes-1, 8);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vsaddu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  ELEM_TYPE maxVal = ~ ELEM_TYPE(0);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 + e2;
          if (dest < e1)
            {
              dest = maxVal;
              csRegs_.write(CsrNumber::VXSAT, PrivilegeMode::Machine, 1);
            }
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVsaddu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vsaddu_vv<uint8_t> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half:   vsaddu_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word:   vsaddu_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2:  vsaddu_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4:  vsaddu_vv<Uint128> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8:  vsaddu_vv<Uint256> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vsaddu_vv<Uint512> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vsaddu_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vsaddu_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;

  ELEM_TYPE maxVal = ~ ELEM_TYPE(0);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 + e2;
          if (dest < e1)
            {
              dest = maxVal;
              csRegs_.write(CsrNumber::VXSAT, PrivilegeMode::Machine, 1);
            }
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVsaddu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vsaddu_vx<uint8_t> (vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Half:   vsaddu_vx<uint16_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word:   vsaddu_vx<uint32_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word2:  vsaddu_vx<uint64_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word4:  vsaddu_vx<Uint128> (vd, vs1, Int128(e2),  group, start, elems, masked); break;
    case EW::Word8:  vsaddu_vx<Uint256> (vd, vs1, Int256(e2),  group, start, elems, masked); break;
    case EW::Word16: vsaddu_vx<Uint512> (vd, vs1, Int512(e2),  group, start, elems, masked); break;
    case EW::Word32: vsaddu_vx<Uint1024>(vd, vs1, Int1024(e2), group, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVsaddu_vi(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(), vs1 = di->op1();
  int32_t imm = di->op2As<int32_t>();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vsaddu_vx<uint8_t> (vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Half:   vsaddu_vx<uint16_t>(vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Word:   vsaddu_vx<uint32_t>(vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Word2:  vsaddu_vx<uint64_t>(vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Word4:  vsaddu_vx<Uint128> (vd, vs1, Int128(imm),  group, start, elems, masked); break;
    case EW::Word8:  vsaddu_vx<Uint256> (vd, vs1, Int256(imm),  group, start, elems, masked); break;
    case EW::Word16: vsaddu_vx<Uint512> (vd, vs1, Int512(imm),  group, start, elems, masked); break;
    case EW::Word32: vsaddu_vx<Uint1024>(vd, vs1, Int1024(imm), group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vsadd_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  uint32_t bitCount = sizeof(ELEM_TYPE)*8;
  ELEM_TYPE minVal = ELEM_TYPE(1) << (bitCount - 1);
  ELEM_TYPE maxVal = minVal - 1;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 + e2;
          bool sameSign = (e1 < 0) == (e2 < 0);
          if (sameSign and ((e1 < 0) != (dest < 0)))
            {
              if (e1 < 0)
                dest = minVal;
              else
                dest = maxVal;
              csRegs_.write(CsrNumber::VXSAT, PrivilegeMode::Machine, 1);
            }
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVsadd_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vsadd_vv<int8_t> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half:   vsadd_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word:   vsadd_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2:  vsadd_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4:  vsadd_vv<Int128> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8:  vsadd_vv<Int256> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vsadd_vv<Int512> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vsadd_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vsadd_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;

  uint32_t bitCount = sizeof(ELEM_TYPE)*8;
  ELEM_TYPE minVal = ELEM_TYPE(1) << (bitCount - 1);
  ELEM_TYPE maxVal = minVal - 1;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 + e2;
          bool sameSign = (e1 < 0) == (e2 < 0);
          if (sameSign and ((e1 < 0) != (dest < 0)))
            {
              if (e1 < 0)
                dest = minVal;
              else
                dest = maxVal;
              csRegs_.write(CsrNumber::VXSAT, PrivilegeMode::Machine, 1);
            }
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVsadd_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vsadd_vx<int8_t> (vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Half:   vsadd_vx<int16_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word:   vsadd_vx<int32_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word2:  vsadd_vx<int64_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word4:  vsadd_vx<Int128> (vd, vs1, Int128(e2),  group, start, elems, masked); break;
    case EW::Word8:  vsadd_vx<Int256> (vd, vs1, Int256(e2),  group, start, elems, masked); break;
    case EW::Word16: vsadd_vx<Int512> (vd, vs1, Int512(e2),  group, start, elems, masked); break;
    case EW::Word32: vsadd_vx<Int1024>(vd, vs1, Int1024(e2), group, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVsadd_vi(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(), vs1 = di->op1();
  int32_t imm = di->op2As<int32_t>();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vsadd_vx<int8_t> (vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Half:   vsadd_vx<int16_t>(vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Word:   vsadd_vx<int32_t>(vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Word2:  vsadd_vx<int64_t>(vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Word4:  vsadd_vx<Int128> (vd, vs1, Int128(imm),  group, start, elems, masked); break;
    case EW::Word8:  vsadd_vx<Int256> (vd, vs1, Int256(imm),  group, start, elems, masked); break;
    case EW::Word16: vsadd_vx<Int512> (vd, vs1, Int512(imm),  group, start, elems, masked); break;
    case EW::Word32: vsadd_vx<Int1024>(vd, vs1, Int1024(imm), group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vssubu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  ELEM_TYPE minVal = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 - e2;
          if (dest > e1)
            {
              dest = minVal;
              csRegs_.write(CsrNumber::VXSAT, PrivilegeMode::Machine, 1);
            }
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVssubu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vssubu_vv<uint8_t> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half:   vssubu_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word:   vssubu_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2:  vssubu_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4:  vssubu_vv<Uint128> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8:  vssubu_vv<Uint256> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vssubu_vv<Uint512> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vssubu_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vssubu_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;

  ELEM_TYPE minVal = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 - e2;
          if (dest > e1)
            {
              dest = minVal;
              csRegs_.write(CsrNumber::VXSAT, PrivilegeMode::Machine, 1);
            }
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVssubu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vssubu_vx<uint8_t> (vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Half:   vssubu_vx<uint16_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word:   vssubu_vx<uint32_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word2:  vssubu_vx<uint64_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word4:  vssubu_vx<Uint128> (vd, vs1, Int128(e2),  group, start, elems, masked); break;
    case EW::Word8:  vssubu_vx<Uint256> (vd, vs1, Int256(e2),  group, start, elems, masked); break;
    case EW::Word16: vssubu_vx<Uint512> (vd, vs1, Int512(e2),  group, start, elems, masked); break;
    case EW::Word32: vssubu_vx<Uint1024>(vd, vs1, Int1024(e2), group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vssub_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  uint32_t bitCount = sizeof(ELEM_TYPE)*8;
  ELEM_TYPE minVal = ELEM_TYPE(1) << (bitCount - 1);
  ELEM_TYPE maxVal = minVal - 1;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1 - e2;
          bool sameSign = (e1 < 0) == (e2 >= 0);
          if (sameSign and ((e1 < 0) != (dest < 0)))
            {
              if (e1 < 0)
                dest = minVal;
              else
                dest = maxVal;
              csRegs_.write(CsrNumber::VXSAT, PrivilegeMode::Machine, 1);
            }
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVssub_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vssub_vv<int8_t> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half:   vssub_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word:   vssub_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2:  vssub_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4:  vssub_vv<Int128> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8:  vssub_vv<Int256> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vssub_vv<Int512> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vssub_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vssub_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, dest = 0;

  uint32_t bitCount = sizeof(ELEM_TYPE)*8;
  ELEM_TYPE minVal = ELEM_TYPE(1) << (bitCount - 1);
  ELEM_TYPE maxVal = minVal - 1;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 - e2;
          bool sameSign = (e1 < 0) == (e2 >= 0);
          if (sameSign and ((e1 < 0) != (dest < 0)))
            {
              if (e1 < 0)
                dest = minVal;
              else
                dest = maxVal;
              csRegs_.write(CsrNumber::VXSAT, PrivilegeMode::Machine, 1);
            }
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVssub_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vssub_vx<int8_t> (vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Half:   vssub_vx<int16_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word:   vssub_vx<int32_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word2:  vssub_vx<int64_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word4:  vssub_vx<Int128> (vd, vs1, Int128(e2),  group, start, elems, masked); break;
    case EW::Word8:  vssub_vx<Int256> (vd, vs1, Int256(e2),  group, start, elems, masked); break;
    case EW::Word16: vssub_vx<Int512> (vd, vs1, Int512(e2),  group, start, elems, masked); break;
    case EW::Word32: vssub_vx<Int1024>(vd, vs1, Int1024(e2), group, start, elems, masked); break;
    }
}


template <typename T>
static
void
roundoff(VecRoundingMode mode, T& value, unsigned d)
{
  if (d == 0)
    return;

  unsigned bit = 0;

  unsigned vd = unsigned((value >> d) & 1);
  unsigned vd_1 = unsigned((value >> (d-1)) & 1);

  switch(mode)
    {
    case VecRoundingMode::NearestUp:
      bit = vd_1;
      break;

    case VecRoundingMode::NearestEven:
      bit = vd_1 & ( (((T(1) << (d-1)) & value) != 0)  |  vd );
      break;

    case VecRoundingMode::Down:
      break;

    case VecRoundingMode::Odd:
      bit = (~vd & 1)  & ( ((T(1) << d) & value) != 0 );
      break;

    default:
      break;
    }

  T extra = bit;
  value = value + (extra << d);
  value >>= d;
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vaadd_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0;

  typedef typename makeDoubleWide<ELEM_TYPE>::type ELEM_TYPE2; // Double wide

  URV rmVal = 0;
  peekCsr(CsrNumber::VXRM, rmVal);
  VecRoundingMode rm = VecRoundingMode(rmVal);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          ELEM_TYPE2 temp = e1;
          temp += e2;
          roundoff(rm, temp, 1);
          ELEM_TYPE dest = ELEM_TYPE(temp);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVaadd_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vaadd_vv<int8_t> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half:   vaadd_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word:   vaadd_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2:  vaadd_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4:  vaadd_vv<Int128> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8:  vaadd_vv<Int256> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vaadd_vv<Int512> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: assert(0 && "1024-bit fixed point not yet implemented"); break;
    }
}


template <typename URV>
void
Hart<URV>::execVaaddu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vaadd_vv<uint8_t> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half:   vaadd_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word:   vaadd_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2:  vaadd_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4:  vaadd_vv<Uint128> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8:  vaadd_vv<Uint256> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vaadd_vv<Uint512> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: assert(0 && "1024-bit fixed point not yet implemented"); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vaadd_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0;

  typedef typename makeDoubleWide<ELEM_TYPE>::type ELEM_TYPE2; // Double wide

  URV rmVal = 0;
  peekCsr(CsrNumber::VXRM, rmVal);
  VecRoundingMode rm = VecRoundingMode(rmVal);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          ELEM_TYPE2 temp = e1;
          temp += e2;
          roundoff(rm, temp, 1);
          ELEM_TYPE dest = ELEM_TYPE(temp);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVaadd_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vaadd_vx<int8_t> (vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Half:   vaadd_vx<int16_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word:   vaadd_vx<int32_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word2:  vaadd_vx<int64_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word4:  vaadd_vx<Int128> (vd, vs1, Int128(e2),  group, start, elems, masked); break;
    case EW::Word8:  vaadd_vx<Int256> (vd, vs1, Int256(e2),  group, start, elems, masked); break;
    case EW::Word16: vaadd_vx<Int512> (vd, vs1, Int512(e2),  group, start, elems, masked); break;
    case EW::Word32: assert(0 && "1024-bit fixed point not yet implemented"); break;
    }
}


template <typename URV>
void
Hart<URV>::execVaaddu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vaadd_vx<uint8_t> (vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Half:   vaadd_vx<uint16_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word:   vaadd_vx<uint32_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word2:  vaadd_vx<uint64_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word4:  vaadd_vx<Uint128> (vd, vs1, Int128(e2),  group, start, elems, masked); break;
    case EW::Word8:  vaadd_vx<Uint256> (vd, vs1, Int256(e2),  group, start, elems, masked); break;
    case EW::Word16: vaadd_vx<Uint512> (vd, vs1, Int512(e2),  group, start, elems, masked); break;
    case EW::Word32: assert(0 && "1024-bit fixed point not yet implemented"); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vasub_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0;

  typedef typename makeDoubleWide<ELEM_TYPE>::type ELEM_TYPE2; // Double wide

  URV rmVal = 0;
  peekCsr(CsrNumber::VXRM, rmVal);
  VecRoundingMode rm = VecRoundingMode(rmVal);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          ELEM_TYPE2 temp = e1;
          temp -= e2;
          roundoff(rm, temp, 1);
          ELEM_TYPE dest = ELEM_TYPE(temp);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVasub_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vasub_vv<int8_t> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half:   vasub_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word:   vasub_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2:  vasub_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4:  vasub_vv<Int128> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8:  vasub_vv<Int256> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vasub_vv<Int512> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: assert(0 && "1024-bit fixed point not yet implemented"); break;
    }
}


template <typename URV>
void
Hart<URV>::execVasubu_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vasub_vv<uint8_t> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half:   vasub_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word:   vasub_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2:  vasub_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4:  vasub_vv<Uint128> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8:  vasub_vv<Uint256> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vasub_vv<Uint512> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: assert(0 && "1024-bit fixed point not yet implemented"); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vasub_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0;

  typedef typename makeDoubleWide<ELEM_TYPE>::type ELEM_TYPE2; // Double wide

  URV rmVal = 0;
  peekCsr(CsrNumber::VXRM, rmVal);
  VecRoundingMode rm = VecRoundingMode(rmVal);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          ELEM_TYPE2 temp = e1;
          temp -= e2;
          roundoff(rm, temp, 1);
          ELEM_TYPE dest = ELEM_TYPE(temp);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVasub_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vasub_vx<int8_t> (vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Half:   vasub_vx<int16_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word:   vasub_vx<int32_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word2:  vasub_vx<int64_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word4:  vasub_vx<Int128> (vd, vs1, Int128(e2),  group, start, elems, masked); break;
    case EW::Word8:  vasub_vx<Int256> (vd, vs1, Int256(e2),  group, start, elems, masked); break;
    case EW::Word16: vasub_vx<Int512> (vd, vs1, Int512(e2),  group, start, elems, masked); break;
    case EW::Word32: assert(0 && "1024-bit fixed point not yet implemented"); break;
    }
}


template <typename URV>
void
Hart<URV>::execVasubu_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vasub_vx<uint8_t> (vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Half:   vasub_vx<uint16_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word:   vasub_vx<uint32_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word2:  vasub_vx<uint64_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word4:  vasub_vx<Uint128> (vd, vs1, Int128(e2),  group, start, elems, masked); break;
    case EW::Word8:  vasub_vx<Uint256> (vd, vs1, Int256(e2),  group, start, elems, masked); break;
    case EW::Word16: vasub_vx<Uint512> (vd, vs1, Int512(e2),  group, start, elems, masked); break;
    case EW::Word32: assert(0 && "1024-bit fixed point not yet implemented"); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vsmul_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0;

  typedef typename makeDoubleWide<ELEM_TYPE>::type ELEM_TYPE2; // Double wide

  URV rmVal = 0;
  peekCsr(CsrNumber::VXRM, rmVal);
  VecRoundingMode rm = VecRoundingMode(rmVal);

  int elemBits = integerWidth<ELEM_TYPE> ();
  ELEM_TYPE minVal = ELEM_TYPE(1) << (elemBits - 1);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          ELEM_TYPE2 dest = 0;
          if (e1 == minVal and e2 == minVal)
            {
              // Result saturates at max positive value.
              dest = minVal - 1;
              csRegs_.write(CsrNumber::VXSAT, PrivilegeMode::Machine, 1);
            }
          else
            {
              ELEM_TYPE2 temp = e1;
              temp *= e2;
              roundoff(rm, temp, sizeof(ELEM_TYPE)*8 - 1);
              dest = temp;
            }

          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVsmul_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vsmul_vv<int8_t> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half:   vsmul_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word:   vsmul_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2:  vsmul_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4:  vsmul_vv<Int128> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8:  vsmul_vv<Int256> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vsmul_vv<Int512> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: assert(0 && "1024-bit fixed point not yet implemented"); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vsmul_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0;

  typedef typename makeDoubleWide<ELEM_TYPE>::type ELEM_TYPE2; // Double wide

  URV rmVal = 0;
  peekCsr(CsrNumber::VXRM, rmVal);
  VecRoundingMode rm = VecRoundingMode(rmVal);

  int elemBits = integerWidth<ELEM_TYPE> ();
  ELEM_TYPE minVal = ELEM_TYPE(1) << (elemBits - 1);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          ELEM_TYPE2 dest = 0;
          if (e1 == minVal and e2 == minVal)
            {
              // Result saturates at max positive value.
              dest = minVal - 1;
              csRegs_.write(CsrNumber::VXSAT, PrivilegeMode::Machine, 1);
            }
          else
            {
              ELEM_TYPE2 temp = e1;
              temp *= e2;
              roundoff(rm, temp, sizeof(ELEM_TYPE)*8 - 1);
              dest = temp;
            }

          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVsmul_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vsmul_vx<int8_t> (vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Half:   vsmul_vx<int16_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word:   vsmul_vx<int32_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word2:  vsmul_vx<int64_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word4:  vsmul_vx<Int128> (vd, vs1, Int128(e2),  group, start, elems, masked); break;
    case EW::Word8:  vsmul_vx<Int256> (vd, vs1, Int256(e2),  group, start, elems, masked); break;
    case EW::Word16: vsmul_vx<Int512> (vd, vs1, Int512(e2),  group, start, elems, masked); break;
    case EW::Word32: assert(0 && "1024-bit fixed point not yet implemented"); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vssr_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0, e2 = 0;

  URV rmVal = 0;
  peekCsr(CsrNumber::VXRM, rmVal);
  VecRoundingMode rm = VecRoundingMode(rmVal);

  unsigned elemBits = integerWidth<ELEM_TYPE> ();
  unsigned mask = elemBits - 1;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          ELEM_TYPE dest = e1;
          unsigned amount = unsigned(e2) & mask;
          roundoff(rm, dest, amount);

          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVssrl_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vssr_vv<uint8_t> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half:   vssr_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word:   vssr_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2:  vssr_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4:  vssr_vv<Uint128> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8:  vssr_vv<Uint256> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vssr_vv<Uint512> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vssr_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vssr_vx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;
  ELEM_TYPE e1 = 0;

  URV rmVal = 0;
  peekCsr(CsrNumber::VXRM, rmVal);
  VecRoundingMode rm = VecRoundingMode(rmVal);

  unsigned elemBits = integerWidth<ELEM_TYPE> ();
  unsigned mask = elemBits - 1;
  unsigned amount = unsigned(e2) & mask;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          ELEM_TYPE dest = e1;
          roundoff(rm, dest, amount);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVssrl_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vssr_vx<uint8_t> (vd, vs1, e2,           group, start, elems, masked); break;
    case EW::Half:   vssr_vx<uint16_t>(vd, vs1, e2,           group, start, elems, masked); break;
    case EW::Word:   vssr_vx<uint32_t>(vd, vs1, e2,           group, start, elems, masked); break;
    case EW::Word2:  vssr_vx<uint64_t>(vd, vs1, e2,           group, start, elems, masked); break;
    case EW::Word4:  vssr_vx<Uint128> (vd, vs1, Uint128(e2),  group, start, elems, masked); break;
    case EW::Word8:  vssr_vx<Uint256> (vd, vs1, Uint256(e2),  group, start, elems, masked); break;
    case EW::Word16: vssr_vx<Uint512> (vd, vs1, Uint512(e2),  group, start, elems, masked); break;
    case EW::Word32: vssr_vx<Uint1024>(vd, vs1, Uint1024(e2), group, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVssrl_vi(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  imm = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vssr_vx<uint8_t> (vd, vs1, imm,           group, start, elems, masked); break;
    case EW::Half:   vssr_vx<uint16_t>(vd, vs1, imm,           group, start, elems, masked); break;
    case EW::Word:   vssr_vx<uint32_t>(vd, vs1, imm,           group, start, elems, masked); break;
    case EW::Word2:  vssr_vx<uint64_t>(vd, vs1, imm,           group, start, elems, masked); break;
    case EW::Word4:  vssr_vx<Uint128> (vd, vs1, Uint128(imm),  group, start, elems, masked); break;
    case EW::Word8:  vssr_vx<Uint256> (vd, vs1, Uint256(imm),  group, start, elems, masked); break;
    case EW::Word16: vssr_vx<Uint512> (vd, vs1, Uint512(imm),  group, start, elems, masked); break;
    case EW::Word32: vssr_vx<Uint1024>(vd, vs1, Uint1024(imm), group, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVssra_vv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vssr_vv<int8_t> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half:   vssr_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word:   vssr_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2:  vssr_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4:  vssr_vv<Int128> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8:  vssr_vv<Int256> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vssr_vv<Int512> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: vssr_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVssra_vx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vssr_vx<int8_t> (vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Half:   vssr_vx<int16_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word:   vssr_vx<int32_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word2:  vssr_vx<int64_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word4:  vssr_vx<Int128> (vd, vs1, Int128(e2),  group, start, elems, masked); break;
    case EW::Word8:  vssr_vx<Int256> (vd, vs1, Int256(e2),  group, start, elems, masked); break;
    case EW::Word16: vssr_vx<Int512> (vd, vs1, Int512(e2),  group, start, elems, masked); break;
    case EW::Word32: vssr_vx<Int1024>(vd, vs1, Int1024(e2), group, start, elems, masked); break;
    }
}


template <typename URV>
void
Hart<URV>::execVssra_vi(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  imm = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vssr_vx<int8_t> (vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Half:   vssr_vx<int16_t>(vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Word:   vssr_vx<int32_t>(vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Word2:  vssr_vx<int64_t>(vd, vs1, imm,          group, start, elems, masked); break;
    case EW::Word4:  vssr_vx<Int128> (vd, vs1, Int128(imm),  group, start, elems, masked); break;
    case EW::Word8:  vssr_vx<Int256> (vd, vs1, Int256(imm),  group, start, elems, masked); break;
    case EW::Word16: vssr_vx<Int512> (vd, vs1, Int512(imm),  group, start, elems, masked); break;
    case EW::Word32: vssr_vx<Int1024>(vd, vs1, Int1024(imm), group, start, elems, masked); break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vnclip_wv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  typedef typename std::make_unsigned<ELEM_TYPE>::type  U_ELEM_TYPE;
  typedef typename makeDoubleWide<ELEM_TYPE>::type ELEM_TYPE2X; // Double wide

  unsigned errors = 0;
  ELEM_TYPE2X e1 = 0;
  ELEM_TYPE e2 = 0;

  URV rmVal = 0;
  peekCsr(CsrNumber::VXRM, rmVal);
  VecRoundingMode rm = VecRoundingMode(rmVal);

  unsigned elemBits = integerWidth<ELEM_TYPE> ();
  unsigned mask = elemBits - 1;
  unsigned group2x = group*2;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group2x, e1) and vecRegs_.read(vs2, ix, group, e2))
        {
          unsigned amount = unsigned(e2) & mask;
          roundoff(rm, e1, amount);

          ELEM_TYPE dest = ELEM_TYPE(e1);
          if (e1 != ELEM_TYPE2X(dest))
            {
              if (std::is_same<ELEM_TYPE, U_ELEM_TYPE>::value)
                dest = maxVal<ELEM_TYPE>();
              else
                dest = (e1 < 0) ? minVal<ELEM_TYPE>() : maxVal<ELEM_TYPE>();
            }

          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVnclipu_wv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vnclip_wv<uint8_t> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half:   vnclip_wv<uint16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word:   vnclip_wv<uint32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2:  vnclip_wv<uint64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4:  vnclip_wv<Uint128> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8:  vnclip_wv<Uint256> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vnclip_wv<Uint512> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vnclip_wx(unsigned vd, unsigned vs1, ELEM_TYPE e2, unsigned group,
                     unsigned start, unsigned elems, bool masked)
{
  typedef typename std::make_unsigned<ELEM_TYPE>::type  U_ELEM_TYPE;
  typedef typename makeDoubleWide<ELEM_TYPE>::type ELEM_TYPE2X; // Double wide

  unsigned errors = 0;
  ELEM_TYPE2X e1 = 0;

  URV rmVal = 0;
  peekCsr(CsrNumber::VXRM, rmVal);
  VecRoundingMode rm = VecRoundingMode(rmVal);

  unsigned elemBits = integerWidth<ELEM_TYPE> ();
  unsigned mask = elemBits - 1;
  unsigned amount = unsigned(e2) & mask;
  unsigned group2x = group*2;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group2x, e1))
        {
          roundoff(rm, e1, amount);

          ELEM_TYPE dest = ELEM_TYPE(e1);
          if (e1 != ELEM_TYPE2X(dest))
            {
              if (std::is_same<ELEM_TYPE, U_ELEM_TYPE>::value)
                dest = maxVal<ELEM_TYPE>();
              else
                dest = (e1 < 0) ? minVal<ELEM_TYPE>() : maxVal<ELEM_TYPE>();
            }

          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVnclipu_wx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vnclip_wx<uint8_t> (vd, vs1, e2,           group, start, elems, masked); break;
    case EW::Half:   vnclip_wx<uint16_t>(vd, vs1, e2,           group, start, elems, masked); break;
    case EW::Word:   vnclip_wx<uint32_t>(vd, vs1, e2,           group, start, elems, masked); break;
    case EW::Word2:  vnclip_wx<uint64_t>(vd, vs1, e2,           group, start, elems, masked); break;
    case EW::Word4:  vnclip_wx<Uint128> (vd, vs1, Uint128(e2),  group, start, elems, masked); break;
    case EW::Word8:  vnclip_wx<Uint256> (vd, vs1, Uint256(e2),  group, start, elems, masked); break;
    case EW::Word16: vnclip_wx<Uint512> (vd, vs1, Uint512(e2),  group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVnclipu_wi(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  imm = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vnclip_wx<uint8_t> (vd, vs1, imm,           group, start, elems, masked); break;
    case EW::Half:   vnclip_wx<uint16_t>(vd, vs1, imm,           group, start, elems, masked); break;
    case EW::Word:   vnclip_wx<uint32_t>(vd, vs1, imm,           group, start, elems, masked); break;
    case EW::Word2:  vnclip_wx<uint64_t>(vd, vs1, imm,           group, start, elems, masked); break;
    case EW::Word4:  vnclip_wx<Uint128> (vd, vs1, Uint128(imm),  group, start, elems, masked); break;
    case EW::Word8:  vnclip_wx<Uint256> (vd, vs1, Uint256(imm),  group, start, elems, masked); break;
    case EW::Word16: vnclip_wx<Uint512> (vd, vs1, Uint512(imm),  group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVnclip_wv(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vnclip_wv<int8_t> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Half:   vnclip_wv<int16_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word:   vnclip_wv<int32_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word2:  vnclip_wv<int64_t>(vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word4:  vnclip_wv<Int128> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word8:  vnclip_wv<Int256> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word16: vnclip_wv<Int512> (vd, vs1, vs2, group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVnclip_wx(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vnclip_wx<int8_t> (vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Half:   vnclip_wx<int16_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word:   vnclip_wx<int32_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word2:  vnclip_wx<int64_t>(vd, vs1, e2,          group, start, elems, masked); break;
    case EW::Word4:  vnclip_wx<Int128> (vd, vs1, Int128(e2),  group, start, elems, masked); break;
    case EW::Word8:  vnclip_wx<Int256> (vd, vs1, Int256(e2),  group, start, elems, masked); break;
    case EW::Word16: vnclip_wx<Int512> (vd, vs1, Int512(e2),  group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
void
Hart<URV>::execVnclip_wi(const DecodedInst* di)
{
  if (not checkMaskableInst(di))
    return;

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  imm = di->op2();
  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  typedef ElementWidth EW;
  switch (sew)
    {
    case EW::Byte:   vnclip_wx<int8_t> (vd, vs1, imm,           group, start, elems, masked); break;
    case EW::Half:   vnclip_wx<int16_t>(vd, vs1, imm,           group, start, elems, masked); break;
    case EW::Word:   vnclip_wx<int32_t>(vd, vs1, imm,           group, start, elems, masked); break;
    case EW::Word2:  vnclip_wx<int64_t>(vd, vs1, imm,           group, start, elems, masked); break;
    case EW::Word4:  vnclip_wx<Int128> (vd, vs1, Int128(imm),  group, start, elems, masked); break;
    case EW::Word8:  vnclip_wx<Int256> (vd, vs1, Int256(imm),  group, start, elems, masked); break;
    case EW::Word16: vnclip_wx<Int512> (vd, vs1, Int512(imm),  group, start, elems, masked); break;
    case EW::Word32: illegalInst(di); break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vectorLoad(const DecodedInst* di, ElementWidth eew, bool faultFirst)
{
  // Compute emul: lmul*eew/sew
  unsigned groupX8 = vecRegs_.groupMultiplierX8();
  groupX8 = groupX8 * vecRegs_.elementWidthInBits(eew) / vecRegs_.elementWidthInBits();
  GroupMultiplier lmul = GroupMultiplier::One;
  bool badConfig = false;
  if (not vecRegs_.groupNumberX8ToSymbol(groupX8, lmul))
    badConfig = true;
  else
    badConfig = not vecRegs_.legalConfig(eew, lmul);

  if (not isVecLegal() or badConfig)
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(), rs1 = di->op1(), errors = 0;
  uint64_t addr = intRegs_.read(rs1);

  unsigned start = vecRegs_.startIndex();
  unsigned elemCount = vecRegs_.elemCount();

  // FIX TODO: check permissions, translate, ....
  for (unsigned ix = start; ix < elemCount; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      auto cause = ExceptionCause::NONE;
      auto secCause = SecondaryCause::NONE;

      ELEM_TYPE elem = 0;
      if constexpr (sizeof(elem) > 8)
        {
          for (unsigned n = 0; n < sizeof(elem); n += 8)
            {
              uint64_t dword = 0;
              cause = determineLoadException(rs1, addr, addr, 8, secCause);
              if (cause != ExceptionCause::NONE)
                break;
              memory_.read(addr + n, dword);
              elem <<= 64;
              elem |= dword;
            }
        }
      else
        {
          if (determineLoadException(rs1, addr, addr, sizeof(elem), secCause) !=
              ExceptionCause::NONE)
            memory_.read(addr, elem);
        }

      if (cause != ExceptionCause::NONE)
        {
          vecRegs_.setStartIndex(ix);
          csRegs_.write(CsrNumber::VSTART, PrivilegeMode::Machine, ix);
          if (ix == 0 or not faultFirst)
            initiateLoadException(cause, addr, secCause);
          break;
        }

      if (not vecRegs_.write(vd, ix, groupX8, elem))
        {
          errors++;
          break;
        }

      addr += sizeof(ELEM_TYPE);
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVle8_v(const DecodedInst* di)
{
  vectorLoad<uint8_t>(di, ElementWidth::Byte, false);
}


template <typename URV>
void
Hart<URV>::execVle16_v(const DecodedInst* di)
{
  vectorLoad<uint16_t>(di, ElementWidth::Half, false);
}


template <typename URV>
void
Hart<URV>::execVle32_v(const DecodedInst* di)
{
  vectorLoad<uint32_t>(di, ElementWidth::Word, false);
}


template <typename URV>
void
Hart<URV>::execVle64_v(const DecodedInst* di)
{
  vectorLoad<uint64_t>(di, ElementWidth::Word2, false);
}


template <typename URV>
void
Hart<URV>::execVle128_v(const DecodedInst* di)
{
  vectorLoad<Uint128>(di, ElementWidth::Word4, false);
}


template <typename URV>
void
Hart<URV>::execVle256_v(const DecodedInst* di)
{
  vectorLoad<Uint256>(di, ElementWidth::Word8, false);
}


template <typename URV>
void
Hart<URV>::execVle512_v(const DecodedInst* di)
{
  vectorLoad<Uint512>(di, ElementWidth::Word16, false);
}


template <typename URV>
void
Hart<URV>::execVle1024_v(const DecodedInst* di)
{
  vectorLoad<Uint1024>(di, ElementWidth::Word32, false);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vectorStore(const DecodedInst* di, ElementWidth eew)
{
  // Compute emul: lmul*eew/sew
  unsigned groupX8 = vecRegs_.groupMultiplierX8();
  groupX8 = groupX8 * vecRegs_.elementWidthInBits(eew) / vecRegs_.elementWidthInBits();
  GroupMultiplier lmul = GroupMultiplier::One;
  bool badConfig = false;
  if (not vecRegs_.groupNumberX8ToSymbol(groupX8, lmul))
    badConfig = true;
  else
    badConfig = not vecRegs_.legalConfig(eew, lmul);

  if (not isVecLegal() or badConfig)
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  uint32_t vd = di->op0(), rs1 = di->op1(), errors = 0;
  uint64_t addr = intRegs_.read(rs1);

  unsigned start = vecRegs_.startIndex();
  unsigned elemCount = vecRegs_.elemCount();

  // TODO check permissions, translate, ....
  for (unsigned ix = start; ix < elemCount; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;
      ELEM_TYPE elem = 0;
      if (not vecRegs_.read(vd, ix, groupX8, elem))
        {
          errors++;
          break;
        }

      auto cause = ExceptionCause::NONE;
      auto secCause = SecondaryCause::NONE;

      if constexpr (sizeof(elem) > 8)
        {
          for (unsigned n = 0; n < sizeof(elem); n += 8)
            {
              uint64_t dword = uint64_t(elem);
              bool forced = false;
              cause = determineStoreException(rs1, URV(addr), addr, dword, secCause, forced);
              if (cause != ExceptionCause::NONE)
                break;

              memory_.write(hartIx_, addr + n, dword);
              elem >>= 64;
            }
        }
      else
        {
          bool forced = false;
          if (determineStoreException(rs1, URV(addr), addr, elem, secCause, forced) !=
              ExceptionCause::NONE)
            memory_.write(hartIx_, addr, elem);
        }

      if (cause != ExceptionCause::NONE)
        {
          vecRegs_.setStartIndex(ix);
          csRegs_.write(CsrNumber::VSTART, PrivilegeMode::Machine, ix);
          initiateStoreException(cause, addr, secCause);
          break;
        }

      addr += sizeof(ELEM_TYPE);
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVse8_v(const DecodedInst* di)
{
  vectorStore<uint8_t>(di, ElementWidth::Byte);
}


template <typename URV>
void
Hart<URV>::execVse16_v(const DecodedInst* di)
{
  vectorStore<uint16_t>(di, ElementWidth::Half);
}


template <typename URV>
void
Hart<URV>::execVse32_v(const DecodedInst* di)
{
  vectorStore<uint32_t>(di, ElementWidth::Word);
}


template <typename URV>
void
Hart<URV>::execVse64_v(const DecodedInst* di)
{
  vectorStore<uint64_t>(di, ElementWidth::Word2);
}


template <typename URV>
void
Hart<URV>::execVse128_v(const DecodedInst* di)
{
  vectorStore<Uint128>(di, ElementWidth::Word4);
}


template <typename URV>
void
Hart<URV>::execVse256_v(const DecodedInst* di)
{
  vectorStore<Uint256>(di, ElementWidth::Word8);
}


template <typename URV>
void
Hart<URV>::execVse512_v(const DecodedInst* di)
{
  vectorStore<Uint512>(di, ElementWidth::Word16);
}


template <typename URV>
void
Hart<URV>::execVse1024_v(const DecodedInst* di)
{
  vectorStore<Uint1024>(di, ElementWidth::Word32);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vectorLoadWholeReg(const DecodedInst* di, ElementWidth eew)
{
  unsigned groupCode = di->op3();
  bool badConfig = groupCode > 3;
  GroupMultiplier gm = GroupMultiplier(groupCode);

  badConfig = badConfig or not vecRegs_.legalConfig(eew, gm);
  if ((not isVecLegal()) or badConfig)
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(), rs1 = di->op1(), errors = 0;
  URV addr = intRegs_.read(rs1);

  unsigned groupX8 = vecRegs_.groupMultiplierX8(gm), start = 0;
  unsigned elemCount = (groupX8*vecRegs_.bytesPerRegister()) / 8;

  // TODO check permissions, translate, ....
  for (unsigned ix = start; ix < elemCount; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;
      bool exception = false;
      ELEM_TYPE elem = 0;
      if constexpr (sizeof(elem) > 8)
        {
          for (unsigned n = 0; n < sizeof(elem) and not exception; n += 8)
            {
              uint64_t dword = 0;
              memory_.read(addr + n, dword);
              elem <<= 64;
              elem |= dword;
            }
        }
      else
        memory_.read(addr, elem);

      if (exception)
        {
          vecRegs_.setStartIndex(ix);
          csRegs_.write(CsrNumber::VSTART, PrivilegeMode::Machine, ix);
          break;
        }

      if (not vecRegs_.write(vd, ix, groupX8, elem))
        {
          errors++;
          break;
        }

      addr += sizeof(ELEM_TYPE);
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVlre8_v(const DecodedInst* di)
{
  vectorLoadWholeReg<uint8_t>(di, ElementWidth::Byte);
}


template <typename URV>
void
Hart<URV>::execVlre16_v(const DecodedInst* di)
{
  vectorLoadWholeReg<uint16_t>(di, ElementWidth::Half);
}


template <typename URV>
void
Hart<URV>::execVlre32_v(const DecodedInst* di)
{
  vectorLoadWholeReg<uint32_t>(di, ElementWidth::Word);
}


template <typename URV>
void
Hart<URV>::execVlre64_v(const DecodedInst* di)
{
  vectorLoadWholeReg<uint64_t>(di, ElementWidth::Word2);
}


template <typename URV>
void
Hart<URV>::execVlre128_v(const DecodedInst* di)
{
  vectorLoadWholeReg<Uint128>(di, ElementWidth::Word4);
}


template <typename URV>
void
Hart<URV>::execVlre256_v(const DecodedInst* di)
{
  vectorLoadWholeReg<Uint256>(di, ElementWidth::Word8);
}


template <typename URV>
void
Hart<URV>::execVlre512_v(const DecodedInst* di)
{
  vectorLoadWholeReg<Uint512>(di, ElementWidth::Word16);
}


template <typename URV>
void
Hart<URV>::execVlre1024_v(const DecodedInst* di)
{
  vectorLoadWholeReg<Uint1024>(di, ElementWidth::Word32);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vectorStoreWholeReg(const DecodedInst* di, ElementWidth eew)
{
  unsigned groupCode = di->op3();
  bool badConfig = groupCode > 3;
  GroupMultiplier gm = GroupMultiplier(groupCode);

  badConfig = badConfig or not vecRegs_.legalConfig(eew, gm);
  if ((not isVecLegal()) or badConfig)
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(), rs1 = di->op1(), errors = 0;
  URV addr = intRegs_.read(rs1);

  unsigned groupX8 = vecRegs_.groupMultiplierX8(gm), start = 0;
  unsigned elemCount = (groupX8*vecRegs_.bytesPerRegister()) / 8;

  // TODO check permissions, translate, ....
  for (unsigned ix = start; ix < elemCount; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;
      bool exception = false;
      ELEM_TYPE elem = 0;
      if (not vecRegs_.read(vd, ix, groupX8, elem))
        {
          errors++;
          break;
        }

      if constexpr (sizeof(elem) > 8)
        {
          for (unsigned n = 0; n < sizeof(elem) and not exception; n += 8)
            {
              uint64_t dword = uint64_t(elem);
              memory_.write(hartIx_, addr + n, dword);
              elem >>= 64;
            }
        }
      else
        memory_.write(hartIx_, addr, elem);

      if (exception)
        {
          vecRegs_.setStartIndex(ix);
          csRegs_.write(CsrNumber::VSTART, PrivilegeMode::Machine, ix);
          break;
        }

      addr += sizeof(ELEM_TYPE);
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVsre8_v(const DecodedInst* di)
{
  vectorStoreWholeReg<uint8_t>(di, ElementWidth::Byte);
}


template <typename URV>
void
Hart<URV>::execVsre16_v(const DecodedInst* di)
{
  vectorStoreWholeReg<uint16_t>(di, ElementWidth::Half);
}


template <typename URV>
void
Hart<URV>::execVsre32_v(const DecodedInst* di)
{
  vectorStoreWholeReg<uint32_t>(di, ElementWidth::Word);
}


template <typename URV>
void
Hart<URV>::execVsre64_v(const DecodedInst* di)
{
  vectorStoreWholeReg<uint64_t>(di, ElementWidth::Word2);
}


template <typename URV>
void
Hart<URV>::execVsre128_v(const DecodedInst* di)
{
  vectorStoreWholeReg<Uint128>(di, ElementWidth::Word4);
}


template <typename URV>
void
Hart<URV>::execVsre256_v(const DecodedInst* di)
{
  vectorStoreWholeReg<Uint256>(di, ElementWidth::Word8);
}


template <typename URV>
void
Hart<URV>::execVsre512_v(const DecodedInst* di)
{
  vectorStoreWholeReg<Uint512>(di, ElementWidth::Word16);
}


template <typename URV>
void
Hart<URV>::execVsre1024_v(const DecodedInst* di)
{
  vectorStoreWholeReg<Uint1024>(di, ElementWidth::Word32);
}


template <typename URV>
void
Hart<URV>::execVle8ff_v(const DecodedInst* di)
{
  vectorLoad<uint8_t>(di, ElementWidth::Byte, true);
}


template <typename URV>
void
Hart<URV>::execVle16ff_v(const DecodedInst* di)
{
  vectorLoad<uint16_t>(di, ElementWidth::Half, true);
}


template <typename URV>
void
Hart<URV>::execVle32ff_v(const DecodedInst* di)
{
  vectorLoad<uint32_t>(di, ElementWidth::Word, true);
}


template <typename URV>
void
Hart<URV>::execVle64ff_v(const DecodedInst* di)
{
  vectorLoad<uint64_t>(di, ElementWidth::Word2, true);
}


template <typename URV>
void
Hart<URV>::execVle128ff_v(const DecodedInst* di)
{
  vectorLoad<Uint128>(di, ElementWidth::Word4, true);
}


template <typename URV>
void
Hart<URV>::execVle256ff_v(const DecodedInst* di)
{
  vectorLoad<Uint256>(di, ElementWidth::Word8, true);
}


template <typename URV>
void
Hart<URV>::execVle512ff_v(const DecodedInst* di)
{
  vectorLoad<Uint512>(di, ElementWidth::Word16, true);
}


template <typename URV>
void
Hart<URV>::execVle1024ff_v(const DecodedInst* di)
{
  vectorLoad<Uint1024>(di, ElementWidth::Word32, true);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vectorLoadStrided(const DecodedInst* di, ElementWidth eew)
{
  // Compute emul: lmul*eew/sew
  unsigned groupX8 = vecRegs_.groupMultiplierX8();
  groupX8 = groupX8 * vecRegs_.elementWidthInBits(eew) / vecRegs_.elementWidthInBits();
  GroupMultiplier lmul = GroupMultiplier::One;
  bool badConfig = false;
  if (not vecRegs_.groupNumberX8ToSymbol(groupX8, lmul))
    badConfig = true;
  else
    badConfig = not vecRegs_.legalConfig(eew, lmul);

  if (not isVecLegal() or badConfig)
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(), rs1 = di->op1(), rs2 = di->op2(), errors = 0;
  uint64_t addr = intRegs_.read(rs1);
  uint64_t stride = intRegs_.read(rs2);

  unsigned start = vecRegs_.startIndex();
  unsigned elemCount = vecRegs_.elemCount();

  // TODO check permissions, translate, ....
  for (unsigned ix = start; ix < elemCount; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      auto cause = ExceptionCause::NONE;
      auto secCause = SecondaryCause::NONE;

      ELEM_TYPE elem = 0;
      if constexpr (sizeof(elem) > 8)
        {
          for (unsigned n = 0; n < sizeof(elem); n += 8)
            {
              uint64_t dword = 0;
              cause = determineLoadException(rs1, addr, addr, 8, secCause);
              if (cause != ExceptionCause::NONE)
                break;
              memory_.read(addr + n, dword);
              elem <<= 64;
              elem |= dword;
            }
        }
      else
        {
          if (determineLoadException(rs1, addr, addr, sizeof(elem), secCause) !=
              ExceptionCause::NONE)
            memory_.read(addr, elem);
        }

      if (cause != ExceptionCause::NONE)
        {
          vecRegs_.setStartIndex(ix);
          csRegs_.write(CsrNumber::VSTART, PrivilegeMode::Machine, ix);
          initiateLoadException(cause, addr, secCause);
          break;
        }

      if (not vecRegs_.write(vd, ix, groupX8, elem))
        {
          errors++;
          break;
        }

      addr += stride;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVlse8_v(const DecodedInst* di)
{
  vectorLoadStrided<uint8_t>(di, ElementWidth::Byte);
}


template <typename URV>
void
Hart<URV>::execVlse16_v(const DecodedInst* di)
{
  vectorLoadStrided<uint16_t>(di, ElementWidth::Half);
}


template <typename URV>
void
Hart<URV>::execVlse32_v(const DecodedInst* di)
{
  vectorLoadStrided<uint32_t>(di, ElementWidth::Word);
}


template <typename URV>
void
Hart<URV>::execVlse64_v(const DecodedInst* di)
{
  vectorLoadStrided<uint64_t>(di, ElementWidth::Word2);
}


template <typename URV>
void
Hart<URV>::execVlse128_v(const DecodedInst* di)
{
  vectorLoadStrided<Uint128>(di, ElementWidth::Word4);
}


template <typename URV>
void
Hart<URV>::execVlse256_v(const DecodedInst* di)
{
  vectorLoadStrided<Uint256>(di, ElementWidth::Word8);
}


template <typename URV>
void
Hart<URV>::execVlse512_v(const DecodedInst* di)
{
  vectorLoadStrided<Uint512>(di, ElementWidth::Word16);
}


template <typename URV>
void
Hart<URV>::execVlse1024_v(const DecodedInst* di)
{
  vectorLoadStrided<Uint1024>(di, ElementWidth::Word32);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vectorStoreStrided(const DecodedInst* di, ElementWidth eew)
{
  // Compute emul: lmul*eew/sew
  unsigned groupX8 = vecRegs_.groupMultiplierX8();
  groupX8 = groupX8 * vecRegs_.elementWidthInBits(eew) / vecRegs_.elementWidthInBits();
  GroupMultiplier lmul = GroupMultiplier::One;
  bool badConfig = false;
  if (not vecRegs_.groupNumberX8ToSymbol(groupX8, lmul))
    badConfig = true;
  else
    badConfig = not vecRegs_.legalConfig(eew, lmul);

  if (not isVecLegal() or badConfig)
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(), rs1 = di->op1(), rs2 = di->op2(), errors = 0;
  uint64_t addr = intRegs_.read(rs1);
  uint64_t stride = intRegs_.read(rs2);

  unsigned start = vecRegs_.startIndex();
  unsigned elemCount = vecRegs_.elemCount();

  // TODO check permissions, translate, ....
  for (unsigned ix = start; ix < elemCount; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;
      bool exception = false;
      ELEM_TYPE elem = 0;
      if (not vecRegs_.read(vd, ix, groupX8, elem))
        {
          errors++;
          break;
        }

      if constexpr (sizeof(elem) > 8)
        {
          for (unsigned n = 0; n < sizeof(elem) and not exception; n += 8)
            {
              uint64_t dword = uint64_t(elem);
              memory_.write(hartIx_, addr + n, dword);
              elem >>= 64;
            }
        }
      else
        memory_.write(hartIx_, addr, elem);

      if (exception)
        {
          vecRegs_.setStartIndex(ix);
          csRegs_.write(CsrNumber::VSTART, PrivilegeMode::Machine, ix);
          break;
        }

      addr += stride;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVsse8_v(const DecodedInst* di)
{
  vectorStoreStrided<uint8_t>(di, ElementWidth::Byte);
}


template <typename URV>
void
Hart<URV>::execVsse16_v(const DecodedInst* di)
{
  vectorStoreStrided<uint16_t>(di, ElementWidth::Half);
}


template <typename URV>
void
Hart<URV>::execVsse32_v(const DecodedInst* di)
{
  vectorStoreStrided<uint32_t>(di, ElementWidth::Word);
}


template <typename URV>
void
Hart<URV>::execVsse64_v(const DecodedInst* di)
{
  vectorStoreStrided<uint64_t>(di, ElementWidth::Word2);
}


template <typename URV>
void
Hart<URV>::execVsse128_v(const DecodedInst* di)
{
  vectorStoreStrided<Uint128>(di, ElementWidth::Word4);
}


template <typename URV>
void
Hart<URV>::execVsse256_v(const DecodedInst* di)
{
  vectorStoreStrided<Uint256>(di, ElementWidth::Word8);
}


template <typename URV>
void
Hart<URV>::execVsse512_v(const DecodedInst* di)
{
  vectorStoreStrided<Uint512>(di, ElementWidth::Word16);
}


template <typename URV>
void
Hart<URV>::execVsse1024_v(const DecodedInst* di)
{
  vectorStoreStrided<Uint1024>(di, ElementWidth::Word32);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vectorLoadIndexed(const DecodedInst* di, ElementWidth offsetEew)
{
  if (not isVecLegal())
    {
      illegalInst(di);
      return;
    }

  uint32_t elemWidth = vecRegs_.elementWidthInBits();
  uint32_t offsetElemWidth = vecRegs_.elementWidthInBits(offsetEew);

  uint32_t groupX8 = vecRegs_.groupMultiplierX8();
  uint32_t offsetGroupX8 = (offsetElemWidth*groupX8)/elemWidth;

  GroupMultiplier offsetGroup{GroupMultiplier::One};
  bool badConfig = vecRegs_.groupNumberX8ToSymbol(offsetGroupX8, offsetGroup);
  if (not badConfig)
    badConfig = not vecRegs_.legalConfig(offsetEew, offsetGroup);
  if (badConfig)
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(), rs1 = di->op1(), vi = di->op2(), errors = 0;
  uint64_t addr = intRegs_.read(rs1);

  unsigned start = vecRegs_.startIndex();
  unsigned elemCount = vecRegs_.elemCount();

  // TODO check permissions, translate, ....
  for (unsigned ix = start; ix < elemCount; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      auto cause = ExceptionCause::NONE;
      auto secCause = SecondaryCause::NONE;

      uint64_t offset = 0;
      if (not vecRegs_.readIndex(vi, ix, offsetEew, offsetGroupX8, offset))
        {
          errors++;
          break;
        }

      ELEM_TYPE elem = 0;
      uint64_t eaddr = addr + offset;

      if constexpr (sizeof(elem) > 8)
        {
          for (unsigned n = 0; n < sizeof(elem); n += 8)
            {
              uint64_t dword = 0;
              cause = determineLoadException(rs1, eaddr, addr, 8, secCause);
              if (cause != ExceptionCause::NONE)
                break;
              memory_.read(eaddr + n, dword);
              elem <<= 64;
              elem |= dword;
            }
        }
      else
        {
          if (determineLoadException(rs1, eaddr, eaddr, sizeof(elem), secCause) !=
              ExceptionCause::NONE)
            memory_.read(eaddr, elem);
        }

      if (cause != ExceptionCause::NONE)
        {
          vecRegs_.setStartIndex(ix);
          csRegs_.write(CsrNumber::VSTART, PrivilegeMode::Machine, ix);
          initiateLoadException(cause, eaddr, secCause);
          break;
        }

      if (not vecRegs_.write(vd, ix, groupX8, elem))
        {
          errors++;
          break;
        }
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVlxei8_v(const DecodedInst* di)
{
  vectorLoadIndexed<uint8_t>(di, ElementWidth::Byte);
}


template <typename URV>
void
Hart<URV>::execVlxei16_v(const DecodedInst* di)
{
  vectorLoadIndexed<uint16_t>(di, ElementWidth::Half);
}


template <typename URV>
void
Hart<URV>::execVlxei32_v(const DecodedInst* di)
{
  vectorLoadIndexed<uint32_t>(di, ElementWidth::Word);
}


template <typename URV>
void
Hart<URV>::execVlxei64_v(const DecodedInst* di)
{
  vectorLoadIndexed<uint64_t>(di, ElementWidth::Word2);
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vectorStoreIndexed(const DecodedInst* di, ElementWidth offsetEew)
{
  if (not isVecLegal())
    {
      illegalInst(di);
      return;
    }

  uint32_t elemWidth = vecRegs_.elementWidthInBits();
  uint32_t offsetElemWidth = vecRegs_.elementWidthInBits(offsetEew);

  uint32_t groupX8 = vecRegs_.groupMultiplierX8();
  uint32_t offsetGroupX8 = (offsetElemWidth*groupX8)/elemWidth;

  GroupMultiplier offsetGroup{GroupMultiplier::One};
  bool badConfig = vecRegs_.groupNumberX8ToSymbol(offsetGroupX8, offsetGroup);
  if (not badConfig)
    badConfig = not vecRegs_.legalConfig(offsetEew, offsetGroup);
  if (badConfig)
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  uint32_t vd = di->op0(), rs1 = di->op1(), vi = di->op2(), errors = 0;
  uint64_t addr = intRegs_.read(rs1);

  unsigned start = vecRegs_.startIndex();
  unsigned elemCount = vecRegs_.elemCount();

  // TODO check permissions, translate, ....
  for (unsigned ix = start; ix < elemCount; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;
      ELEM_TYPE elem = 0;
      if (not vecRegs_.read(vd, ix, groupX8, elem))
        {
          errors++;
          break;
        }

      uint64_t offset = 0;
      if (not vecRegs_.readIndex(vi, ix, offsetEew, offsetGroupX8, offset))
        {
          errors++;
          break;
        }

      uint64_t eaddr = addr + offset;
      auto cause = ExceptionCause::NONE;
      auto secCause = SecondaryCause::NONE;

      if constexpr (sizeof(elem) > 8)
        {
          for (unsigned n = 0; n < sizeof(elem); n += 8)
            {
              uint64_t dword = elem;
              bool forced = false;
              cause = determineStoreException(rs1, URV(eaddr), eaddr, dword, secCause,
                                              forced);
              if (cause != ExceptionCause::NONE)
                break;

              memory_.write(hartIx_, eaddr + n, dword);
              elem >>= 64;
            }
        }
      else
        {
          bool forced = false;
          if (determineStoreException(rs1, URV(eaddr), eaddr, elem, secCause, forced) !=
              ExceptionCause::NONE)
            memory_.write(hartIx_, eaddr, elem);
        }

      if (cause != ExceptionCause::NONE)
        {
          vecRegs_.setStartIndex(ix);
          csRegs_.write(CsrNumber::VSTART, PrivilegeMode::Machine, ix);
          initiateStoreException(cause, eaddr, secCause);
          break;
        }
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVsxei8_v(const DecodedInst* di)
{
  vectorStoreIndexed<uint8_t>(di, ElementWidth::Byte);
}


template <typename URV>
void
Hart<URV>::execVsxei16_v(const DecodedInst* di)
{
  vectorStoreIndexed<uint16_t>(di, ElementWidth::Half);
}


template <typename URV>
void
Hart<URV>::execVsxei32_v(const DecodedInst* di)
{
  vectorStoreIndexed<uint32_t>(di, ElementWidth::Word);
}


template <typename URV>
void
Hart<URV>::execVsxei64_v(const DecodedInst* di)
{
  vectorStoreIndexed<uint64_t>(di, ElementWidth::Word2);
}


template class WdRiscv::Hart<uint32_t>;
template class WdRiscv::Hart<uint64_t>;
