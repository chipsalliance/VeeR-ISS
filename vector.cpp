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
#include "instforms.hpp"
#include "DecodedInst.hpp"
#include "Hart.hpp"


// On pure 32-bit machines, use boost for 128-bit and higher integer types.
#if __x86_64__
  typedef __int128_t  Int128;
  typedef __uint128_t Uint128;
#else
  typedef boost::multiprecision::int128_t  Int128;
  typedef boost::multiprecision::uint128_t Uint128;
#endif

typedef boost::multiprecision::int256_t   Int256;
typedef boost::multiprecision::uint256_t  Uint256;
typedef boost::multiprecision::int512_t   Int512;
typedef boost::multiprecision::uint512_t  Uint512;
typedef boost::multiprecision::int1024_t  Int1024;
typedef boost::multiprecision::uint1024_t Uint1024;


// make_unsigned does not work on boost types -- compensate.
namespace std
{
  template <>
  struct
  make_unsigned<Int128>
  {
    typedef Uint128 type;
  };

  template <>
  struct
  make_unsigned<Int256>
  {
    typedef Uint256 type;
  };

  template <>
  struct
  make_unsigned<Int512>
  {
    typedef Uint512 type;
  };

  template <>
  struct
  make_unsigned<Int1024>
  {
    typedef Uint1024 type;
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
    temp *= b;
    temp >>= bits;
    result = TS(temp);
  }


  /// Specialized mulhsu for 1024-bit unsigned operands.
  template <>
  void mulhsu(const Int1024& a, const Uint1024& b, Int1024& result)
  {
    if (a > 0)
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
}


using namespace WdRiscv;


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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vadd_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vadd_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vadd_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vadd_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vadd_vv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vadd_vv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vadd_vv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vadd_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  unsigned rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  switch (sew)
    {
    case ElementWidth::Byte:
      vadd_vx<int8_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vadd_vx<int16_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vadd_vx<int32_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vadd_vx<int64_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vadd_vx<Int128>(vd, vs1, Int128(e2), group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vadd_vx<Int256>(vd, vs1, Int256(e2), group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vadd_vx<Int512>(vd, vs1, Int512(e2), group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vadd_vx<Int1024>(vd, vs1, Int1024(e2), group, start, elems, masked);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVadd_vi(const DecodedInst* di)
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
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vadd_vx<int8_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vadd_vx<int16_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vadd_vx<int32_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vadd_vx<int64_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vadd_vx<Int128>(vd, vs1, Int128(imm), group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vadd_vx<Int256>(vd, vs1, Int256(imm), group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vadd_vx<Int512>(vd, vs1, Int512(imm), group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vadd_vx<Int1024>(vd, vs1, Int1024(imm), group, start, elems, masked);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vsub_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vsub_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vsub_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vsub_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vsub_vv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vsub_vv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vsub_vv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vsub_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  unsigned rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vsub_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vsub_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vsub_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vsub_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vsub_vx<Int128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vsub_vx<Int256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vsub_vx<Int512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vsub_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  unsigned rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vrsub_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vrsub_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vrsub_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vrsub_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vrsub_vx<Int128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vrsub_vx<Int256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vrsub_vx<Int512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vrsub_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  int32_t imm = di->op2As<int32_t>();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vrsub_vi<int8_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vrsub_vi<int16_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vrsub_vi<int32_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vrsub_vi<int64_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vrsub_vi<Int128>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vrsub_vi<Int256>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vrsub_vi<Int512>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vrsub_vi<Int1024>(vd, vs1, imm, group, start, elems, masked);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
        {
          dest = DWT(e1);
          dest += e2;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vwadd_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vwadd_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vwadd_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vwadd_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vwadd_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vwadd_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vwadd_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVwadd_vv(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vwadd_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vwadd_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vwadd_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vwadd_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vwadd_vv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vwadd_vv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vwadd_vv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
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
          dest += e2;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  URV e2 =intRegs_.read(di->op2());

  switch (sew)
    {
    case ElementWidth::Byte:
      vwadd_vx<uint8_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vwadd_vx<uint16_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vwadd_vx<uint32_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vwadd_vx<uint64_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vwadd_vx<Uint128>(vd, vs1, Uint128(e2), group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vwadd_vx<Uint256>(vd, vs1, Uint256(e2), group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vwadd_vx<Uint512>(vd, vs1, Uint512(e2), group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVwadd_vx(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vwadd_vx<int8_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vwadd_vx<int16_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vwadd_vx<int32_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vwadd_vx<int64_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vwadd_vx<Int128>(vd, vs1, Int128(e2), group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vwadd_vx<Int256>(vd, vs1, Int256(e2), group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vwadd_vx<Int512>(vd, vs1, Int512(e2), group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVwsubu_vx(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vwadd_vx<uint8_t>(vd, vs1, -uint8_t(e2), group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vwadd_vx<uint16_t>(vd, vs1, -uint16_t(e2), group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vwadd_vx<uint32_t>(vd, vs1, -uint32_t(e2), group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vwadd_vx<uint64_t>(vd, vs1, -uint64_t(e2), group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vwadd_vx<Uint128>(vd, vs1, -Uint128(e2), group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vwadd_vx<Uint256>(vd, vs1, Uint256(0)-Uint256(e2), group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vwadd_vx<Uint512>(vd, vs1, Uint256(0)-Uint512(e2), group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVwsub_vx(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vwadd_vx<int8_t>(vd, vs1, -e2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vwadd_vx<int16_t>(vd, vs1, -e2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vwadd_vx<int32_t>(vd, vs1, -e2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vwadd_vx<int64_t>(vd, vs1, -e2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vwadd_vx<Int128>(vd, vs1, -Int128(e2), group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vwadd_vx<Int256>(vd, vs1, -Int256(e2), group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vwadd_vx<Int512>(vd, vs1, -Int512(e2), group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
        {
          dest = DWT(e1);
          dest -= e2;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vwsub_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vwsub_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vwsub_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vwsub_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vwsub_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vwsub_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vwsub_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVwsub_vv(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vwsub_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vwsub_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vwsub_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vwsub_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vwsub_vv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vwsub_vv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vwsub_vv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
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

      if (vecRegs_.read(vs1, ix, doubleGroup, e1) and
          vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1;
          dest += e2;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vwadd_wv<uint8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vwadd_wv<uint16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vwadd_wv<uint32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vwadd_wv<uint64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vwadd_wv<Uint128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vwadd_wv<Uint256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vwadd_wv<Uint512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVwadd_wv(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vwadd_wv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vwadd_wv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vwadd_wv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vwadd_wv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vwadd_wv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vwadd_wv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vwadd_wv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVwaddu_wx(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vadd_vx<uint16_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vadd_vx<uint32_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vadd_vx<uint64_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vadd_vx<Uint128>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vadd_vx<Uint256>(vd, vs1, Uint256(e2), group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vadd_vx<Uint512>(vd, vs1, Uint512(e2), group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vadd_vx<Uint1024>(vd, vs1, Uint1024(e2), group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVwadd_wx(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vadd_vx<int16_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vadd_vx<int32_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vadd_vx<int64_t>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vadd_vx<Int128>(vd, vs1, e2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vadd_vx<Int256>(vd, vs1, Int256(e2), group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vadd_vx<Int512>(vd, vs1, Int512(e2), group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vadd_vx<Int1024>(vd, vs1, Int1024(e2), group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVwsubu_wx(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vadd_vx<uint16_t>(vd, vs1, -e2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vadd_vx<uint32_t>(vd, vs1, -e2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vadd_vx<uint64_t>(vd, vs1, -e2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vadd_vx<Uint128>(vd, vs1, -e2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vadd_vx<Uint256>(vd, vs1, Uint256(0)-Uint256(e2), group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vadd_vx<Uint512>(vd, vs1, Uint512(0)-Uint512(e2), group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vadd_vx<Uint1024>(vd, vs1, Uint1024(0)-Uint1024(e2), group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVwsub_wx(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vadd_vx<int16_t>(vd, vs1, -e2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vadd_vx<int32_t>(vd, vs1, -e2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vadd_vx<int64_t>(vd, vs1, -e2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vadd_vx<Int128>(vd, vs1, -e2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vadd_vx<Int256>(vd, vs1, -Int256(e2), group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vadd_vx<Int512>(vd, vs1, -Int512(e2), group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vadd_vx<Int1024>(vd, vs1, -Int1024(e2), group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
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

      if (vecRegs_.read(vs1, ix, doubleGroup, e1) and
          vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1;
          dest -= e2;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vwsub_wv<uint8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vwsub_wv<uint16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vwsub_wv<uint32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vwsub_wv<uint64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vwsub_wv<Uint128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vwsub_wv<Uint256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vwsub_wv<Uint512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVwsub_wv(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vwsub_wv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vwsub_wv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vwsub_wv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vwsub_wv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vwsub_wv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vwsub_wv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vwsub_wv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      illegalInst(di);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vminu_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vminu_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vminu_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vminu_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vminu_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vminu_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vminu_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vminu_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  unsigned rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vminu_vx<uint8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vminu_vx<uint16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vminu_vx<uint32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vminu_vx<uint64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vminu_vx<Uint128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vminu_vx<Uint256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vminu_vx<Uint512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vminu_vx<Uint1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmin_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vmin_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vmin_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vmin_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vmin_vv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vmin_vv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vmin_vv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vmin_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  unsigned rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmin_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vmin_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vmin_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vmin_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vmin_vx<Int128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vmin_vx<Int256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vmin_vx<Int512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vmin_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmaxu_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vmaxu_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vmaxu_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vmaxu_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vmaxu_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vmaxu_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vmaxu_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vmaxu_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  unsigned rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmaxu_vx<uint8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vmaxu_vx<uint16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vmaxu_vx<uint32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vmaxu_vx<uint64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vmaxu_vx<Uint128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vmaxu_vx<Uint256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vmaxu_vx<Uint512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vmaxu_vx<Uint1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmax_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vmax_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vmax_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vmax_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vmax_vv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vmax_vv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vmax_vv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vmax_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  unsigned rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmax_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vmax_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vmax_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vmax_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vmax_vx<Int128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vmax_vx<Int256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vmax_vx<Int512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vmax_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vand_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vand_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vand_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vand_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vand_vv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vand_vv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vand_vv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vand_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  unsigned rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vand_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vand_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vand_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vand_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vand_vx<Int128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vand_vx<Int256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vand_vx<Int512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vand_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  int32_t imm = di->op2As<int32_t>();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vand_vi<int8_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vand_vi<int16_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vand_vi<int32_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vand_vi<int64_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vand_vi<Int128>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vand_vi<Int256>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vand_vi<Int512>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vand_vi<Int1024>(vd, vs1, imm, group, start, elems, masked);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vor_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vor_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vor_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vor_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vor_vv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vor_vv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vor_vv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vor_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  unsigned rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vor_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vor_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vor_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vor_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vor_vx<Int128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vor_vx<Int256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vor_vx<Int512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vor_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  int32_t imm = di->op2As<int32_t>();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vor_vi<int8_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vor_vi<int16_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vor_vi<int32_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vor_vi<int64_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vor_vi<Int128>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vor_vi<Int256>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vor_vi<Int512>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vor_vi<Int1024>(vd, vs1, imm, group, start, elems, masked);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vxor_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vxor_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vxor_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vxor_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vxor_vv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vxor_vv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vxor_vv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vxor_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  unsigned rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vxor_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vxor_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vxor_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vxor_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vxor_vx<Int128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vxor_vx<Int256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vxor_vx<Int512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vxor_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  int32_t imm = di->op2As<int32_t>();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vxor_vi<int8_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vxor_vi<int16_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vxor_vi<int32_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vxor_vi<int64_t>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vxor_vi<Int128>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vxor_vi<Int256>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vxor_vi<Int512>(vd, vs1, imm, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vxor_vi<Int1024>(vd, vs1, imm, group, start, elems, masked);
      break;
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

  switch (sew)
    {
    case ElementWidth::Byte:
      vrgather_vv<uint8_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vrgather_vv<uint16_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Word:
      vrgather_vv<uint32_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vrgather_vv<uint64_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vrgather_vv<Uint128>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vrgather_vv<Uint256>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vrgather_vv<Uint512>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vrgather_vv<Uint1024>(vd, vs1, vs2, group, start, elems);
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vrgather_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                       unsigned start, unsigned elems)
{
  unsigned errors = 0;

  URV rv2 = intRegs_.read(rs2);
  URV vs1Ix = rv2;    // FIX: if rv2 > VLMAX, rv2 = VLMAX

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

  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  unsigned rs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vrgather_vx<uint8_t>(vd, vs1, rs2, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vrgather_vx<uint16_t>(vd, vs1, rs2, group, start, elems);
      break;

    case ElementWidth::Word:
      vrgather_vx<uint32_t>(vd, vs1, rs2, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vrgather_vx<uint64_t>(vd, vs1, rs2, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vrgather_vx<Uint128>(vd, vs1, rs2, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vrgather_vx<Uint256>(vd, vs1, rs2, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vrgather_vx<Uint512>(vd, vs1, rs2, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vrgather_vx<Uint1024>(vd, vs1, rs2, group, start, elems);
      break;
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

  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  uint32_t imm = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vrgather_vi<uint8_t>(vd, vs1, imm, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vrgather_vi<uint16_t>(vd, vs1, imm, group, start, elems);
      break;

    case ElementWidth::Word:
      vrgather_vi<uint32_t>(vd, vs1, imm, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vrgather_vi<uint64_t>(vd, vs1, imm, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vrgather_vi<Uint128>(vd, vs1, imm, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vrgather_vi<Uint256>(vd, vs1, imm, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vrgather_vi<Uint512>(vd, vs1, imm, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vrgather_vi<Uint1024>(vd, vs1, imm, group, start, elems);
      break;
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

  switch (sew)
    {
    case ElementWidth::Byte:
      vrgatherei16_vv<uint8_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vrgatherei16_vv<uint16_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Word:
      vrgatherei16_vv<uint32_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vrgatherei16_vv<uint64_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vrgatherei16_vv<Uint128>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vrgatherei16_vv<Uint256>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vrgatherei16_vv<Uint512>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vrgatherei16_vv<Uint1024>(vd, vs1, vs2, group, start, elems);
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vcompress_vm(unsigned vd, unsigned vs1, unsigned vs2,
                        unsigned group, unsigned start, unsigned elems)
{
  unsigned errors = 0;

  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  unsigned destIx = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (vecRegs_.read(vs2, ix, group, e2))
        {
          if (e2 & 1)
            {
              if (vecRegs_.read(vs1, ix, group, e1))
                {
                  dest = e1;
                  if (not vecRegs_.write(vd, destIx++, group, dest))
                    errors++;
                }
            }
        }
      else
        errors++;
    }

  assert(errors == 0);
}


template <typename URV>
void
Hart<URV>::execVcompress_vm(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  if (di->isMasked())
    {
      illegalInst(di);
      return;
    }

  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vcompress_vm<uint8_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vcompress_vm<uint16_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Word:
      vcompress_vm<uint32_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vcompress_vm<uint64_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vcompress_vm<Uint128>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vcompress_vm<Uint256>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vcompress_vm<Uint512>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vcompress_vm<Uint1024>(vd, vs1, vs2, group, start, elems);
      break;
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

  switch (sew)
    {
    case ElementWidth::Byte:
      vredsum_vs<int8_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vredsum_vs<int16_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Word:
      vredsum_vs<int32_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vredsum_vs<int64_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vredsum_vs<Int128>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vredsum_vs<Int256>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vredsum_vs<Int512>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vredsum_vs<Int1024>(vd, vs1, vs2, group, start, elems);
      break;
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

  switch (sew)
    {
    case ElementWidth::Byte:
      vredand_vs<uint8_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vredand_vs<uint16_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Word:
      vredand_vs<uint32_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vredand_vs<uint64_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vredand_vs<Uint128>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vredand_vs<Uint256>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vredand_vs<Uint512>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vredand_vs<Uint1024>(vd, vs1, vs2, group, start, elems);
      break;
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

  switch (sew)
    {
    case ElementWidth::Byte:
      vredor_vs<uint8_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vredor_vs<uint16_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Word:
      vredor_vs<uint32_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vredor_vs<uint64_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vredor_vs<Uint128>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vredor_vs<Uint256>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vredor_vs<Uint512>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vredor_vs<Uint1024>(vd, vs1, vs2, group, start, elems);
      break;
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

  switch (sew)
    {
    case ElementWidth::Byte:
      vredxor_vs<uint8_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vredxor_vs<uint16_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Word:
      vredxor_vs<uint32_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vredxor_vs<uint64_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vredxor_vs<Uint128>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vredxor_vs<Uint256>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vredxor_vs<Uint512>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vredxor_vs<Uint1024>(vd, vs1, vs2, group, start, elems);
      break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vredminu_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
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

  switch (sew)
    {
    case ElementWidth::Byte:
      vredminu_vs<uint8_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vredminu_vs<uint16_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Word:
      vredminu_vs<uint32_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vredminu_vs<uint64_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vredminu_vs<Uint128>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vredminu_vs<Uint256>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vredminu_vs<Uint512>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vredminu_vs<Uint1024>(vd, vs1, vs2, group, start, elems);
      break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vredmin_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
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

  switch (sew)
    {
    case ElementWidth::Byte:
      vredmin_vs<int8_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vredmin_vs<int16_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Word:
      vredmin_vs<int32_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vredmin_vs<int64_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vredmin_vs<Int128>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vredmin_vs<Int256>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vredmin_vs<Int512>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vredmin_vs<Int1024>(vd, vs1, vs2, group, start, elems);
      break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vredmaxu_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
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

  switch (sew)
    {
    case ElementWidth::Byte:
      vredmaxu_vs<uint8_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vredmaxu_vs<uint16_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Word:
      vredmaxu_vs<uint32_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vredmaxu_vs<uint64_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vredmaxu_vs<Uint128>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vredmaxu_vs<Uint256>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vredmaxu_vs<Uint512>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vredmaxu_vs<Uint1024>(vd, vs1, vs2, group, start, elems);
      break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
void
Hart<URV>::vredmax_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
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

  switch (sew)
    {
    case ElementWidth::Byte:
      vredmax_vs<int8_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vredmax_vs<int16_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Word:
      vredmax_vs<int32_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vredmax_vs<int64_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vredmax_vs<Int128>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vredmax_vs<Int256>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vredmax_vs<Int512>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vredmax_vs<Int1024>(vd, vs1, vs2, group, start, elems);
      break;
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

  uint8_t* vdData = vecRegs_.getVecBytes(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecBytes(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecBytes(di->op2());
  if (not vs1Data or not vs2Data or not vdData)
    assert(0);

  // If bits indices are byte aligned process bytes
  if ((start & 7) == 0 and (elems & 7) == 0)
    {
      start = start >> 3;
      elems = elems >> 3;
      for (unsigned i = start; i < elems; ++i)
        vdData[i] = vs1Data[i] & vs2Data[i];
      return;
    }

  // Bit indices are not byte aligned.
  for (unsigned i = start; i < elems; ++i)
    {
      unsigned byteIx = i >> 3;
      unsigned bitIx = i & 7; // Bit index in byte
      uint8_t mask = 1 << bitIx;
      vdData[i] = ( (vdData[byteIx] & ~mask) |
                    (vs1Data[byteIx] & vs2Data[byteIx] & mask) );
    }
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

  uint8_t* vdData = vecRegs_.getVecBytes(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecBytes(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecBytes(di->op2());
  if (not vs1Data or not vs2Data or not vdData)
    assert(0);

  // If bits indices are byte aligned process bytes
  if ((start & 7) == 0 and (elems & 7) == 0)
    {
      start = start >> 3;
      elems = elems >> 3;
      for (unsigned i = start; i < elems; ++i)
        vdData[i] = ~ (vs1Data[i] & vs2Data[i]);
      return;
    }

  // Bit indices are not byte aligned.
  for (unsigned i = start; i < elems; ++i)
    {
      unsigned byteIx = i >> 3;
      unsigned bitIx = i & 7; // Bit index in byte
      uint8_t mask = 1 << bitIx;
      vdData[i] = ( (vdData[byteIx] & ~mask) |
                    (~(vs1Data[byteIx] & vs2Data[byteIx]) & mask) );
    }
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

  uint8_t* vdData = vecRegs_.getVecBytes(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecBytes(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecBytes(di->op2());
  if (not vs1Data or not vs2Data or not vdData)
    assert(0);

  // If bits indices are byte aligned process bytes
  if ((start & 7) == 0 and (elems & 7) == 0)
    {
      start = start >> 3;
      elems = elems >> 3;
      for (unsigned i = start; i < elems; ++i)
        vdData[i] = vs1Data[i] & ~vs2Data[i];
      return;
    }

  // Bit indices are not byte aligned.
  for (unsigned i = start; i < elems; ++i)
    {
      unsigned byteIx = i >> 3;
      unsigned bitIx = i & 7; // Bit index in byte
      uint8_t mask = 1 << bitIx;
      vdData[i] = ( (vdData[byteIx] & ~mask) |
                    ((vs1Data[byteIx] & ~vs2Data[byteIx]) & mask) );
    }
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

  uint8_t* vdData = vecRegs_.getVecBytes(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecBytes(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecBytes(di->op2());
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

  // Bit indices are not byte aligned.
  for (unsigned i = start; i < elems; ++i)
    {
      unsigned byteIx = i >> 3;
      unsigned bitIx = i & 7; // Bit index in byte
      uint8_t mask = 1 << bitIx;
      vdData[i] = ( (vdData[byteIx] & ~mask) |
                    ((vs1Data[byteIx] ^ vs2Data[byteIx]) & mask) );
    }
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

  uint8_t* vdData = vecRegs_.getVecBytes(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecBytes(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecBytes(di->op2());
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

  // Bit indices are not byte aligned.
  for (unsigned i = start; i < elems; ++i)
    {
      unsigned byteIx = i >> 3;
      unsigned bitIx = i & 7; // Bit index in byte
      uint8_t mask = 1 << bitIx;
      vdData[i] = ( (vdData[byteIx] & ~mask) |
                    ((vs1Data[byteIx] | vs2Data[byteIx]) & mask) );
    }
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

  uint8_t* vdData = vecRegs_.getVecBytes(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecBytes(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecBytes(di->op2());
  if (not vs1Data or not vs2Data or not vdData)
    assert(0);

  // If bits indices are byte aligned process bytes
  if ((start & 7) == 0 and (elems & 7) == 0)
    {
      start = start >> 3;
      elems = elems >> 3;
      for (unsigned i = start; i < elems; ++i)
        vdData[i] = ~(vs1Data[i] | vs2Data[i]);
      return;
    }

  // Bit indices are not byte aligned.
  for (unsigned i = start; i < elems; ++i)
    {
      unsigned byteIx = i >> 3;
      unsigned bitIx = i & 7; // Bit index in byte
      uint8_t mask = 1 << bitIx;
      vdData[i] = ( (vdData[byteIx] & ~mask) |
                    ( ~(vs1Data[byteIx] | vs2Data[byteIx]) & mask) );
    }
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

  uint8_t* vdData = vecRegs_.getVecBytes(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecBytes(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecBytes(di->op2());
  if (not vs1Data or not vs2Data or not vdData)
    assert(0);

  // If bits indices are byte aligned process bytes
  if ((start & 7) == 0 and (elems & 7) == 0)
    {
      start = start >> 3;
      elems = elems >> 3;
      for (unsigned i = start; i < elems; ++i)
        vdData[i] = vs1Data[i] | ~vs2Data[i];
      return;
    }

  // Bit indices are not byte aligned.
  for (unsigned i = start; i < elems; ++i)
    {
      unsigned byteIx = i >> 3;
      unsigned bitIx = i & 7; // Bit index in byte
      uint8_t mask = 1 << bitIx;
      vdData[i] = ( (vdData[byteIx] & ~mask) |
                    ((vs1Data[byteIx] | ~vs2Data[byteIx]) & mask) );
    }
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

  uint8_t* vdData = vecRegs_.getVecBytes(di->op0());
  uint8_t* vs1Data = vecRegs_.getVecBytes(di->op1());
  uint8_t* vs2Data = vecRegs_.getVecBytes(di->op2());
  if (not vs1Data or not vs2Data or not vdData)
    assert(0);

  // If bits indices are byte aligned process bytes
  if ((start & 7) == 0 and (elems & 7) == 0)
    {
      start = start >> 3;
      elems = elems >> 3;
      for (unsigned i = start; i < elems; ++i)
        vdData[i] = vs1Data[i] ^ ~vs2Data[i];
      return;
    }

  // Bit indices are not byte aligned.
  for (unsigned i = start; i < elems; ++i)
    {
      unsigned byteIx = i >> 3;
      unsigned bitIx = i & 7; // Bit index in byte
      uint8_t mask = 1 << bitIx;
      vdData[i] = ( (vdData[byteIx] & ~mask) |
                    ((vs1Data[byteIx] ^ ~vs2Data[byteIx]) & mask) );
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(), rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vslideup<uint8_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vslideup<uint16_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vslideup<uint32_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vslideup<uint64_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vslideup<Uint128>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vslideup<Uint256>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vslideup<Uint512>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vslideup<Uint1024>(vd, vs1, amount, group, start, elems, masked);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVslideup_vi(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(), imm = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vslideup<uint8_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vslideup<uint16_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vslideup<uint32_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vslideup<uint64_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vslideup<Uint128>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vslideup<Uint256>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vslideup<Uint512>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vslideup<Uint1024>(vd, vs1, amount, group, start, elems, masked);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVslide1up_vx(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(), rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

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

    case ElementWidth::HalfWord:
      vslideup<uint16_t>(vd, vs1, amount, group, start, elems, masked);
      vecRegs_.write(vd, 0, group, int16_t(replacement));
      break;

    case ElementWidth::Word:
      vslideup<uint32_t>(vd, vs1, amount, group, start, elems, masked);
      vecRegs_.write(vd, 0, group, int32_t(replacement));
      break;

    case ElementWidth::DoubleWord:
      vslideup<uint64_t>(vd, vs1, amount, group, start, elems, masked);
      vecRegs_.write(vd, 0, group, int64_t(replacement));
      break;

    case ElementWidth::QuadWord:
      vslideup<Uint128>(vd, vs1, amount, group, start, elems, masked);
      vecRegs_.write(vd, 0, group, Int128(replacement));
      break;

    case ElementWidth::OctWord:
      vslideup<Uint256>(vd, vs1, amount, group, start, elems, masked);
      vecRegs_.write(vd, 0, group, Int256(replacement));
      break;

    case ElementWidth::HalfKbits:
      vslideup<Uint512>(vd, vs1, amount, group, start, elems, masked);
      vecRegs_.write(vd, 0, group, Int512(replacement));
      break;

    case ElementWidth::Kbits:
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(), rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vslidedown<uint8_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vslidedown<uint16_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vslidedown<uint32_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vslidedown<uint64_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vslidedown<Uint128>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vslidedown<Uint256>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vslidedown<Uint512>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vslidedown<Uint1024>(vd, vs1, amount, group, start, elems, masked);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVslidedown_vi(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(), imm = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

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

  switch (sew)
    {
    case ElementWidth::Byte:
      vslidedown<uint8_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vslidedown<uint16_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vslidedown<uint32_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vslidedown<uint64_t>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vslidedown<Uint128>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vslidedown<Uint256>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vslidedown<Uint512>(vd, vs1, amount, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vslidedown<Uint1024>(vd, vs1, amount, group, start, elems, masked);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmul_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vmul_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vmul_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vmul_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vmul_vv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vmul_vv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vmul_vv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vmul_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(), vs1 = di->op1(), rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmul_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vmul_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vmul_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vmul_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vmul_vx<Int128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vmul_vx<Int256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vmul_vx<Int512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vmul_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmulh_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vmulh_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vmulh_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vmulh_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vmulh_vv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vmulh_vv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vmulh_vv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vmul_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked);
      assert(0);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmulh_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vmulh_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vmulh_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vmulh_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vmulh_vx<Int128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vmulh_vx<Int256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vmulh_vx<Int512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vmulh_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVmulhu_vv(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmulh_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vmulh_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vmulh_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vmulh_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vmulh_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vmulh_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vmulh_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vmulh_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmulhu_vx<uint8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vmulhu_vx<uint16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vmulhu_vx<uint32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vmulhu_vx<uint64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vmulhu_vx<Uint128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vmulhu_vx<Uint256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vmulhu_vx<Uint512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vmulhu_vx<Uint1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmulhsu_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vmulhsu_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vmulhsu_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vmulhsu_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vmulhsu_vv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vmulhsu_vv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vmulhsu_vv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vmulhsu_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmulhsu_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vmulhsu_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vmulhsu_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vmulhsu_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vmulhsu_vx<Int128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vmulhsu_vx<Int256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vmulhsu_vx<Int512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vmulhsu_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
        {
          dest = ~ ELEM_TYPE(0); // divide by zero result
          if (e2 != 0)
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vdivu_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vdivu_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vdivu_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vdivu_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vdivu_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vdivu_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vdivu_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vdivu_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
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
          if (e2 != 0)
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vdivu_vv<uint8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vdivu_vv<uint16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vdivu_vv<uint32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vdivu_vv<uint64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vdivu_vv<Uint128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vdivu_vv<Uint256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vdivu_vv<Uint512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vdivu_vv<Uint1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vdiv_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  unsigned elemBits = integerWidth<ELEM_TYPE> ();

  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  ELEM_TYPE minInt = ELEM_TYPE(1) << (elemBits - 1);
  ELEM_TYPE negOne = ELEM_TYPE(-1);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
        {
          dest = negOne; // Divide by zero result
          if (e2 != 0)
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vdiv_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vdiv_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vdiv_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vdiv_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vdiv_vv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vdiv_vv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vdiv_vv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vdiv_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vdiv_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  unsigned elemBits = integerWidth<ELEM_TYPE> ();

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
          if (e2 != 0)
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vdiv_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vdiv_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vdiv_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vdiv_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vdiv_vx<Int128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vdiv_vx<Int256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vdiv_vx<Int512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vdiv_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
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

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1; // divide by zero result
          if (e2 != 0)
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vremu_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vremu_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vremu_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vremu_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vremu_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vremu_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vremu_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vremu_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
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
          if (e2 != 0)
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vremu_vv<uint8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vremu_vv<uint16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vremu_vv<uint32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vremu_vv<uint64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vremu_vv<Uint128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vremu_vv<Uint256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vremu_vv<Uint512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vremu_vv<Uint1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vrem_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  unsigned elemBits = integerWidth<ELEM_TYPE> ();

  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;

  ELEM_TYPE minInt = ELEM_TYPE(1) << (elemBits - 1);
  ELEM_TYPE negOne = ELEM_TYPE(-1);

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
        {
          dest = e1; // Divide by zero remainder
          if (e2 != 0)
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  vs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vrem_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vrem_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vrem_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vrem_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vrem_vv<Int128>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vrem_vv<Int256>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vrem_vv<Int512>(vd, vs1, vs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vrem_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked);
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
void
Hart<URV>::vrem_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  unsigned elemBits = integerWidth<ELEM_TYPE> ();

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
          if (e2 != 0)
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  unsigned vd = di->op0(),  vs1 = di->op1(),  rs2 = di->op2();
  if (masked and vd == 0)  // Dest register cannot overlap mask register v0
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vremu_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfWord:
      vremu_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Word:
      vremu_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vremu_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vremu_vx<Int128>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vremu_vx<Int256>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vremu_vx<Int512>(vd, vs1, rs2, group, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vremu_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked);
      break;
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
Hart<URV>::execVsext_f2(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  ElementWidth sew = vecRegs_.elemWidth();
  ElementWidth eew = sew;  // Effective elem width of source.
  switch (sew)
    {
    case ElementWidth::Byte:       illegalInst(di);                return;
    case ElementWidth::HalfWord:   eew = ElementWidth::Byte;       break;
    case ElementWidth::Word:       eew = ElementWidth::HalfWord;   break;
    case ElementWidth::DoubleWord: eew = ElementWidth::Word;       break;
    case ElementWidth::QuadWord:   eew = ElementWidth::DoubleWord; break;
    case ElementWidth::OctWord:    eew = ElementWidth::QuadWord;   break;
    case ElementWidth::HalfKbits:  eew = ElementWidth::OctWord;    break;
    case ElementWidth::Kbits:      eew = ElementWidth::HalfKbits;  break;
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
    case ElementWidth::Byte:
      illegalInst(di);
      break;

    case ElementWidth::HalfWord:
      vsext<int16_t,int8_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::Word:
      vsext<int32_t, int16_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vsext<int64_t, int32_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vsext<Int128, int64_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vsext<Int256, Int128>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vsext<Int512, Int256>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vsext<Int1024, Int512>(vd, vs1, group, fromGroup, start, elems, masked);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVsext_f4(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  ElementWidth sew = vecRegs_.elemWidth();
  ElementWidth eew = sew;  // Effective elem width of source.
  switch (sew)
    {
    case ElementWidth::Byte:       illegalInst(di);                return;
    case ElementWidth::HalfWord:   illegalInst(di);                return;
    case ElementWidth::Word:       eew = ElementWidth::Byte;       break;
    case ElementWidth::DoubleWord: eew = ElementWidth::HalfWord;   break;
    case ElementWidth::QuadWord:   eew = ElementWidth::Word;       break;
    case ElementWidth::OctWord:    eew = ElementWidth::DoubleWord; break;
    case ElementWidth::HalfKbits:  eew = ElementWidth::QuadWord;   break;
    case ElementWidth::Kbits:      eew = ElementWidth::OctWord;    break;
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
    case ElementWidth::Byte:
    case ElementWidth::HalfWord:
      illegalInst(di);
      break;

    case ElementWidth::Word:
      vsext<int32_t, int8_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vsext<int64_t, int16_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vsext<Int128, int32_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vsext<Int256, int64_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vsext<Int512, Int128>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vsext<Int1024, Int256>(vd, vs1, group, fromGroup, start, elems, masked);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVsext_f8(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  ElementWidth sew = vecRegs_.elemWidth();
  ElementWidth eew = sew;  // Effective elem width of source.
  switch (sew)
    {
    case ElementWidth::Byte:       illegalInst(di);                return;
    case ElementWidth::HalfWord:   illegalInst(di);                return;
    case ElementWidth::Word:       illegalInst(di);                return;
    case ElementWidth::DoubleWord: eew = ElementWidth::Byte;       break;
    case ElementWidth::QuadWord:   eew = ElementWidth::HalfWord;   break;
    case ElementWidth::OctWord:    eew = ElementWidth::Word;       break;
    case ElementWidth::HalfKbits:  eew = ElementWidth::DoubleWord; break;
    case ElementWidth::Kbits:      eew = ElementWidth::QuadWord;   break;
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
    case ElementWidth::Byte:
    case ElementWidth::HalfWord:
    case ElementWidth::Word:
    case ElementWidth::DoubleWord:
      vsext<int64_t, int8_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vsext<Int128, int16_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vsext<Int256, int32_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vsext<Int512, int64_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vsext<Int1024, Int128>(vd, vs1, group, fromGroup, start, elems, masked);
      break;
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
Hart<URV>::execVzext_f2(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  ElementWidth sew = vecRegs_.elemWidth();
  ElementWidth eew = sew;  // Effective elem width of source.
  switch (sew)
    {
    case ElementWidth::Byte:       illegalInst(di);                return;
    case ElementWidth::HalfWord:   eew = ElementWidth::Byte;       break;
    case ElementWidth::Word:       eew = ElementWidth::HalfWord;   break;
    case ElementWidth::DoubleWord: eew = ElementWidth::Word;       break;
    case ElementWidth::QuadWord:   eew = ElementWidth::DoubleWord; break;
    case ElementWidth::OctWord:    eew = ElementWidth::QuadWord;   break;
    case ElementWidth::HalfKbits:  eew = ElementWidth::OctWord;    break;
    case ElementWidth::Kbits:      eew = ElementWidth::HalfKbits;  break;
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
    case ElementWidth::Byte:
      illegalInst(di);
      break;

    case ElementWidth::HalfWord:
      vzext<uint16_t,uint8_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::Word:
      vzext<uint32_t, uint16_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vzext<uint64_t, uint32_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vzext<Uint128, uint64_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vzext<Uint256, Uint128>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vzext<Uint512, Uint256>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vzext<Uint1024, Uint512>(vd, vs1, group, fromGroup, start, elems, masked);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVzext_f4(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  ElementWidth sew = vecRegs_.elemWidth();
  ElementWidth eew = sew;  // Effective elem width of source.
  switch (sew)
    {
    case ElementWidth::Byte:       illegalInst(di);                return;
    case ElementWidth::HalfWord:   illegalInst(di);                return;
    case ElementWidth::Word:       eew = ElementWidth::Byte;       break;
    case ElementWidth::DoubleWord: eew = ElementWidth::HalfWord;   break;
    case ElementWidth::QuadWord:   eew = ElementWidth::Word;       break;
    case ElementWidth::OctWord:    eew = ElementWidth::DoubleWord; break;
    case ElementWidth::HalfKbits:  eew = ElementWidth::QuadWord;   break;
    case ElementWidth::Kbits:      eew = ElementWidth::OctWord;    break;
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
    case ElementWidth::Byte:
    case ElementWidth::HalfWord:
      illegalInst(di);
      break;

    case ElementWidth::Word:
      vzext<uint32_t, uint8_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::DoubleWord:
      vzext<uint64_t, uint16_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vzext<Uint128, uint32_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vzext<Uint256, uint64_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vzext<Uint512, Uint128>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vzext<Uint1024, Uint256>(vd, vs1, group, fromGroup, start, elems, masked);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVzext_f8(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

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

  ElementWidth sew = vecRegs_.elemWidth();
  ElementWidth eew = sew;  // Effective elem width of source.
  switch (sew)
    {
    case ElementWidth::Byte:       illegalInst(di);                return;
    case ElementWidth::HalfWord:   illegalInst(di);                return;
    case ElementWidth::Word:       illegalInst(di);                return;
    case ElementWidth::DoubleWord: eew = ElementWidth::Byte;       break;
    case ElementWidth::QuadWord:   eew = ElementWidth::HalfWord;   break;
    case ElementWidth::OctWord:    eew = ElementWidth::Word;       break;
    case ElementWidth::HalfKbits:  eew = ElementWidth::DoubleWord; break;
    case ElementWidth::Kbits:      eew = ElementWidth::QuadWord;   break;
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
    case ElementWidth::Byte:
    case ElementWidth::HalfWord:
    case ElementWidth::Word:
    case ElementWidth::DoubleWord:
      vzext<uint64_t, uint8_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::QuadWord:
      vzext<Uint128, uint16_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::OctWord:
      vzext<Uint256, uint32_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::HalfKbits:
      vzext<Uint512, uint64_t>(vd, vs1, group, fromGroup, start, elems, masked);
      break;

    case ElementWidth::Kbits:
      vzext<Uint1024, Uint128>(vd, vs1, group, fromGroup, start, elems, masked);
      break;
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
      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
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

  switch (sew)
    {
    case ElementWidth::Byte:
      vmerge_vv<int8_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vmerge_vv<int16_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Word:
      vmerge_vv<int32_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vmerge_vv<int64_t>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vmerge_vv<Int128>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vmerge_vv<Int256>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vmerge_vv<Int512>(vd, vs1, vs2, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vmerge_vv<Int1024>(vd, vs1, vs2, group, start, elems);
      break;
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
  unsigned vd = di->op0();
  unsigned vs1 = di->op1();
  unsigned rs2 = di->op2();
  if (not masked)
    {
      illegalInst(di);
      return;
    }

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  SRV e2 = SRV(intRegs_.read(rs2));

  switch (sew)
    {
    case ElementWidth::Byte:
      vmerge_vx<int8_t>(vd, vs1, e2, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vmerge_vx<int16_t>(vd, vs1, e2, group, start, elems);
      break;

    case ElementWidth::Word:
      vmerge_vx<int32_t>(vd, vs1, e2, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vmerge_vx<int64_t>(vd, vs1, e2, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vmerge_vx<Int128>(vd, vs1, Int128(e2), group, start, elems);
      break;

    case ElementWidth::OctWord:
      vmerge_vx<Int256>(vd, vs1, Int256(e2), group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vmerge_vx<Int512>(vd, vs1, Int512(e2), group, start, elems);
      break;

    case ElementWidth::Kbits:
      vmerge_vx<Int1024>(vd, vs1, Int1024(e2), group, start, elems);
      break;
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

  switch (sew)
    {
    case ElementWidth::Byte:
      vmerge_vx<int8_t>(vd, vs1, imm, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vmerge_vx<int16_t>(vd, vs1, imm, group, start, elems);
      break;

    case ElementWidth::Word:
      vmerge_vx<int32_t>(vd, vs1, imm, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vmerge_vx<int64_t>(vd, vs1, imm, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vmerge_vx<Int128>(vd, vs1, Int128(imm), group, start, elems);
      break;

    case ElementWidth::OctWord:
      vmerge_vx<Int256>(vd, vs1, Int256(imm), group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vmerge_vx<Int512>(vd, vs1, Int512(imm), group, start, elems);
      break;

    case ElementWidth::Kbits:
      vmerge_vx<Int1024>(vd, vs1, Int1024(imm), group, start, elems);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVmv_x_s(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  if (masked)
    {
      illegalInst(di);   // Masked version reserved.
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

    case ElementWidth::HalfWord:
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

    case ElementWidth::DoubleWord:
      {
        int64_t val = 0;
        vecRegs_.read(vs1, 0, groupX8, val);
        intRegs_.write(rd, SRV(val));
      }
      break;

    case ElementWidth::QuadWord:
      {
        Int128 val = 0;
        vecRegs_.read(vs1, 0, groupX8, val);
        intRegs_.write(rd, SRV(val));
      }
      break;

    case ElementWidth::OctWord:
      {
        Int256 val = 0;
        vecRegs_.read(vs1, 0, groupX8, val);
        intRegs_.write(rd, SRV(val));
      }
      break;

    case ElementWidth::HalfKbits:
      {
        Int512 val = 0;
        vecRegs_.read(vs1, 0, groupX8, val);
        intRegs_.write(rd, SRV(val));
      }
      break;

    case ElementWidth::Kbits:
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  if (masked)
    {
      illegalInst(di);    // Masked version reserved.
      return;
    }

  unsigned vd = di->op0(), rs1 = di->op1(), groupX8 = 8;

  ElementWidth sew = vecRegs_.elemWidth();

  SRV val = intRegs_.read(rs1);

  switch (sew)
    {
    case ElementWidth::Byte:
      vecRegs_.write(vd, 0, groupX8, int8_t(val));
      break;

    case ElementWidth::HalfWord:
      vecRegs_.write(vd, 0, groupX8, int16_t(val));
      break;

    case ElementWidth::Word:
      vecRegs_.write(vd, 0, groupX8, int32_t(val));
      break;

    case ElementWidth::DoubleWord:
      vecRegs_.write(vd, 0, groupX8, int64_t(val));
      break;

    case ElementWidth::QuadWord:
      vecRegs_.write(vd, 0, groupX8, Int128(val));
      break;

    case ElementWidth::OctWord:
      vecRegs_.write(vd, 0, groupX8, Int256(val));
      break;

    case ElementWidth::HalfKbits:
      vecRegs_.write(vd, 0, groupX8, Int512(val));
      break;

    case ElementWidth::Kbits:
      vecRegs_.write(vd, 0, groupX8, Int1024(val));
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  assert(not masked);

  unsigned vd = di->op0(),  vs1 = di->op1();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  switch (sew)
    {
    case ElementWidth::Byte:
      vmv_v_v<int8_t>(vd, vs1, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vmv_v_v<int16_t>(vd, vs1, group, start, elems);
      break;

    case ElementWidth::Word:
      vmv_v_v<int32_t>(vd, vs1, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vmv_v_v<int64_t>(vd, vs1, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vmv_v_v<Int128>(vd, vs1, group, start, elems);
      break;

    case ElementWidth::OctWord:
      vmv_v_v<Int256>(vd, vs1, group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vmv_v_v<Int512>(vd, vs1, group, start, elems);
      break;

    case ElementWidth::Kbits:
      vmv_v_v<Int1024>(vd, vs1, group, start, elems);
      break;
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
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  assert(not masked);

  unsigned vd = di->op0();
  unsigned rs1 = di->op1();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  int e1 = SRV(intRegs_.read(rs1));

  switch (sew)
    {
    case ElementWidth::Byte:
      vmv_v_x<int8_t>(vd, e1, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vmv_v_x<int16_t>(vd, e1, group, start, elems);
      break;

    case ElementWidth::Word:
      vmv_v_x<int32_t>(vd, e1, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vmv_v_x<int64_t>(vd, e1, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vmv_v_x<Int128>(vd, Int128(e1), group, start, elems);
      break;

    case ElementWidth::OctWord:
      vmv_v_x<Int256>(vd, Int256(e1), group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vmv_v_x<Int512>(vd, Int512(e1), group, start, elems);
      break;

    case ElementWidth::Kbits:
      vmv_v_x<Int1024>(vd, Int1024(e1), group, start, elems);
      break;
    }
}


template <typename URV>
void
Hart<URV>::execVmv_v_i(const DecodedInst* di)
{
  if (not isVecLegal() or not vecRegs_.legalConfig())
    {
      illegalInst(di);
      return;
    }

  bool masked = di->isMasked();
  assert(not masked);

  unsigned vd = di->op0();

  unsigned group = vecRegs_.groupMultiplierX8(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  int e1 = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      vmv_v_x<int8_t>(vd, e1, group, start, elems);
      break;

    case ElementWidth::HalfWord:
      vmv_v_x<int16_t>(vd, e1, group, start, elems);
      break;

    case ElementWidth::Word:
      vmv_v_x<int32_t>(vd, e1, group, start, elems);
      break;

    case ElementWidth::DoubleWord:
      vmv_v_x<int64_t>(vd, e1, group, start, elems);
      break;

    case ElementWidth::QuadWord:
      vmv_v_x<Int128>(vd, Int128(e1), group, start, elems);
      break;

    case ElementWidth::OctWord:
      vmv_v_x<Int256>(vd, Int256(e1), group, start, elems);
      break;

    case ElementWidth::HalfKbits:
      vmv_v_x<Int512>(vd, Int512(e1), group, start, elems);
      break;

    case ElementWidth::Kbits:
      vmv_v_x<Int1024>(vd, Int1024(e1), group, start, elems);
      break;
    }
}


template class WdRiscv::Hart<uint32_t>;
template class WdRiscv::Hart<uint64_t>;
