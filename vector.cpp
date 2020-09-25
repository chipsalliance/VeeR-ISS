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


using namespace WdRiscv;


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vadd_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vadd_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vadd_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vadd_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vadd_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vadd_vv<Int128>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vadd_vv<Int256>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vadd_vv<Int512>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vadd_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vadd_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                        unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  SRV srv2 = intRegs_.read(rs2);

  ELEM_TYPE val2 = srv2;

  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 + val2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vadd_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vadd_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vadd_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vadd_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vadd_vx<Int128>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vadd_vx<Int256>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vadd_vx<Int512>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vadd_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vadd_vi(unsigned vd, unsigned vs1, int32_t imm, unsigned group,
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
          dest = e1 + val2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vadd_vi<int8_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vadd_vi<int16_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vadd_vi<int32_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vadd_vi<int64_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vadd_vi<Int128>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vadd_vi<Int256>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vadd_vi<Int512>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vadd_vi<Int1024>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vsub_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vsub_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vsub_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vsub_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vsub_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vsub_vv<Int128>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vsub_vv<Int256>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vsub_vv<Int512>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vsub_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vsub_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                        unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  SRV srv2 = intRegs_.read(rs2);

  ELEM_TYPE val2 = srv2;

  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 - val2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vsub_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vsub_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vsub_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vsub_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vsub_vx<Int128>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vsub_vx<Int256>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vsub_vx<Int512>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vsub_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vrsub_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                         unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  SRV srv2 = intRegs_.read(rs2);

  ELEM_TYPE val2 = srv2;

  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = val2 - e1;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vrsub_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vrsub_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vrsub_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vrsub_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vrsub_vx<Int128>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vrsub_vx<Int256>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vrsub_vx<Int512>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vrsub_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vrsub_vi(unsigned vd, unsigned vs1, int32_t imm, unsigned group,
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
          dest = val2 - e1;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vrsub_vi<int8_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vrsub_vi<int16_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vrsub_vi<int32_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vrsub_vi<int64_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vrsub_vi<Int128>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vrsub_vi<Int256>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vrsub_vi<Int512>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vrsub_vi<Int1024>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vminu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vminu_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vminu_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vminu_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vminu_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vminu_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vminu_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vminu_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vminu_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vminu_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                        unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  URV srv2 = intRegs_.read(rs2);  // Spec (sep 24, 2020) says this should be sign extended. We hope they come to their senses.

  ELEM_TYPE val2 = srv2;

  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 < val2 ? e1 : val2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vminu_vx<uint8_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vminu_vx<uint16_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vminu_vx<uint32_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vminu_vx<uint64_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vminu_vx<Uint128>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vminu_vx<Uint256>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vminu_vx<Uint512>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vminu_vx<Uint1024>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vmin_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vmin_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vmin_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vmin_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vmin_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vmin_vv<Int128>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vmin_vv<Int256>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vmin_vv<Int512>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vmin_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vmin_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                        unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  SRV srv2 = intRegs_.read(rs2);

  ELEM_TYPE val2 = srv2;

  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 < val2 ? e1 : val2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vmin_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vmin_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vmin_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vmin_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vmin_vx<Int128>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vmin_vx<Int256>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vmin_vx<Int512>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vmin_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vmaxu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vmaxu_vv<uint8_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vmaxu_vv<uint16_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vmaxu_vv<uint32_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vmaxu_vv<uint64_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vmaxu_vv<Uint128>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vmaxu_vv<Uint256>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vmaxu_vv<Uint512>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vmaxu_vv<Uint1024>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vmaxu_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                         unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  ELEM_TYPE val2 = intRegs_.read(rs2);   // Spec (sep 24, 2020) says this should be sign extended. We hope they come to their senses.

  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 > val2 ? e1 : val2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vmaxu_vx<uint8_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vmaxu_vx<uint16_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vmaxu_vx<uint32_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vmaxu_vx<uint64_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vmaxu_vx<Uint128>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vmaxu_vx<Uint256>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vmaxu_vx<Uint512>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vmaxu_vx<Uint1024>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vmax_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vmax_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vmax_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vmax_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vmax_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vmax_vv<Int128>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vmax_vv<Int256>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vmax_vv<Int512>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vmax_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vmax_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                        unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  SRV srv2 = intRegs_.read(rs2);

  ELEM_TYPE val2 = srv2;

  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 > val2 ? e1 : val2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vmax_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vmax_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vmax_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vmax_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vmax_vx<Int128>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vmax_vx<Int256>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vmax_vx<Int512>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vmax_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vand_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vand_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vand_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vand_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vand_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vand_vv<Int128>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vand_vv<Int256>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vand_vv<Int512>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vand_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vand_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                        unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  SRV srv2 = intRegs_.read(rs2);

  ELEM_TYPE val2 = srv2;

  ELEM_TYPE e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 & val2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vand_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vand_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vand_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vand_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vand_vx<Int128>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vand_vx<Int256>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vand_vx<Int512>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vand_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vand_vi(unsigned vd, unsigned vs1, int32_t imm, unsigned group,
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
          dest = e1 & val2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vand_vi<int8_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vand_vi<int16_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vand_vi<int32_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vand_vi<int64_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vand_vi<Int128>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vand_vi<Int256>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vand_vi<Int512>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vand_vi<Int1024>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vor_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vor_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vor_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vor_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vor_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vor_vv<Int128>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vor_vv<Int256>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vor_vv<Int512>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vor_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vor_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                       unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  SRV srv2 = intRegs_.read(rs2);

  ELEM_TYPE val2 = srv2;

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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vor_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vor_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vor_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vor_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vor_vx<Int128>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vor_vx<Int256>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vor_vx<Int512>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vor_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vor_vi(unsigned vd, unsigned vs1, int32_t imm, unsigned group,
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vor_vi<int8_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vor_vi<int16_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vor_vi<int32_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vor_vi<int64_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vor_vi<Int128>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vor_vi<Int256>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vor_vi<Int512>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vor_vi<Int1024>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vxor_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vxor_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vxor_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vxor_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vxor_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vxor_vv<Int128>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vxor_vv<Int256>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vxor_vv<Int512>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vxor_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vxor_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                        unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  SRV srv2 = intRegs_.read(rs2);

  ELEM_TYPE val2 = srv2;

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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vxor_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vxor_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vxor_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vxor_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vxor_vx<Int128>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vxor_vx<Int256>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vxor_vx<Int512>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vxor_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vxor_vi(unsigned vd, unsigned vs1, int32_t imm, unsigned group,
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vxor_vi<int8_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vxor_vi<int16_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vxor_vi<int32_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vxor_vi<int64_t>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vxor_vi<Int128>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vxor_vi<Int256>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vxor_vi<Int512>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vxor_vi<Int1024>(vd, vs1, imm, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vrgather_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vrgather_vv<uint8_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vrgather_vv<uint16_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vrgather_vv<uint32_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vrgather_vv<uint64_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vrgather_vv<Uint128>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vrgather_vv<Uint256>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vrgather_vv<Uint512>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vrgather_vv<Uint1024>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vrgather_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vrgather_vx<uint8_t>(vd, vs1, rs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vrgather_vx<uint16_t>(vd, vs1, rs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vrgather_vx<uint32_t>(vd, vs1, rs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vrgather_vx<uint64_t>(vd, vs1, rs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vrgather_vx<Uint128>(vd, vs1, rs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vrgather_vx<Uint256>(vd, vs1, rs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vrgather_vx<Uint512>(vd, vs1, rs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vrgather_vx<Uint1024>(vd, vs1, rs2, group, start, elems))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vrgather_vi<uint8_t>(vd, vs1, imm, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vrgather_vi<uint16_t>(vd, vs1, imm, group, start, elems))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vrgather_vi<uint32_t>(vd, vs1, imm, group, start, elems))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vrgather_vi<uint64_t>(vd, vs1, imm, group, start, elems))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vrgather_vi<Uint128>(vd, vs1, imm, group, start, elems))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vrgather_vi<Uint256>(vd, vs1, imm, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vrgather_vi<Uint512>(vd, vs1, imm, group, start, elems))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vrgather_vi<Uint1024>(vd, vs1, imm, group, start, elems))
        errors++;
      break;
    }
}




template <typename URV>
template <typename ELEM_TYPE>
bool Hart<URV>::vrgatherei16_vv(unsigned vd, unsigned vs1, unsigned vs2,
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vrgatherei16_vv<uint8_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vrgatherei16_vv<uint16_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vrgatherei16_vv<uint32_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vrgatherei16_vv<uint64_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vrgatherei16_vv<Uint128>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vrgatherei16_vv<Uint256>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vrgatherei16_vv<Uint512>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vrgatherei16_vv<Uint1024>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vcompress_vm<uint8_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vcompress_vm<uint16_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vcompress_vm<uint32_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vcompress_vm<uint64_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vcompress_vm<Uint128>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vcompress_vm<Uint256>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vcompress_vm<Uint512>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vcompress_vm<Uint1024>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;
    }
}



template <typename URV>
template<typename ELEM_TYPE>
bool
Hart<URV>::vredsum_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                      unsigned start, unsigned elems)
{
  unsigned errors = 0;

  ELEM_TYPE e2 = 0;
  if (vecRegs_.read(vs2, 0, 1, e2))
    return false;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result += e1;
    else
      errors++;

  if (errors == 0)
    if (not vecRegs_.write(vd, 0, group, result))
      errors++;

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vredsum_vs<int8_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vredsum_vs<int16_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vredsum_vs<int32_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vredsum_vs<int64_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vredsum_vs<Int128>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vredsum_vs<Int256>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vredsum_vs<Int512>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vredsum_vs<Int1024>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
bool
Hart<URV>::vredand_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                      unsigned start, unsigned elems)
{
  unsigned errors = 0;

  ELEM_TYPE e2 = 0;
  if (vecRegs_.read(vs2, 0, 1, e2))
    return false;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result = result & e1;
    else
      errors++;

  if (errors == 0)
    if (not vecRegs_.write(vd, 0, group, result))
      errors++;

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vredand_vs<uint8_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vredand_vs<uint16_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vredand_vs<uint32_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vredand_vs<uint64_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vredand_vs<Uint128>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vredand_vs<Uint256>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vredand_vs<Uint512>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vredand_vs<Uint1024>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
bool
Hart<URV>::vredor_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                     unsigned start, unsigned elems)
{
  unsigned errors = 0;

  ELEM_TYPE e2 = 0;
  if (vecRegs_.read(vs2, 0, 1, e2))
    return false;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result = result | e1;
    else
      errors++;

  if (errors == 0)
    if (not vecRegs_.write(vd, 0, group, result))
      errors++;

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vredor_vs<uint8_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vredor_vs<uint16_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vredor_vs<uint32_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vredor_vs<uint64_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vredor_vs<Uint128>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vredor_vs<Uint256>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vredor_vs<Uint512>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vredor_vs<Uint1024>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
bool
Hart<URV>::vredxor_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                      unsigned start, unsigned elems)
{
  unsigned errors = 0;

  ELEM_TYPE e2 = 0;
  if (vecRegs_.read(vs2, 0, 1, e2))
    return false;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result = result ^ e1;
    else
      errors++;

  if (errors == 0)
    if (not vecRegs_.write(vd, 0, group, result))
      errors++;

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vredxor_vs<uint8_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vredxor_vs<uint16_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vredxor_vs<uint32_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vredxor_vs<uint64_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vredxor_vs<Uint128>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vredxor_vs<Uint256>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vredxor_vs<Uint512>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vredxor_vs<Uint1024>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
bool
Hart<URV>::vredminu_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                       unsigned start, unsigned elems)
{
  unsigned errors = 0;

  ELEM_TYPE e2 = 0;
  if (vecRegs_.read(vs2, 0, 1, e2))
    return false;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result = result < e1 ? result : e1;
    else
      errors++;

  if (errors == 0)
    if (not vecRegs_.write(vd, 0, group, result))
      errors++;

  return errors == 0;

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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vredminu_vs<uint8_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vredminu_vs<uint16_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vredminu_vs<uint32_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vredminu_vs<uint64_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vredminu_vs<Uint128>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vredminu_vs<Uint256>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vredminu_vs<Uint512>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vredminu_vs<Uint1024>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
bool
Hart<URV>::vredmin_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                      unsigned start, unsigned elems)
{
  unsigned errors = 0;

  ELEM_TYPE e2 = 0;
  if (vecRegs_.read(vs2, 0, 1, e2))
    return false;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result = result < e1 ? result : e1;
    else
      errors++;

  if (errors == 0)
    if (not vecRegs_.write(vd, 0, group, result))
      errors++;

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vredmin_vs<int8_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vredmin_vs<int16_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vredmin_vs<int32_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vredmin_vs<int64_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vredmin_vs<Int128>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vredmin_vs<Int256>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vredmin_vs<Int512>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vredmin_vs<Int1024>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
bool
Hart<URV>::vredmaxu_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                       unsigned start, unsigned elems)
{
  unsigned errors = 0;

  ELEM_TYPE e2 = 0;
  if (vecRegs_.read(vs2, 0, 1, e2))
    return false;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result = result > e1 ? result : e1;
    else
      errors++;

  if (errors == 0)
    if (not vecRegs_.write(vd, 0, group, result))
      errors++;

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vredmaxu_vs<uint8_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vredmaxu_vs<uint16_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vredmaxu_vs<uint32_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vredmaxu_vs<uint64_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vredmaxu_vs<Uint128>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vredmaxu_vs<Uint256>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vredmaxu_vs<Uint512>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vredmaxu_vs<Uint1024>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;
    }
}


template <typename URV>
template<typename ELEM_TYPE>
bool
Hart<URV>::vredmax_vs(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                      unsigned start, unsigned elems)
{
  unsigned errors = 0;

  ELEM_TYPE e2 = 0;
  if (vecRegs_.read(vs2, 0, 1, e2))
    return false;
  
  ELEM_TYPE e1 = 0, result = e2;

  for (unsigned ix = start; ix < elems; ++ix)
    if (vecRegs_.read(vs1, ix, group, e1))
      result = result > e1 ? result : e1;
    else
      errors++;

  if (errors == 0)
    if (not vecRegs_.write(vd, 0, group, result))
      errors++;

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(), start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vredmax_vs<int8_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vredmax_vs<int16_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vredmax_vs<int32_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vredmax_vs<int64_t>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vredmax_vs<Int128>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vredmax_vs<Int256>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vredmax_vs<Int512>(vd, vs1, vs2, group, start, elems))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vredmax_vs<Int1024>(vd, vs1, vs2, group, start, elems))
        errors++;
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
bool
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned dist = vd > vs1 ? vd - vs1 : vs1 - vd;
  if (dist > group)
    {
      illegalInst(di);  // Source/dest vecs cannot overlap
      return;
    }

  URV amount = intRegs_.read(rs2);

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vslideup<uint8_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vslideup<uint16_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vslideup<uint32_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vslideup<uint64_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vslideup<Uint128>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vslideup<Uint256>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vslideup<Uint512>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vslideup<Uint1024>(vd, vs1, amount, group, start, elems, masked))
        errors++;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned dist = vd > vs1 ? vd - vs1 : vs1 - vd;
  if (dist > group)
    {
      illegalInst(di);  // Source/dest vecs cannot overlap
      return;
    }

  URV amount = imm;

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vslideup<uint8_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vslideup<uint16_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vslideup<uint32_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vslideup<uint64_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vslideup<Uint128>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vslideup<Uint256>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vslideup<Uint512>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vslideup<Uint1024>(vd, vs1, amount, group, start, elems, masked))
        errors++;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned dist = vd > vs1 ? vd - vs1 : vs1 - vd;
  if (dist > group)
    {
      illegalInst(di);  // Source/dest vecs cannot overlap
      return;
    }

  URV amount = 1,  replacement = intRegs_.read(rs2);

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vslideup<uint8_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      else
        vecRegs_.write(vd, 0, group, uint8_t(replacement));
      break;

    case ElementWidth::HalfWord:
      if (not vslideup<uint16_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      else
        vecRegs_.write(vd, 0, group, uint16_t(replacement));
      break;

    case ElementWidth::Word:
      if (not vslideup<uint32_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      else
        vecRegs_.write(vd, 0, group, uint32_t(replacement));
      break;

    case ElementWidth::DoubleWord:
      if (not vslideup<uint64_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      else
        vecRegs_.write(vd, 0, group, uint64_t(replacement));
      break;

    case ElementWidth::QuadWord:
      if (not vslideup<Uint128>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      else
        vecRegs_.write(vd, 0, group, Uint128(replacement));
      break;

    case ElementWidth::OctWord:
      if (not vslideup<Uint256>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      else
        vecRegs_.write(vd, 0, group, Uint256(replacement));
      break;

    case ElementWidth::HalfKbits:
      if (not vslideup<Uint512>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      else
        vecRegs_.write(vd, 0, group, Uint512(replacement));
      break;

    case ElementWidth::Kbits:
      if (not vslideup<Uint1024>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      else
        vecRegs_.write(vd, 0, group, Uint1024(replacement));
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned dist = vd > vs1 ? vd - vs1 : vs1 - vd;
  if (dist > group)
    {
      illegalInst(di);  // Source/dest vecs cannot overlap
      return;
    }

  URV amount = intRegs_.read(rs2);

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vslidedown<uint8_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vslidedown<uint16_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vslidedown<uint32_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vslidedown<uint64_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vslidedown<Uint128>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vslidedown<Uint256>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vslidedown<Uint512>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vslidedown<Uint1024>(vd, vs1, amount, group, start, elems, masked))
        errors++;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned dist = vd > vs1 ? vd - vs1 : vs1 - vd;
  if (dist > group)
    {
      illegalInst(di);  // Source/dest vecs cannot overlap
      return;
    }

  URV amount = imm;

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vslidedown<uint8_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vslidedown<uint16_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vslidedown<uint32_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vslidedown<uint64_t>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vslidedown<Uint128>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vslidedown<Uint256>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vslidedown<Uint512>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vslidedown<Uint1024>(vd, vs1, amount, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool
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

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vmul_vv<int8_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vmul_vv<int16_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vmul_vv<int32_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vmul_vv<int64_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vmul_vv<Int128>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vmul_vv<Int256>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vmul_vv<Int512>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vmul_vv<Int1024>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool
Hart<URV>::vmul_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  SRV srv2 = intRegs_.read(rs2);

  ELEM_TYPE val2 = srv2, e1 = 0, dest = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          dest = e1 * val2;
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vmul_vx<int8_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vmul_vx<int16_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vmul_vx<int32_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vmul_vx<int64_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vmul_vx<Int128>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vmul_vx<Int256>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vmul_vx<Int512>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      if (not vmul_vx<Int1024>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE, typename ELEM_TYPE_X2>
bool
Hart<URV>::vmulh_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  unsigned elemBits = 8*sizeof(ELEM_TYPE);
  if constexpr (std::is_same<Int128, ELEM_TYPE>::value)
    elemBits = 128;
  if constexpr (std::is_same<Int256, ELEM_TYPE>::value)
    elemBits = 256;
  if constexpr (std::is_same<Int512, ELEM_TYPE>::value)
    elemBits = 512;
  if constexpr (std::is_same<Int1024, ELEM_TYPE>::value)
    elemBits = 1024;

  ELEM_TYPE e1 = 0, e2 = 0, dest = 0;
  ELEM_TYPE_X2 temp = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
        {
          temp = e1;
          temp = temp * e2;
          temp = temp >> elemBits;
          dest = ELEM_TYPE(temp);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vmulh_vv<int8_t,int16_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vmulh_vv<int16_t,int32_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vmulh_vv<int32_t,int64_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vmulh_vv<int64_t,Int128>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vmulh_vv<Int128,Int256>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vmulh_vv<Int256,Int512>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vmulh_vv<Int512,Int1024>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      std::cerr << "vmulh_vv not yet supported for SEW=1024\n";
      errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE, typename ELEM_TYPE_X2>
bool
Hart<URV>::vmulh_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  unsigned elemBits = 8*sizeof(ELEM_TYPE);
  if constexpr (std::is_same<Int128, ELEM_TYPE>::value)
    elemBits = 128;
  if constexpr (std::is_same<Int256, ELEM_TYPE>::value)
    elemBits = 256;
  if constexpr (std::is_same<Int512, ELEM_TYPE>::value)
    elemBits = 512;
  if constexpr (std::is_same<Int1024, ELEM_TYPE>::value)
    elemBits = 1024;

  ELEM_TYPE e1 = 0, e2 = SRV(intRegs_.read(rs2)), dest = 0;
  ELEM_TYPE_X2 temp = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          temp = e1;
          temp = temp * e2;
          temp = temp >> elemBits;
          dest = ELEM_TYPE(temp);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vmulh_vx<int8_t,int16_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vmulh_vx<int16_t,int32_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vmulh_vx<int32_t,int64_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vmulh_vx<int64_t,Int128>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vmulh_vx<Int128,Int256>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vmulh_vx<Int256,Int512>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vmulh_vx<Int512,Int1024>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      std::cerr << "vmulh_vx not yet supported for SEW=1024\n";
      errors++;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vmulh_vv<uint8_t,uint16_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vmulh_vv<uint16_t,uint32_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vmulh_vv<uint32_t,uint64_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vmulh_vv<uint64_t,Uint128>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vmulh_vv<Uint128,Uint256>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vmulh_vv<Uint256,Uint512>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vmulh_vv<Uint512,Uint1024>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      std::cerr << "vmulhu_vv not yet supported for SEW=1024\n";
      errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE, typename ELEM_TYPE_X2>
bool
Hart<URV>::vmulhu_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  unsigned errors = 0;

  unsigned elemBits = 8*sizeof(ELEM_TYPE);
  if constexpr (std::is_same<Int128, ELEM_TYPE>::value)
    elemBits = 128;
  if constexpr (std::is_same<Int256, ELEM_TYPE>::value)
    elemBits = 256;
  if constexpr (std::is_same<Int512, ELEM_TYPE>::value)
    elemBits = 512;
  if constexpr (std::is_same<Int1024, ELEM_TYPE>::value)
    elemBits = 1024;

  ELEM_TYPE e1 = 0, e2 = intRegs_.read(rs2), dest = 0;
  ELEM_TYPE_X2 temp = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          temp = e1;
          temp = temp * e2;
          temp = temp >> elemBits;
          dest = ELEM_TYPE(temp);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vmulhu_vx<uint8_t,uint16_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vmulhu_vx<uint16_t,uint32_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vmulhu_vx<uint32_t,uint64_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vmulhu_vx<uint64_t,Uint128>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vmulhu_vx<Uint128,Uint256>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vmulhu_vx<Uint256,Uint512>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vmulhu_vx<Uint512,Uint1024>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      std::cerr << "vmulhu_vx not yet supported for SEW=1024\n";
      errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE, typename ELEM_TYPE_X2>
bool
Hart<URV>::vmulhsu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                      unsigned start, unsigned elems, bool masked)
{
  typedef typename std::make_unsigned<ELEM_TYPE>::type U_ELEM_TYPE;
  unsigned errors = 0;

  unsigned elemBits = 8*sizeof(ELEM_TYPE);
  if constexpr (std::is_same<Int128, ELEM_TYPE>::value)
    elemBits = 128;
  if constexpr (std::is_same<Int256, ELEM_TYPE>::value)
    elemBits = 256;
  if constexpr (std::is_same<Int512, ELEM_TYPE>::value)
    elemBits = 512;
  if constexpr (std::is_same<Int1024, ELEM_TYPE>::value)
    elemBits = 1024;

  ELEM_TYPE e1 = 0, dest = 0;
  U_ELEM_TYPE e2 = 0;
  ELEM_TYPE_X2 temp = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1) and
          vecRegs_.read(vs2, ix, group, e2))
        {
          temp = e1;
          temp = temp * e2;
          temp = temp >> elemBits;
          dest = ELEM_TYPE(temp);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vmulhsu_vv<int8_t,int16_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vmulhsu_vv<int16_t,int32_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vmulhsu_vv<int32_t,int64_t>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vmulhsu_vv<int64_t,Int128>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vmulhsu_vv<Int128,Int256>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vmulhsu_vv<Int256,Int512>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vmulhsu_vv<Int512,Int1024>(vd, vs1, vs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      std::cerr << "vmulhsu_vv not yet supported for SEW=1024\n";
      errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE, typename ELEM_TYPE_X2>
bool
Hart<URV>::vmulhsu_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                      unsigned start, unsigned elems, bool masked)
{
  typedef typename std::make_unsigned<ELEM_TYPE>::type U_ELEM_TYPE;
  unsigned errors = 0;

  unsigned elemBits = 8*sizeof(ELEM_TYPE);
  if constexpr (std::is_same<Int128, ELEM_TYPE>::value)
    elemBits = 128;
  if constexpr (std::is_same<Int256, ELEM_TYPE>::value)
    elemBits = 256;
  if constexpr (std::is_same<Int512, ELEM_TYPE>::value)
    elemBits = 512;
  if constexpr (std::is_same<Int1024, ELEM_TYPE>::value)
    elemBits = 1024;

  ELEM_TYPE e1 = 0, dest = 0;
  U_ELEM_TYPE e2 = intRegs_.read(rs2);
  ELEM_TYPE_X2 temp = 0;

  for (unsigned ix = start; ix < elems; ++ix)
    {
      if (masked and not vecRegs_.isActive(0, ix))
        continue;

      if (vecRegs_.read(vs1, ix, group, e1))
        {
          temp = e1;
          temp = temp * e2;
          temp = temp >> elemBits;
          dest = ELEM_TYPE(temp);
          if (not vecRegs_.write(vd, ix, group, dest))
            errors++;
        }
      else
        errors++;
    }

  return errors == 0;
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

  unsigned group = vecRegs_.groupMultiplier(),  start = vecRegs_.startIndex();
  unsigned elems = vecRegs_.elemCount();
  ElementWidth sew = vecRegs_.elemWidth();

  unsigned errors = 0;

  switch (sew)
    {
    case ElementWidth::Byte:
      if (not vmulhsu_vx<int8_t,int16_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfWord:
      if (not vmulhsu_vx<int16_t,int32_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Word:
      if (not vmulhsu_vx<int32_t,int64_t>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::DoubleWord:
      if (not vmulhsu_vx<int64_t,Int128>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::QuadWord:
      if (not vmulhsu_vx<Int128,Int256>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::OctWord:
      if (not vmulhsu_vx<Int256,Int512>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::HalfKbits:
      if (not vmulhsu_vx<Int512,Int1024>(vd, vs1, rs2, group, start, elems, masked))
        errors++;
      break;

    case ElementWidth::Kbits:
      std::cerr << "vmulhsu_vx not yet supported for SEW=1024\n";
      errors++;
      break;
    }
}


template <typename URV>
template <typename ELEM_TYPE>
bool
Hart<URV>::vdivu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  return false;
}


template <typename URV>
void
Hart<URV>::execVdivu_vv(const DecodedInst* di)
{
}


template <typename URV>
template <typename ELEM_TYPE>
bool
Hart<URV>::vdivu_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  return false;
}


template <typename URV>
void
Hart<URV>::execVdivu_vx(const DecodedInst* di)
{
}


template <typename URV>
template <typename ELEM_TYPE>
bool
Hart<URV>::vdiv_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  return false;
}


template <typename URV>
void
Hart<URV>::execVdiv_vv(const DecodedInst* di)
{
}


template <typename URV>
template <typename ELEM_TYPE>
bool
Hart<URV>::vdiv_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  return false;
}


template <typename URV>
void
Hart<URV>::execVdiv_vx(const DecodedInst* di)
{
}


template <typename URV>
template <typename ELEM_TYPE>
bool
Hart<URV>::vremu_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  return false;
}


template <typename URV>
void
Hart<URV>::execVremu_vv(const DecodedInst* di)
{
}


template <typename URV>
template <typename ELEM_TYPE>
bool
Hart<URV>::vremu_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                    unsigned start, unsigned elems, bool masked)
{
  return false;
}


template <typename URV>
void
Hart<URV>::execVremu_vx(const DecodedInst* di)
{
}


template <typename URV>
template <typename ELEM_TYPE>
bool
Hart<URV>::vrem_vv(unsigned vd, unsigned vs1, unsigned vs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  return false;
}


template <typename URV>
void
Hart<URV>::execVrem_vv(const DecodedInst* di)
{
}


template <typename URV>
template <typename ELEM_TYPE>
bool
Hart<URV>::vrem_vx(unsigned vd, unsigned vs1, unsigned rs2, unsigned group,
                   unsigned start, unsigned elems, bool masked)
{
  return false;
}


template <typename URV>
void
Hart<URV>::execVrem_vx(const DecodedInst* di)
{
}


template class WdRiscv::Hart<uint32_t>;
template class WdRiscv::Hart<uint64_t>;
