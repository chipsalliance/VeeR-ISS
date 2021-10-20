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

#include "FpRegs.hpp"
#include "softfloat-util.hpp"

using namespace WdRiscv;


float
Float16::toFloat() const
{
#ifdef SOFT_FLOAT

  auto sf16 = nativeToSoft(*this);
  auto sf32 = f16_to_f32(sf16);
  return softToNative(sf32);

#else

  bool sign = (i16 >> 15) & 1;
  if (isInf())
    {
      float x = std::numeric_limits<float>::infinity();
      return sign? -x : x;
    }

  if (isSnan())
    {
      float x = std::numeric_limits<float>::signaling_NaN();
      return sign? -x : x;
    }

  if (isNan())
    {
      float x = std::numeric_limits<float>::quiet_NaN();
      return sign? -x : x;
    }

  if (isZero())
    return sign? -0.0f : 0.0f;

  if (isSubnormal())
    {
      // Subnormal in half precision would be normal in float.
      // Renormalize.
      uint32_t sig = sigBits();
      assert(sig != 0);
      uint32_t exp = expBits();
      unsigned mssb = __builtin_clz(sig);  // Most sig set bit
      assert(mssb <= 9);
      unsigned shift = 10 - mssb;
      sig = sig & ~(uint32_t(1) << shift);  // Clear most sig bit
      sig = sig << shift;
      exp = exp - shift;
      exp = exp - 15 + 127;  // Update bias
      uint32_t val = sign? 1 : 0;
      val = (val << 31) | (exp << 23) | sig;
      Uint32FloatUnion uf{val};
      return uf.f;
    }

  // Normalized number. Update exponent for float bias.
  uint32_t sig = sigBits();
  uint32_t exp = expBits();
  exp = exp - 15 + 127;
  uint32_t val = sign? 1 : 0;
  val = (val << 31) | (exp << 23) | sig;
  Uint32FloatUnion uf{val};
  return uf.f;

#endif
}


Float16
Float16::fromFloat(float val)
{
#ifdef SOFT_FLOAT

  auto sf32 = nativeToSoft(val);
  auto sf16 = f32_to_f16(sf32);
  return softToNative(sf16);

#else

  bool sign = std::signbit(val);
  if (std::isinf(val))
    {
      Float16 x = Float16::infinity();
      return sign? -x : x;
    }

  if (std::isnan(val))
    {
      Float16 x = WdRiscv::isSnan(val)? Float16::signalingNan() : Float16::quietNan();
      return sign? -x : x;
    }

  if (val == 0 or not std::isnormal(val))
    return sign? -Float16{} : Float16{};

  // Normalized number. Update exponent for float16 bias.
  Uint32FloatUnion uf{val};

  uint32_t sig = (uf.u << 9) >> 9;
  int exp = ((uf.u) >> 23) & 0xff;
  exp = exp - 127 + 15;
  if (exp < -10)
    return sign? -Float16{} : Float16{};
  if (exp < 0)
    {
      assert(0);
    }
  if (exp >= 0x1f)
    return sign? -Float16::infinity() : Float16::infinity();

  uint16_t res = sign? 1 : 0;
  res = (res << 15) | uint16_t(exp << 10) | uint16_t(sig >> 13);

  return Float16::fromBits(res);

#endif
}


FpRegs::FpRegs(unsigned regCount)
  : regs_(regCount, 0)
{
  numberToName_.resize(32);

  for (unsigned ix = 0; ix < 32; ++ix)
    {
      std::string name = "f" + std::to_string(ix);
      nameToNumber_[name] = FpRegNumber(ix);
      numberToName_.at(ix) = name;
    }

  numberToAbiName_ = { "ft0", "ft1", "ft2", "ft3", "ft4", "ft5", "ft6", "ft7",
		       "fs0", "fs1", "fa0", "fa1", "fa2", "fa3", "fa4", "fa5",
		       "fa6", "fa7", "fs2", "fs3", "fs4", "fs5", "fs6", "fs7",
		       "fs8", "fs9", "fs10", "fs11", "ft8", "ft9", "ft10", "ft11" };

  for (unsigned ix = 0; ix < 32; ++ix)
    {
      std::string abiName = numberToAbiName_.at(ix);
      nameToNumber_[abiName] = FpRegNumber(ix);
    }
}


bool
FpRegs::findReg(const std::string& name, unsigned& ix) const
{
  const auto iter = nameToNumber_.find(name);
  if (iter == nameToNumber_.end())
    return false;

  ix = iter->second;
  return true;
}


void
FpRegs::reset(bool hasHalf, bool hasSingle, bool hasDouble)
{
  hasHalf_ = hasHalf_;
  hasSingle_ = hasSingle_;
  hasDouble_ = hasDouble_;

  if (hasDouble)
    {
      for (auto& reg : regs_)
	reg = 0;
    }
  else if (hasSingle)
    {
      // F extension present without D. Reset to NAN-boxed
      // single-precision zeros.
      for (size_t i = 0; i < regs_.size(); ++i)
	writeSingle(i, 0);
    }
  else if (hasHalf)
    {
      // F16 extension present without F or D. Reset to NAN-boxed
      // half-precision zeros.
      for (size_t i = 0; i < regs_.size(); ++i)
	writeHalf(i, Float16());
    }

  clearLastWrittenReg();
}
