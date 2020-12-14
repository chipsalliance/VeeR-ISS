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


// This file contains the methods of the class Hart supporting
// the RISCV floating point instructions.

#include <cfenv>
#include <cmath>
#include <emmintrin.h>
#include <array>

#ifdef SOFT_FLOAT
extern "C" {
#include <softfloat.h>
}
#endif

#include "Hart.hpp"
#include "instforms.hpp"


using namespace WdRiscv;


/// Unsigned-float union: reinterpret bits as uint32_t or float
union Uint32FloatUnion
{
  Uint32FloatUnion(uint32_t u) : u(u)
  { }

  Uint32FloatUnion(float f) : f(f)
  { }

  uint32_t u = 0;
  float f;
};


/// Unsigned-float union: reinterpret bits as uint64_t or double
union Uint64DoubleUnion
{
  Uint64DoubleUnion(uint64_t u) : u(u)
  { }

  Uint64DoubleUnion(double d) : d(d)
  { }

  uint64_t u = 0;
  double d;
};


template <typename URV>
void
Hart<URV>::resetFloat()
{
  // Enable FP if f/d extension and linux/newlib.
  if ((isRvf() or isRvd()) and (newlib_ or linux_))
    {
      URV val = csRegs_.peekMstatus();
      MstatusFields<URV> fields(val);
      fields.bits_.FS = unsigned(FpFs::Initial);
      csRegs_.write(CsrNumber::MSTATUS, PrivilegeMode::Machine, fields.value_);
    }

  if (isRvf() or isRvd())
    fpRegs_.reset(isRvd());

  #ifdef SOFT_FLOAT

  softfloat_exceptionFlags = 0;
  softfloat_detectTininess = softfloat_tininess_afterRounding;

  #endif
}


template <typename URV>
void
Hart<URV>::setFcsrFlags(FpFlags flag)
{
  auto prev = getFpFlags();
  auto val = prev | unsigned(flag);
  if (val != prev)
    {
      setFpFlags(val);
      recordCsrWrite(CsrNumber::FCSR);
    }
}


template <typename URV>
inline
RoundingMode
Hart<URV>::effectiveRoundingMode(RoundingMode instMode)
{
  if (instMode != RoundingMode::Dynamic)
    return instMode;

  return getFpRoundingMode();
}


/// Return the exponent bits of the given floating point value.
unsigned
spExponentBits(float sp)
{
  Uint32FloatUnion uf(sp);
  return (uf.u << 1) >> 24;
}


/// Return the exponent bits of the given double precision value.
unsigned
dpExponentBits(double dp)
{
  Uint64DoubleUnion ud(dp);
  return (ud.u << 1) >> 53;
}


#ifdef FAST_SLOPPY

template <typename URV>
inline
void
Hart<URV>::updateAccruedFpBits(float, bool)
{
}


template <typename URV>
inline
void
Hart<URV>::updateAccruedFpBits(double, bool)
{
}


template <typename URV>
inline
void
Hart<URV>::markFsDirty()
{
}


#else

template <typename URV>
inline
void
Hart<URV>::updateAccruedFpBits(float res, bool invalid)
{
  URV val = getFpFlags();
  URV prev = val;

#ifdef SOFT_FLOAT
  res = res; // Passify compiler.
  int flags = softfloat_exceptionFlags;
  if (flags)
    {
      if (flags & softfloat_flag_inexact)   val |= URV(FpFlags::Inexact);
      if (flags & softfloat_flag_underflow) val |= URV(FpFlags::Underflow);
      if (flags & softfloat_flag_overflow)  val |= URV(FpFlags::Overflow);
      if (flags & softfloat_flag_infinite)  val |= URV(FpFlags::DivByZero);
      if (flags & softfloat_flag_invalid)   val |= URV(FpFlags::Invalid);
    }      
#else
  int flags = fetestexcept(FE_ALL_EXCEPT);
  if (flags)
    {
      if (flags & FE_INEXACT)
        val |= URV(FpFlags::Inexact);
      if ((flags & FE_UNDERFLOW) and (spExponentBits(res) == 0))
        val |= URV(FpFlags::Underflow);
      if (flags & FE_OVERFLOW)
        val |= URV(FpFlags::Overflow);
      if (flags & FE_DIVBYZERO)
        val |= URV(FpFlags::DivByZero);
      if (flags & FE_INVALID)
        val |= URV(FpFlags::Invalid);
    }
#endif

  if (invalid)
    val |= URV(FpFlags::Invalid);

  if (val != prev)
    {
      setFpFlags(val);
      recordCsrWrite(CsrNumber::FCSR);
    }
}


template <typename URV>
inline
void
Hart<URV>::updateAccruedFpBits(double res, bool invalid)
{
  URV val = getFpFlags();
  URV prev = val;

#ifdef SOFT_FLOAT
  res = res; // Passify compiler.
  int flags = softfloat_exceptionFlags;
  if (flags)
    {
      if (flags & softfloat_flag_inexact)   val |= URV(FpFlags::Inexact);
      if (flags & softfloat_flag_underflow) val |= URV(FpFlags::Underflow);
      if (flags & softfloat_flag_overflow)  val |= URV(FpFlags::Overflow);
      if (flags & softfloat_flag_infinite)  val |= URV(FpFlags::DivByZero);
      if (flags & softfloat_flag_invalid)   val |= URV(FpFlags::Invalid);
    }      
#else
  int flags = fetestexcept(FE_ALL_EXCEPT);
  if (flags)
    {
      if (flags & FE_INEXACT)
        val |= URV(FpFlags::Inexact);
      if ((flags & FE_UNDERFLOW) and (spExponentBits(res) == 0))
        val |= URV(FpFlags::Underflow);
      if (flags & FE_OVERFLOW)
        val |= URV(FpFlags::Overflow);
      if (flags & FE_DIVBYZERO)
        val |= URV(FpFlags::DivByZero);
      if (flags & FE_INVALID)
        val |= URV(FpFlags::Invalid);
    }
#endif

  if (invalid)
    val |= URV(FpFlags::Invalid);

  if (val != prev)
    {
      setFpFlags(val);
      recordCsrWrite(CsrNumber::FCSR);
    }
}


template <typename URV>
inline
void
Hart<URV>::markFsDirty()
{
  if (mstatusFs_ == FpFs::Dirty)
    return;

  URV val = csRegs_.peekMstatus();
  MstatusFields<URV> fields(val);
  fields.bits_.FS = unsigned(FpFs::Dirty);

  csRegs_.poke(CsrNumber::MSTATUS, fields.value_);

  URV newVal = csRegs_.peekMstatus();
  if (val != newVal)
    recordCsrWrite(CsrNumber::MSTATUS);

  updateCachedMstatusFields();
}

#endif


#ifdef SOFT_FLOAT
/// Map a RISCV rounding mode to a soft-float constant.
static std::array<int, 5> riscvRoungingModeToSoftFloat =
  {
   softfloat_round_near_even,   // NearsetEven
   softfloat_round_minMag,      // Zero
   softfloat_round_min,         // Down
   softfloat_round_max,         // Up
   softfloat_round_near_maxMag  // NearestMax
  };
#endif


#ifdef SOFT_FLOAT

static
inline
int
mapRiscvRoundingModeToSoftFloat(RoundingMode mode)
{
  uint32_t ix = uint32_t(mode);
  return riscvRoungingModeToSoftFloat.at(ix);
}


static
inline
int
setSimulatorRoundingMode(RoundingMode mode)
{
  int previous = softfloat_roundingMode;
  int next = mapRiscvRoundingModeToSoftFloat(mode);

  softfloat_roundingMode = next;

  return previous;
}


/// Clear the floating point flags in the machine running this
/// simulator. Do nothing in the simuated RISCV machine.
static
inline
void
clearSimulatorFpFlags()
{
  softfloat_exceptionFlags = 0;
}

#else

/// Map a RISCV rounding mode to an fetsetround constant.
static std::array<int, 5> riscvRoungingModeToFe =
  {
   FE_TONEAREST,  // NearsetEven
   FE_TOWARDZERO, // Zero
   FE_DOWNWARD,   // Down
   FE_UPWARD,     // Up
   FE_TONEAREST   // NearestMax
  };


static
inline
int
mapRiscvRoundingModeToFe(RoundingMode mode)
{
  uint32_t ix = uint32_t(mode);
  return riscvRoungingModeToFe.at(ix);
}
  

static
inline
int
setSimulatorRoundingMode(RoundingMode mode)
{
  int previous = std::fegetround();
  int next = mapRiscvRoundingModeToFe(mode);

  if (next != previous)
    std::fesetround(next);

  return previous;
}


/// Clear the floating point flags in the machine running this
/// simulator. Do nothing in the simuated RISCV machine.
static
inline
void
clearSimulatorFpFlags()
{
  uint32_t val = _mm_getcsr();
  val &= ~uint32_t(0x3f);
  _mm_setcsr(val);
  // std::feclearexcept(FE_ALL_EXCEPT);
}

#endif


template <typename URV>
inline
bool
Hart<URV>::checkRoundingModeSp(const DecodedInst* di)
{
if (not isFpLegal())
    {
      illegalInst(di);
      return false;
    }

  RoundingMode riscvMode = effectiveRoundingMode(di->roundingMode());
  if (riscvMode >= RoundingMode::Invalid1)
    {
      illegalInst(di);
      return false;
    }

  clearSimulatorFpFlags();
  setSimulatorRoundingMode(riscvMode);
  return true;
}


template <typename URV>
inline
bool
Hart<URV>::checkRoundingModeDp(const DecodedInst* di)
{
  if (not isDpLegal())
    {
      illegalInst(di);
      return false;
    }

  RoundingMode riscvMode = effectiveRoundingMode(di->roundingMode());
  if (riscvMode >= RoundingMode::Invalid1)
    {
      illegalInst(di);
      return false;
    }

  clearSimulatorFpFlags();
  setSimulatorRoundingMode(riscvMode);
  return true;
}


template <typename URV>
void
Hart<URV>::execFlw(const DecodedInst* di)
{
  if (not isFpLegal())
    {
      illegalInst(di);
      return;
    }

  uint32_t rd = di->op0(), rs1 = di->op1();
  SRV imm = di->op2As<SRV>();

  URV base = intRegs_.read(rs1);
  URV virtAddr = base + imm;

  ldStAddr_ = virtAddr;   // For reporting load addr in trace-mode.
  ldStAddrValid_ = true;  // For reporting load addr in trace-mode.
  uint64_t addr = virtAddr;
  unsigned ldSize = 4;

  auto secCause = SecondaryCause::NONE;
  auto cause = ExceptionCause::NONE;

#ifndef FAST_SLOPPY
  if (loadQueueEnabled_)
    removeFromLoadQueue(rs1, false);

  if (hasActiveTrigger())
    {
      if (ldStAddrTriggerHit(virtAddr, TriggerTiming::Before, true /*isLoad*/,
			     privMode_, isInterruptEnabled()))
	triggerTripped_ = true;
      if (triggerTripped_)
	return;
    }

  cause = determineLoadException(rs1, base, addr, ldSize, secCause);
  if (cause != ExceptionCause::NONE)
    {
      if (not triggerTripped_)
        initiateLoadException(cause, virtAddr, secCause);
      return;
    }
#endif

  uint32_t word = 0;
  if (memory_.read(addr, word))
    {
      if (loadQueueEnabled_)
        {
          uint64_t prevRdVal = 0;
          peekFpReg(rd, prevRdVal);
          putInLoadQueue(ldSize, addr, rd, prevRdVal, false /*wide*/, true /*fp*/);
        }
      Uint32FloatUnion ufu(word);
      fpRegs_.writeSingle(rd, ufu.f);
      markFsDirty();
      return;
    }

  cause = ExceptionCause::LOAD_ACC_FAULT;
  secCause = SecondaryCause::LOAD_ACC_MEM_PROTECTION;
  if (isAddrMemMapped(addr))
    secCause = SecondaryCause::LOAD_ACC_PIC;

  initiateLoadException(cause, virtAddr, secCause);
}


template <typename URV>
void
Hart<URV>::execFsw(const DecodedInst* di)
{
if (not isFpLegal())
    {
      illegalInst(di);
      return;
    }

  uint32_t rs1 = di->op1(), rs2 = di->op0();
  SRV imm = di->op2As<SRV>();

  URV base = intRegs_.read(rs1);
  URV addr = base + imm;

  // This operation does not check for proper NAN boxing. We read raw bits.
  uint64_t val = fpRegs_.readBitsRaw(rs2);

  store<uint32_t>(rs1, base, addr, uint32_t(val));
}


#ifdef SOFT_FLOAT

/// Convert softfloat float32_t to a native float.
inline float
f32ToFloat(float32_t f32)
{
  Uint32FloatUnion tmp(f32.v);
  return tmp.f;
}


/// Convert softfloat float64_t to a native double.
inline double
f64ToDouble(float64_t f64)
{
  Uint64DoubleUnion tmp(f64.v);
  return tmp.d;
}


/// Convert a native float to a softfloat float32_t
inline float32_t
floatToF32(float x)
{
  Uint32FloatUnion tmp(x);
  return float32_t{tmp.u};
}


/// Convert a native double to a softfloat float64_t
inline float64_t
doubleToF64(double x)
{
  Uint64DoubleUnion tmp(x);
  return float64_t{tmp.u};
}

#endif


/// Use fused mutiply-add to perform x*y + z.
/// Set invalid to true if x and y are zero and infinity or
/// vice versa since RISCV consider that as an invalid operation.
static
float
fusedMultiplyAdd(float x, float y, float z, bool& invalid)
{
#ifndef SOFT_FLOAT
  #ifdef __FP_FAST_FMA
  float res = x*y + z;
  #else
  float res = std::fma(x, y, z);
  #endif
#else
  float32_t tmp = f32_mulAdd(floatToF32(x), floatToF32(y), floatToF32(z));
  float res = f32ToFloat(tmp);
#endif

  invalid = (std::isinf(x) and y == 0) or (x == 0 and std::isinf(y));
  return res;
}


/// Use fused mutiply-add to perform x*y + z.
static
double
fusedMultiplyAdd(double x, double y, double z, bool& invalid)
{
#ifndef SOFT_FLOAT
  #ifdef __FP_FAST_FMA
  double res = x*y + z;
  #else
  double res = std::fma(x, y, z);
  #endif
#else
  float64_t tmp = f64_mulAdd(doubleToF64(x), doubleToF64(y), doubleToF64(z));
  double res = f64ToDouble(tmp);
#endif

  invalid = (std::isinf(x) and y == 0) or (x == 0 and std::isinf(y));
  return res;
}


template <typename URV>
void
Hart<URV>::execFmadd_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float f3 = fpRegs_.readSingle(di->op3());

  bool invalid = false;
  float res = fusedMultiplyAdd(f1, f2, f3, invalid);
  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  updateAccruedFpBits(res, invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmsub_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float f3 = -fpRegs_.readSingle(di->op3());

  bool invalid = false;
  float res = fusedMultiplyAdd(f1, f2, f3, invalid);
  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  updateAccruedFpBits(res, invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFnmsub_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = -fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float f3 = fpRegs_.readSingle(di->op3());

  bool invalid = false;
  float res = fusedMultiplyAdd(f1, f2, f3, invalid);
  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  updateAccruedFpBits(res, invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFnmadd_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  // we want -(f[op1] * f[op2]) - f[op3]

  float f1 = -fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float f3 = -fpRegs_.readSingle(di->op3());

  bool invalid = false;
  float res = fusedMultiplyAdd(f1, f2, f3, invalid);
  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  updateAccruedFpBits(res, invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFadd_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());

#ifdef SOFT_FLOAT
  float res = f32ToFloat(f32_add(floatToF32(f1), floatToF32(f2)));
#else
  float res = f1 + f2;
#endif

  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsub_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());

#ifdef SOFT_FLOAT
  float res = f32ToFloat(f32_sub(floatToF32(f1), floatToF32(f2)));
#else
  float res = f1 - f2;
#endif

  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmul_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());

#ifdef SOFT_FLOAT
  float res = f32ToFloat(f32_mul(floatToF32(f1), floatToF32(f2)));
#else
  float res = f1 * f2;
#endif

  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFdiv_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());

#ifdef SOFT_FLOAT
  float res = f32ToFloat(f32_div(floatToF32(f1), floatToF32(f2)));
#else
  float res = f1 / f2;
#endif

  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsqrt_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());

#ifdef SOFT_FLOAT
  float res = f32ToFloat(f32_sqrt(floatToF32(f1)));
#else
  float res = std::sqrt(f1);
#endif

  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsgnj_s(const DecodedInst* di)
{
  if (not isFpLegal())
    {
      illegalInst(di);
      return;
    }

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float res = std::copysignf(f1, f2);  // Magnitude of rs1 and sign of rs2
  fpRegs_.writeSingle(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsgnjn_s(const DecodedInst* di)
{
if (not isFpLegal())
    {
      illegalInst(di);
      return;
    }

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());
  float res = std::copysignf(f1, f2);  // Magnitude of rs1 and sign of rs2
  res = -res;  // Magnitude of rs1 and negative the sign of rs2
  fpRegs_.writeSingle(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsgnjx_s(const DecodedInst* di)
{
if (not isFpLegal())
    {
      illegalInst(di);
      return;
    }

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());

  int sign1 = (std::signbit(f1) == 0) ? 0 : 1;
  int sign2 = (std::signbit(f2) == 0) ? 0 : 1;
  int sign = sign1 ^ sign2;

  float x = sign? -1 : 1;

  float res = std::copysignf(f1, x);  // Magnitude of rs1 and sign of x
  fpRegs_.writeSingle(di->op0(), res);

  markFsDirty();
}


/// Return true if given float is a signaling not-a-number.
static
bool
issnan(float f)
{
  if (std::isnan(f))
    {
      Uint32FloatUnion ufu(f);

      // Most sig bit of significant must be zero.
      return ((ufu.u >> 22) & 1) == 0;
    }
  return false;
}


/// Return true if given double is a signaling not-a-number.
static
bool
issnan(double d)
{
  if (std::isnan(d))
    {
      Uint64DoubleUnion udu(d);

      // Most sig bit of significant must be zero.
      return ((udu.u >> 51) & 1) == 0;
    }
  return false;
}


template <typename URV>
void
Hart<URV>::execFmin_s(const DecodedInst* di)
{
if (not isFpLegal())
    {
      illegalInst(di);
      return;
    }

  float in1 = fpRegs_.readSingle(di->op1());
  float in2 = fpRegs_.readSingle(di->op2());
  float res = 0;

  bool isNan1 = std::isnan(in1), isNan2 = std::isnan(in2);
  if (isNan1 and isNan2)
    res = std::numeric_limits<float>::quiet_NaN();
  else if (isNan1)
    res = in2;
  else if (isNan2)
    res = in1;
  else
    res = std::fminf(in1, in2);

  if (issnan(in1) or issnan(in2))
    setFcsrFlags(FpFlags::Invalid);
  else if (std::signbit(in1) != std::signbit(in2) and in1 == in2)
    res = std::copysign(res, -1.0F);  // Make sure min(-0, +0) is -0.

  fpRegs_.writeSingle(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmax_s(const DecodedInst* di)
{
if (not isFpLegal())
    {
      illegalInst(di);
      return;
    }

  float in1 = fpRegs_.readSingle(di->op1());
  float in2 = fpRegs_.readSingle(di->op2());
  float res = 0;

  bool isNan1 = std::isnan(in1), isNan2 = std::isnan(in2);
  if (isNan1 and isNan2)
    res = std::numeric_limits<float>::quiet_NaN();
  else if (isNan1)
    res = in2;
  else if (isNan2)
    res = in1;
  else
    res = std::fmaxf(in1, in2);

  if (issnan(in1) or issnan(in2))
    setFcsrFlags(FpFlags::Invalid);
  else if (std::signbit(in1) != std::signbit(in2) and in1 == in2)
    res = std::copysign(res, 1.0F);  // Make sure max(-0, +0) is +0.

  fpRegs_.writeSingle(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_w_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  SRV result = 0;
  bool valid = false;

#ifdef SOFT_FLOAT
  result = f32_to_i32(floatToF32(f1), softfloat_roundingMode, true);
  valid = true;  // We get invalid from softfloat library.
#else

  int32_t minInt = int32_t(1) << 31;
  int32_t maxInt = (~uint32_t(0)) >> 1;

  unsigned signBit = std::signbit(f1);
  if (std::isinf(f1))
    result = signBit ? minInt : maxInt;
  else if (std::isnan(f1))
    result = maxInt;
  else
    {
      float near = std::nearbyint(f1);
      if (near >= float(maxInt))
	result = maxInt;
      else if (near < float(minInt))
	result = SRV(minInt);
      else
	{
	  valid = true;
          result = int32_t(std::lrintf(f1));
	}
    }

#endif

  intRegs_.write(di->op0(), result);

  updateAccruedFpBits(0.0, not valid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_wu_s(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  SRV result = 0;

#ifdef SOFT_FLOAT

  // In 64-bit mode, we sign extend the result to 64-bits.
  result = SRV(int32_t(f32_to_ui32(floatToF32(f1), softfloat_roundingMode, true)));
  if (not std::isnan(f1) and f1 < 0)
    {
      softfloat_exceptionFlags &= ~softfloat_flag_inexact;  // Should not set inexact.
      softfloat_exceptionFlags |= softfloat_flag_invalid;
    }
  updateAccruedFpBits(0.0f, false);

#else

  bool valid = false;
  bool exact = true;

  uint32_t maxUint32 = ~uint32_t(0);
  if (std::isnan(f1))
    {
      result = ~URV(0);
    }
  else if (std::signbit(f1) and f1 != 0)
    {
      result = 0;
    }
  else
    {
      double near = std::nearbyint(f1);
      if (near > double(maxUint32))
        {
          result = ~URV(0);
        }
      else if (near < 0)
        {
          result = 0;
        }
      else if (near == 0)
        {
          result = 0;
          valid = true;
          exact = near == f1;
        }
      else
        {
          result = SRV(int32_t(std::lrint(f1)));
          valid = true;
          exact = near == f1;
        }
    }

  if (not valid)
    setFcsrFlags(FpFlags::Invalid);
  if (not exact)
    setFcsrFlags(FpFlags::Inexact);

#endif

  intRegs_.write(di->op0(), result);
  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmv_x_w(const DecodedInst* di)
{
if (not isFpLegal())
    {
      illegalInst(di);
      return;
    }

  // This operation does not check for proper NAN boxing. We read raw bits.
  uint64_t v1 = fpRegs_.readBitsRaw(di->op1());
  int32_t s1 = v1;  // Keep lower 32 bits

  SRV value = SRV(s1); // Sign extend.

  intRegs_.write(di->op0(), value);
}

 
template <typename URV>
void
Hart<URV>::execFeq_s(const DecodedInst* di)
{
if (not isFpLegal())
    {
      illegalInst(di);
      return;
    }

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());

  URV res = 0;

  if (std::isnan(f1) or std::isnan(f2))
    {
      if (issnan(f1) or issnan(f2))
	setFcsrFlags(FpFlags::Invalid);
    }
  else
    res = (f1 == f2)? 1 : 0;

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFlt_s(const DecodedInst* di)
{
if (not isFpLegal())
    {
      illegalInst(di);
      return;
    }

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());

  URV res = 0;

  if (std::isnan(f1) or std::isnan(f2))
    setFcsrFlags(FpFlags::Invalid);
  else
    res = (f1 < f2)? 1 : 0;
    
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFle_s(const DecodedInst* di)
{
if (not isFpLegal())
    {
      illegalInst(di);
      return;
    }

  float f1 = fpRegs_.readSingle(di->op1());
  float f2 = fpRegs_.readSingle(di->op2());

  URV res = 0;

  if (std::isnan(f1) or std::isnan(f2))
    setFcsrFlags(FpFlags::Invalid);
  else
    res = (f1 <= f2)? 1 : 0;

  intRegs_.write(di->op0(), res);
}


bool
mostSignificantFractionBit(float x)
{
  Uint32FloatUnion ufu(x);
  return (ufu.u >> 22) & 1;
}


bool
mostSignificantFractionBit(double x)
{
  Uint64DoubleUnion udu(x);
  return (udu.u >> 51) & 1;
}



template <typename URV>
void
Hart<URV>::execFclass_s(const DecodedInst* di)
{
if (not isFpLegal())
    {
      illegalInst(di);
      return;
    }

  float f1 = fpRegs_.readSingle(di->op1());
  URV result = 0;

  bool pos = not std::signbit(f1);
  int type = std::fpclassify(f1);

  if (type == FP_INFINITE)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosInfinity);
      else
	result |= URV(FpClassifyMasks::NegInfinity);
    }
  else if (type == FP_NORMAL)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosNormal);
      else
	result |= URV(FpClassifyMasks::NegNormal);
    }
  else if (type == FP_SUBNORMAL)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosSubnormal);
      else
	result |= URV(FpClassifyMasks::NegSubnormal);
    }
  else if (type == FP_ZERO)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosZero);
      else
	result |= URV(FpClassifyMasks::NegZero);
    }
  else if (type == FP_NAN)
    {
      bool quiet = mostSignificantFractionBit(f1);
      if (quiet)
	result |= URV(FpClassifyMasks::QuietNan);
      else
	result |= URV(FpClassifyMasks::SignalingNan);
    }

  intRegs_.write(di->op0(), result);
}


template <typename URV>
void
Hart<URV>::execFcvt_s_w(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  int32_t i1 = intRegs_.read(di->op1());

#ifdef SOFT_FLOAT
  float res = f32ToFloat(i32_to_f32(i1));
#else
  float res = float(i1);
#endif

  fpRegs_.writeSingle(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_s_wu(const DecodedInst* di)
{
  if (not checkRoundingModeSp(di))
    return;

  uint32_t u1 = intRegs_.read(di->op1());

#ifdef SOFT_FLOAT
  float res = f32ToFloat(ui32_to_f32(u1));
#else
  float res = float(u1);
#endif

  fpRegs_.writeSingle(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmv_w_x(const DecodedInst* di)
{
  if (not isFpLegal())
    {
      illegalInst(di);
      return;
    }

  uint32_t u1 = intRegs_.read(di->op1());

  Uint32FloatUnion ufu(u1);
  fpRegs_.writeSingle(di->op0(), ufu.f);

  markFsDirty();
}


template <>
void
Hart<uint32_t>::execFcvt_l_s(const DecodedInst* di)
{
  illegalInst(di);  // fcvt.l.s is not an RV32 instruction.
}


template <>
void
Hart<uint64_t>::execFcvt_l_s(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  SRV result = 0;
  bool valid = false;

#ifdef SOFT_FLOAT
  result = f32_to_i64(floatToF32(f1), softfloat_roundingMode, true);
  valid = true;  // We get invalid from softfloat library.
#else

  int64_t maxInt = (~uint64_t(0)) >> 1;
  int64_t minInt = int64_t(1) << 63;

  unsigned signBit = std::signbit(f1);
  if (std::isinf(f1))
    {
      if (signBit)
	result = minInt;
      else
	result = maxInt;
    }
  else if (std::isnan(f1))
    result = maxInt;
  else
    {
      double near = std::nearbyint(double(f1));
      if (near >= double(maxInt))
	result = maxInt;
      else if (near < double(minInt))
	result = minInt;
      else
	{
	  valid = true;
          result = std::lrint(f1);
	}
    }

#endif

  intRegs_.write(di->op0(), result);

  updateAccruedFpBits(0.0, not valid);

  markFsDirty();
}


template <>
void
Hart<uint32_t>::execFcvt_lu_s(const DecodedInst* di)
{
  illegalInst(di);  // RV32 does not have fcvt.lu.s
}


template <>
void
Hart<uint64_t>::execFcvt_lu_s(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeSp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  uint64_t result = 0;

#ifdef SOFT_FLOAT

  result = f32_to_ui64(floatToF32(f1), softfloat_roundingMode, true);
  if (not std::isnan(f1) and f1 < 0)
    {
      softfloat_exceptionFlags &= ~softfloat_flag_inexact;  // Should not set inexact.
      softfloat_exceptionFlags |= softfloat_flag_invalid;
    }
  updateAccruedFpBits(0.0f, false);

#else

  bool valid = false;
  bool exact = true;

  uint64_t maxUint = ~uint64_t(0);

  unsigned signBit = std::signbit(f1);
  if (std::isinf(f1))
    {
      if (signBit)
	result = 0;
      else
	result = maxUint;
    }
  else if (std::isnan(f1))
    result = maxUint;
  else if (std::signbit(f1) and f1 != 0)
    result = 0;
  else
    {
      double near = std::nearbyint(double(f1));
      if (near == 0)
        {
          result = 0;
          valid = true;
          exact = near == f1;
        }
      else if (near < 0)
        {
          result = 0;
        }
      else
        {
          // Using "near > maxUint" will not work beacuse of rounding.
          if (near >= 2*double(uint64_t(1)<<63))
            result = maxUint;
          else
            {
              // std::lprint will produce an overflow if most sig bit
              // of result is 1 (it thinks there's an overflow).  We
              // compensate with the divide multiply by 2.
              if (f1 < (uint64_t(1) << 63))
                result = std::llrint(f1);
              else
                {
                  result = std::llrint(f1/2);
                  result *= 2;
                }
              valid = true;
              exact = near == f1;
            }
        }
    }


  if (not valid)
    setFcsrFlags(FpFlags::Invalid);
  if (not exact)
    setFcsrFlags(FpFlags::Inexact);

#endif

  intRegs_.write(di->op0(), result);
  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_s_l(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeSp(di))
    return;

  SRV i1 = intRegs_.read(di->op1());

#ifdef SOFT_FLOAT
  float res = f32ToFloat(i64_to_f32(i1));
#else
  float res = float(i1);
#endif

  fpRegs_.writeSingle(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_s_lu(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeSp(di))
    return;

  URV i1 = intRegs_.read(di->op1());

#ifdef SOFT_FLOAT
  float res = f32ToFloat(ui64_to_f32(i1));
#else
  float res = float(i1);
#endif

  fpRegs_.writeSingle(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFld(const DecodedInst* di)
{
  if (not isDpLegal())
    {
      illegalInst(di);
      return;
    }

  uint32_t rd = di->op0(), rs1 = di->op1();
  SRV imm = di->op2As<SRV>();

  URV base = intRegs_.read(rs1);
  URV virtAddr = base + imm;

  ldStAddr_ = virtAddr;   // For reporting load addr in trace-mode.
  ldStAddrValid_ = true;  // For reporting load addr in trace-mode.
  unsigned ldSize = 8;
  uint64_t addr = virtAddr;

  auto secCause = SecondaryCause::NONE;
  auto cause = ExceptionCause::NONE;

#ifndef FAST_SLOPPY
  if (loadQueueEnabled_)
    removeFromLoadQueue(rs1, false);

  if (hasActiveTrigger())
    {
      if (ldStAddrTriggerHit(virtAddr, TriggerTiming::Before, true /*isLoad*/,
			     privMode_, isInterruptEnabled()))
	triggerTripped_ = true;
      if (triggerTripped_)
	return;
    }

  cause = determineLoadException(rs1, base, addr, ldSize, secCause);
  if (cause != ExceptionCause::NONE)
    {
      if (not triggerTripped_)
        initiateLoadException(cause, virtAddr, secCause);
      return;
    }
#endif

  union UDU  // Unsigned double union: reinterpret bits as unsigned or double
  {
    uint64_t u;
    double d;
  };

  uint64_t val64 = 0;
  if (memory_.read(addr, val64))
    {
      if (loadQueueEnabled_)
        {
          uint64_t prevRdVal = 0;
          peekFpReg(rd, prevRdVal);
          putInLoadQueue(ldSize, addr, rd, prevRdVal, false /*wide*/, true /*fp*/);
        }

      UDU udu;
      udu.u = val64;
      fpRegs_.write(di->op0(), udu.d);

      markFsDirty();
      return;
    }

  cause = ExceptionCause::LOAD_ACC_FAULT;
  secCause = SecondaryCause::LOAD_ACC_MEM_PROTECTION;
  if (isAddrMemMapped(addr))
    secCause = SecondaryCause::LOAD_ACC_PIC;

  initiateLoadException(cause, virtAddr, secCause);
}


template <typename URV>
void
Hart<URV>::execFsd(const DecodedInst* di)
{
  if (not isDpLegal())
    {
      illegalInst(di);
      return;
    }

  uint32_t rs1 = di->op1();
  uint32_t rs2 = di->op0();

  URV base = intRegs_.read(rs1);
  URV addr = base + di->op2As<SRV>();
  double val = fpRegs_.read(rs2);

  union UDU  // Unsigned double union: reinterpret bits as unsigned or double
  {
    uint64_t u;
    double d;
  };

  UDU udu;
  udu.d = val;

  store<uint64_t>(rs1, base, addr, udu.u);
}


template <typename URV>
void
Hart<URV>::execFmadd_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double f1 = fpRegs_.read(di->op1());
  double f2 = fpRegs_.read(di->op2());
  double f3 = fpRegs_.read(di->op3());

  bool invalid = false;
  double res = fusedMultiplyAdd(f1, f2, f3, invalid);
  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  updateAccruedFpBits(res, invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmsub_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double f1 = fpRegs_.read(di->op1());
  double f2 = fpRegs_.read(di->op2());
  double f3 = -fpRegs_.read(di->op3());

  bool invalid = false;
  double res = fusedMultiplyAdd(f1, f2, f3, invalid);

  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  updateAccruedFpBits(res, invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFnmsub_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double f1 = -fpRegs_.read(di->op1());
  double f2 = fpRegs_.read(di->op2());
  double f3 = fpRegs_.read(di->op3());

  bool invalid = false;
  double res = fusedMultiplyAdd(f1, f2, f3, invalid);
  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  updateAccruedFpBits(res, invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFnmadd_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  // we want -(f[op1] * f[op2]) - f[op3]

  double f1 = -fpRegs_.read(di->op1());
  double f2 = fpRegs_.read(di->op2());
  double f3 = -fpRegs_.read(di->op3());

  bool invalid = false;
  double res = fusedMultiplyAdd(f1, f2, f3, invalid);
  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  updateAccruedFpBits(res, invalid);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFadd_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());

#ifdef SOFT_FLOAT
  double res = f64ToDouble(f64_add(doubleToF64(d1), doubleToF64(d2)));
#else
  double res = d1 + d2;
#endif

  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsub_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());

#ifdef SOFT_FLOAT
  double res = f64ToDouble(f64_sub(doubleToF64(d1), doubleToF64(d2)));
#else
  double res = d1 - d2;
#endif

  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmul_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());

#ifdef SOFT_FLOAT
  double res = f64ToDouble(f64_mul(doubleToF64(d1), doubleToF64(d2)));
#else
  double res = d1 * d2;
#endif

  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFdiv_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());

#ifdef SOFT_FLOAT
  double res = f64ToDouble(f64_div(doubleToF64(d1), doubleToF64(d2)));
#else
  double res = d1 / d2;
#endif

  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsgnj_d(const DecodedInst* di)
{
  if (not isDpLegal())
    {
      illegalInst(di);
      return;
    }

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());
  double res = copysign(d1, d2);  // Magnitude of rs1 and sign of rs2
  fpRegs_.write(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsgnjn_d(const DecodedInst* di)
{
  if (not isDpLegal())
    {
      illegalInst(di);
      return;
    }

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());
  double res = copysign(d1, d2);  // Magnitude of rs1 and sign of rs2
  res = -res;  // Magnitude of rs1 and negative the sign of rs2
  fpRegs_.write(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsgnjx_d(const DecodedInst* di)
{
  if (not isDpLegal())
    {
      illegalInst(di);
      return;
    }

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());

  int sign1 = (std::signbit(d1) == 0) ? 0 : 1;
  int sign2 = (std::signbit(d2) == 0) ? 0 : 1;
  int sign = sign1 ^ sign2;

  double x = sign? -1 : 1;

  double res = copysign(d1, x);  // Magnitude of rs1 and sign of x
  fpRegs_.write(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmin_d(const DecodedInst* di)
{
  if (not isDpLegal())
    {
      illegalInst(di);
      return;
    }

  double in1 = fpRegs_.read(di->op1());
  double in2 = fpRegs_.read(di->op2());
  double res = 0;

  bool isNan1 = std::isnan(in1), isNan2 = std::isnan(in2);
  if (isNan1 and isNan2)
    res = std::numeric_limits<double>::quiet_NaN();
  else if (isNan1)
    res = in2;
  else if (isNan2)
    res = in1;
  else
    res = fmin(in1, in2);

  if (issnan(in1) or issnan(in2))
    setFcsrFlags(FpFlags::Invalid);
  else if (std::signbit(in1) != std::signbit(in2) and in1 == in2)
    res = std::copysign(res, -1.0);  // Make sure min(-0, +0) is -0.

  fpRegs_.write(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmax_d(const DecodedInst* di)
{
  if (not isDpLegal())
    {
      illegalInst(di);
      return;
    }

  double in1 = fpRegs_.read(di->op1());
  double in2 = fpRegs_.read(di->op2());
  double res = 0;

  bool isNan1 = std::isnan(in1), isNan2 = std::isnan(in2);
  if (isNan1 and isNan2)
    res = std::numeric_limits<double>::quiet_NaN();
  else if (isNan1)
    res = in2;
  else if (isNan2)
    res = in1;
  else
    res = std::fmax(in1, in2);

  if (issnan(in1) or issnan(in2))
    setFcsrFlags(FpFlags::Invalid);
  else if (std::signbit(in1) != std::signbit(in2) and in1 == in2)
    res = std::copysign(res, 1.0);  // Make sure max(-0, +0) is +0.

  fpRegs_.write(di->op0(), res);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_d_s(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  float f1 = fpRegs_.readSingle(di->op1());
  double res = f1;
  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_s_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());

#ifdef SOFT_FLOAT
  float res = f32ToFloat(f64_to_f32(doubleToF64(d1)));
#else
  float res = float(d1);
#endif

  if (std::isnan(res))
    res = std::numeric_limits<float>::quiet_NaN();

  fpRegs_.writeSingle(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFsqrt_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());

#ifdef SOFT_FLOAT
  double res = f64ToDouble(f64_sqrt(doubleToF64(d1)));
#else
  double res = std::sqrt(d1);
#endif

  if (std::isnan(res))
    res = std::numeric_limits<double>::quiet_NaN();

  fpRegs_.write(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFle_d(const DecodedInst* di)
{
  if (not isDpLegal())
    {
      illegalInst(di);
      return;
    }

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());

  URV res = 0;

  if (std::isnan(d1) or std::isnan(d2))
    setFcsrFlags(FpFlags::Invalid);
  else
    res = (d1 <= d2)? 1 : 0;

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFlt_d(const DecodedInst* di)
{
  if (not isDpLegal())
    {
      illegalInst(di);
      return;
    }

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());

  URV res = 0;

  if (std::isnan(d1) or std::isnan(d2))
    setFcsrFlags(FpFlags::Invalid);
  else
    res = (d1 < d2)? 1 : 0;

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFeq_d(const DecodedInst* di)
{
  if (not isDpLegal())
    {
      illegalInst(di);
      return;
    }

  double d1 = fpRegs_.read(di->op1());
  double d2 = fpRegs_.read(di->op2());

  URV res = 0;

  if (std::isnan(d1) or std::isnan(d2))
    {
      if (issnan(d1) or issnan(d2))
	setFcsrFlags(FpFlags::Invalid);
    }
  else
    res = (d1 == d2)? 1 : 0;

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFcvt_w_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());
  SRV result = 0;
  bool valid = false;

#ifdef SOFT_FLOAT
  result = f64_to_i32(doubleToF64(d1), softfloat_roundingMode, true);
  valid = true;  // We get invalid from softfloat library.
#else

  int32_t minInt = int32_t(1) << 31;
  int32_t maxInt = (~uint32_t(0)) >> 1;

  unsigned signBit = std::signbit(d1);
  if (std::isinf(d1))
    result = signBit ? minInt : maxInt;
  else if (std::isnan(d1))
    result = maxInt;
  else
    {
      double near = std::nearbyint(d1);
      if (near > double(maxInt))
	result = maxInt;
      else if (near < double(minInt))
	result = minInt;
      else
	{
	  valid = true;
	  result = SRV(int32_t(std::lrint(d1)));
	}
    }

#endif

  intRegs_.write(di->op0(), result);

  updateAccruedFpBits(0.0, not valid);

  markFsDirty();
}



template <typename URV>
void
Hart<URV>::execFcvt_wu_d(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  double d1 = fpRegs_.read(di->op1());
  SRV result = 0;

#ifdef SOFT_FLOAT

  result = SRV(int32_t(f64_to_ui32(doubleToF64(d1), softfloat_roundingMode, true)));
  if (not std::isnan(d1) and d1 < 0)
    {
      softfloat_exceptionFlags &= ~softfloat_flag_inexact;  // Should not set inexact.
      softfloat_exceptionFlags |= softfloat_flag_invalid;
    }
  updateAccruedFpBits(0.0f, false);

#else


  bool valid = false;
  bool exact = true;

  uint32_t maxUint32 = ~uint32_t(0);
  if (std::isnan(d1))
    {
      result = ~URV(0);
    }
  else if (std::signbit(d1) and d1 != 0)
    {
      result = 0;
    }
  else
    {
      double near = std::nearbyint(d1);
      if (near > double(maxUint32))
        {
          result = ~URV(0);
        }
      else if (near == 0)
        {
          result = 0;
          valid = true;
          exact = near == d1;
        }
      else if (near < 0)
        {
          result = 0;
        }
      else
        {
          result = SRV(int32_t(std::lrint(d1)));
          valid = true;
          exact = near == d1;
        }
    }

  if (not valid)
    setFcsrFlags(FpFlags::Invalid);
  if (not exact)
    setFcsrFlags(FpFlags::Inexact);

#endif

  intRegs_.write(di->op0(), result);
  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_d_w(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  int32_t i1 = intRegs_.read(di->op1());
  double res = i1;
  fpRegs_.write(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_d_wu(const DecodedInst* di)
{
  if (not checkRoundingModeDp(di))
    return;

  uint32_t i1 = intRegs_.read(di->op1());
  double res = i1;
  fpRegs_.write(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFclass_d(const DecodedInst* di)
{
  if (not isDpLegal())
    {
      illegalInst(di);
      return;
    }

  double d1 = fpRegs_.read(di->op1());
  URV result = 0;

  bool pos = not std::signbit(d1);
  int type = std::fpclassify(d1);

  if (type == FP_INFINITE)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosInfinity);
      else
	result |= URV(FpClassifyMasks::NegInfinity);
    }
  else if (type == FP_NORMAL)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosNormal);
      else
	result |= URV(FpClassifyMasks::NegNormal);
    }
  else if (type == FP_SUBNORMAL)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosSubnormal);
      else
	result |= URV(FpClassifyMasks::NegSubnormal);
    }
  else if (type == FP_ZERO)
    {
      if (pos)
	result |= URV(FpClassifyMasks::PosZero);
      else
	result |= URV(FpClassifyMasks::NegZero);
    }
  else if(type == FP_NAN)
    {
      bool quiet = mostSignificantFractionBit(d1);
      if (quiet)
	result |= URV(FpClassifyMasks::QuietNan);
      else
	result |= URV(FpClassifyMasks::SignalingNan);
    }

  intRegs_.write(di->op0(), result);
}


template <>
void
Hart<uint32_t>::execFcvt_l_d(const DecodedInst* di)
{
  illegalInst(di);  // fcvt.l.d not available in RV32
}


template <>
void
Hart<uint64_t>::execFcvt_l_d(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeDp(di))
    return;

  double f1 = fpRegs_.read(di->op1());
  SRV result = 0;
  bool valid = false;

  int64_t maxInt = (~uint64_t(0)) >> 1;
  int64_t minInt = int64_t(1) << 63;

  unsigned signBit = std::signbit(f1);
  if (std::isinf(f1))
    {
      if (signBit)
	result = minInt;
      else
	result = maxInt;
    }
  else if (std::isnan(f1))
    result = maxInt;
  else
    {
      double near = std::nearbyint(f1);

      // Note "near > double(maxInt)" will not work because of
      // rounding.
      if (near >= double(uint64_t(1) << 63))
	result = maxInt;
      else if (near < double(minInt))
	result = minInt;
      else
	{
	  valid = true;
	  result = std::lrint(f1);
	}
    }

  intRegs_.write(di->op0(), result);

  updateAccruedFpBits(0.0, not valid);

  markFsDirty();
}


template <>
void
Hart<uint32_t>::execFcvt_lu_d(const DecodedInst* di)
{
  illegalInst(di);  /// fcvt.lu.d is not available in RV32.
}


template <>
void
Hart<uint64_t>::execFcvt_lu_d(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeDp(di))
    return;

  double f1 = fpRegs_.read(di->op1());
  uint64_t result = 0;

#ifdef SOFT_FLOAT

  result = f64_to_ui64(doubleToF64(f1), softfloat_roundingMode, true);
  if (not std::isnan(f1) and f1 < 0)
    {
      softfloat_exceptionFlags &= ~softfloat_flag_inexact;  // Should not set inexact.
      softfloat_exceptionFlags |= softfloat_flag_invalid;
    }
  updateAccruedFpBits(0.0f, false);

#else

  bool valid = false;
  bool exact = true;

  uint64_t maxUint = ~uint64_t(0);

  unsigned signBit = std::signbit(f1);
  if (std::isinf(f1))
    {
      if (signBit)
	result = 0;
      else
	result = maxUint;
    }
  else if (std::isnan(f1))
    result = maxUint;
  else if (std::signbit(f1) and f1 != 0)
    result = 0;
  else
    {
      double near = std::nearbyint(f1);
      if (near == 0)
        {
          result = 0;
          valid = true;
          exact = near == f1;
        }
      else if (near < 0)
        {
          result = 0;
        }
      else
        {
          // Using "near > maxUint" will not work beacuse of rounding.
          if (near >= 2*double(uint64_t(1)<<63))
            result = maxUint;
          else
            {
              // std::llrint will produce an overflow if most sig bit
              // of result is 1 (it thinks there's an overflow).  We
              // compensate with the divide multiply by 2.
              if (f1 < (uint64_t(1) << 63))
                result = std::llrint(f1);
              else
                {
                  result = std::llrint(f1/2);
                  result *= 2;
                }
              valid = true;
              exact = near == f1;
            }
        }
    }


  if (not valid)
    setFcsrFlags(FpFlags::Invalid);
  if (not exact)
    setFcsrFlags(FpFlags::Inexact);

#endif

  intRegs_.write(di->op0(), result);
  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_d_l(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeDp(di))
    return;

  SRV i1 = intRegs_.read(di->op1());
  double res = double(i1);
  fpRegs_.write(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFcvt_d_lu(const DecodedInst* di)
{
  if (not isRv64())
    {
      illegalInst(di);
      return;
    }

  if (not checkRoundingModeDp(di))
    return;

  URV i1 = intRegs_.read(di->op1());
  double res = double(i1);
  fpRegs_.write(di->op0(), res);

  updateAccruedFpBits(res, false /*invalid*/);

  markFsDirty();
}


template <typename URV>
void
Hart<URV>::execFmv_d_x(const DecodedInst* di)
{
  if (not isRv64() or not isDpLegal())
    {
      illegalInst(di);
      return;
    }

  uint64_t u1 = intRegs_.read(di->op1());

  union UDU  // Unsigned double union: reinterpret bits as unsigned or double
  {
    uint64_t u;
    double d;
  };

  UDU udu;
  udu.u = u1;

  fpRegs_.write(di->op0(), udu.d);

  markFsDirty();
}


// In 32-bit harts, fmv_x_d is an illegal instruction.
template <>
void
Hart<uint32_t>::execFmv_x_d(const DecodedInst* di)
{
  illegalInst(di);
}


template <>
void
Hart<uint64_t>::execFmv_x_d(const DecodedInst* di)
{
  if (not isRv64() or not isDpLegal())
    {
      illegalInst(di);
      return;
    }

  uint64_t v1 = fpRegs_.readBitsRaw(di->op1());
  intRegs_.write(di->op0(), v1);
}


template class WdRiscv::Hart<uint32_t>;
template class WdRiscv::Hart<uint64_t>;
