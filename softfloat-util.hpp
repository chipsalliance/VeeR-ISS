#pragma once

#ifdef SOFT_FLOAT
extern "C" {
#include <softfloat.h>
}


namespace WdRiscv
{
  
  /// Convert softfloat float16_t type to Float16.
  inline Float16 softToNative(float16_t f16)
  {
    return Float16::fromBits(f16.v);
  }


  /// Convert softfloat float32_t type to float.
  inline float softToNative(float32_t f32)
  {
    Uint32FloatUnion tmp(f32.v);
    return tmp.f;
  }


  /// Convert softfloat float64_t to double.
  inline double softToNative(float64_t f64)
  {
    Uint64DoubleUnion tmp(f64.v);
    return tmp.d;
  }


  /// Convert Float16 to a softfloat float16_t
  inline float16_t
  nativeToSoft(Float16 x)
  {
    return float16_t{x.bits()};
  }


  /// Convert a native float to a softfloat float32_t
  inline float32_t
  nativeToSoft(float x)
  {
    Uint32FloatUnion tmp(x);
    return float32_t{tmp.u};
  }


  /// Convert a native double to a softfloat float64_t
  inline float64_t
  nativeToSoft(double x)
  {
    Uint64DoubleUnion tmp(x);
    return float64_t{tmp.u};
  }


  /// Perform a floating point add using the softfloat library.
  inline float
  softAdd(float a, float b)
  {
    float res = softToNative(f32_add(nativeToSoft(a), nativeToSoft(b)));
    return res;
  }


  /// Perform a floating point add using the softfloat library.
  inline double
  softAdd(double a, double b)
  {
    double res = softToNative(f64_add(nativeToSoft(a), nativeToSoft(b)));
    return res;
  }


  /// Perform a floating point add using the softfloat library.
  inline Float16
  softAdd(Float16 a, Float16 b)
  {
    Float16 res = softToNative(f16_add(nativeToSoft(a), nativeToSoft(b)));
    return res;
  }


  /// Perform a floating point mul using the softfloat library.
  inline float
  softMul(float a, float b)
  {
    float res = softToNative(f32_mul(nativeToSoft(a), nativeToSoft(b)));
    return res;
  }


  /// Perform a floating point mul using the softfloat library.
  inline double
  softMul(double a, double b)
  {
    double res = softToNative(f64_mul(nativeToSoft(a), nativeToSoft(b)));
    return res;
  }


  /// Perform a floating point mul using the softfloat library.
  inline Float16
  softMul(Float16 a, Float16 b)
  {
    Float16 res = softToNative(f16_mul(nativeToSoft(a), nativeToSoft(b)));
    return res;
  }


  /// Perform a floating point divide using the softfloat library.
  inline float
  softDiv(float a, float b)
  {
    float res = softToNative(f32_div(nativeToSoft(a), nativeToSoft(b)));
    return res;
  }


  /// Perform a floating point divide using the softfloat library.
  inline double
  softDiv(double a, double b)
  {
    double res = softToNative(f64_div(nativeToSoft(a), nativeToSoft(b)));
    return res;
  }


  /// Perform a floating point divide using the softfloat library.
  inline Float16
  softDiv(Float16 a, Float16 b)
  {
    Float16 res = softToNative(f16_div(nativeToSoft(a), nativeToSoft(b)));
    return res;
  }


  /// Perform a floating point sqrt using the softfloat library.
  inline float
  softSqrt(float a)
  {
    float res = softToNative(f32_sqrt(nativeToSoft(a)));
    return res;
  }


  /// Perform a floating point sqrt using the softfloat library.
  inline double
  softSqrt(double a)
  {
    double res = softToNative(f64_sqrt(nativeToSoft(a)));
    return res;
  }


  /// Perform a floating point sqrt using the softfloat library.
  inline Float16
  softSqrt(Float16 a)
  {
    Float16 res = softToNative(f16_sqrt(nativeToSoft(a)));
    return res;
  }
}

#endif

