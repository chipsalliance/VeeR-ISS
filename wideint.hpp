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


/// Provide fixed size integer types wider than 64-bit.
/// The reason we do not use boost is that we want 2's complement representation
/// of signed integers.


#include <type_traits>
#include <cstdint>


namespace WdRiscv
{

  class Int128;

  /// Unsigned 128-bit integer.
  class Uint128
  {
  public:

    typedef Uint128  SelfType;
    typedef Int128   SignedType;
    typedef uint64_t HalfType;
    typedef uint32_t QuarterType;

    static int width()     { return 8*sizeof(SelfType); }
    static int halfWidth() { return 8*sizeof(HalfType); }

    /// Default constructor.
    Uint128()
    { }

    /// Copy constructor.
    Uint128(const Uint128& x)
      : low_(x.low_), high_(x.high_)
    { }

    /// Construct from a 128-bit int: Copy bits.
    Uint128(const Int128& x);

    /// Construct from a 64-bit unsigned int.
    Uint128(HalfType x)
      : low_(x), high_(0)
    { }

    /// Construct from a pair of half-type ints.
    Uint128(HalfType high, HalfType low)
      : low_(low), high_(high)
    { }

    /// Assignment constructor.
    SelfType& operator = (const SelfType& x)
    { low_ = x.low_; high_ = x.high_; return *this; }

    /// Return least sig half.
    HalfType low() const
    { return low_; }

    /// Return most sig half.
    HalfType high() const
    { return high_; }

    /// Convert to a built-in integral type.
    template <typename INT,
              std::enable_if_t<std::is_integral<INT>::value, int> = 0>
    explicit operator INT() const
    { return low_; }

    SelfType& operator += (const SelfType& x)
    {
      HalfType prevLow = low_;
      low_ += x.low_;
      high_ += x.high_;
      if (low_ < prevLow)
        high_++;
      return *this;
    }

    SelfType& operator -= (const SelfType& x)
    {
      HalfType prevLow = low_;
      low_ -= x.low_;
      high_ -= x.high_;
      if (low_ > prevLow)
        high_--;
      return *this;
    }

    SelfType& operator |= (const SelfType& x)
    { low_ |= x.low_; high_ |= x.high_; return *this; }

    SelfType& operator &= (const SelfType& x)
    { low_ &= x.low_; high_ &= x.high_; return *this; }

    SelfType& operator ^= (const SelfType& x)
    { low_ ^= x.low_; high_ ^= x.high_; return *this; }

    SelfType& operator ++ ()
    { *this += 1; return *this; }

    SelfType operator ++ (int)
    { SelfType temp = *this; temp += 1; return temp; }

    SelfType& operator -- ()
    { *this -= 1; return *this; }

    SelfType operator -- (int)
    { SelfType temp = *this; temp -= 1; return temp; }

    SelfType operator ~ () const
    { SelfType temp( ~high_, ~low_); return temp; }

    SelfType& operator *= (const SelfType& x);

    SelfType& operator /= (const SelfType& x);

    SelfType& operator %= (const SelfType& x);

    SelfType& operator >>= (int n);

    SelfType& operator <<= (int n);

    bool operator == (const SelfType& x) const
    { return high_ == x.high_ and low_ == x.low_; }

    bool operator != (const SelfType& x) const
    { return not (*this == x); }

    bool operator < (const SelfType& x) const
    { return high_ < x.high_ or (high_ == x.high_ and low_ < x.low_); }

    bool operator > (const SelfType& x) const
    { return high_ > x.high_ or (high_ == x.high_ and low_ > x.low_); }

    bool operator <= (const SelfType& x) const
    { return high_ < x.high_ or (high_ == x.high_ and low_ <= x.low_); }

    bool operator >= (const SelfType& x) const
    { return high_ > x.high_ or (high_ == x.high_ and low_ >= x.low_); }

  protected:

    HalfType low_ = 0;
    HalfType high_ = 0;
  };


  /// Signed 128-bit integer.
  class Int128
  {
  public:

    typedef Int128   SelfType;
    typedef Uint128  UnsignedType;
    typedef int64_t  HalfType;
    typedef uint64_t HalfUnsigned;

    static int width()     { return 8*sizeof(SelfType); }
    static int halfWidth() { return 8*sizeof(HalfType); }

    /// Copy constructor.
    Int128(const Int128& x)
      : low_(x.low_), high_(x.high_)
    { }

    /// Construct from an unsigned 128-bit int: Copy bits.
    Int128(Uint128 x)
      : low_(x.low()), high_(x.high())
    { }

    /// Construct from a 64-bit int.
    Int128(int64_t x)
      : low_(x), high_(x < 0? ~int64_t(0) : 0)
    { }

    /// Construct from a 64-bit unsigned.
    Int128(uint64_t x)
      : low_(x), high_(0)
    { }

    /// Construct from a 32-bit int.
    Int128(int32_t x)
      : low_(int64_t(x)), high_(0)
    { }

    /// Construct from a 64-bit unsigned.
    Int128(uint32_t x)
      : low_(x), high_(0)
    { }

    /// Construct from a pair of half-type ints.
    Int128(HalfType high, HalfType low)
      : low_(low), high_(high)
    { }

    /// Assignment constructor.
    SelfType& operator = (const SelfType& x)
    { low_ = x.low_; high_ = x.high_; return *this; }

    /// Return least sig half.
    HalfType low() const
    { return low_; }

    /// Return most sig half.
    HalfType high() const
    { return high_; }

    /// Convert to a built-in integral type.
    template <typename INT,
              std::enable_if_t<std::is_integral<INT>::value, int> = 0>
    explicit operator INT() const
    { return low_; }

    SelfType& operator += (const SelfType& x)
    {
      UnsignedType* a = reinterpret_cast<UnsignedType*> (this);
      const UnsignedType* b = reinterpret_cast<const UnsignedType*> (&x);
      *a += *b;
      return *this;
    }

    SelfType& operator -= (const SelfType& x)
    {
      UnsignedType* a = reinterpret_cast<UnsignedType*> (this);
      const UnsignedType* b = reinterpret_cast<const UnsignedType*> (&x);
      *a -= *b;
      return *this;
    }

    SelfType& operator |= (const SelfType& x)
    { low_ |= x.low_; high_ |= x.high_; return *this; }

    SelfType& operator &= (const SelfType& x)
    { low_ &= x.low_; high_ &= x.high_; return *this; }

    SelfType& operator ^= (const SelfType& x)
    { low_ ^= x.low_; high_ ^= x.high_; return *this; }

    SelfType& operator ++ ()
    { *this += 1; return *this; }

    SelfType operator ++ (int)
    { SelfType temp = *this; temp += 1; return temp; }

    SelfType& operator -- ()
    { *this -= 1; return *this; }

    SelfType operator -- (int)
    { SelfType temp = *this; temp -= 1; return temp; }

    SelfType operator ~ () const
    { SelfType temp( ~high_, ~low_); return temp; }

    SelfType& operator *= (const SelfType& x);

    SelfType& operator /= (const SelfType& x);

    SelfType& operator %= (const SelfType& x);

    SelfType& operator >>= (int n);

    SelfType& operator <<= (int n);

    bool operator == (const SelfType& x) const
    { return high_ == x.high_ and low_ == x.low_; }

    bool operator != (const SelfType& x) const
    { return not (*this == x); }

    bool operator < (const SelfType& x) const
    {
      return (high_ < x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) < HalfUnsigned(x.low_)));
    }

    bool operator > (const SelfType& x) const
    {
      return (high_ > x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) > HalfUnsigned(x.low_)));
    }

    bool operator <= (const SelfType& x) const
    {
      return (high_ < x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) <= HalfUnsigned(x.low_)));
    }

    bool operator >= (const SelfType& x) const
    {
      return (high_ > x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) >= HalfUnsigned(x.low_)));
    }

  protected:

    HalfType low_ = 0;
    HalfType high_ = 0;
  };


  class Int256;

  /// Unsigned 256-bit integer.
  class Uint256
  {
  public:

    typedef Uint256  SelfType;
    typedef Int256   SignedType;
    typedef Uint128  HalfType;
    typedef uint64_t QuarterType;

    static int width()     { return 8*sizeof(SelfType); }
    static int halfWidth() { return 8*sizeof(HalfType); }

    /// Default constructor.
    Uint256()
    { }

    /// Copy constructor.
    Uint256(const Uint256& x)
      : low_(x.low_), high_(x.high_)
    { }

    /// Construct from a 256-bit int: Copy bits.
    Uint256(const Int256& x);

    /// Construct from a 128-bit unsigned int.
    Uint256(const HalfType& x)
      : low_(x), high_(0)
    { }

    /// Construct from a 64-bit unsigned int.
    Uint256(uint64_t x)
      : low_(x), high_(0)
    { }

    /// Construct from a pair of half-type ints.
    Uint256(HalfType high, HalfType low)
      : low_(low), high_(high)
    { }

    /// Assignment constructor.
    SelfType& operator = (const SelfType& x)
    { low_ = x.low_; high_ = x.high_; return *this; }

    /// Return least sig half.
    HalfType low() const
    { return low_; }

    /// Return most sig half.
    HalfType high() const
    { return high_; }

    /// Convert to a half-type.
    operator HalfType() const
    { return low_; }

    /// Convert to a built-in integral type.
    template <typename INT,
              std::enable_if_t<std::is_integral<INT>::value, int> = 0>
    explicit operator INT() const
    { return INT(low_); }

    SelfType& operator += (const SelfType& x)
    {
      HalfType prevLow = low_;
      low_ += x.low_;
      high_ += x.high_;
      if (low_ < prevLow)
        high_++;
      return *this;
    }

    SelfType& operator -= (const SelfType& x)
    {
      HalfType prevLow = low_;
      low_ -= x.low_;
      high_ -= x.high_;
      if (low_ > prevLow)
        high_--;
      return *this;
    }

    SelfType& operator |= (const SelfType& x)
    { low_ |= x.low_; high_ |= x.high_; return *this; }

    SelfType& operator &= (const SelfType& x)
    { low_ &= x.low_; high_ &= x.high_; return *this; }

    SelfType& operator ^= (const SelfType& x)
    { low_ ^= x.low_; high_ ^= x.high_; return *this; }

    SelfType& operator ++ ()
    { *this += 1; return *this; }

    SelfType operator ++ (int)
    { SelfType temp = *this; temp += 1; return temp; }

    SelfType& operator -- ()
    { *this -= 1; return *this; }

    SelfType operator -- (int)
    { SelfType temp = *this; temp -= 1; return temp; }

    SelfType operator ~ () const
    { SelfType temp( ~high_, ~low_); return temp; }

    SelfType& operator *= (const SelfType& x);

    SelfType& operator /= (const SelfType& x);

    SelfType& operator %= (const SelfType& x);

    SelfType& operator >>= (int n);

    SelfType& operator <<= (int n);

    bool operator == (const SelfType& x) const
    { return high_ == x.high_ and low_ == x.low_; }

    bool operator != (const SelfType& x) const
    { return not (*this == x); }

    bool operator < (const SelfType& x) const
    { return high_ < x.high_ or (high_ == x.high_ and low_ < x.low_); }

    bool operator > (const SelfType& x) const
    { return high_ > x.high_ or (high_ == x.high_ and low_ > x.low_); }

    bool operator <= (const SelfType& x) const
    { return high_ < x.high_ or (high_ == x.high_ and low_ <= x.low_); }

    bool operator >= (const SelfType& x) const
    { return high_ > x.high_ or (high_ == x.high_ and low_ >= x.low_); }

  protected:

    HalfType low_ = 0;
    HalfType high_ = 0;
  };


  /// Signed 256-bit integer.
  class Int256
  {
  public:

    typedef Int256  SelfType;
    typedef Uint256 UnsignedType;
    typedef Int128  HalfType;
    typedef Uint128 HalfUnsigned;

    static int width()     { return 8*sizeof(SelfType); }
    static int halfWidth() { return 8*sizeof(HalfType); }

    /// Default constructor.
    Int256()
    { }

    /// Copy constructor.
    Int256(const Int256& x)
      : low_(x.low_), high_(x.high_)
    { }

    /// Construct from an unsigned 256-bit int: Copy bits.
    Int256(const Uint256& x)
      : low_(x.low()), high_(x.high())
    { }

    /// Construct from a 128-bit int.
    Int256(const HalfType& x)
      : low_(x), high_(x < HalfType(0)? ~HalfType(0) : HalfType(0))
    { }

    /// Construct from a 128-bit unsigned int.
    Int256(const HalfUnsigned& x)
      : low_(x), high_(0)
    { }

    /// Construct from a 64-bit int.
    Int256(int64_t x)
      : low_(x), high_(x < 0? ~HalfType(0) : HalfType(0))
    { }

    /// Construct from a pair of half-type ints.
    Int256(HalfType high, HalfType low)
      : low_(low), high_(high)
    { }

    /// Assignment constructor.
    SelfType& operator = (const SelfType& x)
    { low_ = x.low_; high_ = x.high_; return *this; }

    /// Return least sig half.
    HalfType low() const
    { return low_; }

    /// Return most sig half.
    HalfType high() const
    { return high_; }

    /// Convert to a half.
    operator HalfType() const
    { return low_; }

    /// Convert to a built-in integral type.
    template <typename INT,
              std::enable_if_t<std::is_integral<INT>::value, int> = 0>
    explicit operator INT() const
    { return INT(low_); }

    SelfType& operator += (const SelfType& x)
    {
      UnsignedType* a = reinterpret_cast<UnsignedType*> (this);
      const UnsignedType* b = reinterpret_cast<const UnsignedType*> (&x);
      *a += *b;
      return *this;
    }

    SelfType& operator -= (const SelfType& x)
    {
      UnsignedType* a = reinterpret_cast<UnsignedType*> (this);
      const UnsignedType* b = reinterpret_cast<const UnsignedType*> (&x);
      *a -= *b;
      return *this;
    }

    SelfType& operator |= (const SelfType& x)
    { low_ |= x.low_; high_ |= x.high_; return *this; }

    SelfType& operator &= (const SelfType& x)
    { low_ &= x.low_; high_ &= x.high_; return *this; }

    SelfType& operator ^= (const SelfType& x)
    { low_ ^= x.low_; high_ ^= x.high_; return *this; }

    SelfType& operator ++ ()
    { *this += 1; return *this; }

    SelfType operator ++ (int)
    { SelfType temp = *this; temp += 1; return temp; }

    SelfType& operator -- ()
    { *this -= 1; return *this; }

    SelfType operator -- (int)
    { SelfType temp = *this; temp -= 1; return temp; }

    SelfType operator ~ () const
    { SelfType temp( ~high_, ~low_); return temp; }

    SelfType& operator *= (const SelfType& x);

    SelfType& operator /= (const SelfType& x);

    SelfType& operator %= (const SelfType& x);

    SelfType& operator >>= (int n);

    SelfType& operator <<= (int n);

    bool operator == (const SelfType& x) const
    { return high_ == x.high_ and low_ == x.low_; }

    bool operator != (const SelfType& x) const
    { return not (*this == x); }

    bool operator < (const SelfType& x) const
    {
      return (high_ < x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) < HalfUnsigned(x.low_)));
    }

    bool operator > (const SelfType& x) const
    {
      return (high_ > x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) > HalfUnsigned(x.low_)));
    }

    bool operator <= (const SelfType& x) const
    {
      return (high_ < x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) <= HalfUnsigned(x.low_)));
    }

    bool operator >= (const SelfType& x) const
    {
      return (high_ > x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) >= HalfUnsigned(x.low_)));
    }

  protected:

    HalfType low_ = 0;
    HalfType high_ = 0;
  };


  class Int512;

  /// Unsigned 512-bit integer.
  class Uint512
  {
  public:

    typedef Uint512  SelfType;
    typedef Int512   SignedType;
    typedef Uint256  HalfType;
    typedef Uint128  QuarterType;

    static int width()     { return 8*sizeof(SelfType); }
    static int halfWidth() { return 8*sizeof(HalfType); }

    /// Default constructor.
    Uint512()
    { }

    /// Copy constructor.
    Uint512(const Uint512& x)
      : low_(x.low_), high_(x.high_)
    { }

    /// Construct from a 512-bit int: Copy bits.
    Uint512(const Int512& x);

    /// Construct from a 256-bit unsigned int.
    Uint512(HalfType x)
      : low_(x), high_(0)
    { }

    /// Construct from a 128-bit unsigned int.
    Uint512(Uint128 x)
      : low_(x), high_(0)
    { }

    /// Construct from a 64-bit unsigned int.
    Uint512(uint64_t x)
      : low_(x), high_(0)
    { }

    /// Construct from a pair of 256-bit unsigned ints.
    Uint512(Uint256 high, Uint256 low)
      : low_(low), high_(high)
    { }

    /// Assignment constructor.
    SelfType& operator = (const SelfType& x)
    { low_ = x.low_; high_ = x.high_; return *this; }

    /// Return least sig half.
    HalfType low() const
    { return low_; }

    /// Return most sig half.
    HalfType high() const
    { return high_; }

    /// Convert to a half-type.
    operator HalfType() const
    { return low_; }

    /// Convert to a built-in integral type.
    template <typename INT,
              std::enable_if_t<std::is_integral<INT>::value, int> = 0>
    explicit operator INT() const
    { return INT(low_); }

    SelfType& operator += (const SelfType& x)
    {
      HalfType prevLow = low_;
      low_ += x.low_;
      high_ += x.high_;
      if (low_ < prevLow)
        high_++;
      return *this;
    }

    SelfType& operator -= (const SelfType& x)
    {
      HalfType prevLow = low_;
      low_ -= x.low_;
      high_ -= x.high_;
      if (low_ > prevLow)
        high_--;
      return *this;
    }

    SelfType& operator |= (const SelfType& x)
    { low_ |= x.low_; high_ |= x.high_; return *this; }

    SelfType& operator &= (const SelfType& x)
    { low_ &= x.low_; high_ &= x.high_; return *this; }

    SelfType& operator ^= (const SelfType& x)
    { low_ ^= x.low_; high_ ^= x.high_; return *this; }

    SelfType& operator ++ ()
    { *this += 1; return *this; }

    SelfType operator ++ (int)
    { SelfType temp = *this; temp += 1; return temp; }

    SelfType& operator -- ()
    { *this -= 1; return *this; }

    SelfType operator -- (int)
    { SelfType temp = *this; temp -= 1; return temp; }

    SelfType operator ~ () const
    { SelfType temp( ~high_, ~low_); return temp; }

    SelfType& operator *= (const SelfType& x);

    SelfType& operator /= (const SelfType& x);

    SelfType& operator %= (const SelfType& x);

    SelfType& operator >>= (int n);

    SelfType& operator <<= (int n);

    bool operator == (const SelfType& x) const
    { return high_ == x.high_ and low_ == x.low_; }

    bool operator != (const SelfType& x) const
    { return not (*this == x); }

    bool operator < (const SelfType& x) const
    { return high_ < x.high_ or (high_ == x.high_ and low_ < x.low_); }

    bool operator > (const SelfType& x) const
    { return high_ > x.high_ or (high_ == x.high_ and low_ > x.low_); }

    bool operator <= (const SelfType& x) const
    { return high_ < x.high_ or (high_ == x.high_ and low_ <= x.low_); }

    bool operator >= (const SelfType& x) const
    { return high_ > x.high_ or (high_ == x.high_ and low_ >= x.low_); }

  protected:

    HalfType low_ = 0;
    HalfType high_ = 0;
  };


  /// Signed 512-bit integer.
  class Int512
  {
  public:

    typedef Int512  SelfType;
    typedef Uint512 UnsignedType;
    typedef Int256  HalfType;
    typedef Uint256 HalfUnsigned;

    static int width()     { return 8*sizeof(SelfType); }
    static int halfWidth() { return 8*sizeof(HalfType); }

    /// Default constructor.
    Int512()
    { }

    /// Copy constructor.
    Int512(const Int512& x)
      : low_(x.low_), high_(x.high_)
    { }

    /// Construct from an unsigned 512-bit int: Copy bits.
    Int512(Uint512 x)
      : low_(x.low()), high_(x.high())
    { }

    /// Construct from a 256-bit int.
    Int512(HalfType x)
      : low_(x), high_(x < HalfType(0)? ~HalfType(0) : HalfType(0))
    { }

    /// Construct from a 256-bit unsigned int.
    Int512(HalfUnsigned x)
      : low_(x), high_(0)
    { }

    /// Construct from a 128-bit int.
    Int512(Int128 x)
      : low_(x), high_(x < HalfType(0)? ~HalfType(0) : HalfType(0))
    { }

    /// Construct from a 64-bit int.
    Int512(int64_t x)
      : low_(x), high_(x < 0? ~HalfType(0) : HalfType(0))
    { }

    /// Construct from a pair of 256-bit ints.
    Int512(HalfType high, HalfType low)
      : low_(low), high_(high)
    { }

    /// Assignment constructor.
    SelfType& operator = (const SelfType& x)
    { low_ = x.low_; high_ = x.high_; return *this; }

    /// Return least sig half.
    HalfType low() const
    { return low_; }

    /// Return most sig half.
    HalfType high() const
    { return high_; }

    /// Convert to a half.
    operator HalfType() const
    { return low_; }

    /// Convert to a built-in integral type.
    template <typename INT,
              std::enable_if_t<std::is_integral<INT>::value, int> = 0>
    explicit operator INT() const
    { return INT(low_); }

    SelfType& operator += (const SelfType& x)
    {
      UnsignedType* a = reinterpret_cast<UnsignedType*> (this);
      const UnsignedType* b = reinterpret_cast<const UnsignedType*> (&x);
      *a += *b;
      return *this;
    }

    SelfType& operator -= (const SelfType& x)
    {
      UnsignedType* a = reinterpret_cast<UnsignedType*> (this);
      const UnsignedType* b = reinterpret_cast<const UnsignedType*> (&x);
      *a -= *b;
      return *this;
    }

    SelfType& operator |= (const SelfType& x)
    { low_ |= x.low_; high_ |= x.high_; return *this; }

    SelfType& operator &= (const SelfType& x)
    { low_ &= x.low_; high_ &= x.high_; return *this; }

    SelfType& operator ^= (const SelfType& x)
    { low_ ^= x.low_; high_ ^= x.high_; return *this; }

    SelfType& operator ++ ()
    { *this += 1; return *this; }

    SelfType operator ++ (int)
    { SelfType temp = *this; temp += 1; return temp; }

    SelfType& operator -- ()
    { *this -= 1; return *this; }

    SelfType operator -- (int)
    { SelfType temp = *this; temp -= 1; return temp; }

    SelfType operator ~ () const
    { SelfType temp( ~high_, ~low_); return temp; }

    SelfType& operator *= (const SelfType& x);

    SelfType& operator /= (const SelfType& x);

    SelfType& operator %= (const SelfType& x);

    SelfType& operator >>= (int n);

    SelfType& operator <<= (int n);

    bool operator == (const SelfType& x) const
    { return high_ == x.high_ and low_ == x.low_; }

    bool operator != (const SelfType& x) const
    { return not (*this == x); }

    bool operator < (const SelfType& x) const
    {
      return (high_ < x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) < HalfUnsigned(x.low_)));
    }

    bool operator > (const SelfType& x) const
    {
      return (high_ > x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) > HalfUnsigned(x.low_)));
    }

    bool operator <= (const SelfType& x) const
    {
      return (high_ < x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) <= HalfUnsigned(x.low_)));
    }

    bool operator >= (const SelfType& x) const
    {
      return (high_ > x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) >= HalfUnsigned(x.low_)));
    }

  protected:

    HalfType low_ = 0;
    HalfType high_ = 0;
  };


  class Int1024;

  /// Unsigned 1024-bit integer.
  class Uint1024
  {
  public:

    typedef Uint1024  SelfType;
    typedef Int1024   SignedType;
    typedef Uint512   HalfType;
    typedef Uint256   QuarterType;

    static int width()     { return 8*sizeof(SelfType); }
    static int halfWidth() { return 8*sizeof(HalfType); }

    /// Default constructor.
    Uint1024()
    { }

    /// Copy constructor.
    Uint1024(const Uint1024& x)
      : low_(x.low_), high_(x.high_)
    { }

    /// Construct from a 1024-bit int: Copy bits.
    Uint1024(const Int1024& x);

    /// Construct from a 512-bit unsigned int.
    Uint1024(HalfType x)
      : low_(x), high_(0)
    { }

    /// Construct from a 256-bit unsigned int.
    Uint1024(Uint256 x)
      : low_(x), high_(0)
    { }

    /// Construct from a 128-bit unsigned int.
    Uint1024(Uint128 x)
      : low_(x), high_(0)
    { }

    /// Construct from a 64-bit unsigned int.
    Uint1024(uint64_t x)
      : low_(x), high_(0)
    { }

    /// Construct from a pair of 512-bit unsigned ints.
    Uint1024(Uint512 high, Uint512 low)
      : low_(low), high_(high)
    { }

    /// Assignment constructor.
    SelfType& operator = (const SelfType& x)
    { low_ = x.low_; high_ = x.high_; return *this; }

    /// Return least sig half.
    HalfType low() const
    { return low_; }

    /// Return most sig half.
    HalfType high() const
    { return high_; }

    /// Convert to a half-type.
    operator HalfType() const
    { return low_; }

    /// Convert to a built-in integral type.
    template <typename INT,
              std::enable_if_t<std::is_integral<INT>::value, int> = 0>
    explicit operator INT() const
    { return INT(low_); }

    SelfType& operator += (const SelfType& x)
    {
      HalfType prevLow = low_;
      low_ += x.low_;
      high_ += x.high_;
      if (low_ < prevLow)
        high_++;
      return *this;
    }

    SelfType& operator -= (const SelfType& x)
    {
      HalfType prevLow = low_;
      low_ -= x.low_;
      high_ -= x.high_;
      if (low_ > prevLow)
        high_--;
      return *this;
    }

    SelfType& operator |= (const SelfType& x)
    { low_ |= x.low_; high_ |= x.high_; return *this; }

    SelfType& operator &= (const SelfType& x)
    { low_ &= x.low_; high_ &= x.high_; return *this; }

    SelfType& operator ^= (const SelfType& x)
    { low_ ^= x.low_; high_ ^= x.high_; return *this; }

    SelfType& operator ++ ()
    { *this += 1; return *this; }

    SelfType operator ++ (int)
    { SelfType temp = *this; temp += 1; return temp; }

    SelfType& operator -- ()
    { *this -= 1; return *this; }

    SelfType operator -- (int)
    { SelfType temp = *this; temp -= 1; return temp; }

    SelfType operator ~ () const
    { SelfType temp( ~high_, ~low_); return temp; }

    SelfType& operator *= (const SelfType& x);

    SelfType& operator /= (const SelfType& x);

    SelfType& operator %= (const SelfType& x);

    SelfType& operator >>= (int n);

    SelfType& operator <<= (int n);

    bool operator == (const SelfType& x) const
    { return high_ == x.high_ and low_ == x.low_; }

    bool operator != (const SelfType& x) const
    { return not (*this == x); }

    bool operator < (const SelfType& x) const
    { return high_ < x.high_ or (high_ == x.high_ and low_ < x.low_); }

    bool operator > (const SelfType& x) const
    { return high_ > x.high_ or (high_ == x.high_ and low_ > x.low_); }

    bool operator <= (const SelfType& x) const
    { return high_ < x.high_ or (high_ == x.high_ and low_ <= x.low_); }

    bool operator >= (const SelfType& x) const
    { return high_ > x.high_ or (high_ == x.high_ and low_ >= x.low_); }

  protected:

    HalfType low_ = 0;
    HalfType high_ = 0;
  };


  /// Signed 1024-bit integer.
  class Int1024
  {
  public:

    typedef Int1024  SelfType;
    typedef Uint1024 UnsignedType;
    typedef Int512   HalfType;
    typedef Uint512  HalfUnsigned;

    static int width()     { return 8*sizeof(SelfType); }
    static int halfWidth() { return 8*sizeof(HalfType); }

    /// Default constructor.
    Int1024()
    { }

    /// Copy constructor.
    Int1024(const Int1024& x)
      : low_(x.low_), high_(x.high_)
    { }

    /// Construct from an unsigned 1024-bit int: Copy bits.
    Int1024(Uint1024 x)
      : low_(x.low()), high_(x.high())
    { }

    /// Construct from a 512-bit int.
    Int1024(HalfType x)
      : low_(x), high_(x < HalfType(0)? ~HalfType(0) : HalfType(0))
    { }

    /// Construct from a 512-bit unsigned int.
    Int1024(HalfUnsigned x)
      : low_(x), high_(0)
    { }

    /// Construct from a 256-bit int.
    Int1024(Int256 x)
      : low_(x), high_(x < HalfType(0)? ~HalfType(0) : HalfType(0))
    { }

    /// Construct from a 128-bit int.
    Int1024(Int128 x)
      : low_(x), high_(x < 0? ~HalfType(0) : HalfType(0))
    { }

    /// Construct from a 64-bit int.
    Int1024(int64_t x)
      : low_(x), high_(x < 0? ~HalfType(0) : HalfType(0))
    { }

    /// Construct from a pair of 512-bit ints.
    Int1024(HalfType high, HalfType low)
      : low_(low), high_(high)
    { }

    /// Assignment constructor.
    SelfType& operator = (const SelfType& x)
    { low_ = x.low_; high_ = x.high_; return *this; }

    /// Return least sig half.
    HalfType low() const
    { return low_; }

    /// Return most sig half.
    HalfType high() const
    { return high_; }

    /// Convert to a half.
    operator HalfType() const
    { return low_; }

    /// Convert to a built-in integral type.
    template <typename INT,
              std::enable_if_t<std::is_integral<INT>::value, int> = 0>
    explicit operator INT() const
    { return INT(low_); }

    SelfType& operator += (const SelfType& x)
    {
      UnsignedType* a = reinterpret_cast<UnsignedType*> (this);
      const UnsignedType* b = reinterpret_cast<const UnsignedType*> (&x);
      *a += *b;
      return *this;
    }

    SelfType& operator -= (const SelfType& x)
    {
      UnsignedType* a = reinterpret_cast<UnsignedType*> (this);
      const UnsignedType* b = reinterpret_cast<const UnsignedType*> (&x);
      *a -= *b;
      return *this;
    }

    SelfType& operator |= (const SelfType& x)
    { low_ |= x.low_; high_ |= x.high_; return *this; }

    SelfType& operator &= (const SelfType& x)
    { low_ &= x.low_; high_ &= x.high_; return *this; }

    SelfType& operator ^= (const SelfType& x)
    { low_ ^= x.low_; high_ ^= x.high_; return *this; }

    SelfType& operator ++ ()
    { *this += 1; return *this; }

    SelfType operator ++ (int)
    { SelfType temp = *this; temp += 1; return temp; }

    SelfType& operator -- ()
    { *this -= 1; return *this; }

    SelfType operator -- (int)
    { SelfType temp = *this; temp -= 1; return temp; }

    SelfType operator ~ () const
    { SelfType temp( ~high_, ~low_); return temp; }

    SelfType& operator *= (const SelfType& x);

    SelfType& operator /= (const SelfType& x);

    SelfType& operator %= (const SelfType& x);

    SelfType& operator >>= (int n);

    SelfType& operator <<= (int n);

    bool operator == (const SelfType& x) const
    { return high_ == x.high_ and low_ == x.low_; }

    bool operator != (const SelfType& x) const
    { return not (*this == x); }

    bool operator < (const SelfType& x) const
    {
      return (high_ < x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) < HalfUnsigned(x.low_)));
    }

    bool operator > (const SelfType& x) const
    {
      return (high_ > x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) > HalfUnsigned(x.low_)));
    }

    bool operator <= (const SelfType& x) const
    {
      return (high_ < x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) <= HalfUnsigned(x.low_)));
    }

    bool operator >= (const SelfType& x) const
    {
      return (high_ > x.high_ or
              (high_ == x.high_ and HalfUnsigned(low_) >= HalfUnsigned(x.low_)));
    }

  protected:

    HalfType low_ = 0;
    HalfType high_ = 0;
  };


  inline
  Uint128::Uint128(const Int128& x)
    : low_(x.low()), high_(x.high())
  { }

  inline Uint128 operator + (Uint128 a, Uint128 b)
  { a += b; return a; }

  inline Uint128 operator - (Uint128 a, Uint128 b)
  { a -= b; return a; }

  inline Uint128 operator * (Uint128 a, Uint128 b)
  { a *= b; return a; }

  inline Uint128 operator / (Uint128 a, Uint128 b)
  { a /= b; return a; }

  inline Uint128 operator % (Uint128 a, Uint128 b)
  { a %= b; return a; }

  inline Uint128 operator - (Uint128 a)
  { Uint128 c = 0; c -= a; return c; }

  inline Uint128 operator >> (Uint128 x, int n)
  { x >>= n; return x; }

  inline Uint128 operator << (Uint128 x, int n)
  { x <<= n; return x; }

  inline Uint128 operator | (Uint128 a, Uint128 b)
  { a |= b; return a; }

  inline Uint128 operator & (Uint128 a, Uint128 b)
  { a &= b; return a; }

  inline Uint128 operator ^ (Uint128 a, Uint128 b)
  { a ^= b; return a; }


  inline Int128 operator + (Int128 a, Int128 b)
  { a += b; return a; }

  inline Int128 operator - (Int128 a, Int128 b)
  { a -= b; return a; }

  inline Int128 operator * (Int128 a, Int128 b)
  { a *= b; return a; }

  inline Int128 operator / (Int128 a, Int128 b)
  { a /= b; return a; }

  inline Int128 operator % (Int128 a, Int128 b)
  { a %= b; return a; }

  inline Int128 operator - (Int128 a)
  { Int128 c = 0; c -= a; return c; }

  inline Int128 operator >> (Int128 x, int n)
  { x >>= n; return x; }

  inline Int128 operator << (Int128 x, int n)
  { x <<= n; return x; }

  inline Int128 operator | (Int128 a, Int128 b)
  { a |= b; return a; }

  inline Int128 operator & (Int128 a, Int128 b)
  { a &= b; return a; }

  inline Int128 operator ^ (Int128 a, Int128 b)
  { a ^= b; return a; }



  inline
  Uint256::Uint256(const Int256& x)
    : low_(x.low()), high_(x.high())
  { }

  inline Uint256 operator + (Uint256 a, Uint256 b)
  { a += b; return a; }

  inline Uint256 operator - (Uint256 a, Uint256 b)
  { a -= b; return a; }

  inline Uint256 operator * (Uint256 a, Uint256 b)
  { a *= b; return a; }

  inline Uint256 operator / (Uint256 a, Uint256 b)
  { a /= b; return a; }

  inline Uint256 operator % (Uint256 a, Uint256 b)
  { a %= b; return a; }

  inline Uint256 operator - (Uint256 a)
  { Uint256 c = 0UL; c -= a; return c; }

  inline Uint256 operator >> (Uint256 x, int n)
  { x >>= n; return x; }

  inline Uint256 operator << (Uint256 x, int n)
  { x <<= n; return x; }

  inline Uint256 operator | (Uint256 a, Uint256 b)
  { a |= b; return a; }

  inline Uint256 operator & (Uint256 a, Uint256 b)
  { a &= b; return a; }

  inline Uint256 operator ^ (Uint256 a, Uint256 b)
  { a ^= b; return a; }


  inline Int256 operator + (Int256 a, Int256 b)
  { a += b; return a; }

  inline Int256 operator - (Int256 a, Int256 b)
  { a -= b; return a; }

  inline Int256 operator * (Int256 a, Int256 b)
  { a *= b; return a; }

  inline Int256 operator / (Int256 a, Int256 b)
  { a /= b; return a; }

  inline Int256 operator % (Int256 a, Int256 b)
  { a %= b; return a; }

  inline Int256 operator - (Int256 a)
  { Int256 c = 0L; c -= a; return c; }

  inline Int256 operator >> (Int256 x, int n)
  { x >>= n; return x; }

  inline Int256 operator << (Int256 x, int n)
  { x <<= n; return x; }

  inline Int256 operator | (Int256 a, Int256 b)
  { a |= b; return a; }

  inline Int256 operator & (Int256 a, Int256 b)
  { a &= b; return a; }

  inline Int256 operator ^ (Int256 a, Int256 b)
  { a ^= b; return a; }



  inline
  Uint512::Uint512(const Int512& x)
    : low_(x.low()), high_(x.high())
  { }

  inline Uint512 operator + (Uint512 a, Uint512 b)
  { a += b; return a; }

  inline Uint512 operator - (Uint512 a, Uint512 b)
  { a -= b; return a; }

  inline Uint512 operator * (Uint512 a, Uint512 b)
  { a *= b; return a; }

  inline Uint512 operator / (Uint512 a, Uint512 b)
  { a /= b; return a; }

  inline Uint512 operator % (Uint512 a, Uint512 b)
  { a %= b; return a; }

  inline Uint512 operator - (Uint512 a)
  { Uint512 c = 0UL; c -= a; return c; }

  inline Uint512 operator >> (Uint512 x, int n)
  { x >>= n; return x; }

  inline Uint512 operator << (Uint512 x, int n)
  { x <<= n; return x; }

  inline Uint512 operator | (Uint512 a, Uint512 b)
  { a |= b; return a; }

  inline Uint512 operator & (Uint512 a, Uint512 b)
  { a &= b; return a; }

  inline Uint512 operator ^ (Uint512 a, Uint512 b)
  { a ^= b; return a; }


  inline Int512 operator + (Int512 a, Int512 b)
  { a += b; return a; }

  inline Int512 operator - (Int512 a, Int512 b)
  { a -= b; return a; }

  inline Int512 operator * (Int512 a, Int512 b)
  { a *= b; return a; }

  inline Int512 operator / (Int512 a, Int512 b)
  { a /= b; return a; }

  inline Int512 operator % (Int512 a, Int512 b)
  { a %= b; return a; }

  inline Int512 operator - (Int512 a)
  { Int512 c = 0L; c -= a; return c; }

  inline Int512 operator >> (Int512 x, int n)
  { x >>= n; return x; }

  inline Int512 operator << (Int512 x, int n)
  { x <<= n; return x; }

  inline Int512 operator | (Int512 a, Int512 b)
  { a |= b; return a; }

  inline Int512 operator & (Int512 a, Int512 b)
  { a &= b; return a; }

  inline Int512 operator ^ (Int512 a, Int512 b)
  { a ^= b; return a; }



  inline
  Uint1024::Uint1024(const Int1024& x)
    : low_(x.low()), high_(x.high())
  { }

  inline Uint1024 operator + (Uint1024 a, Uint1024 b)
  { a += b; return a; }

  inline Uint1024 operator - (Uint1024 a, Uint1024 b)
  { a -= b; return a; }

  inline Uint1024 operator * (Uint1024 a, Uint1024 b)
  { a *= b; return a; }

  inline Uint1024 operator / (Uint1024 a, Uint1024 b)
  { a /= b; return a; }

  inline Uint1024 operator % (Uint1024 a, Uint1024 b)
  { a %= b; return a; }

  inline Uint1024 operator - (Uint1024 a)
  { Uint1024 c = 0UL; c -= a; return c; }

  inline Uint1024 operator >> (Uint1024 x, int n)
  { x >>= n; return x; }

  inline Uint1024 operator << (Uint1024 x, int n)
  { x <<= n; return x; }

  inline Uint1024 operator | (Uint1024 a, Uint1024 b)
  { a |= b; return a; }

  inline Uint1024 operator & (Uint1024 a, Uint1024 b)
  { a &= b; return a; }

  inline Uint1024 operator ^ (Uint1024 a, Uint1024 b)
  { a ^= b; return a; }


  inline Int1024 operator + (Int1024 a, Int1024 b)
  { a += b; return a; }

  inline Int1024 operator - (Int1024 a, Int1024 b)
  { a -= b; return a; }

  inline Int1024 operator * (Int1024 a, Int1024 b)
  { a *= b; return a; }

  inline Int1024 operator / (Int1024 a, Int1024 b)
  { a /= b; return a; }

  inline Int1024 operator % (Int1024 a, Int1024 b)
  { a %= b; return a; }

  inline Int1024 operator - (Int1024 a)
  { Int1024 c = 0L; c -= a; return c; }

  inline Int1024 operator >> (Int1024 x, int n)
  { x >>= n; return x; }

  inline Int1024 operator << (Int1024 x, int n)
  { x <<= n; return x; }

  inline Int1024 operator | (Int1024 a, Int1024 b)
  { a |= b; return a; }

  inline Int1024 operator & (Int1024 a, Int1024 b)
  { a &= b; return a; }

  inline Int1024 operator ^ (Int1024 a, Int1024 b)
  { a ^= b; return a; }

}
