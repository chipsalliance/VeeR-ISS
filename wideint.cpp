#include <cassert>
#include <iostream>
#include "wideint.hpp"

using namespace WdRiscv;


Uint128&
Uint128::operator *= (const Uint128& x)
{
  unsigned qw = width() / 4;

  QuarterType a0 = (low_ << qw) >> qw, a1 = low_ >> qw;
  QuarterType a2 = (high_ << qw) >> qw, a3 = high_ >> qw;
  QuarterType b0 = (x.low_ << qw) >> qw, b1 = x.low_ >> qw;
  QuarterType b2 = (x.high_ << qw) >> qw, b3 = x.high_ >> qw;

  *this = 0;
  SelfType shifted = 0;

  HalfType prod = a0; prod *= b0; *this += prod;
  prod = a0; prod *= b1; shifted = prod; shifted <<= qw; *this += shifted;
  prod = a0; prod *= b2; shifted = prod; shifted <<= 2*qw; *this += shifted;
  prod = a0; prod *= b3; shifted = prod; shifted <<= 3*qw; *this += shifted;

  prod = a1; prod *= b0; shifted = prod; shifted <<= qw; *this += shifted;
  prod = a1; prod *= b1; shifted = prod; shifted <<= 2*qw; *this += shifted;
  prod = a1; prod *= b2; shifted = prod; shifted <<= 3*qw; *this += shifted;
  prod = a1; prod *= b3; shifted = prod; shifted <<= 4*qw; *this += shifted;

  prod = a2; prod *= b0; shifted = prod; shifted <<= 2*qw; *this += shifted;
  prod = a2; prod *= b1; shifted = prod; shifted <<= 3*qw; *this += shifted;
  prod = a2; prod *= b2; shifted = prod; shifted <<= 4*qw; *this += shifted;
  prod = a2; prod *= b3; shifted = prod; shifted <<= 5*qw; *this += shifted;

  prod = a3; prod *= b0; shifted = prod; shifted <<= 3*qw; *this += shifted;
  prod = a3; prod *= b1; shifted = prod; shifted <<= 4*qw; *this += shifted;
  prod = a3; prod *= b2; shifted = prod; shifted <<= 5*qw; *this += shifted;
  prod = a3; prod *= b3; shifted = prod; shifted <<= 6*qw; *this += shifted;

  return *this;
}


Uint128&
Uint128::operator /= (const Uint128& x)
{
  unsigned n = width();
  SelfType rem(0), result(0);

  SelfType y = *this;  // Dividend

  uint8_t* remLow = reinterpret_cast<uint8_t*> (&rem);  // Least sig byte of rem
  uint8_t* resultLow = reinterpret_cast<uint8_t*> (&result);  // Least sig byte of result
  uint8_t* yHigh = (reinterpret_cast<uint8_t*> (&y)) + sizeof(y) - 1; // Most sig byte of dividend

  for (unsigned i = 0; i < n; ++i)
    {
      uint8_t yMsb = *yHigh >> 7; // Most sig bit of dividend
      rem <<= 1;
      result <<= 1;
      y <<= 1;
      *remLow |= yMsb;
      if (x <= rem)
        {
          *resultLow |= 1;
          rem -= x;
        }
    }

  *this = result;
  return *this;
}


Uint128&
Uint128::operator %= (const Uint128& x)
{
  unsigned n = width();
  SelfType rem(0), result(0);

  SelfType y = *this;  // Dividend

  uint8_t* remLow = reinterpret_cast<uint8_t*> (&rem);  // Least sig byte of rem
  uint8_t* resultLow = reinterpret_cast<uint8_t*> (&result);  // Least sig byte of result
  uint8_t* yHigh = (reinterpret_cast<uint8_t*> (&y)) + sizeof(y) - 1; // Most sig byte of dividend

  for (unsigned i = 0; i < n; ++i)
    {
      uint8_t yMsb = *yHigh >> 7; // Most sig bit of dividend
      rem <<= 1;
      result <<= 1;
      y <<= 1;
      *remLow |= yMsb;
      if (x <= rem)
        {
          *resultLow |= 1;
          rem -= x;
        }
    }

  *this = rem;
  return *this;
}


Uint128&
Uint128::operator >>= (int n)
{
  int halfw = halfWidth();

  if (n >= width())
    low_ = high_ = HalfType(0);
  else if (n >= halfw)
    {
      low_ = high_;
      high_ = HalfType(0);
      low_ >>= (n - halfw);
    }
  else
    {
      HalfType temp = high_ << (halfw - n);
      high_ >>= n;
      low_ >>= n;
      low_ |= temp;
    }
  return *this;
}


Uint128&
Uint128::operator <<= (int n)
{
  int halfw = halfWidth();

  if (n >= width())
    low_ = high_ = HalfType(0);
  else if (n >= halfw)
    {
      high_ = low_;
      low_ = HalfType(0);
      high_ <<= (n - halfw);
    }
  else
    {
      HalfType temp = low_ >> (halfw - n);
      high_ <<= n;
      low_ <<= n;
      high_ |= temp;
    }
  return *this;
}


//


Int128&
Int128::operator *= (const Int128& xx)
{
  if (*this >= SelfType(0) and xx >= SelfType(0))
    {
      UnsignedType a = *this, b = xx;
      a *= b;
      *this = a;
      return *this;
    }

  SelfType minInt(1);
  minInt <<= width() - 1;

  if (*this == minInt or xx == minInt)
    {
      *this = SelfType(0);
      return *this;
    }

  bool neg = false;

  SelfType b = xx;

  if (*this < SelfType(0) and xx >= SelfType(0))
    {
      neg = true;
      *this = - *this;
    }
  else if (*this >= SelfType(0) and xx < SelfType(0))
    {
      neg = true;
      b = -xx;
    }
  else
    {
      *this = - *this;
      b = -xx;
    }
      
  SelfType a = *this;
  a *= b;
  if (neg)
    a = -a;

  *this = a;

  return *this;
}


Int128&
Int128::operator /= (const Int128& xx)
{
  if (*this == xx)
    {
      *this = 1;
      return *this;
    }

  SelfType minInt(1);
  minInt <<= width() - 1;

  if (xx == minInt)
    {
      *this = 0;
      return *this;
    }

  bool neg = false;

  UnsignedType bb = xx;

  if (*this < 0 and xx >= 0)
    {
      neg = true;
      *this = - *this;
    }
  else if (*this >= 0 and xx < 0)
    {
      neg = true;
      bb = - xx;
    }
  else if (*this < 0 and xx < 0)
    {
      *this = - *this;
      bb = - xx;
    }

  UnsignedType aa = *this;
  aa /= bb;
  *this = aa;
  if (neg)
    *this = - *this;

  return *this;
}


Int128&
Int128::operator %= (const Int128& xx)
{
  if (*this == xx)
    {
      *this = 0;
      return *this;
    }

  SelfType minInt(1);
  minInt <<= width() - 1;

  if (xx == minInt)
    return *this;

  bool neg = false;
  UnsignedType bb = xx;

  if (*this < 0)
    {
      neg = true;
      *this = - *this;
    }
  if (xx < 0)
    bb = - xx;

  UnsignedType aa = *this;
  aa %= bb;
  *this = aa;
  if (neg)
    *this = - *this;

  return *this;
}


Int128&
Int128::operator >>= (int n)
{
  bool neg = high_ < HalfType(0);

  int halfw = halfWidth();

  if (n >= width())
    {
      if (neg)
        low_ = high_ = ~HalfType(0);
      else
        low_ = high_ = HalfType(0);
    }
  else if (n >= halfw)
    {
      low_ = high_;
      if (neg)
        low_ |= HalfType(1) << (sizeof(low_)*8 - 1);
      high_ = neg? ~HalfType(0) : HalfType(0);
      low_ >>= (n - halfw);
    }
  else
    {
      HalfUnsigned temp = high_ << (halfw - n);
      high_ >>= n;
      low_ = HalfUnsigned(low_) >> n;
      low_ |= temp;
    }

  return *this;
}


Int128&
Int128::operator <<= (int n)
{
  int halfw = halfWidth();

  if (n >= width())
    low_ = high_ = HalfType(0);
  else if (n >= halfw)
    {
      high_ = low_;
      low_ = HalfType(0);
      high_ <<= (n - halfw);
    }
  else
    {
      HalfUnsigned temp = low_ >> (halfw - n);
      high_ <<= n;
      low_ <<= n;
      high_ |= temp;
    }
  return *this;
}



// 256



Uint256&
Uint256::operator *= (const Uint256& x)
{
  unsigned qw = width() / 4;

  QuarterType a0 = QuarterType((low_ << qw) >> qw), a1 = QuarterType(low_ >> qw);
  QuarterType a2 = QuarterType((high_ << qw) >> qw), a3 = QuarterType(high_ >> qw);
  QuarterType b0 = QuarterType((x.low_ << qw) >> qw), b1 = QuarterType(x.low_ >> qw);
  QuarterType b2 = QuarterType((x.high_ << qw) >> qw), b3 = QuarterType(x.high_ >> qw);

  *this = 0;
  SelfType shifted = 0;

  HalfType prod = a0; prod *= b0; *this += prod;
  prod = a0; prod *= b1; shifted = prod; shifted <<= qw; *this += shifted;
  prod = a0; prod *= b2; shifted = prod; shifted <<= 2*qw; *this += shifted;
  prod = a0; prod *= b3; shifted = prod; shifted <<= 3*qw; *this += shifted;

  prod = a1; prod *= b0; shifted = prod; shifted <<= qw; *this += shifted;
  prod = a1; prod *= b1; shifted = prod; shifted <<= 2*qw; *this += shifted;
  prod = a1; prod *= b2; shifted = prod; shifted <<= 3*qw; *this += shifted;
  prod = a1; prod *= b3; shifted = prod; shifted <<= 4*qw; *this += shifted;

  prod = a2; prod *= b0; shifted = prod; shifted <<= 2*qw; *this += shifted;
  prod = a2; prod *= b1; shifted = prod; shifted <<= 3*qw; *this += shifted;
  prod = a2; prod *= b2; shifted = prod; shifted <<= 4*qw; *this += shifted;
  prod = a2; prod *= b3; shifted = prod; shifted <<= 5*qw; *this += shifted;

  prod = a3; prod *= b0; shifted = prod; shifted <<= 3*qw; *this += shifted;
  prod = a3; prod *= b1; shifted = prod; shifted <<= 4*qw; *this += shifted;
  prod = a3; prod *= b2; shifted = prod; shifted <<= 5*qw; *this += shifted;
  prod = a3; prod *= b3; shifted = prod; shifted <<= 6*qw; *this += shifted;

  return *this;
}


Uint256&
Uint256::operator /= (const Uint256& x)
{
  unsigned n = width();
  SelfType rem(0), result(0);

  SelfType y = *this;  // Dividend

  uint8_t* remLow = reinterpret_cast<uint8_t*> (&rem);  // Least sig byte of rem
  uint8_t* resultLow = reinterpret_cast<uint8_t*> (&result);  // Least sig byte of result
  uint8_t* yHigh = (reinterpret_cast<uint8_t*> (&y)) + sizeof(y) - 1; // Most sig byte of dividend

  for (unsigned i = 0; i < n; ++i)
    {
      uint8_t yMsb = *yHigh >> 7; // Most sig bit of dividend
      rem <<= 1;
      result <<= 1;
      y <<= 1;
      *remLow |= yMsb;
      if (x <= rem)
        {
          *resultLow |= 1;
          rem -= x;
        }
    }

  *this = result;
  return *this;
}


Uint256&
Uint256::operator %= (const Uint256& x)
{
  unsigned n = width();
  SelfType rem(0), result(0);

  SelfType y = *this;  // Dividend

  uint8_t* remLow = reinterpret_cast<uint8_t*> (&rem);  // Least sig byte of rem
  uint8_t* resultLow = reinterpret_cast<uint8_t*> (&result);  // Least sig byte of result
  uint8_t* yHigh = (reinterpret_cast<uint8_t*> (&y)) + sizeof(y) - 1; // Most sig byte of dividend

  for (unsigned i = 0; i < n; ++i)
    {
      uint8_t yMsb = *yHigh >> 7; // Most sig bit of dividend
      rem <<= 1;
      result <<= 1;
      y <<= 1;
      *remLow |= yMsb;
      if (x <= rem)
        {
          *resultLow |= 1;
          rem -= x;
        }
    }

  *this = rem;
  return *this;
}


Uint256&
Uint256::operator >>= (int n)
{
  int halfw = halfWidth();

  if (n >= width())
    low_ = high_ = HalfType(0);
  else if (n >= halfw)
    {
      low_ = high_;
      high_ = HalfType(0);
      low_ >>= (n - halfw);
    }
  else
    {
      HalfType temp = high_ << (halfw - n);
      high_ >>= n;
      low_ >>= n;
      low_ |= temp;
    }
  return *this;
}


Uint256&
Uint256::operator <<= (int n)
{
  int halfw = halfWidth();

  if (n >= width())
    low_ = high_ = HalfType(0);
  else if (n >= halfw)
    {
      high_ = low_;
      low_ = HalfType(0);
      high_ <<= (n - halfw);
    }
  else
    {
      HalfType temp = low_ >> (halfw - n);
      high_ <<= n;
      low_ <<= n;
      high_ |= temp;
    }
  return *this;
}


//


Int256&
Int256::operator *= (const Int256& xx)
{
  if (*this >= SelfType(0) and xx >= SelfType(0))
    {
      UnsignedType a = *this, b = xx;
      a *= b;
      *this = a;
      return *this;
    }

  SelfType minInt(1);
  minInt <<= width() - 1;

  if (*this == minInt or xx == minInt)
    {
      *this = SelfType(0);
      return *this;
    }

  bool neg = false;

  SelfType b = xx;

  if (*this < SelfType(0) and xx >= SelfType(0))
    {
      neg = true;
      *this = - *this;
    }
  else if (*this >= SelfType(0) and xx < SelfType(0))
    {
      neg = true;
      b = -xx;
    }
  else
    {
      *this = - *this;
      b = -xx;
    }
      
  SelfType a = *this;
  a *= b;
  if (neg)
    a = -a;

  *this = a;

  return *this;
}


Int256&
Int256::operator /= (const Int256& xx)
{
  if (*this == xx)
    {
      *this = 1;
      return *this;
    }

  SelfType minInt(1);
  minInt <<= width() - 1;

  if (xx == minInt)
    {
      *this = 0;
      return *this;
    }

  bool neg = false;

  UnsignedType bb = xx;

  if (*this < 0 and xx >= 0)
    {
      neg = true;
      *this = - *this;
    }
  else if (*this >= 0 and xx < 0)
    {
      neg = true;
      bb = - xx;
    }
  else if (*this < 0 and xx < 0)
    {
      *this = - *this;
      bb = - xx;
    }

  UnsignedType aa = *this;
  aa /= bb;
  *this = aa;
  if (neg)
    *this = - *this;

  return *this;
}


Int256&
Int256::operator %= (const Int256& xx)
{
  if (*this == xx)
    {
      *this = 0;
      return *this;
    }

  SelfType minInt(1);
  minInt <<= width() - 1;

  if (xx == minInt)
    return *this;

  bool neg = false;
  UnsignedType bb = xx;

  if (*this < 0)
    {
      neg = true;
      *this = - *this;
    }
  if (xx < 0)
    bb = - xx;

  UnsignedType aa = *this;
  aa %= bb;
  *this = aa;
  if (neg)
    *this = - *this;

  return *this;
}


Int256&
Int256::operator >>= (int n)
{
  bool neg = high_ < HalfType(0);

  int halfw = halfWidth();

  if (n >= width())
    {
      if (neg)
        low_ = high_ = ~HalfType(0);
      else
        low_ = high_ = HalfType(0);
    }
  else if (n >= halfw)
    {
      low_ = high_;
      if (neg)
        low_ |= HalfType(1) << int(sizeof(low_)*8 - 1);
      high_ = neg? ~HalfType(0) : HalfType(0);
      low_ >>= (n - halfw);
    }
  else
    {
      HalfUnsigned temp = high_ << (halfw - n);
      high_ >>= n;
      low_ = HalfUnsigned(low_) >> n;
      low_ |= temp;
    }

  return *this;
}


Int256&
Int256::operator <<= (int n)
{
  int halfw = halfWidth();

  if (n >= width())
    low_ = high_ = HalfType(0);
  else if (n >= halfw)
    {
      high_ = low_;
      low_ = HalfType(0);
      high_ <<= (n - halfw);
    }
  else
    {
      HalfUnsigned temp = low_ >> (halfw - n);
      high_ <<= n;
      low_ <<= n;
      high_ |= temp;
    }
  return *this;
}



// 512



Uint512&
Uint512::operator *= (const Uint512& x)
{
  unsigned qw = width() / 4;

  QuarterType a0 = (low_ << qw) >> qw, a1 = low_ >> qw;
  QuarterType a2 = (high_ << qw) >> qw, a3 = high_ >> qw;
  QuarterType b0 = (x.low_ << qw) >> qw, b1 = x.low_ >> qw;
  QuarterType b2 = (x.high_ << qw) >> qw, b3 = x.high_ >> qw;

  *this = 0;
  SelfType shifted = 0;

  HalfType prod = a0; prod *= b0; *this += prod;
  prod = a0; prod *= b1; shifted = prod; shifted <<= qw; *this += shifted;
  prod = a0; prod *= b2; shifted = prod; shifted <<= 2*qw; *this += shifted;
  prod = a0; prod *= b3; shifted = prod; shifted <<= 3*qw; *this += shifted;

  prod = a1; prod *= b0; shifted = prod; shifted <<= qw; *this += shifted;
  prod = a1; prod *= b1; shifted = prod; shifted <<= 2*qw; *this += shifted;
  prod = a1; prod *= b2; shifted = prod; shifted <<= 3*qw; *this += shifted;
  prod = a1; prod *= b3; shifted = prod; shifted <<= 4*qw; *this += shifted;

  prod = a2; prod *= b0; shifted = prod; shifted <<= 2*qw; *this += shifted;
  prod = a2; prod *= b1; shifted = prod; shifted <<= 3*qw; *this += shifted;
  prod = a2; prod *= b2; shifted = prod; shifted <<= 4*qw; *this += shifted;
  prod = a2; prod *= b3; shifted = prod; shifted <<= 5*qw; *this += shifted;

  prod = a3; prod *= b0; shifted = prod; shifted <<= 3*qw; *this += shifted;
  prod = a3; prod *= b1; shifted = prod; shifted <<= 4*qw; *this += shifted;
  prod = a3; prod *= b2; shifted = prod; shifted <<= 5*qw; *this += shifted;
  prod = a3; prod *= b3; shifted = prod; shifted <<= 6*qw; *this += shifted;

  return *this;
}


Uint512&
Uint512::operator /= (const Uint512& x)
{
  unsigned n = width();
  SelfType rem(0), result(0);

  SelfType y = *this;  // Dividend

  uint8_t* remLow = reinterpret_cast<uint8_t*> (&rem);  // Least sig byte of rem
  uint8_t* resultLow = reinterpret_cast<uint8_t*> (&result);  // Least sig byte of result
  uint8_t* yHigh = (reinterpret_cast<uint8_t*> (&y)) + sizeof(y) - 1; // Most sig byte of dividend

  for (unsigned i = 0; i < n; ++i)
    {
      uint8_t yMsb = *yHigh >> 7; // Most sig bit of dividend
      rem <<= 1;
      result <<= 1;
      y <<= 1;
      *remLow |= yMsb;
      if (x <= rem)
        {
          *resultLow |= 1;
          rem -= x;
        }
    }

  *this = result;
  return *this;
}


Uint512&
Uint512::operator %= (const Uint512& x)
{
  unsigned n = width();
  SelfType rem(0), result(0);

  SelfType y = *this;  // Dividend

  uint8_t* remLow = reinterpret_cast<uint8_t*> (&rem);  // Least sig byte of rem
  uint8_t* resultLow = reinterpret_cast<uint8_t*> (&result);  // Least sig byte of result
  uint8_t* yHigh = (reinterpret_cast<uint8_t*> (&y)) + sizeof(y) - 1; // Most sig byte of dividend

  for (unsigned i = 0; i < n; ++i)
    {
      uint8_t yMsb = *yHigh >> 7; // Most sig bit of dividend
      rem <<= 1;
      result <<= 1;
      y <<= 1;
      *remLow |= yMsb;
      if (x <= rem)
        {
          *resultLow |= 1;
          rem -= x;
        }
    }

  *this = rem;
  return *this;
}


Uint512&
Uint512::operator >>= (int n)
{
  int halfw = halfWidth();

  if (n >= width())
    low_ = high_ = HalfType(0);
  else if (n >= halfw)
    {
      low_ = high_;
      high_ = HalfType(0);
      low_ >>= (n - halfw);
    }
  else
    {
      HalfType temp = high_ << (halfw - n);
      high_ >>= n;
      low_ >>= n;
      low_ |= temp;
    }
  return *this;
}


Uint512&
Uint512::operator <<= (int n)
{
  int halfw = halfWidth();

  if (n >= width())
    low_ = high_ = HalfType(0);
  else if (n >= halfw)
    {
      high_ = low_;
      low_ = HalfType(0);
      high_ <<= (n - halfw);
    }
  else
    {
      HalfType temp = low_ >> (halfw - n);
      high_ <<= n;
      low_ <<= n;
      high_ |= temp;
    }
  return *this;
}


//


Int512&
Int512::operator *= (const Int512& xx)
{
  if (*this >= SelfType(0) and xx >= SelfType(0))
    {
      UnsignedType a = *this, b = xx;
      a *= b;
      *this = a;
      return *this;
    }

  SelfType minInt(1);
  minInt <<= width() - 1;

  if (*this == minInt or xx == minInt)
    {
      *this = SelfType(0);
      return *this;
    }

  bool neg = false;

  SelfType b = xx;

  if (*this < SelfType(0) and xx >= SelfType(0))
    {
      neg = true;
      *this = - *this;
    }
  else if (*this >= SelfType(0) and xx < SelfType(0))
    {
      neg = true;
      b = -xx;
    }
  else
    {
      *this = - *this;
      b = -xx;
    }
      
  SelfType a = *this;
  a *= b;
  if (neg)
    a = -a;

  *this = a;

  return *this;
}


Int512&
Int512::operator /= (const Int512& xx)
{
  if (*this == xx)
    {
      *this = 1;
      return *this;
    }

  SelfType minInt(1);
  minInt <<= width() - 1;

  if (xx == minInt)
    {
      *this = 0;
      return *this;
    }

  bool neg = false;

  UnsignedType bb = xx;

  if (*this < 0 and xx >= 0)
    {
      neg = true;
      *this = - *this;
    }
  else if (*this >= 0 and xx < 0)
    {
      neg = true;
      bb = - xx;
    }
  else if (*this < 0 and xx < 0)
    {
      *this = - *this;
      bb = - xx;
    }

  UnsignedType aa = *this;
  aa /= bb;
  *this = aa;
  if (neg)
    *this = - *this;

  return *this;
}


Int512&
Int512::operator %= (const Int512& xx)
{
  if (*this == xx)
    {
      *this = 0;
      return *this;
    }

  SelfType minInt(1);
  minInt <<= width() - 1;

  if (xx == minInt)
    return *this;

  bool neg = false;
  UnsignedType bb = xx;

  if (*this < 0)
    {
      neg = true;
      *this = - *this;
    }
  if (xx < 0)
    bb = - xx;

  UnsignedType aa = *this;
  aa %= bb;
  *this = aa;
  if (neg)
    *this = - *this;

  return *this;
}


Int512&
Int512::operator >>= (int n)
{
  bool neg = high_ < HalfType(0);

  int halfw = halfWidth();

  if (n >= width())
    {
      if (neg)
        low_ = high_ = ~HalfType(0);
      else
        low_ = high_ = HalfType(0);
    }
  else if (n >= halfw)
    {
      low_ = high_;
      if (neg)
        low_ |= HalfType(1) << int(sizeof(low_)*8 - 1);
      high_ = neg? ~HalfType(0) : HalfType(0);
      low_ >>= (n - halfw);
    }
  else
    {
      HalfUnsigned temp = high_ << (halfw - n);
      high_ >>= n;
      low_ = HalfUnsigned(low_) >> n;
      low_ |= temp;
    }

  return *this;
}


Int512&
Int512::operator <<= (int n)
{
  int halfw = halfWidth();

  if (n >= width())
    low_ = high_ = HalfType(0);
  else if (n >= halfw)
    {
      high_ = low_;
      low_ = HalfType(0);
      high_ <<= (n - halfw);
    }
  else
    {
      HalfUnsigned temp = low_ >> (halfw - n);
      high_ <<= n;
      low_ <<= n;
      high_ |= temp;
    }
  return *this;
}



// 1024



Uint1024&
Uint1024::operator *= (const Uint1024& x)
{
  unsigned qw = width() / 4;

  QuarterType a0 = (low_ << qw) >> qw, a1 = low_ >> qw;
  QuarterType a2 = (high_ << qw) >> qw, a3 = high_ >> qw;
  QuarterType b0 = (x.low_ << qw) >> qw, b1 = x.low_ >> qw;
  QuarterType b2 = (x.high_ << qw) >> qw, b3 = x.high_ >> qw;

  *this = 0;
  SelfType shifted = 0;

  HalfType prod = a0; prod *= b0; *this += prod;
  prod = a0; prod *= b1; shifted = prod; shifted <<= qw; *this += shifted;
  prod = a0; prod *= b2; shifted = prod; shifted <<= 2*qw; *this += shifted;
  prod = a0; prod *= b3; shifted = prod; shifted <<= 3*qw; *this += shifted;

  prod = a1; prod *= b0; shifted = prod; shifted <<= qw; *this += shifted;
  prod = a1; prod *= b1; shifted = prod; shifted <<= 2*qw; *this += shifted;
  prod = a1; prod *= b2; shifted = prod; shifted <<= 3*qw; *this += shifted;
  prod = a1; prod *= b3; shifted = prod; shifted <<= 4*qw; *this += shifted;

  prod = a2; prod *= b0; shifted = prod; shifted <<= 2*qw; *this += shifted;
  prod = a2; prod *= b1; shifted = prod; shifted <<= 3*qw; *this += shifted;
  prod = a2; prod *= b2; shifted = prod; shifted <<= 4*qw; *this += shifted;
  prod = a2; prod *= b3; shifted = prod; shifted <<= 5*qw; *this += shifted;

  prod = a3; prod *= b0; shifted = prod; shifted <<= 3*qw; *this += shifted;
  prod = a3; prod *= b1; shifted = prod; shifted <<= 4*qw; *this += shifted;
  prod = a3; prod *= b2; shifted = prod; shifted <<= 5*qw; *this += shifted;
  prod = a3; prod *= b3; shifted = prod; shifted <<= 6*qw; *this += shifted;

  return *this;
}


Uint1024&
Uint1024::operator /= (const Uint1024& x)
{
  unsigned n = width();
  SelfType rem(0), result(0);

  SelfType y = *this;  // Dividend

  uint8_t* remLow = reinterpret_cast<uint8_t*> (&rem);  // Least sig byte of rem
  uint8_t* resultLow = reinterpret_cast<uint8_t*> (&result);  // Least sig byte of result
  uint8_t* yHigh = (reinterpret_cast<uint8_t*> (&y)) + sizeof(y) - 1; // Most sig byte of dividend

  for (unsigned i = 0; i < n; ++i)
    {
      uint8_t yMsb = *yHigh >> 7; // Most sig bit of dividend
      rem <<= 1;
      result <<= 1;
      y <<= 1;
      *remLow |= yMsb;
      if (x <= rem)
        {
          *resultLow |= 1;
          rem -= x;
        }
    }

  *this = result;
  return *this;
}


Uint1024&
Uint1024::operator %= (const Uint1024& x)
{
  unsigned n = width();
  SelfType rem(0), result(0);

  SelfType y = *this;  // Dividend

  uint8_t* remLow = reinterpret_cast<uint8_t*> (&rem);  // Least sig byte of rem
  uint8_t* resultLow = reinterpret_cast<uint8_t*> (&result);  // Least sig byte of result
  uint8_t* yHigh = (reinterpret_cast<uint8_t*> (&y)) + sizeof(y) - 1; // Most sig byte of dividend

  for (unsigned i = 0; i < n; ++i)
    {
      uint8_t yMsb = *yHigh >> 7; // Most sig bit of dividend
      rem <<= 1;
      result <<= 1;
      y <<= 1;
      *remLow |= yMsb;
      if (x <= rem)
        {
          *resultLow |= 1;
          rem -= x;
        }
    }

  *this = rem;
  return *this;
}


Uint1024&
Uint1024::operator >>= (int n)
{
  int halfw = halfWidth();

  if (n >= width())
    low_ = high_ = HalfType(0);
  else if (n >= halfw)
    {
      low_ = high_;
      high_ = HalfType(0);
      low_ >>= (n - halfw);
    }
  else
    {
      HalfType temp = high_ << (halfw - n);
      high_ >>= n;
      low_ >>= n;
      low_ |= temp;
    }
  return *this;
}


Uint1024&
Uint1024::operator <<= (int n)
{
  int halfw = halfWidth();

  if (n >= width())
    low_ = high_ = HalfType(0);
  else if (n >= halfw)
    {
      high_ = low_;
      low_ = HalfType(0);
      high_ <<= (n - halfw);
    }
  else
    {
      HalfType temp = low_ >> (halfw - n);
      high_ <<= n;
      low_ <<= n;
      high_ |= temp;
    }
  return *this;
}


//


Int1024&
Int1024::operator *= (const Int1024& xx)
{
  if (*this >= SelfType(0) and xx >= SelfType(0))
    {
      UnsignedType a = *this, b = xx;
      a *= b;
      *this = a;
      return *this;
    }

  SelfType minInt(1);
  minInt <<= width() - 1;

  if (*this == minInt or xx == minInt)
    {
      *this = SelfType(0);
      return *this;
    }

  bool neg = false;

  SelfType b = xx;

  if (*this < SelfType(0) and xx >= SelfType(0))
    {
      neg = true;
      *this = - *this;
    }
  else if (*this >= SelfType(0) and xx < SelfType(0))
    {
      neg = true;
      b = -xx;
    }
  else
    {
      *this = - *this;
      b = -xx;
    }
      
  SelfType a = *this;
  a *= b;
  if (neg)
    a = -a;

  *this = a;

  return *this;
}


Int1024&
Int1024::operator /= (const Int1024& xx)
{
  if (*this == xx)
    {
      *this = 1;
      return *this;
    }

  SelfType minInt(1);
  minInt <<= width() - 1;

  if (xx == minInt)
    {
      *this = 0;
      return *this;
    }

  bool neg = false;

  UnsignedType bb = xx;

  if (*this < 0 and xx >= 0)
    {
      neg = true;
      *this = - *this;
    }
  else if (*this >= 0 and xx < 0)
    {
      neg = true;
      bb = - xx;
    }
  else if (*this < 0 and xx < 0)
    {
      *this = - *this;
      bb = - xx;
    }

  UnsignedType aa = *this;
  aa /= bb;
  *this = aa;
  if (neg)
    *this = - *this;

  return *this;
}


Int1024&
Int1024::operator %= (const Int1024& xx)
{
  if (*this == xx)
    {
      *this = 0;
      return *this;
    }

  SelfType minInt(1);
  minInt <<= width() - 1;

  if (xx == minInt)
    return *this;

  bool neg = false;
  UnsignedType bb = xx;

  if (*this < 0)
    {
      neg = true;
      *this = - *this;
    }
  if (xx < 0)
    bb = - xx;

  UnsignedType aa = *this;
  aa %= bb;
  *this = aa;
  if (neg)
    *this = - *this;

  return *this;
}


Int1024&
Int1024::operator >>= (int n)
{
  bool neg = high_ < HalfType(0);

  int halfw = halfWidth();

  if (n >= width())
    {
      if (neg)
        low_ = high_ = ~HalfType(0);
      else
        low_ = high_ = HalfType(0);
    }
  else if (n >= halfw)
    {
      low_ = high_;
      if (neg)
        low_ |= HalfType(1) << int(sizeof(low_)*8 - 1);
      high_ = neg? ~HalfType(0) : HalfType(0);
      low_ >>= (n - halfw);
    }
  else
    {
      HalfUnsigned temp = high_ << (halfw - n);
      high_ >>= n;
      low_ = HalfUnsigned(low_) >> n;
      low_ |= temp;
    }

  return *this;
}


Int1024&
Int1024::operator <<= (int n)
{
  int halfw = halfWidth();

  if (n >= width())
    low_ = high_ = HalfType(0);
  else if (n >= halfw)
    {
      high_ = low_;
      low_ = HalfType(0);
      high_ <<= (n - halfw);
    }
  else
    {
      HalfUnsigned temp = low_ >> (halfw - n);
      high_ <<= n;
      low_ <<= n;
      high_ |= temp;
    }
  return *this;
}


# if 0

#include <random>

int
main(int argc, char* argv[])
{
  int64_t i64 = uint32_t(0xffffffff);
  Int128 i128 = uint64_t(0xffffffffffffffffLL);
  Uint128 u128 = 0; u128 = ~u128;
  unsigned bit = unsigned((u128 >> 127) & 1);
  u128 = u128 >> 127;
  Int256 i256 = u128;
  i256 = -1;
  Int512 i512{i256};
  bool d = i512 != 0;
  d = u128 != 0;
  
  uint64_t a = ~uint64_t(0);
  __uint128_t c = ~ __uint128_t(0);
  c *= a;

  Uint128 cc(a, a);
  cc *= a;

  assert(cc.high() == uint64_t(c >> 64));
  assert(cc.low() == (c << 64) >> 64);

  c = ~ __uint128_t(0);
  c *= c;

  cc = Uint128(a, a);
  cc *= cc;

  if ((cc.high() != c >> 64) or (cc.low() != (c << 64) >> 64))
    std::cerr << "Failed\n";

  std::mt19937 gen;
  std::uniform_int_distribution<uint64_t> distrib(0, ~uint64_t(0));

  for (unsigned i = 0; i < 100000000; ++i)
    {
      if ((i % 1000000) == 0)
        std::cerr << std::dec << i << '\n';

      uint64_t low1 = distrib(gen), high1 = distrib(gen);
      __uint128_t x1 = high1;
      x1 = (x1 << 64) | low1;
      Uint128 y1(high1, low1);

      uint64_t low2 = distrib(gen), high2 = distrib(gen);
      __uint128_t x2 = high2;
      x2 = (x2 << 64) | low2;
      Uint128 y2(high2, low2);

      __uint128_t z1 = x1 * x2;
      Uint128 z2 = y1 * y2;

      if ((z2.high() != z1 >> 64) or (z2.low() != (z1 << 64) >> 64))
        {
          std::cerr << std::dec << "Failed mulu " << i << std::hex << " on 0x"
                    << high1 << " 0x" << low1
                    << " * 0x" << high2 << " 0x " << low2 << '\n';
          return 1;
        }

      __int128_t sx1 = x1, sx2 = x2;
      __int128_t sz1 = sx1 * sx2;

      Int128 sy1 = y1, sy2 = y2;
      Int128 sz2 = sy1 * sy2;

      if ((sz2.high() != sz1 >> 64 or z2.low() != (z1 << 64) >> 64))
        {
          std::cerr << std::dec << "Failed mul " << i << std::hex << " on 0x"
                    << high1 << " 0x" << low1
                    << " * 0x" << high2 << " 0x " << low2 << '\n';
          return 1;
        }

      if (x2 == 0)
        continue;

      z1 = x1 / x2;
      z2 = y1 / y2;

      if ((z2.high() != z1 >> 64) or (z2.low() != (z1 << 64) >> 64))
        {
          std::cerr << std::dec << "Failed divu" << i << std::hex << " on 0x"
                    << high1 << " 0x" << low1
                    << " * 0x" << high2 << " 0x " << low2 << '\n';
          return 1;
        }

      sz1 = sx1 / sx2;
      sz2 = sy1 / sy2;

      if ((sz2.high() != sz1 >> 64) or (sz2.low() != (sz1 << 64) >> 64))
        {
          std::cerr << std::dec << "Failed div" << i << std::hex << " on 0x"
                    << high1 << " 0x" << low1
                    << " * 0x" << high2 << " 0x " << low2 << '\n';
          return 1;
        }

    }      

  return 0;
}

#endif
