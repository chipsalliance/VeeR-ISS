#include <cassert>
#include "wideint.hpp"

using namespace WdRiscv;


Uint128&
Uint128::operator *= (const Uint128& x)
{
  QuarterType* multiplicand = (QuarterType*) &low_;
  QuarterType* multiplier = (QuarterType*) &x.low_;

  SelfType prod[2];
  QuarterType* acc = (QuarterType*) &prod[0];

  for (unsigned i = 0; i < 4; ++i, ++acc)
    {
      QuarterType a = multiplier[i];
      QuarterType carry = 0;
      QuarterType* row = acc;
      for (unsigned j = 0; j < 4; ++j, ++row)
        {
          HalfType b = multiplicand[j];
          b *= a;
          HalfType* sum = (HalfType*) row;
          HalfType prev = *sum;
          *sum += b + HalfType(carry);
          carry = (*sum < prev)? QuarterType(1) : QuarterType(0);
        }
    }
  *this = prod[0];
  return *this;
}


Uint128&
Uint128::operator /= (const Uint128& x)
{
  unsigned n = width();
  SelfType rem(0), result(0);

  SelfType y = *this;  // Dividend

  uint8_t* remLow = (uint8_t*) &rem;  // Least sig byte of rem
  uint8_t* resultLow = (uint8_t*) &result;  // Least sig byte of result
  uint8_t* yHigh = ((uint8_t*) &y) + sizeof(y) - 1; // Most sig byte of dividend

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

  uint8_t* remLow = (uint8_t*) &rem;  // Least sig byte of rem
  uint8_t* resultLow = (uint8_t*) &result;  // Least sig byte of result
  uint8_t* yHigh = ((uint8_t*) &y) + sizeof(y) - 1; // Most sig byte of dividend

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
  QuarterType* multiplicand = (QuarterType*) &low_;
  QuarterType* multiplier = (QuarterType*) &x.low_;

  SelfType prod[2];
  QuarterType* acc = (QuarterType*) &prod[0];

  for (unsigned i = 0; i < 4; ++i, ++acc)
    {
      QuarterType a = multiplier[i];
      QuarterType carry = 0;
      QuarterType* row = acc;
      for (unsigned j = 0; j < 4; ++j, ++row)
        {
          HalfType b = multiplicand[j];
          b *= a;
          HalfType* sum = (HalfType*) row;
          HalfType prev = *sum;
          *sum += b + HalfType(carry);
          carry = (*sum < prev)? QuarterType(1) : QuarterType(0);
        }
    }
  *this = prod[0];
  return *this;
}


Uint256&
Uint256::operator /= (const Uint256& x)
{
  unsigned n = width();
  SelfType rem(0), result(0);

  SelfType y = *this;  // Dividend

  uint8_t* remLow = (uint8_t*) &rem;  // Least sig byte of rem
  uint8_t* resultLow = (uint8_t*) &result;  // Least sig byte of result
  uint8_t* yHigh = ((uint8_t*) &y) + sizeof(y) - 1; // Most sig byte of dividend

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

  uint8_t* remLow = (uint8_t*) &rem;  // Least sig byte of rem
  uint8_t* resultLow = (uint8_t*) &result;  // Least sig byte of result
  uint8_t* yHigh = ((uint8_t*) &y) + sizeof(y) - 1; // Most sig byte of dividend

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
  QuarterType* multiplicand = (QuarterType*) &low_;
  QuarterType* multiplier = (QuarterType*) &x.low_;

  SelfType prod[2];
  QuarterType* acc = (QuarterType*) &prod[0];

  for (unsigned i = 0; i < 4; ++i, ++acc)
    {
      QuarterType a = multiplier[i];
      QuarterType carry = 0;
      QuarterType* row = acc;
      for (unsigned j = 0; j < 4; ++j, ++row)
        {
          HalfType b = multiplicand[j];
          b *= a;
          HalfType* sum = (HalfType*) row;
          HalfType prev = *sum;
          *sum += b + HalfType(carry);
          carry = (*sum < prev)? QuarterType(1) : QuarterType(0);
        }
    }
  *this = prod[0];
  return *this;
}


Uint512&
Uint512::operator /= (const Uint512& x)
{
  unsigned n = width();
  SelfType rem(0), result(0);

  SelfType y = *this;  // Dividend

  uint8_t* remLow = (uint8_t*) &rem;  // Least sig byte of rem
  uint8_t* resultLow = (uint8_t*) &result;  // Least sig byte of result
  uint8_t* yHigh = ((uint8_t*) &y) + sizeof(y) - 1; // Most sig byte of dividend

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

  uint8_t* remLow = (uint8_t*) &rem;  // Least sig byte of rem
  uint8_t* resultLow = (uint8_t*) &result;  // Least sig byte of result
  uint8_t* yHigh = ((uint8_t*) &y) + sizeof(y) - 1; // Most sig byte of dividend

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
  QuarterType* multiplicand = (QuarterType*) &low_;
  QuarterType* multiplier = (QuarterType*) &x.low_;

  SelfType prod[2];
  QuarterType* acc = (QuarterType*) &prod[0];

  for (unsigned i = 0; i < 4; ++i, ++acc)
    {
      QuarterType a = multiplier[i];
      QuarterType carry = 0;
      QuarterType* row = acc;
      for (unsigned j = 0; j < 4; ++j, ++row)
        {
          HalfType b = multiplicand[j];
          b *= a;
          HalfType* sum = (HalfType*) row;
          HalfType prev = *sum;
          *sum += b + HalfType(carry);
          carry = (*sum < prev)? QuarterType(1) : QuarterType(0);
        }
    }
  *this = prod[0];
  return *this;
}


Uint1024&
Uint1024::operator /= (const Uint1024& x)
{
  unsigned n = width();
  SelfType rem(0), result(0);

  SelfType y = *this;  // Dividend

  uint8_t* remLow = (uint8_t*) &rem;  // Least sig byte of rem
  uint8_t* resultLow = (uint8_t*) &result;  // Least sig byte of result
  uint8_t* yHigh = ((uint8_t*) &y) + sizeof(y) - 1; // Most sig byte of dividend

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

  uint8_t* remLow = (uint8_t*) &rem;  // Least sig byte of rem
  uint8_t* resultLow = (uint8_t*) &result;  // Least sig byte of result
  uint8_t* yHigh = ((uint8_t*) &y) + sizeof(y) - 1; // Most sig byte of dividend

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
