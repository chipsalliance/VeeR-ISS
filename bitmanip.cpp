#include "Hart.hpp"

using namespace WdRiscv;


template <typename URV>
void
Hart<URV>::execClz(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());

  if (v1 == 0)
    v1 = mxlen_;
  else
    {
      if (mxlen_ == 32)
        v1 = __builtin_clz(v1);
      else
        v1 = __builtin_clzl(v1);
    }

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execCtz(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());

  if (mxlen_ == 32)
    v1 = __builtin_ctz(v1);
  else
    v1 = __builtin_ctzl(v1);

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execCpop(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());

  if (mxlen_ == 32)
    v1 = __builtin_popcount(v1);
  else
    v1 = __builtin_popcountl(v1);

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execClzw(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  uint32_t v1 = intRegs_.read(di->op1());

  if (v1 == 0)
    v1 = 32;
  else
    v1 = __builtin_clz(v1);

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execCtzw(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  uint32_t v1 = intRegs_.read(di->op1());
  v1 = __builtin_ctz(v1);

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execCpopw(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  uint32_t v1 = intRegs_.read(di->op1());
  URV res = __builtin_popcount(v1);
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execAndn(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV res = v1 & ~v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execOrn(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV res = v1 | ~v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execXnor(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV res = v1 ^ ~v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execMin(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  SRV v1 = intRegs_.read(di->op1());
  SRV v2 = intRegs_.read(di->op2());
  SRV res = v1 < v2? v1 : v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execMax(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  SRV v1 = intRegs_.read(di->op1());
  SRV v2 = intRegs_.read(di->op2());
  SRV res = v1 > v2? v1 : v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execMinu(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV res = v1 < v2? v1 : v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execMaxu(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV res = v1 > v2? v1 : v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execRol(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV mask = shiftMask();
  URV rot = intRegs_.read(di->op2()) & mask;  // Rotate amount

  URV v1 = intRegs_.read(di->op1());
  URV res = (v1 << rot) | (v1 >> ((mxlen_ - rot) & mask));

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execRor(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV mask = shiftMask();
  URV rot = intRegs_.read(di->op2()) & mask;  // Rotate amount

  URV v1 = intRegs_.read(di->op1());
  URV res = (v1 >> rot) | (v1 << ((mxlen_ - rot) & mask));

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execRori(const DecodedInst* di)
{
  if (not isRvzbb() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV rot = di->op2();
  if (not checkShiftImmediate(di, rot))
    return;

  URV v1 = intRegs_.read(di->op1());
  URV res = (v1 >> rot) | (v1 << ((mxlen_ - rot) & shiftMask()));

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execRolw(const DecodedInst* di)
{
  if (not isRv64() or (not isRvzbb() and not isRvzbp()))
    {
      illegalInst(di);
      return;
    }

  unsigned len = 32;
  unsigned mask = len - 1;
  unsigned rot = intRegs_.read(di->op2()) & mask;  // Rotate amount

  uint32_t v1 = intRegs_.read(di->op1());
  uint32_t res32 = (v1 << rot) | (v1 >> ((len - rot) & mask));

  uint64_t res64 = int32_t(res32);  // Sign extend to 64-bits.

  intRegs_.write(di->op0(), res64);
}


template <typename URV>
void
Hart<URV>::execRorw(const DecodedInst* di)
{
  if (not isRv64() or (not isRvzbb() and not isRvzbp()))
    {
      illegalInst(di);
      return;
    }

  unsigned len = 32;
  unsigned mask = len - 1;
  unsigned rot = intRegs_.read(di->op2()) & mask;  // Rotate amount

  uint32_t v1 = intRegs_.read(di->op1());
  uint32_t res32 = (v1 >> rot) | (v1 << ((len - rot) & mask));

  uint64_t res64 = int32_t(res32);  // Sign extend to 64-bits.

  intRegs_.write(di->op0(), res64);
}


template <typename URV>
void
Hart<URV>::execRoriw(const DecodedInst* di)
{
  if (not isRv64() or (not isRvzbb() and not isRvzbp()))
    {
      illegalInst(di);
      return;
    }

  unsigned len = 32;
  unsigned mask = len - 1;
  unsigned rot = di->op2();

  uint32_t v1 = intRegs_.read(di->op1());
  uint32_t res32 = (v1 >> rot) | (v1 << ((len - rot) & mask));

  uint64_t res64 = int32_t(res32);  // Sign extend to 64-bits.

  intRegs_.write(di->op0(), res64);
}


template <typename URV>
void
Hart<URV>::execRev8(const DecodedInst* di)
{
  // Byte swap.

  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());

  if (mxlen_ == 32)
    v1 = __builtin_bswap32(v1);
  else
    v1 = __builtin_bswap64(v1);

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execSext_b(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  int8_t byte = intRegs_.read(di->op1());
  SRV value = byte;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execSext_h(const DecodedInst* di)
{
  if (not isRvzbb())
    {
      illegalInst(di);
      return;
    }

  int16_t half = intRegs_.read(di->op1());
  SRV value = half;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execPack(const DecodedInst* di)
{
  // Zext.h is a zbb pseudo-inst that maps to pack: pack rd, rs1, zero.
  bool zext_h = (di->op2() == 0);

  bool legal = isRvzbe() or isRvzbf() or isRvzbp();
  legal = legal or (isRv64() and isRvzbm());
  if (zext_h)
    legal = legal or isRvzbb();

  if (not legal)
    {
      illegalInst(di);
      return;
    }

  unsigned halfXlen = mxlen_ >> 1;
  URV lower = (intRegs_.read(di->op1()) << halfXlen) >> halfXlen;
  URV upper = intRegs_.read(di->op2()) << halfXlen;
  URV res = upper | lower;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execSlli_uw(const DecodedInst* di)
{
  if (not isRv64() or not isRvzba())
    {
      illegalInst(di);
      return;
    }

  uint32_t amount(di->op2());

  if (amount > 0x3f)
    {
      illegalInst(di);   // Bits 6 of immediate must be zero.
      return;
    }

  URV v1 = uint32_t(intRegs_.read(di->op1()));
  URV value = v1 << amount;

  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execPackh(const DecodedInst* di)
{
  if (not isRvzbe() and not isRvzbf() and not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV lower = intRegs_.read(di->op1()) & 0xff;
  URV upper = (intRegs_.read(di->op2()) & 0xff) << 8;
  URV value = lower | upper;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execPacku(const DecodedInst* di)
{
  bool legal = isRvzbp() or (isRv64() and isRvzbm());
  if (not legal)
    {
      illegalInst(di);
      return;
    }

  unsigned halfXlen = mxlen_ >> 1;

  URV lower = intRegs_.read(di->op1()) >> halfXlen;
  URV upper = (intRegs_.read(di->op2()) >> halfXlen) << halfXlen;
  URV value = lower | upper;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execPackw(const DecodedInst* di)
{
  // zext.h is an alias for packw and is part of zbb.
  if (not isRv64() or (not isRvzbe() and not isRvzbf() and not isRvzbp()))
    {
      illegalInst(di);
      return;
    }

  uint32_t v1 = intRegs_.read(di->op1());
  uint32_t v2 = intRegs_.read(di->op2());

  uint32_t lower = (v1 << 16) >> 16;
  uint32_t upper = v2 << 16;
  uint32_t value = lower | upper;

  int64_t res = int32_t(value);  // Sign extend.
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execPackuw(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  uint32_t v1 = intRegs_.read(di->op1());
  uint32_t v2 = intRegs_.read(di->op2());

  uint32_t lower = v1 >> 16;
  uint32_t upper = (v2 >> 16) << 16;
  uint32_t value = lower | upper;

  int64_t res = int32_t(value);  // sign extend.
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execGrev(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  if (mxlen_ == 32)
    {
      unsigned shamt = v2 & 31;
      if (shamt & 1)
        v1 = ((v1 & 0x55555555) << 1)  | ((v1 & 0xaaaaaaaa) >> 1);
      if (shamt & 2)
        v1 = ((v1 & 0x33333333) << 2)  | ((v1 & 0xcccccccc) >> 2);
      if (shamt & 4)
        v1 = ((v1 & 0x0f0f0f0f) << 4)  | ((v1 & 0xf0f0f0f0) >> 4);
      if (shamt & 8)
        v1 = ((v1 & 0x00ff00ff) << 8)  | ((v1 & 0xff00ff00) >> 8);
      if (shamt & 16)
        v1 = ((v1 & 0x0000ffff) << 16) | ((v1 & 0xffff0000) >> 16);
    }
  else
    {
      int shamt = v2 & 63;
      if (shamt & 1)
        v1 = ((v1 & 0x5555555555555555ll) << 1)  | ((v1 & 0xaaaaaaaaaaaaaaaall) >> 1);
      if (shamt & 2)
        v1 = ((v1 & 0x3333333333333333ll) << 2)  | ((v1 & 0xccccccccccccccccll) >> 2);
      if (shamt & 4)
        v1 = ((v1 & 0x0f0f0f0f0f0f0f0fll) << 4)  | ((v1 & 0xf0f0f0f0f0f0f0f0ll) >> 4);
      if (shamt & 8)
        v1 = ((v1 & 0x00ff00ff00ff00ffll) << 8)  | ((v1 & 0xff00ff00ff00ff00ll) >> 8);
      if (shamt & 16)
        v1 = ((v1 & 0x0000ffff0000ffffll) << 16) | ((v1 & 0xffff0000ffff0000ll) >> 16);
      if (shamt & 32)
        v1 = ((v1 & 0x00000000ffffffffll) << 32) | ((v1 & 0xffffffff00000000ll) >> 32);
    }

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execGrevi(const DecodedInst* di)
{
  URV shamt = di->op2();

  bool zbb = false;  // True if variant is also a zbb instruction.
  if (isRv64())
    zbb = shamt == 0x38;  // rev8 is also in zbb
  else
    zbb = shamt == 0x18;  // rev8 is also in zbb

  bool illegal = not isRvzbp();
  if (zbb)
    illegal = not isRvzbb() and not isRvzbp();

  if (illegal)
    {
      illegalInst(di);
      return;
    }

  if (not checkShiftImmediate(di, shamt))
    return;

  URV v1 = intRegs_.read(di->op1());

  if (mxlen_ == 32)
    {
      if (shamt & 1)
        v1 = ((v1 & 0x55555555) << 1)  | ((v1 & 0xaaaaaaaa) >> 1);
      if (shamt & 2)
        v1 = ((v1 & 0x33333333) << 2)  | ((v1 & 0xcccccccc) >> 2);
      if (shamt & 4)
        v1 = ((v1 & 0x0f0f0f0f) << 4)  | ((v1 & 0xf0f0f0f0) >> 4);
      if (shamt & 8)
        v1 = ((v1 & 0x00ff00ff) << 8)  | ((v1 & 0xff00ff00) >> 8);
      if (shamt & 16)
        v1 = ((v1 & 0x0000ffff) << 16) | ((v1 & 0xffff0000) >> 16);
    }
  else
    {
      if (shamt & 1)
        v1 = ((v1 & 0x5555555555555555ll) << 1)  | ((v1 & 0xaaaaaaaaaaaaaaaall) >> 1);
      if (shamt & 2)
        v1 = ((v1 & 0x3333333333333333ll) << 2)  | ((v1 & 0xccccccccccccccccll) >> 2);
      if (shamt & 4)
        v1 = ((v1 & 0x0f0f0f0f0f0f0f0fll) << 4)  | ((v1 & 0xf0f0f0f0f0f0f0f0ll) >> 4);
      if (shamt & 8)
        v1 = ((v1 & 0x00ff00ff00ff00ffll) << 8)  | ((v1 & 0xff00ff00ff00ff00ll) >> 8);
      if (shamt & 16)
        v1 = ((v1 & 0x0000ffff0000ffffll) << 16) | ((v1 & 0xffff0000ffff0000ll) >> 16);
      if (shamt & 32)
        v1 = ((v1 & 0x00000000ffffffffll) << 32) | ((v1 & 0xffffffff00000000ll) >> 32);
    }

  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execGrevw(const DecodedInst* di)
{
  if (not isRvzbp() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  uint32_t v2 = intRegs_.read(di->op2());

  unsigned shamt = v2 & 31;
  if (shamt & 1)
    v1 = ((v1 & 0x55555555) << 1)  | ((v1 & 0xaaaaaaaa) >> 1);
  if (shamt & 2)
    v1 = ((v1 & 0x33333333) << 2)  | ((v1 & 0xcccccccc) >> 2);
  if (shamt & 4)
    v1 = ((v1 & 0x0f0f0f0f) << 4)  | ((v1 & 0xf0f0f0f0) >> 4);
  if (shamt & 8)
    v1 = ((v1 & 0x00ff00ff) << 8)  | ((v1 & 0xff00ff00) >> 8);
  if (shamt & 16)
    v1 = ((v1 & 0x0000ffff) << 16) | ((v1 & 0xffff0000) >> 16);

  int64_t res = int32_t(v1);  // Sign extend 32-bit to 64-bit.

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execGreviw(const DecodedInst* di)
{
  URV shamt = di->op2();

  bool zbb = false;  // True if variant is also a zbb instruction.
  if (isRv64())
    zbb = shamt == 0x38;  // rev8 is also in zbb
  else
    zbb = shamt == 0x18;  // rev8 is also in zbb

  bool illegal = not isRvzbp();
  if (zbb)
    illegal = not isRvzbb() and not isRvzbp();

  if (illegal)
    {
      illegalInst(di);
      return;
    }

  if (not checkShiftImmediate(di, shamt))
    return;

  unsigned v1 = intRegs_.read(di->op1());

  if (shamt & 1)
    v1 = ((v1 & 0x55555555) << 1)  | ((v1 & 0xaaaaaaaa) >> 1);
  if (shamt & 2)
    v1 = ((v1 & 0x33333333) << 2)  | ((v1 & 0xcccccccc) >> 2);
  if (shamt & 4)
    v1 = ((v1 & 0x0f0f0f0f) << 4)  | ((v1 & 0xf0f0f0f0) >> 4);
  if (shamt & 8)
    v1 = ((v1 & 0x00ff00ff) << 8)  | ((v1 & 0xff00ff00) >> 8);
  if (shamt & 16)
    v1 = ((v1 & 0x0000ffff) << 16) | ((v1 & 0xffff0000) >> 16);

  int64_t res = int32_t(v1);  // Sign extend 32-bit to 64-bit.

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execGorci(const DecodedInst* di)
{
  URV shamt = di->op2();

  bool orc_b = (shamt == 0x7);  // orc.b is also in zbb

  bool legal = isRvzbp();
  if (orc_b)
    legal = legal or isRvzbb();

  if (not legal)
    {
      illegalInst(di);
      return;
    }

  if (not checkShiftImmediate(di, shamt))
    return;

  URV v1 = intRegs_.read(di->op1());

  if (mxlen_ == 32)
    {
      if (shamt & 1)
        v1 |= ((v1 & 0xaaaaaaaa) >>  1) | ((v1 & 0x55555555) <<  1);
      if (shamt & 2)
        v1 |= ((v1 & 0xcccccccc) >>  2) | ((v1 & 0x33333333) <<  2);
      if (shamt & 4)
        v1 |= ((v1 & 0xf0f0f0f0) >>  4) | ((v1 & 0x0f0f0f0f) <<  4);
      if (shamt & 8)
        v1 |= ((v1 & 0xff00ff00) >>  8) | ((v1 & 0x00ff00ff) <<  8);
      if (shamt & 16)
        v1 |= ((v1 & 0xffff0000) >> 16) | ((v1 & 0x0000ffff) << 16);

      intRegs_.write(di->op0(), v1);
      return;
    }

  if constexpr (sizeof(URV) == 8)
    {
      if (shamt & 1)
        v1 |= ((v1 & 0xaaaaaaaaaaaaaaaa) >>  1) | ((v1 & 0x5555555555555555) <<  1);
      if (shamt & 2)
        v1 |= ((v1 & 0xcccccccccccccccc) >>  2) | ((v1 & 0x3333333333333333) <<  2);
      if (shamt & 4)
        v1 |= ((v1 & 0xf0f0f0f0f0f0f0f0) >>  4) | ((v1 & 0x0f0f0f0f0f0f0f0f) <<  4);
      if (shamt & 8)
        v1 |= ((v1 & 0xff00ff00ff00ff00) >>  8) | ((v1 & 0x00ff00ff00ff00ff) <<  8);
      if (shamt & 16)
        v1 |= ((v1 & 0xffff0000ffff0000) >> 16) | ((v1 & 0x0000ffff0000ffff) << 16);
      if (shamt & 32)
        v1 |= ((v1 & 0xffffffff00000000) >> 32) | ((v1 & 0x00000000ffffffff) << 32);

      intRegs_.write(di->op0(), v1);
    }
  else
    assert(0 and "mxlen is 64 in rv32");
}


template <typename URV>
void
Hart<URV>::execGorciw(const DecodedInst* di)
{
  URV shamt = di->op2();

  bool orc_b = (shamt == 0x7);  // orc.b is also in zbb

  bool legal = isRvzbp() and isRv64();
  if (orc_b)
    legal = legal or isRvzbb();

  if (not legal)
    {
      illegalInst(di);
      return;
    }

  if (not checkShiftImmediate(di, shamt))
    return;

  uint32_t v1 = intRegs_.read(di->op1());

  if (shamt & 1)
    v1 |= ((v1 & 0xaaaaaaaa) >>  1) | ((v1 & 0x55555555) <<  1);
  if (shamt & 2)
    v1 |= ((v1 & 0xcccccccc) >>  2) | ((v1 & 0x33333333) <<  2);
  if (shamt & 4)
    v1 |= ((v1 & 0xf0f0f0f0) >>  4) | ((v1 & 0x0f0f0f0f) <<  4);
  if (shamt & 8)
    v1 |= ((v1 & 0xff00ff00) >>  8) | ((v1 & 0x00ff00ff) <<  8);
  if (shamt & 16)
    v1 |= ((v1 & 0xffff0000) >> 16) | ((v1 & 0x0000ffff) << 16);

  int64_t res = int32_t(v1);  // Sign extend 32-bit result to 64-bits.
  intRegs_.write(di->op0(), res);
}


static
uint32_t
shuffleStage32(uint32_t src, uint32_t maskL, uint32_t maskR, unsigned n)
{
  uint32_t x = src & ~(maskL | maskR);
  x |= ((src << n) & maskL) | ((src >> n) & maskR);
  return x;
}


static
uint32_t
shuffle32(uint32_t x, unsigned shamt)
{
  if (shamt & 8)
    x = shuffleStage32(x, 0x00ff0000, 0x0000ff00, 8);
  if (shamt & 4)
    x = shuffleStage32(x, 0x0f000f00, 0x00f000f0, 4);
  if (shamt & 2)
    x = shuffleStage32(x, 0x30303030, 0x0c0c0c0c, 2);
  if (shamt & 1)
    x = shuffleStage32(x, 0x44444444, 0x22222222, 1);

  return x;
}


static
uint32_t
unshuffle32(uint32_t x, unsigned shamt)
{
  if (shamt & 1)
    x = shuffleStage32(x, 0x44444444, 0x22222222, 1);
  if (shamt & 2)
    x = shuffleStage32(x, 0x30303030, 0x0c0c0c0c, 2);
  if (shamt & 4)
    x = shuffleStage32(x, 0x0f000f00, 0x00f000f0, 4);
  if (shamt & 8)
    x = shuffleStage32(x, 0x00ff0000, 0x0000ff00, 8);

  return x;
}


static
uint64_t
shuffleStage64(uint64_t src, uint64_t maskL, uint64_t maskR, unsigned n)
{
  uint64_t x = src & ~(maskL | maskR);
  x |= ((src << n) & maskL) | ((src >> n) & maskR);
  return x;
}


static
uint64_t
shuffle64(uint64_t x, unsigned shamt)
{
  if (shamt & 16)
    x = shuffleStage64(x, 0x0000ffff00000000LL, 0x00000000ffff0000LL, 16);
  if (shamt & 8)
    x = shuffleStage64(x, 0x00ff000000ff0000LL, 0x0000ff000000ff00LL, 8);
  if (shamt & 4)
    x = shuffleStage64(x, 0x0f000f000f000f00LL, 0x00f000f000f000f0LL, 4);
  if (shamt & 2)
    x = shuffleStage64(x, 0x3030303030303030LL, 0x0c0c0c0c0c0c0c0cLL, 2);
  if (shamt & 1)
    x = shuffleStage64(x, 0x4444444444444444LL, 0x2222222222222222LL, 1);

  return x;
}


static
uint64_t
unshuffle64(uint64_t x, unsigned shamt)
{
  if (shamt & 1)
    x = shuffleStage64(x, 0x4444444444444444LL, 0x2222222222222222LL, 1);
  if (shamt & 2)
    x = shuffleStage64(x, 0x3030303030303030LL, 0x0c0c0c0c0c0c0c0cLL, 2);
  if (shamt & 4)
    x = shuffleStage64(x, 0x0f000f000f000f00LL, 0x00f000f000f000f0LL, 4);
  if (shamt & 8)
    x = shuffleStage64(x, 0x00ff000000ff0000LL, 0x0000ff000000ff00LL, 8);
  if (shamt & 16)
    x = shuffleStage64(x, 0x0000ffff00000000LL, 0x00000000ffff0000LL, 16);

  return x;
}


template <typename URV>
void
Hart<URV>::execBset(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV mask = shiftMask();
  unsigned bitIx = intRegs_.read(di->op2()) & mask;

  URV value = intRegs_.read(di->op1()) | (URV(1) << bitIx);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execBclr(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV mask = shiftMask();
  unsigned bitIx = intRegs_.read(di->op2()) & mask;

  URV value = intRegs_.read(di->op1()) & ~(URV(1) << bitIx);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execBinv(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV mask = shiftMask();
  unsigned bitIx = intRegs_.read(di->op2()) & mask;

  URV value = intRegs_.read(di->op1()) ^ (URV(1) << bitIx);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execBext(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV mask = shiftMask();
  unsigned bitIx = intRegs_.read(di->op2()) & mask;

  URV value = (intRegs_.read(di->op1()) >> bitIx) & 1;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execBseti(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV bitIx = di->op2();
  if (not checkShiftImmediate(di, bitIx))
    return;

  URV value = intRegs_.read(di->op1()) | (URV(1) << bitIx);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execBclri(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV bitIx = di->op2();
  if (not checkShiftImmediate(di, bitIx))
    return;

  URV value = intRegs_.read(di->op1()) & ~(URV(1) << bitIx);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execBinvi(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV bitIx = di->op2();
  if (not checkShiftImmediate(di, bitIx))
    return;

  URV value = intRegs_.read(di->op1()) ^ (URV(1) << bitIx);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execBexti(const DecodedInst* di)
{
  if (not isRvzbs())
    {
      illegalInst(di);
      return;
    }

  URV bitIx = di->op2();
  if (not checkShiftImmediate(di, bitIx))
    return;

  URV value = (intRegs_.read(di->op1()) >> bitIx) & 1;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execBcompress(const DecodedInst* di)
{
  if (not isRvzbe())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV res = 0;
  for (unsigned i = 0, j = 0; i < mxlen_; ++i)
    if ((v2 >> i) & 1)
      {
        if ((v1 >> i) & 1)
          res |= URV(1) << j;
        ++j;
      }

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execBdecompress(const DecodedInst* di)
{
  if (not isRvzbe())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV res = 0;
  for (unsigned i = 0, j = 0; i < mxlen_; ++i)
    if ((v2 >> i) & 1)
      {
        if ((v1 >> j) & 1)
          res |= URV(1) << i;
        j++;
      }

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execBcompressw(const DecodedInst* di)
{
  if (not isRvzbe())
    {
      illegalInst(di);
      return;
    }

  uint32_t v1 = intRegs_.read(di->op1());
  uint32_t v2 = intRegs_.read(di->op2());

  uint32_t res = 0;
  for (unsigned i = 0, j = 0; i < 32; ++i)
    if ((v2 >> i) & 1)
      {
        if ((v1 >> i) & 1)
          res |= URV(1) << j;
        ++j;
      }

  int64_t val = int32_t(res);  // sign extend
  intRegs_.write(di->op0(), val);
}


template <typename URV>
void
Hart<URV>::execBdecompressw(const DecodedInst* di)
{
  if (not isRvzbe())
    {
      illegalInst(di);
      return;
    }

  uint32_t v1 = intRegs_.read(di->op1());
  uint32_t v2 = intRegs_.read(di->op2());

  uint32_t res = 0;
  for (unsigned i = 0, j = 0; i < 32; ++i)
    if ((v2 >> i) & 1)
      {
        if ((v1 >> j) & 1)
          res |= URV(1) << i;
        j++;
      }

  int64_t val = int32_t(res);  // sign extend.
  intRegs_.write(di->op0(), val);
}


template <typename URV>
void
Hart<URV>::execBfp(const DecodedInst* di)
{
  if (not isRvzbf())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV cfg = v2 >> (mxlen_ / 2);
  if ((cfg >> 30) == 2)
    cfg = cfg >> 16;

  unsigned len = (cfg >> 8) & (mxlen_ / 2 - 1);
  unsigned off = cfg & (mxlen_ - 1);
  len = len ? len : mxlen_ / 2;
  URV mask = ~(~URV(0) << len);
  mask = mask << off;
  URV data = v2 << off;

  URV res = (data & mask) | (v1 & ~mask);
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execBfpw(const DecodedInst* di)
{
  if (not isRvzbf() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint32_t v1 = intRegs_.read(di->op1());
  uint32_t v2 = intRegs_.read(di->op2());

  uint32_t cfg = v2 >> 16;
  if ((cfg >> 30) == 2)
    cfg = cfg >> 16;

  unsigned xlen = 32;
  uint32_t len = (cfg >> 8) & (xlen / 2 - 1);
  uint32_t off = cfg & (xlen - 1);
  len = len ? len : xlen / 2;
  uint32_t mask = ~(~uint32_t(0) << len);
  mask = mask << off;
  uint32_t data = v2 << off;

  uint32_t value = (data & mask) | (v1 & ~mask);
  int64_t res = int32_t(value);  // Sign extend.
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execGorc(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  uint32_t mask = mxlen_ - 1;  // 0x1f in rv32 and 0x3f in rv64
  uint32_t shamt = intRegs_.read(di->op2()) & mask;

  if (mxlen_ == 32)
    {
      if (shamt & 1)
        v1 |= ((v1 & 0xaaaaaaaa) >>  1) | ((v1 & 0x55555555) <<  1);
      if (shamt & 2)
        v1 |= ((v1 & 0xcccccccc) >>  2) | ((v1 & 0x33333333) <<  2);
      if (shamt & 4)
        v1 |= ((v1 & 0xf0f0f0f0) >>  4) | ((v1 & 0x0f0f0f0f) <<  4);
      if (shamt & 8)
      v1 |= ((v1 & 0xff00ff00) >>  8) | ((v1 & 0x00ff00ff) <<  8);
      if (shamt & 16)
        v1 |= ((v1 & 0xffff0000) >> 16) | ((v1 & 0x0000ffff) << 16);

      intRegs_.write(di->op0(), v1);
      return;
    }

  if constexpr (sizeof(URV) == 8)
    {
      if (shamt & 1)
        v1 |= ((v1 & 0xaaaaaaaaaaaaaaaa) >>  1) | ((v1 & 0x5555555555555555) <<  1);
      if (shamt & 2)
        v1 |= ((v1 & 0xcccccccccccccccc) >>  2) | ((v1 & 0x3333333333333333) <<  2);
      if (shamt & 4)
        v1 |= ((v1 & 0xf0f0f0f0f0f0f0f0) >>  4) | ((v1 & 0x0f0f0f0f0f0f0f0f) <<  4);
      if (shamt & 8)
        v1 |= ((v1 & 0xff00ff00ff00ff00) >>  8) | ((v1 & 0x00ff00ff00ff00ff) <<  8);
      if (shamt & 16)
        v1 |= ((v1 & 0xffff0000ffff0000) >> 16) | ((v1 & 0x0000ffff0000ffff) << 16);
      if (shamt & 32)
        v1 |= ((v1 & 0xffffffff00000000) >> 32) | ((v1 & 0x00000000ffffffff) << 32);

      intRegs_.write(di->op0(), v1);
    }
  else
    assert(0 and "mxlen is 64 in rv32");
}


template <typename URV>
void
Hart<URV>::execGorcw(const DecodedInst* di)
{
  if (not isRvzbp() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  URV v1 = uint32_t(intRegs_.read(di->op1()));  // Clear most sig 32 bits
  uint32_t shamt = intRegs_.read(di->op2()) & 0x1f;

  if (shamt & 1)
    v1 |= ((v1 & 0xaaaaaaaa) >>  1) | ((v1 & 0x55555555) <<  1);
  if (shamt & 2)
    v1 |= ((v1 & 0xcccccccc) >>  2) | ((v1 & 0x33333333) <<  2);
  if (shamt & 4)
    v1 |= ((v1 & 0xf0f0f0f0) >>  4) | ((v1 & 0x0f0f0f0f) <<  4);
  if (shamt & 8)
    v1 |= ((v1 & 0xff00ff00) >>  8) | ((v1 & 0x00ff00ff) <<  8);
  if (shamt & 16)
    v1 |= ((v1 & 0xffff0000) >> 16) | ((v1 & 0x0000ffff) << 16);

  v1 = SRV(int32_t(v1));  // Sign extend least sig 32-bits to 64-bits.
  intRegs_.write(di->op0(), v1);
}


template <typename URV>
void
Hart<URV>::execClmul(const DecodedInst* di)
{
  if (not isRvzbc())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV x = 0;
  for (unsigned i = 0; i < mxlen_; ++i)
    if ((v2 >> i) & 1)
      x ^= v1 << i;

  intRegs_.write(di->op0(), x);
}


template <typename URV>
void
Hart<URV>::execClmulh(const DecodedInst* di)
{
  if (not isRvzbc())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV x = 0;
  for (unsigned i = 1; i < mxlen_; ++i)
    if ((v2 >> i) & 1)
      x ^= v1 >> (mxlen_ - i);

  intRegs_.write(di->op0(), x);
}


template <typename URV>
void
Hart<URV>::execShfl(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV val = 0;

  if (mxlen_ == 32)
    {
      unsigned shamt = v2 & 15;
      val = shuffle32(v1, shamt);
    }
  else
    {
      unsigned shamt = v2 & 31;
      val = shuffle64(v1, shamt);
    }

  intRegs_.write(di->op0(), val);
}


template <typename URV>
void
Hart<URV>::execShflw(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  uint32_t v1 = intRegs_.read(di->op1());
  unsigned shamt = intRegs_.read(di->op2()) & 0xf;

  uint32_t res32 = shuffle32(v1, shamt);

  int64_t res = int32_t(res32);   // Sign extend to 64-bits.
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execShfli(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV amt = di->op2();
  URV val = 0;

  if (mxlen_ == 32)
    {
      if (amt > 15)
        {
          illegalInst(di);
          return;
        }
      val = shuffle32(v1, amt);
    }
  else
    {
      if (amt > 31)
        {
          illegalInst(di);
          return;
        }
      val = shuffle64(v1, amt);
    }

  intRegs_.write(di->op0(), val);
}


template <typename URV>
void
Hart<URV>::execUnshfl(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV val = 0;

  if (mxlen_ == 32)
    {
      unsigned shamt = v2 & 15;
      val = unshuffle32(v1, shamt);
    }
  else
    {
      unsigned shamt = v2 & 31;
      val = unshuffle64(v1, shamt);
    }

  intRegs_.write(di->op0(), val);
}


template <typename URV>
void
Hart<URV>::execUnshfli(const DecodedInst* di)
{
  URV amt = di->op2();

  bool legal = isRvzbp();

  // Instructions zip8 and unzip16 are aliases to unshfli for
  // immediate values of 0x18 and 0x10 and are legal when zbm is on.
  legal = legal or (isRvzbm() and (amt == 0x18 or amt == 0x10));

  if (not legal)
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV val = 0;

  if (mxlen_ == 32)
    {
      if (amt > 15)
        {
          illegalInst(di);
          return;
        }
      val = unshuffle32(v1, amt);
    }
  else
    {
      if (amt > 31)
        {
          illegalInst(di);
          return;
        }
      val = unshuffle64(v1, amt);
    }

  intRegs_.write(di->op0(), val);
}


template <typename URV>
void
Hart<URV>::execUnshflw(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  uint32_t v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  unsigned shamt = v2 & 15;
  uint32_t val = unshuffle32(v1, shamt);

  uint64_t res = int32_t(val);  // Sign extend to 64-bits.

  intRegs_.write(di->op0(), res);
}


static
uint32_t
xperm32(uint32_t v1, uint32_t v2, unsigned log2Width)
{
  uint32_t res = 0;
  uint32_t width = 1u << log2Width;
  uint32_t mask = (1u << width) - 1;
  for (unsigned i = 0; i < 32; i += width)
    {
      uint32_t pos = ((v2 >> i) & mask) << log2Width;
      if (pos < 32)
        res |= ((v1 >> pos) & mask) << i;
    }
  return res;
}

  
static
uint64_t
xperm64(uint64_t v1, uint64_t v2, unsigned log2Width)
{
  uint64_t res = 0;
  uint32_t width = 1u << log2Width;
  uint64_t mask = (uint64_t(1) << width) - 1;
  for (unsigned i = 0; i < 64; i += width)
    {
      uint64_t pos = ((v2 >> i) & mask) << log2Width;
      if (pos < 64)
        res |= ((v1 >> pos) & mask) << i;
    }
  return res;
}


template <typename URV>
void
Hart<URV>::execXperm_n(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV res = 0;

  if (mxlen_ == 32)
    res = xperm32(v1, v2, 2);
  else
    res = xperm64(v1, v2, 2);

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execXperm_b(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV res = 0;

  if (mxlen_ == 32)
    res = xperm32(v1, v2, 3);
  else
    res = xperm64(v1, v2, 3);

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execXperm_h(const DecodedInst* di)
{
  if (not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV res = 0;

  if (mxlen_ == 32)
    res = xperm32(v1, v2, 4);
  else
    res = xperm64(v1, v2, 4);

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execXperm_w(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbp())
    {
      illegalInst(di);
      return;
    }

  uint64_t v1 = intRegs_.read(di->op1());
  uint64_t v2 = intRegs_.read(di->op2());
  uint64_t res = xperm64(v1, v2, 5);

  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execClmulr(const DecodedInst* di)
{
  if (not isRvzbc())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV x = 0;
  for (unsigned i = 0; i < mxlen_; ++i)
    if ((v2 >> i) & 1)
      x ^= v1 >> (mxlen_ - i - 1);

  intRegs_.write(di->op0(), x);
}


template <typename URV>
void
Hart<URV>::execSh1add(const DecodedInst* di)
{
  if (not isRvzba())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV res = (v1 << 1) + v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execSh2add(const DecodedInst* di)
{
  if (not isRvzba())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV res = (v1 << 2) + v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execSh3add(const DecodedInst* di)
{
  if (not isRvzba())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());

  URV res = (v1 << 3) + v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execSh1add_uw(const DecodedInst* di)
{
  if (not isRv64() or not isRvzba())
    {
      illegalInst(di);
      return;
    }

  URV v1 = uint32_t(intRegs_.read(di->op1()));
  URV v2 = intRegs_.read(di->op2());

  URV res = (v1 << 1) + v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execSh2add_uw(const DecodedInst* di)
{
  if (not isRv64() or not isRvzba())
    {
      illegalInst(di);
      return;
    }

  URV v1 = uint32_t(intRegs_.read(di->op1()));
  URV v2 = intRegs_.read(di->op2());

  URV res = (v1 << 2) + v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execSh3add_uw(const DecodedInst* di)
{
  if (not isRv64() or not isRvzba())
    {
      illegalInst(di);
      return;
    }

  URV v1 = uint32_t(intRegs_.read(di->op1()));
  URV v2 = intRegs_.read(di->op2());

  URV res = (v1 << 3) + v2;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execAdd_uw(const DecodedInst* di)
{
  if (not isRv64() or not isRvzba())
    {
      illegalInst(di);
      return;
    }

  URV v1 = uint32_t(intRegs_.read(di->op1()));
  URV v2 = intRegs_.read(di->op2());

  URV value = v1 + v2;
  intRegs_.write(di->op0(), value);
}


template <typename URV>
static
URV
crc32(URV x, unsigned nbits)
{
  for (unsigned i = 0; i < nbits; ++i)
    x = (x >> 1) ^ (0xedb88320 & ~((x & 1) - 1));
  return x;
}


template <typename URV>
static
URV
crc32c(URV x, unsigned nbits)
{
  for (unsigned i = 0; i < nbits; ++i)
    x = (x >> 1) ^ (0x82F63B78  & ~((x & 1) - 1));
  return x;
}



template <typename URV>
void
Hart<URV>::execCrc32_b(const DecodedInst* di)
{
  if (not isRvzbr())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32(intRegs_.read(di->op1()), 8);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execCrc32_h(const DecodedInst* di)
{
  if (not isRvzbr())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32(intRegs_.read(di->op1()), 16);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execCrc32_w(const DecodedInst* di)
{
  if (not isRvzbr())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32(intRegs_.read(di->op1()), 32);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execCrc32_d(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbr())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32(intRegs_.read(di->op1()), 64);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execCrc32c_b(const DecodedInst* di)
{
  if (not isRvzbr())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32c(intRegs_.read(di->op1()), 8);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execCrc32c_h(const DecodedInst* di)
{
  if (not isRvzbr())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32c(intRegs_.read(di->op1()), 16);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execCrc32c_w(const DecodedInst* di)
{
  if (not isRvzbr())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32c(intRegs_.read(di->op1()), 32);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execCrc32c_d(const DecodedInst* di)
{
  if (not isRv64() or not isRvzbr())
    {
      illegalInst(di);
      return;
    }
  URV value = crc32c(intRegs_.read(di->op1()), 64);
  intRegs_.write(di->op0(), value);
}


template <typename URV>
void
Hart<URV>::execBmator(const DecodedInst* di)
{
  if (not isRvzbm() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint8_t u[8]; // rows of rs1
  uint8_t v[8]; // cols of rs2

  uint64_t rs1 = intRegs_.read(di->op1());
  uint64_t rs2 = intRegs_.read(di->op2());

  uint64_t rs2t = rs2;
  rs2t = shuffle64(rs2t, 31);
  rs2t = shuffle64(rs2t, 31);
  rs2t = shuffle64(rs2t, 31);

  for (int i = 0; i < 8; i++)
    {
      u[i] = rs1 >> (i*8);
      v[i] = rs2t >> (i*8);
    }

  uint64_t x = 0;
  for (int i = 0; i < 64; i++)
    {
      if ((u[i / 8] & v[i % 8]) != 0)
        x |= 1LL << i;
    }

  intRegs_.write(di->op0(), x);
}


template <typename URV>
void
Hart<URV>::execBmatxor(const DecodedInst* di)
{
  if (not isRvzbm() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint8_t u[8]; // rows of rs1
  uint8_t v[8]; // cols of rs2

  uint64_t rs1 = intRegs_.read(di->op1());
  uint64_t rs2 = intRegs_.read(di->op2());

  uint64_t rs2t = rs2;
  rs2t = shuffle64(rs2t, 31);
  rs2t = shuffle64(rs2t, 31);
  rs2t = shuffle64(rs2t, 31);

  for (int i = 0; i < 8; i++)
    {
      u[i] = rs1 >> (i*8);
      v[i] = rs2t >> (i*8);
    }

  uint64_t x = 0;
  for (int i = 0; i < 64; i++)
    {
      if (__builtin_popcount(u[i / 8] & v[i % 8]) & 1)
        x ^= 1LL << i;
    }

  intRegs_.write(di->op0(), x);
}


template <typename URV>
void
Hart<URV>::execBmatflip(const DecodedInst* di)
{
  if (not isRvzbm() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint64_t rs1 = intRegs_.read(di->op1());
  rs1 = shuffle64(rs1, 31);
  rs1 = shuffle64(rs1, 31);
  rs1 = shuffle64(rs1, 31);

  intRegs_.write(di->op0(), rs1);
}


template <typename URV>
void
Hart<URV>::execCmov(const DecodedInst* di)
{
  if (not isRvzbt())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV v3 = intRegs_.read(di->op3());

  URV res = v2 ? v1 : v3;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execCmix(const DecodedInst* di)
{
  if (not isRvzbt())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV v3 = intRegs_.read(di->op3());

  URV res = (v1 & v2) | (v3 & ~v2);
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFsl(const DecodedInst* di)
{
  if (not isRvzbt())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV v3 = intRegs_.read(di->op3());

  unsigned shamt = v2 & (2*mxlen_ - 1);

  URV aa = v1, bb = v3;

  if (shamt >= mxlen_)
    {
      shamt -= mxlen_;
      aa = v3;
      bb = v1;
    }

  URV res = shamt ? (aa << shamt) | (bb >> (mxlen_ - shamt)) : aa;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFsr(const DecodedInst* di)
{
  if (not isRvzbt())
    {
      illegalInst(di);
      return;
    }

  URV v1 = intRegs_.read(di->op1());
  URV v2 = intRegs_.read(di->op2());
  URV v3 = intRegs_.read(di->op3());

  unsigned shamt = v2 & (2*mxlen_ - 1);

  URV aa = v1, bb = v3;

  if (shamt >= mxlen_)
    {
      shamt -= mxlen_;
      aa = v3;
      bb = v1;
    }

  URV res = shamt ? (aa >> shamt) | (bb << (mxlen_ - shamt)) : aa;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFsri(const DecodedInst* di)
{
  if (not isRvzbt())
    {
      illegalInst(di);
      return;
    }

  URV aa = intRegs_.read(di->op1());
  URV bb = intRegs_.read(di->op2());
  URV imm = intRegs_.read(di->op3());

  unsigned shamt = imm & (2*mxlen_ - 1);

  if (shamt >= mxlen_)
    {
      shamt -= mxlen_;
      std::swap(aa, bb);
    }

  URV res = shamt ? (aa >> shamt) | (bb << (mxlen_ - shamt)) : aa;
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFslw(const DecodedInst* di)
{
  if (not isRvzbt() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint32_t v1 = intRegs_.read(di->op1());
  uint32_t v2 = intRegs_.read(di->op2());
  uint32_t v3 = intRegs_.read(di->op3());

  uint32_t len = 32;

  unsigned shamt = v2 & (2*len - 1);

  URV aa = v1, bb = v3;

  if (shamt >= len)
    {
      shamt -= len;
      aa = v3;
      bb = v1;
    }

  uint32_t val = shamt ? (aa << shamt) | (bb >> (len - shamt)) : aa;

  int64_t res = int32_t(val);  // sign extend
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFsrw(const DecodedInst* di)
{
  if (not isRvzbt() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint32_t v1 = intRegs_.read(di->op1());
  uint32_t v2 = intRegs_.read(di->op2());
  uint32_t v3 = intRegs_.read(di->op3());

  uint32_t len = 32;

  unsigned shamt = v2 & (2*len - 1);

  URV aa = v1, bb = v3;

  if (shamt >= len)
    {
      shamt -= len;
      aa = v3;
      bb = v1;
    }

  uint32_t val = shamt ? (aa >> shamt) | (bb << (len - shamt)) : aa;

  int64_t res = int32_t(val);  // sign extend
  intRegs_.write(di->op0(), res);
}


template <typename URV>
void
Hart<URV>::execFsriw(const DecodedInst* di)
{
  if (not isRvzbt() or not isRv64())
    {
      illegalInst(di);
      return;
    }

  uint32_t aa = intRegs_.read(di->op1());
  uint32_t bb = intRegs_.read(di->op2());
  uint32_t imm = intRegs_.read(di->op3());

  unsigned len = 32;

  unsigned shamt = imm & (2*len - 1);

  if (shamt >= len)
    {
      shamt -= len;
      std::swap(aa, bb);
    }

  uint32_t val = shamt ? (aa >> shamt) | (bb << (len - shamt)) : aa;

  int64_t res = int32_t(val);  // sign extend
  intRegs_.write(di->op0(), res);
}


template class WdRiscv::Hart<uint32_t>;
template class WdRiscv::Hart<uint64_t>;
