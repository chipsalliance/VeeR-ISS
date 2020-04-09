//
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2018-2019-2020 Western Digital Corporation or its affiliates.
// 
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
// 
// You should have received a copy of the GNU General Public License along with
// this program. If not, see <https://www.gnu.org/licenses/>.
//

#include <iostream>
#include <cmath>
#include <cassert>
#include "PmaManager.hpp"

using namespace WdRiscv;


PmaManager::PmaManager(uint64_t memSize, uint64_t pageSize)
  : memSize_(memSize), pageSize_(pageSize)
{
  assert(memSize >= pageSize);
  assert(pageSize >= 64);

  uint64_t logPageSize = static_cast<uint64_t>(std::log2(pageSize_));
  uint64_t p2PageSize = uint64_t(1) << logPageSize;
  assert(p2PageSize == pageSize_);
  pageShift_ = logPageSize;

  uint64_t pageCount = memSize_ / pageSize_;
  assert(pageCount * pageSize_ == memSize_);

  // Whole memory is intially set for instruction/data/atomic access.
  // No iccm/dccm/mmr/io.
  pagePmas_.resize(pageCount, Pma::Default);
}


void
PmaManager::enable(uint64_t a0, uint64_t a1, Pma::Attrib attrib)
{
  a0 = (a0 >> 2) << 2;   // Make word aligned.
  a1 = (a1 >> 2) << 2;   // Make word aligned.

  unsigned mask = attrib;

  while (a0 <= a1)
    {
      uint64_t p0 = getPageStartAddr(a0);
      uint64_t p1 = getPageStartAddr(a1);
      bool doWord = (a0 != p0) or (p0 == p1 and (a1 - a0 + 4 != pageSize_));
      Pma pma = getPma(a0);
      doWord = doWord or pma.word_;
      if (doWord)
        {
          fracture(a0);
          uint64_t last = std::min(a1, p0 + pageSize_ - 4);
          for ( ; a0 <= last; a0 += 4)
            {
              pma = getPma(a0);
              pma.attrib_ = pma.attrib_ | mask;
              wordPmas_[a0>>2] = pma;
            }
        }
      else
        {
          pma.attrib_ = pma.attrib_ | mask;
          uint64_t pageIx = getPageIx(p0);
          pagePmas_[pageIx] = pma;
          a0 += pageSize_;
        }
    }
}


void
PmaManager::disable(uint64_t a0, uint64_t a1, Pma::Attrib attrib)
{
  a0 = (a0 >> 2) << 2;   // Make word aligned.
  a1 = (a1 >> 2) << 2;   // Make word aligned.

  unsigned mask = ~attrib;

  while (a0 <= a1)
    {
      uint64_t p0 = getPageStartAddr(a0);
      uint64_t p1 = getPageStartAddr(a1);
      bool doWord = (a0 != p0) or (p0 == p1 and (a1 - a0 + 4 != pageSize_));
      Pma pma = getPma(a0);
      doWord = doWord or pma.word_;
      if (doWord)
        {
          fracture(a0);
          uint64_t last = std::min(a1, p0 + pageSize_ - 4);
          for ( ; a0 <= last; a0 += 4)
            {
              pma = getPma(a0);
              pma.attrib_ = pma.attrib_ & mask;
              wordPmas_[a0>>2] = pma;
            }
        }
      else
        {
          pma.attrib_ = pma.attrib_ & mask;
          uint64_t pageIx = getPageIx(p0);
          pagePmas_[pageIx] = pma;
          a0 += pageSize_;
        }
    }
}


void
PmaManager::setAttribute(uint64_t a0, uint64_t a1, Pma::Attrib attrib)
{
  a0 = (a0 >> 2) << 2;   // Make word aligned.
  a1 = (a1 >> 2) << 2;   // Make word aligned.

  while (a0 <= a1)
    {
      uint64_t p0 = getPageStartAddr(a0);
      uint64_t p1 = getPageStartAddr(a1);
      bool doWord = (a0 != p0) or (p0 == p1 and (a1 - a0 + 4 != pageSize_));
      Pma pma = getPma(a0);
      doWord = doWord or pma.word_;
      if (doWord)
        {
          fracture(a0);
          Pma pma(attrib);
          pma.word_ = true;
          uint64_t last = std::min(a1, p0 + pageSize_ - 4);
          for ( ; a0 <= last; a0 += 4)
            wordPmas_[a0>>2] = pma;
        }
      else
        {
          uint64_t pageIx = getPageIx(p0);
          pagePmas_[pageIx] = Pma(attrib);
          a0 += pageSize_;
        }
    }
}


void
PmaManager::setMemMappedMask(uint64_t addr, uint32_t mask)
{
  assert(addr >= memMappedBase_);
  uint64_t wordIx = (addr - memMappedBase_) / 4;
  assert(wordIx < memMappedMasks_.size());
  memMappedMasks_.at(wordIx) = mask;
}


uint32_t
PmaManager::getMemMappedMask(uint64_t addr) const
{
  if (addr < memMappedBase_)
    return 0xffffffff;
  uint64_t wordIx = (addr - memMappedBase_) / 4;
  if (wordIx < memMappedMasks_.size())
    return memMappedMasks_[wordIx];
  return 0xffffffff;
}


bool
PmaManager::defineMemMappedArea(uint64_t base, uint64_t size)
{
  memMappedBase_ = base;
  memMappedSize_ = (size >> 2) << 2;
  uint64_t wordCount = memMappedSize_ / 4;
  memMappedMasks_.resize(wordCount);
  memMappedRegs_.resize(wordCount);
  return size == memMappedSize_;
}


bool
PmaManager::readRegister(uint64_t addr, uint32_t& value) const
{
  if ((addr & 3) != 0)
    return false;  // Address must be workd-aligned.
  if (addr < memMappedBase_)
    return false;
  uint64_t wordIx = (addr - memMappedBase_) / 4;
  if (wordIx >= memMappedRegs_.size())
    return false;
  value = memMappedRegs_[wordIx];
  return true;
}


bool
PmaManager::writeRegister(uint64_t addr, uint32_t value)
{
  if ((addr & 3) != 0)
    return false;  // Address must be workd-aligned.
  if (addr < memMappedBase_)
    return false;
  uint64_t wordIx = (addr - memMappedBase_) / 4;
  if (wordIx >= memMappedRegs_.size())
    return false;
  uint32_t mask = memMappedMasks_.at(wordIx);
  memMappedRegs_.at(wordIx) = value & mask;
  return true;
}


bool
PmaManager::writeRegisterNoMask(uint64_t addr, uint32_t value)
{
  if ((addr & 3) != 0)
    return false;  // Address must be workd-aligned.
  if (addr < memMappedBase_)
    return false;
  uint64_t wordIx = (addr - memMappedBase_) / 4;
  if (wordIx >= memMappedRegs_.size())
    return false;
  memMappedRegs_.at(wordIx) = value;
  return true;
}


bool
PmaManager::writeRegisterByte(uint64_t addr, uint8_t value)
{
  if (addr < memMappedBase_)
    return false;
  uint64_t wordIx = (addr - memMappedBase_) / 4;
  if (wordIx >= memMappedRegs_.size())
    return false;
  unsigned byteIx = addr & 3;
  unsigned shift = byteIx * 8;
  uint32_t byteMask = 0xff << shift;
  uint32_t vv = (uint32_t(value) << shift) & memMappedMasks_.at(wordIx);
  memMappedRegs_.at(wordIx) = memMappedRegs_.at(wordIx) & ~byteMask;
  memMappedRegs_.at(wordIx) = memMappedRegs_.at(wordIx) | vv;
  return true;
}


bool
PmaManager::changeMemMappedBase(uint64_t newBase)
{
  if (newBase == memMappedBase_)
    return true;

  if (memSize_ - memMappedSize_ < newBase)
    return false;

  // Makr old area as non-memory-mapped.
  disable(memMappedBase_, memMappedBase_ + memMappedSize_ - 1, Pma::MemMapped);

  // Mark new area as read/write/memory-mapped.
  Pma::Attrib attrib = Pma::Attrib(Pma::Read | Pma::Write | Pma::MemMapped);
  setAttribute(newBase, newBase + memMappedSize_ - 1, attrib);

  memMappedBase_ = newBase;
  return true;
}


void
PmaManager::fracture(uint64_t addr)
{
  uint64_t pageIx = getPageIx(addr);
  Pma pma = pagePmas_.at(pageIx);
  if (pma.word_)
    return;
  pma.word_= true;

  uint64_t words = pageSize_ / 4;
  uint64_t wordIx = (pageIx*pageSize_) >> 2;
  for (uint64_t i = 0; i < words; ++i, wordIx++)
    wordPmas_[wordIx] = pma;
}
