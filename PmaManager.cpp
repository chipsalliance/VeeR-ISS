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
#include <cmath>
#include <cassert>
#include "PmaManager.hpp"

using namespace WdRiscv;


PmaManager::PmaManager(uint64_t memSize, uint64_t sectionSize)
  : memSize_(memSize), sectionSize_(sectionSize)
{
  assert(memSize >= sectionSize);
  assert(sectionSize >= 64);

  uint64_t logSectionSize = static_cast<uint64_t>(std::log2(sectionSize_));
  uint64_t p2SectionSize = uint64_t(1) << logSectionSize;
  assert(p2SectionSize == sectionSize_);
  sectionShift_ = logSectionSize;

  uint64_t sectionCount = memSize_ / sectionSize_;
  assert(sectionCount * sectionSize_ == memSize_);

  // Whole memory is intially set for instruction/data/atomic access.
  // No iccm/dccm/mmr/io.
  sectionPmas_.resize(sectionCount, Pma::Default);
}


void
PmaManager::enable(uint64_t a0, uint64_t a1, Pma::Attrib attrib)
{
  a0 = (a0 >> 2) << 2;   // Make word aligned.
  a1 = (a1 >> 2) << 2;   // Make word aligned.

  unsigned mask = attrib;

  while (a0 <= a1)
    {
      uint64_t p0 = getSectionStartAddr(a0);
      uint64_t p1 = getSectionStartAddr(a1);
      bool doWord = (a0 != p0) or (p0 == p1 and (a1 - a0 + 4 != sectionSize_));
      Pma pma = getPma(a0);
      doWord = doWord or pma.word_;
      if (doWord)
        {
          fracture(a0);
          uint64_t last = std::min(a1, p0 + sectionSize_ - 4);
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
          uint64_t sectionIx = getSectionIx(p0);
          sectionPmas_[sectionIx] = pma;
          a0 += sectionSize_;
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
      uint64_t p0 = getSectionStartAddr(a0);
      uint64_t p1 = getSectionStartAddr(a1);
      bool doWord = (a0 != p0) or (p0 == p1 and (a1 - a0 + 4 != sectionSize_));
      Pma pma = getPma(a0);
      doWord = doWord or pma.word_;
      if (doWord)
        {
          fracture(a0);
          uint64_t last = std::min(a1, p0 + sectionSize_ - 4);
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
          uint64_t sectionIx = getSectionIx(p0);
          sectionPmas_[sectionIx] = pma;
          a0 += sectionSize_;
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
      uint64_t p0 = getSectionStartAddr(a0);
      uint64_t p1 = getSectionStartAddr(a1);
      bool doWord = (a0 != p0) or (p0 == p1 and (a1 - a0 + 4 != sectionSize_));
      Pma pma = getPma(a0);
      doWord = doWord or pma.word_;
      if (doWord)
        {
          fracture(a0);
          Pma pma(attrib);
          pma.word_ = true;
          uint64_t last = std::min(a1, p0 + sectionSize_ - 4);
          for ( ; a0 <= last; a0 += 4)
            wordPmas_[a0>>2] = pma;
        }
      else
        {
          uint64_t sectionIx = getSectionIx(p0);
          sectionPmas_[sectionIx] = Pma(attrib);
          a0 += sectionSize_;
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

  // Mark old area as non-memory-mapped.
  disable(memMappedBase_, memMappedBase_ + memMappedSize_ - 1, Pma::MemMapped);

  // Mark new area as memory-mapped.
  enable(newBase, newBase + memMappedSize_ - 1, Pma::MemMapped);

  memMappedBase_ = newBase;
  return true;
}


void
PmaManager::fracture(uint64_t addr)
{
  uint64_t sectionIx = getSectionIx(addr);
  Pma pma = sectionPmas_.at(sectionIx);
  if (pma.word_)
    return;
  pma.word_= true;
  sectionPmas_.at(sectionIx) = pma;

  uint64_t words = sectionSize_ / 4;
  uint64_t wordIx = (sectionIx*sectionSize_) >> 2;
  for (uint64_t i = 0; i < words; ++i, wordIx++)
    wordPmas_[wordIx] = pma;
}
