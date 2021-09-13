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


bool
PmaManager::changeMemMappedBase(uint64_t newBase)
{
	if(memMappedSections_.size() != 1)
		return false;
	auto mms = memMappedSections_[0];
	uint64_t currBase = mms.base;
	if (newBase == currBase)
		return true;
	if (memSize_ - mms.size < newBase)
	    return false;
	std::unordered_map<uint64_t, MemMappedRegister> tmp;
	tmp.insert(memMappedRegs_.begin(), memMappedRegs_.end());
	memMappedRegs_.clear();
	for(auto& r: tmp)
		memMappedRegs_.insert(std::make_pair(uint64_t(newBase+currBase-r.first), r.second));

	// Mark old area as non-memory-mapped.
	disable(mms.base, mms.base+mms.size - 1, Pma::MemMapped);

	// Mark new area as memory-mapped.
	enable(newBase, newBase + mms.size - 1, Pma::MemMapped);

	mms.base = newBase;
	return true;
}


void
PmaManager::setMemMappedMask(uint64_t addr, uint32_t mask, uint8_t size)
{
	assert((addr & (size-1)) == 0);
	assert(getMemMappedSection(addr)>=0);
	memMappedRegs_.insert(std::make_pair(addr, MemMappedRegister(mask, size)));

}


uint64_t
PmaManager::getMemMappedMask(uint64_t addr, uint64_t& size) const
{
	auto t = memMappedRegs_.find(addr);
	if(t != memMappedRegs_.end()) {
		size = t->second.size;
		return t->second.mask;
	}
	size =  4;
	return 0;
}


bool
PmaManager::defineMemMappedArea(uint64_t base, uint64_t size, bool isInternal)
{
	if(size & 3) return false;
	memMappedSections_.push_back(MemMappedSection(base, size, isInternal));
	return true;
}



bool
PmaManager::writeRegisterByte(uint64_t addr, uint8_t value)
{
    uint64_t wordAddr;
    if(auto r = getRegister(addr, wordAddr)) {
		unsigned shift = (addr & (r->size-1))*8;
		uint64_t mask =  uint64_t(0xff) << shift;
		r->data = ((int64_t(value)<<shift) | (r->data & ~mask)) & r->mask;
		return true;
    }
	return getMemMappedSection(addr) >= 0;
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
