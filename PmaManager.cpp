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
#include "PmaManager.hpp"

using namespace WdRiscv;


PmaManager::PmaManager(uint64_t memSize, uint64_t pageSize)
  : memSize_(memSize), pageSize_(pageSize)
{
  if (pageSize_ >= 64*1024*1024)
    {
      pageSize_ = 64*1024*1024;
      std::cerr << "Page size (" << pageSize << ") too large: Using "
                << pageSize_ << '\n';
    }

  if ((memSize & 4) != 0)
    {
      memSize_ = (memSize_ >> 2) << 2;
      std::cerr << "Memory size (" << memSize << ") is not a multiple of 4."
                << " Using " << memSize_ << '\n';
    }

  uint64_t logPageSize = static_cast<uint64_t>(std::log2(pageSize_));
  uint64_t p2PageSize = uint64_t(1) << logPageSize;
  if (p2PageSize != pageSize_)
    {
      std::cerr << "Memory page size (0x" << std::hex << pageSize_ << ") "
		<< "is not a power of 2 -- using 0x" << p2PageSize << '\n'
		<< std::dec;
      pageSize_ = p2PageSize;
    }
  pageShift_ = logPageSize;

  if (memSize_ < pageSize_)
    {
      std::cerr << "Unreasonably small memory size (less than 0x "
		<< std::hex << pageSize_ << ") -- using 0x" << pageSize_
		<< '\n' << std::dec;
      memSize_ = pageSize_;
    }

  uint64_t pageCount = memSize_ / pageSize_;
  if (pageCount * pageSize_ != memSize_)
    {
      pageCount++;
      uint64_t newSize = pageCount * pageSize_;
      std::cerr << "Memory size (0x" << std::hex << memSize_ << ") is not a "
		<< "multiple of page size (0x" << pageSize_ << ") -- "
		<< "using 0x" << newSize << '\n' << std::dec;

      memSize_ = newSize;
    }

  // Whole memory is intially set for instruction/data/atomic access.
  // No iccm/dccm/mmr/io.
  pagePmas_.resize(pageCount);
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
