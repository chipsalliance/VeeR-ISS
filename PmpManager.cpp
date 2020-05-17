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
#include <fstream>
#include "PmpManager.hpp"

using namespace WdRiscv;


PmpManager::PmpManager(uint64_t memSize, uint64_t pageSize)
  : memSize_(memSize), pageSize_(pageSize), accessCount_(16),
    typeCount_(Pmp::Type::_Count)
{
  assert(memSize >= pageSize);
  assert(pageSize >= 64);

  uint64_t logPageSize = static_cast<uint64_t>(std::log2(pageSize_));
  uint64_t p2PageSize = uint64_t(1) << logPageSize;
  assert(p2PageSize == pageSize_);
  pageShift_ = logPageSize;

  uint64_t pageCount = memSize_ / pageSize_;
  assert(pageCount * pageSize_ == memSize_);

  // Mark memory as no access (machine mode still has access because
  // it is not checked).
  pagePmps_.resize(pageCount, Pmp::None);
}


PmpManager::~PmpManager()
{
}


void
PmpManager::reset()
{
  for (auto& entry : pagePmps_)
    entry = Pmp(Pmp::None);
  wordPmps_.clear();
}


void
PmpManager::setMode(uint64_t a0, uint64_t a1, Pmp::Type type, Pmp::Mode mode,
                    unsigned pmpIx, bool lock)
{
  if (a1 >= memSize_)
    a1 = memSize_ - 1;

  a0 = (a0 >> 2) << 2;   // Make word aligned.
  a1 = (a1 >> 2) << 2;   // Make word aligned.

  while (a0 <= a1)
    {
      uint64_t p0 = getPageStartAddr(a0);
      uint64_t p1 = getPageStartAddr(a1);

      bool doWord = (a0 != p0) or (p0 == p1 and (a1 - a0 + 4 != pageSize_));
      Pmp prev = getPmp(a0);
      doWord = doWord or prev.word_;

      Pmp pmp(mode, pmpIx, lock, type);

      if (doWord)
        {
          fracture(a0);
          pmp.word_ = true;
          uint64_t last = std::min(a1, p0 + pageSize_ - 4);
          for ( ; a0 <= last; a0 += 4)
            wordPmps_[a0>>2] = pmp;
        }
      else
        {
          uint64_t pageIx = getPageIx(p0);
          pagePmps_.at(pageIx) = pmp;
          a0 += pageSize_;
        }
    }
}


bool
PmpManager::printStats(std::ostream& out) const
{
  out << "PMP entry access count:\n";
  for (size_t i = 0; i < accessCount_.size(); ++i)
    out << "  entry " << i << ": " << accessCount_.at(i) << '\n';

  out << "PMP type count:\n";
  for (size_t i = 0; i < typeCount_.size(); ++i)
    {
      out << "  ";
      switch(i)
        {
        case Pmp::Type::Off:   out << "OFF:   "; break;
        case Pmp::Type::Tor:   out << "TOR:   "; break;
        case Pmp::Type::Na4:   out << "NA4:   "; break;
        case Pmp::Type::Napot: out << "NAPOT: "; break;
        default:               out << "??:    "; break;
        }
      out << typeCount_.at(i) << '\n';
    }

  return true;
}


bool
PmpManager::printStats(const std::string& path) const
{
  std::ofstream out(path);
  if (not out.good())
    {
      std::cerr << "Failed to open file '" << path << "' for output\n";
      return false;
    }

  return printStats(out);
}
