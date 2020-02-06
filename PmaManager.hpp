//
// SPED-License-Identifier: GPL-3.0-or-later
// Copyright 2018-2020 Western Digital Corporation or its affiliates.
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

#pragma once

#include <cstdint>
#include <vector>
#include <unordered_map>

namespace WdRiscv
{

  /// Physical memory attribute. An instance of this is usually
  /// associated with a memory page. For sub-page attribution, an
  /// instance is associated with a word-aligned memory word.
  class Pma
  {
  public:

    enum Attrib { None = 0, Inst = 1, Data = 2, Mapped = 3,
                  Idempotent = 4, Atomic = 8, Iccm = 16,
                  Dccm = 32, MemMapped = 64, Cached = 128 };

    friend class PmaManager;

    /// Default constructor: All access allowed. No-dccm, no-iccm,
    /// no-mmr, atomic.
    Pma()
      : attrib_(uint16_t(Attrib::Inst) | uint16_t(Attrib::Data) |
                uint16_t(Attrib::Idempotent) | uint16_t(Attrib::Atomic))
    { }

    bool isMapped() const        { return attrib_ & Attrib::Mapped; }
    bool isIccm() const          { return attrib_ & Attrib::Iccm; }
    bool isDccm() const          { return attrib_ & Attrib::Dccm; }
    bool isMemMappedReg() const  { return attrib_ & Attrib::MemMapped; }
    bool isIdempotent() const    { return attrib_ & Attrib::Idempotent; }
    bool isCacheable() const     { return attrib_ & Attrib::Cached; }

    /// True if instruction fech is allowed at word-aligned word
    /// covering given address.
    bool hasInst() const         { return attrib_ & Attrib::Inst; }

    /// True if load/store is allowed at word-aligned word covering
    /// given address.
    bool hasData() const         { return attrib_ & Attrib::Data; }

    /// True if atomic instructions allowed at word-aligned word
    /// covering given address.
    bool hasAtomic() const       { return attrib_ & Attrib::Atomic; }

  private:

    uint16_t attrib_ = 0;
    bool word_ = false;     // True if word granularity otherwise page.
  };


  /// Physical memory attribute manager. One per memory. Shared
  /// between cores and harts. Physical memory attributes apply to
  /// word-aligned regions as small as 1 word (but are expected to be
  /// applied to a few number of large regions).
  class PmaManager
  {
  public:
    PmaManager(uint64_t memorySize, uint64_t pageSize);

    /// Return the physical memory attribute associated with the
    /// word-aligned word designated by the given address. Return an
    /// unmapped attribute if the given address is out of memory
    /// range.
    Pma getPma(uint64_t addr) const
    {
      uint64_t ix = getPageIx(addr);
      if (ix >= pagePmas_.size())
        return Pma();
      Pma pma = pagePmas_[ix];
      if (pma.word_)
        {
          addr = (addr >> 2);  // Get word index.
          pma = wordPmas_.at(addr);
        }
      return pma;
    }

    /// Enable given attribute in word-aligned words overlapping given
    /// region.
    void enable(uint64_t addr0, uint64_t addr1, Pma::Attrib attrib);

    /// Disable given attribute in word-aligned words overlapping given
    /// region.
    void disable(uint64_t addr0, uint64_t addr1, Pma::Attrib attrib);

    /// Return start address of page containing given address.
    size_t getPageStartAddr(size_t addr) const
    { return (addr >> pageShift_) << pageShift_; }

  private:

    /// Fracture attribute of page overlapping given address into word
    /// attributes.
    void fracture(uint64_t addr)
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

    uint64_t getPageIx(uint64_t addr) const
    { return addr >> pageShift_; }

  private:

    std::vector<Pma> pagePmas_;
    std::unordered_map<size_t, Pma> wordPmas_; // Map word index to pma.
    uint64_t memSize_;
    uint64_t pageSize_ = 4*1024;
    unsigned pageShift_ = 12;
  };
}
