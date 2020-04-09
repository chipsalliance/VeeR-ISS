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

    friend class PmaManager;

    enum Attrib
      {
       None = 0, Read = 1, Write = 2, Exec = 4,
       Idempotent = 8, Atomic = 16, Iccm = 32,
       Dccm = 64, MemMapped = 128, Cached = 256,
       Aligned = 1024,
       ReadWrite = Read | Write,
       Mapped = Exec | Read | Write,
       Default = Mapped | Idempotent | Atomic
      };

    /// Default constructor: No access allowed. No-dccm, no-iccm,
    /// no-mmr, no-atomic.
    Pma(Attrib a = None)
      : attrib_(a), word_(false)
    { }

    /// Return true if mapped.
    bool isMapped() const
    { return attrib_ & Mapped; }

    /// Return true if in ICCM region (instruction closely coupled
    /// memory).
    bool isIccm() const
    { return attrib_ & Iccm; }

    /// Return true if in DCCM region (instruction closely coupled
    /// memory).
    bool isDccm() const
    { return attrib_ & Dccm; }

    /// Return true if in memory-mapped-register region.
    bool isMemMappedReg() const
    { return attrib_ & MemMapped; }

    /// Return true if in idempotent region.
    bool isIdempotent() const
    { return attrib_ & Idempotent; }

    /// Return true if in cacheable region.
    bool isCacheable() const     { return attrib_ & Attrib::Cached; }

    /// Return true if in readable (ld instructions allowed) region.
    bool isRead() const
    { return attrib_ & Read; }

    /// Return true if in writeable (st instructions allowed) region.
    bool isWrite() const
    { return attrib_ & Write; }

    /// Return true if in executable (fetch allowed) region.
    bool isExec() const
    { return attrib_ & Exec; }

    /// Return true in region where access must be aligned.
    bool isAligned() const
    { return attrib_ & Aligned; }

    /// Return true in region where atomic instructions are allowed.
    bool isAtomic() const
    { return attrib_ & Atomic; }

    /// Return true in cached region.
    bool isCached() const
    { return attrib_ & Cached; }

    /// Return true if this object has the same attributes as the
    /// given object.
    bool operator== (const Pma& other) const
    { return attrib_ == other.attrib_; }

    /// Return true if this object has different attributes from those
    /// of the given object.
    bool operator!= (const Pma& other) const
    { return attrib_ != other.attrib_; }

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

    friend class Memory;

    PmaManager(uint64_t memorySize, uint64_t pageSize);

    /// Return the physical memory attribute associated with the
    /// word-aligned word designated by the given address. Return an
    /// unmapped attribute if the given address is out of memory
    /// range. Internally we associate a pma object with each page of
    /// a regions where the first/last address is aligned with the
    /// first/last address of a page. For a region where the first/last
    /// address is not page-aligned we associate a pma object with
    /// each word before/after the first/last page aligned address.
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

    /// Set attribute of word-aligned words overlapping given region.
    void setAttribute(uint64_t addr0, uint64_t addr1, Pma::Attrib attrib);

    /// Return start address of page containing given address.
    uint64_t getPageStartAddr(uint64_t addr) const
    { return (addr >> pageShift_) << pageShift_; }

    /// Associate a mask with the word-aligned word at the given address.
    void setMemMappedMask(uint64_t addr, uint32_t mask);

    /// Return mask associated with the word-aligned word at the given
    /// address.  Return 0xffffffff if no mask was ever associated
    /// with given address.
    uint32_t getMemMappedMask(uint64_t addr) const;

    /// Return true if page of given address is in data closed coupled
    /// memory.
    bool isAddrInDccm(size_t addr) const
    { Pma pma = getPma(addr); return pma.isDccm(); }

    /// Return true if given address is in memory-mapped register region.
    bool isAddrMemMapped(size_t addr) const
    { Pma pma = getPma(addr); return pma.isMemMappedReg(); }

    /// Change the base address of the memory mapped register area to
    /// the given new base. Return true on success and false if new
    /// base cannot accomodate the size of the memory mapped area (new
    /// base too close to end of memory).
    bool changeMemMappedBase(uint64_t newBase);

  protected:

    /// Reset (to zero) all memory mapped registers.
    void resetMemMapped()
    { memMappedRegs_.assign(memMappedRegs_.size(), 0); }

    /// Bool an aread for memory mapped registers. Size is in bytes.
    /// Size must be a multiple of 4.
    bool defineMemMappedArea(uint64_t base, uint64_t size);

    /// Set value to the value of the memory mapped regiser at addr
    /// returning true if addr is valid. Return false if addr is not a
    /// memory mapped reg leaving vlaue unmodified.
    bool readRegister(uint64_t addr, uint32_t& value) const;

    /// Set the value of the memory mapped regiser at addr to the
    /// given value returning true if addr is valid. Return false if
    /// addr is not a memory mapped reg leaving vlaue unmodified.
    bool writeRegister(uint64_t addr, uint32_t value);

    /// Similar to writeRgister but no masking is applied to value.
    bool writeRegisterNoMask(uint64_t addr, uint32_t value);

    /// Set the value of the memory mapped regiser byte at addr to the
    /// given value returning true if addr is valid. Return false if
    /// addr is not a memory mapped reg leaving vlaue unmodified.
    bool writeRegisterByte(uint64_t addr, uint8_t value);

  private:

    /// Fracture attribute of page overlapping given address into word
    /// attributes.
    void fracture(uint64_t addr);

    /// Return the page number corresponding to the given address.
    uint64_t getPageIx(uint64_t addr) const
    { return addr >> pageShift_; }

  private:

    std::vector<Pma> pagePmas_;
    std::unordered_map<uint64_t, Pma> wordPmas_; // Map word index to pma.
    uint64_t memSize_;
    uint64_t pageSize_ = 4*1024;
    unsigned pageShift_ = 12;

    uint64_t memMappedBase_ = 0;
    uint64_t memMappedSize_ = 0;
    std::vector<uint32_t> memMappedMasks_;  // One entry per register.
    std::vector<uint32_t> memMappedRegs_;  // One entry per register.
  };
}
