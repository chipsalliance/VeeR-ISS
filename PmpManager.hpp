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

  /// Physical memory protection. An instance of this is usually
  /// associated with a memory page. For sub-page attribution, an
  /// instance is associated with a word-aligned memory word.
  class Pmp
  {
  public:

    friend class PmpManager;

    enum Mode
      {
       None = 0, Exec = 1, Read = 2, Write = 4, ReadWrite = Read | Write,
       Default = Exec | Read | Write
      };

    /// Default constructor: No access allowed.
    Pmp(Mode m = None)
      : mode_(m), word_(false)
    { }

    /// Return true if read (i.e. load instructions) access allowed 
    bool isRead() const
    { return mode_ & Read; }

    /// Return true if write (i.e. store instructions) access allowed.
    bool isWrite() const
    { return mode_ & Write; }

    /// Return true if instruction fecth is allowed.
    bool isExec() const
    { return mode_ & Exec; }

    /// Return true if this object has the mode attributes as the
    /// given object.
    bool operator== (const Pmp& other) const
    { return mode_ == other.mode_; }

    /// Return true if this object has different attributes from those
    /// of the given object.
    bool operator!= (const Pmp& other) const
    { return mode_ != other.mode_; }

  private:

    uint8_t mode_ = 0;
    bool word_ = false;     // True if word granularity otherwise page.
  } __attribute__((packed));


  /// Physical memory protection manager. One per hart.  Physical
  /// memory attributes apply to word-aligned regions as small as 1
  /// word (but are expected to be applied to a few number of large
  /// regions).
  class PmpManager
  {
  public:

    friend class Memory;

    PmpManager(uint64_t memorySize, uint64_t pageSize);

    /// Return the physical memory protection object (pmm) associated
    /// with the word-aligned word designated by the given
    /// address. Return a no-access object if the given address is out
    /// of memory range. Internally we associate a pmp object with
    /// each page of a region where the first/last address is aligned
    /// with the first/last address of a page. For a region where the
    /// first/last address is not page-aligned we associate a pma
    /// object with each word before/after the first/last page aligned
    /// address.
    Pmp getPmp(uint64_t addr) const
    {
      uint64_t ix = getPageIx(addr);
      if (ix >= pagePmps_.size())
        return Pmp();
      Pmp pmp = pagePmps_[ix];
      if (pmp.word_)
        {
          addr = (addr >> 2);  // Get word index.
          pmp = wordPmps_.at(addr);
        }
      return pmp;
    }

    /// Enable given mode in word-aligned words overlapping given
    /// region.
    void enable(uint64_t addr0, uint64_t addr1, Pmp::Mode mode);

    /// Disable given mode in word-aligned words overlapping given
    /// region.
    void disable(uint64_t addr0, uint64_t addr1, Pmp::Mode mode);

    /// Set access mode of word-aligned words overlapping given region.
    void setMode(uint64_t addr0, uint64_t addr1, Pmp::Mode mode);

    /// Return start address of page containing given address.
    uint64_t getPageStartAddr(uint64_t addr) const
    { return (addr >> pageShift_) << pageShift_; }

  private:

    /// Fracture attribute of page overlapping given address into word
    /// attributes.
    void fracture(uint64_t addr)
    {
      uint64_t pageIx = getPageIx(addr);
      Pmp pmp = pagePmps_.at(pageIx);
      if (pmp.word_)
        return;
      pmp.word_= true;

      uint64_t words = pageSize_ / 4;
      uint64_t wordIx = (pageIx*pageSize_) >> 2;
      for (uint64_t i = 0; i < words; ++i, wordIx++)
        wordPmps_[wordIx] = pmp;
    }

    uint64_t getPageIx(uint64_t addr) const
    { return addr >> pageShift_; }

  private:

    std::vector<Pmp> pagePmps_;
    std::unordered_map<uint64_t, Pmp> wordPmps_; // Map word index to pma.
    uint64_t memSize_;
    uint64_t pageSize_ = 4*1024;
    unsigned pageShift_ = 12;
  };
}
