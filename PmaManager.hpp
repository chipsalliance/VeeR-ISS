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

#pragma once

#include <cstdint>
#include <vector>
#include <unordered_map>

namespace WdRiscv
{

  /// Physical memory attribute. An instance of this is typically
  /// associated with a section of the address space. The address
  /// space is evenly divided into contiguous, equally sized sections,
  /// aligned to the section size.
  /// For sub-section attribution, an instance is associated with a
  /// word-aligned memory word. To reduce footprint of the PmaMgr
  /// object, we typically use a section size of 8 or more pages.
  class Pma
  {
  public:

    friend class PmaManager;

    enum Attrib
      {
       None = 0, Read = 1, Write = 2, Exec = 4,
       Idempotent = 8, Atomic = 16, Iccm = 32,
       Dccm = 64, MemMapped = 128,
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
    { return attrib_ & (Mapped | MemMapped); }

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

    /// Return true if in readable (ld instructions allowed) region.
    bool isRead() const
    { return attrib_ & (Read | MemMapped); }

    /// Return true if in writeable (st instructions allowed) region.
    bool isWrite() const
    { return attrib_ & (Write | MemMapped); }

    /// Return true if in executable (fetch allowed) region.
    bool isExec() const
    { return attrib_ & Exec; }

    /// Return true in region where atomic instructions are allowed.
    bool isAtomic() const
    { return attrib_ & Atomic; }

    /// Return true if this object has the same attributes as the
    /// given object.
    bool operator== (const Pma& other) const
    { return attrib_ == other.attrib_; }

    /// Return true if this object has different attributes from those
    /// of the given object.
    bool operator!= (const Pma& other) const
    { return attrib_ != other.attrib_; }

  private:

    unsigned attrib_ : 8;
    bool word_       : 8;     // True if word granularity otherwise section.
  } __attribute__((packed));


  /// Physical memory attribute manager. One per memory. Shared
  /// among cores and harts. Physical memory attributes apply to
  /// word-aligned regions as small as 1 word (but are expected to be
  /// applied to a few number of large regions).
  class PmaManager
  {
  public:

    friend class Memory;

    PmaManager(uint64_t memorySize, uint64_t sectionSize = 32*1024);

    /// Return the physical memory attribute associated with the
    /// word-aligned word designated by the given address. Return an
    /// unmapped attribute if the given address is out of memory
    /// range.
    Pma getPma(uint64_t addr) const
    {
      uint64_t ix = getSectionIx(addr);
      if (ix >= sectionPmas_.size())
        return Pma();
      Pma pma = sectionPmas_[ix];
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

    /// Associate a mask with the word-aligned word at the given address.
    void setMemMappedMask(uint64_t addr, uint32_t mask);

    /// Return mask associated with the word-aligned word at the given
    /// address.  Return 0xffffffff if no mask was ever associated
    /// with given address.
    uint32_t getMemMappedMask(uint64_t addr) const;

    /// Return true if the word-algined word containing given address
    /// is in data closed coupled memory.
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

    /// Internally, for a user specified region, we associate a pma
    /// object with each section of that region where the first/last
    /// address is aligned with the first/last address of a
    /// section. For a region where the first/last address is not
    /// section-aligned we associate a pma object with each word
    /// before/after the first/last section aligned address.

    /// Reset (to zero) all memory mapped registers.
    void resetMemMapped()
    { memMappedRegs_.assign(memMappedRegs_.size(), 0); }

    /// Bool an aread for memory mapped registers. Size is in bytes.
    /// Size must be a multiple of 4.
    bool defineMemMappedArea(uint64_t base, uint64_t size);

    /// Set value to the value of the memory mapped regiser at addr
    /// returning true if addr is valid. Return false if addr is not word
    /// aligned or is outside of the memory-mapped-regiser area.
    bool readRegister(uint64_t addr, uint32_t& value) const;

    /// Set value to the value of the memory mapped regiser byte at
    /// addr returning true if addr is valid. Return false if addr is
    /// is outside of the memory-mapped-regiser area.
    bool readRegisterByte(uint64_t addr, uint8_t& value) const
    {
      uint32_t word = 0;
      uint64_t wordAddr = (addr >> 2) << 2;
      if (not readRegister(wordAddr, word))
        return false;
      unsigned byteInWord = addr & 3;
      value = (word >> (byteInWord*8)) & 0xff;
      return true;
    }

    /// Set the value of the memory mapped regiser at addr to the
    /// given value returning true if addr is valid. Return false if
    /// addr is not a memory mapped reg leaving vlaue unmodified.
    bool writeRegister(uint64_t addr, uint32_t value);

    /// Similar to writeRgister but no masking is applied to value.
    bool writeRegisterNoMask(uint64_t addr, uint32_t value);

    /// Set the value of the memory mapped regiser byte at addr to the
    /// given value applying masking and returning true if addr is
    /// valid. Return false if addr is not a memory mapped reg leaving
    /// vlaue unmodified.
    bool writeRegisterByte(uint64_t addr, uint8_t value);

    /// Return start address of section containing given address.
    uint64_t getSectionStartAddr(uint64_t addr) const
    { return (addr >> sectionShift_) << sectionShift_; }

  private:

    /// Fracture attribute of section overlapping given address into word
    /// attributes.
    void fracture(uint64_t addr);

    /// Return the section number corresponding to the given address.
    uint64_t getSectionIx(uint64_t addr) const
    { return addr >> sectionShift_; }

  private:

    std::vector<Pma> sectionPmas_;
    std::unordered_map<uint64_t, Pma> wordPmas_; // Map word index to pma.
    uint64_t memSize_;
    uint64_t sectionSize_ = 32*1024;
    unsigned sectionShift_ = 15;

    uint64_t memMappedBase_ = 0;
    uint64_t memMappedSize_ = 0;
    std::vector<uint32_t> memMappedMasks_;  // One entry per register.
    std::vector<uint32_t> memMappedRegs_;  // One entry per register.
  };
}
