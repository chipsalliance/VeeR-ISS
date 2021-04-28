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
#include <string>
#include <unordered_map>
#include "CsRegs.hpp"

namespace WdRiscv
{

  /// Physical memory protection. An instance of this is usually
  /// associated with a section of the address space. The address
  /// space is evenly divided into contiguous, equally sized sections,
  /// aligned to the section size.
  /// For sub-section protection, an instance is associated with a
  /// word-aligned memory word. To reduce footprint of the PmaMgr
  /// object, we typically use a section size of 8 or more pages.
  class Pmp
  {
  public:

    friend class PmpManager;

    enum Type { Off = 0, Tor = 1, Na4 = 2, Napot = 3, _Count = 4 };

    enum Mode
      {
       None = 0, Read = 1, Write = 2, Exec = 4, ReadWrite = Read | Write,
       Default = Read | Write | Exec
      };

    /// Default constructor: No access allowed.
    Pmp(Mode m = None, unsigned pmpIx = 0, bool locked = false,
        Type type = Type::Off)
      : mode_(m), type_(type), locked_(locked), pmpFree_(false),
        pmpIx_(pmpIx), word_(false)
    { }

    /// Return true if read (i.e. load instructions) access allowed 
    bool isRead(PrivilegeMode mode, PrivilegeMode prevMode, bool mprv) const
    {
      if (mprv)
        mode = prevMode;
      bool check = (mode != PrivilegeMode::Machine) or locked_ or pmpFree_;
      return check ? mode_ & Read : true;
    }

    /// Return true if write (i.e. store instructions) access allowed.
    bool isWrite(PrivilegeMode mode, PrivilegeMode prevMode, bool mprv) const
    {
      bool check = (mode != PrivilegeMode::Machine or locked_ or
                    (mprv and prevMode != PrivilegeMode::Machine));
      return check ? mode_ & Write : true;
    }

    /// Return true if instruction fecth is allowed.
    bool isExec(PrivilegeMode mode, PrivilegeMode prevMode, bool mprv) const
    {
      bool check = (mode != PrivilegeMode::Machine or locked_ or
                    (mprv and prevMode != PrivilegeMode::Machine));
      return check ? mode_ & Exec : true;
    }

    /// Return true if this object has the mode attributes as the
    /// given object.
    bool operator== (const Pmp& other) const
    { return mode_ == other.mode_ and pmpIx_ == other.pmpIx_; }

    /// Return true if this object has different attributes from those
    /// of the given object.
    bool operator!= (const Pmp& other) const
    { return mode_ != other.pmpIx_ or pmpIx_ != other.pmpIx_; }

    /// Return string representation of the given PMP type.
    static std::string toString(Type t);

    /// Return string representation of the given PMP mode.
    static std::string toString(Mode m);
      

  protected:

    /// Return the index of the PMP entry from which this object was
    /// created.
    unsigned pmpIndex() const
    { return pmpIx_; }

  private:

    uint8_t mode_ = 0;
    Type type_      : 8;
    bool locked_    : 1;
    bool pmpFree_   : 1;  // Not covered by any pmp register.
    unsigned pmpIx_ : 5;  // Index of corresponding pmp register.
    bool word_      : 1;
  } __attribute__((packed));


  /// Physical memory protection manager. One per hart.  rotection
  /// applies to word-aligned regions as small as 1 word but are
  /// expected to be applied to a small number (less than or equal 16)
  /// of large regions.
  class PmpManager
  {
  public:

    friend class Memory;

    /// Constructor. Mark all memory as no access to user/supervisor
    /// (machine mode does access since it is not checked).
    PmpManager(uint64_t memorySize, uint64_t sectionSize = 32*1024);

    /// Destructor.
    ~PmpManager();

    /// Reset: Mark all memory as no access to user/supervisor
    /// (machine mode does have access because it is not checked).
    void reset();

    /// Return the physical memory protection object (pmp) associated
    /// with the word-aligned word designated by the given
    /// address. Return a no-access object if the given address is out
    /// of memory range.
    Pmp getPmp(uint64_t addr) const
    {
      uint64_t ix = getSectionIx(addr);
      if (ix >= sectionPmps_.size())
        return Pmp();
      Pmp pmp = sectionPmps_[ix];
      if (pmp.word_)
        {
          addr = (addr >> 2);  // Get word index.
          pmp = wordPmps_.at(addr);
        }
      return pmp;
    }

    /// Similar to getPmp but it also updates the access count associated with
    /// each PMP entry.
    inline Pmp accessPmp(uint64_t addr) const
    {
      uint64_t ix = getSectionIx(addr);
      if (ix >= sectionPmps_.size())
        return Pmp();
      Pmp pmp = sectionPmps_[ix];
      if (pmp.word_)
        {
          addr = (addr >> 2);  // Get word index.
          pmp = wordPmps_.at(addr);
        }
      accessCount_.at(pmp.pmpIndex())++;
      typeCount_.at(pmp.type_)++;
      return pmp;
    }

    /// Set access mode of word-aligned words overlapping given region
    /// for user/supervisor.
    void setMode(uint64_t addr0, uint64_t addr1, Pmp::Type type, Pmp::Mode mode,
                 unsigned pmpIx, bool locked);

    /// Return start address of section containing given address.
    uint64_t getSectionStartAddr(uint64_t addr) const
    { return (addr >> sectionShift_) << sectionShift_; }

    /// Print statistics on the given stream.
    bool printStats(std::ostream& out) const;

    /// Print statistics on the given file.
    bool printStats(const std::string& path) const;

  private:

    /// Internally, for a user specified region, we associate a pmp
    /// object with each section of that region where the first/last
    /// address is aligned with the first/last address of a
    /// section. For a region where the first/last address is not
    /// section-aligned we associate a pmp object with each word
    /// before/after the first/last section aligned address.

    /// Fracture attribute of section overlapping given address into
    /// word attributes.
    void fracture(uint64_t addr)
    {
      uint64_t sectionIx = getSectionIx(addr);
      if (sectionIx > sectionPmps_.size())
        return;

      Pmp pmp = sectionPmps_.at(sectionIx);
      if (pmp.word_)
        return;
      pmp.word_= true;
      sectionPmps_.at(sectionIx) = pmp;

      uint64_t words = sectionSize_ / 4;
      uint64_t wordIx = (sectionIx*sectionSize_) >> 2;
      for (uint64_t i = 0; i < words; ++i, wordIx++)
        wordPmps_[wordIx] = pmp;
    }

    uint64_t getSectionIx(uint64_t addr) const
    { return addr >> sectionShift_; }

  private:

    std::vector<Pmp> sectionPmps_;
    std::unordered_map<uint64_t, Pmp> wordPmps_; // Map word index to pmp.
    uint64_t memSize_;
    uint64_t sectionSize_ = 32*1024;
    unsigned sectionShift_ = 15;
    mutable std::vector<uint64_t> accessCount_;  // PMP entry access count.
    mutable std::vector<uint64_t> typeCount_;  // PMP type access count.
  };
}
