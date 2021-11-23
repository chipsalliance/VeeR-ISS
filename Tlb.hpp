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

namespace WdRiscv
{

  /// Translation lookaside buffer entry.
  struct TlbEntry
  {
    uint64_t virtPageNum_ = 0;
    uint64_t physPageNum_ = 0;
    uint64_t time_ = 0;      // Access time (we use order to approximate time).
    uint32_t asid_ = 0;      // Address space identifier.
    bool valid_ = false;
    bool global_ = false;    // 
    bool user_ = false;      // User-mode entry if true.
    bool read_ = false;      // Has read access.
    bool write_ = false;     // Write access.
    bool exec_ = false;      // Execute Access.
    bool accessed_ = false;
    bool dirty_ = false;
    uint8_t levels_ = 3;
  };


  /// Translation lookaside buffer.
  class Tlb
  {
  public:

    /// Define a a TLB with the given size (number of entries).
    Tlb(unsigned size);

    /// Return pointer to TLB entry associated with given virtual page
    /// number and address space identifier.  Return nullptr if no
    /// such entry.
    TlbEntry* findEntry(uint64_t pageNum, uint32_t asid)
    {
      if(useFullTlb_) {
    	  auto ait = fullTlb_[0].find((pageNum<<16) | uint64_t(asid));
    	  if(ait!=fullTlb_[0].end())
    		  return &ait->second;
    	  auto it = fullTlb_[1].find(pageNum<<16);
    	  if(it!=fullTlb_[0].end())
			  return &it->second;
    	  return nullptr;
      }
      for (auto& entry : entries_)
        if (entry.valid_ and entry.virtPageNum_ == pageNum)
          if (entry.global_ or entry.asid_ == asid)
            {
              entry.time_ = time_++;
              return &entry;
            }
      return nullptr;
    }
    /// print TLB content
    void printTlb(std::ostream& ost) const;
    /// print TLB entry
    void printEntry(std::ostream& ost, const TlbEntry& te) const;
    /// Insert a TLB entry for the given translation parameters. If TLB is full
    /// the contents of the  least recently accessed slot are replaced by the
    /// given parameters.
    void insertEntry(uint64_t virtPageNum, uint64_t phyPageNUm,
                     uint32_t asid, bool global, bool isUser, bool read,
                     bool write, bool exec);

    /// Insert copy of given entry.
    void insertEntry(const TlbEntry& entry);

    void invalidate()
    { for (auto& entry : entries_) entry.valid_ = false; }

  protected:

  private:
    const bool useFullTlb_;
    std::vector<TlbEntry> entries_;
    uint64_t time_ = 0;  // Access time (we use access order as approximation).
    std::unordered_map<uint64_t, TlbEntry> fullTlb_[2];
  };
}

    
