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
    uint64_t time_ = 0;   // Access time (we use order to approximate time).
    uint32_t asid_ = 0;   // Address space identifier.
    bool user_ = false;   // User-mode entry if trye.
    bool valid_ = false;
    bool read_ = false;   // Has read access.
    bool write_ = false;  // Write access.
    bool exec_ = false;   // Execute Access.
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
    const TlbEntry* findEntry(uint64_t pageNum, uint32_t asid)
    {
      for (auto& entry : entries_)
        if (entry.valid_ and entry.asid_ == asid and entry.virtPageNum_ == pageNum)
          {
            entry.time_ = time_++;
            return &entry;
          }
      return nullptr;
    }

    /// Insert a TLB entry for the given translation parameters. If TLB is full
    /// the contents of the  least recently accessed slot are replaced by the
    /// given parameters.
    void insertEntry(uint64_t virtPageNum, uint64_t phyPageNUm,
                     uint32_t asid, bool isUser, bool read, bool write, bool exec);

  private:

    std::vector<TlbEntry> entries_;
    uint64_t time_ = 0;  // Access time (we use access order as approximation).
  };
}

    
