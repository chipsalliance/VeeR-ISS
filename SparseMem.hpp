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
#include <cstring>
#include <unordered_map>


namespace WdRiscv
{

  /// Memory model. Host machine memory is conserved by allocating
  /// pages only for addresses refereneced.
  class SparseMem
  {
  public:

    SparseMem();

    ~SparseMem();

    /// Read an unsigned item of the given size (1, 2, 4 or 8 bytes)
    /// from the given target-machine address placing the item bits
    /// in value. Return true on success and false on failure (out of
    /// bounds address or size).
    bool read(uint64_t addr, unsigned size, uint64_t& value);

    /// Write an unsigned item of the given size (1, 2, 4 or 8 bytes)
    /// to the given target-machine address gettng the item bits from
    /// the least sig bits of value. Return true on success and false
    /// on failure (out of bounds address or size).
    bool write(uint64_t addr, unsigned size, uint64_t value);

    /// Write the contents of the memory to a verilog hex file. Return
    /// true on success and false on failure.
    bool writeHexFile(const std::string& path) const;

  protected:

    /// Read from given target-machine address an item of type U
    /// (uint8_t, uint16_t, uint32_t or uin64_t) and place it in
    /// value. Address must be aligned.
    template <typename U>
    bool
    read(uint64_t addr, uint64_t& value)
    {
      uint64_t pageRank = getPageRank(addr);
      uint8_t* page = findOrCreatePage(pageRank);
      if (not page)
        return false;
      unsigned offset = addr & pageMask_;
      value = *( reinterpret_cast<U*>(page + offset) );
      return true;
    }

    /// Write to the given target-machine address an item of type U
    /// (uint8_t, uint16_t, uint32_t or uin64_t). Address must be
    /// aligned. Item bits are in the least significant bits of
    /// value.
    template <typename U>
    bool
    write(uint64_t addr, uint64_t value)
    {
      uint64_t pageRank = getPageRank(addr);
      uint8_t* page = findOrCreatePage(pageRank);
      if (not page)
        return false;
      unsigned offset = addr & pageMask_;
      *( reinterpret_cast<U*>(page + offset) ) = value;
      return true;
    }

    /// Return the page number of the page containing the byte at the
    /// given address.
    uint64_t getPageRank(uint64_t addr)
    { return addr >> pageShift_; }

    /// Return host-machine address of the target-machine page with
    /// the given page number creating such a page (and zeroing it) if
    /// it has never been accessed before.
    uint8_t* findOrCreatePage(uint64_t pageRank)
    {
      auto iter = pageMap_.find(pageRank);
      if (iter != pageMap_.end())
        return iter->second;
      uint8_t* page = new uint8_t[pageSize_];
      memset(page, 0, pageSize_);
      pageMap_[pageRank] = page;
      return page;
    }

  private:

    size_t pageSize_ = 4*1024;
    unsigned pageShift_ = 12;
    unsigned pageMask_ = 0xfff;

    std::unordered_map<uint64_t, uint8_t*> pageMap_;  // Map address to page
  };
}
