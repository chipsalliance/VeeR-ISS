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
#include <map>
#include <unordered_map>
#include <vector>
#include <cassert>


namespace WdRiscv
{

  /// Model a cache. This is for the support of the performance model.
  /// We keep track of the addresses of the lines in the cache.  We
  /// do not keep track of the data.
  class Cache
  {
  public:

    /// Replacement policy.
    enum Policy { RANDOM, LRU };

    /// Define a cache with the given total data size and line size
    /// (all sizes in bytes and refer to the data part of the cache
    /// and not the tags) and set-associativity. The total-size must
    /// be a power of 2 and must be a multiple of the line-size. The
    /// line-size and set-size must also be powers of 2. Peformance
    /// will degrade significanlty (quadratic cost) if set size is
    /// larger than 64.
    ///
    /// Typical total-size: 2*1024*1024  (2 MB)
    /// Typical line-size: 64 bytes
    /// Typical set-size: 16 (16-way set associative)
    Cache(uint64_t totalSize, unsigned lineSize, unsigned setCount);

    ~Cache();

    /// Insert line overlapping given address into the cahce.
    void insert(uint64_t addr)
    {
      uint64_t lineNumber = getLineNumber(addr);
      uint64_t setIndex = getSetIndex(lineNumber);
      auto& lines = linesPerSet_.at(setIndex);
      accesses_++;

      // Find line number or oldest entry.
      size_t bestIx = 0;
      for (size_t ix = 0; ix < lines.size(); ++ix)
        {
          auto& entry = lines[ix];
          if (entry.tag_ == lineNumber)
            {
              hits_++;
              bestIx = ix;
              break;
            }
          if (not entry.valid() or entry.time_ < lines[bestIx].time_)
            bestIx = ix;
        }
      lines[bestIx].tag_ = lineNumber;
      lines[bestIx].time_ = time_++;
    }

    /// Invalidate line overlapping given address.
    void invalidate(uint64_t addr)
    {
      uint64_t lineNumber = getLineNumber(addr);
      uint64_t setIndex = getSetIndex(lineNumber);

      auto& lines = linesPerSet_.at(setIndex);
      for (auto& entry : lines)
        if (entry.tag_ == lineNumber)
          entry.invalidate();
    }

    /// Return true if line overlapping given address is present in
    /// the cache updating the access time of the line; return flase
    /// otherwise. A line is present after it is inserted and until it
    /// is evicted or invalidated.
    bool access(uint64_t addr)
    {
      uint64_t lineNumber = getLineNumber(addr);
      uint64_t setIndex = getSetIndex(lineNumber);
      auto& lines = linesPerSet_.at(setIndex);
      accesses_++;

      for (auto& entry : lines)
        if (entry.tag_ == lineNumber)
          {
            hits_++;
            entry.time_ = time_++;
            return true;
          }
      return false;
    }

    /// Fill the given vector (cleared on entry) with the addresses of
    /// the lines curently in the cache in descending order (oldest
    /// one first) by age.
    void getLineAddresses(std::vector<uint64_t>& result) const;

    /// Take a snapshot of the cache tags into the given file. Return
    /// true on success or false on failure
    bool saveSnapshot(const std::string& path);

    /// Load the cache tags from the snapshot file. Return true on
    /// success and fase on failure.
    bool loadSnapshot(const std::string& path);

  protected:

    /// Return the line number corresponding to the given address.
    uint64_t getLineNumber(uint64_t addr) const
    { return addr >> lineNumberShift_; }

    /// Cache is organized as an array of sets. Return the index of the
    /// set corresponding to the given line number.
    uint64_t getSetIndex(uint64_t lineNumber) const
    { return lineNumber & setIndexMask_; }

  private:

    struct Entry
    {
      uint64_t tag_ = ~uint64_t(0);
      uint64_t time_ = 0;

      bool valid() const { return tag_ != ~uint64_t(0); }
      void invalidate() { tag_ = ~uint64_t(0); }
    };

    /// Lines in set: Map a line address to an access time.
    //typedef std::unordered_map<uint64_t, uint64_t> LineToTime;
    typedef std::vector<Entry> LinesInSet;

    /// Map a set index (line-address modulo secCount) to a line-to-time
    /// map.
    std::vector<LinesInSet> linesPerSet_;

    uint64_t size_ = 0;
    uint64_t time_ = 0;
    unsigned lineSize_ = 0;
    unsigned setSize_ = 0;
    unsigned lineNumberShift_ = 0;
    unsigned setIndexMask_ = 0;

    uint64_t hits_ = 0;
    uint64_t accesses_ = 0;
  };
}
