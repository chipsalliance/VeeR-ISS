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
    /// and not the tags). The total-size must be a power of 2 and
    /// must be a multiple of the line-size. The line-size must be a
    /// power of 2. To have a fully associated cache, arrange for the
    /// set-count to be the same as the line count (which is the
    /// total-size divided by the line-size).
    ///
    /// Typical line-size: 64
    /// Typical set-count: 16
    /// Typical total-size: 2*1024*1024  (2 MB)
    Cache(uint64_t totalSize, unsigned lineSize, unsigned setCount);

    /// Insert line overlapping given address into the cahce.
    void insert(uint64_t addr)
    {
      uint64_t lineNumber = getLineNumber(addr);
      uint64_t setIndex = getSetIndex(lineNumber);

      auto& lines = linesPerSet_.at(setIndex);
      auto& times = timesPerSet_.at(setIndex);

      // If line present update access time.
      auto lineIter = lines.find(lineNumber);
      bool present = lineIter != lines.end();
      if (present)
        {
          auto time = lineIter->second;
          auto timeIter = times.find(time);
          assert(timeIter != times.end() and timeIter->second == lineNumber);
          times.erase(timeIter);
          uint64_t nextTime = time_++;
          times[nextTime] = lineNumber;
          lines[lineNumber] = nextTime;
          return;
        }

      // If set is not full, insert line.
      if (lines.size() < setCount_)
        {
          assert(times.size() < setCount_);
          uint64_t nextTime = time_++;
          lines[lineNumber] = nextTime;
          times[nextTime] = lineNumber;
          return;
        }

      // Set is full. Evict line oldest line (one with smallest time).
      assert(lines.size() == setCount_);
      assert(times.size() == setCount_);

      auto timeIter = times.begin();
      auto evictedLine = timeIter->second;
      lineIter = lines.find(evictedLine);
      assert(lineIter != lines.end() and lineIter->second == timeIter->first);
      times.erase(timeIter);
      lines.erase(lineIter);

      uint64_t nextTime = time_++;
      lines[lineNumber] = nextTime;
      times[nextTime] = lineNumber;
    }

    /// Invalidate line overlapping given address.
    void invalidate(uint64_t addr)
    {
      uint64_t lineNumber = getLineNumber(addr);
      uint64_t setIndex = getSetIndex(lineNumber);

      auto& lines = linesPerSet_.at(setIndex);
      auto lineIter = lines.find(lineNumber);
      if (lineIter == lines.end())
        return; // Line not in cache
      auto time = lineIter->second;

      auto& times = timesPerSet_.at(setIndex);
      auto timeIter = times.find(time);
      assert(timeIter != times.end());
      assert(timeIter->second == lineNumber);

      times.erase(timeIter);
      lines.erase(lineIter);
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
      auto lineIter = lines.find(lineNumber);
      bool present = lineIter != lines.end();
      if (present)
        {
          uint64_t time = lineIter->second;
          auto& times = timesPerSet_.at(setIndex);
          auto timeIter = times.find(time);
          assert(timeIter != times.end() and timeIter->second == lineNumber);
          times.erase(timeIter);
          uint64_t nextTime = time_++;
          times[nextTime] = lineNumber;
          lines[lineNumber] = nextTime;
        }
      return present;
    }

  protected:

    /// Return the line number corresponding to the given address.
    uint64_t getLineNumber(uint64_t addr) const
    { return addr >> lineNumberShift_; }

    /// Cache is organized as an array of sets. Return the index of the
    /// set corresponding to the given line number.
    uint64_t getSetIndex(uint64_t lineNumber) const
    { return lineNumber >> setIndexShift_; }

  private:

    /// Map a line address to an access time.
    typedef std::map<uint64_t, uint64_t> LineToTime;

    /// Map access-time to a line address
    typedef std::map<uint64_t, uint64_t> TimeToLine;

    /// Map a set index (line-address modulo secCount) to an time-to-line
    /// map.
    std::vector<LineToTime> timesPerSet_;

    /// Map a set index (line-address modulo secCount) to a line-to-time
    /// map.
    std::vector<TimeToLine> linesPerSet_;

    uint64_t size_ = 0;
    uint64_t time_ = 0;
    unsigned lineSize_ = 0;
    unsigned setCount_ = 0;
    unsigned lineNumberShift_ = 0;
    unsigned setIndexShift_ = 0;
  };
}
