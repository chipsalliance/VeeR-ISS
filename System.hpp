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

#include <memory>               // For shared_ptr
#include <functional>
#include "Memory.hpp"
#include "SparseMem.hpp"

namespace WdRiscv
{

  template <typename URV>
  class Hart;

  template <typename URV>
  class Core;


  /// Model a system consisting of n cores with m-harts per core and a
  /// memory. The harts in the system are indexed from 0 to n*m -
  /// 1. The type URV (unsigned register value) is that of the integer
  /// register and is either uint32_t or uint64_t.
  template <typename URV>
  class System
  {
  public:

    typedef Hart<URV> HartClass;
    typedef Core<URV> CoreClass;

    /// Constructor: Construct a system with n (coreCount) cores each
    /// consisting of m (hartsPerCore) harts. The harts in this system
    /// are indexed with 0 to n*m - 1.  Each core is assigned a
    /// hart-id start from the sequence 0, hartIdOffset,
    /// 2*hartIdOffset, ...  Harts in a core are assigned consecutive
    /// hart-ids (values of MHARTID CSRs) starting with the start if
    /// od the core.
    System(unsigned coreCount, unsigned hartsPerCore, unsigned hartIdOffset,
           size_t memSize, size_t pageSize, size_t regionSize);

    ~System();

    /// Return count of cores in this system.
    unsigned coreCount() const
    { return cores_.size(); }

    /// Return the number of harts per core.
    unsigned hartsPerCore() const
    { return hartsPerCore_; }

    /// Return count of harts (coreCount * hartsPerCore) in this
    /// system.
    unsigned hartCount() const
    { return hartCount_; }

    /// Return pointer to the ith hart in the system or null if i is
    /// out of bounds. A hart index is valid if it is less than the
    /// value returned by the hartCount method.
    std::shared_ptr<HartClass> ithHart(unsigned i) const
    {
      if (i >= sysHarts_.size())
	return std::shared_ptr<HartClass>();
      return sysHarts_.at(i);
    }

    /// Return pointer to this system hart having the given value as
    /// its hart-id (value of MHARTID CSR) or null if no such hart.
    std::shared_ptr<HartClass> findHartByHartId(URV hartId) const
    {
      const auto& iter = hartIdToIndex_.find(hartId);
      if (iter != hartIdToIndex_.end())
        return ithHart(iter->second);
      return std::shared_ptr<HartClass>();
    }

    /// Return pointer to the ith core in the system or null if i is
    /// out of bounds.
    std::shared_ptr<CoreClass> ithCore(unsigned i)
    {
      if (i >= cores_.size())
	return std::shared_ptr<CoreClass>();
      return cores_.at(i);
    }

    /// With a true flag, when loading ELF files, error out if ELF
    /// file refers to unmapped memory. With a false flag, ignore
    /// unmapped memory in the ELF file.
    void checkUnmappedElf(bool flag);

    /// Define read memory callback. This (along with
    /// defineWriteMemoryCallback) allows the caller to bypass the
    /// memory model with their own.
    void defineReadMemoryCallback(
         std::function<bool(uint64_t, unsigned, uint64_t&)> callback )
    {
      memory_->defineReadMemoryCallback(callback);
    }

    /// Define write memory callback. This (along with
    /// defineReadMemoryCallback) allows the caller to bypass the
    /// memory model with their own.
    void defineWriteMemoryCallback(
         std::function<bool(uint64_t, unsigned, uint64_t)> callback )
    {
      memory_->defineWriteMemoryCallback(callback);
    }

    /// Break a hart-index-in-system into a core-index and a
    /// hart-index in core. Return true if successful and false if
    /// igven hart-index-in-system is out of bounds.
    bool unpackSystemHartIx(unsigned hartIxInSys, unsigned& coreIx,
                            unsigned& hartIxInCore)
    {
      if (hartIxInSys >= sysHarts_.size())
        return false;
      coreIx = hartIxInSys / hartsPerCore_;
      hartIxInCore = hartIxInSys % hartsPerCore_;
      return true;
    }

    /// Write contents of memory accessed by current run in verilog
    /// hex format to the file at the given path. Return true on
    /// success and false on failure. Currently this will write the
    /// contents of accessed pages.
    bool writeAccessedMemory(const std::string& path) const;

  private:

    unsigned hartCount_;
    unsigned hartsPerCore_;

    std::vector< std::shared_ptr<CoreClass> > cores_;
    std::vector< std::shared_ptr<HartClass> > sysHarts_; // All harts in system.
    std::unordered_map<URV, unsigned> hartIdToIndex_;
    std::shared_ptr<Memory> memory_ = nullptr;
    SparseMem* sparseMem_ = nullptr;
  };
}
