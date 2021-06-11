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
#include "Memory.hpp"


namespace WdRiscv
{

  template <typename URV>
  class Hart;

  template <typename URV>
  class Core
  {
  public:

    typedef Hart<URV> HartClass;

    /// Constructor: construct a core with n (hartsPerCore) harts
    /// assigning to them hart-ids (value in mhartid CSR) idBase to
    /// idBase + n - 1. CoreIx the index of this core in the system
    /// (cores are indexed 0 to m-1 where m is the number of cores).
    Core(URV idBase, unsigned coreIx, unsigned hartsPerCore, Memory& memory);

    ~Core();

    /// Return number of hrats in this core.
    unsigned hartCount() const
    { return harts_.size(); }

    /// Return pointer to ith hart in this core or null if i is
    /// greater-than or equal-to hartCount().
    std::shared_ptr<HartClass> ithHart(unsigned i)
    { return i < harts_.size() ? harts_.at(i) : std::shared_ptr<HartClass>(); }

  private:
    
    std::vector< std::shared_ptr<HartClass> > harts_;
  };

}
