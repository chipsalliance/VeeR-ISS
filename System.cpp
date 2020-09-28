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

#include "Hart.hpp"
#include "Core.hpp"
#include "System.hpp"

using namespace WdRiscv;


template <typename URV>
System<URV>::System(unsigned coreCount, unsigned hartsPerCore, size_t memSize,
                    size_t pageSize)
  : hartCount_(coreCount * hartsPerCore), hartsPerCore_(hartsPerCore)
{
  cores_.resize(coreCount);

  memory_ = std::make_shared<Memory>(memSize, pageSize);

  Memory& mem = *(memory_.get());
  mem.setHartCount(hartCount_);

  for (unsigned ix = 0; ix < coreCount; ++ix)
    {
      URV hartIdBase = ix * hartsPerCore;
      cores_.at(ix) = std::make_shared<CoreClass>(hartIdBase, hartsPerCore, mem);

      // Maintain a vector of all the harts in the system.
      auto core = cores_.at(ix);
      for (unsigned i = 0; i < hartsPerCore; ++i)
        {
          auto hart = core->ithHart(i);
          sysHarts_.push_back(hart);
        }
    }
}


template <typename URV>
System<URV>::~System()
{
}


template <typename URV>
void
System<URV>::checkUnmappedElf(bool flag)
{
  if (memory_)
    memory_->checkUnmappedElf(flag);
}


template class WdRiscv::System<uint32_t>;
template class WdRiscv::System<uint64_t>;
