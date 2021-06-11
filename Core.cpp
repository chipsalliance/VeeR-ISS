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

using namespace WdRiscv;


template <typename URV>
Core<URV>::Core(URV hartIdBase, unsigned coreIx, unsigned hartsPerCore, Memory& memory)
{
  harts_.resize(hartsPerCore);

  for (unsigned ix = 0; ix < hartsPerCore; ++ix)
    {
      URV hartId = hartIdBase + ix;  // Value in MHARTID of hart.
      unsigned hartIx = coreIx * hartsPerCore + ix;  // Rank of hart in system.
      harts_.at(ix) = std::make_shared<HartClass>(hartIx, hartId, memory);
    }
}


template <typename URV>
Core<URV>::~Core()
{
}


template class WdRiscv::Core<uint32_t>;
template class WdRiscv::Core<uint64_t>;
