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
System<URV>::System(unsigned coreCount, unsigned hartsPerCore, Memory& memory)
  : hartCount_(coreCount * hartsPerCore), hartsPerCore_(hartsPerCore)
{
  cores_.resize(coreCount);

  for (unsigned ix = 0; ix < coreCount; ++ix)
    {
      URV hartIdBase = ix * hartsPerCore;
      cores_.at(ix) = std::make_shared<CoreClass>(hartIdBase, hartsPerCore,
						  memory);
    }
}


template <typename URV>
System<URV>::~System()
{
}


template class WdRiscv::System<uint32_t>;
template class WdRiscv::System<uint64_t>;
