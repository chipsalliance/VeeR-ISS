//
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2018 Western Digital Corporation or its affiliates.
// 
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
// 
// You should have received a copy of the GNU General Public License along with
// this program. If not, see <https://www.gnu.org/licenses/>.
//

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
