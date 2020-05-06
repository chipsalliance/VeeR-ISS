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

using namespace WdRiscv;


template <typename URV>
Core<URV>::Core(URV idBase, unsigned hartsPerCore, Memory& memory)
{
  harts_.resize(hartsPerCore);

  for (unsigned ix = 0; ix < hartsPerCore; ++ix)
    {
      URV hartId = idBase + ix;
      harts_.at(ix) = std::make_shared<HartClass>(hartId, memory);
    }
}


template <typename URV>
Core<URV>::~Core()
{
}


template class WdRiscv::Core<uint32_t>;
template class WdRiscv::Core<uint64_t>;
