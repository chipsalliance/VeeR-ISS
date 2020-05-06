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


#pragma once

#include <memory>               // For shapre_ptr
#include <Memory.hpp>


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
    /// assigning to them hart-ids idBase to idBase + n - 1.
    Core(URV idBase, unsigned hartsPerCore, Memory& memory);

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
