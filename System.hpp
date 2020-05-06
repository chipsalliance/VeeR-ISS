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

#include <Memory.hpp>


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
    /// are indexed with 0 to n*m - 1.
    System(unsigned coreCount, unsigned hartsPerCore, Memory& memory);

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
    /// out of bounds.
    std::shared_ptr<HartClass> ithHart(unsigned i)
    {
      unsigned coreIx = i / hartsPerCore_;
      if (coreIx >= cores_.size())
	return std::shared_ptr<HartClass>();
      unsigned hartInCore = i % hartsPerCore_;
      return cores_.at(coreIx)->ithHart(hartInCore);
    }

  private:

    unsigned hartCount_;
    unsigned hartsPerCore_;

    std::vector< std::shared_ptr<CoreClass> > cores_;
  };
}
