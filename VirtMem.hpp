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

#include "trapEnums.hpp"
#include "Memory.hpp"


namespace WdRiscv
{

  /// Structure to unpack the fields of a 32-bit page table entry.
  struct Pte32Bits
  {
    unsigned valid_    : 1;
    unsigned read_     : 1;
    unsigned write_    : 1;
    unsigned exec_     : 1;
    unsigned user_     : 1;
    unsigned global_   : 1;
    unsigned accessed_ : 1;
    unsigned dirty_    : 1;
    unsigned rsw_      : 2;   // Reserved for supervisor.
    unsigned ppn0_     : 10;  // Physical page num
    unsigned ppn1_     : 12;  // Physical page num
  };


  /// 32-bit page table entry.
  union Pte32
  {
    Pte32Bits bits_;
    uint32_t data_ = 0;

    Pte32(uint32_t word) : data_(word)
    { }

    bool valid() const      { return bits_.valid_; }

    bool read() const       { return bits_.read_; }

    bool write() const      { return bits_.write_; }

    bool exec() const       { return bits_.exec_; }

    bool user() const       { return bits_.user_; }

    bool accessed() const   { return bits_.accessed_; }

    bool dirty() const      { return bits_.dirty_; }

    uint32_t ppn() const    { return ppn0() | (ppn1() << 10); }

    uint32_t ppn0() const   { return bits_.ppn0_; }

    uint32_t ppn1() const   { return bits_.ppn1_; }

    uint32_t levels() const { return 2; }

    uint32_t size() const   { return sizeof(data_); }

    uint32_t ppn(int i) const
    {
      if (i == 0) return ppn0();
      if (i == 1) return ppn1();
      assert(0); return 0;
    }

    uint32_t paPpnShift(int i) const
    {
      if (i == 0) return 12;
      if (i == 1) return 22;
      assert(0); return 0;
    }
  };


  /// Struct to unpack the fields of a page table entry for Sv39.
  struct Pte39Bits
  {
    unsigned valid_    : 1;
    unsigned read_     : 1;
    unsigned write_    : 1;
    unsigned exec_     : 1;
    unsigned user_     : 1;
    unsigned global_   : 1;
    unsigned accessed_ : 1;
    unsigned dirty_    : 1;
    unsigned rsw_      : 2;   // Reserved for supervisor.
    unsigned ppn0_     : 9;   // Physical page num
    unsigned ppn1_     : 9;   // Physical page num
    unsigned ppn2_     : 26;  // Physical page num
    unsigned res_      : 10;  // Reserved
  };


  /// Page table entry for Sv39
  union Pte39
  {
    Pte39Bits bits_;
    uint64_t data_ = 0;

    Pte39(uint64_t data) : data_(data)
    { }

    bool valid() const      { return bits_.valid_; }

    bool read() const       { return bits_.read_; }

    bool write() const      { return bits_.write_; }

    bool exec() const       { return bits_.exec_; }

    bool user() const       { return bits_.user_; }

    bool accessed() const   { return bits_.accessed_; }

    bool dirty() const      { return bits_.dirty_; }

    uint32_t ppn() const    { return ppn0() | (ppn1() << 9) | (ppn2() << 18); }

    uint32_t ppn0() const   { return bits_.ppn0_; }

    uint32_t ppn1() const   { return bits_.ppn1_; }

    uint32_t ppn2() const   { return bits_.ppn2_; }

    uint32_t levels() const { return 3; }

    uint32_t size() const   { return sizeof(data_); }

    uint32_t ppn(int i) const
    {
      if (i == 0) { return ppn0(); }
      if (i == 1) { return ppn1(); }
      if (i == 2) { return ppn2(); }
      assert(0);
      return 0;
    }

    uint32_t paPpnShift(int i) const
    {
      if (i == 0) { return 12; }
      if (i == 1) { return 21; }
      if (i == 2) { return 30; }
      assert(0);
      return 0;
    }
  };


  struct Pte48Bits
  {
    unsigned valid_    : 1;
    unsigned read_     : 1;
    unsigned write_    : 1;
    unsigned exec_     : 1;
    unsigned user_     : 1;
    unsigned global_   : 1;
    unsigned accessed_ : 1;
    unsigned dirty_    : 1;
    unsigned rsw_      : 2;   // Reserved for supervisor.
    unsigned ppn0_     : 9;   // Physical page num
    unsigned ppn1_     : 9;   // Physical page num
    unsigned ppn2_     : 9;   // Physical page num
    unsigned ppn3_     : 17;  // Physical page num
    unsigned res_      : 10;  // Reserved
  };


  /// Page table entry for Sv48
  union Pte48
  {
    Pte48Bits bits_;
    uint64_t data_ = 0;

    Pte48(uint64_t data) : data_(data)
    { }

    bool valid() const      { return bits_.valid_; }

    bool read() const       { return bits_.read_; }

    bool write() const      { return bits_.write_; }

    bool exec() const       { return bits_.exec_; }

    bool user() const       { return bits_.user_; }

    bool accessed() const   { return bits_.accessed_; }

    bool dirty() const      { return bits_.dirty_; }

    uint32_t ppn() const    { return ppn0() | (ppn1() << 9) | (ppn2() << 18) | (ppn3() << 27); }

    uint32_t ppn0() const   { return bits_.ppn0_; }

    uint32_t ppn1() const   { return bits_.ppn1_; }

    uint32_t ppn2() const   { return bits_.ppn2_; }

    uint32_t ppn3() const   { return bits_.ppn3_; }

    uint32_t levels() const { return 4; }

    uint32_t size() const   { return sizeof(data_); }

    uint32_t ppn(int i) const
    {
      if (i == 0) { return ppn0(); }
      if (i == 1) { return ppn1(); }
      if (i == 2) { return ppn2(); }
      if (i == 3) { return ppn3(); }
      assert(0);
      return 0;
    }

    uint32_t paPpnShift(int i) const
    {
      if (i == 0) { return 12; }
      if (i == 1) { return 21; }
      if (i == 2) { return 30; }
      if (i == 3) { return 39; }
      assert(0);
      return 0;
    }
  };



  /// Structure to unpack the fields of 32-bit virtual address.
  struct Va32Bits
  {
    unsigned offset_ : 12;
    unsigned vpn0_   : 10;
    unsigned vpn1_   : 10;
  };


  /// 32-bit virtual address.
  union Va32
  {
    Va32Bits bits_;
    uint32_t data_ = 0;

    Va32(uint32_t word) : data_(word)
    { }

    uint32_t offset() const { return bits_.offset_; }

    uint32_t vpn0() const   { return bits_.vpn0_; }

    uint32_t vpn1() const   { return bits_.vpn1_; }

    uint32_t vpn(int i) const
    {
      if (i == 0) return vpn0();
      if (i == 1) return vpn1();
      assert(0);
      return 0;
    }
  };


  /// Structure to unpack the fields of Sv39 virtual address.
  struct Va39Bits
  {
    unsigned offset_ : 12;
    unsigned vpn0_   : 9;
    unsigned vpn1_   : 9;
    unsigned vpn2_   : 9;
  };


  /// 39-bit virtual address.
  union Va39
  {
    Va39Bits bits_;
    uint64_t data_ = 0;

    Va39(uint64_t data) : data_(data)
    { }

    uint32_t offset() const { return bits_.offset_; }

    uint32_t vpn0() const   { return bits_.vpn0_; }

    uint32_t vpn1() const   { return bits_.vpn1_; }

    uint32_t vpn2() const   { return bits_.vpn2_; }

    uint32_t vpn(int i) const
    {
      if (i == 0) return vpn0();
      if (i == 1) return vpn1();
      if (i == 2) return vpn2();
      assert(0);
      return 0;
    }
  };


  /// Structure to unpack the fields of Sv48 virtual address.
  struct Va48Bits
  {
    unsigned offset_ : 12;
    unsigned vpn0_   : 9;
    unsigned vpn1_   : 9;
    unsigned vpn2_   : 9;
    unsigned vpn3_   : 9;
  };


  /// 48-bit virtual address.
  union Va48
  {
    Va48Bits bits_;
    uint64_t data_ = 0;

    Va48(uint64_t data) : data_(data)
    { }

    uint32_t offset() const { return bits_.offset_; }

    uint32_t vpn0() const   { return bits_.vpn0_; }

    uint32_t vpn1() const   { return bits_.vpn1_; }

    uint32_t vpn2() const   { return bits_.vpn2_; }

    uint32_t vpn3() const   { return bits_.vpn2_; }

    uint32_t vpn(int i) const
    {
      if (i == 0) return vpn0();
      if (i == 1) return vpn1();
      if (i == 2) return vpn2();
      if (i == 3) return vpn3();
      assert(0);
      return 0;
    }
  };


  class VirtMem
  {
  public:

    enum Mode { Bare, Sv32, Sv39, Sv48, Sv57, Sv64 };

    VirtMem(Memory& memory, unsigned pageSize);

    /// Perform virtual to physical memory address translation.
    /// Return encoutered exception on failure or ExceptionType::NONE
    /// on success.
    ExceptionCause translate(size_t va, PrivilegeMode pm, bool read,
                             bool write, bool exec, size_t& pa);

    template <typename PTE, typename VA>
    ExceptionCause translate_(size_t va, PrivilegeMode pm, bool read,
                              bool write, bool exec, size_t& pa);

  protected:
    void setPageTableRoot(uint64_t root)
    { pageTableRoot_ = root; }

  private:

    Memory& memory_;
    uint64_t pageTableRoot_ = 0;
    Mode mode_ = Bare;
    unsigned pageSize_ = 4096;
    unsigned pageBits_ = 12;

    // Cached mstatus bits
    bool execReadable_ = false;  // MXR bit
    bool supervisorOk_ = false;  // SUM bit
  };

}

