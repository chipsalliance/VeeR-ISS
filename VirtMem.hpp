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

#include <iosfwd>
#include "trapEnums.hpp"
#include "Memory.hpp"
#include "Tlb.hpp"


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
  } __attribute__((packed));


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

    bool global() const     { return bits_.global_; }

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
  } __attribute((packed));


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

    bool global() const     { return bits_.global_; }

    bool accessed() const   { return bits_.accessed_; }

    bool dirty() const      { return bits_.dirty_; }

    uint64_t ppn() const    { return ppn0() | (ppn1() << 9) | (ppn2() << 18); }

    uint64_t ppn0() const   { return bits_.ppn0_; }

    uint64_t ppn1() const   { return bits_.ppn1_; }

    uint64_t ppn2() const   { return bits_.ppn2_; }

    uint32_t levels() const { return 3; }

    uint32_t size() const   { return sizeof(data_); }

    uint64_t ppn(int i) const
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
  } __attribute__((packed));


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

    bool global() const     { return bits_.global_; }

    bool accessed() const   { return bits_.accessed_; }

    bool dirty() const      { return bits_.dirty_; }

    uint64_t ppn() const    { return ppn0() | (ppn1() << 9) | (ppn2() << 18) | (ppn3() << 27); }

    uint64_t ppn0() const   { return bits_.ppn0_; }

    uint64_t ppn1() const   { return bits_.ppn1_; }

    uint64_t ppn2() const   { return bits_.ppn2_; }

    uint64_t ppn3() const   { return bits_.ppn3_; }

    uint32_t levels() const { return 4; }

    uint32_t size() const   { return sizeof(data_); }

    uint64_t ppn(int i) const
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
  } __attribute__((packed));


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
  } __attribute__((packed));


  /// 39-bit virtual address.
  union Va39
  {
    Va39Bits bits_;
    uint64_t data_ = 0;

    Va39(uint64_t data) : data_(data)
    { }

    uint64_t offset() const { return bits_.offset_; }

    uint64_t vpn0() const   { return bits_.vpn0_; }

    uint64_t vpn1() const   { return bits_.vpn1_; }

    uint64_t vpn2() const   { return bits_.vpn2_; }

    uint64_t vpn(int i) const
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
  } __attribute__((packed));


  /// 48-bit virtual address.
  union Va48
  {
    Va48Bits bits_;
    uint64_t data_ = 0;

    Va48(uint64_t data) : data_(data)
    { }

    uint64_t offset() const { return bits_.offset_; }

    uint64_t vpn0() const   { return bits_.vpn0_; }

    uint64_t vpn1() const   { return bits_.vpn1_; }

    uint64_t vpn2() const   { return bits_.vpn2_; }

    uint64_t vpn3() const   { return bits_.vpn2_; }

    uint64_t vpn(int i) const
    {
      if (i == 0) return vpn0();
      if (i == 1) return vpn1();
      if (i == 2) return vpn2();
      if (i == 3) return vpn3();
      assert(0);
      return 0;
    }
  };


  template <typename URV>
  class Hart;

  class PmpManager;

  class VirtMem
  {
  public:

    friend class Hart<uint32_t>;
    friend class Hart<uint64_t>;

    enum Mode { Bare = 0, Sv32 = 1, Sv39 = 8, Sv48 = 9, Sv57 = 10, Sv64 = 11 };

    VirtMem(unsigned hartIx, Memory& memory, unsigned pageSize,
            PmpManager& pmpMgr, unsigned tlbSize);

    /// Perform virtual to physical memory address translation and
    /// check for read access if the read flag is true (similary also
    /// check for write access is the write flag is true ...).  Return
    /// encoutered exception on failure or ExceptionType::NONE on
    /// success.
    ExceptionCause translate(uint64_t va, PrivilegeMode pm, bool read,
                             bool write, bool exec, uint64_t& pa);

    /// Same as translate but only check for execute access.
    ExceptionCause translateForFetch(uint64_t va, PrivilegeMode pm, uint64_t& pa);

    /// Same as translate but only check for read access.
    ExceptionCause translateForLoad(uint64_t va, PrivilegeMode pm, uint64_t& pa);

    /// Same as translate but only check for write access.
    ExceptionCause translateForStore(uint64_t va, PrivilegeMode pm, uint64_t& pa);

    /// Return page size.
    unsigned pageSize() const
    { return pageSize_; }

    /// Return the address of the first byte in the page containing
    /// the given address.
    uint64_t pageStartAddress(uint64_t address) const
    { return (address >> pageBits_) << pageBits_; }

    /// Debug method: Print all the entries in the page table.
    void printPageTable(std::ostream& os) const;

    /// Print all the page table entries at or below the page table
    /// page rooted at the given address. This is a helper to
    /// printPageTable.
    template <typename PTE, typename VA>
    void printEntries(std::ostream& os, uint64_t addr, std::string path) const;

  protected:

    /// Helper to translate method.
    template <typename PTE, typename VA>
    ExceptionCause pageTableWalk(uint64_t va, PrivilegeMode pm, bool read, bool write,
                                 bool exec, uint64_t& pa, TlbEntry& tlbEntry);

    /// Helper to translate method.
    ExceptionCause pageTableWalkUpdateTlb(uint64_t va, PrivilegeMode pm, bool read,
                                          bool write, bool exec, uint64_t& pa);

    /// Set the page table root page: The root page is placed in
    /// physical memory at address root * page_size
    void setPageTableRootPage(uint64_t root)
    { pageTableRootPage_ = root; }

    // Change the translation mode to m.  Page size is reset to 4096.
    void setMode(Mode m)
    {
      mode_ = m;
      pageSize_ = 4096;
      pageBits_ = 12;
    }

    /// Set the address space id (asid).
    void setAddressSpace(uint32_t asid)
    { asid_ = asid; }

    /// Make executable pages also readable (supports MXR bit in MSTATUS).
    void setExecReadable(bool flag)
    { execReadable_ = flag; }

    /// Allow supervisor-mode code to access user-mode pages (supports SUM
    /// bit in MSTATUS).
    void setSupervisorAccessUser(bool flag)
    { supervisorOk_ = flag; }

    /// Return true if successful and false if page size is not supported.
    bool setPageSize(uint64_t size);

    /// Return current address space id.
    uint32_t addressSpace() const
    { return asid_; }

  private:

    Memory& memory_;
    uint64_t pageTableRootPage_ = 0;
    Mode mode_ = Bare;
    uint32_t asid_ = 0;
    unsigned pageSize_ = 4096;
    unsigned pageBits_ = 12;
    uint64_t pageMask_ = 0xfff;
    unsigned hartIx_ = 0;

    uint64_t time_ = 0;  //  Access order

    // Cached mstatus bits
    bool execReadable_ = false;  // MXR bit
    bool supervisorOk_ = false;  // SUM bit
    bool faultOnFirstAccess_ = true;  // Make this configurable.

    PmpManager& pmpMgr_;
    Tlb tlb_;
  };

}
