#include <cmath>
#include "VirtMem.hpp"

using namespace WdRiscv;


VirtMem::VirtMem(unsigned hartIx, Memory& memory, unsigned pageSize)
  : memory_(memory), mode_(Sv32), pageSize_(pageSize), hartIx_(hartIx)
{
  pageBits_ = static_cast<unsigned>(std::log2(pageSize_));
  unsigned p2PageSize =  unsigned(1) << pageBits_;

  assert(p2PageSize == pageSize);
  assert(pageSize >= 64);
}


ExceptionCause
VirtMem::translate(size_t va, PrivilegeMode pm, bool read, bool write,
                   bool exec, size_t& pa)
{
  if (mode_ == Bare)
    {
      pa = va;
      return ExceptionCause::NONE;
    }

  if (mode_ == Sv32)
    return translate_<Pte32, Va32>(va, pm, read, write, exec, pa);

  if (mode_ == Sv39)
    return translate_<Pte39, Va39>(va, pm, read, write, exec, pa);

  if (mode_ == Sv48)
    return translate_<Pte48, Va48>(va, pm, read, write, exec, pa);

  assert(0 and "Unspupported virtual memory mode.");
  return ExceptionCause::LOAD_PAGE_FAULT;
}


inline
ExceptionCause
pageFaultType(bool read, bool write, bool exec)
{
  if (exec)  return ExceptionCause::INST_PAGE_FAULT;
  if (read)  return ExceptionCause::LOAD_PAGE_FAULT;
  if (write) return ExceptionCause::STORE_PAGE_FAULT;
  assert(0);
  return ExceptionCause::STORE_PAGE_FAULT;
}


template<typename PTE, typename VA>
ExceptionCause
VirtMem::translate_(size_t address, PrivilegeMode privMode,
                    bool read, bool write, bool exec,
                    size_t& pa)
{
  // TBD check xlen against valen.

  PTE pte(0);

  const unsigned levels = pte.levels();
  const unsigned pteSize = pte.size();

  VA va(address);

  uint64_t root = pageTableRoot_;
  uint64_t pteAddr = 0;

  // 2.
  int ii = levels - 1;

  while (true)
    {
      // 3.

      uint32_t vpn = va.vpn(ii);
      uint64_t pteAddr = root + vpn*pteSize;

      // TBD: Check pmp
      if (! memory_.read(pteAddr, pte.data_))
        return pageFaultType(read, write, exec);

      // 4.
      if (not pte.valid() or (not pte.read() and pte.write()))
        return pageFaultType(read, write, exec);

      // 5.
      if (not pte.read() and not pte.exec())
        {
          ii = ii - 1;
          if (ii < 0)
            return pageFaultType(read, write, exec);
          root = pte.ppn();
          // goto 3.
        }
      else
        break;  // goto 6.
    }

  // 6.  pte.read_ or pte.exec_ : leaf pte
  if (privMode == PrivilegeMode::User and not pte.user())
    return pageFaultType(read, write, exec);
  if (privMode == PrivilegeMode::Supervisor and pte.user() and
      not supervisorOk_)
    return pageFaultType(read, write, exec);

  bool pteRead = pte.read() or (execReadable_ and pte.exec());
  if ((read and not pteRead) or (write and not pte.write()) or
      (exec and pte.exec()))
        return pageFaultType(read, write, exec);

  // 7.
  if (ii > 0 and pte.ppn1() != 0 and pte.ppn0() != 0)
    return pageFaultType(read, write, exec);

  // 8.
  if (not pte.accessed() or (write and not pte.dirty()))
    {
      // We have a choice:
      // A. Page fault
      if (faultOnFirstAccess_)
        return pageFaultType(read, write, exec);  // A

      // Or B
      // B1. Set pte->accessed_ to 1 and, if a write, set pte->dirty_ to 1.
      // B2. Access fault if PMP violation.
      pte.bits_.accessed_ = 1;
      if (write)
        pte.bits_.dirty_ = 1;

      // TBD: Check pmp
      if (not memory_.write(hartIx_, pteAddr, pte.data_))
        return pageFaultType(read, write, exec);
    }

  // 9.
  pa = va.offset();
  if (ii > 0)
    for (int i = ii - 1; i >= 0; --i)
      pa = pa | (va.vpn(i) << pte.paPpnShift(i));

  for (int i = levels - 1; i >= ii; --i)
    pa = pa | pte.ppn(i) << pte.paPpnShift(i);

  return ExceptionCause::NONE;
}
