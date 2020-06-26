#include <cmath>
#include "VirtMem.hpp"

using namespace WdRiscv;


Tlb::Tlb(unsigned size)
  : entries_(size)
{
}


void
Tlb::insertEntry(uint64_t virtPageNum, uint64_t physPageNum, uint32_t asid,
                 bool global, bool isUser, bool read, bool write, bool exec)
{
  TlbEntry* best = nullptr;
  for (size_t i = 0; i < entries_.size(); ++i)
    {
      auto& entry = entries_[i];
      if (not entry.valid_ or (best and entry.time_ < best->time_))
        best = &entry;
    }
  best->valid_ = true;
  best->virtPageNum_ = virtPageNum;
  best->physPageNum_ = physPageNum;
  best->time_ = time_++;
  best->asid_ = asid;
  best->global_ = global;
  best->user_ = isUser;
  best->read_ = read;
  best->write_ = write;
  best->exec_ = exec;
}


void
Tlb::insertEntry(const TlbEntry& te)
{
  TlbEntry* best = nullptr;
  for (size_t i = 0; i < entries_.size(); ++i)
    {
      auto& entry = entries_[i];
      if (not entry.valid_)
        {
          best = &entry;
          break;
        }
      if (not best or entry.time_ < best->time_)
        best = &entry;
    }

  *best = te;
  best->time_ = time_++;
}
