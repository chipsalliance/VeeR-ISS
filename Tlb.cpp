#include <cmath>
#include "VirtMem.hpp"
#include <bitset>
#include <sstream>
#include <iomanip>

using namespace WdRiscv;


Tlb::Tlb(unsigned size)
  : useFullTlb_(size>1024), entries_((size<=1024)?size:0)
{
}


void
Tlb::insertEntry(uint64_t virtPageNum, uint64_t physPageNum, uint32_t asid,
                 bool global, bool isUser, bool read, bool write, bool exec)
{
  TlbEntry* best = nullptr;
  if(useFullTlb_) {
	  best = &fullTlb_[global].insert(std::make_pair(((virtPageNum<<16) | (asid*global)), TlbEntry())).first->second;
  }
  else {
	  for (size_t i = 0; i < entries_.size(); ++i)
		{
		  auto& entry = entries_[i];
		  if (not entry.valid_ or (best and entry.time_ < best->time_))
			best = &entry;
		}
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

void Tlb::printTlb(std::ostream& ost) const {
	for(auto& glb: {true,false})
		for(auto& it: fullTlb_[glb]) printEntry(ost, it.second);
	for(auto&te: entries_) printEntry(ost, te);
}

void Tlb::printEntry(std::ostream& ost, const TlbEntry& te) const {
	if(not te.valid_) return;
	ost << std::hex << std::setfill('0') << std::setw(10) << te.virtPageNum_ << ",";
	if(te.global_) ost << "***";
	else ost << std::setfill('0') << std::setw(4) << te.asid_;
	ost << " -> " << std::setfill('0') << std::setw(10) << te.physPageNum_;
	ost << " P:" << (te.read_?"r":"-") << (te.write_?"w":"-") << (te.exec_?"x":"-");
	ost << " A:" << (te.accessed_ ? "a" : "-") << (te.dirty_?"d":"-");
	ost << " S:" << (te.levels_==3?"1G" : te.levels_==2?"2M" : "4K") << "\n";
}

void
Tlb::insertEntry(const TlbEntry& te)
{
  if(useFullTlb_) {
	  fullTlb_[te.global_].insert(std::make_pair(((te.virtPageNum_<<16) | (te.asid_*te.global_)), TlbEntry()));
	  return;
  }
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
