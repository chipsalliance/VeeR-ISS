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

#include <iostream>
#include <cmath>
#include <string.h>
#include "VecRegs.hpp"


using namespace WdRiscv;


VecRegs::VecRegs()
{
  // Intialize structure (vector of vectors) defining legal
  // element-width/group-multiplier combinations to all false: No
  // combination is supported. Configuration code will later change
  // this.
  legalConfigs_.resize(size_t(VecEnums::WidthLimit));
  for (auto& groupFlags : legalConfigs_)
    groupFlags.resize(size_t(VecEnums::GroupLimit));

  // Start by making all combinations of width/grouping legal. This
  // gets adjusted when config method is called.
  for (auto& groupFlags : legalConfigs_)
    groupFlags.assign(groupFlags.size(), true);

  // Operands effective group multiplier is used for tracing/logging.
  // Worst case we have 3 vector operands.
  opsEmul_.resize(3);
  opsEmul_.assign(opsEmul_.size(), 1);
}


VecRegs::~VecRegs()
{
  regCount_ = 0;
  bytesPerReg_ = 0;
  bytesInRegFile_ = 0;
  delete [] data_;
  data_ = nullptr;
}


void
VecRegs::config(unsigned bytesPerReg, unsigned minBytesPerElem,
		unsigned maxBytesPerElem)
{
  if (bytesPerReg > 4096)
    {
      std::cerr << "VecRegs::configure: bytes-per-register too large (" << bytesPerReg
                << ") -- using 4096\n";
      bytesPerReg = 4096;
    }

  if (bytesPerReg <= 4)
    {
      std::cerr << "VecRegs::configure: bytes-per-register too small (" << bytesPerReg
                << ") -- using 4\n";
      bytesPerReg = 4;
    }

  unsigned l2BytesPerReg = std::log2(bytesPerReg);
  unsigned p2BytesPerReg = uint32_t(1) << l2BytesPerReg;
  if (p2BytesPerReg != bytesPerReg)
    {
      std::cerr << "VecRegs::configure: bytes-per-register (" << bytesPerReg
                << ") not a power of 2 -- using " << p2BytesPerReg << "\n";
      bytesPerReg = p2BytesPerReg;
    }

  if (minBytesPerElem < 1)
    {
      std:: cerr << "VecRegd::configure: zero min-bytes-per-element -- using 1\n";
      minBytesPerElem = 1;
    }

  if (maxBytesPerElem < 1)
    {
      std:: cerr << "VecRegd::configure: zero max-bytes-per-element -- using 1\n";
      maxBytesPerElem = 1;
    }

  if (minBytesPerElem > maxBytesPerElem)
    {
      std:: cerr << "VecRegd::configure: min-bytes-per-elem larger than max -- using max\n";
      minBytesPerElem = maxBytesPerElem;
    }

  unsigned l2BytesPerElem = std::log2(maxBytesPerElem);
  unsigned p2BytesPerElem = uint32_t(1) << l2BytesPerElem;
  if (p2BytesPerElem != maxBytesPerElem)
    {
      std::cerr << "VecRegs::configure: max-bytes-per-element (" << maxBytesPerElem
                << ") not a power of 2 -- using " << p2BytesPerElem << "\n";
      maxBytesPerElem = p2BytesPerElem;
    }

  if (maxBytesPerElem > bytesPerReg)
    {
      std::cerr << "VecRegs::configure: max-bytes-per-element (" << maxBytesPerElem
                << ") is greater than the bytes-per-register (" << bytesPerReg
                << " -- using " << bytesPerReg << "\n";
      maxBytesPerElem = bytesPerReg;
    }

  l2BytesPerElem = std::log2(minBytesPerElem);
  p2BytesPerElem = uint32_t(1) << l2BytesPerElem;
  if (p2BytesPerElem != minBytesPerElem)
    {
      std::cerr << "VecRegs::configure: min-bytes-per-element (" << minBytesPerElem
                << ") not a power of 2 -- using " << p2BytesPerElem << "\n";
      minBytesPerElem = p2BytesPerElem;
    }

  if (minBytesPerElem > bytesPerReg)
    {
      std::cerr << "VecRegs::configure: min-bytes-per-element (" << minBytesPerElem
                << ") is greater than the bytes-per-register (" << bytesPerReg
                << " -- using " << bytesPerReg << "\n";
      minBytesPerElem = bytesPerReg;
    }

  regCount_ = 32;
  bytesPerReg_ = bytesPerReg;
  minBytesPerElem_ = minBytesPerElem;
  maxBytesPerElem_ = maxBytesPerElem;
  bytesInRegFile_ = regCount_ * bytesPerReg_;

  // Make illegal all group entries for element-widths greater than
  // the max-element-width (which is in bytesPerElem_) or smaller
  // than the min-element-width.
  for (unsigned i = 0; i <= unsigned(ElementWidth::Word32); ++i)
    {
      ElementWidth ew = ElementWidth(i);
      unsigned bytes = VecRegs::elementWidthInBytes(ew);
      if (bytes > maxBytesPerElem_ or bytes < minBytesPerElem_ )
	{
	  auto& groupFlags = legalConfigs_.at(size_t(ew));
	  groupFlags.assign(groupFlags.size(), false);
	}
    }

  delete [] data_;
  data_ = new uint8_t[bytesInRegFile_];
  memset(data_, 0, bytesInRegFile_);
}


void
VecRegs::reset()
{
  if (data_)
    memset(data_, 0, bytesInRegFile_);
  lastWrittenReg_ = -1;
  lastGroupX8_ = 8;
}
