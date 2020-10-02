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

  // Temporary, for testing, make all combinations legal. FIX.
  for (auto& groupFlags : legalConfigs_)
    groupFlags.assign(groupFlags.size(), true);

  // Temporary, for testing. FIX.
  sew_ = ElementWidth::DoubleWord;
  sewInBits_ = 64;

  // Temporary, for testing. FIX.
  elems_ = 2;

  std::cerr << "VecRegs::VecRegs: Remove test code\n";
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
VecRegs::config(unsigned bytesPerReg, unsigned bytesPerElem)
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

  if (bytesPerElem < 1)
    {
      std:: cerr << "VecRegd::configure: zero max-bytes-per-element -- using 1\n";
      bytesPerElem = 1;
    }

  unsigned l2BytesPerElem = std::log2(bytesPerElem);
  unsigned p2BytesPerElem = uint32_t(1) << l2BytesPerElem;
  if (p2BytesPerElem != bytesPerElem)
    {
      std::cerr << "VecRegs::configure: max-bytes-per-element (" << bytesPerElem
                << ") not a power of 2 -- using " << p2BytesPerElem << "\n";
      bytesPerElem = p2BytesPerElem;
    }

  if (bytesPerElem > bytesPerReg)
    {
      std::cerr << "VecRegs::configure: max-bytes-per-element (" << bytesPerElem
                << ") is greater than the bytes-per-regidter (" << bytesPerReg
                << " -- using " << p2BytesPerReg << "\n";
      bytesPerElem = bytesPerReg;
    }

  regCount_ = 32;
  bytesPerReg_ = bytesPerReg;
  bytesPerElem_ = bytesPerElem;
  bytesInRegFile_ = regCount_ * bytesPerReg_;

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
  lastElemIx_ = 0;
  lastElemWidth_ = 0;
}
