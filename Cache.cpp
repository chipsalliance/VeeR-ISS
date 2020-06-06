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
#include "Cache.hpp"

using namespace WdRiscv;


Cache::Cache(uint64_t totalSize, unsigned lineSize, unsigned setSize)
  : size_(totalSize), lineSize_(lineSize), setSize_(setSize)
{
  unsigned logSize = static_cast<unsigned>(std::log2(totalSize));
  uint64_t p2Size = uint64_t(1) << logSize;
  assert(p2Size == totalSize);

  unsigned logSetCount = static_cast<unsigned>(std::log2(setSize));
  unsigned p2SetCount = unsigned(1) << logSetCount;
  assert(p2SetCount == setSize);

  unsigned logLineSize = static_cast<unsigned>(std::log2(lineSize));
  unsigned p2LineSize = unsigned(1) << logLineSize;
  assert(p2LineSize == lineSize);

  lineNumberShift_ = logLineSize;

  assert(totalSize >= lineSize);

  uint64_t lineCount = totalSize / lineSize;
  assert(lineCount >= setSize);

  uint64_t count = lineCount / setSize;
  unsigned logCount = static_cast<unsigned>(std::log2(count));
  uint64_t p2Count = uint64_t(1) << logCount;
  assert(p2Count == count);

  setIndexMask_ = count - 1;

  timesPerSet_.resize(count);
  linesPerSet_.resize(count);

  for (auto& lines : linesPerSet_)
    lines.resize(setSize_);
}


Cache::~Cache()
{
  std::cerr << "Cache access: " << accesses_ << '\n';
  std::cerr << "Cache hits: " << hits_ << '\n';

  double ratio = accesses_ == 0? 0. : double(hits_)/double(accesses_);
  std::cerr << "Hit ratio: " << ratio << '\n';
}
