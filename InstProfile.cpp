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

#include <algorithm>
#include "InstProfile.hpp"

using namespace WdRiscv;


void
InstProfiles::configure()
{
  unsigned instCount = unsigned(InstId::maxId) + 1;

  // For vector instructions, we track each combination of
  // instruction id and vector element size (for size = 8, 16, 32,
  // and 64). To do so we create 4 entries per instruction (one for
  // each size). First copy is used for regular instructions and
  // for vector element size 8. Remaining are for vec elem sizes
  // 16, 32, and 64.
  unsigned sizeCount = unsigned(ElementWidth::Word2) + 1;
  vec_.resize(instCount*sizeCount);

  unsigned regCount = 32;

  for (unsigned sizeIx = 0; sizeIx < sizeCount; ++sizeIx)
    {
      ElementWidth ew{sizeIx};
      for (unsigned instIx = 0; instIx < instCount; ++instIx)
	{
	  InstId id{instIx};
	  size_t vecIx = instCount*sizeIx + instIx;

	  auto& inst = vec_.at(vecIx);

	  inst.id_ = id;
	  inst.elemWidth_ = ew;
	  inst.destRegFreq_.resize(regCount);
	  inst.srcRegFreq_.resize(3);  // Up to 3 source operands
	  for (auto& vec : inst.srcRegFreq_)
	    vec.resize(regCount);

	  inst.srcHisto_.resize(3);  // Up to 3 source historgrams
	  for (auto& vec : inst.srcHisto_)
	    vec.resize(13);  // FIX: avoid magic 13
	}
    }
}


void
InstProfiles::sort(std::vector<size_t>& indices) const
{
  struct CompareFreq
  {
    CompareFreq(const std::vector<InstProfile>& profileVec)
      : profileVec(profileVec)
    { }

    bool operator()(size_t a, size_t b) const
    { return profileVec.at(a).freq_ < profileVec.at(b).freq_; }

    const std::vector<InstProfile>& profileVec;
  };

  indices.clear();
  indices.resize(vec_.size());
  for (size_t i = 0; i < indices.size(); ++i)
    indices.at(i) = i;
  std::sort(indices.begin(), indices.end(), CompareFreq(vec_));
}
