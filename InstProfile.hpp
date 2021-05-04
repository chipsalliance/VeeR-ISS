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

#include <vector>
#include <unordered_map>
#include "InstId.hpp"


namespace WdRiscv
{

  struct InstProfile
  {
    InstId id_ = InstId::illegal;
    uint64_t freq_ = 0;       // Number of times instruction was executed
    uint64_t user_ = 0;       // Number of times exeuted in user mode.
    uint64_t supervisor_ = 0; // Number of times exeuted in supervisor mode.
    uint64_t machine_ = 0;    // Number of times exeuted in machine mode.

    // One entry per integer register: Count of times register was used
    // by instruction as destination register.
    std::vector<uint64_t> destRegFreq_;

    // SrcRegFreq_[0] corresponds to the 1st source register operand,
    // srcRegFreq_[1] to the second and srcRegFreq_[3] to the
    // 3rd/immediate operand. Each entry of srcRegFreq_[i] is a vector
    // with one entry per integer/fp register indicating the count of
    // times that register was used by the instruction as a source
    // operand.
    std::vector< std::vector<uint64_t> > srcRegFreq_;

    std::vector< std::vector<uint64_t> > srcHisto_;

    bool hasImm_ = false;
    int32_t minImm_ = 0;  // Minimum immediate operand value.
    int32_t maxImm_ = 0;  // Maximum immediate operand value.
  };
}


    
