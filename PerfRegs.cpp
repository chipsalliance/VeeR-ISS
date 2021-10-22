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
#include "PerfRegs.hpp"


using namespace WdRiscv;


PerfRegs::PerfRegs(unsigned numCounters)
{
  // 29 counters: MHPMCOUNTER3 to MHPMCOUNTER31
  counters_.resize(29);

  config(numCounters);
}


void
PerfRegs::config(unsigned numCounters)
{
  assert(numCounters < counters_.size());

  eventOfCounter_.resize(numCounters);
  enableUser_.resize(numCounters);
  enableMachine_.resize(numCounters);
}


bool
PerfRegs::applyPerfEventAssign()
{
  if (not hasPending_)
    return false;

  hasPending_ = false;

  if (pendingCounter_ >= eventOfCounter_.size())
    return false;

  eventOfCounter_.at(pendingCounter_) = pendingEvent_;
  enableUser_.at(pendingCounter_) = pendingUser_;
  enableMachine_.at(pendingCounter_) = pendingMachine_;

  return true;
}


void
PerfRegs::reset()
{
  eventOfCounter_.assign(eventOfCounter_.size(), EventNumber::None);
}
