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

#include "IntRegs.hpp"


using namespace WdRiscv;


template <typename URV>
IntRegs<URV>::IntRegs(unsigned regCount)
  : regs_(regCount, 0)
{
  numberToName_.resize(32);

  for (unsigned ix = 0; ix < 32; ++ix)
    {
      std::string name = "x" + std::to_string(ix);
      nameToNumber_[name] = IntRegNumber(ix);
      numberToName_[ix] = name;
    }

  numberToAbiName_ = { "zero", "ra", "sp", "gp", "tp", "t0", "t1", "t2",
		       "s0", "s1", "a0", "a1", "a2", "a3", "a4", "a5",
		       "a6", "a7", "s2", "s3", "s4", "s5", "s6", "s7",
		       "s8", "s9", "s10", "s11", "t3", "t4", "t5", "t6" };

  for (unsigned ix = 0; ix < 32; ++ix)
    {
      std::string abiName = numberToAbiName_.at(ix);
      nameToNumber_[abiName] = IntRegNumber(ix);
    }

  nameToNumber_["fp"] = RegX8;   // Fp, s0 and x8 name the same reg.
}


template <typename URV>
bool
IntRegs<URV>::findReg(const std::string& name, unsigned& ix) const
{
  const auto iter = nameToNumber_.find(name);
  if (iter == nameToNumber_.end())
    return false;

  unsigned num = iter->second;
  if (num >= regs_.size())
    return false;

  ix = num;
  return true;
}


template class WdRiscv::IntRegs<uint32_t>;
template class WdRiscv::IntRegs<uint64_t>;
