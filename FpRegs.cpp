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

#include "FpRegs.hpp"


using namespace WdRiscv;


FpRegs::FpRegs(unsigned regCount)
  : regs_(regCount, 0)
{
  numberToName_.resize(32);

  for (unsigned ix = 0; ix < 32; ++ix)
    {
      std::string name = "f" + std::to_string(ix);
      nameToNumber_[name] = FpRegNumber(ix);
      numberToName_[ix] = name;
    }

  numberToAbiName_ = { "ft0", "ft1", "ft2", "ft3", "ft4", "ft5", "ft6", "ft7",
		       "fs0", "fs1", "fa0", "fa1", "fa2", "fa3", "fa4", "fa5",
		       "fa6", "fa7", "fs2", "fs3", "fs4", "fs5", "fs6", "fs7",
		       "fs8", "fs9", "fs10", "fs11", "ft8", "ft9", "ft0", "ft11" };

  for (unsigned ix = 0; ix < 32; ++ix)
    {
      std::string abiName = numberToAbiName_.at(ix);
      nameToNumber_[abiName] = FpRegNumber(ix);
    }
}


bool
FpRegs::findReg(const std::string& name, unsigned& ix) const
{
  const auto iter = nameToNumber_.find(name);
  if (iter == nameToNumber_.end())
    return false;

  ix = iter->second;
  return true;
}
