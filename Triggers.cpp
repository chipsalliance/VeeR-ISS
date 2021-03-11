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

#include "Triggers.hpp"


using namespace WdRiscv;


template <typename URV>
Triggers<URV>::Triggers(unsigned count)
  : triggers_(count)
{
  // Define each triggers as a single-element chain.
  for (unsigned i = 0; i < count; ++i)
    triggers_.at(i).setChainBounds(i, i+1);
}


template <typename URV>
bool
Triggers<URV>::readData1(URV trigger, URV& value) const
{
  if (trigger >= triggers_.size())
    return false;

  value = triggers_.at(trigger).readData1();
  return true;
}


template <typename URV>
bool
Triggers<URV>::readData2(URV trigger, URV& value) const
{
  if (trigger >= triggers_.size())
    return false;

  value = triggers_.at(trigger).readData2();
  return true;
}


template <typename URV>
bool
Triggers<URV>::readData3(URV trigger, URV& value) const
{
  if (trigger >= triggers_.size())
    return false;

  value = triggers_.at(trigger).readData3();
  return true;
}


template <typename URV>
bool
Triggers<URV>::writeData1(URV trigIx, bool debugMode, URV value)
{
  if (trigIx >= triggers_.size())
    return false;

  auto& trig = triggers_.at(trigIx);

  Data1Bits d1bits(value); // Unpack value of attempted write.

  // If next trigger has a dmode of 1 (debugger only), then clear
  // chain bit in attempted wirte if write would also set dmode to 0.
  // Otherwise we would have a chain with different dmodes.
  if (trigIx + 1 < triggers_.size())
    {
      auto& nextTrig = triggers_.at(trigIx + 1);
      if (nextTrig.isDebugModeOnly() and not d1bits.dmodeOnly())
        {
          d1bits.mcontrol_.chain_ = 0;
          value = d1bits.value_;
        }
    }

  // Write is ignored if it would set dmode and previous trigger has
  // both dmode=0 and chain=1. Otherwise, we would have a chain with
  // different dmodes.
  if (d1bits.dmodeOnly() and trigIx > 0)
    {
      auto& prevTrig = triggers_.at(trigIx - 1);
      if (prevTrig.getChain() and not prevTrig.isDebugModeOnly())
        return false;
    }

  bool oldChain = trig.getChain();

  if (not trig.writeData1(debugMode, value))
    return false;

  bool newChain = trig.getChain();
  if (oldChain != newChain)
    defineChainBounds();

  return true;
}


template <typename URV>
bool
Triggers<URV>::writeData2(URV trigger, bool debugMode, URV value)
{
  if (trigger >= triggers_.size())
    return false;

  return triggers_.at(trigger).writeData2(debugMode, value);
}


template <typename URV>
bool
Triggers<URV>::writeData3(URV trigger, bool /*debugMode*/, URV /*value*/)
{
  if (trigger >= triggers_.size())
    return false;

  return false;
}


template <typename URV>
bool
Triggers<URV>::updateChainHitBit(Trigger<URV>& trigger)
{
  bool chainHit = true;  // True if all items in chain hit
  TriggerTiming  timing = trigger.getTiming();
  bool uniformTiming = true;

  size_t beginChain = 0, endChain = 0;
  trigger.getChainBounds(beginChain, endChain);

  for (size_t i = beginChain; i < endChain; ++i)
    {
      auto& trig = triggers_.at(i);
      chainHit = chainHit and trig.getLocalHit();
      uniformTiming = uniformTiming and (timing == trig.getTiming());
      if (chainHit)
        trig.setHit(true);
    }

  if (not chainHit or not uniformTiming)
    return false;

  for (size_t i = beginChain; i < endChain; ++i)
    triggers_.at(i).setTripped(true);

  return true;
}


template <typename URV>
bool
Triggers<URV>::ldStAddrTriggerHit(URV address, TriggerTiming timing,
				  bool isLoad, PrivilegeMode mode,
                                  bool interruptEnabled)
{
  // Check if we should skip tripping because we are running in
  // machine mode and interrupts disabled.
  bool skip = mode == PrivilegeMode::Machine and not interruptEnabled;

  bool chainHit = false;  // Chain hit.
  for (auto& trigger : triggers_)
    {
      if (not trigger.isEnterDebugOnHit() and skip)
	continue;

      if (not trigger.matchLdStAddr(address, timing, isLoad, mode))
	continue;

      trigger.setLocalHit(true);

      if (updateChainHitBit(trigger))
	chainHit = true;
    }
  return chainHit;
}


template <typename URV>
bool
Triggers<URV>::ldStDataTriggerHit(URV value, TriggerTiming timing, bool isLoad,
				  PrivilegeMode mode, bool interruptEnabled)
{
  // Check if we should skip tripping because we are running in
  // machine mode and interrupts disabled.
  bool skip = mode == PrivilegeMode::Machine and not interruptEnabled;

  bool chainHit = false;  // Chain hit.
  for (auto& trigger : triggers_)
    {
      if (not trigger.isEnterDebugOnHit() and skip)
	continue;

      if (not trigger.matchLdStData(value, timing, isLoad, mode))
	continue;

      trigger.setLocalHit(true);

      if (updateChainHitBit(trigger))
	chainHit = true;
    }

  return chainHit;
}


template <typename URV>
bool
Triggers<URV>::instAddrTriggerHit(URV address, TriggerTiming timing,
                                  PrivilegeMode mode, bool interruptEnabled)
{
  // Check if we should skip tripping because we are running in
  // machine mode and interrupts disabled.
  bool skip = mode == PrivilegeMode::Machine and not interruptEnabled;

  bool chainHit = false;  // Chain hit.
  for (auto& trigger : triggers_)
    {
      if (not trigger.isEnterDebugOnHit() and skip)
	continue;

      if (not trigger.matchInstAddr(address, timing, mode))
	continue;

      trigger.setLocalHit(true);

      if (updateChainHitBit(trigger))
	chainHit = true;
    }

  return chainHit;
}


template <typename URV>
bool
Triggers<URV>::instOpcodeTriggerHit(URV opcode, TriggerTiming timing,
                                    PrivilegeMode mode,
				    bool interruptEnabled)
{
  // Check if we should skip tripping because we are running in
  // machine mode and interrupts disabled.
  bool skip = mode == PrivilegeMode::Machine and not interruptEnabled;

  bool hit = false;
  for (auto& trigger : triggers_)
    {
      if (not trigger.isEnterDebugOnHit() and skip)
	continue;

      if (not trigger.matchInstOpcode(opcode, timing, mode))
	continue;

      trigger.setLocalHit(true);

      if (updateChainHitBit(trigger))
	hit = true;
    }

  return hit;
}


template <typename URV>
bool
Triggers<URV>::icountTriggerHit(PrivilegeMode mode, bool interruptEnabled)
{
  // Check if we should skip tripping because we are running in
  // machine mode and interrupts disabled.
  bool skip = mode == PrivilegeMode::Machine and not interruptEnabled;

  bool hit = false;

  for (auto& trig : triggers_)
    {
      if (not trig.isEnterDebugOnHit() and skip)
	continue;

      if (trig.isModified())
	continue; // Trigger was written by current instruction.

      if (not trig.instCountdown(mode))
	continue;

      hit = true;
      trig.setHit(true);
      trig.setLocalHit(true);
    }
  return hit;
}


template <typename URV>
bool
Triggers<URV>::config(unsigned trigger,
                      uint64_t rv1, uint64_t rv2, uint64_t rv3,
		      uint64_t wm1, uint64_t wm2, uint64_t wm3,
		      uint64_t pm1, uint64_t pm2, uint64_t pm3)
{
  if (trigger <= triggers_.size())
    triggers_.resize(trigger + 1);

  triggers_.at(trigger).configData1(rv1, wm1, pm1);
  triggers_.at(trigger).configData2(rv2, wm2, pm2);
  triggers_.at(trigger).configData3(rv3, wm3, pm3);

  triggers_.at(trigger).writeData2(true, rv2);  // Define compare mask.

  defineChainBounds();

  return true;
}


template <typename URV>
void
Triggers<URV>::reset()
{
  for (auto& trigger : triggers_)
    trigger.reset();
  defineChainBounds();
}



template <typename URV>
bool
Triggers<URV>::peek(unsigned trigger, uint64_t& data1, uint64_t& data2,
                    uint64_t& data3) const
{
  if (trigger >= triggers_.size())
    return false;

  return triggers_.at(trigger).peek(data1, data2, data3);
}


template <typename URV>
bool
Triggers<URV>::peek(unsigned trigger,
                    uint64_t& data1, uint64_t& data2, uint64_t& data3,
		    uint64_t& wm1, uint64_t& wm2, uint64_t& wm3,
		    uint64_t& pm1, uint64_t& pm2, uint64_t& pm3) const
{
  if (trigger >= triggers_.size())
    return false;

  const Trigger<URV>& trig = triggers_.at(trigger);
  return trig.peek(data1, data2, data3, wm1, wm2, wm3, pm1, pm2, pm3);
}


template <typename URV>
bool
Triggers<URV>::poke(URV trigger, URV v1, URV v2, URV v3)
{
  if (trigger >= triggers_.size())
    return false;

  Trigger<URV>& trig = triggers_.at(trigger);

  trig.pokeData1(v1);
  trig.pokeData2(v2);
  trig.pokeData3(v3);

  return true;
}


template <typename URV>
bool
Triggers<URV>::pokeData1(URV trigIx, URV value)
{
  if (trigIx >= triggers_.size())
    return false;

  auto& trig = triggers_.at(trigIx);

  Data1Bits d1bits(value); // Unpack value of attempted write.

  // If next trigger has a dmode of 1 (debugger only), then clear
  // chain bit in attempted wirte if write would also set dmode to 0.
  // Otherwise we would have a chain with different dmodes.
  if (trigIx + 1 < triggers_.size())
    {
      auto& nextTrig = triggers_.at(trigIx + 1);
      if (nextTrig.isDebugModeOnly() and not d1bits.dmodeOnly())
        {
          d1bits.mcontrol_.chain_ = 0;
          value = d1bits.value_;
        }
    }

  // Write is ignored if it would set dmode and previous trigger has
  // both dmode=0 and chain=1. Otherwise, we would have a chain with
  // different dmodes.
  if (d1bits.dmodeOnly() and trigIx > 0)
    {
      auto& prevTrig = triggers_.at(trigIx - 1);
      if (prevTrig.getChain() and not prevTrig.isDebugModeOnly())
        return false;
    }

  bool oldChain = trig.getChain();

  trig.pokeData1(value);

  bool newChain = trig.getChain();
  if (oldChain != newChain)
    defineChainBounds();

  return true;
}


template <typename URV>
bool
Triggers<URV>::pokeData2(URV trigger, URV val)
{
  if (trigger >= triggers_.size())
    return false;

  Trigger<URV>& trig = triggers_.at(trigger);

  trig.pokeData2(val);
  return true;
}


template <typename URV>
bool
Triggers<URV>::pokeData3(URV trigger, URV val)
{
  if (trigger >= triggers_.size())
    return false;

  Trigger<URV>& trig = triggers_.at(trigger);

  trig.pokeData3(val);
  return true;
}


template <typename URV>
void
Triggers<URV>::defineChainBounds()
{
  if (chainPairs_)
    {
      // Reset each trigger to a chain of length 1.
      for  (size_t i = 0; i < triggers_.size(); ++i)
	triggers_.at(i).setChainBounds(i, i+1);

      // Only chain consecutive even/odd pairs if chain bit of even is set.
      for (size_t i = 0; i < triggers_.size(); i += 2)
	{
	  if (triggers_.at(i).getChain() and i + 1 < triggers_.size())
	    {
	      triggers_.at(i).setChainBounds(i, i + 2);
	      triggers_.at(i+1).setChainBounds(i, i + 2);
	    }
	}
      return;
    }

  size_t begin = 0, end = 0;

  for (size_t i = 0; i < triggers_.size(); ++i)
    {
      if (not triggers_.at(i).getChain())
	{
	  end = i + 1;
	  for (size_t j = begin; j < end; j++)
	    triggers_.at(j).setChainBounds(begin, end);
	  begin = end;
	}
    }

  end = triggers_.size();
  for  (size_t i = begin; i < end; ++i)
    triggers_.at(i).setChainBounds(begin, end);
}


template <typename URV>
bool
Trigger<URV>::matchLdStAddr(URV address, TriggerTiming timing, bool isLoad,
                            PrivilegeMode mode) const
{
  if (not data1_.isAddrData())
    return false;  // Not an address trigger.

  const Mcontrol<URV>& ctl = data1_.mcontrol_;

  if (mode == PrivilegeMode::Machine and not ctl.m_)
    return false;  // Not enabled;

  if (mode == PrivilegeMode::Supervisor and not ctl.s_)
    return false;  // Not enabled;

  if (mode == PrivilegeMode::User and not ctl.u_)
    return false;  // Not enabled;

  if (mode == PrivilegeMode::Reserved)
    return false;  // Not enabled;

  bool isStore = not isLoad;
  bool clearBit0 = false;

  if (TriggerTiming(ctl.timing_) == timing and
      Select(ctl.select_) == Select::MatchAddress and
      ((isLoad and ctl.load_) or (isStore and ctl.store_)))
    return doMatch(address, clearBit0);

  return false;
}


template <typename URV>
bool
Trigger<URV>::matchLdStData(URV value, TriggerTiming timing, bool isLoad,
                            PrivilegeMode mode) const
{
  if (not data1_.isAddrData())
    return false;  // Not an address trigger.

  const Mcontrol<URV>& ctl = data1_.mcontrol_;

  if (mode == PrivilegeMode::Machine and not ctl.m_)
    return false;  // Not enabled;

  if (mode == PrivilegeMode::Supervisor and not ctl.s_)
    return false;  // Not enabled;

  if (mode == PrivilegeMode::User and not ctl.u_)
    return false;  // Not enabled;

  if (mode == PrivilegeMode::Reserved)
    return false;  // Not enabled;

  bool isStore = not isLoad;
  bool clearBit0 = false;

  if (TriggerTiming(ctl.timing_) == timing and
      Select(ctl.select_) == Select::MatchData and
      ((isLoad and ctl.load_) or (isStore and ctl.store_)))
    return doMatch(value, clearBit0);

  return false;
}


template <typename URV>
bool
Trigger<URV>::doMatch(URV item, bool clearBit0) const
{
  URV data2 = data2_;
  if (clearBit0)
    {
      data2 = (data2 >> 1) << 1;
      item = (item >> 1) << 1;
    }

  switch (Match(data1_.mcontrol_.match_))
    {
    case Match::Equal:
      return item == data2;

    case Match::Masked:
      return (item & data2CompareMask_) == (data2 & data2CompareMask_);

    case Match::GE:
      return item >= data2;

    case Match::LT:
      return item < data2;

    case Match::MaskHighEqualLow:
      {
	unsigned halfBitCount = 4*sizeof(URV);
	// Mask low half of item with data2_ high half
	item = item & (data2 >> halfBitCount);
	// Compare low half
	return (item << halfBitCount) == (data2 << halfBitCount);
      }

    case Match::MaskLowEqualHigh:
      {
	unsigned halfBitCount = 4*sizeof(URV);
	// Mask high half of item with data2_ low half
	item = item & (data2 << halfBitCount);
	// Compare high half
	return (item >> halfBitCount) == (data2 >> halfBitCount);
      }
    }

  return false;
}


template <typename URV>
bool
Trigger<URV>::matchInstAddr(URV address, TriggerTiming timing,
                            PrivilegeMode mode) const
{
  if (not data1_.isAddrData())
    return false;  // Not an address trigger.

  const Mcontrol<URV>& ctl = data1_.mcontrol_;

  if (mode == PrivilegeMode::Machine and not ctl.m_)
    return false;  // Not enabled;

  if (mode == PrivilegeMode::Supervisor and not ctl.s_)
    return false;  // Not enabled;

  if (mode == PrivilegeMode::User and not ctl.u_)
    return false;  // Not enabled;

  if (mode == PrivilegeMode::Reserved)
    return false;  // Not enabled;

  bool clearBit0 = true;  // Clear bit0 of address before matching.
  if (TriggerTiming(ctl.timing_) == timing and
      Select(ctl.select_) == Select::MatchAddress and
      ctl.execute_)
    return doMatch(address, clearBit0);

  return false;
}


template <typename URV>
bool
Trigger<URV>::matchInstOpcode(URV opcode, TriggerTiming timing,
                              PrivilegeMode mode) const
{
  if (not data1_.isAddrData())
    return false;  // Not an address trigger.

  const Mcontrol<URV>& ctl = data1_.mcontrol_;

  if (mode == PrivilegeMode::Machine and not ctl.m_)
    return false;  // Not enabled;

  if (mode == PrivilegeMode::Supervisor and not ctl.s_)
    return false;  // Not enabled;

  if (mode == PrivilegeMode::User and not ctl.u_)
    return false;  // Not enabled;

  if (mode == PrivilegeMode::Reserved)
    return false;  // Not enabled;

  bool clearBit0 = false;

  if (TriggerTiming(ctl.timing_) == timing and
      Select(ctl.select_) == Select::MatchData and
      ctl.execute_)
    return doMatch(opcode, clearBit0);

  return false;
}


template class WdRiscv::Trigger<uint32_t>;
template class WdRiscv::Trigger<uint64_t>;

template class WdRiscv::Triggers<uint32_t>;
template class WdRiscv::Triggers<uint64_t>;
