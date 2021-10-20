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

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include "HartConfig.hpp"
#include "System.hpp"
#include "Core.hpp"
#include "Hart.hpp"


using namespace WdRiscv;


HartConfig::HartConfig()
{
  config_ = new nlohmann::json();
}


HartConfig::~HartConfig()
{
  delete config_;
  config_ = nullptr;
}


bool
HartConfig::loadConfigFile(const std::string& filePath)
{
  std::ifstream ifs(filePath);
  if (not ifs.good())
    {
      std::cerr << "Failed to open config file '" << filePath
		<< "' for input.\n";
      return false;
    }

  try
    {
      ifs >> *config_;
    }
  catch (std::exception& e)
    {
      std::cerr << e.what() << "\n";
      return false;
    }
  catch (...)
    {
      std::cerr << "Caught unknown exception while parsing "
		<< " config file '" << filePath << "'\n";
      return false;
    }

  return true;
}


namespace WdRiscv
{
  
  /// Convert given json entry to an unsigned integer value honoring
  /// hexadecimal prefix (0x) if any. Return true on succes and false
  /// if given entry does not represent an integer.
  template <typename URV>
  bool
  getJsonUnsigned(const std::string& tag, const nlohmann::json& js, URV& value)
  {
    value = 0;

    if (js.is_number())
      {
        value = js.get<URV>();
        return true;
      }

    if (js.is_string())
      {
	char *end = nullptr;
	std::string str = js.get<std::string>();
	uint64_t u64 = strtoull(str.c_str(), &end, 0);
	if (end and *end)
          {
            std::cerr << "Invalid config file value for '" << tag << "': "
                      << str << '\n';
            return false;
          }
	value = static_cast<URV>(u64);
	if (value != u64)
          {
            std::cerr << "Overflow in config file value for '" << tag << "': "
                      << str << '\n';
            return false;
          }

	return true;
      }

    std::cerr << "Config file entry '" << tag << "' must contain a number\n";
    return false;
  }


  /// Convert given json array value to a vector of unsigned integers
  /// honoring any hexadecimal prefix (0x) if any. Return true on
  /// sucess an false on failure.
  template <typename URV>
  bool
  getJsonUnsignedVec(const std::string& tag, const nlohmann::json& js,
                     std::vector<URV>& vec)
  {
    vec.clear();

    if (not js.is_array())
      {
	std::cerr << "Invalid config file value for '" << tag << "'"
		  << " -- expecting array of numbers\n";
	return false;
      }

    unsigned errors = 0;

    for (const auto& item :js)
      {
	if (item.is_number())
	  vec.push_back(item.get<unsigned>());
	else if (item.is_string())
	  {
	    char *end = nullptr;
	    std::string str = item.get<std::string>();
	    uint64_t u64 = strtoull(str.c_str(), &end, 0);
	    if (end and *end)
	      {
		std::cerr << "Invalid config file value for '" << tag << "': "
			  << str << '\n';
                errors++;
		continue;
	      }

	    URV val = static_cast<URV>(u64);
	    if (val != u64)
              {
                std::cerr << "Overflow in config file value for '" << tag
                          << "': " << str << '\n';
                errors++;
                continue;
              }

	    vec.push_back(val);
	  }
	else
          {
            std::cerr << "Invalid config file value for '" << tag << "'"
                      << " -- expecting array of number\n";
            errors++;
          }
      }

    return errors == 0;
  }


  /// Convert given json entry to a boolean value. Return ture on
  /// success and false on failure.
  bool
  getJsonBoolean(const std::string& tag, const nlohmann::json& js, bool& value)
  {
    value = false;

    if (js.is_boolean())
      {
        value = js.get<bool>();
        return true;
      }

    if (js.is_number())
      {
        value = js.get<unsigned>() != 0;
        return true;
      }

    if (js.is_string())
      {
	std::string str = js.get<std::string>();
	if (str == "0" or str == "false" or str == "False")
          value = false;
	else if (str == "1" or str == "true" or str == "True")
          value = true;
        else
          {
            std::cerr << "Invalid config file value for '" << tag << "': "
                      << str << '\n';
            return false;
          }
        return true;
      }

    std::cerr << "Config file entry '" << tag << "' must contain a bool\n";
    return false;
  }

}


static
bool
validateStackChecker(const nlohmann::json& csrs)
{
  // If any of the stack checker CSRs is present then all must be
  // present.
  auto tags = { "mspcba", "mspcta", "mspcc" };
  std::string present;
  unsigned count = 0;
  for (const auto& tag : tags)
    if (csrs.count(tag))
      {
	present = tag;
	count++;
      }

  if (count == 0)
    return true;

  if (count != tags.size())
    {
      std::cerr << "Error: Not all stack checker CSRs are defined:\n";
      std::cerr << "  Defined: ";
      std::string sep = "";
      for (const auto& tag: tags)
	if (csrs.count(tag))
	  {
	    std::cerr << sep << tag;
	    sep = ", ";
	  }

      sep.clear();
      std::cerr << "  Missing: ";
      for (const auto& tag: tags)
	if (not csrs.count(tag))
	  {
	    std::cerr << sep << tag;
	    sep = ", ";
	  }

      return false;
    }

  return true;
}


template <typename URV>
static
bool
applyCsrConfig(Hart<URV>& hart, const nlohmann::json& config, bool verbose)
{
  if (not config.count("csr"))
    return true;  // Nothing to apply

  const auto& csrs = config.at("csr");
  if (not csrs.is_object())
    {
      std::cerr << "Invalid csr entry in config file (expecting an object)\n";
      return false;
    }

  unsigned errors = 0;
  for (auto it = csrs.begin(); it != csrs.end(); ++it)
    {
      const std::string& csrName = it.key();
      const auto& conf = it.value();

      URV reset = 0, mask = 0, pokeMask = 0;
      bool isDebug = false, exists = true, shared = false;

      Csr<URV>* csr = hart.findCsr(csrName);
      if (csr)
	{
	  reset = csr->getResetValue();
	  mask = csr->getWriteMask();
	  pokeMask = csr->getPokeMask();
	  isDebug = csr->isDebug();
	}

      if (conf.count("reset"))
        if (not getJsonUnsigned(csrName + ".reset", conf.at("reset"), reset))
          errors++;

      if (conf.count("mask"))
	{
	  if (not getJsonUnsigned(csrName + ".mask", conf.at("mask"), mask))
            errors++;

	  // If defining a non-standard CSR (as popposed to
	  // configuring an existing CSR) then default the poke-mask
	  // to the write-mask.
	  if (not csr)
	    pokeMask = mask;
	}

      if (conf.count("poke_mask"))
	if (not getJsonUnsigned(csrName + ".poke_mask", conf.at("poke_mask"), pokeMask))
          errors++;

      if (conf.count("debug"))
	if (not getJsonBoolean(csrName + ".bool", conf.at("debug"), isDebug))
          errors++;

      if (conf.count("exists"))
	if (not getJsonBoolean(csrName + ".bool", conf.at("exists"), exists))
          errors++;

      if (conf.count("shared"))
        if (not getJsonBoolean(csrName + ".bool", conf.at("shared"), shared))
          errors++;

      // If number present and csr is not defined, then define a new
      // CSR; otherwise, configure.
      if (conf.count("number"))
	{
          unsigned number = 0;
	  if (not getJsonUnsigned<unsigned>(csrName + ".number", conf.at("number"), number))
            errors++;
          else
            {
              if (csr)
                {
                  if (csr->getNumber() != CsrNumber(number))
                    {
                      std::cerr << "Invalid config file entry for CSR "
                                << csrName << ": Number (0x" << std::hex << number
                                << ") does not match that of previous definition ("
                                << "0x" << unsigned(csr->getNumber())
                                << ")\n" << std::dec;
                      errors++;
                      continue;
                    }
                  // If number matches we configure below
                }
              else if (hart.defineCsr(csrName, CsrNumber(number), exists,
                                      reset, mask, pokeMask, isDebug))
                {
                  csr = hart.findCsr(csrName);
                  assert(csr);
                }
              else
                {
                  std::cerr << "Invalid config file CSR definition with name "
                            << csrName << " and number 0x" << std::hex << number
                            << ": Number already in use\n" << std::dec;
                  errors++;
                  continue;
                }
            }
        }

      if (not csr)
	{
	  std::cerr << "A CSR number must be provided in configuration of non-standard CSR "
		    << csrName << '\n';
	  errors++;
	  continue;
	}
      bool exists0 = csr->isImplemented(), isDebug0 = csr->isDebug();
      bool shared0 = csr->isShared();
      URV reset0 = csr->getResetValue(), mask0 = csr->getWriteMask();
      URV pokeMask0 = csr->getPokeMask();

      if (csrName == "mhartstart")
        if (hart.sysHartIndex() == 0 and (reset & 1) == 0)
          std::cerr << "Warning: Bit corresponding to hart 0 is cleared "
                    << "in reset value of mhartstart CSR -- Bit is ignored\n";

      if (csrName == "mhartid")
        {
          std::cerr << "CSR mhartid cannot be configured.\n";
          std::cerr << "Ignoring mhartid CSR configuration in config file.\n";
          continue;
        }

      if (not hart.configCsr(csrName, exists, reset, mask, pokeMask,
			     isDebug, shared))
	{
	  std::cerr << "Invalid CSR (" << csrName << ") in config file.\n";
	  errors++;
	}
      else if (verbose)
	{
	  if (exists0 != exists or isDebug0 != isDebug or reset0 != reset or
	      mask0 != mask or pokeMask0 != pokeMask)
	    {
	      std::cerr << "Configuration of CSR (" << csrName <<
		") changed in config file:\n";

	      if (exists0 != exists)
		std::cerr << "  implemented: " << exists0 << " to "
			  << exists << '\n';

	      if (isDebug0 != isDebug)
		std::cerr << "  debug: " << isDebug0 << " to "
			  << isDebug << '\n';

	      if (shared0 != shared)
		std::cerr << "  shred: " << shared0 << " to "
			  << shared << '\n';

	      if (reset0 != reset)
		std::cerr << "  reset: 0x" << std::hex << reset0
			  << " to 0x" << reset << '\n' << std::dec;

	      if (mask0 != mask)
		std::cerr << "  mask: 0x" << std::hex << mask0
			  << " to 0x" << mask << '\n' << std::dec;

	      if (pokeMask0 != pokeMask)
		std::cerr << "  poke_mask: " << std::hex << pokeMask0
			  << " to 0x" << pokeMask << '\n' << std::dec;
	    }
	}
    }

  // Stack checker.
  if (not validateStackChecker(csrs))
    errors++;

  return errors == 0;
}


template <typename URV>
static
bool
applyMemMappedRegConfig(Hart<URV>& hart, const nlohmann::json& config)
{
  if (not config.count("memory_mapped_registers"))
    return true;  // Nothing to apply.

  const auto& mmr = config.at("memory_mapped_registers");

  unsigned errors = 0;


  // Define memory-mapped-register region.
  uint64_t addr = 0, size = 0;
  if (getJsonUnsigned("address", mmr.at("address"), addr) and
      getJsonUnsigned("size", mmr.at("size"), size))
    {
      if (not hart.defineMemoryMappedRegisterArea(addr, size))
        return false;
    }
  else
    errors++;

  // Start by giving all registers in region a default mask.
  size_t possibleRegCount = size / 4;
  if (mmr.count("default_mask"))
    {
      uint32_t mask = 0;
      if (getJsonUnsigned("default_mask", mmr.at("default_mask"), mask))
        for (size_t ix = 0; ix < possibleRegCount; ++ix)
          hart.defineMemoryMappedRegisterWriteMask(addr + ix*4, mask);
      else
        errors++;
    }

  if (not mmr.count("registers"))
    return true;

  const auto& regs = mmr.at("registers");
  if (not regs.is_object())
    {
      std::cerr << "Invalid memory_mapped_registers.registers entry "
                << "in config file (expecting an object)\n";
      return false;
    }

  for (auto it = regs.begin(); it != regs.end(); ++it)
    {
      const std::string& name = it.key();
      const auto& conf = it.value();

      if (not conf.count("count") or not conf.count("address") or
          not conf.count("mask"))
        {
          std::cerr << "Register entry \"" << name << "\"under "
                    << "memory_mapped_registers must have count/address/mask\n";
          errors++;
          continue;
        }
      URV count = 0, mask = 0, addr = 0;
      if (getJsonUnsigned("count", conf.at("count"), count) and
          getJsonUnsigned("mask", conf.at("mask"), mask) and
          getJsonUnsigned("address", conf.at("address"), addr))
        {
          for (URV ix = 0; ix < count; ++ix)
            if (not hart.defineMemoryMappedRegisterWriteMask(addr + ix*4, mask))
              errors++;
        }
      else
        errors++;
    }

  return errors == 0;
}


template <typename URV>
static
bool
applyIccmConfig(Hart<URV>& hart, const nlohmann::json& config)
{
  if (not config.count("iccm"))
    return true;

  const auto& iccm = config.at("iccm");
  if (iccm.count("region") and iccm.count("size") and iccm.count("offset"))
    {
      size_t region = 0, size = 0, offset = 0;
      if (getJsonUnsigned("iccm.region", iccm.at("region"), region) and
          getJsonUnsigned("iccm.size",   iccm.at("size"), size) and
          getJsonUnsigned("iccm.offset", iccm.at("offset"), offset))
        {
          size_t addr = region*hart.regionSize() + offset;
          return hart.defineIccm(addr, size);
        }
      else
        return false;
    }

  std::cerr << "The ICCM entry in the configuration file must contain "
            << "a region, offset and a size entry.\n";
  return false;
}


template <typename URV>
static
bool
applyDccmConfig(Hart<URV>& hart, const nlohmann::json& config)
{
  if (not config.count("dccm"))
    return true;

  const auto& dccm = config.at("dccm");
  if (dccm.count("region") and dccm.count("size") and dccm.count("offset"))
    {
      size_t region = 0, size = 0, offset = 0;
      if (getJsonUnsigned("dccm.region", dccm.at("region"), region) and
          getJsonUnsigned("dccm.size",   dccm.at("size"), size) and
          getJsonUnsigned("dccm.offset", dccm.at("offset"), offset))
        {
          size_t addr = region*hart.regionSize() + offset;
          return hart.defineDccm(addr, size);
        }
      return false;
    }

  std::cerr << "The DCCM entry in the configuration file must contain "
            << "a region, offset and a size entry.\n";
  return false;
}


template <typename URV>
static
bool
applyTriggerConfig(Hart<URV>& hart, const nlohmann::json& config)
{
  if (not config.count("triggers"))
    return true;  // Nothing to apply

  const auto& triggers = config.at("triggers");
  if (not triggers.is_array())
    {
      std::cerr << "Invalid triggers entry in config file (expecting an array)\n";
      return false;
    }

  unsigned errors = 0;
  unsigned ix = 0;
  for (auto it = triggers.begin(); it != triggers.end(); ++it, ++ix)
    {
      const auto& trig = *it;
      std::string name = std::string("trigger") + std::to_string(ix);
      if (not trig.is_object())
	{
	  std::cerr << "Invalid trigger in config file triggers array "
		    << "(expecting an object at index " << ix << ")\n";
	  ++errors;
	  break;
	}
      bool ok = true;
      for (const auto& tag : {"reset", "mask", "poke_mask"})
	if (not trig.count(tag))
	  {
	    std::cerr << "Trigger " << name << " has no '" << tag
		      << "' entry in config file\n";
	    ok = false;
	  }
      if (not ok)
	{
	  errors++;
	  continue;
	}

      std::vector<URV> resets, masks, pokeMasks;
      ok = (getJsonUnsignedVec(name + ".reset", trig.at("reset"), resets) and
            getJsonUnsignedVec(name + ".mask", trig.at("mask"), masks) and
            getJsonUnsignedVec(name + ".poke_mask", trig.at("poke_mask"), pokeMasks));
      if (not ok)
        {
          errors++;
          continue;
        }

      if (resets.size() != 3)
	{
	  std::cerr << "Trigger " << name << ": Bad item count (" << resets.size()
		    << ") for 'reset' field in config file. Expecting 3.\n";
	  ok = false;
	}

      if (masks.size() != 3)
	{
	  std::cerr << "Trigger " << name << ": Bad item count (" << masks.size()
		    << ") for 'mask' field in config file. Expecting 3.\n";
	  ok = false;
	}

      if (pokeMasks.size() != 3)
	{
	  std::cerr << "Trigger " << name << ": Bad item count (" << pokeMasks.size()
		    << ") for 'poke_mask' field in config file. Expecting 3.\n";
	  ok = false;
	}

      if (not ok)
	{
	  errors++;
	  continue;
	}
      if (not hart.configTrigger(ix, resets.at(0), resets.at(1), resets.at(2),
				 masks.at(0), masks.at(1), masks.at(2),
				 pokeMasks.at(0), pokeMasks.at(1),
                                 pokeMasks.at(2)))
	{
	  std::cerr << "Failed to configure trigger " << std::dec << ix << '\n';
	  ++errors;
	}
    }

  return errors == 0;
}


template <typename URV>
static
bool
applyPerfEvents(Hart<URV>& hart, const nlohmann::json& config,
                bool userMode, bool /*verbose*/)
{
  unsigned errors = 0;

  std::string tag = "num_mmode_perf_regs";
  if (config.count(tag))
    {
      unsigned count = 0;
      if (not getJsonUnsigned<unsigned>(tag, config.at(tag), count))
        errors++;
      else
        {
          if (not hart.configMachineModePerfCounters(count))
            errors++;
          if (userMode)
            if (not hart.configUserModePerfCounters(count))
              errors++;
        }
    }

  unsigned maxPerfId = 0;
  tag = "max_mmode_perf_event";
  if (config.count(tag))
    {
      if (not getJsonUnsigned<unsigned>(tag, config.at(tag), maxPerfId))
        errors++;
      else
        {
          unsigned limit = 16*1024;
          if (maxPerfId > limit)
            {
              std::cerr << "Config file max_mmode_perf_event too large -- Using "
                        << limit << '\n';
              maxPerfId = limit;
            }
          hart.configMachineModeMaxPerfEvent(maxPerfId);
        }
    }

  tag = "mmode_perf_events";
  if (config.count(tag))
    {
      std::vector<unsigned> eventsVec;

      const auto& events = config.at(tag);
      if (not events.is_array())
        {
          std::cerr << "Invalid mmode_perf_events entry in config file (expecting an array)\n";
          errors++;
        }
      else
        {
          unsigned ix = 0;
          for (auto it = events.begin(); it != events.end(); ++it, ++ix)
            {
              const auto& event = *it;
              std::string elemTag = tag + "element " + std::to_string(ix);
              unsigned eventId = 0;
              if (not getJsonUnsigned<unsigned>(elemTag, event, eventId))
                errors++;
              else
                eventsVec.push_back(eventId);
            }
        }
      hart.configPerfEvents(eventsVec);
    }

  return errors == 0;
}


template <typename URV>
static
bool
applyVectorConfig(Hart<URV>& hart, const nlohmann::json& config)
{
  if (not config.count("vector"))
    return true;  // Nothing to apply

  unsigned errors = 0;
  const auto& vconf = config.at("vector");

  unsigned bytesPerVec = 0;
  std::string tag = "bytes_per_vec";
  if (not vconf.count(tag))
    {
      std::cerr << "Error: Missing " << tag << " tag in vector section of config file\n";
      errors++;
    }
  else
    {
      if (not getJsonUnsigned(tag, vconf.at(tag), bytesPerVec))
        errors++;
      else if (bytesPerVec == 0 or bytesPerVec > 4096)
        {
          std::cerr << "Error: Invalid config file bytes_per_vec number: "
                    << bytesPerVec << '\n';
          errors++;
        }
      else
        {
          unsigned l2BytesPerVec = std::log2(bytesPerVec);
          unsigned p2BytesPerVec = uint32_t(1) << l2BytesPerVec;
          if (p2BytesPerVec != bytesPerVec)
            {
              std::cerr << "Error: Config file bytes_per_vec ("
                        << bytesPerVec << ") is not a power of 2\n";
              errors++;
            }
        }
    }

  unsigned bytesPerElem = 0;
  tag = "max_bytes_per_elem";
  if (not vconf.count(tag))
    {
      std::cerr << "Error: Missing " << tag << " tag in vector section of config file\n";
      errors++;
    }
  else
    {
      if (not getJsonUnsigned(tag, vconf.at(tag), bytesPerElem))
        errors++;
      else if (bytesPerElem == 0 or bytesPerElem > bytesPerVec)
        {
          std::cerr << "Error: Invalid config file max_bytes_per_elem number: "
                    << bytesPerElem << '\n';
          errors++;
        }
      else
        {
          unsigned l2BytesPerElem = std::log2(bytesPerElem);
          unsigned p2BytesPerElem = uint32_t(1) << l2BytesPerElem;
          if (p2BytesPerElem != bytesPerElem)
            {
              std::cerr << "Error: Config file max_bytes_per_elem ("
                        << bytesPerElem << ") is not a power of 2\n";
              errors++;
            }
        }
    }

  if (errors == 0)
    hart.configVector(bytesPerVec, bytesPerElem);

  return errors == 0;
}


template <typename URV>
static
bool
applyInstMemConfig(Hart<URV>& hart, const nlohmann::json& config)
{
  using std::cerr;

  if (not config.is_array())
    {
      cerr << "Invalid inst entry in config file memmap (execpting an array)\n";
      return false;
    }

  // Pairs of addresses in which inst fetch is allowed.
  std::vector<std::pair<URV,URV>> windows;

  unsigned errors = 0;
  unsigned ix = 0;
  for (auto it = config.begin(); it != config.end(); ++it, ++ix)
    {
      const auto& addrPair = *it;
      if (not addrPair.is_array() or addrPair.size() != 2)
	{
	  cerr << "Invalid address range in config file memap inst ("
	       << "expecting an array of 2 numbers at index " << ix << ")\n";
	  ++errors;
	  break;
	}

      std::vector<URV> vec;
      if (not getJsonUnsignedVec("memmap.inst", addrPair, vec))
        errors++;
      else if (vec.size() != 2)
	errors++;
      else
	{
	  auto p = std::make_pair(vec.at(0), vec.at(1));
	  windows.push_back(p);
	}
    }

  if (not errors and not windows.empty())
    if (not hart.configMemoryFetch(windows))
      errors++;

  return errors == 0;
}


template <typename URV>
static
bool
applyDataMemConfig(Hart<URV>& hart, const nlohmann::json& config)
{
  using std::cerr;

  if (not config.is_array())
    {
      cerr << "Invalid data entry in config file memmap (execpting an array)\n";
      return false;
    }

  // Pairs of addresses in which data access is allowed.
  std::vector<std::pair<URV,URV>> windows;

  unsigned errors = 0;
  unsigned ix = 0;
  for (auto it = config.begin(); it != config.end(); ++it, ++ix)
    {
      const auto& addrPair = *it;
      if (not addrPair.is_array() or addrPair.size() != 2)
	{
	  cerr << "Invalid address range in config file memap data ("
	       << "expecting an array of 2 numbers at index " << ix << ")\n";
	  ++errors;
	  break;
	}

      std::vector<URV> vec;
      if (not getJsonUnsignedVec("memmap.data", addrPair, vec))
        errors++;
      if (vec.size() != 2)
	errors++;
      else
	{
	  auto p = std::make_pair(vec.at(0), vec.at(1));
	  windows.push_back(p);
	}
    }

  if (not errors and not windows.empty())
    if (not hart.configMemoryDataAccess(windows))
      errors++;

  return errors == 0;
}


template<typename URV>
bool
HartConfig::applyMemoryConfig(Hart<URV>& hart, bool iccmRw, bool /*verbose*/) const
{
  unsigned errors = 0;
  if (not applyIccmConfig(hart, *config_))
    errors++;

  if (not applyDccmConfig(hart, *config_))
    errors++;

  if (not applyMemMappedRegConfig(hart, *config_))
    errors++;

  if (config_ -> count("iccmrw"))
    if (not getJsonBoolean("iccmrw", config_ ->at("iccmrw"), iccmRw))
      errors++;

  hart.finishCcmConfig(iccmRw);

  if (config_ -> count("memmap"))
    {
      // Apply memory protection windows.
      const auto& memmap = config_ -> at("memmap");
      std::string tag = "inst";
      if (memmap.count(tag))
	if (not applyInstMemConfig(hart, memmap.at(tag)))
	  errors++;
      tag = "data";
      if (memmap.count(tag))
	if (not applyDataMemConfig(hart, memmap.at(tag)))
	  errors++;
    }

  return errors == 0;
}


template<typename URV>
bool
HartConfig::applyConfig(Hart<URV>& hart, bool userMode, bool verbose) const
{
  unsigned errors = 0;

  // Define PC value after reset.
  std::string tag = "reset_vec";
  if (config_ -> count(tag))
    {
      URV resetPc = 0;
      if (getJsonUnsigned(tag, config_ -> at(tag), resetPc))
        hart.defineResetPc(resetPc);
      else
        errors++;
    }

  // Define non-maskable-interrupt pc
  tag = "nmi_vec";
  if (config_ -> count(tag))
    {
      URV nmiPc = 0;
      if (getJsonUnsigned(tag, config_ -> at(tag), nmiPc))
        hart.defineNmiPc(nmiPc);
      else
        errors++;
    }

  // Use ABI register names (e.g. sp instead of x2).
  tag = "abi_names";
  if (config_ -> count(tag))
    {
      bool abiNames = false;
      if (getJsonBoolean(tag, config_ -> at(tag), abiNames))
        hart.enableAbiNames(abiNames);
      else
        errors++;
    }

  // Print memory address of load/store instruction in trace log.
  tag = "print_load_store_address";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.setTraceLoadStore(flag);
      else
        errors++;
    }

  // Atomic instructions illegal outside of DCCM.
  tag = "amo_illegal_outside_dccm";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ ->at(tag), flag))
        hart.setAmoInDccmOnly(flag);
      else
        errors++;
    }

  // Atomic instructions illegal outside cacheable regions.
  tag = "amo_illegal_outside_cacheable";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ ->at(tag), flag))
        hart.setAmoInCacheableOnly(flag);
      else
        errors++;
    }

  // Wide (64-bit) load/store. WDC special.
  tag = "enable_wide_load_store";
  if (config_ -> count(tag))
    {
      bool flag = true;
      if (getJsonBoolean(tag, config_ ->at(tag), flag))
        hart.enableWideLoadStore(flag);
      else
        errors++;
    }

  // Bus barrier. WDC special.
  tag = "enable_bus_barrier";
  if (config_ -> count(tag))
    {
      bool flag = true;
      if (getJsonBoolean(tag, config_ ->at(tag), flag))
        hart.enableBusBarrier(flag);
      else
        errors++;
    }

  // Ld/st instructions trigger misaligned exception if base address
  // (value in rs1) and effective address refer to regions of
  // different types.
  tag = "effective_address_compatible_with_base";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ ->at(tag), flag))
        hart.setEaCompatibleWithBase(flag);
      else
        errors++;
    }

  // Enable debug triggers.
  tag = "enable_triggers";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ ->at(tag), flag))
        hart.enableTriggers(flag);
      else
        errors++;
    }

  // Enable performance counters.
  tag = "enable_performance_counters";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ ->at(tag), flag))
        hart.enablePerformanceCounters(flag);
      else
        errors++;
    }

  // Enable rollback of memory on store error.
  tag = "store_error_rollback";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.enableStoreErrorRollback(flag);
      else
        errors++;
    }

  // Enable rollback of register on load error.
  tag = "load_error_rollback";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        {
          hart.enableLoadErrorRollback(flag);
          hart.enableBenchLoadExceptions(flag);
        }
      else
        errors++;
    }

  // Enable fast interrupts.
  tag = "fast_interrupt_redirect";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.enableFastInterrupts(flag);
      else
        errors++;
    }

  tag = "enable_per_mode_counter_control";
  if (config_ -> count(tag))
    {
      bool flag = true;
      if (getJsonBoolean(tag, config_ ->at(tag), flag))
        hart.enablePerModeCounterControl(flag);
      else
        errors++;
    }

  // Enable zbb.
  tag = "enable_zbmini";
  if (config_ -> count(tag))
    {
      std::cerr << "Config file tag \"enable_zbmini\" deprecated: "
		<< "Using \"enable_zbb\" and \"enable_zbs\"\n";
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        {
          hart.enableRvzbb(flag);
          hart.enableRvzbs(flag);
        }
      else
        errors++;

    }

  tag = "enable_zba";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.enableRvzba(flag);
      else
        errors++;
    }

  tag = "enable_zbb";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.enableRvzbb(flag);
      else
        errors++;
    }

  tag = "enable_zbc";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.enableRvzbc(flag);
      else
        errors++;
    }

  tag = "enable_zbe";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.enableRvzbe(flag);
      else
        errors++;
    }

  tag = "enable_zbf";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.enableRvzbf(flag);
      else
        errors++;
    }

  tag = "enable_zbm";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.enableRvzbm(flag);
      else
        errors++;
    }

  tag = "enable_zbp";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.enableRvzbp(flag);
      else
        errors++;
    }

  tag = "enable_zbr";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.enableRvzbr(flag);
      else
        errors++;
    }

  tag = "enable_zbs";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.enableRvzbs(flag);
      else
        errors++;
    }

  tag = "enable_zbt";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.enableRvzbt(flag);
      else
        errors++;
    }

  tag = "enable_zfh";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.enableZfh(flag);
      else
        errors++;
    }

  tag = "load_queue_size";
  if (config_ -> count(tag))
    {
      unsigned lqs = 0;
      if (getJsonUnsigned(tag, config_ -> at(tag), lqs))
        {
          if (lqs > 64)
            {
              std::cerr << "Config file load queue size (" << lqs << ") too large"
                        << " -- using 64.\n";
              lqs = 64;
            }
          hart.setLoadQueueSize(lqs);
        }
      else
        errors++;
    }

  tag = "even_odd_trigger_chains";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.configEvenOddTriggerChaining(flag);
      else
        errors++;
    }

  if (not applyPerfEvents(hart, *config_, userMode, verbose))
    errors++;

  if (not applyCsrConfig(hart, *config_, verbose))
    errors++;

  if (not applyTriggerConfig(hart, *config_))
    errors++;

  if (not applyVectorConfig(hart, *config_))
    errors++;

  tag = "load_data_trigger";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.configLoadDataTrigger(flag);
      else
        errors++;
    }

  tag = "exec_opcode_trigger";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (getJsonBoolean(tag, config_ -> at(tag), flag))
        hart.configExecOpcodeTrigger(flag);
      else
        errors++;
    }

  tag = "memmap";
  if (config_ -> count(tag))
    {
      const auto& memmap = config_ -> at(tag);
      tag = "consoleio";
      if (memmap.count(tag))
	{
          URV io = 0;
	  if (getJsonUnsigned("memmap.consoleio", memmap.at(tag), io))
            hart.setConsoleIo(io);
          else
            errors++;
	}
    }

  tag = "syscall_slam_area";
  if (config_ -> count(tag))
    {
      uint64_t addr = 0;
      if (getJsonUnsigned(tag, config_ -> at(tag), addr))
        {
          hart.defineSyscallSlam(addr);
          hart.enableLinux(true);
        }
      else
        errors++;
    }

  tag = "physical_memory_protection_grain";
  if (config_ -> count(tag))
    {
      uint64_t size = 0;
      if (getJsonUnsigned<uint64_t>(tag, config_ -> at(tag), size))
        hart.configMemoryProtectionGrain(size);
      else
        errors++;
    }

  tag = "enable_misaligned_data";
  if (config_ -> count(tag))
    {
      bool misal = true;
      if (getJsonBoolean(tag, config_ ->at(tag), misal))
        hart.enableMisalignedData(misal);
      else
        errors++;
    }

  tag = "force_rounding_mode";
  if (config_ -> count(tag))
    {
      std::string str = config_->at(tag).get<std::string>();
      if (str == "rne")
	hart.forceRoundingMode(RoundingMode::NearestEven);
      else if (str == "rtz")
	hart.forceRoundingMode(RoundingMode::Zero);
      else if (str == "rdn")
	hart.forceRoundingMode(RoundingMode::Down);
      else if (str == "rup")
	hart.forceRoundingMode(RoundingMode::Up);
      else if (str == "rmm")
	hart.forceRoundingMode(RoundingMode::NearestMax);
      else
	{
	  std::cerr << "Invalid force_rounding_mode config: " << str << '\n';
	  errors++;
	}
    }

  tag = "enable_csv_log";
  if (config_ -> count(tag))
    {
      bool flag = false;
      if (not getJsonBoolean(tag, config_ -> at(tag), flag))
	errors++;
      else
	hart.enableCsvLog(flag);
    }

  return errors == 0;
}


template<typename URV>
bool
HartConfig::configHarts(System<URV>& system, bool userMode,
                        bool verbose) const
{
  userMode = userMode or this->userModeEnabled();

  // Apply JSON configuration.
  for (unsigned i = 0; i < system.hartCount(); ++i)
    if (not applyConfig(*system.ithHart(i), userMode, verbose))
      return false;

  return finalizeCsrConfig(system);
}


template<typename URV>
bool
HartConfig::configMemory(System<URV>& system, bool iccmRw, bool unmappedElfOk,
                         bool verbose) const
{
  system.checkUnmappedElf(not unmappedElfOk);

  auto& hart0 = *system.ithHart(0);
  if (not applyMemoryConfig(hart0, iccmRw, verbose))
    return false;

  for (unsigned i = 1; i < system.hartCount(); ++i)
    system.ithHart(i)->copyMemRegionConfig(hart0);

  return true;
}


bool
HartConfig::getXlen(unsigned& xlen) const
{
  if (config_ -> count("xlen"))
    return getJsonUnsigned("xlen", config_ -> at("xlen"), xlen);
  return false;
}


bool
HartConfig::getCoreCount(unsigned& count) const
{
  if (config_ -> count("cores"))
    return getJsonUnsigned("cores", config_ -> at("cores"), count);
  return false;
}


bool
HartConfig::getHartsPerCore(unsigned& count) const
{
  if (config_ -> count("harts"))
    return getJsonUnsigned("harts", config_ -> at("harts"), count);
  return false;
}


bool
HartConfig::getPageSize(size_t& pageSize) const
{
  if (not config_ -> count("memmap"))
    return false;

  auto& mem = config_ -> at("memmap");
  if (not mem.count("page_size"))
    return false;

  return getJsonUnsigned("memmap.page_size", mem.at("page_size"), pageSize);
}


bool
HartConfig::getHartIdOffset(unsigned& offset) const
{
  std::string tag = "core_hart_id_offset";
  if (not config_ -> count(tag))
    return false;

  return getJsonUnsigned(tag, config_->at(tag), offset);
}


bool
HartConfig::getMemorySize(size_t& memSize) const
{
  if (not config_ -> count("memmap"))
    return false;

  auto& mem = config_ -> at("memmap");
  if (not mem.count("size"))
    return false;

  return getJsonUnsigned("memmap.size", mem.at("size"), memSize);
}


bool
HartConfig::userModeEnabled() const
{
  uint64_t resetVal = 0;
  if (not getMisaReset(resetVal))
    return false;

  // User mode enabled if bit corresponding to U extension is set.
  return ((uint64_t(1) << ('u' - 'a')) & resetVal) != 0;
}


bool
HartConfig::supervisorModeEnabled() const
{
  uint64_t resetVal = 0;
  if (not getMisaReset(resetVal))
    return false;

  // User mode enabled if bit corresponding to S extension is set.
  return ((uint64_t(1) << ('s' - 'a')) & resetVal) != 0;
}


void
HartConfig::clear()
{
  config_ -> clear();
}


// Meaning of the maco bits depend on the desing.
// For 32-bit design we have Maco32Masks; Otherwise Maco64Masks.
enum class Maco32Masks : uint32_t
  {
   Enable       = 1,
   Lock         = 2,
   SideEffect   = 4,
   PreciseStore = 8
  };

enum class Maco64Masks : uint32_t
  {
   Enable       = 1,
   Lock         = 2,
   Cacheable    = 4,
   SideEffect   = 8
  };


template <typename URV>
void
unpackMacoValue(URV value, URV mask, bool rv32, uint64_t& start, uint64_t& end,
                bool& idempotent, bool& cacheable)
{
  start = end = 0;
  bool enable = false;

  if (rv32)
    {
      idempotent = (value & URV(Maco32Masks::SideEffect)) == 0;
      cacheable = false;
      enable = value & URV(Maco32Masks::Enable);
    }
  else
    {
      idempotent = (value & URV(Maco64Masks::SideEffect)) == 0;
      cacheable = (value & URV(Maco64Masks::Cacheable));
      enable = value & URV(Maco64Masks::Enable);
    }

  if (not enable)
    return;  // Disabled: start same as end

  // We want first zero starting from bit 7 and going towards most sig bit.
  value |= 0x7f;  // Set least sig 7 bits to 1
  if (((value & mask) >> 7) == ((mask | 0x7f) >> 7))
    {
      start = end = 0;   // Illegal setup: Address bits all set.
      return;
    }

  unsigned rzi = __builtin_ctzl(~value);  // Rightmost-zero-bit-index.
  start = (value >> rzi) << rzi; // Clear bits below rightmost zero bit.
  uint64_t sizeM1 = (uint64_t(1) << (rzi + 1)) - 1;  // Size minus 1.
  end = start + sizeM1;
}


/// Associate callbacks with write/poke of the maco0-maco4 CSRs.
template <typename URV>
void
defineMacoSideEffects(System<URV>& system)
{
  bool rv32 = sizeof(URV) == 4;

  for (unsigned i = 0; i < system.hartCount(); ++i)
    {
      auto hart = system.ithHart(i);

      // Maco registers, if present, start at maco0 and are defined
      // sequentially.  We allow up to 32.
      unsigned macoIx = 0;
      for ( ; macoIx < 32; ++macoIx)
        {
          auto name = std::string("maco") + std::to_string(macoIx);
          auto csrPtr = hart->findCsr(name);
          if (not csrPtr)
            break;

          // Define pre-write/pre-poke callback: If lock bit is set,
          // then preserve CSR value
          auto pre = [hart, rv32] (Csr<URV>& csr, URV& val) -> void {
                       URV previous = 0;
                       hart->peekCsr(csr.getNumber(), previous);
                       if (previous & URV(Maco32Masks::Lock))
                         val = previous;  // Locked: keep previous value
                       else if (not rv32)
                         {
                           // Combination side-effect/cacable not allowed
                           bool side = val & URV(Maco64Masks::SideEffect);
                           bool cache = val & URV(Maco64Masks::Cacheable);
                           if (cache and side)
                             val &= ~URV(Maco64Masks::Cacheable);
                         }
                      };

          // Define post-write/post-poke callback. Upddate the idempotent
          // regions of the hart.
          auto post = [hart, macoIx, rv32] (Csr<URV>& csr, URV val) -> void {
                        uint64_t start = 0, end = 0;
                        bool idempotent = false, cacheable = false;
                        URV mask = csr.getWriteMask();
                        unpackMacoValue(val, mask, rv32, start, end, idempotent, cacheable);
                        hart->definePmaOverride(macoIx, start, end, idempotent, cacheable);
                      };

          csrPtr->registerPrePoke(pre);
          csrPtr->registerPreWrite(pre);

          csrPtr->registerPostPoke(post);
          csrPtr->registerPostWrite(post);
        }

      hart->definePmaOverrideRegions(macoIx);
    }
}


/// Associate callbacks with write/poke of the mrac CSR. Mrac is
/// memory region accss control CSR which defines which regions are
/// cacahble/idempotent.
template <typename URV>
void
defineMracSideEffects(System<URV>& system)
{
  unsigned count = 0; // Count of mrac definitions.

  for (unsigned i = 0; i < system.hartCount(); ++i)
    {
      auto hart = system.ithHart(i);
      auto csrPtr = hart->findCsr("mrac");
      if (not csrPtr)
        continue;

      count++;

      auto pre = [] (Csr<URV>&, URV& val) -> void {
                   // A value of 0b11 (io/cacheable) for the ith
                   // region is invalid: Make it 0b10
                   // (io/non-cacheable).
                   URV mask = 0b11;
                   unsigned xlen = sizeof(URV) * 8; // FIX.
                   for (unsigned i = 0; i < xlen; i += 2)
                     {
                       if ((val & mask) == mask)
                         val = (val & ~mask) | (0b10 << i);
                       mask = mask << 2;
                     }
                 };

      // Mrac is shared between harts. If one hart writes it, all
      // harts must be updated.
      auto post = [&system] (Csr<URV>&, URV val) -> void {
                    for (unsigned hix = 0; hix < system.hartCount(); ++hix)
                      {
                        auto ht = system.ithHart(hix);
                        for (unsigned region = 0; region < 16; ++region)
                          {
                            unsigned bit = (val >> (region*2 + 1)) & 1;
                            ht->markRegionIdempotent(region, bit == 0);
                          }
                      }
                  };

      csrPtr->registerPrePoke(pre);
      csrPtr->registerPreWrite(pre);
      csrPtr->registerPostPoke(post);
      csrPtr->registerPostWrite(post);
    }

  if (count and sizeof(URV) == 8)
    {
      std::cerr << "Warning: mrac CSR is defined with rv32. This is likely "
                << "a mistake. Only least significant 32 bits will be used.\n";
    }
}


/// Associate callbacks with write/poke of the mdac CSR. This is the
/// default access control CSR used alongside the maco CSRs to define
/// cachable/idempotent regions.
template <typename URV>
void
defineMdacSideEffects(System<URV>& system)
{
  for (unsigned i = 0; i < system.hartCount(); ++i)
    {
      auto hart = system.ithHart(i);
      auto csrPtr = hart->findCsr("mdac");
      if (not csrPtr)
        continue;

      auto pre = [] (Csr<URV>&, URV& val) -> void {
                   // A value of 0b11 (io/cacheable) is invalid: Make it 0b10
                   // (io/non-cacheable).
                   URV mask = 0b11;
                   if ((val & mask) == mask)
                     val = (val & ~mask) | 0b10;
                 };

      // Mdac is not shared between harts.
      auto post = [hart] (Csr<URV>&, URV val) -> void {
                    bool idempotent = (val & 2) == 0;
                    bool cacheable = (val & 1);
                    hart->setDefaultIdempotent(idempotent);
                    hart->setDefaultCacheable(cacheable);
                  };

      auto reset = [hart] (Csr<URV>& csr) -> void {
          URV val = csr.read();
          bool idempotent = (val & 2) == 0;
          bool cacheable = (val & 1);
          hart->setDefaultIdempotent(idempotent);
          hart->setDefaultCacheable(cacheable);
        };

      csrPtr->registerPrePoke(pre);
      csrPtr->registerPreWrite(pre);

      csrPtr->registerPostPoke(post);
      csrPtr->registerPostWrite(post);

      csrPtr->registerPostReset(reset);
    }
}


/// Associate callbacks with write/poke of the mpicbaddr (pic base
/// address csr).
template <typename URV>
void
defineMpicbaddrSideEffects(System<URV>& system)
{
  for (unsigned i = 0; i < system.hartCount(); ++i)
    {
      auto hart = system.ithHart(i);
      auto csrPtr = hart->findCsr("mpicbaddr");
      if (not csrPtr)
        continue;

      auto post = [hart] (Csr<URV>&, URV val) -> void {
                    hart->changeMemMappedBase(val);
                  };

      csrPtr->registerPostPoke(post);
      csrPtr->registerPostWrite(post);
    }
}


/// Associate callbacks with write/poke of mhartstart to start harts
/// when corresponding bits are set in that CSR.
template <typename URV>
void
defineMhartstartSideEffects(System<URV>& system)
{
  for (unsigned i = 0; i < system.hartCount(); ++i)
    {
      auto hart = system.ithHart(i);
      auto csrPtr = hart->findCsr("mhartstart");
      if (not csrPtr)
        continue;
      auto csrNum = csrPtr->getNumber();

      auto post = [&system] (Csr<URV>&, URV val) -> void {
                    // Start harts corresponding to set bits
                    for (unsigned ii = 0; ii < system.hartCount(); ++ii)
                      {
                        auto ht = system.ithHart(ii);
                        URV id = ht->sysHartIndex();
                        if (val & (URV(1) << id))
                          ht->setStarted(true);
                      }
                  };

      auto pre = [&system, csrNum] (Csr<URV>&, URV& val) -> void {
                   // Implement write-one semantics. We let hart 0 do
                   // the shared CSR value change.
                   auto ht = system.ithHart(0);
                   URV prev = 0;
                   ht->peekCsr(csrNum, prev);
                   val |= prev;
                 };

      csrPtr->registerPostPoke(post);
      csrPtr->registerPostWrite(post);
      csrPtr->registerPrePoke(pre);
      csrPtr->registerPreWrite(pre);
    }
}


/// Associate callback with write/poke of mnmipdel to deletage
/// non-maskable-interrupts to harts.
template <typename URV>
void
defineMnmipdelSideEffects(System<URV>& system)
{
  for (unsigned ix = 0; ix < system.hartCount(); ++ix)
    {
      auto hart = system.ithHart(ix);
      auto csrPtr = hart->findCsr("mnmipdel");
      if (not csrPtr)
        continue;

      // Enable NMI for harts corresponding to set bits in mnmipdel.
      auto post = [&system] (Csr<URV>& csr, URV val) -> void {
                    if ((val & csr.getWriteMask()) == 0)
                      return;
                    for (unsigned i = 0; i < system.hartCount(); ++i)
                      {
                        auto ht = system.ithHart(i);
                        bool enable = (val & (URV(1) << i)) != 0;
                        ht->enableNmi(enable);
                      }
                  };

      // If an attempt to change writeable bits to all-zero, keep
      // previous value.
      auto pre = [] (Csr<URV>& csr, URV& val) -> void {
                   URV prev = csr.read();
                   if ((val & csr.getWriteMask()) == 0)
                     val = prev;
                 };

      // On reset, enable NMI in harts according to the bits of mnmipdel
      auto reset = [hart] (Csr<URV>& csr) -> void {
                   URV val = csr.read();
                   URV id = hart->sysHartIndex();
                   bool flag = (val & (URV(1) << id)) != 0;
                   hart->enableNmi(flag);
                 };

      csrPtr->registerPostPoke(post);
      csrPtr->registerPostWrite(post);

      csrPtr->registerPrePoke(pre);
      csrPtr->registerPreWrite(pre);

      csrPtr->registerPostReset(reset);
    }
}


/// Associate callback with write/poke of mpmpc
template <typename URV>
void
defineMpmcSideEffects(System<URV>& system)
{
  for (unsigned i = 0; i < system.hartCount(); ++i)
    {
      auto hart = system.ithHart(i);
      auto csrPtr = hart->findCsr("mpmc");
      if (not csrPtr)
        continue;

      // Writing 3 to pmpc enables external interrupts unless in debug mode.
      auto prePoke = [hart] (Csr<URV>& csr, URV& val) -> void {
                       if (hart->inDebugMode() or (val & 3) != 3 or
                           (val & csr.getPokeMask()) == 0)
                         return;
                       URV mval = 0;
                       if (not hart->peekCsr(CsrNumber::MSTATUS, mval))
                         return;
                       MstatusFields<URV> fields(mval);
                       fields.bits_.MIE = 1;
                       hart->pokeCsr(CsrNumber::MSTATUS, fields.value_);
                     };

      auto preWrite = [hart] (Csr<URV>& csr, URV& val) -> void {
                        if (hart->inDebugMode() or (val & 3) != 3 or
                           (val & csr.getWriteMask()) == 0)
                          return;
                        URV mval = 0;
                        if (not hart->peekCsr(CsrNumber::MSTATUS, mval))
                          return;
                       MstatusFields<URV> fields(mval);
                       fields.bits_.MIE = 1;
                       hart->pokeCsr(CsrNumber::MSTATUS, fields.value_);
                       hart->recordCsrWrite(CsrNumber::MSTATUS);
                     };


      csrPtr->registerPrePoke(prePoke);
      csrPtr->registerPreWrite(preWrite);
    }
}


// Define callback to react to write/poke to mgpmc CSR. This is a
// mon-standard WD CSR that is only in ehx1 and that controls the
// performance counters. Mgpmc will not be present if mcountinhibit is
// present.
template <typename URV>
void
defineMgpmcSideEffects(System<URV>& system)
{
  for (unsigned i = 0; i < system.hartCount(); ++i)
    {
      auto hart = system.ithHart(i);
      auto csrPtr = hart->findCsr("mgpmc");
      if (not csrPtr)
        continue;

      // For poke, the effect takes place immediately (next instruction
      // will see the new control).
      auto postPoke = [hart] (Csr<URV>&, URV val) -> void {
                        bool enable = (val & 1) == 1;
                        URV mask = enable? ~URV(0) : 0;
                        mask |= 7; // cycle/time/instret not controlled by mgpmc
                        hart->setPerformanceCounterControl(mask);
                        hart->setPerformanceCounterControl(mask);
                      };

      // For write (invoked from current instruction), the effect
      // takes place on the following instruction.
      auto postWrite = [hart] (Csr<URV>&, URV val) -> void {
                        bool enable = (val & 1) == 1;
                        URV mask = enable? ~URV(0) : 0;
                        mask |= 7; // cycle/time/instret not controlled by mgpmc
                        hart->setPerformanceCounterControl(mask);
                       };

      csrPtr->registerPostPoke(postPoke);
      csrPtr->registerPostWrite(postWrite);
    }
}


/// Associate callback with write/poke of mcounthinibit
template <typename URV>
void
defineMcountinhibitSideEffects(System<URV>& system)
{
  for (unsigned i = 0; i < system.hartCount(); ++i)
    {
      auto hart = system.ithHart(i);
      auto csrPtr = hart->findCsr("mcountinhibit");
      if (not csrPtr)
        continue;

      // For poke, the effect takes place immediately (next instruction
      // will see the new control).
      auto postPoke = [hart] (Csr<URV>&, URV val) -> void {
                        hart->setPerformanceCounterControl(~val);
                        hart->setPerformanceCounterControl(~val);
                      };

      // For write (invoked from current instruction), the effect
      // takes place on the following instruction.
      auto postWrite = [hart] (Csr<URV>&, URV val) -> void {
                         hart->setPerformanceCounterControl(~val);
                       };

      csrPtr->registerPostPoke(postPoke);
      csrPtr->registerPostWrite(postWrite);
    }
}


template <typename URV>
bool
HartConfig::finalizeCsrConfig(System<URV>& system) const
{
  if (system.hartCount() == 0)
    return false;

  // Make shared CSRs in each hart except first one in core point to
  // the corresponding values in the first hart zero.
  for (unsigned ci = 0; ci < system.coreCount(); ++ci)
    {
      auto corePtr = system.ithCore(ci);
      if (not corePtr)
        continue;

      auto hart0 = corePtr->ithHart(0);
      if (not hart0)
        continue;

      for (unsigned hi = 1; hi < corePtr->hartCount(); ++hi)
        {
          auto hartPtr = corePtr->ithHart(hi);
          if (hartPtr)
            hartPtr->tieSharedCsrsTo(*hart0);
        }
    }

  // The following are WD non-standard CSRs. We implement their
  // actions by associating callbacks with the write/poke CSR methods.
  defineMhartstartSideEffects(system);
  defineMnmipdelSideEffects(system);
  defineMpicbaddrSideEffects(system);
  defineMpmcSideEffects(system);
  defineMgpmcSideEffects(system);
  defineMacoSideEffects(system);
  defineMracSideEffects(system);
  defineMdacSideEffects(system);

  // Define callback to react to write/poke to mcountinhibit CSR.
  defineMcountinhibitSideEffects(system);

  return true;
}


bool
HartConfig::getMisaReset(uint64_t& val) const
{
  val = 0;

  if (not config_ -> count("csr"))
    return false;  // No csr section

  const auto& csrs = config_ -> at("csr");
  if (not csrs.is_object())
    return false;  // No csr section in this config.

  if (not csrs.count("misa"))
    return false;  // CSR misa not present in csr section

  const auto& misa = csrs.at("misa");
  if (not misa.is_object())
    return false;

  if (not misa.count("reset"))
    return false;  // No reset entry under misa

  uint64_t resetVal = 0;
  if (not getJsonUnsigned("csr.misa.reset", misa.at("reset"), resetVal))
    return false;

  val = resetVal;
  return true;
}


// Instantiate tempate member functions

template bool
HartConfig::applyConfig<uint32_t>(Hart<uint32_t>&, bool, bool) const;

template bool
HartConfig::applyConfig<uint64_t>(Hart<uint64_t>&, bool, bool) const;

template bool
HartConfig::configHarts<uint32_t>(System<uint32_t>&, bool, bool) const;

template bool
HartConfig::configHarts<uint64_t>(System<uint64_t>&, bool, bool) const;

template bool
HartConfig::configMemory(System<uint32_t>&, bool, bool, bool) const;

template bool
HartConfig::configMemory(System<uint64_t>&, bool, bool, bool) const;


template bool
HartConfig::applyMemoryConfig<uint32_t>(Hart<uint32_t>&, bool, bool) const;

template bool
HartConfig::applyMemoryConfig<uint64_t>(Hart<uint64_t>&, bool, bool) const;

template bool
HartConfig::finalizeCsrConfig<uint32_t>(System<uint32_t>&) const;

template bool
HartConfig::finalizeCsrConfig<uint64_t>(System<uint64_t>&) const;


template
void
unpackMacoValue<uint32_t>(uint32_t value, uint32_t mask, bool rv32,
                          uint64_t& start, uint64_t& end, bool& idempotent,
                          bool& cacheable);

template
void
unpackMacoValue<uint64_t>(uint64_t value, uint64_t mask, bool rv32,
                          uint64_t& start, uint64_t& end, bool& idempotent,
                          bool& cacheable);
