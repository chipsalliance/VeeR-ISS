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
#include <cstdio>
#include "SparseMem.hpp"


using namespace WdRiscv;


SparseMem::SparseMem()
{
}


SparseMem::~SparseMem()
{
  for (auto kv : pageMap_)
    delete [] kv.second;
}


bool
SparseMem::read(uint64_t addr, unsigned size, uint64_t& value)
{
  if (size == 1)
    return read<uint8_t>(addr, value);
  bool aligned = (addr & (size-1)) == 0;

  if (size == 2)
    {
      if (aligned)
        return read<uint16_t>(addr, value);
      uint64_t low, high;
      if (read<uint8_t>(addr, low) and read<uint8_t>(addr+1, high))
        {
          value = (high << 8) | low;
          return true;
        }
      return false;
    }

  if (size == 4)
    {
      if (aligned)
        return read<uint32_t>(addr, value);
      uint64_t a0, a1, a2, a3;
      if (read<uint8_t>(addr, a0) and read<uint8_t>(addr+1, a1) and
          read<uint8_t>(addr+2, a2) and read<uint8_t>(addr+3, a3))
        {
          value = (a3 << 24) | (a2 << 16) | (a1 << 8) | a0;
          return true;
        }
      return false;
    }

  if (size == 8)
    {
      if (aligned)
        return read<uint64_t>(addr, value);
      uint64_t a0, a1, a2, a3, a4, a5, a6, a7;
      if (read<uint8_t>(addr, a0) and read<uint8_t>(addr+1, a1) and
          read<uint8_t>(addr+2, a2) and read<uint8_t>(addr+3, a3) and
          read<uint8_t>(addr+4, a4) and read<uint8_t>(addr+5, a5) and
          read<uint8_t>(addr+6, a6) and read<uint8_t>(addr+7, a7))
        {
          value = ( (a7 << 56) | (a6 << 48) | (a5 << 40) | (a4 << 32) |
                    (a3 << 24) | (a2 << 16) | (a1 << 8) | a0 );
          return true;
        }
      return false;
    }

  return false;
}


bool
SparseMem::write(uint64_t addr, unsigned size, uint64_t val)
{
  if (size == 1)
    return write<uint8_t>(addr, val);
  bool aligned = (addr & (size-1)) == 0;

  if (size == 2)
    {
      if (aligned)
        return write<uint16_t>(addr, val);
      return write<uint8_t>(addr, val) and write<uint8_t>(addr+1, val >> 8);
    }

  if (size == 4)
    {
      if (aligned)
        return write<uint32_t>(addr, val);

      return (write<uint8_t>(addr, val) and write<uint8_t>(addr+1, val >> 8) and
              write<uint8_t>(addr+2, val >> 16) and write<uint8_t>(addr+3, val >> 24));
    }

  if (size == 8)
    {
      if (aligned)
        return write<uint64_t>(addr, val);

      return (write<uint8_t>(addr,   val)       and
              write<uint8_t>(addr+1, val >> 8)  and
              write<uint8_t>(addr+2, val >> 16) and
              write<uint8_t>(addr+3, val >> 24) and
              write<uint8_t>(addr+4, val >> 32) and
              write<uint8_t>(addr+5, val >> 40) and
              write<uint8_t>(addr+6, val >> 48) and
              write<uint8_t>(addr+7, val >> 56));
    }

  return false;
}


bool
SparseMem::writeHexFile(const std::string& path) const
{
  FILE* out = fopen(path.c_str(), "w");
  if (not out)
    {
      std::cerr << "SparseMem::writeHexFile failed - cannot open "
                << path << " for write\n";
      return false;
    }

  bool ok = true;

  for (const auto& kv : pageMap_)
    {
      uint64_t addr = kv.first * pageSize_;    // Page address
      uint8_t* data = kv.second;   // Page data
      if (fprintf(out, "@%0lx\n", addr) < 0)
        {
          ok = false;
          break;
        }

      size_t remain = pageSize_;
      while (remain and ok)
        {
          size_t chunk = std::min(remain, size_t(16));
          const char* sep = "";
          for (size_t i = 0; i < chunk; ++i)
            {
              if (fprintf(out, "%s%02x", sep, *data++) < 0)
                ok = false;
              sep = " ";
            }
          if (fprintf(out, "\n") < 0)
            ok = false;
          remain -= chunk;
        }
    }

  fclose(out);

  return ok;
}
