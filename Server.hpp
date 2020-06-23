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

#include "System.hpp"


namespace WdRiscv
{

  /// Manage server mode.
  template <typename URV>
  class Server
  {
  public:

    /// Constructor.
    Server(System<URV>& system);

    /// Server mode poke command.
    bool pokeCommand(const WhisperMessage& req, WhisperMessage& reply);

    /// Server mode peek command.
    bool peekCommand(const WhisperMessage& req, WhisperMessage& reply);

    // Server mode disassemble command.
    void disassembleAnnotateInst(Hart<URV>& hart,
                                 uint32_t inst, bool interrupted,
				 bool hasPreTrigger, bool hasPostTrigger,
				 std::string& text);

    /// Server mode step command.
    bool stepCommand(const WhisperMessage& req, 
		     std::vector<WhisperMessage>& pendingChanges,
		     WhisperMessage& reply,
		     FILE* traceFile);

    /// Server mode exception command.
    bool exceptionCommand(const WhisperMessage& req, WhisperMessage& reply,
			  std::string& text);

    /// Server mode loop: Receive command and send reply till a quit
    /// command is received. Return true on successful termination (quit
    /// received). Return false otherwise.
    bool interact(int soc, FILE* traceFile, FILE* commandLog);

  protected:

    /// Process changes of a single-step command. Put the changes in the
    /// pendingChanges vector (which is cleared on entry). Put the
    /// number of change record in the reply parameter along with the
    /// instruction address, opcode and assembly text. Use hasPre
    /// (instruction tripped a "before" trigger), hasPost (tripped an
    /// "after" trigger) and interrupted (instruction encountered an
    /// external interrupt) to annotate the assembly text.
    void processStepCahnges(Hart<URV>&, uint32_t inst,
			    std::vector<WhisperMessage>& pendingChanges,
			    bool interrupted, bool hasPre, bool hasPost,
			    WhisperMessage& reply);

  private:

    /// Check if target hart id is valid. Return true if it is, and
    /// false otherwise setting reply to invalid.
    bool checkHartId(const WhisperMessage& reg, WhisperMessage& reply);

    /// Check if target hart is valid and is started. Return true if
    /// it is, and false otherwise setting reply to invalid.  Complain
    /// about command receiven in non-started state.
    bool checkHart(const WhisperMessage& reg, const std::string& command,
                   WhisperMessage& reply);

    System<URV>& system_;
  };

}
