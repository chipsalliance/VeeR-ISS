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
#include <fstream>
#include <sstream>
#include <thread>
#include <atomic>
#if defined(__cpp_lib_filesystem)
  #include <filesystem>
  namespace FileSystem = std::filesystem;
#else
  #include <experimental/filesystem>
  namespace FileSystem = std::experimental::filesystem;
#endif

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#ifdef __MINGW64__
#include <winsock2.h>
typedef int socklen_t;
#define close(s)          closesocket((s))
#define setlinebuf(f)     setvbuf((f),NULL,_IOLBF,0)
#define strerror_r(a,b,c) strerror((a))
#else
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#endif

#include <csignal>
#include "HartConfig.hpp"
#include "WhisperMessage.h"
#include "Hart.hpp"
#include "Core.hpp"
#include "System.hpp"
#include "Server.hpp"
#include "Interactive.hpp"


using namespace WdRiscv;


/// Return format string suitable for printing an integer of type URV
/// in hexadecimal form.
template <typename URV>
static
const char*
getHexForm()
{
  if (sizeof(URV) == 4)
    return "0x%08x";
  if (sizeof(URV) == 8)
    return "0x%016x";
  if (sizeof(URV) == 16)
    return "0x%032x";
  return "0x%x";
}


/// Convert the command line string numberStr to a number using
/// strotull and a base of zero (prefixes 0 and 0x are
/// honored). Return true on success and false on failure (string does
/// not represent a number). TYPE is an integer type (e.g
/// uint32_t). Option is the command line option associated with the
/// string and is used for diagnostic messages.
template <typename TYPE>
static
bool
parseCmdLineNumber(const std::string& option,
		   const std::string& numberStr,
		   TYPE& number)
{
  std::string str = numberStr;
  bool good = not str.empty();
  uint64_t scale = 1;
  if (good)
    {
      char suffix = str.back();
      if (suffix == 'k')
        scale = 1024;
      else if (suffix == 'm')
        scale = 1024*1024;
      else if (suffix == 'g')
        scale = 1024*1024*1024;
      if (scale != 1)
        {
          str = str.substr(0, str.length() - 1);
          if (str.empty())
            good = false;
        }
    }

  if (good)
    {
      typedef typename std::make_signed_t<TYPE> STYPE;

      char* end = nullptr;
      
      bool bad = false;

      if (std::is_same<TYPE, STYPE>::value)
        {
          int64_t val = strtoll(str.c_str(), &end, 0) * scale;
          number = static_cast<TYPE>(val);
          bad = val != number;
        }
      else
        {
          uint64_t val = strtoull(str.c_str(), &end, 0) * scale;
          number = static_cast<TYPE>(val);
          bad = val != number;
        }

      if (bad)
	{
	  std::cerr << "parseCmdLineNumber: Number too large: " << numberStr
		    << '\n';
	  return false;
	}
      if (end and *end)
	good = false;  // Part of the string are non parseable.
    }

  if (not good)
    std::cerr << "Invalid command line " << option << " value: " << numberStr
	      << '\n';
  return good;
}


/// Aapter for the parseCmdLineNumber for optionals.
template <typename TYPE>
static
bool
parseCmdLineNumber(const std::string& option,
		   const std::string& numberStr,
		   std::optional<TYPE>& number)
{
  TYPE n;
  if (not parseCmdLineNumber(option, numberStr, n))
    return false;
  number = n;
  return true;
}


typedef std::vector<std::string> StringVec;


/// Hold values provided on the command line.
struct Args
{
  StringVec   hexFiles;        // Hex files to be loaded into simulator memory.
  std::string traceFile;       // Log of state change after each instruction.
  std::string commandLogFile;  // Log of interactive or socket commands.
  std::string consoleOutFile;  // Console io output file.
  std::string serverFile;      // File in which to write server host and port.
  std::string instFreqFile;    // Instruction frequency file.
  std::string configFile;      // Configuration (JSON) file.
  std::string bblockFile;      // Basci block file.
  std::string isa;
  std::string snapshotDir = "snapshot"; // Dir prefix for saving snapshots
  std::string loadFrom;        // Directory for loading a snapshot
  std::string stdoutFile;      // Redirect target program stdout to this.
  std::string stderrFile;      // Redirect target program stderr to this. 
  StringVec   zisa;
  StringVec   regInits;        // Initial values of regs
  StringVec   targets;         // Target (ELF file) programs and associated
                               // program options to be loaded into simulator
                               // memory. Each target plus args is one string.
  StringVec   isaVec;          // Isa string split around _ with rv32/rv64 prefix removed.
  std::string targetSep = " "; // Target program argument separator.

  std::optional<std::string> toHostSym;
  std::optional<std::string> consoleIoSym;

  // Ith item is a vector of strings representing ith target and its args.
  std::vector<StringVec> expandedTargets;

  std::optional<uint64_t> startPc;
  std::optional<uint64_t> endPc;
  std::optional<uint64_t> toHost;
  std::optional<uint64_t> consoleIo;
  std::optional<uint64_t> instCountLim;
  std::optional<uint64_t> memorySize;
  std::optional<uint64_t> snapshotPeriod;
  std::optional<uint64_t> alarmInterval;
  std::optional<uint64_t> swInterrupt;  // Sotware interrupt mem mapped address
  std::optional<uint64_t> clint;  // Clint mem mapped address
  std::optional<uint64_t> syscallSlam;

  unsigned regWidth = 32;
  unsigned harts = 1;
  unsigned cores = 1;
  unsigned pageSize = 4*1024;
  uint64_t bblockInsts = ~uint64_t(0);

  bool help = false;
  bool hasRegWidth = false;
  bool hasHarts = false;
  bool hasCores = false;
  bool trace = false;
  bool interactive = false;
  bool verbose = false;
  bool version = false;
  bool traceLdSt = false;  // Trace ld/st data address if true.
  bool csv = false;        // Log files in CSV format when true.
  bool triggers = false;   // Enable debug triggers when true.
  bool counters = false;   // Enable performance counters when true.
  bool gdb = false;        // Enable gdb mode when true.
  std::vector<unsigned> gdbTcpPort;        // Enable gdb mode over TCP when port is positive.
  bool abiNames = false;   // Use ABI register names in inst disassembly.
  bool newlib = false;     // True if target program linked with newlib.
  bool linux = false;      // True if target program linked with Linux C-lib.
  bool raw = false;        // True if bare-metal program (no linux no newlib).
  bool elfisa = false;     // Use ELF file RISCV architecture tags to set MISA if true.
  bool fastExt = false;    // True if fast external interrupt dispatch enabled.
  bool unmappedElfOk = false;
  bool iccmRw = false;
  bool quitOnAnyHart = false;    // True if run quits when any hart finishes.
  bool noConInput = false;       // If true console io address is not used for input (ld).
  bool relativeInstCount = false;

  // Expand each target program string into program name and args.
  void expandTargets();
};


void
Args::expandTargets()
{
  this->expandedTargets.clear();
  for (const auto& target : this->targets)
    {
      StringVec tokens;
      boost::split(tokens, target, boost::is_any_of(this->targetSep),
		   boost::token_compress_on);
      this->expandedTargets.push_back(tokens);
    }
}


static
void
printVersion()
{
  unsigned version = 1;
  unsigned subversion = 739;
  std::cout << "Version " << version << "." << subversion << " compiled on "
	    << __DATE__ << " at " << __TIME__ << '\n';
}


static
bool
collectCommandLineValues(const boost::program_options::variables_map& varMap,
			 Args& args)
{
  bool ok = true;

  if (varMap.count("startpc"))
    {
      auto numStr = varMap["startpc"].as<std::string>();
      if (not parseCmdLineNumber("startpc", numStr, args.startPc))
	ok = false;
    }

  if (varMap.count("endpc"))
    {
      auto numStr = varMap["endpc"].as<std::string>();
      if (not parseCmdLineNumber("endpc", numStr, args.endPc))
	ok = false;
    }

  if (varMap.count("tohost"))
    {
      auto numStr = varMap["tohost"].as<std::string>();
      if (not parseCmdLineNumber("tohost", numStr, args.toHost))
	ok = false;
    }

  if (varMap.count("consoleio"))
    {
      auto numStr = varMap["consoleio"].as<std::string>();
      if (not parseCmdLineNumber("consoleio", numStr, args.consoleIo))
	ok = false;
    }

  if (varMap.count("maxinst"))
    {
      auto numStr = varMap["maxinst"].as<std::string>();
      if (not parseCmdLineNumber("maxinst", numStr, args.instCountLim))
	ok = false;
      args.relativeInstCount = not numStr.empty() and numStr.at(0) == '+';
    }

  if (varMap.count("memorysize"))
    {
      auto numStr = varMap["memorysize"].as<std::string>();
      if (not parseCmdLineNumber("memorysize", numStr, args.memorySize))
        ok = false;
    }

  if (varMap.count("snapshotperiod"))
    {
      auto numStr = varMap["snapshotperiod"].as<std::string>();
      if (not parseCmdLineNumber("snapshotperiod", numStr, args.snapshotPeriod))
        ok = false;
      else if (*args.snapshotPeriod == 0)
        std::cerr << "Warning: Zero snapshot period ignored.\n";
    }

  if (varMap.count("tohostsym"))
    args.toHostSym = varMap["tohostsym"].as<std::string>();

  if (varMap.count("consoleiosym"))
    args.consoleIoSym = varMap["consoleiosym"].as<std::string>();

  if (varMap.count("xlen"))
    args.hasRegWidth = true;

  if (varMap.count("cores"))
    args.hasCores = true;

  if (varMap.count("harts"))
    args.hasHarts = true;

  if (varMap.count("alarm"))
    {
      auto numStr = varMap["alarm"].as<std::string>();
      if (not parseCmdLineNumber("alarm", numStr, args.alarmInterval))
        ok = false;
      else if (*args.alarmInterval == 0)
        std::cerr << "Warning: Zero alarm period ignored.\n";
    }

  if (varMap.count("clint"))
    {
      auto numStr = varMap["clint"].as<std::string>();
      if (not parseCmdLineNumber("clint", numStr, args.clint))
        ok = false;
      else if ((*args.clint & 7) != 0)
        {
          std::cerr << "Error: clint address must be a multiple of 8\n";
          ok = false;
        }
    }

  if (varMap.count("softinterrupt"))
    {
      auto numStr = varMap["softinterrupt"].as<std::string>();
      if (not parseCmdLineNumber("softinterrupt", numStr, args.swInterrupt))
        ok = false;
      else if ((*args.swInterrupt & 3) != 0)
        {
          std::cerr << "Error: softinterrupt address must be a multiple of 4\n";
          ok = false;
        }
    }

  if (varMap.count("syscallslam"))
    {
      auto numStr = varMap["syscallslam"].as<std::string>();
      if (not parseCmdLineNumber("syscallslam", numStr, args.syscallSlam))
        ok = false;
    }

  if (args.interactive)
    args.trace = true;  // Enable instruction tracing in interactive mode.

  return ok;
}


/// Parse command line arguments. Place option values in args.
/// Return true on success and false on failure. Exists program
/// if --help is used.
static
bool
parseCmdLineArgs(int argc, char* argv[], Args& args)
{
  try
    {
      // Define command line options.
      namespace po = boost::program_options;
      po::options_description desc("options");
      desc.add_options()
	("help,h", po::bool_switch(&args.help),
	 "Produce this message.")
	("log,l", po::bool_switch(&args.trace),
	 "Enable tracing to standard output of executed instructions.")
	("isa", po::value(&args.isa),
	 "Specify instruction set extensions to enable. Supported extensions "
	 "are a, c, d, f, i, m, s and u. Default is imc.")
	("zisa", po::value(&args.zisa)->multitoken(),
	 "Specify instruction set z-extension to enable. Only z-extensions "
	 "currently supported are zbb and zbs (Exammple --zisa zbb)")
	("xlen", po::value(&args.regWidth),
	 "Specify register width (32 or 64), defaults to 32")
	("harts", po::value(&args.harts),
	 "Specify number of hardware threads per core (default=1).")
	("cores", po::value(&args.cores),
	 "Specify number of core per system (default=1).")
	("pagesize", po::value(&args.pageSize),
	 "Specify memory page size.")
	("target,t", po::value(&args.targets)->multitoken(),
	 "Target program (ELF file) to load into simulator memory. In "
	 "newlib/Linux emulation mode, program options may follow program name.")
	("targetsep", po::value(&args.targetSep),
	 "Target program argument separator.")
	("hex,x", po::value(&args.hexFiles)->multitoken(),
	 "HEX file to load into simulator memory.")
	("logfile,f", po::value(&args.traceFile),
	 "Enable tracing to given file of executed instructions. Output is compressed (with /usr/bin/gzip) if file name ends with \".gz\".")
	("csvlog", po::bool_switch(&args.csv),
	 "Enable CSV format for log file.")
	("consoleoutfile", po::value(&args.consoleOutFile),
	 "Redirect console output to given file.")
	("commandlog", po::value(&args.commandLogFile),
	 "Enable logging of interactive/socket commands to the given file.")
	("server", po::value(&args.serverFile),
	 "Interactive server mode. Put server hostname and port in file.")
	("startpc,s", po::value<std::string>(),
	 "Set program entry point. If not specified, use entry point of the "
	 "most recently loaded ELF file.")
	("endpc,e", po::value<std::string>(),
	 "Set stop program counter. Simulator will stop once instruction at "
	 "the stop program counter is executed.")
	("tohost", po::value<std::string>(),
	 "Memory address to which a write stops simulator.")
	("tohostsym", po::value<std::string>(),
	 "ELF symbol to use for setting tohost from ELF file (in the case "
	 "where tohost is not specified on the command line). Default: "
	 "\"tohost\".")
	("consoleio", po::value<std::string>(),
	 "Memory address corresponding to console io. Reading/writing "
	 "(lw/lh/lb sw/sh/sb) from given address reads/writes a byte from the "
         "console.")
	("consoleiosym", po::value<std::string>(),
	 "ELF symbol to use as console-io address (in the case where "
         "consoleio is not specified on the command line). Deafult: "
         "\"__whisper_console_io\".")
	("maxinst,m", po::value<std::string>(),
	 "Limit executed instruction count to arg. With a leading plus sign interpret the count as relative to the loaded (from a snapshot) instruction count.")
	("memorysize", po::value<std::string>(),
	 "Memory size (must be a multiple of 4096).")
	("interactive,i", po::bool_switch(&args.interactive),
	 "Enable interactive mode.")
	("traceload", po::bool_switch(&args.traceLdSt),
	 "Enable tracing of load/store instruction data address.")
	("triggers", po::bool_switch(&args.triggers),
	 "Enable debug triggers (triggers are on in interactive and server modes)")
	("counters", po::bool_switch(&args.counters),
	 "Enable performance counters")
	("gdb", po::bool_switch(&args.gdb),
	 "Run in gdb mode enabling remote debugging from gdb (this requires gdb version"
         "8.2 or higher).")
	("gdb-tcp-port", po::value(&args.gdbTcpPort)->multitoken(),
	 	 "TCP port number for gdb; If port num is negative,"
			" gdb will work with stdio (default -1).")
	("profileinst", po::value(&args.instFreqFile),
	 "Report instruction frequency to file.")
	("setreg", po::value(&args.regInits)->multitoken(),
	 "Initialize registers. Apply to all harts unless specific prefix "
	 "present (hart is 1 in 1:x3=0xabc). Example: --setreg x1=4 x2=0xff "
	 "1:x3=0xabc")
	("configfile", po::value(&args.configFile),
	 "Configuration file (JSON file defining system features).")
	("bblockfile", po::value(&args.bblockFile),
	 "Basic blocks output stats file.")
	("bblockinterval", po::value(&args.bblockInsts),
	 "Basic block stats are reported even mulitples of given instruction counts and once at end of run.")
	("snapshotdir", po::value(&args.snapshotDir),
	 "Directory prefix for saving snapshots.")
	("snapshotperiod", po::value<std::string>(),
	 "Snapshot period: Save snapshot using snapshotdir every so many instructions.")
	("loadfrom", po::value(&args.loadFrom),
	 "Snapshot directory from which to restore a previously saved (snapshot) state.")
	("stdout", po::value(&args.stdoutFile),
	 "Redirect standard output of newlib/Linux target program to this.")
	("stderr", po::value(&args.stderrFile),
	 "Redirect standard error of newlib/Linux target program to this.")
	("abinames", po::bool_switch(&args.abiNames),
	 "Use ABI register names (e.g. sp instead of x2) in instruction disassembly.")
	("newlib", po::bool_switch(&args.newlib),
	 "Emulate (some) newlib system calls. Done automatically if newlib "
         "symbols are detected in the target ELF file.")
	("linux", po::bool_switch(&args.linux),
	 "Emulate (some) Linux system calls. Done automatically if Linux "
         "symbols are detected in the target ELF file.")
	("raw", po::bool_switch(&args.raw),
	 "Bare metal mode: Disble emulation of Linux/newlib system call emulation "
         "even if Linux/newlib symbols detected in the target ELF file.")
	("elfisa", po::bool_switch(&args.elfisa),
	 "Confiure reset value of MISA according to the RISCV architecture tag(s) "
         "encoded into the laoded ELF file(s) if any.")
	("fastext", po::bool_switch(&args.fastExt),
	 "Enable fast external interrupt dispatch.")
	("unmappedelfok", po::bool_switch(&args.unmappedElfOk),
	 "Do not flag as error ELF file sections targeting unmapped "
         " memory.")
	("alarm", po::value<std::string>(),
	 "External interrupt period in micro-seconds: Convert arg to an "
         "instruction count, n, assuming a 1ghz clock, and force an external "
         " interrupt every n instructions. No-op if arg is zero.")
        ("softinterrupt", po::value<std::string>(),
         "Address of memory mapped word(s) controlling software interrupts. In "
         "an n-hart system, words at addresses a, a+4, ... a+(n-1)*4 "
         "are associated with the n harts (\"a\" being the address "
         "specified by this option and must be a multiple of "
         "4). Writing 0/1 to one of these addresses (using sw) "
         "clear/sets the software interrupt bit in the the MIP (machine "
         "interrupt pending) CSR of the corresponding hart. If a "
         "software interrupt is taken, it is up to interrupt handler to "
         "write zero to the same location to clear the corresponding "
         "bit in MIP. Writing values besides 0/1 will not affect the "
         "MIP bit and neither will writing using sb/sh/sd or writing to "
         "non-multiple-of-4 addresses.")
        ("clint", po::value<std::string>(),
         "Define address, a, of memory mapped area for clint (core local "
         "interruptor). In an n-hart system, words at addresses a, a+4, ... "
         "a+(n-1)*4, are  associated with the n harts. Store a 0/1 to one of "
         "these locations clears/sets the software interrupt bit in the MIP CSR "
         "of the corresponding hart. Similary, addresses b, b+8, ... b+(n-1)*8, "
         "where b is a+0x4000, are associated with the n harts. Writing to one "
         "of these double words sets the timer-limit of the corresponding hart. "
         "A timer interrupt in such a hart becomes pending when the timer value "
         "equals or exceeds the timer limit.")
        ("syscallslam", po::value<std::string>(),
         "Define address, a, of a non-cached memory area in which the "
         "memory changes of an emulated system call will be slammed. This "
         "is used in server mode to relay the effects of a system call "
         "to the RTL simulator. The memory area at location a will be filled "
         "with a sequence of pairs of double words designating addresses and "
         "corresponding values. A zero/zero pair will indicate the end of "
         "sequence.")
        ("iccmrw", po::bool_switch(&args.iccmRw),
         "Temporary switch to make ICCM region available to ld/st isntructions.")
        ("quitany", po::bool_switch(&args.quitOnAnyHart),
         "Terminate multi-threaded run when any hart finishes (default is to wait "
         "for all harts.)")
        ("noconinput", po::bool_switch(&args.noConInput),
         "Do not use console IO address for input. Loads from the cosole io address "
         "simply return last value stored there.")
	("verbose,v", po::bool_switch(&args.verbose),
	 "Be verbose.")
	("version", po::bool_switch(&args.version),
	 "Print version.");

      // Define positional options.
      po::positional_options_description pdesc;
      pdesc.add("target", -1);

      // Parse command line options.
      po::variables_map varMap;
      po::command_line_parser parser(argc, argv);
      auto parsed = parser.options(desc).positional(pdesc).run();
      po::store(parsed, varMap);
      po::notify(varMap);

      // auto unparsed = po::collect_unrecognized(parsed.options, po::include_positional);

      if (args.version)
        printVersion();

      if (args.help)
	{
	  std::cout <<
	    "Simulate a RISCV system running the program specified by the given ELF\n"
	    "and/or HEX file. With --newlib/--linux, the ELF file is a newlib/linux linked\n"
	    "program and may be followed by corresponding command line arguments.\n"
	    "All numeric arguments are interpreted as hexadecimal numbers when prefixed"
	    " with 0x."
	    "Examples:\n"
	    "  whisper --target prog --log\n"
	    "  whisper --target prog --setreg sp=0xffffff00\n"
	    "  whisper --newlib --log --target \"prog -x -y\"\n"
	    "  whisper --linux --log --targetsep ':' --target \"prog:-x:-y\"\n\n";
	  std::cout << desc;
	  return true;
	}

      if (not collectCommandLineValues(varMap, args))
	return false;
    }

  catch (std::exception& exp)
    {
      std::cerr << "Failed to parse command line args: " << exp.what() << '\n';
      return false;
    }

  return true;
}


/// Apply register initializations specified on the command line.
template<typename URV>
static
bool
applyCmdLineRegInit(const Args& args, Hart<URV>& hart)
{
  bool ok = true;

  URV hartId = hart.sysHartIndex();

  for (const auto& regInit : args.regInits)
    {
      // Each register initialization is a string of the form reg=val
      // or hart:reg=val
      std::vector<std::string> tokens;
      boost::split(tokens, regInit, boost::is_any_of("="),
		   boost::token_compress_on);
      if (tokens.size() != 2)
	{
	  std::cerr << "Invalid command line register initialization: "
		    << regInit << '\n';
	  ok = false;
	  continue;
	}

      std::string regName = tokens.at(0);
      const std::string& regVal = tokens.at(1);

      bool specificHart = false;
      unsigned id = 0;
      size_t colonIx = regName.find(':');
      if (colonIx != std::string::npos)
	{
	  std::string hartStr = regName.substr(0, colonIx);
	  regName = regName.substr(colonIx + 1);
	  if (not parseCmdLineNumber("hart", hartStr, id))
	    {
	      std::cerr << "Invalid command line register initialization: "
			<< regInit << '\n';
	      ok = false;
	      continue;
	    }
	  specificHart = true;
	}

      URV val = 0;
      if (not parseCmdLineNumber("register", regVal, val))
	{
	  ok = false;
	  continue;
	}

      if (specificHart and id != hartId)
	continue;

      if (unsigned reg = 0; hart.findIntReg(regName, reg))
	{
	  if (args.verbose)
	    std::cerr << "Setting register " << regName << " to command line "
		      << "value 0x" << std::hex << val << std::dec << '\n';
	  hart.pokeIntReg(reg, val);
	  continue;
	}

      if (unsigned reg = 0; hart.findFpReg(regName, reg))
	{
	  if (args.verbose)
	    std::cerr << "Setting register " << regName << " to command line "
		      << "value 0x" << std::hex << val << std::dec << '\n';
	  hart.pokeFpReg(reg, val);
	  continue;
	}

      auto csr = hart.findCsr(regName);
      if (csr)
	{
	  if (args.verbose)
	    std::cerr << "Setting register " << regName << " to command line "
		      << "value 0x" << std::hex << val << std::dec << '\n';
	  hart.pokeCsr(csr->getNumber(), val);
	  continue;
	}

      std::cerr << "No such RISCV register: " << regName << '\n';
      ok = false;
    }

  return ok;
}


template<typename URV>
static
bool
applyZisaString(const std::string& zisa, Hart<URV>& hart)
{
  if (zisa.empty())
    return true;

  std::string ext = zisa;

  if (ext.at(0) == 'z')
    ext = ext.substr(1);

  if (boost::starts_with(ext, "ba"))
    hart.enableRvzba(true);
  else if (boost::starts_with(ext, "bb"))
    hart.enableRvzbb(true);
  else if (boost::starts_with(ext, "bc"))
    hart.enableRvzbc(true);
  else if (boost::starts_with(ext, "be"))
    hart.enableRvzbe(true);
  else if (boost::starts_with(ext, "bf"))
    hart.enableRvzbf(true);
  else if (boost::starts_with(ext, "bm"))
    hart.enableRvzbm(true);
  else if (boost::starts_with(ext, "bp"))
    hart.enableRvzbp(true);
  else if (boost::starts_with(ext, "br"))
    hart.enableRvzbr(true);
  else if (boost::starts_with(ext, "bs"))
    hart.enableRvzbs(true);
  else if (boost::starts_with(ext, "bt"))
    hart.enableRvzbt(true);
  else if (boost::starts_with(ext, "bmini"))
    {
      hart.enableRvzbb(true);
      hart.enableRvzbs(true);
      std::cerr << "ISA option zbmini is deprecated. Using zbb and zbs.\n";
    }
  else if (boost::starts_with(ext, "fh"))
    hart.enableZfh(true);
  else
    {
      std::cerr << "No such Z extension: " << zisa << '\n';
      return false;
    }

  return true;
}


template<typename URV>
static
bool
applyZisaStrings(const StringVec& zisa, Hart<URV>& hart)
{
  unsigned errors = 0;

  for (const auto& ext : zisa)
    {
      if (not applyZisaString(ext, hart))
	errors++;
    }

  return errors == 0;
}


template<typename URV>
static
bool
applyIsaStrings(const StringVec& isaStrings, Hart<URV>& hart)
{
  if (isaStrings.empty())
    return true;

  URV isa = 0;
  unsigned errors = 0;

  for (const auto& isaStr : isaStrings)
    {
      if (isaStr.empty())
	continue;

      char c = isaStr.at(0);
      switch(c)
	{
	case 'a':
	case 'c':
	case 'd':
        case 'e':
	case 'f':
	case 'i':
	case 'm':
	case 's':
	case 'u':
        case 'v':
	  isa |= URV(1) << (c -  'a');
	  break;

	case 'b':
	  std::cerr << "Warning: Extension \"" << c << "\" is not supported.\n";
	  break;

        case 'g':  // Enable a, d, f, and m
          isa |= 0x1 | 0x8 | 0x20 | 0x1000;
          break;

	case 'z':
	  if (not applyZisaString(isaStr, hart))
	    errors++;
	  break;
	  
	default:
	  std::cerr << "Error: Extension \"" << c << "\" is not supported.\n";
	  errors++;
	  break;
	}
    }

  if (not (isa & (URV(1) << ('i' - 'a'))))
    {
      std::cerr << "Extension \"i\" implicitly enabled\n";
      isa |= URV(1) << ('i' -  'a');
    }

  if (isa & (URV(1) << ('d' - 'a')))
    if (not (isa & (URV(1) << ('f' - 'a'))))
      std::cerr << "Warning Extension \"d\" enabled without \"f\"\n";

  // Set the xlen bits: 1 for 32-bits and 2 for 64.
  URV xlen = sizeof(URV) == 4? 1 : 2;
  isa |= xlen << (8*sizeof(URV) - 2);

  bool resetMemoryMappedRegs = false;

  URV mask = 0, pokeMask = 0;
  bool implemented = true, isDebug = false, shared = true;
  if (not hart.configCsr("misa", implemented, isa, mask, pokeMask, isDebug,
                         shared))
    {
      std::cerr << "Failed to configure MISA CSR\n";
      errors++;
    }
  else
    hart.reset(resetMemoryMappedRegs); // Apply effects of new misa value.

  return errors == 0;
}


/// Enable linux or newlib based on the symbols in the ELF files.
/// Return true if either is enabled.
template<typename URV>
static
bool
enableNewlibOrLinuxFromElf(const Args& args, Hart<URV>& hart)
{
  bool newlib = args.newlib, linux = args.linux;
  if (args.raw)
    {
      if (newlib or linux)
	std::cerr << "Raw mode not comptible with newlib/linux. Sticking"
		  << " with raw mode.\n";
      return false;
    }

  if (linux or newlib)
    ;  // Emulation preference already set by user.
  else
    {
      // At this point ELF files have not been loaded: Cannot use
      // hart.findElfSymbol.
      for (auto target : args.expandedTargets)
	{
	  auto elfPath = target.at(0);
	  if (not linux)
	    linux = Memory::isSymbolInElfFile(elfPath, "__libc_csu_init");

	  if (not newlib)
	    newlib = Memory::isSymbolInElfFile(elfPath, "__call_exitprocs");
	}

      if (args.verbose and linux)
	std::cerr << "Deteced linux symbol in ELF\n";

      if (args.verbose and newlib)
	std::cerr << "Deteced newlib symbol in ELF\n";

      if (newlib and linux)
	{
	  std::cerr << "Fishy: Both newlib and linux symbols present in "
		    << "ELF file(s). Doing linux emulation.\n";
	  newlib = false;
	}
    }

  hart.enableNewlib(newlib);
  hart.enableLinux(linux);

  return newlib or linux;
}


/// Set stack pointer to a reasonable value for linux/newlib.
template<typename URV>
static
void
sanitizeStackPointer(Hart<URV>& hart, bool verbose)
{
  // Set stack pointer to the 128 bytes below end of memory.
  size_t memSize = hart.getMemorySize();
  if (memSize > 128)
    {
      size_t spValue = memSize - 128;
      if (verbose)
	std::cerr << "Setting stack pointer to 0x" << std::hex << spValue
		  << std::dec << " for newlib/linux\n";
      hart.pokeIntReg(IntRegNumber::RegSp, spValue);
    }
}


/// Load register and memory state from snapshot previously saved
/// in the given directory. Return true on success and false on
/// failure.
template <typename URV>
static
bool
loadSnapshot(Hart<URV>& hart, const std::string& snapDir)
{
  using std::cerr;

  if (not FileSystem::is_directory(snapDir))
    {
      cerr << "Error: Path is not a snapshot directory: " << snapDir << '\n';
      return false;
    }

  FileSystem::path path(snapDir);
  FileSystem::path regPath = path / "registers";
  if (not FileSystem::is_regular_file(regPath))
    {
      cerr << "Error: Snapshot file does not exists: " << regPath << '\n';
      return false;
    }

  FileSystem::path memPath = path / "memory";
  if (not FileSystem::is_regular_file(regPath))
    {
      cerr << "Error: Snapshot file does not exists: " << memPath << '\n';
      return false;
    }

  if (not hart.loadSnapshot(path))
    {
      cerr << "Error: Failed to load sanpshot from dir " << snapDir << '\n';
      return false;
    }

  return true;
}



template<typename URV>
static
void
configureClint(Hart<URV>& hart, System<URV>& system, uint64_t clintStart,
               uint64_t clintLimit, uint64_t timerAddr)
{
  // Define callback to associate a memory mapped software interrupt
  // location to its corresponding hart so that when such a location
  // is written the software interrupt bit is set/cleared in the MIP
  // register of that hart.
  uint64_t swAddr = clintStart;
  auto swAddrToHart = [swAddr, &system](URV addr) -> Hart<URV>* {
    uint64_t addr2 = swAddr + system.hartCount()*4; // 1 word per hart
    if (addr >= swAddr and addr < addr2)
      {
        size_t ix = (addr - swAddr) / 4;
        return system.ithHart(ix).get();
      }
    return nullptr;
  };

  // Same for timer limit addresses.
  auto timerAddrToHart = [timerAddr, &system](URV addr) -> Hart<URV>* {
    uint64_t addr2 = timerAddr + system.hartCount()*8; // 1 double word per hart
    if (addr >= timerAddr and addr < addr2)
      {
        size_t ix = (addr - timerAddr) / 8;
        return system.ithHart(ix).get();
      }
    return nullptr;
  };

  hart.configClint(clintStart, clintLimit, swAddrToHart, timerAddrToHart);
}


static
bool
getElfFilesIsaString(const Args& args, std::string& isaString)
{
  std::vector<std::string> archTags;

  unsigned errors = 0;
  
  for (const auto& target : args.expandedTargets)
    {
      const auto& elfFile = target.front();
      if (not Memory::collectElfRiscvTags(elfFile, archTags))
        errors++;
    }

  if (archTags.empty())
    return errors == 0;

  const std::string& ref = archTags.front();

  for (const auto& tag : archTags)
    if (tag != ref)
      std::cerr << "Warning differen ELF files have different ISA strings: "
		<< tag << " and " << ref << '\n';

  isaString = ref;

  if (args.verbose)
    std::cerr << "ISA string from ELF file(s): " << isaString << '\n';

  return errors == 0;
}


/// Return the string representing the current contents of the MISA CSR.
template<typename URV>
static
std::string
getIsaStringFromCsr(const Hart<URV>& hart)
{
  std::string res;

  URV val;
  if (not hart.peekCsr(CsrNumber::MISA, val))
    return res;

  URV mask = 1;
  for (char c = 'a'; c <= 'z'; ++c, mask <<= 1)
    if (val & mask)
      res += c;

  return res;
}

/// Apply command line arguments: Load ELF and HEX files, set
/// start/end/tohost. Return true on success and false on failure.
template<typename URV>
static
bool
applyCmdLineArgs(const Args& args, StringVec isaVec, Hart<URV>& hart, System<URV>& system)
{
  unsigned errors = 0;

  // Handle linux/newlib adjusting stack if needed.
  bool clib = enableNewlibOrLinuxFromElf(args, hart);

  if (clib and isaVec.empty() and not args.raw)
    {
      if (args.verbose)
        std::cerr << "Adding a/c/m/f/d extensions for newlib/linux\n";
      std::string isa = getIsaStringFromCsr(hart);
      for (auto c : std::string("icmafd"))
	if (isa.find(c) == std::string::npos)
	  isa += c;
      for (auto c : isa)
	isaVec.push_back(std::string(1, c));
    }

  if (not applyIsaStrings(isaVec, hart))
      errors++;

  if (not applyZisaStrings(args.zisa, hart))
    errors++;

  if (clib)  // Linux or newlib enabled.
    sanitizeStackPointer(hart, args.verbose);

  if (args.toHostSym)
    hart.setTohostSymbol(*args.toHostSym);

  if (args.consoleIoSym)
    hart.setConsoleIoSymbol(*args.consoleIoSym);

  // Load ELF files.
  for (const auto& target : args.expandedTargets)
    {
      const auto& elfFile = target.front();
      if (args.verbose)
	std::cerr << "Loading ELF file " << elfFile << '\n';
      size_t entryPoint = 0;
      if (hart.loadElfFile(elfFile, entryPoint))
	hart.pokePc(URV(entryPoint));
      else
	errors++;
    }

  // Load HEX files.
  for (const auto& hexFile : args.hexFiles)
    {
      if (args.verbose)
	std::cerr << "Loading HEX file " << hexFile << '\n';
      if (not hart.loadHexFile(hexFile))
	errors++;
    }

  if (not args.instFreqFile.empty())
    hart.enableInstructionFrequency(true);

  if (not args.loadFrom.empty())
    if (not loadSnapshot(hart, args.loadFrom))
      errors++;

  if (not args.stdoutFile.empty())
    if (not hart.redirectOutputDescriptor(STDOUT_FILENO, args.stdoutFile))
      errors++;

  if (not args.stderrFile.empty())
    if (not hart.redirectOutputDescriptor(STDERR_FILENO, args.stderrFile))
      errors++;

  // Command line to-host overrides that of ELF and config file.
  if (args.toHost)
    hart.setToHostAddress(*args.toHost);

  // Command-line entry point overrides that of ELF.
  if (args.startPc)
    hart.pokePc(URV(*args.startPc));

  // Command-line exit point overrides that of ELF.
  if (args.endPc)
    hart.setStopAddress(URV(*args.endPc));

  // Command-line console io address overrides config file.
  if (args.consoleIo)
    hart.setConsoleIo(URV(*args.consoleIo));

  hart.enableConsoleInput(! args.noConInput);

  if (args.clint and args.swInterrupt)
    std::cerr << "Ignoring --sontinterrupt: incompatible with --clint.\n";

  if (args.clint)
    {
      uint64_t swAddr = *args.clint;
      uint64_t timerAddr = swAddr + 0x4000;
      uint64_t clintLimit = swAddr + 0x40000000 - 1;
      configureClint(hart, system, swAddr, clintLimit, timerAddr);
    }
  else if (args.swInterrupt)
    {
      uint64_t swAddr = *args.swInterrupt;
      uint64_t timerAddr = swAddr + 0x4000;
      uint64_t clintLimit = swAddr + system.hartCount() * 4 - 1;
      configureClint(hart, system, swAddr, clintLimit, timerAddr);
    }

  if (args.syscallSlam)
    hart.defineSyscallSlam(*args.syscallSlam);

  // Set instruction count limit.
  if (args.instCountLim)
    {
      uint64_t count = args.relativeInstCount? hart.getInstructionCount() : 0;
      count += *args.instCountLim;
      hart.setInstructionCountLimit(count);
    }

  // Print load-instruction data-address when tracing instructions.
  if (args.traceLdSt)
    hart.setTraceLoadStore(args.traceLdSt);

  // Setup periodic external interrupts.
  if (args.alarmInterval)
    {
      // Convert from micro-seconds to processor ticks. Assume a 1
      // ghz-processor.
      uint64_t ticks = (*args.alarmInterval)*1000;
      hart.setupPeriodicTimerInterrupts(ticks);
    }

  if (args.triggers)
    hart.enableTriggers(args.triggers);
  hart.enableGdb(args.gdb);
  if (args.gdbTcpPort.size()>hart.sysHartIndex())
    hart.setGdbTcpPort(args.gdbTcpPort[hart.sysHartIndex()]);
  if (args.counters)
    hart.enablePerformanceCounters(args.counters);
  if (args.abiNames)
    hart.enableAbiNames(args.abiNames);

  if (args.fastExt)
    hart.enableFastInterrupts(args.fastExt);

  // Apply register initialization.
  if (not applyCmdLineRegInit(args, hart))
    errors++;

  if (args.expandedTargets.empty())
    return errors == 0;

  // Setup target program arguments.
  if (clib)
    {
      if (args.loadFrom.empty())
        if (not hart.setTargetProgramArgs(args.expandedTargets.front()))
          {
            size_t memSize = hart.memorySize();
            size_t suggestedStack = memSize - 4;

            std::cerr << "Failed to setup target program arguments -- stack "
                      << "is not writable\n"
                      << "Try using --setreg sp=<val> to set the stack pointer "
                      << "to a\nwritable region of memory (e.g. --setreg "
                      << "sp=0x" << std::hex << suggestedStack << '\n'
                      << std::dec;
            errors++;
          }
    }
  else if (args.expandedTargets.front().size() > 1)
    {
      std::cerr << "Warning: Target program options present which requires\n"
		<< "         the use of --newlib/--linux. Options ignored.\n";
    }

  if (args.csv)
    hart.enableCsvLog(args.csv);

  return errors == 0;
}


/// Open a server socket and put opened socket information (hostname
/// and port number) in the given server file. Wait for one
/// connection. Service connection. Return true on success and false
/// on failure.
template <typename URV>
static
bool
runServer(System<URV>& system, const std::string& serverFile,
	  FILE* traceFile, FILE* commandLog)
{
  char hostName[1024];
  if (gethostname(hostName, sizeof(hostName)) != 0)
    {
      std::cerr << "Failed to obtain name of this computer\n";
      return false;
    }

  int soc = socket(AF_INET, SOCK_STREAM, 0);
  if (soc < 0)
    {
      char buffer[512];
      char* p = buffer;
#ifdef __APPLE__
      strerror_r(errno, buffer, 512);
#else
      p = strerror_r(errno, buffer, 512);
#endif
      std::cerr << "Failed to create socket: " << p << '\n';
      return -1;
    }

  sockaddr_in serverAddr;
  memset(&serverAddr, 0, sizeof(serverAddr));
  serverAddr.sin_family = AF_INET;
  serverAddr.sin_addr.s_addr = htonl(INADDR_ANY);
  serverAddr.sin_port = htons(0);

  if (bind(soc, (sockaddr*) &serverAddr, sizeof(serverAddr)) < 0)
    {
      perror("Socket bind failed");
      return false;
    }

  if (listen(soc, 1) < 0)
    {
      perror("Socket listen failed");
      return false;
    }

  sockaddr_in socAddr;
  socklen_t socAddrSize = sizeof(socAddr);
  socAddr.sin_family = AF_INET;
  socAddr.sin_port = 0;
  if (getsockname(soc, (sockaddr*) &socAddr,  &socAddrSize) == -1)
    {
      perror("Failed to obtain socket information");
      return false;
    }

  {
    std::ofstream out(serverFile);
    if (not out.good())
      {
	std::cerr << "Failed to open file '" << serverFile << "' for output\n";
	return false;
      }
    out << hostName << ' ' << ntohs(socAddr.sin_port) << std::endl;
  }

  sockaddr_in clientAddr;
  socklen_t clientAddrSize = sizeof(clientAddr);
  int newSoc = accept(soc, (sockaddr*) & clientAddr, &clientAddrSize);
  if (newSoc < 0)
    {
      perror("Socket accept failed");
      return false;
    }

  bool ok = true;

  try
    {
      Server<URV> server(system);
      ok = server.interact(newSoc, traceFile, commandLog);
    }
  catch(...)
    {
      ok = false;
    }

  close(newSoc);
  close(soc);

  return ok;
}


template <typename URV>
static
bool
reportInstructionFrequency(Hart<URV>& hart, const std::string& outPath)
{
  FILE* outFile = fopen(outPath.c_str(), "w");
  if (not outFile)
    {
      std::cerr << "Failed to open instruction frequency file '" << outPath
		<< "' for output.\n";
      return false;
    }
  hart.reportInstructionFrequency(outFile);
  hart.reportTrapStat(outFile);
  fprintf(outFile, "\n");
  hart.reportPmpStat(outFile);
  fprintf(outFile, "\n");
  hart.reportLrScStat(outFile);

  fclose(outFile);
  return true;
}


/// Open the trace-file, command-log and console-output files
/// specified on the command line. Return true if successful or false
/// if any specified file fails to open.
static
bool
openUserFiles(const Args& args, FILE*& traceFile, FILE*& commandLog,
	      FILE*& consoleOut, FILE*& bblockFile)
{
  size_t len = args.traceFile.size();
  bool doGzip = len > 3 and args.traceFile.substr(len-3) == ".gz";

  if (not args.traceFile.empty())
    {
      if (doGzip)
	{
	  std::string cmd = "/usr/bin/gzip -c > ";
	  cmd += args.traceFile;
	  traceFile = popen(cmd.c_str(), "w");
	}
      else
	traceFile = fopen(args.traceFile.c_str(), "w");
      if (not traceFile)
	{
	  std::cerr << "Failed to open trace file '" << args.traceFile
		    << "' for output\n";
	  return false;
	}
    }

  if (args.trace and traceFile == NULL)
    traceFile = stdout;
  if (traceFile and not doGzip)
    setlinebuf(traceFile);  // Make line-buffered.

  if (not args.commandLogFile.empty())
    {
      commandLog = fopen(args.commandLogFile.c_str(), "w");
      if (not commandLog)
	{
	  std::cerr << "Failed to open command log file '"
		    << args.commandLogFile << "' for output\n";
	  return false;
	}
      setlinebuf(commandLog);  // Make line-buffered.
    }

  if (not args.consoleOutFile.empty())
    {
      consoleOut = fopen(args.consoleOutFile.c_str(), "w");
      if (not consoleOut)
	{
	  std::cerr << "Failed to open console output file '"
		    << args.consoleOutFile << "' for output\n";
	  return false;
	}
    }

  if (not args.bblockFile.empty())
    {
      bblockFile = fopen(args.bblockFile.c_str(), "w");
      if (not bblockFile)
	{
	  std::cerr << "Failed to open basic block file '"
		    << args.bblockFile << "' for output\n";
	  return false;
	}
    }

  return true;
}


/// Counterpart to openUserFiles: Close any open user file.
static
void
closeUserFiles(const Args& args, FILE*& traceFile, FILE*& commandLog,
	       FILE*& consoleOut, FILE*& bblockFile)
{
  if (consoleOut and consoleOut != stdout)
    fclose(consoleOut);
  consoleOut = nullptr;

  if (traceFile and traceFile != stdout)
    {
      size_t len = args.traceFile.size();
      bool doGzip = len > 3 and args.traceFile.substr(len-3) == ".gz";
      if (doGzip)
	pclose(traceFile);
      else
	fclose(traceFile);
    }
  traceFile = nullptr;

  if (commandLog and commandLog != stdout)
    fclose(commandLog);
  commandLog = nullptr;

  if (bblockFile and bblockFile != stdout)
    fclose(bblockFile);
  bblockFile = nullptr;
}


// In interactive mode, keyboard interrupts (typically control-c) are
// ignored.
static void
kbdInterruptHandler(int)
{
  std::cerr << "keyboard interrupt\n";
}


template <typename URV>
static bool
batchRun(System<URV>& system, FILE* traceFile, bool waitAll)
{
  if (system.hartCount() == 0)
    return true;

  if (system.hartCount() == 1)
    {
      auto& hart = *system.ithHart(0);
      bool ok = hart.run(traceFile);
#ifdef FAST_SLOPPY
      hart.reportOpenedFiles(std::cout);
#endif
      return ok;
    }

  // Run each hart in its own thread.

  std::vector<std::thread> threadVec;

  std::atomic<bool> result = true;
  std::atomic<unsigned> finished = 0;  // Count of finished threads. 
  std::atomic<bool> hart0Done = false;

  auto threadFunc = [&traceFile, &result, &finished, &hart0Done] (Hart<URV>* hart) {
                      // In multi-hart system, wait till hart is started by hart0.
                      while (not hart->isStarted())
                        if (hart0Done)
                          return;  // We are not going to be started.
		      bool r = hart->run(traceFile);
		      result = result and r;
                      finished++;
                      if (hart->sysHartIndex() == 0)
                        hart0Done = true;
		    };

  for (unsigned i = 0; i < system.hartCount(); ++i)
    {
      Hart<URV>* hart = system.ithHart(i).get();
      threadVec.emplace_back(std::thread(threadFunc, hart));
    }                             

  if (waitAll)
    {
      for (auto& t : threadVec)
        t.join();
    }
  else
    {
      // First thread to finish terminates run.
      while (finished == 0)
        ;

      extern void forceUserStop(int);
      forceUserStop(0);
      
      for (auto& t : threadVec)
        t.join();
    }

  return result;
}


/// Run producing a snapshot after each snapPeriod instructions. Each
/// snapshot goes into its own directory names <dir><n> where <dir> is
/// the string in snapDir and <n> is a sequential integer starting at
/// 0. Return true on success and false on failure.
template <typename URV>
static
bool
snapshotRun(System<URV>& system, FILE* traceFile,
            const std::string& snapDir, uint64_t snapPeriod)
{
  if (not snapPeriod)
    {
      std::cerr << "Warning: Zero snap period ignored.\n";
      return batchRun(system, traceFile, true /* waitAll */);
    }

  assert(system.hartCount() == 1);
  Hart<URV>& hart = *(system.ithHart(0));

  bool done = false;
  uint64_t globalLimit = hart.getInstructionCountLimit();

  while (not done)
    {
      uint64_t nextLimit = hart.getInstructionCount() +  snapPeriod;
      if (nextLimit >= globalLimit)
        done = true;
      nextLimit = std::min(nextLimit, globalLimit);
      hart.setInstructionCountLimit(nextLimit);
      hart.run(traceFile);
      if (hart.hasTargetProgramFinished())
        done = true;
      if (not done)
        {
          unsigned index = hart.snapshotIndex();
          FileSystem::path path(snapDir + std::to_string(index));
          if (not FileSystem::is_directory(path))
            if (not FileSystem::create_directories(path))
              {
                std::cerr << "Error: Failed to create snapshot directory " << path << '\n';
                return false;
              }
          hart.setSnapshotIndex(index + 1);
          if (not hart.saveSnapshot(path))
            {
              std::cerr << "Error: Failed to save a snapshot\n";
              return false;
            }
        }
    }

#ifdef FAST_SLOPPY
  hart.reportOpenedFiles(std::cout);
#endif

  return true;
}


static
bool
extractOptVersion(const std::string& isa, size_t& i, std::string& opt)
{
  while (i+1 < isa.size() and std::isdigit(isa.at(i+1)))
    opt.push_back(isa.at(++i));
  if (std::isdigit(opt.back()))
    {
      if (i+1 < isa.size() and isa.at(i+1) == 'p')
	{
	  opt.push_back(isa.at(++i));
	  while (i+1 < isa.size() and std::isdigit(isa.at(i+1)))
	    opt.push_back(isa.at(++i));
	  if (opt.back() == 'p')
	    return false;
	}
    }
  return true;
}


static
bool
determineIsa(const Args& args, StringVec& isaVec)
{
  isaVec.clear();

  if (not args.isa.empty() and args.elfisa)
    std::cerr << "Warning: Both --isa and --elfisa present: Using --isa\n";

  std::string isa = args.isa;

  if (isa.empty() and args.elfisa)
    if (not getElfFilesIsaString(args, isa))
      return false;

  if (isa.empty())
    return true;

  // Xlen part of isa is processed in determineRegisterWidth
  if (boost::starts_with(isa, "rv32") or boost::starts_with(isa, "rv64"))
    isa = isa.substr(4);
  else if (boost::starts_with(isa, "rv") and isa.size() > 2 and std::isdigit(isa.at(2)))
    {
      std::cerr << "Unsupported ISA: " << isa << '\n';
      return false;
    }

  bool hasZ = false, good = true;
  
  for (size_t i = 0; i < isa.size() and good; ++i)
    {
      char c = isa.at(i);
      if (c == '_')  { good = i > 0; continue; }
      if (c == 'z')
	{
	  if (i == 0)
	    good = false;
	  else if (hasZ and i > 0 and isa.at(i-1) != '_')
	    good = false;
	  else if (i+3 > isa.size())
	    good = false;
	  else if (std::isdigit(i+1) or std::isdigit(i+2))
	    good = false;
	  else
	    {
	      hasZ = true;
	      std::string opt = isa.substr(i, 3);
	      i += 2;
	      good = extractOptVersion(isa, i, opt);
	      if (good)
		isaVec.push_back(opt);
	    }
	}
      else if (c >= 'a' and c < 'z')
	{
	  if (hasZ)
	    good = false;
	  else
	    {
	      std::string opt = isa.substr(i, 1);
	      good = extractOptVersion(isa, i, opt);
	      if (good)
		isaVec.push_back(opt);
	    }
	}
      else
	good = false;
    }

  if (not good)
    std::cerr << "Invalid ISA: " << isa << '\n';

  return good;
}


/// Depending on command line args, start a server, run in interactive
/// mode, or initiate a batch run.
template <typename URV>
static
bool
sessionRun(System<URV>& system, const Args& args, FILE* traceFile, FILE* cmdLog)
{
  StringVec isaVec;
  if (not determineIsa(args, isaVec))
    return false;

  for (unsigned i = 0; i < system.hartCount(); ++i)
    if (not applyCmdLineArgs(args, isaVec, *system.ithHart(i), system))
      if (not args.interactive)
	return false;

  // In server/interactive modes: enable triggers and performance counters.
  bool serverMode = not args.serverFile.empty();
  if (serverMode or args.interactive)
    {
      for (unsigned i = 0; i < system.hartCount(); ++i)
        {
          auto& hart = *system.ithHart(i);
          hart.enableTriggers(true);
          hart.enablePerformanceCounters(true);
        }
    }
  else
    {
      // Load error rollback is an annoyance if not in server/interactive mode
      for (unsigned i = 0; i < system.hartCount(); ++i)
        {
          auto& hart = *system.ithHart(i);
          hart.enableLoadErrorRollback(false);
          hart.enableBenchLoadExceptions(false);
        }
    }

  if (serverMode)
    return runServer(system, args.serverFile, traceFile, cmdLog);

  if (args.interactive)
    {
      // Ignore keyboard interrupt for most commands. Long running
      // commands will enable keyboard interrupts while they run.
#ifdef __MINGW64__
      signal(SIGINT, kbdInterruptHandler);
#else
      struct sigaction newAction;
      sigemptyset(&newAction.sa_mask);
      newAction.sa_flags = 0;
      newAction.sa_handler = kbdInterruptHandler;
      sigaction(SIGINT, &newAction, nullptr);
#endif

      Interactive interactive(system);
      return interactive.interact(traceFile, cmdLog);
    }

  if (args.snapshotPeriod and *args.snapshotPeriod)
    {
      uint64_t period = *args.snapshotPeriod;
      std::string dir = args.snapshotDir;
      if (system.hartCount() == 1)
        return snapshotRun(system, traceFile, dir, period);
      std::cerr << "Warning: Snapshots not supported for multi-thread runs\n";
    }

  bool waitAll = not args.quitOnAnyHart;
  return batchRun(system, traceFile, waitAll);
}


/// Santize memory parameters. Page/region sizes must be greater and
/// equal to 4 and must be powers of 2.
/// Region size must be a multiple of page size.
/// Memory size must be a multiple of region size.
/// Return true if given parameters are good. False if any parameters
/// is changed to meet expectation.
static
bool
checkAndRepairMemoryParams(size_t& memSize, size_t& pageSize,
                           size_t& regionSize)
{
  bool ok = true;

  unsigned logPageSize = static_cast<unsigned>(std::log2(pageSize));
  size_t p2PageSize = size_t(1) << logPageSize;
  if (p2PageSize != pageSize)
    {
      std::cerr << "Memory page size (0x" << std::hex << pageSize << ") "
		<< "is not a power of 2 -- using 0x" << p2PageSize << '\n'
		<< std::dec;
      pageSize = p2PageSize;
      ok = false;
    }

  if (pageSize < 64)
    {
      std::cerr << "Page size (" << pageSize << ") is less than 64. Using 64.\n";
      pageSize = 64;
      ok = false;
    }

  size_t logRegionSize = static_cast<size_t>(std::log2(regionSize));
  size_t p2RegionSize = size_t(1) << logRegionSize;
  if (p2RegionSize != regionSize)
    {
      std::cerr << "Memory region size (0x" << std::hex << regionSize << ") "
		<< "is not a power of 2 -- using 0x" << p2RegionSize << '\n'
		<< std::dec;
      regionSize = p2RegionSize;
      ok = false;
    }

  if (regionSize < pageSize)
    {
      std::cerr << "Memory region size (0x" << std::hex << regionSize << ") "
		<< "smaller than page size (0x" << pageSize << ") -- "
		<< "using page size\n" << std::dec;
      regionSize = pageSize;
      ok = false;
    }

  size_t pagesInRegion = regionSize / pageSize;
  size_t multiple = pagesInRegion * pageSize;
  if (multiple != regionSize)
    {
      std::cerr << "Memory region size (0x" << std::hex << regionSize << ") "
		<< "is not a multiple of page size (0x" << pageSize << ") -- "
		<< "using 0x" << multiple << " as region size\n" << std::dec;
      regionSize = multiple;
      ok = false;
    }

  if (memSize < pageSize)
    {
      std::cerr << "Memory size (0x" << std::hex << memSize << ") "
		<< "smaller than page size (0x" << pageSize << ") -- "
                << "using 0x" << pageSize << " as memory size\n" << std::dec;
      memSize = pageSize;
      ok = false;
    }

  size_t pageCount = memSize / pageSize;
  if (pageCount * pageSize != memSize)
    {
      pageCount++;
      size_t newSize = pageCount * pageSize;
      std::cerr << "Memory size (0x" << std::hex << memSize << ") is not a "
		<< "multiple of page size (0x" << pageSize << ") -- "
		<< "using 0x" << newSize << '\n' << std::dec;
      memSize = newSize;
      ok = false;
    }

  return ok;
}


static
bool
getPrimaryConfigParameters(const Args& args, const HartConfig& config,
                           unsigned& hartsPerCore, unsigned& coreCount,
                           size_t& pageSize, size_t& memorySize)
{
  config.getHartsPerCore(hartsPerCore);
  if (args.hasHarts)
    hartsPerCore = args.harts;
  if (hartsPerCore == 0 or hartsPerCore > 16)
    {
      std::cerr << "Unsupported hart count: " << hartsPerCore;
      std::cerr << " (1 to 16 currently suppored)\n";
      return false;
    }

  config.getCoreCount(coreCount);
  if (args.hasCores)
    coreCount = args.cores;
  if (coreCount == 0 or coreCount > 16)
    {
      std::cerr << "Unsupported core count: " << coreCount;
      std::cerr << " (1 to 16 currently suppored)\n";
      return false;
    }

  // Determine simulated memory size. Default to 4 gigs.
  // If running a 32-bit machine (pointer size = 32 bits), try 2 gigs.
  if (memorySize == 0)
    memorySize = size_t(1) << 31;  // 2 gigs
  config.getMemorySize(memorySize);
  if (args.memorySize)
    memorySize = *args.memorySize;

  if (not config.getPageSize(pageSize))
    pageSize = args.pageSize;

  return true;
}


template <typename URV>
static
bool
session(const Args& args, const HartConfig& config)
{
  // Collect primary configuration paramters.
  unsigned hartsPerCore = 1;
  unsigned coreCount = 1;
  size_t pageSize = 4*1024;
  size_t regionSize = 256*1024*1024;
  size_t memorySize = size_t(1) << 32;  // 4 gigs

  if (not getPrimaryConfigParameters(args, config, hartsPerCore, coreCount,
                                     pageSize, memorySize))
    return false;

  checkAndRepairMemoryParams(memorySize, pageSize, regionSize);

  // Create cores & harts.
  unsigned hartIdOffset = hartsPerCore;
  config.getHartIdOffset(hartIdOffset);
  if (hartIdOffset < hartsPerCore)
    {
      std::cerr << "Invalid core_hart_id_offset: " << hartIdOffset
                << ",  must be greater than harts_per_core: " << hartsPerCore << '\n';
      return false;
    }
  System<URV> system(coreCount, hartsPerCore, hartIdOffset, memorySize, pageSize);
  assert(system.hartCount() == coreCount*hartsPerCore);
  assert(system.hartCount() > 0);

  // Configure harts. Define callbacks for non-standard CSRs.
  bool userMode = args.isa.find_first_of("uU") != std::string::npos;
  if (not config.configHarts(system, userMode, args.verbose))
    if (not args.interactive)
      return false;

  // Configure memory.
  if (not config.configMemory(system, args.iccmRw, args.unmappedElfOk, args.verbose))
    return false;

  if (args.hexFiles.empty() and args.expandedTargets.empty()
      and not args.interactive)
    {
      std::cerr << "No program file specified.\n";
      return false;
    }

  FILE* traceFile = nullptr;
  FILE* commandLog = nullptr;
  FILE* consoleOut = stdout;
  FILE* bblockFile = nullptr;
  if (not openUserFiles(args, traceFile, commandLog, consoleOut, bblockFile))
    return false;

  for (unsigned i = 0; i < system.hartCount(); ++i)
    {
      auto& hart = *system.ithHart(i);
      hart.setConsoleOutput(consoleOut);
      if (bblockFile)
	hart.enableBasicBlocks(bblockFile, args.bblockInsts);
      hart.reset();
    }

  bool result = sessionRun(system, args, traceFile, commandLog);

  auto& hart0 = *system.ithHart(0);
  if (not args.instFreqFile.empty())
    result = reportInstructionFrequency(hart0, args.instFreqFile) and result;

  closeUserFiles(args, traceFile, commandLog, consoleOut, bblockFile);

  return result;
}


/// Determine regiser width (xlen) from ELF file.  Return true if
/// successful and false otherwise (xlen is left unmodified).
static
bool
getXlenFromElfFile(const Args& args, unsigned& xlen)
{
  if (args.expandedTargets.empty())
    return false;

  // Get the length from the first target.
  auto& elfPath = args.expandedTargets.front().front();
  bool is32 = false, is64 = false, isRiscv = false;
  if (not Memory::checkElfFile(elfPath, is32, is64, isRiscv))
    return false;  // ELF does not exist.

  if (not is32 and not is64)
    return false;

  if (is32 and is64)
    {
      std::cerr << "Error: ELF file '" << elfPath << "' has both"
		<< " 32  and 64-bit calss\n";
      return false;
    }

  if (is32)
    xlen = 32;
  else
    xlen = 64;

  if (args.verbose)
    std::cerr << "Setting xlen to " << xlen << " based on ELF file "
	      <<  elfPath << '\n';
  return true;
}


/// Obtain integer-register width (xlen). Command line has top
/// priority, then config file, then ELF file.
static
unsigned
determineRegisterWidth(const Args& args, const HartConfig& config)
{
  unsigned isaLen = 0;
  if (not args.isa.empty())
    {
      if (boost::starts_with(args.isa, "rv32"))
	isaLen = 32;
      else if (boost::starts_with(args.isa, "rv64"))
	isaLen = 64;
    }

  // 1. If --xlen used, go with that.
  if (args.hasRegWidth)
    {
      unsigned xlen = args.regWidth;
      if (args.verbose)
	std::cerr << "Setting xlen from --xlen: " << xlen << "\n";
      if (isaLen and xlen != isaLen)
	{
	  std::cerr << "Xlen value from --xlen (" << xlen
		    << ") different from --isa (" << args.isa
		    << "), using: " << xlen << "\n";
	}
      return xlen;
    }

  // 2. If --isa specifies xlen, go with that.
  if (isaLen)
    {
      if (args.verbose)
        std::cerr << "Setting xlen from --isa: " << isaLen << "\n";
      return isaLen;
    }

  // 3. If config file specifies xlen, go with that.
  unsigned xlen = 32;
  if (config.getXlen(xlen))
    {
      if (args.verbose)
	std::cerr << "Setting xlen from config file: " << xlen << "\n";
      return xlen;
    }

  // 4. Get xlen from ELF file.
  if (getXlenFromElfFile(args, xlen))
    {
      if (args.verbose)
	std::cerr << "Setting xlen from ELF file: " << xlen << "\n";
    }
  else if (args.verbose)
    std::cerr << "Using default for xlen: " << xlen << "\n";
  
  return xlen;
}


int
main(int argc, char* argv[])
{
  Args args;
  if (not parseCmdLineArgs(argc, argv, args))
    return 1;

  if (args.help)
    return 0;

  // Expand each target program string into program name and args.
  args.expandTargets();

  // Load configuration file.
  HartConfig config;
  if (not args.configFile.empty())
    if (not config.loadConfigFile(args.configFile))
      return 1;

  unsigned regWidth = determineRegisterWidth(args, config);

  bool ok = true;

  try
    {
      if (regWidth == 32)
	ok = session<uint32_t>(args, config);
      else if (regWidth == 64)
	ok = session<uint64_t>(args, config);
      else
	{
	  std::cerr << "Invalid register width: " << regWidth;
	  std::cerr << " -- expecting 32 or 64\n";
	  ok = false;
	}
    }
  catch (std::exception& e)
    {
      std::cerr << e.what() << '\n';
      ok = false;
    }
	
  return ok? 0 : 1;
}
