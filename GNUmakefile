INSTALL_DIR := .

PROJECT := whisper

# For Static Linking to Boost Library
# STATIC_LINK := 1
# For Dynamic linking to Boost Library change the following line to
# STATIC_LINK := 0
# or
# make STATIC_LINK=0
STATIC_LINK := 1

# For non-default compiler toolchain uncomment and change the following variables
#CC := gcc-8
#CXX := g++-8
#AR := gcc-ar-8
# Or run make with these options
# $ make CC=gcc-8 CXX=g++-8 AR=gcc-ar-8

# We use boost 1.67.
# Set the BOOST_ROOT environment variable to point to the base install
# location of the Boost Libraries
BOOST_ROOT := /wdc/apps/utilities/boost-1.67
BOOST_DIR := $(BOOST_ROOT)
# For Various Installation types of Boost Library
BOOST_INC := $(wildcard $(BOOST_DIR) $(BOOST_DIR)/include)

# These boost libraries must be compiled with: "g++ -std=c++14" or "g++ -std=c++17"
# For Various Installation types of Boost Library
BOOST_LIB_DIR := $(wildcard $(BOOST_DIR)/stage/lib $(BOOST_DIR)/lib)

# Specify only the basename of the Boost libraries
BOOST_LIBS := boost_program_options

# Add extra dependency libraries here
ifeq (CYGWIN_NT-10.0,$(shell uname -s))
EXTRA_LIBS := -lpthread -lz -lstdc++fs
else
EXTRA_LIBS := -lpthread -lz -static-libstdc++
endif

ifeq (Linux,$(shell uname -s))
EXTRA_LIBS += -lstdc++fs
endif

ifeq (mingw,$(findstring mingw,$(shell $(CXX) -v 2>&1 | grep Target | cut -d' ' -f2)))
EXTRA_LIBS += -lws2_32
endif

ifdef SOFT_FLOAT
override CPPFLAGS += -I$(PWD)/softfloat/source/include
override CPPFLAGS += -DSOFT_FLOAT
soft_float_build := $(wildcard $(PWD)/softfloat/build/RISCV-GCC)
soft_float_lib := $(soft_float_build)/softfloat.a
endif

ifdef MEM_CALLBACKS
override CPPFLAGS += -DMEM_CALLBACKS
endif


# Add External Library location paths here
LINK_DIRS := $(addprefix -L,$(BOOST_LIB_DIR))

# Generating the Linker options for dependent libraries
ifeq ($(STATIC_LINK), 1)
  LINK_LIBS := $(addprefix -l:lib, $(addsuffix .a, $(BOOST_LIBS))) $(EXTRA_LIBS)
else
  COMMA := ,
  LINK_DIRS += $(addprefix -Wl$(COMMA)-rpath=, $(BOOST_LIB_DIR))
  LINK_LIBS := $(addprefix -l, $(BOOST_LIBS)) $(EXTRA_LIBS)
endif

ifeq (Darwin,$(shell uname))
   LINK_LIBS := $(BOOST_LIB_DIR)/lib$(BOOST_LIBS).a $(EXTRA_LIBS)
endif

# For out of source build
BUILD_DIR := build-$(shell uname -s)
MKDIR_P ?= mkdir -p
RM := rm -rf
# Optimization flags.  Use -g for debug.
OFLAGS := -O3

# Include paths.
IFLAGS := $(addprefix -isystem ,$(BOOST_INC)) -I.

# Command to compile .cpp files.
ifeq (CYGWIN_NT-10.0,$(shell uname -s))
override CXXFLAGS += -MMD -MP -mfma -std=c++17 -D_GNU_SOURCE $(OFLAGS) $(IFLAGS) -pedantic -Wall -Wextra
else
override CXXFLAGS += -MMD -MP -mfma -std=c++17 $(OFLAGS) $(IFLAGS) -fPIC -pedantic -Wall -Wextra
endif

# Rule to make a .o from a .cpp file.
$(BUILD_DIR)/%.cpp.o:  %.cpp
	@if [ ! -d "$(dir $@)" ]; then $(MKDIR_P) $(dir $@); fi
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

# Rule to make a .o from a .c file.
$(BUILD_DIR)/%.c.o:  %.c
	@if [ ! -d "$(dir $@)" ]; then $(MKDIR_P) $(dir $@); fi
	$(CC) $(CFLAGS) -c -o $@ $<

# Main target.(only linking)
$(BUILD_DIR)/$(PROJECT): $(BUILD_DIR)/whisper.cpp.o \
                         $(BUILD_DIR)/librvcore.a \
			 $(soft_float_lib)
	$(CXX) -o $@ $^ $(LINK_DIRS) $(LINK_LIBS)

# List of all CPP sources needed for librvcore.a
RVCORE_SRCS := IntRegs.cpp CsRegs.cpp FpRegs.cpp instforms.cpp \
            Memory.cpp Hart.cpp InstEntry.cpp Triggers.cpp \
            PerfRegs.cpp gdb.cpp HartConfig.cpp \
            Server.cpp Interactive.cpp decode.cpp disas.cpp \
	    Syscall.cpp PmaManager.cpp DecodedInst.cpp snapshot.cpp \
	    PmpManager.cpp VirtMem.cpp Core.cpp System.cpp Cache.cpp \
	    Tlb.cpp VecRegs.cpp vector.cpp wideint.cpp float.cpp bitmanip.cpp \
	    amo.cpp SparseMem.cpp

# List of All CPP Sources for the project
SRCS_CXX += $(RVCORE_SRCS) whisper.cpp

# List of All C Sources for the project
SRCS_C :=

# List of all object files for the project
OBJS_GEN := $(SRCS_CXX:%=$(BUILD_DIR)/%.o) $(SRCS_C:%=$(BUILD_DIR)/%.o)

# List of all auto-genreated dependency files.
DEPS_FILES := $(OBJS_GEN:.o=.d)

# Include Generated Dependency files if available.
-include $(DEPS_FILES)

# Object files needed for librvcore.a
OBJS := $(RVCORE_SRCS:%=$(BUILD_DIR)/%.o) $(SRCS_C:%=$(BUILD_DIR)/%.o)

$(BUILD_DIR)/librvcore.a: $(OBJS)
	$(AR) cr $@ $^

$(soft_float_lib):
	$(MAKE) -C $(soft_float_build)

install: $(BUILD_DIR)/$(PROJECT)
	@if test "." -ef "$(INSTALL_DIR)" -o "" == "$(INSTALL_DIR)" ; \
         then echo "INSTALL_DIR is not set or is same as current dir" ; \
         else echo cp $^ $(INSTALL_DIR); cp $^ $(INSTALL_DIR); \
         fi

clean:
	$(RM) $(BUILD_DIR)/$(PROJECT) $(OBJS_GEN) $(BUILD_DIR)/librvcore.a $(DEPS_FILES)

help:
	@echo "Possible targets: $(BUILD_DIR)/$(PROJECT) install clean"
	@echo "To compile for debug: make OFLAGS=-g"
	@echo "To install: make INSTALL_DIR=<target> install"
	@echo "To browse source code: make cscope"

cscope:
	( find . \( -name \*.cpp -or -name \*.hpp -or -name \*.c -or -name \*.h \) -print | xargs cscope -b ) && cscope -d && $(RM) cscope.out

.PHONY: install clean help cscope

