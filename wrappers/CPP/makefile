# ============================================================================
# Name        : Makefile
# Author      : Rod Persky (rodney.persky@connect.qut.edu.au) 
# Version     : 0.0.3
# Copyright   : Use and modify at your own risk.
# Description : Makefile for CoolProp.
# ============================================================================
# The installation procedure should be as follows:
# 1) make
# 2) sudo make install
 
# If you want a static library then add LIBTYPE=a to *ALL* of the
#  make and install commands. Also whilst you're adding commands may as well
#  include -j $(nproc) to increase build speed. So for a static library use:
#  make LIBTYPE=a -j $(nproc)
#  sudo make LIBTYPE=a install
#
# If you need a clean install use:
#  sudo make uninstall install
#
# You can simply uninstall by using:
#  sudo make uninstall
#
# When using the library in Eclipse in your source files add:
#  #include <CoolProp/CoolProp.h>
#
# And in the linker Librarys add (dl is only needed for the static library):
#  CoolProp
#  dl
#
# Finally, in the Library search path add:
#  /usr/include/CoolProp
# ============================================================================
# General commands:
# These commands are faily generic accross all platforms
 
CC    := g++
RM    := rm -f
CP    := cp 
CH    := chmod 0644 
MK    := mkdir -p -m 655
LD    := ldconfig 
LN    := ln -sf 
LT    := libtool 
AR    := ar rvs
MV    := mv
SHELL := /bin/sh
LIB   := /usr/lib
 
# What library type are we going to build by default (shared)
LIBTYPE = so
 
# We're going to install into /usr/include because then when
#  using the library you don't get build warnings from the library
DESTDIR     = /usr/include/$(BASENAME)
COOLPROPDIR = ../../CoolProp/
BINDIR      = ./build
 
# Optimisation flags, change to your wants and desires
USER_OPTFLAGS = -O0 #-funroll-loops -ffast-math -ffloat-store
 
# Debug flags
DBGFLAGS = -O0 -g3
WRNFLAGS = -Wall -Wpedantic -Wextra
 
 
# =============================
# ~ END OF CUSTOMISABLE CODE ~
# =============================
# Sort out exactly what we want to do,
# shared or static, debug or standard
# build. Get the final names of the execs
# =============================
 
BASENAME=CoolProp
VFILE = $(COOLPROPDIR)/version.h
VERSION=$(shell cat $(VFILE) | grep -e [0-9\.]* -o)
MAJOR=$(shell cat $(VFILE) | grep -e \"[0-9] -o | grep [0-9] -o)
MINOR=$(shell cat $(VFILE) | grep -e \\.[0-9]\\. -o | grep [0-9] -o)
 
ifeq ($(LIBTYPE),a)
  LIBTYPECLEAN=a
else
  LIBTYPECLEAN=so
endif
 
ifeq ($(DEBUG),debug)
  DEBUGCLEAN=_dbg
endif
 

EXECUTABLE=lib$(BASENAME)$(DEBUGCLEAN).$(LIBTYPECLEAN).$(VERSION)
SONAME=lib$(BASENAME)$(DEBUGCLEAN).$(LIBTYPECLEAN).$(MAJOR).$(MINOR)
EXECUTABLE_MAJOR=lib$(BASENAME)$(DEBUGCLEAN).$(LIBTYPECLEAN).$(MAJOR)
EXECUTABLE_BASE=lib$(BASENAME)$(DEBUGCLEAN).$(LIBTYPECLEAN)
 
 
# ======================
# Compile all the flags
# ======================
 
# The standard flags are meant to represent the safest compilation,
#  not fastest resulting library speed.
 
_LIBRARIES = -ldl
_INCLUDES = -I$(COOLPROPDIR)
 
 
STD_OPTFLAGS = -fPIC
_CPPFLAGS=$(_LIBRARIES) $(STD_OPTFLAGS) $(USER_OPTFLAGS)
 
ifeq ($(LIBTYPE),a)
  #AR doesn't take flags
else
  _LIBFLAGS=-shared -Wl,-soname,$(SONAME)
endif
 
ifeq ($(DEBUG),debug)
  _CPPFLAGS=$(STD_OPTFLAGS) $(DBGFLAGS) $(WRNFLAGS)
endif
 
#Get sources
COOLCPP_FILES=$(shell echo $(COOLPROPDIR)*.cpp)
COOLOBJ_FILES  := $(addprefix $(BINDIR)/,$(notdir $(COOLCPP_FILES:.cpp=.o)))
 
# Path for prerequisite files
vpath %.cpp $(COOLPROPDIR)

# ========================
# Compile coolprop sources
# ========================
 
.PHONY: all
all: $(COOLOBJ_FILES) $(BINDIR)/$(EXECUTABLE)

$(COOLOBJ_FILES): | $(BINDIR)

$(BINDIR):
	$(MK) -m 755 $(BINDIR)

# Make the .o files (for both static and shared)
$(BINDIR)/%.o : %.cpp
	$(CC) $(_CPPFLAGS) -c $(_INCLUDES) $< -o $@

# Link shared library
$(BINDIR)/lib$(BASENAME)$(DEBUGCLEAN).so.$(VERSION): $(COOLOBJ_FILES) 
	$(CC) $(_LIBFLAGS) $(COOLOBJ_FILES) -o $@

# Link static library
$(BINDIR)/lib$(BASENAME)$(DEBUGCLEAN).a.$(VERSION): $(COOLOBJ_FILES)
	$(AR) $@ $^

# ========================
# Install and Uninstall
# ========================

.PHONY : install
install:
	$(MK) $(DESTDIR)
	$(CP) -Ru $(COOLPROPDIR)/* $(DESTDIR)
	find $(DESTDIR) -type f ! -name "*.h" -delete
	$(CP) $(BINDIR)/$(EXECUTABLE) $(DESTDIR)
	$(LN) $(DESTDIR)/$(EXECUTABLE) $(DESTDIR)/$(SONAME)
	$(LN) $(DESTDIR)/$(EXECUTABLE) $(DESTDIR)/$(EXECUTABLE_MAJOR)
	$(LN) $(DESTDIR)/$(EXECUTABLE) $(DESTDIR)/$(EXECUTABLE_BASE)
	$(LN) $(DESTDIR)/$(EXECUTABLE) $(LIB)/$(SONAME)
	$(LD) -n $(DESTDIR)
	$(LD)

.PHONY : uninstall 
uninstall:
	rm -rf $(DESTDIR)
	rm -rf /usr/lib/lib$(BASENAME)$(DEBUGCLEAN).so.$(MAJOR).$(MINOR)
	rm -rf /usr/lib/lib$(BASENAME)$(DEBUGCLEAN).a.$(MAJOR).$(MINOR)
	$(LD)
	
# ========================
# Accessory functions
# ========================

.PHONY : clean
clean:
	$(RM) $(BINDIR)/*
 
.PHONY : check_flags
check_flags:
	@printf "\nBuild settings:\n\
	  EXECUTABLE=$(EXECUTABLE)\n\
	  SONAME=$(SONAME)\n\
	  LIBFLAGS=$(_LIBFLAGS)\n\
	  CPPFLAGS=$(_CPPFLAGS)\n\
	  INCLUDES=$(_INCLUDES)\n\n\
	  COOLCPP_FILES:\n$(COOLCPP_FILES)\n\n\
	  COOLOBJ_FILES:\n$(COOLOBJ_FILES)\n\n"
