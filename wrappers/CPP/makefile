# ============================================================================
# Name        : Makefile
# Author      : Rod Persky (rodney.persky@connect.qut.edu.au) 
# Version     : 0
# Copyright   : Use and modify at your own risk.
# Description : Makefile for CoolProp.
# ============================================================================
# The installation procedure should be as follows:
# 1) make
# 2) sudo make install
 
# Decide if you want a static library then add LIBTYPE=a to the make command
# whilst you're there may as well add -j $(nproc) to increase build speed:
#   make LIBTYPE=a -j $(nproc)
#
# If you need a clean install replace sudo make install with:
#   sudo make uninstall install
#
# You can simply uninstall by using:
# 3) sudo make uninstall
#
# When using the library (in Eclipse say) in your .c files add
#    #include <CoolProp/CoolProp.h>
# and in the linker settings add CoolProp to the libraries and
#    /usr/include/CoolProp to the library search path.
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
 
# Optimisation flags, taylor to your wants and desires
USER_OPTFLAGS = -funroll-loops -ffast-math #-ffloat-store
 
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
 
# define exec name here because we have ensured it's either .a or .so
# in the prior ifeq statement. The SONAME assumes there is no
# drastic changes between minor versions. It also assumes that
# estimations of properties will get better and that a application
# targeting an older version can still get the same data in the
# same way (ID's don't change).
 
# ToDo: Check the above assumption otherwise add minor to the SONAME
 
EXECUTABLE=lib$(BASENAME)$(DEBUGCLEAN).$(LIBTYPECLEAN).$(VERSION)
SONAME=lib$(BASENAME)$(DEBUGCLEAN).$(LIBTYPECLEAN).$(MAJOR).$(MINOR)
EXECUTABLE_MAJOR=lib$(BASENAME)$(DEBUGCLEAN).$(LIBTYPECLEAN).$(MAJOR)
EXECUTABLE_BASE=lib$(BASENAME)$(DEBUGCLEAN).$(LIBTYPECLEAN)
 
 
# ======================
# Compile all the flags
# ======================
 
#CPPFLAGS= $(OPTFLAGS) $(EXTFLAGS) -fPIC -c -Wall -I$(shell pwd)
#LDFLAGS=-shared -ldl -Wl,-soname,$(EXECUTABLE)
 
# The standard flags are meant to represent the safest compilation,
#  not fastest resulting library speed.
 
_LIBRARIES = -ldl
_INCLUDES = -I$(COOLPROPDIR)
 
 
STD_OPTFLAGS = -fPIC
_CPPFLAGS=$(STD_OPTFLAGS) $(USER_OPTFLAGS)
 
ifeq ($(LIBTYPE),a)
  _LIBFLAGS=-static $(_LIBRARIES)
else
  _LIBFLAGS=-shared -Wl,-soname,$(SONAME) $(_LIBRARIES)
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

$(BINDIR)/$(EXECUTABLE): $(COOLOBJ_FILES) 
	$(CC) $(_LIBFLAGS) $(COOLOBJ_FILES) -o $@
 
$(BINDIR)/%.o : %.cpp
	$(CC) $(_CPPFLAGS) -c $(_INCLUDES) $< -o $@

# ========================
# Install and Uninstall
# ========================

.PHONY : install
install:
	$(MK) $(DESTDIR)
	$(CP) -Rn $(COOLPROPDIR)/* $(DESTDIR)
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
	rm -rf /usr/lib/$(SONAME)
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
