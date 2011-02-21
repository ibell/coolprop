#
#  Simple Makefile
#  For compiling C extension modules for 
#  python.
#

###################################################
#
# Basic Swig Procedure :~)
# 
# // First make the wrapper from the interface file
# swig -python utility.i
# // Then compile both the utility_wrap.c and the utility.c code.
# gcc -c utility.c utility_wrap.c -I/usr/include/python2.2
# // Then link :~)
# ld -shared utility.o utility_wrap.o -o _utility.so
#
# NOTE:
# With a little tweaking this makefile should work on OS X as well.
#
# You currently need to make sure that it matches your include paths
# and python version at the top of the file.
##################################################

# ------------------------- You'll need to set this stuff ---------------------------------------------------
# Set this as well:
python = python2.6

### (old) -I/usr/include/${python} -I/usr/local/include/ -I/usr/lib/${python}
EPD= /Library/Frameworks/Python.framework/Versions/6.1
### -I${EPD}/include/ -I${EPD}/lib/ -I/usr/local/include
linuxincludes =   -I/usr/include/${python}
finkincludes = -I/sw/include/${python} -I/sw/lib/${python}
# Change this line to suit your platform.  I should really get autoconfigure setup :~)
includes = ${linuxincludes}
#includes = ${finkincludes}

#---------------------- This you shouldn't have to mess with -----------------------------------------

interface = *.i
source = *.c
objects = *.o
library = _BellProp.so
junk = *.so *.o *_wrap.c *_wrap.doc
linkopts = -shared -o
ccopts = -c

# Useful Utility commands, just like Uname 
swiglibrary = $(shell swig -python -ldflags)
pythonlib = $(shell python-config)

UNAME=$(shell uname -s)
ifeq ($(UNAME),Linux)
	linkopts = ${swiglibrary} -shared -o
	ccopts = ${swiglibrary} -Wcast-qual -c
endif
ifeq ($(UNAME),Darwin)
	#### -dylib 
	linkopts = ${swiglibrary} -bundle  -flat_namespace -undefined suppress -arch i686 -Wall -o
	ccopts = ${swiglibrary} -arch i686 -Wall  -c
endif

${library} : ${objects} #pointer
	g++ ${linkopts} ${library} ${objects}

${objects} : interface ${source} 
	g++ ${ccopts} ${source} ${includes}

interface : ${interface}
	swig -python ${interface}

clean :
	-rm ${junk}

