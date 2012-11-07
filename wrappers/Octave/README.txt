This CoolProp.oct module is built using SWIG, and is built for octave for windows built with Visual Studio.  Separate versions are available for 3.6.1 and 3.6.2

Developer Notes:
===============

A. 3.6.1 needs to use VS 2008 to build; 3.6.2 needs to use VS 2010

B.
In the win32 distro of 3.6.2, the hard-coded path in OCTAVE/bin/include/math.h around line 73 might need to be changed to 

/* Include VC++ original math.h */

#include <c:/Program Files (x86)/Microsoft Visual Studio 10.0/VC/include/math.h>

depending on where your VS is installed