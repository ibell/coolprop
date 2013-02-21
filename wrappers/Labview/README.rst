Labview wrapper of CoolProp
============================

by Ian Bell and Arnaud Legros

University of Liege and Bell Thermal Consultants

February 2013

To Install
----------
1. Copy the files CoolProp.vi and CoolProp.dll from this folder to somewhere you want
2. Add CoolProp module to your code

To Use
------
Call it using the same sorts of input parameters as Props function : http://coolprop.sourceforge.net/apidoc/CoolProp.html#CoolProp.CoolProp.Props

Notes
-----
Wrapper currently in its infancy, absolutely no error checking is carried out, use at your own risk.

For Developers
--------------

1. To regenerate DLL, run the build script (wrappers/Labview/BuildDLL.bat).  You will need to have Visual Studio 2010 installed, or change the path to vcvarsall.bat in build script
2. Uses __cdecl in combination with extern "C" to make Labview happy
