The *.mexw64 files in this folder are wrappers of CoolProp for MATLAB on 64-bit Windows.

Place them somewhere on the MATLAB path

For Users
=========
There is an example file called MATLAB_sample.m that when run will demonstrate calling both
Props (for fluid properties) and HAProps(for humid air properties)

For Developers
==============
To build the mex file on other platforms, you should be able to just type

    matlab -r MATLABBuilder

at the command prompt.  Please send an email to ian.h.bell@gmail.com if it works or if you have problems

(c) Ian Bell, 2012