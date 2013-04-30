The *.mexw64 files in this folder are wrappers of CoolProp for 64-bit MATLAB
The *.mexw32 files in this folder are wrappers of CoolProp for 32-bit MATLAB
You can tell what type of MATLAB you have when you start up MATLAB.  It will say on the splash screen

Place them somewhere on the MATLAB path

For Users
=========
There is an example file called MATLAB_sample.m that when run will demonstrate calling both
Props (for fluid properties) and HAProps(for humid air properties)

For Developers
==============
To build the mex files on windows, you should enter these commands at the command prompt::

    svn checkout svn://svn.code.sf.net/p/coolprop/code/trunk coolprop-code
    cd coolprop-code/wrappers/MATLAB
    matlab -r MATLABBuilder.m

You will need to have subversion installed (google it).  You will also need a compiler installed (Visual studio express works).

On OSX/Linux, the same idea.  Do this::

    svn checkout svn://svn.code.sf.net/p/coolprop/code/trunk coolprop-code
    cd coolprop-code/wrappers/MATLAB
    # Change line 6 MATLABBuilder_OSX.m according to instructions in file.
    matlab -r MATLABBuilder_OSX.m

Please send an email to ian.h.bell@gmail.com if it works or if you have problems

(c) Ian Bell, 2012-