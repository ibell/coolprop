The *.mexw64 files in this folder are wrappers of CoolProp for MATLAB on 64-bit Windows.

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

You will need to have subversion installed (google it)

On OSX/Linux, the same idea.  Do this::

    svn checkout svn://svn.code.sf.net/p/coolprop/code/trunk coolprop-code
    cd coolprop-code/wrappers/MATLAB
    matlab -r MATLABBuilder_OSX.m

Please send an email to ian.h.bell@gmail.com if it works or if you have problems

(c) Ian Bell, 2012