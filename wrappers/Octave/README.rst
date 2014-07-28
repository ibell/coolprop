A wrapper of CoolProp for the Octave programming language (an open-source version of MATLAB)

This CoolProp.oct module is built using SWIG, and is built for octave for windows built with MinGW.  

Installation
============
A version of the .oct file is available for Octave 3.6.4 for mingw on windows. Version 3.8.0 is currently aviable in Ubuntu 14.04 repositories. 

Put the .oct file for the version you have in somewhere in the octave path, or use the ``addpath`` function to add the folder that contains CoolProp.oct to the Octave path

On Linux systems you can put .oct file in
"/usr/share/octave/octave.version.number/m" folder. You will need superuser
privileges to do this.

If you place .oct file somewhere outside octave path, you have to use
"addpath" function at begining of your code.

Example: adding the folder that contains CoolProp.oct file to the Octave path:
    "addpath('/home/?user_name?/Some_folder/CoolProp')"

Developer Notes:
===============


Building on Linux (Ubuntu and derivatives)
---------------------------------------
1. You will need to run 
      sudo apt-get update
      sudo apt-get install octave liboctave-dev swig
   to install the necesary dependencies.  The install of octave might not be necessary but it cant hurt
2. Check out the full source for coolprop from github
      git clone https://github.com/ibell/coolprop
3. Change into the coolprop-code/wrappers/Octave folder
      cd coolprop/wrappers/Octave
4. Call
      octave _OctaveBuilder_Linux.m
5. Call
      octave sample_code.m
   to run the sample


Building on Linux (openSUSE)
---------------------------------------
1. In YaST program manager select for installation the following libraries: 
      octave
      octave-devel
      octave-doc
      octave-mathgl
      plplot-octave
      swig
   Accept the install of the additional necessary dependencies. The install of octave might not be necessary but it cant hurt
2. Check out the full source for coolprop from github
      git clone https://github.com/ibell/coolprop
3. Change into the coolprop-code/wrappers/Octave folder
      cd coolprop/wrappers/Octave
4. Call
      octave _OctaveBuilder_Linux.m
5. Call
      octave sample_code.m
   to run the sample

   
Building on Raspberry PI
------------------------
1. You will need to run
      sudo aptitude update
      sudo aptitude install octave liboctave-dev swig
2. Download all the sources from subversion using
      git clone https://github.com/ibell/coolprop
3. Change into the folder
      cd coolprop/wrappers/Octave
4. Run the build script
      octave _OctaveBuilder_Linux.m
5. Call 
      octave sample_code.m
    to run the sample
    

