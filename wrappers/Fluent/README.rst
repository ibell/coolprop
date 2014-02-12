CoolProp wrapper for FLUENT
===========================

Contributors
------------
Primary CoolProp Developer: Ian Bell, University of Liege, Belgium (ian.h.bell@gmail.com)
FLUENT experts: Joris Degroote and Iva Papes, University of Gent, Belgium
First release: October 17, 2013

Requirements
------------
A linux version of FLUENT
g++

To Build
--------

Let us call the directory the directory where the Fluent wrapper is (coolprop/wrappers/Fluent) as CUSTOM_DIRECTORY.

1. Open compile.sh shell file and edit variables SOLVER and FLUENT_BIN_FOLDER.
   a. The SOLVER variable should specify what kind of solver you want to run (2d, 2ddp, 3d, 3ddp).
   b. The FLUENT_BIN_FOLDER variable should specify the path to the bin folder in your ANSYS installation. For example:
      FLUENT_BIN_FOLDER="/home/ansys_inc/v145/fluent/bin"
      In case you can run Fluent from terminal command line, you can leave it as it is ("NULL").
   c. Do not modify anything else.
2. Make sure your case file is in the same directory as compile.sh and UDF.c.
   a. If they are not, move your case/mesh files to CUSTOM_DIRECTORY. Do NOT move compile.sh and UDF.c.
3. Run the script compile.sh (sh compile.sh), this should generate the libudf folder.
   a. Several warnings may show up, those should not be a problem.
4. Run Fluent.
   a. Make sure it runs the same solver as specified in the SOLVER variable
5. Open your case/mesh.
6. Compile the UDF.c.
   a. Using the Fluent interface, go to Define > User-Defined > Functions > Compiled
   b. Select the CUSTOM_DIRECTORY/UDF.c file
   c. Make sure "libudf" is written in the Library field
   d. Hit Load
7. The default UDF.c provided by the Coolprop wrapper is an EXECUTE_ON_DEMAND file.
   a. To check if it is working: Define > User-defined > Execute on Demand
   b. Select "call_coolprop::libudf" and hit execute
  
Warning
-------
Absolutely no guarantee of utility or accuracy can be made, although we have done our best to ensure useful and accurate results.  Caveat emptor!
