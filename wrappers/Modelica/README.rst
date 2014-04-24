
Notice
======

Please note that CoolProp support is about to be added to ExternalMedia, which can 
be considered the standard package for pure and pseudo-pure fluid property calculations 
in the Modelica language. You can find the latest version at 
https://github.com/modelica/ExternalMedia. 

As of March 2014, CoolProp is fully integrated with the developement version of ExternalMedia, 
but there is no official release with CoolProp support, yet. Some more changes to the basic 
structure of ExternalMedia are going to be implemented soon.

While we wait for the next release of ExternalMedia, please feel free to use the wrapper 
provided here.


Wrapper of CoolProp for Modelica
================================

About
-----
This is the first fully-open-source fluid property package for Modelica

Installation
------------
Make the visual studio static library for your compiler by running one of the scripts in this folder.
Note that you have to download the CoolProp source code as well. Place the static library you want to 
use in the C:\\Program Files (x86)\\Dymola 2013\\bin\\lib folder (or the relevant path for your 
installation) and put the CoolPropLib.h file in the C:\\Program Files (x86)\\Dymola 2013\\Source folder 
(or the relevant path for your installation).

Please visit https://github.com/thermocycle/CoolProp2Modelica-library to download a copy of the Modelica 
files required to access the CoolProp database.

Usage
-----
See the examples folder and the general CoolProp documentation.

Credits
-------
The CoolProp developers and the developers of the initial ExternalMedia release (this wrapper is based on it).
