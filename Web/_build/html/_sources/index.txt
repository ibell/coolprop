.. CoolProp documentation file


*********************
Welcome to CoolProp
*********************

.. toctree::
   :maxdepth: 1

   Installation.rst

What it is:
===========

An open-source database of refrigerant properties, for use in a wide range of programming environments.  The current refrigerants available at this time are

* Argon
* Nitrogen
* R134a
* R290 (Propane)
* R32 (No transport properties)
* R404A (Pseudo-pure properties)
* R407C (Pseudo-pure properties)
* R410A (Pseudo-pure properties)
* R507A (Pseudo-pure properties)
* R717 (Ammonia)
* R744 (Carbon Dioxide)

Originally developed for use in the programming language c, wrappers have been written for Python, and a DLL that can be called from any other language.  In addition, so called flooded properties are included for a mixture of liquids and refrigerants.

Who made it?
============

Ian Bell (ian.h.bell .at. gmail.com), during his time as a Ph.D. student in the Herrick Laboratories of Purdue University.
    
How to get it?
==============

To grab a zip file with the CoolProp DLL, head to https://sourceforge.net/projects/coolprop/files/ and download the most recent version.  Or to pull CoolProp source code or Python installer files.

CoolProp is an open-source project, and is actively looking for developers.  The project is hosted in a subversion repository at http://coolprop.svn.sourceforge.net/viewvc/coolprop/trunk/ .  In order to download the sources, click on the "Download GNU tarball" at the bottom of the page.  To unpack the tarball, `7-zip <http://www.7-zip.org>`_  is recommended.

How do I use it (5 seconds introduction)?
=========================================


MATLAB and languages other than Python:
---------------------------------------

The source folder includes a folder unimaginatively called "Examples".  In there is an example using MATLAB, copied verbatim: 

.. literalinclude:: ../Examples/MATLABExample.m

which yields::

	Functions in library CoolPropDLL:

	Help_dll
	[double, cstring] Props_dll(int8, int8, double, int8, double, cstring)
	[double, cstring] T_hp_dll(cstring, double, double, double)
	[double, cstring] Tcrit_dll(cstring)
	[double, cstring] Tsat_dll(cstring, double, double, double)
	[int32, cstring] errCode_dll(cstring)
	[double, cstring] h_sp_dll(cstring, double, double, double)
	[double, cstring] pcrit_dll(cstring)


	The critical temperature of R410A is:
	  344.4940

	The saturated vapor enthalpy of Propane at 275 K is:
	  576.7495

	The density of nitrogen at STP is:
	    1.1458

This example assumes that the compiled CoolpropDLL.dll and CoolProp_dll.h are located either on the MATLAB path somewhere, or in this case in particular, that the DLL is located one folder up and in the CoolPropDLL folder (that's what the relative path '..' means).  The protocol should remain basically the same for other programming languages other than Python.

This example demonstrates the two main types of calls to Coolprop

* Saturation temperature known in the two phase region
* Temperature and pressure known away from two-phase dome

See the rest of the documentation for further information.

Python:
-------

Python is very well suited to the use of code programmed in C.  With the use of SWIG, VERY minimal modifications are required to call c-code from Python.  With CoolProp installed properly, the same example as for MATLAB shows the Python syntax:

.. literalinclude:: ../Examples/PythonExample.py




