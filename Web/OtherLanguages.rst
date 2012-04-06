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

This example assumes that the compiled CoolpropDLL.dll and CoolProp_dll.h are located either on the MATLAB path somewhere, or in this case in particular, that the DLL is located one folder up and in the CoolPropDLL folder (that's what the relative path '..' means).  The protocol should remain basically the same for other programming languages other than MATLAB.

This example demonstrates the two main types of calls to Coolprop

* Saturation temperature known in the two phase region
* Temperature and pressure known away from two-phase dome

See the rest of the documentation for further information.