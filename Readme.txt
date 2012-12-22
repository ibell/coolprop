Changelog:

2.4 (revision 240)
Added the fluids R1234ze, DME, R143a, n-Pentane, n-Hexane, n-Octane, n-Heptane, CycleHexane, 1-Butene,trans-2-Butene, cis-2-Butene,IsoButene, MethylLinoleate, MethylLinolenate, MethylOleate, MethylPalmitate, MethylStearate
Added C# wrappers (built for Windows) (r240)
Added Phase_Trho() and Phase_Tp() functions (r240)
Cleanup of the build process.  svnrevision is saved to a file that is built in.  Can access the svn revision through the functions get_svnrevision() and get_version()
Added a genetic algorithm to build ancillaries to dev folder (r226)
Added third partial derivatives of all the Helmholtz Energy terms (r238)
Bugfixes:
    Fixed Q(T,rho) (r237) (https://sourceforge.net/p/coolprop/tickets/42/)
    dhdT and dhdrho added back (r232)
    Surface tension now properly has the units of N/m as specified in the docs (r228)
    Fixed bug from Reiner with V and Vda (r227)
    Added a Brent solver to fix the solution for the saturation around the critical point (r220)(https://sourceforge.net/p/coolprop/tickets/38/)
    Repaired saturation LUT (r214-r216)
    Fixed bugs in IsFluidType as well as fixed bugs in Brine entropy calculations (r213)

2.3 (revision 212)
Added updated correlations for brines and subcooled liquids from Melinder 2010 (r207)
Added aliases to docs and python and DLL (r211)
Excel wrapper updated to catch errors and output them to a message box
Big speed update to p,Q as inputs (as fast as REFPROP now) (r202)
Doxygen now gets updated as well (r200)
Bugfixes:
    Updated inputs for brines (order doesn't matter) (r208)
    Fixed REFPROP with single-input props (r206)
    Fixed Manifest file for source distro (r206)
    Fixed bug with REFPROP mixtures not being properly parsed (r205 & r212)
    Added a backup Brent method for HAProps when solving at low humidity ratio: closed https://sourceforge.net/p/coolprop/tickets/32/ (r204)
    Added an example to show how to get version of CoolProp: closed https://sourceforge.net/p/coolprop/tickets/34/ (r204)
    Closed the bugs/issues in https://sourceforge.net/p/coolprop/tickets/35/ (r203)
    Resolved memory leak with ECS (r201)

2.2.5 (revision 199)
P,h and p,s as inputs solve for almost all fluids under almost all conditions
Octave modules for 3.6.1 and 3.6.2 now build and run properly for VS build on Windows
Builds properly on Linux now
Bugfixes:
    REFPROP.cpp bug with mixtures (r195)
    fixes around critical point (r198)
    Ancillaries for R134a updated in the vicinity of critical point

2.2.4 (revision 192)
Does not die if pseudo-pure T,P are in the two-phase region
Fixed bug with dewpoint as an input for dewpoints below 0C
Added a CoolPropStateClass for elegantly handling inputs - internal codebase will soon transition to this entirely
Fixed derivatives of drhodp|h and drhodh|p in two-phase region
Improved ancillary equations for Siloxanes (were terrible!)
Improved ancillary equations for Ethanol
Improved ancillary equations for SES36
Tmin is now an option for CoolProp and REFPROP fluids - ex: Props("REFPROP-MDM","Tmin") or Props("MDM","Tmin")
T_hp is now faster than REFPROP 
Added Excel 2003 Add-in for CoolProp - not clear it is working though
Improved the Distro builder


2.2.3 (revision 172)
Added Ethylene, SF6, Ethanol, Methane, Ethane, n-Butane, Isobutane
x(h,p) is much faster due to the avoidance of a lot of saturation routine calls
x(p,Q) is about 200 times faster!!
Added Quality 'Q' as an output
Fixed properties for Air
Fixed ancillaries for Siloxanes

2.2.2 (revision 169)
Added MATLAB wrappers and compiled versions on Windows to batch
Added plots to check solvers for (T,p) and (h,p) in subcooled liquid and superheated vapor regions

2.2.1 (revision 166)
Added the fluid SES36
HAProps added to CoolProp wrapper and added to Excel addin
When using pseudo-pure fluid, saturation density are calculated based on solving for density given T,P and guess value given by ancillary for density 
Improved saturated vapor ancillary for SES36
Changed default names: R717 -> Ammonia, R744 -> CarbonDioxide, R290 -> Propane

2.2.0 (revision 164)
Added the Siloxanes (MM,MDM,MD2M,MD3M,MD4M,D4,D5,D6)
Added a script that will build all the parts (Excel DLL, Python, MATLAB, etc.) and upload to Sourceforge
Very-alpha code for use of CoolProp in Modelica
Enthalpy and pressure are valid inputs for Brine fluids
Added support for quantities package in Python code (If you provide quantities.Quantity instance to CoolProp.CoolProp.Props, the units will be converted to the default units for CoolProp; Default units can be obtained by calling get_index_units(iParam) as a std::string; If a string for the desired output units is passed to Props the units will be converted to the output units)
Internals of CoolProp changed again, added a function called IProps that uses the integer indices for the input terms as well as the fluids - significant speedup.  This is mostly for use with CoolProp.State.State in Python although the same principle can be used elsewhere
Bug fixes for ECS

2.1.0 (revision 154)
Added the fluids Hydrogen, Oxygen, and Helium
Added the output term 'accentric' to get the accentric factor of the fluid
Checking of input temperature now yields errors for bad temperatures below fluid min temp
Fixed T(h,p) and T(s,p) in two-phase region 
Fixed Units on surface tension to N/m

2.0.6 (revision 147)
Fixed entropy of humid air at above-atmospheric pressure (Typo in RP-1485)
Added specific heat of humid air
Changes to setup.py so that it will not build if cython version < 0.17 which is a requirement due to the use of STL containers
Changes to plot module to allow for showing right after plot

2.0.5 (revision 143)
Fixed wetbulb and dewpoint calculations - works correctly now
Added wrappers for MATLAB and Octave to subversion code - not included in source distro

2.0.4 (revision 132)
Fixed density for subcooled liquid
Fixed setup.py for OSX (I think)
Using cython for wrapping of CoolProp module
CoolProp module - removed T_hp and h_sp - use Props instead
Added IceProps function to HumidAirProps
Added and fixed CO2 transport properties

2.0.1 (revision 122)
Implemented the method of Akasaka to calculate the saturation state (works great).  H/T to FPROPS for the recommendation
Fixed the calculations for T(h,p) up to a subcooling of 50 K, works fine in superheated vapor
Added the ideal-gas specific heat with key of C0

2.0.0 (revision 107)
MAJOR revision to the internals of CoolProp
Entropy added for humid air (Only fully validated at atmospheric pressure)
Added the fluids R22, R1234yf and the 20 industrial fluids from Lemmon, 2000
Added ECS model for calculation of transport properties (somewhat experimental)
Added surface tension for all fluids.  Property key is 'I' for surface tension
Some functions have been removed in order to better handle errors at the C++ level.  
    Tcrit(), Tsat() and pcrit() are gone, in Python call Props('R134a','Tcrit') for instance to get Tcrit
Many other bug fixes.
Documentation to follow.

1.4.0 (revision 75)
Internal codebase rewritten in C++ to allow for better exception handling and function overloading
All work now happens in CoolProp.cpp (inspired by FPROPS)
Added 2-D lookup table (temperature and pressure) directly in CoolProp.  Enable by calling UseSinglePhaseLUT(1) to turn on, UseSinglePhaseLUT(0) to turn off
Compiled with the -builtin compilation flag
Documentation updated for UseSinglePhaseLUT

1.3.2 (revision 49)
Added functions to use Isothermal compressibility correlation UseIsothermCompressCorrelation and ideal gas compressibility UseIdealGasEnthalpyCorrelations

1.3.1 (revision 48)
Updated documentation
Added ability to use virial term correlations for Humid air by call to UseVirialCorrelation(1)

1.3 (revision 41):
Added pseudo-pure fluid Air using EOS from Lemmon
Added EOS for ice from IAPWS
Updated Humid Air Thermo Props to use analysis from ASHRAE RP-1845, though IAPWS-1995 is used throughout for water vapor
Enable the use of lookup tables for refrigerant saturation properties[ call UseSaturationLUT(1) to turn on, and UseSaturationLUT(0) to turn off]  Speed up is very significant!

1.2.2 (revision 35): 
Added some simple cycles for comparison of different working fluids
Fixed quality calculations to agree with REFPROP