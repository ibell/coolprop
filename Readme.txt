Changelog:

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