CoolProp is a property database and wrappers for a selection of programming environments

Changelog:
3.1 (revision 534)
Added the fluids Propylene, Cyclopentane, R236FA, R236EA, R227EA, R123, R152A, R227EA, R365MFC, R161, HFE143M, Benzene, R11, Undecane, R125, Cyclopropane, Neon, R124
Added the viscosity and conductivity correlations for a lot of fluids
Added surface tension, Lennard-Jones parameters for a lot of fluids
Added enthalpy, entropy as inputs
Added pressure, density as inputs
CoolProp builds on Raspberry PI
CoolProp works in MATLAB on OSX
Python unit tests have been added in wrappers/CoolProp/CoolProp/tests - a work in progress

3.0 (revision 325)
Added Tabular Taylor Series Expansion (see documentation)
All the way to the critical point for almost all fluids
Support added for Modelica, Python 3.x and Labview

2.5 (revision 247)
Added EES wrapper (r245-r247)
Saturation derivatives dhdp and d2hdp2 (r244)
Caching of Helholtz derivatives in CPState.cpp (r243)
Added Xylenes and EthylBenzene (r242)
Added n-Dodecane, R23, DMC (r241)


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

