Changelog:

1.3 (revision 40):
Added pseudo-pure fluid Air using EOS from Lemmon
Added EOS for ice from IAPWS
Updated Humid Air Thermo Props to use analysis from ASHRAE RP-1845, though IAPWS-1995 is used throughout for water vapor
Enable the use of lookup tables for refrigerant saturation properties[ call UseSaturationLUT(1) to turn on, and UseSaturationLUT(0) to turn off]  Speed up is very significant!

1.2.2 (revision 35): 
Added some simple cycles for comparison of different working fluids
Fixed quality calculations to agree with REFPROP