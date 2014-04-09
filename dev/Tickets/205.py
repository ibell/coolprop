import CoolProp.CoolProp as CP                                                                                                                                              
fluid = "n-Pentane"                                                                                                                                                                                
CP.enable_TTSE_LUT(fluid)                                                                                                                                                                                
print CP.PropsSI('H','T',477.8,'D',31.77945717536664,fluid)
print CP.set_debug_level(1000)
print CP.PropsSI('H','T',477.8,'P',1.4700000000000000e+06,fluid)
