addpath('..')
loadlibrary('CoolPropDLL','CoolProp_dll.h')
disp 'The functions available are:'
libfunctions CoolPropDLL -full
disp 'The critical temperature of R410A is:'
res=calllib('CoolPropDLL','Tcrit_dll','R410A');
disp(res)
disp 'The saturated vapor enthalpy of Propane at 275 K is:'
res=calllib('CoolPropDLL','Props_dll',int8('H'),int8('T'),275.0,int8('Q'),1.0,'R290')
disp(res)
unloadlibrary('CoolPropDLL')