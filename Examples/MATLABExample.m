addpath('../CoolPropDLL')
loadlibrary('CoolPropDLL','CoolProp_dll.h')
libfunctions CoolPropDLL -full
disp 'The critical temperature of R410A is:'
res=calllib('CoolPropDLL','Tcrit_dll','R410A');
disp(res)
disp 'The saturated vapor enthalpy of Propane at 275 K is:'
h=calllib('CoolPropDLL','Props_dll',int8('H'),int8('T'),275.0,int8('Q'),1.0,'R290');
disp(h)
disp 'The density of nitrogen at STP is:'
rho=calllib('CoolPropDLL','Props_dll',int8('D'),int8('T'),298.0,int8('P'),101.325,'Nitrogen');
disp(rho)
unloadlibrary('CoolPropDLL')