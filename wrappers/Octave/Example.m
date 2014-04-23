% Example of CoolProp for Octave
% Ian Bell, 2013

CoolProp
disp(['CoolProp version: ', CoolProp.get_global_param_string('version')])
disp(['CoolProp gitrevision: ', CoolProp.get_global_param_string('gitrevision')])
disp(['CoolProp fluids: ', CoolProp.get_global_param_string('FluidsList')])

disp(' ')
disp('************ USING EOS *************')
disp(' ')
disp('FLUID STATE INDEPENDENT INPUTS')
sprintf('Critical Density Propane: %g kg/m^3', CoolProp.Props1SI('Propane','rhocrit'))
disp(['TWO PHASE INPUTS (Pressure)'])
sprintf('Density of saturated liquid Propane at 101325 Pa: %g kg/m^3', CoolProp.PropsSI('D','P',101325,'Q',0,'Propane'))
sprintf('Density of saturated vapor R290 at 101325 Pa: %g kg/m^3', CoolProp.PropsSI('D','P',101325,'Q',1,'R290'))
disp(['TWO PHASE INPUTS (Temperature)'])
sprintf('Density of saturated liquid Propane at 300 K: %g kg/m^3', CoolProp.PropsSI('D','T',300,'Q',0,'Propane'))
sprintf('Density of saturated vapor R290 at 300 K: %g kg/m^3', CoolProp.PropsSI('D','T',300,'Q',1,'R290'))
disp(['SINGLE PHASE CYCLE (propane)'])
p = CoolProp.PropsSI('P','T',300,'D',1,'Propane'); 
h = CoolProp.PropsSI('H','T',300,'D',1,'Propane');
sprintf('T,D -> P,H %g,%g --> %g,%g', 300, 1, p, h)
T = CoolProp.PropsSI('T','P',p,'H',h,'Propane');
D = CoolProp.PropsSI('D','P',p,'H',h,'Propane');
sprintf('P,H -> T,D %g,%g --> %g,%g', p, h, T, D)
 
disp([' '])
disp(['************ USING TTSE ***************'])
disp([' '])
CoolProp.enable_TTSE_LUT('Propane');
disp(['TWO PHASE INPUTS (Pressure)'])
sprintf('Density of saturated liquid Propane at 101325 Pa: %g kg/m^3', CoolProp.PropsSI('D','P',101325,'Q',0,'Propane'))
sprintf('Density of saturated vapor R290 at 101325 Pa: %g kg/m^3', CoolProp.PropsSI('D','P',101325,'Q',1,'R290'))
disp(['TWO PHASE INPUTS (Temperature)'])
sprintf('Density of saturated liquid Propane at 300 K: %g kg/m^3', CoolProp.PropsSI('D','T',300,'Q',0,'Propane'))
sprintf('Density of saturated vapor R290 at 300 K: %g kg/m^3', CoolProp.PropsSI('D','T',300,'Q',1,'R290'))
disp(['SINGLE PHASE CYCLE (propane)'])
p = CoolProp.PropsSI('P','T',300,'D',1,'Propane'); 
h = CoolProp.PropsSI('H','T',300,'D',1,'Propane');
sprintf('T,D -> P,H %g,%g --> %g,%g', 300, 1, p, h)
T = CoolProp.PropsSI('T','P',p,'H',h,'Propane');
D = CoolProp.PropsSI('D','P',p,'H',h,'Propane');
sprintf('P,H -> T,D %g,%g --> %g,%g', p, h, T, D)
CoolProp.disable_TTSE_LUT('Propane');

try
    disp(' ')
    disp('************ USING REFPROP ***************')
    disp(' ')
    disp('FLUID STATE INDEPENDENT INPUTS')
    disp(['Critical Density Propane:', CoolProp.PropsSI('REFPROP-Propane','rhocrit'), 'kg/m^3'])
    disp(['TWO PHASE INPUTS (Pressure)'])
    sprintf('Density of saturated liquid Propane at 101325 Pa: %g kg/m^3', CoolProp.PropsSI('D','P',101325,'Q',0,'REFPROP-Propane'))
    sprintf('Density of saturated vapor R290 at 101325 Pa: %g kg/m^3', CoolProp.PropsSI('D','P',101325,'Q',1,'REFPROP-Propane'))
    disp(['TWO PHASE INPUTS (Temperature)'])
    sprintf('Density of saturated liquid Propane at 300 K: %g kg/m^3', CoolProp.PropsSI('D','T',300,'Q',0,'REFPROP-Propane'))
    sprintf('Density of saturated vapor R290 at 300 K: %g kg/m^3', CoolProp.PropsSI('D','T',300,'Q',1,'REFPROP-Propane'))
    disp(['SINGLE PHASE CYCLE (propane)'])
    p = CoolProp.PropsSI('P','T',300,'D',1,'REFPROP-Propane'); 
    h = CoolProp.PropsSI('H','T',300,'D',1,'REFPROP-Propane');
    sprintf('T,D -> P,H %g,%g --> %g,%g', 300, 1, p, h)
    T = CoolProp.PropsSI('T','P',p,'H',h,'REFPROP-Propane');
    D = CoolProp.PropsSI('D','P',p,'H',h,'REFPROP-Propane');
    sprintf('P,H -> T,D %g,%g --> %g,%g', p, h, T, D)
catch
    disp(' ')
    disp('************ CANT USE REFPROP ************')
    disp(' ')
end
 
disp(' ')
disp('************ BRINES AND SECONDARY WORKING FLUIDS *************')
disp(' ')
sprintf('Density of 50%% (mass) ethylene glycol/water at 300 K, 101325 Pa: %g kg/m^3', CoolProp.PropsSI('D','T',300,'P',101325,'EG-50%'))
sprintf('Viscosity of Therminol D12 at 350 K, 101325 Pa: %g Pa-s', CoolProp.PropsSI('V', 'T', 350, 'P', 101325, 'TD12'))

disp(' ')
disp('************ HUMID AIR PROPERTIES *************')
disp(' ')
sprintf('Humidity ratio of 50%% rel. hum. air at 300 K, 101325 kPa: %g kg_w/kg_da', CoolProp.HAProps('W','T',300,'P',101.325,'R',0.5))
sprintf('Relative humidity from last calculation: %g (fractional)', CoolProp.HAProps('R','T',300,'P',101.325,'W',HAProps('W','T',300,'P',101.325,'R',0.5)))
