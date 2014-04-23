% Example of CoolProp for MATLAB
% Ian Bell, 2013

disp(['CoolProp version: ', PropsSI('version')])
disp(['CoolProp gitrevision: ', PropsSI('gitrevision')])
disp(['CoolProp fluids: ', PropsSI('FluidsList')])

disp(' ')
disp('************ USING EOS *************')
disp(' ')
disp('FLUID STATE INDEPENDENT INPUTS')
disp(['Critical Density Propane: ', num2str(PropsSI('Propane','rhocrit')), ' kg/m^3'])
disp(['TWO PHASE INPUTS (Pressure)'])
disp(['Density of saturated liquid Propane at 101325 Pa: ', num2str(PropsSI('D','P',101325,'Q',0,'Propane')), ' kg/m^3'])
disp(['Density of saturated vapor R290 at 101325 Pa: ', num2str(PropsSI('D','P',101325,'Q',1,'R290')), ' kg/m^3'])
disp(['TWO PHASE INPUTS (Temperature)'])
disp(['Density of saturated liquid Propane at 300 K: ', num2str(PropsSI('D','T',300,'Q',0,'Propane')), ' kg/m^3'])
disp(['Density of saturated vapor R290 at 300 K: ', num2str(PropsSI('D','T',300,'Q',1,'R290')), ' kg/m^3'])
disp(['SINGLE PHASE CYCLE (propane)'])
p = PropsSI('P','T',300,'D',1,'Propane'); 
h = PropsSI('H','T',300,'D',1,'Propane');
disp(['T,D -> P,H ', num2str(300),',',num2str(1), ' --> ',num2str(p),',',num2str(h)])
T = PropsSI('T','P',p,'H',h,'Propane'); 
D = PropsSI('D','P',p,'H',h,'Propane');
disp(['P,H -> T,D', num2str(p),',',num2str(h),'-->',num2str(T),',',num2str(D)])
 
disp([' '])
disp(['************ USING TTSE ***************'])
disp([' '])
PropsSI('Propane','enable_TTSE')
disp(['TWO PHASE INPUTS (Pressure)'])
disp(['Density of saturated liquid Propane at 101325 Pa: ', num2str(PropsSI('D','P',101325,'Q',0,'Propane')), ' kg/m^3'])
disp(['Density of saturated vapor R290 at 101325 Pa: ', num2str(PropsSI('D','P',101325,'Q',1,'R290')), ' kg/m^3'])
disp(['TWO PHASE INPUTS (Temperature)'])
disp(['Density of saturated liquid Propane at 300 K: ', num2str(PropsSI('D','T',300,'Q',0,'Propane')), ' kg/m^3'])
disp(['Density of saturated vapor R290 at 300 K: ', num2str(PropsSI('D','T',300,'Q',1,'R290')), ' kg/m^3'])
disp(['SINGLE PHASE CYCLE (propane)'])
p = PropsSI('P','T',300,'D',1,'Propane'); 
h = PropsSI('H','T',300,'D',1,'Propane');
disp(['T,D -> P,H ', num2str(300),',',num2str(1), ' --> ',num2str(p),',',num2str(h)])
T = PropsSI('T','P',p,'H',h,'Propane'); 
D = PropsSI('D','P',p,'H',h,'Propane');
disp(['P,H -> T,D ', num2str(p),',',num2str(h),' --> ',num2str(T),',',num2str(D)])
PropsSI('Propane','disable_TTSE')

try
    disp(' ')
    disp('************ USING REFPROP ***************')
    disp(' ')
    disp('FLUID STATE INDEPENDENT INPUTS')
    disp(['Critical Density Propane:', num2str(PropsSI('REFPROP-Propane','rhocrit')), ' kg/m^3'])
    disp(['TWO PHASE INPUTS (Pressure)'])
    disp(['Density of saturated liquid Propane at 101325 Pa: ', num2str(PropsSI('D','P',101325,'Q',0,'REFPROP-Propane')), ' kg/m^3'])
    disp(['Density of saturated vapor R290 at 101325 Pa: ', num2str(PropsSI('D','P',101325,'Q',1,'REFPROP-Propane')), ' kg/m^3'])
    disp(['TWO PHASE INPUTS (Temperature)'])
    disp(['Density of saturated liquid Propane at 300 K: ', num2str(PropsSI('D','T',300,'Q',0,'REFPROP-Propane')), ' kg/m^3'])
    disp(['Density of saturated vapor R290 at 300 K: ', num2str(PropsSI('D','T',300,'Q',1,'REFPROP-Propane')), ' kg/m^3'])
    disp(['SINGLE PHASE CYCLE (propane)'])
    p = PropsSI('P','T',300,'D',1,'REFPROP-Propane');
    h = PropsSI('H','T',300,'D',1,'REFPROP-Propane');
    disp(['T,D -> P,H ', num2str(300),',',num2str(1), ' --> ',num2str(p),',',num2str(h)])
    T = PropsSI('T','P',p,'H',h,'REFPROP-Propane'); 
    D = PropsSI('D','P',p,'H',h,'REFPROP-Propane');
    disp(['P,H -> T,D ', num2str(p),',',num2str(h),' --> ',num2str(T),',',num2str(D)])
catch
    disp(' ')
    disp('************ CANT USE REFPROP ************')
    disp(' ')
end
 
disp(' ')
disp('************ BRINES AND SECONDARY WORKING FLUIDS *************')
disp(' ')
disp(['Density of 50% (mass) ethylene glycol/water at 300 K, 101325 Pa: ', num2str(PropsSI('D','T',300,'P',101325,'EG-50%')), 'kg/m^3'])
disp(['Viscosity of Therminol D12 at 350 K, 101325 Pa: ', num2str(PropsSI('V', 'T', 350, 'P', 101325, 'TD12')), 'Pa-s'])

disp(' ')
disp('************ HUMID AIR PROPERTIES *************')
disp(' ')
disp(['Humidity ratio of 50% rel. hum. air at 300 K, 101.325 kPa: ', num2str(HAProps('W','T',300,'P',101.325,'R',0.5)), ' kg_w/kg_da'])
disp(['Relative humidity from last calculation: ', num2str(HAProps('R','T',300,'P',101.325,'W',HAProps('W','T',300,'P',101.325,'R',0.5))), ' (fractional)'])
