# Example of CoolProp for Python
# Ian Bell, 2013

from __future__ import print_function

import CoolProp
import CoolProp.CoolProp as CP

print('CoolProp version: ', CoolProp.__version__)
print('CoolProp gitrevision: ', CoolProp.__gitrevision__)
print('CoolProp fluids: ', CoolProp.__fluids__)

print(' ')
print('************ USING EOS *************')
print(' ')
print('FLUID STATE INDEPENDENT INPUTS')
print('Critical Density Propane:', CP.PropsSI('Propane', 'rhocrit'), 'kg/m^3')
print('TWO PHASE INPUTS (Pressure)')
print('Density of saturated liquid Propane at 101325 Pa:',
       CP.PropsSI('D', 'P', 101325, 'Q', 0, 'Propane'), 'kg/m^3')
print('Density of saturated vapor R290 at 101325 Pa:',
       CP.PropsSI('D', 'P', 101325, 'Q', 1, 'R290'), 'kg/m^3')
print('TWO PHASE INPUTS (Temperature)')
print('Density of saturated liquid Propane at 300 K:',
       CP.PropsSI('D', 'T', 300, 'Q', 0, 'Propane'), 'kg/m^3')
print('Density of saturated vapor R290 at 300 K:',
       CP.PropsSI('D', 'T', 300, 'Q', 1, 'R290'), 'kg/m^3')

p = CP.PropsSI('P', 'T', 300, 'D', 1, 'Propane')
h = CP.PropsSI('H', 'T', 300, 'D', 1, 'Propane')
T = CP.PropsSI('T', 'P', p, 'H', h, 'Propane')
D = CP.PropsSI('D', 'P', p, 'H', h, 'Propane')
print('SINGLE PHASE CYCLE (propane)')
print('T,D -> P,H', 300, ',', 1, '-->', p, ',', h)
print('P,H -> T,D', p, ',', h, '-->', T, ',', D)


CP.enable_TTSE_LUT('Propane')
print(' ')
print('************ USING TTSE ***************')
print(' ')
print('TWO PHASE INPUTS (Pressure)')
print('Density of saturated liquid Propane at 101325 Pa:',
       CP.PropsSI('D', 'P', 101325, 'Q', 0, 'Propane'), 'kg/m^3')
print('Density of saturated vapor R290 at 101325 Pa:',
       CP.PropsSI('D', 'P', 101325, 'Q', 1, 'R290'), 'kg/m^3')
print('TWO PHASE INPUTS (Temperature)')
print('Density of saturated liquid Propane at 300 K:',
       CP.PropsSI('D', 'T', 300, 'Q', 0, 'Propane'), 'kg/m^3')
print('Density of saturated vapor R290 at 300 K:',
       CP.PropsSI('D', 'T', 300, 'Q', 1, 'R290'), 'kg/m^3')

p = CP.PropsSI('P', 'T', 300, 'D', 1, 'Propane')
h = CP.PropsSI('H', 'T', 300, 'D', 1, 'Propane')
T = CP.PropsSI('T', 'P', p, 'H', h, 'Propane')
D = CP.PropsSI('D', 'P', p, 'H', h, 'Propane')
print('SINGLE PHASE CYCLE (propane)')
print('T,D -> P,H', 300, ',', 1, '-->', p, ',', h)
print('P,H -> T,D', p, ',', h, '-->', T, ',', D)
CP.disable_TTSE_LUT('Propane')

try:
    print(' ')
    print('************ USING REFPROP ***************')
    print(' ')
    print('TWO PHASE INPUTS (Pressure)')
    print('Density of saturated liquid Propane at 101325 Pa:',
           CP.PropsSI('D', 'P', 101325, 'Q', 0, 'REFPROP-Propane'), 'kg/m^3')
    print('Density of saturated vapor Propane at 101325 Pa:',
           CP.PropsSI('D', 'P', 101325, 'Q', 1, 'REFPROP-propane'), 'kg/m^3')
    print('TWO PHASE INPUTS (Temperature)')
    print('Density of saturated liquid Propane at 300 K:',
           CP.PropsSI('D', 'T', 300, 'Q', 0, 'REFPROP-propane'), 'kg/m^3')
    print('Density of saturated vapor Propane at 300 K:',
           CP.PropsSI('D', 'T', 300, 'Q', 1, 'REFPROP-propane'), 'kg/m^3')

    p = CP.PropsSI('P', 'T', 300, 'D', 1, 'Propane')
    h = CP.PropsSI('H', 'T', 300, 'D', 1, 'Propane')
    T = CP.PropsSI('T', 'P', p, 'H', h, 'Propane')
    D = CP.PropsSI('D', 'P', p, 'H', h, 'Propane')
    print('SINGLE PHASE CYCLE (propane)')
    print('T,D -> P,H', 300, ',', 1, '-->', p, ',', h)
    print('P,H -> T,D', p, ',', h, '-->', T, ',', D)
except:
    print(' ')
    print('************ CANT USE REFPROP ************')
    print(' ')

print(' ')
print('************ BRINES AND SECONDARY WORKING FLUIDS *************')
print(' ')
print('Density of 50% (mass) ethylene glycol/water at 300 K, 101325 Pa:',
       CP.PropsSI('D', 'T', 300, 'P', 101325, 'EG-50%'), 'kg/m^3')
print('Viscosity of Therminol D12 at 350 K, 101325 Pa:',
       CP.PropsSI('V', 'T', 350, 'P', 101325, 'TD12'), 'Pa-s')

print(' ')
print('************ HUMID AIR PROPERTIES *************')
print(' ')
print('Humidity ratio of 50% rel. hum. air at 300 K, 101.325 kPa:',
       CP.HAProps('W', 'T', 300, 'P', 101.325, 'R', 0.5), 'kg_w/kg_da')
print('Relative humidity from last calculation:',
       CP.HAProps('R', 'T', 300, 'P', 101.325, 'W',
                   CP.HAProps('W', 'T', 300, 'P', 101.325, 'R', 0.5)),
       '(fractional)')

