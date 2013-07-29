###########################
Fluid Information
###########################

CoolProp Fluid Information
==========================
.. toctree::
    :maxdepth: 1

    PseudoPureFluids.rst
    PureFluids.rst
    
There are 4 basic classes of fluids that can be used in CoolProp, an example is provided for each one.

Pure Fluids, Pseudo-Pure Fluids
-------------------------------
All the fluids listed in Pure-Fluids and Pseudo-Pure-Fluids sections above can be used.  To use one of these fluids, do something like

.. ipython::

    In [1]: from CoolProp.CoolProp import Props
    
    #Density of dry air at 1 atm. and 25C
    In [1]: Props('D','T',298,'P',101.325,'Air')
    
You can also use any of the aliases of the fluids that are the listed on the fluid page.  For instance, R717 is the refrigerant number for ammonia

.. ipython::

    In [1]: from CoolProp.CoolProp import Props
    
    #Density of saturated ammonia vapor at 1 atm.
    In [1]: Props('D','Q',1,'P',101.325,'R717')
    
    #Density of saturated ammonia vapor at 1 atm.
    In [1]: Props('D','Q',1,'P',101.325,'Ammonia')
    

REFPROP Fluids and Mixtures
---------------------------
If you are on Windows and have REFPROP installed, you can use it with CoolProp.  REFPROP needs to be installed in c:\\Program Files\\REFPROP.  If it is somewhere else, just copy it to this location.

It is also possible to use REFPROP on Linux.  Please follow the instructions from https://github.com/jowr/librefprop.so to install the library from Fortran sources.  Additionally, you also need to copy the fluid and mixture files to /opt/refprop. 

All the pure fluids in REFPROP are used just like the CoolProp fluids except that "REFPROP-" is added at the beginning of the fluid name.  You can use any fluid that is included in REFPROP, but you must use the REFPROP fluid file name.  For CoolProp Fluids, you can use the ``get_REFPROPName()`` function to get the REFPROP name for the fluid.

.. ipython::

    In [1]: from CoolProp.CoolProp import Props
    
    #Saturated isobutane vapor density at 1 atmosphere
    In [1]: Props('D','Q',1,'P',101.325,'REFPROP-ISOBUTAN')
    
You can also use mixtures in REFPROP, there is a special format for the fluid name.  The fluid name is set up like this: ``"REFPROP-MIX:R32[0.697615]&R125[0.302385]"`` -  this is R410A.  The numbers within the brackets are the mole fractions of the components.  They must add up to 1.0
    
.. ipython::

    In [1]: from CoolProp.CoolProp import Props
    
    #Saturated R410 vapor density at 1 atmosphere using the mixture properties
    In [1]: Props('D','Q',1,'P',101.325,'REFPROP-MIX:R32[0.697615]&R125[0.302385]')

Brines
------
A number of aqueous solutions are implemented using the coefficients from Melinder (2010).  The list of diluents implemented are

==========================   ===================================================
Fluid Name                   Description
==========================   ===================================================
``EG``                       Ethylene Glycol
``PG``                       Propylene Glycol
``EA``                       Ethyl Alcohol (Ethanol)
``MA``                       Methyl Alcohol (Methanol)
``Glycerol``                 Glycerol
``NH3``                      Ammonia
``K2CO3``                    Potassium Carbonate
``CaCl2``                    Calcium Chloride
``MgCl2``                    Magnesium Chloride
``NaCl``                     Sodium Chloride
``KAC``                      Potassium Acetate
``KFO``                      Potassium Formate
``LiCl``                     Lithium Chloride
==========================   ===================================================

To use them the fluid name is something like ``"EG-20%"`` which is a 20% by mass ethylene glycol solution

.. ipython::

    In [1]: from CoolProp.CoolProp import Props
    
    #Specific heat 20% mass ethylene glycol solution at 300 K and 1 atm.
    In [1]: Props('C','T',300,'P',101.325,'EG-20%')
    
Incompressible Liquids
----------------------
There is also a selection of incompressible liquids implemented.  These only allow for calls with temperature and pressure as input and provide only a subset of thermophysical properties, namely: density, heat capacity, internal energy, enthalpy, entropy, viscosity and thermal conductivity.

.. ipython::

    In [1]: from CoolProp.CoolProp import Props
    
    #Density of HFE-7100 at 300 K and 1 atm.
    In [1]: Props('D','T',300,'P',101.325,'HFE')
 

For refrigeration applications, 8 fluids were implemented from Melinder 2010 and coefficients are obtained from a fit between -80 and +100 degrees Celsius. The reference point for 0 entropy and 0 internal energy is 25 degrees Celsius. 

==========================   ===================================================
Fluid Name                   Description
==========================   ===================================================
``DEB``                      Diethyl Benzene
``HCM``                      Hydrocarbon Mixture (Therminol D12 Solutia)
``HFE``                      Hydrofluoroether HFE-7100
``PMS1``                     Polydimethylsiloxan 1.
``PMS2``                     Polydimethylsiloxan 2.
``SAB``                      Synthetic alkyl benzene
``HCB``                      Hydrocarbon blend (Dynalene MV)
``TCO``                      Terpene from citrus oils
==========================   ===================================================

There are also a few high temperature heat transfer fluids with individual temperature ranges. Please refer to the file IncompLiquid.h for a complete overview. For these fluids, information from commercial data sheets was used to obtain coefficients.

==========================   ===================================================
Fluid Name                   Description
==========================   ===================================================
``TD12``                     Therminol D12 (-85 to +230 C)
``TVP1``                     Therminol VP-1 (+12 to +397 C)
``T72``                      Therminol 72 (-10 to +380 C)
``DowJ``                     Dowtherm J (-80 to +345 C)
``DowQ``                     Dowtherm Q (-35 to +360 C)
``TX22``                     Texatherm 22 (+0 to +350 C)
``NaK``                      Nitrate Salt Blend (+300 to +600 C)
==========================   ===================================================

All fluids are implemented with polynomials for density and heat capacity with typically 4 coefficients C and hence a third order polynomial. Thermal conductivity is a second order polynomial and viscosity and vapour pressure are exponential functions. 

.. math::

    \rho   = \sum_{i=0}^n C_{\rho}[i] \cdot T^i
    
    c_p    = \sum_{i=0}^n C_{c_p}[i] \cdot T^i
    
    h      = \sum_{i=0}^n \frac{1}{i+1} \cdot C_{c_p}[i] 
             \cdot \left( T^{i+1} - T_0^{i+1} \right)
             
    s      = C_{c_p}[0] \cdot \ln\left(\frac{T}{T_0}\right) 
             + \sum_{i=1}^n \frac{1}{i+1} \cdot C_{c_p}[i] 
             \cdot \left( T^{i+1} - T_0^{i+1} \right)
             
    \lambda= \sum_{i=0}^n C_{\lambda}[i] \cdot T^i
    
    \mu    = \exp\left( \frac{C_{\mu}[0]}{T+C_{\mu}[1]} - C_{\mu}[2] \right)
    
    p_{sat}= \exp\left( \frac{C_{sat}[0]}{T+C_{sat}[1]} - C_{sat}[2] \right)

