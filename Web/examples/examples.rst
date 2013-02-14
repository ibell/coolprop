Examples
========
 
Sample Props Code
-------------------
To use the Props function, import it and do some calls, do something like this

.. ipython::

    #import the things you need 
    In [1]: from CoolProp.CoolProp import Props
    
    # print the currently used version of coolprop
    In [1]: import CoolProp; print CoolProp.__version__
    
    @suppress
    In [1]: from CoolProp.CoolProp import UseSinglePhaseLUT; UseSinglePhaseLUT(0)
    
    #Density of carbon dioxide (R744) at 100 bar and 25C
    In [2]: Props('D','T',298.15,'P',10000,'R744')
    
    #Saturated vapor enthalpy [kJ/kg] of R134a at STP
    In [2]: Props('H','T',298.15,'Q',1,'R134a')

Or go to the :ref:`Fluid-Properties` documentation.

All the possible input and output parameters are listed in the :py:func:`CoolProp.CoolProp.Props` documentation

Sample HAProps Code
-------------------
To use the HAProps function, import it and do some calls, do something like this

.. ipython::

    #import the things you need 
    In [1]: from CoolProp.HumidAirProp import HAProps, HAProps_Aux
    
    #Enthalpy (kJ per kg dry air) as a function of temperature, pressure, 
    #    and relative humidity at STP
    In [2]: h=HAProps('H','T',298.15,'P',101.325,'R',0.5); print h
    
    #Temperature of saturated air at the previous enthalpy
    In [2]: T=HAProps('T','P',101.325,'H',h,'R',1.0); print T
    
    #Temperature of saturated air - order of inputs doesn't matter
    In [2]: T=HAProps('T','H',h,'R',1.0,'P',101.325); print T
    
    
Or go to the :ref:`Humid-Air` documentation.

Plotting
--------

To make a Temperature-entropy plot for propane (R290), you can just do this

.. plot::
    :include-source:
    
    from CoolProp.Plots import Ts
    Ts('R290')
    
or for a pressure-enthalpy plot of R410A

.. plot::
    :include-source:
    
    from CoolProp.Plots import Ph
    Ph('R410A')
    
and to overlay a simple four-component cycle on a R410A P-h plot

.. plot::
    :include-source:
    
    from CoolProp.Plots import Ph,SimpleCycle
    Ph('R410A')
    SimpleCycle('R410A',250,300,5,5,0.7)
    
