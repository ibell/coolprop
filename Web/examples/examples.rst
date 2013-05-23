Examples
========
 
The following examples are written in Python to demonstrate some of the 
functionalities of CoolProp.  Similar calling conventions are used in the wrappers
for other programming languages

Sample Props Code
-------------------
To use the Props function, import it and do some calls, do something like this

.. ipython::

    #import the things you need 
    In [1]: from CoolProp.CoolProp import Props
    
    # print the currently used version of coolprop
    In [1]: import CoolProp; print(CoolProp.__version__)
    
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
    
    from CoolProp.Plots.Plots import Ts
    Ts('R290')
    
or for a pressure-enthalpy plot of R410A

.. plot::
    :include-source:
    
    from CoolProp.Plots.Plots import Ph
    Ph('R410A')
    
and to overlay a simple four-component cycle on a R410A P-h plot.

.. plot::
    :include-source:
    
    from CoolProp.Plots.Plots import Ph,SimpleCycle
    Ph('R410A')
    SimpleCycle('R410A',250,300,5,5,0.7)

A more advanced example using built-in functions to draw lines of constant properties
is given below. Note the different ways to invoke drawIsoLines:
    
.. plot::
    :include-source:
     
    from CoolProp.Plots.Plots import Ts,drawIsoLines
    Ref = 'n-Pentane'
    ax = Ts(Ref)
    ax.set_xlim([-0.5,1.5])
    ax.set_ylim([300,530])
    quality    = drawIsoLines(Ref, 'Ts', 'Q', [0.3,  0.5, 0.7, 0.8] , axis=ax)
    isobars    = drawIsoLines(Ref, 'Ts', 'P', [100, 2000]    , num=5, axis=ax)
    isochores  = drawIsoLines(Ref, 'Ts', 'D', [2,    600]    , num=7, axis=ax)
    #    isenthalps = drawIsoLines(Ref, 'Ts', 'H', [100, 300], num=5, axis=ax)
