
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
    
    from CoolProp.Plots.Plots import Ph
    from CoolProp.Plots.SimpleCycles import SimpleCycle
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
