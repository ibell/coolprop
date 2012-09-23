
********************
Air
********************

Equation of State Reference
===========================
Lemmon, E.W., Jacobsen, R.T, Penoncello, S.G., and Friend, D.G.,"Thermodynamic Properties of Air and Mixtures of Nitrogen, Argon, and Oxygen from 60 to 2000 K at Pressures to 2000 MPa, J. Phys. Chem. Ref. Data, 29(3):331-385, 2000.

Transport Properties Information
================================
E.W. Lemmon and R. T. Jacobsen, "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air", International Journal of Thermophysics, Vol. 25, No. 1, January 2004


Fluid Data
==========

Fluid Parameters

=========================  ==============================
Mole Mass [kg/kmol]        28.96546
Triple Point [K]           59.750
=========================  ==============================

Critical Parameters

==========================  ==============================
Temperature [K]             132.53
Density [kg/m\ :sup:`3`\ ]   342.684564
Pressure [kPa]              3786.00000
==========================  ==============================


Saturated Vapor Deviations
==========================

.. plot::

    Fluid = "Air"
    RPFluid = "REFPROP-Air"

    #Saturated Vapor
    from CoolProp.CoolProp import Props
    from numpy import linspace,array,abs
    import matplotlib.pyplot as plt

    Tt = Props(Fluid,'Ttriple')
    Tc = Props(Fluid,'Tcrit')
    Tv = linspace(Tt+0.01,0.95*Tc,20)

    #All the CoolProp calculations
    p = array([Props('P','T',T,'Q',1,Fluid) for T in Tv])
    rho = array([Props('D','T',T,'Q',1,Fluid) for T in Tv])
    cp = array([Props('C','T',T,'Q',1,Fluid) for T in Tv])
    cv = array([Props('O','T',T,'Q',1,Fluid) for T in Tv])
    sigma = array([Props('I','T',T,'Q',1,Fluid) for T in Tv])

    Rp = array([Props('P','T',T,'Q',1,RPFluid) for T in Tv])
    Rrho = array([Props('D','T',T,'Q',1,RPFluid) for T in Tv])
    Rcp = array([Props('C','T',T,'Q',1,RPFluid) for T in Tv])
    Rcv = array([Props('O','T',T,'Q',1,RPFluid) for T in Tv])
    Rsigma = array([Props('I','T',T,'Q',1,RPFluid) for T in Tv])

    fig = plt.figure()

    ax = fig.add_axes((0.15,0.15,0.8,0.8))
    ax.semilogy(Tv/Tc,abs(p/Rp-1)*100,'o',label='Pressure')
    ax.semilogy(Tv/Tc,abs(rho/Rrho-1)*100,'o',label='Density')
    ax.semilogy(Tv/Tc,abs(cp/Rcp-1)*100,'o',label='Specific heat (cp)')
    ax.semilogy(Tv/Tc,abs(cv/Rcv-1)*100,'o',label='Specific heat (cv)')
    ax.semilogy(Tv/Tc,abs(sigma/Rsigma-1)*100,'o',label='Surface tension')
    ax.set_ylim(1e-16,100)
    ax.set_title('Saturated Vapor Deviations from REFPROP 9.0')
    ax.set_xlabel('Reduced temperature T/Tc')
    ax.set_ylabel('Absolute deviation [%]')
    ax.legend(numpoints=1,loc='best')
    plt.show()

Saturated Liquid Deviations
===========================

.. plot::

    Fluid = "Air"
    RPFluid = "REFPROP-Air"

    #Saturated Liquid
    from CoolProp.CoolProp import Props
    from numpy import linspace,array,abs
    import matplotlib.pyplot as plt

    Tt = Props(Fluid,'Ttriple')
    Tc = Props(Fluid,'Tcrit')
    Tv = linspace(Tt+0.01,0.95*Tc,20)

    #All the CoolProp calculations
    p = array([Props('P','T',T,'Q',0,Fluid) for T in Tv])
    rho = array([Props('D','T',T,'Q',0,Fluid) for T in Tv])
    cp = array([Props('C','T',T,'Q',0,Fluid) for T in Tv])
    cv = array([Props('O','T',T,'Q',0,Fluid) for T in Tv])
    sigma = array([Props('I','T',T,'Q',0,Fluid) for T in Tv])

    Rp = array([Props('P','T',T,'Q',0,RPFluid) for T in Tv])
    Rrho = array([Props('D','T',T,'Q',0,RPFluid) for T in Tv])
    Rcp = array([Props('C','T',T,'Q',0,RPFluid) for T in Tv])
    Rcv = array([Props('O','T',T,'Q',0,RPFluid) for T in Tv])
    Rsigma = array([Props('I','T',T,'Q',0,RPFluid) for T in Tv])

    fig = plt.figure()

    ax = fig.add_axes((0.15,0.15,0.8,0.8))
    ax.semilogy(Tv/Tc,abs(p/Rp-1)*100,'o',label='Pressure')
    ax.semilogy(Tv/Tc,abs(rho/Rrho-1)*100,'o',label='Density')
    ax.semilogy(Tv/Tc,abs(cp/Rcp-1)*100,'o',label='Specific heat (cp)')
    ax.semilogy(Tv/Tc,abs(cv/Rcv-1)*100,'o',label='Specific heat (cv)')
    ax.semilogy(Tv/Tc,abs(sigma/Rsigma-1)*100,'o',label='Surface tension')
    ax.set_ylim(1e-16,100)
    ax.set_title('Saturated Liquid Deviations from REFPROP 9.0')
    ax.set_xlabel('Reduced temperature T/Tc')
    ax.set_ylabel('Absolute deviation [%]')
    ax.legend(numpoints=1,loc='best')
    plt.show()

Along the critical isotherm where T=T\ :sub:`c`
================================================
.. plot::

    Fluid = "Air"
    RPFluid = "REFPROP-Air"

    #Critical isotherm
    from CoolProp.CoolProp import Props
    from numpy import linspace,array,abs
    import matplotlib.pyplot as plt

    Tc = Props(Fluid,'Tcrit')
    rhoc = Props(Fluid,'rhocrit')
    rhov = linspace(1e-12,2*rhoc)

    #All the CoolProp calculations
    p = array([Props('P','T',Tc,'D',D,Fluid) for D in rhov])
    rho = array([Props('D','T',Tc,'D',D,Fluid) for D in rhov])
    cp = array([Props('C','T',Tc,'D',D,Fluid) for D in rhov])
    cv = array([Props('O','T',Tc,'D',D,Fluid) for D in rhov])
    visc = array([Props('V','T',Tc,'D',D,Fluid) for D in rhov])

    Rp = array([Props('P','T',Tc,'D',D,RPFluid) for D in rhov])
    Rrho = array([Props('D','T',Tc,'D',D,RPFluid) for D in rhov])
    Rcp = array([Props('C','T',Tc,'D',D,RPFluid) for D in rhov])
    Rcv = array([Props('O','T',Tc,'D',D,RPFluid) for D in rhov])
    Rvisc = array([Props('V','T',Tc,'D',D,RPFluid) for D in rhov])

    fig = plt.figure()

    ax = fig.add_axes((0.15,0.15,0.8,0.8))
    ax.semilogy(rhov/rhoc,abs(p/Rp-1)*100,'o',label='Pressure')
    ax.semilogy(rhov/rhoc,abs(cp/Rcp-1)*100,'o',label='Specific heat (cp)')
    ax.semilogy(rhov/rhoc,abs(cv/Rcv-1)*100,'o',label='Specific heat (cv)')
    ax.semilogy(rhov/rhoc,abs(visc/Rvisc-1)*100,'o',label='Viscosity')
    ax.set_ylim(1e-16,100)
    ax.set_title('Critical isotherm Deviations from REFPROP 9.0')
    ax.set_xlabel(r'Reduced density $\rho/\rho_c$')
    ax.set_ylabel('Absolute deviation [%]')
    ax.legend(numpoints=1,loc='best')
    plt.show()
