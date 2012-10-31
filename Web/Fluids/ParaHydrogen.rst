
********************
ParaHydrogen
********************

Equation of State Reference
===========================
"Fundamental Equations of State for Parahydrogen, Normal Hydrogen, and Orthohydrogen"by J.W. Leachman and R.T. Jacobsen and S.G. Penoncello and E.W. Lemmon, J. Phys. Chem. Ref. Data, Vol. 38, No. 3, 2009, pp 721-748

Transport Properties Information
================================
Viscosity & Surface Tension: McCarty, R.D. and Weber, L.A., "Thermophysical properties of parahydrogen from the freezing liquid line to 5000 R for pressures to 10,000 psia," Natl. Bur. Stand., Tech. Note 617, 1972.


Fluid Data
==========

Fluid Parameters

=========================  ==============================
Mole Mass [kg/kmol]        2.01588
Triple Point [K]           13.803
=========================  ==============================

Critical Parameters

==========================  ==============================
Temperature [K]             32.94
Density [kg/m\ :sup:`3`\ ]   31.322743
Pressure [kPa]              1285.80000
==========================  ==============================


Saturated Vapor Deviations
==========================

.. plot::

    Fluid = "ParaHydrogen"
    RPFluid = "REFPROP-PARAHYD"

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
    h0 = Props('H','T',(Tt+Tc)/2.0,'Q',1,Fluid)
    h = array([Props('H','T',T,'Q',1,Fluid)-h0 for T in Tv])
    s0 = Props('S','T',(Tt+Tc)/2.0,'Q',1,Fluid)
    s = array([Props('S','T',T,'Q',1,Fluid)-s0 for T in Tv])   
    sigma = array([Props('I','T',T,'Q',1,Fluid) for T in Tv])

    Rp = array([Props('P','T',T,'Q',1,RPFluid) for T in Tv])
    Rrho = array([Props('D','T',T,'Q',1,RPFluid) for T in Tv])
    Rcp = array([Props('C','T',T,'Q',1,RPFluid) for T in Tv])
    Rcv = array([Props('O','T',T,'Q',1,RPFluid) for T in Tv])
    Rh0 = Props('H','T',(Tt+Tc)/2.0,'Q',1,RPFluid)
    Rh = array([Props('H','T',T,'Q',1,RPFluid)-Rh0 for T in Tv])
    Rs0 = Props('S','T',(Tt+Tc)/2.0,'Q',1,RPFluid)
    Rs = array([Props('S','T',T,'Q',1,RPFluid)-Rs0 for T in Tv])
    Rsigma = array([Props('I','T',T,'Q',1,RPFluid) for T in Tv])

    fig = plt.figure()

    ax = fig.add_axes((0.15,0.15,0.8,0.8))
    ax.semilogy(Tv/Tc,abs(p/Rp-1)*100,'o',label='Pressure')
    ax.semilogy(Tv/Tc,abs(rho/Rrho-1)*100,'o',label='Density')
    ax.semilogy(Tv/Tc,abs(cp/Rcp-1)*100,'o',label='Specific heat (cp)')
    ax.semilogy(Tv/Tc,abs(cv/Rcv-1)*100,'o',label='Specific heat (cv)')
    ax.semilogy(Tv/Tc,abs(h/Rh-1)*100,'o',label='Enthalpy')
    ax.semilogy(Tv/Tc,abs(s/Rs-1)*100,'o',label='Entropy')  
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

    Fluid = "ParaHydrogen"
    RPFluid = "REFPROP-PARAHYD"

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
    h0 = Props('H','T',(Tt+Tc)/2.0,'Q',0,Fluid)
    h = array([Props('H','T',T,'Q',0,Fluid)-h0 for T in Tv])
    s0 = Props('S','T',(Tt+Tc)/2.0,'Q',0,Fluid)
    s = array([Props('S','T',T,'Q',0,Fluid)-s0 for T in Tv])    
    sigma = array([Props('I','T',T,'Q',0,Fluid) for T in Tv])

    Rp = array([Props('P','T',T,'Q',0,RPFluid) for T in Tv])
    Rrho = array([Props('D','T',T,'Q',0,RPFluid) for T in Tv])
    Rcp = array([Props('C','T',T,'Q',0,RPFluid) for T in Tv])
    Rcv = array([Props('O','T',T,'Q',0,RPFluid) for T in Tv])
    Rh0 = Props('H','T',(Tt+Tc)/2.0,'Q',0,RPFluid)
    Rh = array([Props('H','T',T,'Q',0,RPFluid)-Rh0 for T in Tv])
    Rs0 = Props('S','T',(Tt+Tc)/2.0,'Q',0,RPFluid)
    Rs = array([Props('S','T',T,'Q',0,RPFluid)-Rs0 for T in Tv])
    Rsigma = array([Props('I','T',T,'Q',0,RPFluid) for T in Tv])

    fig = plt.figure()

    ax = fig.add_axes((0.15,0.15,0.8,0.8))
    ax.semilogy(Tv/Tc,abs(p/Rp-1)*100,'o',label='Pressure')
    ax.semilogy(Tv/Tc,abs(rho/Rrho-1)*100,'o',label='Density')
    ax.semilogy(Tv/Tc,abs(cp/Rcp-1)*100,'o',label='Specific heat (cp)')
    ax.semilogy(Tv/Tc,abs(cv/Rcv-1)*100,'o',label='Specific heat (cv)')
    ax.semilogy(Tv/Tc,abs(h/Rh-1)*100,'o',label='Enthalpy')
    ax.semilogy(Tv/Tc,abs(s/Rs-1)*100,'o',label='Entropy')
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

    Fluid = "ParaHydrogen"
    RPFluid = "REFPROP-PARAHYD"

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
    h0 = Props('H','T',0.95*Tc,'Q',1,Fluid)
    h = array([Props('H','T',Tc,'D',D,Fluid)-h0 for D in rhov])
    s0 = Props('S','T',0.95*Tc,'Q',1,Fluid)
    s = array([Props('S','T',Tc,'D',D,Fluid)-s0 for D in rhov])
    visc = array([Props('V','T',Tc,'D',D,Fluid) for D in rhov])

    Rp = array([Props('P','T',Tc,'D',D,RPFluid) for D in rhov])
    Rrho = array([Props('D','T',Tc,'D',D,RPFluid) for D in rhov])
    Rcp = array([Props('C','T',Tc,'D',D,RPFluid) for D in rhov])
    Rcv = array([Props('O','T',Tc,'D',D,RPFluid) for D in rhov])
    Rh0 = Props('H','T',0.95*Tc,'Q',1,RPFluid)
    Rh = array([Props('H','T',Tc,'D',D,RPFluid)-Rh0 for D in rhov])
    Rs0 = Props('S','T',0.95*Tc,'Q',1,RPFluid)
    Rs = array([Props('S','T',Tc,'D',D,RPFluid)-Rs0 for D in rhov])
    Rvisc = array([Props('V','T',Tc,'D',D,RPFluid) for D in rhov])

    fig = plt.figure()

    ax = fig.add_axes((0.15,0.15,0.8,0.8))
    ax.semilogy(rhov/rhoc,abs(p/Rp-1)*100,'o',label='Pressure')
    ax.semilogy(rhov/rhoc,abs(cp/Rcp-1)*100,'o',label='Specific heat (cp)')
    ax.semilogy(rhov/rhoc,abs(cv/Rcv-1)*100,'o',label='Specific heat (cv)')
    ax.semilogy(rhov/rhoc,abs(h/Rh-1)*100,'o',label='Enthalpy')
    ax.semilogy(rhov/rhoc,abs(s/Rs-1)*100,'o',label='Entropy') 
    ax.semilogy(rhov/rhoc,abs(visc/Rvisc-1)*100,'o',label='Viscosity')
    ax.set_ylim(1e-16,100)
    ax.set_title('Critical isotherm Deviations from REFPROP 9.0')
    ax.set_xlabel(r'Reduced density $\rho/\rho_c$')
    ax.set_ylabel('Absolute deviation [%]')
    ax.legend(numpoints=1,loc='best')
    plt.show()

Along the isobar corresponding to :math:`T_{sat}=(T_t+T_c)/2`
=================================================================
.. plot::

    Fluid = "ParaHydrogen"

    from CoolProp.CoolProp import Props
    from numpy import linspace,array,abs
    import matplotlib.pyplot as plt

    Tt = Props(Fluid,'Ttriple')
    Tc = Props(Fluid,'Tcrit')
    Tm = (Tt+Tc)/2.0
    p = Props('P','T',Tm,'Q',0,Fluid)

    fig = plt.figure()

    Tmin = Tm - 30 if Tm - 30 > Tt else Tt+0.1
    Tv = linspace(Tmin, Tm - 0.001, 20)
    #Start with T&P as inputs
    rho = array([Props('D','T',T,'P',p,Fluid) for T in Tv])
    h = array([Props('H','T',T,'P',p,Fluid) for T in Tv])
    s = array([Props('S','T',T,'P',p,Fluid) for T in Tv])
    p_Trho = array([Props('P','T',T,'D',D,Fluid) for T,D in zip(Tv,rho)])
    T_hp = array([Props('T','H',H,'P',p,Fluid) for H in h])
    ax1 = fig.add_subplot(121)
    ax1.semilogy(Tv,abs(p/p_Trho-1)*100,'o',label='p & p(T,rho)')
    ax1.semilogy(Tv,abs(T_hp/Tv-1)*100,'o',label='T & T(h,p)')
    ax1.set_ylim(1e-16,100)
    ax1.set_title('Subcooled liquid')
    ax1.set_xlabel('Temperature [K]')
    ax1.set_ylabel('Absolute deviation [%]')
    ax1.legend(numpoints=1,loc='best')

    Tv = linspace(Tm+0.001, Tm + 30, 20)
    #Start with T&P as inputs
    rho = array([Props('D','T',T,'P',p,Fluid) for T in Tv])
    h = array([Props('H','T',T,'P',p,Fluid) for T in Tv])
    s = array([Props('S','T',T,'P',p,Fluid) for T in Tv])
    p_Trho = array([Props('P','T',T,'D',D,Fluid) for T,D in zip(Tv,rho)])
    T_hp = array([Props('T','H',H,'P',p,Fluid) for H in h])
    ax2 = fig.add_subplot(122)
    ax2.semilogy(Tv,abs(p/p_Trho-1)*100,'o',label='p & p(T,rho)')
    ax2.semilogy(Tv,abs(T_hp/Tv-1)*100,'o',label='T & T(h,p)')
    ax2.set_ylim(1e-16,100)
    ax2.set_title('Superheated vapor')
    ax2.set_xlabel('Temperature [K]')
    ax2.set_ylabel('Absolute deviation [%]')
    ax2.legend(numpoints=1,loc='best')

    fig.tight_layout()
    plt.show()
