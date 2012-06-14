import CoolProp.CoolProp as CP
import CoolProp
import textwrap
import os

def CriticalIsotherm(Fluid):
    
    return textwrap.dedent(
    """
    Along the critical isotherm where T=T\ :sub:`c`
    ================================================
    .. plot::
    
        Fluid = "{Fluid:s}"
        RPFluid = "{RPFluid:s}"
        
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
            
        Rp = array([Props('P','T',Tc,'D',D,RPFluid) for D in rhov])
        Rrho = array([Props('D','T',Tc,'D',D,RPFluid) for D in rhov])
        Rcp = array([Props('C','T',Tc,'D',D,RPFluid) for D in rhov])
        Rcv = array([Props('O','T',Tc,'D',D,RPFluid) for D in rhov])
            
        fig = plt.figure()

        ax = fig.add_axes((0.15,0.15,0.8,0.8))
        ax.semilogy(rhov/rhoc,abs(p/Rp-1)*100,'o',label='Pressure')
        ax.semilogy(rhov/rhoc,abs(cp/Rcp-1)*100,'o',label='Specific heat (cp)')
        ax.semilogy(rhov/rhoc,abs(cv/Rcv-1)*100,'o',label='Specific heat (cv)')
        ax.set_ylim(1e-16,100)
        ax.set_title('Critical isotherm Deviations from REFPROP 9.0')
        ax.set_xlabel(r'Reduced density $\\rho/\\rho_c$')
        ax.set_ylabel('Absolute deviation [%]')
        ax.legend(numpoints=1,loc='best')
        plt.show()
    """.format(Fluid=Fluid,RPFluid = 'REFPROP-'+CP.get_REFPROPname(Fluid)))

def SatVaporParity(Fluid):
    
    return textwrap.dedent(
    """
    Saturated Vapor Deviations
    ==========================
    
    .. plot::
    
        Fluid = "{Fluid:s}"
        RPFluid = "{RPFluid:s}"
        
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
    """.format(Fluid=Fluid,RPFluid = 'REFPROP-'+CP.get_REFPROPname(Fluid)))

def SatLiquidParity(Fluid):
    
    return textwrap.dedent(
    """
    Saturated Liquid Deviations
    ===========================
    
    .. plot::
    
        Fluid = "{Fluid:s}"
        RPFluid = "{RPFluid:s}"
        
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
    """.format(Fluid=Fluid,RPFluid = 'REFPROP-'+CP.get_REFPROPname(Fluid)))
    
def params_table(Fluid):
    params = dict(mm = CP.Props(Fluid,'molemass'),
                  Tt = CP.Props(Fluid,'Ttriple'),
                  pc = CP.Props(Fluid,'pcrit'),
                  )
    
    return textwrap.dedent(
    """
    Fluid Data
    ==========
        
    Fluid Parameters
    
    =========================  ==============================
    Mole Mass [kg/kmol]        {mm:0.5f}
    Triple Point [K]           {Tt:0.3f}
    =========================  ==============================
    """.format(**params))
    
def critical_table(Fluid):
    params = dict(Tc = CP.Props(Fluid,'Tcrit'),
                  rhoc = CP.Props(Fluid,'rhocrit'),
                  pc = CP.Props(Fluid,'pcrit'),
                  )
    
    return textwrap.dedent(
    """
    Critical Parameters
    
    ==========================  ==============================
    Temperature [K]             {Tc:0.2f}
    Density [kg/m\ :sup:`3`\ ]   {rhoc:0.6f}
    Pressure [kPa]              {pc:0.5f}
    ==========================  ==============================
    
    """.format(**params))

def buildrst(Fluid):
    pass
    
def index_file(Fluids):
    Fluids=[Fluid+'.rst' for Fluid in Fluids]
    FluidList = '\n    '.join(Fluids)
    return textwrap.dedent( 
"""
.. toctree::
    :maxdepth: 1

    {:s}
""".format(FluidList))
    
def fluid_header(Fluid):
    return textwrap.dedent(
"""
********************
{Fluid:s}
********************

Equation of State Reference
===========================
{Reference:s}

Transport Properties Information
================================
{Transport:s}
    
""".format(Fluid=Fluid,
               Reference=CP.get_EOSReference(Fluid),
               Transport=CP.get_TransportReference(Fluid)))
    
fp=open(os.path.join('Fluids','FluidInformation.rst'),'w')
fp.write('###########################\nFluid Information\n###########################\n')
fp.write(
"""
.. toctree::
    :maxdepth: 1

    PseudoPureFluids.rst
    PureFluids.rst
"""
)
fp.close()

pseudo_pure_fluids = ['Air','R404A','R410A','R407C','R507A']
with open(os.path.join('Fluids','PseudoPureFluids.rst'),'w') as fp:
    fp.write('#######################\nPseudo-Pure Fluids\n#######################\n')
    fp.write(index_file(pseudo_pure_fluids))

for Fluid in pseudo_pure_fluids:
    fp=open(os.path.join('Fluids',Fluid+'.rst'),'w')
    s = fluid_header(Fluid) + params_table(Fluid) + critical_table(Fluid) + SatVaporParity(Fluid) + SatLiquidParity(Fluid) + CriticalIsotherm(Fluid)
    fp.write(s)
    fp.close()

pure_fluids = sorted([fluid for fluid in CoolProp.__fluids__ if fluid not in pseudo_pure_fluids])
with open(os.path.join('Fluids','PureFluids.rst'),'w') as fp:
    fp.write('#######################\nPure Fluids\n#######################\n')
    fp.write(index_file(sorted(pure_fluids)))
    
for Fluid in sorted(pure_fluids):
    fp=open(os.path.join('Fluids',Fluid+'.rst'),'w')
    s = fluid_header(Fluid) + params_table(Fluid) + critical_table(Fluid) + SatVaporParity(Fluid) + SatLiquidParity(Fluid) + CriticalIsotherm(Fluid)
    fp.write(s)
    fp.close()

    