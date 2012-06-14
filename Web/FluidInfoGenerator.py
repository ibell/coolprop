import CoolProp.CoolProp as CP
import CoolProp
import textwrap
import os

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
    s = fluid_header(Fluid) + params_table(Fluid) + critical_table(Fluid)
    fp.write(s)
    fp.close()

pure_fluids = sorted([fluid for fluid in CoolProp.__fluids__ if fluid not in pseudo_pure_fluids])
with open(os.path.join('Fluids','PureFluids.rst'),'w') as fp:
    fp.write('#######################\nPure Fluids\n#######################\n')
    fp.write(index_file(sorted(pure_fluids)))
    
for Fluid in sorted(pure_fluids):
    fp=open(os.path.join('Fluids',Fluid+'.rst'),'w')
    s = fluid_header(Fluid) + params_table(Fluid) + critical_table(Fluid)
    fp.write(s)
    fp.close()

    