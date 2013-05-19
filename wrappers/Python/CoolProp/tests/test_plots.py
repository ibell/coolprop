import numpy as np
import pylab
from CoolProp.Plots.Plots import getIsoLines, drawIsoLines, Ts

#def test_Trho():
#    for Fluid in CoolProp.__fluids__:
#        for T in np.linspace(Props(Fluid,'Tmin')+1e-5,Props(Fluid,'Tcrit')+100,20):
#            for rho in np.linspace(1e-10,Props(Fluid,'rhocrit')*3,20):
#                yield check_rho,Fluid,T,rho
#
#def check_rho(Fluid,T,rho):
#    p = Props('P','T',T,'D',rho,Fluid)
    
def test_namedPlots():
    return True

def test_functions():
    
    Ref = 'n-Pentane'
        
    # Version A, get lines and do the plotting here
    fig, ax = pylab.subplots(1, 1)
    # Set limits, be sure to use internal CoolProp Units!
    ax.set_xlim([-0.5,1.5])
    ax.set_ylim([300,530])
    lines = []
    lines.extend(getIsoLines(Ref, 'Ts', 'Q', [0.0, 0.6,  0.75, 0.775, 0.8, 1.0], axis=ax))
    lines.extend(getIsoLines(Ref, 'Ts', 'P', [100, 2000], num=5, axis=ax))
    lines.extend(getIsoLines(Ref, 'Ts', 'D', [2,    600], num=5, axis=ax))
#    lines.extend(getIsoLines(Ref, 'Ts', 'H', 100,  300, 5, ax))
    for line in lines:
        ax.plot(line['x'],np.array(line['y'])-273.15,**line['opts'])        
    # Adjust the T limits to Celsius
    ax.set_ylim([25,250])
    
    # Version B, use built-in drawing functions and receiver line objects
    fig, ax = pylab.subplots(1, 1)
    ax = Ts(Ref)
    ax.set_xlim([-0.5,1.5])
    ax.set_ylim([300,530])
    quality    = drawIsoLines(Ref, 'Ts', 'Q', [0.3,  0.75, 0.775, 0.8], axis=ax)
    isobars    = drawIsoLines(Ref, 'Ts', 'P', [100, 2000]             , num=5, axis=ax)
    isochores  = drawIsoLines(Ref, 'Ts', 'D', [2,    600]             , num=7, axis=ax)
#    isenthalps = drawIsoLines(Ref, 'Ts', 'H', 100,  300, 5, ax)

    
    #pylab.show()
    
    
if __name__=='__main__':
    import nose
    nose.runmodule()