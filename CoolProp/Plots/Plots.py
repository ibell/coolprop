import pylab, numpy as np, CoolProp.CoolProp as cp

#Turn off lookup for sure
cp.UseSaturationLUT(0) 

def show():
    """
    A convenience function to call pylab.show()
    """
    pylab.show()
    
def Ts(Ref,Tmin=220,axis=None,show=False, **kwargs):
    """
    Make a temperature- entropy plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """
    
    if axis==None:
        ax=pylab.gca()
    else:
        ax=axis

    Tsat = np.linspace(Tmin,cp.Props(Ref,"Tcrit")-0.001,100)
    (ssatL,psatL,ssatV,psatV)=(0.0*Tsat,0.0*Tsat,0.0*Tsat,0.0*Tsat)
    for i in range(len(Tsat)):
        ssatL[i] = cp.Props('S','T',Tsat[i],'Q',0,Ref)
        ssatV[i] = cp.Props('S','T',Tsat[i],'Q',1,Ref)
        
    ax.plot(ssatL,Tsat,'k')
    ax.plot(ssatV,Tsat,'k')
    ax.plot(np.r_[ssatL[-1],ssatV[-1]],np.r_[Tsat[-1],Tsat[-1]],'k')

    ax.set_xlabel('Entropy [kJ/kg$\cdot$K]')
    ax.set_ylabel('Temperature [K]')
    ax.autoscale(enable=True)
    if show:
        pylab.show()

def Ph(Ref,axis=None,Tmin=220,show = False, **kwargs):
    
    """
    Make a pressure-enthalpy plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """
    if axis==None:
        ax=pylab.gca()
    else:
        ax=axis
    Tsat = np.linspace(Tmin,cp.Props(Ref,"Tcrit")-0.001,1000)
    (hsatL,psatL,hsatV,psatV)=(0.0*Tsat,0.0*Tsat,0.0*Tsat,0.0*Tsat)
    for i in range(len(Tsat)):
        hsatL[i] = cp.Props('H','T',Tsat[i],'Q',0,Ref)
        hsatV[i] = cp.Props('H','T',Tsat[i],'Q',1,Ref)
        psatL[i] = cp.Props('P','T',Tsat[i],'Q',0,Ref)
        psatV[i] = cp.Props('P','T',Tsat[i],'Q',1,Ref)


    ax.plot(hsatL,psatL,'k')
    ax.plot(hsatV,psatV,'k')
    ax.plot(np.r_[hsatL[-1],hsatV[-1]],np.r_[psatL[-1],psatV[-1]],'k')
    
    ax.set_xlabel('Enthalpy [kJ/kg]')
    ax.set_ylabel('Pressure [kPa]')
    ax.autoscale(enable=True)
    if show:
        pylab.show()
    
    
if __name__=='__main__':
    Ph('R290', show = True)
    Ts('R290', show = True)
