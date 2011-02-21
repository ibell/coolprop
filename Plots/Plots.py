import pylab
import numpy as np
import CoolProp.CoolProp as cp
def Ts(Ref,**kwargs):

    if 'Tmin' in kwargs:
        Tmin=float(kwargs['Tmin'])
    else:
        Tmin=220.
        
    if 'axis' in kwargs:
        ax=kwargs['axis']
    else:
        ax=pylab.gca()
        
    if 'bounds' in kwargs:
        if kwargs['bounds']=='R410A':
            smin=0
            smax=3
            Tmin=200
            Tmax=450
            sticks=(1,2,3)
        if kwargs['bounds']=='R744' or kwargs['bounds']=='CO2':
            smin=-2.0
            smax=0.0
            sticks=(-1.8,-1.4,-1.0,-0.6,-0.2)
            Tmin=200
            Tmax=450
        if kwargs['bounds']=='Nitrogen':
            smin=3
            smax=10
            sticks=(4,6,8)
            Tmin=100
            Tmax=600
        if kwargs['bounds']=='Argon':
            smin=-3
            smax=0
            sticks=(-3.0,-2.0,-1.0,0.0)
            Tmin=100
            Tmax=600
    else:
        smin=-2.0
        smax=0.0
        sticks=(-1.8,-1.4,-1.0,-0.6,-0.2)
        Tmin=200
        Tmax=450
        
    if 'sbounds' in kwargs:
        smin=kwargs['sbounds'][0]
        smax=kwargs['sbounds'][1]
    if 'Tbounds' in kwargs:
        Tmin=kwargs['Tbounds'][0]
        Tmax=kwargs['Tbounds'][1]

    Tsat = np.linspace(Tmin,cp.Tcrit(Ref)-0.0001,1000)
    (ssatL,psatL,ssatV,psatV)=(0.0*Tsat,0.0*Tsat,0.0*Tsat,0.0*Tsat)
    for i in range(len(Tsat)):
        ssatL[i] = cp.Props('S','T',Tsat[i],'Q',0,Ref)
        ssatV[i] = cp.Props('S','T',Tsat[i],'Q',1,Ref)
        
    ax.plot(ssatL,Tsat,'k')
    ax.plot(ssatV,Tsat,'k')
    
    ax.set_xticks(sticks)
    ax.set_xlim((smin,smax))
    ax.set_ylim((Tmin,Tmax))
    ax.set_xlabel('Entropy [kJ/kg$\cdot$K]')
    ax.set_ylabel('Temperature [K]')

def Ph(Ref,**kwargs):
    if 'Tmin' in kwargs:
        Tmin=float(kwargs['Tmin'])
    else:
        Tmin=220.
        
    if 'axis' in kwargs:
        ax=kwargs['axis']
    else:
        ax=pylab.gca()
        
    if 'bounds' in kwargs:
        if kwargs['bounds']=='R410A':
            hmin=100
            hmax=600
            pmin=0
            pmax=2500
            hticks=(100,200,300,400,500,600)
        if kwargs['bounds']=='R744' or kwargs['bounds']=='CO2':
            hmin=-500
            hmax=200
            hticks=(-400,-200,0,200)
            pmin=0
            pmax=13000
        if kwargs['bounds']=='Nitrogen':
            hmin=-100
            hmax=600
            pmin=0
            pmax=20000
            hticks=(0,200,400,600)
        if kwargs['bounds']=='Argon':
            hmin=-300
            hmax=200
            pmin=0
            pmax=20000
            hticks=(-300,-200,-100,0,100,200)    
    else:
        hmin=-500
        hmax=200
        pmin=0
        pmax=8000
        hticks=(-400,-200,0,200)
        
    if 'hbounds' in kwargs:
        hmin=kwargs['hbounds'][0]
        hmax=kwargs['hbounds'][1]
    if 'pbounds' in kwargs:
        pmin=kwargs['pbounds'][0]
        pmax=kwargs['pbounds'][1]

    Tsat = np.linspace(Tmin,cp.Tcrit(Ref)-0.0001,1000)
    (hsatL,psatL,hsatV,psatV)=(0.0*Tsat,0.0*Tsat,0.0*Tsat,0.0*Tsat)
    for i in range(len(Tsat)):
        hsatL[i] = cp.Props('H','T',Tsat[i],'Q',0,Ref)
        hsatV[i] = cp.Props('H','T',Tsat[i],'Q',1,Ref)
        psatL[i] = cp.Props('P','T',Tsat[i],'Q',0,Ref)
        psatV[i] = cp.Props('P','T',Tsat[i],'Q',1,Ref)

    ax.plot(hsatL,psatL,'k')
    ax.plot(hsatV,psatV,'k')
    ax.plot(np.r_[hsatL[-1],hsatV[-1]],np.r_[psatL[-1],psatV[-1]],'k')
    
    epsh=10
    epsp=-100
    
    ax.set_xticks(hticks)
    ax.set_xlim((hmin,hmax))
    ax.set_ylim((pmin,pmax))

    ax.set_xlabel('Enthalpy [kJ/kg]')
    ax.set_ylabel('Pressure [kPa]')