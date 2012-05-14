from CoolProp.CoolProp import Props,Ttriple,Tcrit,Help,UseSaturationLUT
import numpy as np,pylab

def SaturationValidationPlot(Ref,REFPROPRef):
    
    T=np.linspace(Props(Ref,'Ttriple')+0.1,Props(Ref,'Tcrit')-0.1,700)
    hL_REFPROP=np.zeros_like(T)
    hV_REFPROP=np.zeros_like(T)
    uL_REFPROP=np.zeros_like(T)
    uV_REFPROP=np.zeros_like(T)
    sL_REFPROP=np.zeros_like(T)
    sV_REFPROP=np.zeros_like(T)
    rhoL_REFPROP=np.zeros_like(T)
    rhoV_REFPROP=np.zeros_like(T)
    viscL_REFPROP=np.zeros_like(T)
    viscV_REFPROP=np.zeros_like(T)
    kL_REFPROP=np.zeros_like(T)
    kV_REFPROP=np.zeros_like(T)
    pL_REFPROP=np.zeros_like(T)
    pV_REFPROP=np.zeros_like(T)
    
    hL=np.zeros_like(T)
    hV=np.zeros_like(T)
    uL=np.zeros_like(T)
    uV=np.zeros_like(T)
    sL=np.zeros_like(T)
    sV=np.zeros_like(T)
    rhoL=np.zeros_like(T)
    rhoV=np.zeros_like(T)
    viscL=np.zeros_like(T)
    viscV=np.zeros_like(T)
    kL=np.zeros_like(T)
    kV=np.zeros_like(T)
    pL=np.zeros_like(T)
    pV=np.zeros_like(T)
    
    rho0=Props('D','T',0.9*Tcrit(Ref),'Q',0,REFPROPRef)
    h0=Props('H','T',0.9*Tcrit(Ref),'D',rho0,Ref)
    s0=Props('S','T',0.9*Tcrit(Ref),'D',rho0,Ref)
    u0=Props('U','T',0.9*Tcrit(Ref),'D',rho0,Ref)
    h0_REFPROP=Props('H','T',0.9*Tcrit(Ref),'Q',0,REFPROPRef)
    s0_REFPROP=Props('S','T',0.9*Tcrit(Ref),'Q',0,REFPROPRef)
    u0_REFPROP=Props('U','T',0.9*Tcrit(Ref),'Q',0,REFPROPRef)
    for i in range(len(T)-1,-1,-1):
        hL_REFPROP[i]=Props('H','T',T[i],'Q',0,REFPROPRef)-h0_REFPROP
        hV_REFPROP[i]=Props('H','T',T[i],'Q',1,REFPROPRef)-h0_REFPROP
        uL_REFPROP[i]=Props('U','T',T[i],'Q',0,REFPROPRef)-u0_REFPROP
        uV_REFPROP[i]=Props('U','T',T[i],'Q',1,REFPROPRef)-u0_REFPROP
        sL_REFPROP[i]=Props('S','T',T[i],'Q',0,REFPROPRef)-s0_REFPROP
        sV_REFPROP[i]=Props('S','T',T[i],'Q',1,REFPROPRef)-s0_REFPROP
        rhoL_REFPROP[i]=Props('D','T',T[i],'Q',0,REFPROPRef)
        rhoV_REFPROP[i]=Props('D','T',T[i],'Q',1,REFPROPRef)
        viscL_REFPROP[i]=Props('V','T',T[i],'Q',0,REFPROPRef)
        viscV_REFPROP[i]=Props('V','T',T[i],'Q',1,REFPROPRef)
        kL_REFPROP[i]=Props('L','T',T[i],'Q',0,REFPROPRef)
        kV_REFPROP[i]=Props('L','T',T[i],'Q',1,REFPROPRef)
        pL_REFPROP[i]=Props('P','T',T[i],'D',rhoL_REFPROP[i],REFPROPRef)
        pV_REFPROP[i]=Props('P','T',T[i],'D',rhoV_REFPROP[i],REFPROPRef)
        rhoL[i]=Props('D','T',T[i],'Q',0.0,Ref)
        rhoV[i]=Props('D','T',T[i],'Q',1.0,Ref)
        hL[i]=Props('H','T',T[i],'D',rhoL[i],Ref)-h0
        hV[i]=Props('H','T',T[i],'D',rhoV[i],Ref)-h0
        uL[i]=Props('U','T',T[i],'D',rhoL[i],Ref)-u0
        uV[i]=Props('U','T',T[i],'D',rhoV[i],Ref)-u0
        sL[i]=Props('S','T',T[i],'D',rhoL[i],Ref)-s0
        sV[i]=Props('S','T',T[i],'D',rhoV[i],Ref)-s0
        viscL[i]=Props('V','T',T[i],'D',rhoL[i],Ref)
        viscV[i]=Props('V','T',T[i],'D',rhoV[i],Ref)
        kL[i]=Props('L','T',T[i],'D',rhoL[i],Ref)
        kV[i]=Props('L','T',T[i],'D',rhoV[i],Ref)
        pL[i]=Props('P','T',T[i],'D',rhoL[i],Ref)
        pV[i]=Props('P','T',T[i],'D',rhoV[i],Ref)
        
    
    nR=7
    nC=2
    pylab.figure(figsize=(8,12))
    pylab.suptitle(Ref)
    
    ax=pylab.subplot(nR,nC,1)
    ax.set_title('Sat. Liquid')
    ax.semilogy(T,np.abs(rhoL_REFPROP/rhoL-1)*100)
    ax.set_ylabel(r'Error: $\rho$ [kg/m$^3$]')
    
    ax=pylab.subplot(nR,nC,2)
    ax.set_title('Sat. Vapor')
    ax.semilogy(T,np.abs(rhoV_REFPROP/rhoV-1)*100)
    
    ax=pylab.subplot(nR,nC,3)
    ax.semilogy(T,np.abs(pL_REFPROP/pL-1)*100)
    ax.set_ylabel('Error: p [$\%$]')
    
    ax=pylab.subplot(nR,nC,4)
    ax.semilogy(T,np.abs(pV_REFPROP/pV-1)*100)
    
    ax=pylab.subplot(nR,nC,5)
    ax.semilogy(T,np.abs(hL_REFPROP/hL-1)*100)
    ax.set_ylabel('Error: h [$\%$]')
    
    ax=pylab.subplot(nR,nC,6)
    ax.semilogy(T,np.abs(hV_REFPROP/hV-1)*100)
    
    ax=pylab.subplot(nR,nC,7)
    ax.semilogy(T,np.abs(uL_REFPROP/uL-1)*100)
    ax.set_ylabel('Error: u [$\%$]')
    
    ax=pylab.subplot(nR,nC,8)
    ax.semilogy(T,np.abs(uV_REFPROP/uV-1)*100)
    
    ax=pylab.subplot(nR,nC,9)
    ax.semilogy(T,np.abs(sL_REFPROP/sL-1)*100)
    ax.set_ylabel('Error: s [$\%$]')
    
    ax=pylab.subplot(nR,nC,10)
    ax.semilogy(T,np.abs(sV_REFPROP/sV-1)*100)
    
    ax=pylab.subplot(nR,nC,11)
    ax.semilogy(T,np.abs(viscL_REFPROP/viscL-1)*100)
    ax.set_ylabel('Error: $\mu$ [$\%$]')
    
    ax=pylab.subplot(nR,nC,12)
    ax.semilogy(T,np.abs(viscV_REFPROP/viscV-1)*100)
    
    ax=pylab.subplot(nR,nC,13)
    ax.semilogy(T,np.abs(kL_REFPROP/kL-1)*100)
    ax.set_ylabel('Error: k [$\%$]')
    
    ax=pylab.subplot(nR,nC,14)
    ax.semilogy(T,np.abs(kV_REFPROP/kV-1)*100)
    
    pylab.show()
    
if __name__=='__main__':
##     SaturationValidationPlot('R290','REFPROP-Propane')
##     SaturationValidationPlot('Water','REFPROP-water')
##     SaturationValidationPlot('Air','REFPROP-Air')
    SaturationValidationPlot('R717','REFPROP-ammonia')
##     SaturationValidationPlot('R744','REFPROP-CO2')

##   These ones are fully validated
    UseSaturationLUT(0)
    SaturationValidationPlot('R134a','REFPROP-R134a')
##     SaturationValidationPlot('Nitrogen','REFPROP-Nitrogen')
##     SaturationValidationPlot('Argon','REFPROP-Argon')
    pass