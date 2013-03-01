import CoolProp
from CoolProp.CoolProp import Props,get_REFPROPname
import numpy as np

print len(CoolProp.__fluids__)

def test_visc():
    for Fluid in CoolProp.__fluids__:
        for Tsat in np.linspace(Props(Fluid,"Tcrit")-0.1,Props(Fluid,"Ttriple")+0.1,10):
            try:
                mu1 = Props("P",'T',Tsat,'Q',0,Fluid)
                mu2 = Props("P",'T',Tsat,'Q',0,'REFPROP-'+get_REFPROPname(Fluid))
                err = (mu1 - mu2)/mu2*100
                if abs(err)>1e-1:
                    print Fluid,Tsat,mu1,mu2,err,'%'
            except:
                Fluid, Tsat
                
if __name__=='__main__':
    test_visc()