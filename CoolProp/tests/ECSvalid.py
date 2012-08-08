from CoolProp.CoolProp import Props
import numpy as np

Ref = 'R32'
RRef = 'REFPROP-'+Ref
for Tsat in np.linspace(Props(Ref,'Ttriple')+0.1,Props(Ref,'Tcrit')-0.1,10):
    mu1 = Props('L','T',Tsat,'Q',1.0,Ref)
    mu2 = Props('L','T',Tsat,'Q',1.0,RRef)
    err = (mu1-mu2)/mu2*100
    print err,"%"
    