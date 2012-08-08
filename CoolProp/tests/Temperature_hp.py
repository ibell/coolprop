from CoolProp.CoolProp import Props
import numpy as np

def test_THP():
    for Tsat in np.linspace(Props("R134a","Tcrit")-0.1,Props("R134a","Ttriple")+0.1):
        p = Props("P",'T',Tsat,'Q',0,"R134a")
        for DT in np.linspace(0.0001,50,5):
            T = Tsat+DT
            h = Props("H",'T',T,'P',p,"R134a")
            Tchk = Props("T",'H',h,'P',p,"R134a")
            err = Tchk - T
            if abs(err)>1e-8:
                print Tsat,DT,T,h,Tchk,err
                
if __name__=='__main__':
    test_THP()