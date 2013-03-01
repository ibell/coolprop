from CoolProp.CoolProp import Props
import numpy as np

def test_TSP():
    for Tsat in np.linspace(Props("R134a","Tcrit")-0.1,Props("R134a","Ttriple")+0.1):
        p = Props("P",'T',Tsat,'Q',0,"R134a")
        for DT in np.linspace(0.0001,50,5):
            T = Tsat+DT
            try:
                rho = Props("D",'T',T,'P',p,"R134a")
                s = Props("S",'T',T,'D',rho,"R134a")
                Tchk = Props("T",'S',s,'P',p,"R134a")
                err = Tchk - T
                if abs(err)>1e-8:
                    print Tsat,DT,T,s,p,Tchk,err
            except ValueError as e:
                print DT,Tsat
                print e
                raise
                
if __name__=='__main__':
    test_TSP()