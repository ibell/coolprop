from CoolProp import CoolProp as CP
from PDSim.misc.datatypes import Collector
import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import *
import textwrap

#
#fluid = 'Propane'
#Rfluid = 'REFPROP-propane'
#e_k = 263.88
#sigma = 0.49748
#
fluid = 'R245fa'
Rfluid = 'REFPROP-R245fa'
e_k = 329.72
sigma = 0.5529
molemass = CP.Props(fluid,'molemass')

Ttriple = CP.Props(fluid,'Ttriple')
Tcrit = CP.Props(fluid,'Tcrit')
rhocrit = CP.Props(fluid,'rhocrit')
n = 6
m = 3
NP = 1
Nb = 0
N = (n-1)*(m+1)+3+Nb

mu,mu_dilute,RHO,TTT = Collector(),Collector(),Collector(),Collector()

rhomax = CP.Props('D','T',Ttriple,'Q',0,'R245fa')
#Build a database of "experimental" data
for T in np.linspace(Ttriple,Tcrit+30,400):
    for rho in np.linspace(1e-10,rhomax,400):
        muval = CP.Props('V','T',T,'D',rho,Rfluid)
        mudilute = CP.viscosity_dilute(fluid,T,rho,e_k,sigma)
        
        #Want positive value, and single-phase
        if (muval > 0 and T > Tcrit or rho > CP.rhosatL_anc(fluid,T) or rho < CP.rhosatV_anc(fluid,T)):
            mu << muval
            mu_dilute << mudilute
            TTT << T
            RHO << rho

from CoolProp.Plots.Plots import Trho
Trho(fluid)
plt.plot(RHO.vec,TTT.vec,'.')
plt.show()

#tau = np.array(TTT.vec)/Tcrit
tau = np.array(TTT.vec)/Tcrit
delta = np.array(RHO.vec)/rhocrit
Tstar = np.array(TTT.vec)/e_k

#Define the objective function
def OBJECTIVE_fit(c,x):
    tau = x[0,:]
    delta = x[1,:]
    #Unpack the inputs into e matrix and f vector
    e = np.zeros((n+1,m+1))
    
    sum = 0
    k = 0
    for i in range(2,n+1):
        for j in range(0,m+1):
            e[i][j] = c[k]
            sum += e[i][j]*delta**i/tau**j
            k += 1
            
    for o in range(0,NP):
        f1 = c[k+o*3]
        g1 = c[k+1+o*3]
        g2 = c[k+2+o*3]
    
    delta_0 = g1*(1+g2*tau**0.5)
    sum += f1*(delta/(delta_0-delta)-delta/delta_0)
    return sum + np.array(mu_dilute.vec)
    
print 'starting fit'
XXX = np.r_[np.array(tau,ndmin = 2), np.array(delta,ndmin=2)]
mod = Model(OBJECTIVE_fit)
mydata = Data(XXX.copy(), np.array(mu.vec))
beta0  = [1 for _ in range(N)]
myodr = ODR(mydata, mod, beta0=beta0)
myoutput = myodr.run()
myoutput.pprint()
print myoutput.sum_square
YFIT = OBJECTIVE_fit(myoutput.beta,XXX)
plt.plot(np.array(mu.vec),YFIT/np.array(mu.vec),'o')
plt.show()

rel_error = (YFIT)/np.array(mu.vec)-1
MAE = np.mean(np.abs(rel_error))*100
SSE = np.sum(np.power(YFIT-np.array(mu.vec),2))
print SSE

def write_output(c):
    e = np.zeros((n+1,m+1))
    k = 0
    edata = ''
    for i in range(2,n+1):
        erow = ''
        for j in range(0,m+1):
            e[i][j] = c[k]
            erow += 'e[{i:d}][{j:d}] = {val:0.16g}; '.format(val = e[i][j],i=i,j=j)
            k += 1
        edata += erow + '\n'
    f1 = c[k]
    g1 = c[k+1]
    g2 = c[k+2]
        
    template = textwrap.dedent(
    """
    double {name:s}Class::viscosity_Trho(double T, double rho)
    {{
        // This function was generated by fitting REFPROP ECS data
        // to the functional form of Vogel, 1998 (propane viscosity)
        // The script entitled dev/fit_avoid_ECS.py was used to make this 
        // function.  The mean absolute error of the fit is equal to
        // {MAE:g} % 
        
        double delta_0, sum, DELTA_H_eta, e_k, sigma, tau, delta;
        double e[{n:d}+1][{m:d}+1];
        
        tau = T/reduce.T; //[Opposite to normal definition]
        delta = rho/reduce.rho;
        
        // Load the coefficients
        double f1 = {f1:0.16g}, g1 = {g1:0.16g}, g2 = {g2:0.16g};
        for (int i=0;i<={n:d};i++){{ for(int j=0;j<={m:d};j++){{ e[i][j]=0.0; }} }}
        {edata:s}
        delta_0=g1*(1+g2*sqrt(tau)); //[no units]
        sum=0;
        for (int i=2;i<={n:d};i++){{
            for (int j=0;j<={m:d};j++){{
                sum += e[i][j]*pow(delta,i)/pow(tau,j);
            }}
        }}
        DELTA_H_eta = sum + f1*(delta/(delta_0-delta)-delta/delta_0); //[Pa-s]
        
        try{{
            // Get the ECS params for the fluid if it has them
            ECSParams(&e_k,&sigma);
        }}
        catch(NotImplementedError)
        {{
            throw ValueError(format("Your fluid does not implement ECSParams"));
        }}
        
        return viscosity_dilute(T,e_k,sigma) + DELTA_H_eta;
    }} 
    """
    )
    
    values = dict(f1 = f1,
                  g1 = g1,
                  g2 = g2,
                  n = n,
                  m = m,
                  name = fluid,
                  edata = edata,
                  MAE = MAE)
    
    print template.format(**values)

write_output(myoutput.beta)