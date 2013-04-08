import numpy as np
from scipy.odr import *
import scipy.optimize, random

def example1():
    #The simplest example, fit slope and intercept to data that is a line
    def f(B, x):
        ''' Linear function y = m*x + b '''
        return B[0]*x + B[1]

        # B is a vector of the parameters.
        # x is an array of the current x values.
        # x is same format as the x passed to Data or RealData.

        # Return an array in the same format as y passed to Data or RealData.
        
    linear = Model(f)

    x = np.linspace(0.01,1000,300)
    y = 1.5*x+3
    mydata = Data(x, y)

    myodr = ODR(mydata, linear, beta0=[1., 2.])

    myoutput = myodr.run()

    myoutput.pprint()
    
    
def example2(Ref, ClassName,  N = 4, addTr = True, Llog = True):
    
    from CoolProp.CoolProp import Props
    
    Tc = Props(Ref,'Tcrit')
    pc = Props(Ref,'pcrit')
    rhoc = Props(Ref,'rhocrit')
    Tmin = Props(Ref,'Tmin')
    
    TT = np.linspace(Tmin+1e-6, Tc-0.00001,300)
    p = [Props('P','T',T,'Q',0,Ref) for T in TT]
    rhoL = [Props('D','T',T,'Q',0,Ref) for T in TT]
    rhoV = [Props('D','T',T,'Q',1,Ref) for T in TT]
    logppc = (np.log(p)-np.log(pc))*TT/Tc
    logrhoLrhoc = np.log(rhoL)-np.log(rhoc)
    logrhoVrhoc = np.log(rhoV)-np.log(rhoc)
        
    def f_p(B, x):
        # B is a vector of the parameters.
        # x is an array of the current x values.
        return B[0]*x + B[1]*x**1.5 + B[2]*x**2.3 + B[3]*x**3.6 + B[4]*x**5.2 + B[5]*x**7.3 + B[6]*x**9
    x = 1.0-TT/Tc
    y = logppc
    
    linear = Model(f_p)
    mydata = Data(x, y)
    myodr = ODR(mydata, linear, beta0=[0, 0, 0, 0, 0, 0, 0])
    myoutput = myodr.run()
    
    psat_N = str(list(myoutput.beta)).lstrip('[').rstrip(']')
    p_fit = np.exp(f_p(myoutput.beta,x)*Tc/TT)*pc
    max_abserror = np.max(np.abs((p_fit/p)-1)*100)
    print max_abserror
    psat_error = max_abserror
    
    x = 1.0-TT/Tc
            
    t = 0
    max_abserror = 0

    if addTr:
        Tr = TT/Tc
    else:
        Tr = 1
        
    y = np.log(np.array(rhoV)/rhoc)*Tr
    
    def f_coeffs(n, x):
        global t, max_abserror

        def f_rho(B, x):
            return sum([_B*x**_n for _B,_n in zip(B,n)])
        linear = Model(f_rho)
        mydata = Data(x, y)
        myodr = ODR(mydata, linear, beta0=[1]*N)
        myoutput = myodr.run()
    
        rho_fit = np.exp(f_rho(myoutput.beta,x)/Tr)*rhoc
        max_abserror = np.max(np.abs((rho_fit/rhoV)-1)*100)
        
        t = myoutput.beta
        print Ref, 'rhoV SSE =', myoutput.sum_square, max_abserror,'%', list(myoutput.beta), list(n)
        return myoutput.sum_square
        
    n0 = [0.345]+range(1,N)
    assert(len(n0) == N)
    n0 = scipy.optimize.minimize(f_coeffs,n0, args=(x,), method = 'Powell').x
    max_abserror = globals()['max_abserror']
    t = globals()['t']
    rhosatV_t = str(list(t)).lstrip('[').rstrip(']')
    rhosatV_N = str(list(n0)).lstrip('[').rstrip(']')
    rhosatV_error = max_abserror
    
    if Llog:
        y = np.log(np.array(rhoL)/rhoc)
    else:
        y = np.array(rhoL)/rhoc-1
        
    beta0 = [1]*N
    def f_coeffs(n,x):
        global t, max_abserror, beta0
        def f_rho(B, x):
            return sum([_B*x**_n for _B,_n in zip(B,n)])
            
        linear = Model(f_rho)
        mydata = Data(x, y)
        myodr = ODR(mydata, linear, beta0=[1]*N)
        myoutput = myodr.run()
        
        if Llog:
            rho_fit = np.exp(f_rho(myoutput.beta,x))*rhoc
        else:
            rho_fit = (f_rho(myoutput.beta,x)+1)*rhoc
            
        max_abserror = np.max(np.abs((rho_fit/rhoL)-1)*100)
        
        t = myoutput.beta
        print Ref, 'rhoL SSE =', myoutput.sum_square, max_abserror,'%', list(myoutput.beta), list(n), max_abserror,'%'
        return myoutput.sum_square
    
    n0 = [0.345]+list(0.75*np.array(range(1,N)))
    assert(len(n0) == N)
    n0 = scipy.optimize.minimize(f_coeffs,n0, args=(x,), method = 'Powell').x
    max_abserror = globals()['max_abserror']
    t = globals()['t']
    rhosatL_t = str(list(t)).lstrip('[').rstrip(']')
    rhosatL_N = str(list(n0)).lstrip('[').rstrip(']')
    rhosatL_error = max_abserror
    
    import textwrap
    template = textwrap.dedent(
    """
    double {name:s}Class::psat(double T)
    {{
        // Maximum absolute error is {psat_error:f} % between {Tmin:f} K and {Tmax:f} K
        const double ti[]={{0,1.0,1.5,2.3,3.6,5.2,7.3,9}};
        const double Ni[]={{0,{psat_N:s} }};
        double summer=0,theta;
        int i;
        theta=1-T/reduce.T;
        for (i=1;i<=6;i++)
        {{
            summer=summer+Ni[i]*pow(theta,ti[i]);
        }}
        return reduce.p*exp(reduce.T/T*summer);
    }}
    double {name:s}Class::rhosatL(double T)
    {{
        // Maximum absolute error is {rhosatL_error:f} % between {Tmin:f} K and {Tmax:f} K
        const double ti[]={{0,{rhosatL_N:s}}};
        const double Ni[]={{0,{rhosatL_t:s}}};
        double summer=0;
        int i;
        double theta;
        theta=1-T/reduce.T;
        for (i=1;i<={N:d};i++)
        {{
            summer+=Ni[i]*pow(theta,ti[i]);
        }}
        return reduce.rho*exp(summer);
    }}
    double {name:s}Class::rhosatV(double T)
    {{
        // Maximum absolute error is {rhosatV_error:f} % between {Tmin:f} K and {Tmax:f} K
        const double ti[]={{0,{rhosatV_N:s}}};
        const double Ni[]={{0,{rhosatV_t:s}}};
        double summer=0,theta;
        int i;
        theta=1.0-T/reduce.T;
        for (i=1;i<={N:d};i++)
        {{
            summer=summer+Ni[i]*pow(theta,ti[i]);
        }}
        return reduce.rho*exp(crit.T/T*summer);
    }}
    """)
    
    f = open('anc.txt','a')
    f.write(template.format(N = N,
                            psat_N = psat_N,
                            rhosatL_t = rhosatL_t,
                            rhosatL_N = rhosatL_N,
                            rhosatV_t = rhosatV_t,
                            rhosatV_N = rhosatV_N,
                            name = ClassName,
                            Tmin = TT[0],
                            Tmax = TT[-1],
                            psat_error = psat_error,
                            rhosatL_error = rhosatL_error,
                            rhosatV_error = rhosatV_error
                            ))
    f.close()
                            
    
#example1()

## example2('REFPROP-pentane','nPentane',N = 4)
## example2('REFPROP-HEXANE','nHexane',N = 5)
## example2('REFPROP-T2Butene','Trans2Butene',N = 5)
## example2('REFPROP-1Butene','1Butene',N = 5)
## example2('REFPROP-C2BUTENE','Cis2Butene',N = 5)
## example2('REFPROP-IBUTENE','IsoButene',N = 5)
## example2('REFPROP-HEPTANE','nHeptane',N = 5)
#example2('REFPROP-CYCLOHEX','Cyclohexane', N = 5)
#example2('REFPROP-SES36','SES36', N = 5)
#example2('REFPROP-R236EA','R236EA', N = 7)
#example2('REFPROP-R227EA','R227EA', N = 8)
#example2('REFPROP-R365MFC','R365MFC', N = 8)
example2('REFPROP-R123','R123', N = 6)
#example2('REFPROP-PROPYLEN','Propylene', N = 5)

#example2('REFPROP-DMC','DMC',N = 5)
#example2('REFPROP-MXYLENE','mXylene',N = 5)
#example2('REFPROP-PXYLENE','pXylene',N = 5)
#example2('REFPROP-EBENZENE','EthylBenzene',N = 5)
