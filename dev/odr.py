from __future__ import division
import numpy as np
from scipy.odr import *
import scipy.optimize, random
import matplotlib.pyplot as plt
import textwrap

def rsquared(x, y):
    """ 
    Return R^2 where x and y are array-like.
    
    from http://stackoverflow.com/questions/893657/how-do-i-calculate-r-squared-using-python-and-numpy
    """
    import scipy.stats
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    return r_value**2

def saturation_density(Ref, ClassName, form = 'A', LV = 'L', perc_error_allowed = 0.3):
    """
    
    Parameters
    ----------
    Ref : string
        The fluid name for the fluid that will be used to generate the saturation data
    ClassName : The name of the class that will be used in the C++ code
    form : string
        If ``'A'``, use a term of the form 
    """
    
    from CoolProp.CoolProp import Props
    
    Tc = Props(Ref,'Tcrit')
    pc = Props(Ref,'pcrit')
    rhoc = Props(Ref,'rhocrit')
    Tmin = Props(Ref,'Tmin')
    
    TT = np.linspace(Tmin+1e-6, Tc-0.00001, 1000)
    p = np.array([Props('P','T',T,'Q',0,Ref) for T in TT])
    rhoL = np.array([Props('D','T',T,'Q',0,Ref) for T in TT])
    rhoV = np.array([Props('D','T',T,'Q',1,Ref) for T in TT])
    
    # Start with a large library of potential powers
    n = [i/6.0 for i in range(1,200)]#+[0.35+i/200.0 for i in range(1,70)]+[0.05+0.01*i for i in range(1,70)]
    
    x = 1.0-TT/Tc
    if form == 'A' and LV == 'L':
        y = np.array(rhoL)/rhoc-1
    elif form == 'A' and LV == 'V':
        y = np.array(rhoV)/rhoc-1
    elif form == 'B' and LV == 'L':
        y = (np.log(rhoL)-np.log(rhoc))*TT/Tc
    elif form == 'B' and LV == 'V':
        y = (np.log(rhoV)-np.log(rhoc))*TT/Tc
    else:
        raise ValueError
    
    plt.plot(x,y)
    plt.show()
        
    max_abserror = 0
    while len(n) > 3:
        print max_abserror, len(n)
    
        def f_p(B, x):
            # B is a vector of the parameters.
            # x is an array of the current x values.
            return sum([_B*x**(_n) for _B,_n in zip(B,n)])
        
        linear = Model(f_p)
        mydata = Data(x, y)
        myodr = ODR(mydata, linear, beta0=[0]*len(n))
        myoutput = myodr.run()
        
        beta = myoutput.beta
        sd = myoutput.sd_beta
        
        if form == 'A':
            rho_fit = (f_p(myoutput.beta,x)+1)*rhoc
        elif form == 'B':
            rho_fit = np.exp(f_p(myoutput.beta,x)*Tc/TT)*rhoc
        else:
            raise ValueError
        
        if LV == 'L':
            max_abserror = np.max(np.abs(rho_fit/rhoL-1))*100
        else:
            max_abserror = np.max(np.abs(rho_fit/rhoV-1))*100
            
        dropped_indices = [i for i in range(len(n)//2) if abs(sd[i])<1e-15 ]
        if dropped_indices:
            for i in reversed(sorted(dropped_indices)):
                n.pop(i)
            print 'popping...', len(n), 'terms remaining'
            continue
        
        if max_abserror > perc_error_allowed:
            break # The last good run will be used
        else:
            print max_abserror
            Ncoeffs = str(list(myoutput.beta)).lstrip('[').rstrip(']')
            tcoeffs = str(n).lstrip('[').rstrip(']')
            maxerror = max_abserror
            if form == 'A':
                code_template = textwrap.dedent(
                """
                for (int i=1; i<={count:d}; i++)
                {{
                    summer += N[i]*pow(theta,t[i]);
                }}
                return reduce.rho*(summer+1);
                """.format(count = len(n))
                )
            elif form == 'B':
                code_template = textwrap.dedent(
                """
                for (int i=1; i<={count:d}; i++)
                {{
                    summer += N[i]*pow(theta,t[i]);
                }}
                return reduce.rho*exp(reduce.T/T*summer);
                """.format(count = len(n))
                )
            else:
                raise ValueError
            
        # Find the least significant entry (the one with the largest relative standard error)
        # and remove it
        n.pop(np.argmax(np.abs(sd/beta)))
        
        #Remove elements that are not 
    template = textwrap.dedent(
    """
    double {name:s}Class::rhosat{LV:s}(double T)
    {{
        // Maximum absolute error is {error:f} % between {Tmin:f} K and {Tmax:f} K
        const double t[] = {{0, {tcoeffs:s}}};
        const double N[] = {{0, {Ncoeffs:s}}};
        double summer=0,theta;
        theta=1-T/reduce.T;
        \t{code:s}
    }}
    """)
    f = open('anc.txt','a')
    f.write(template.format(tcoeffs = tcoeffs,
                            Ncoeffs = Ncoeffs,
                            name = ClassName,
                            Tmin = TT[0],
                            Tmax = TT[-1],
                            error = maxerror,
                            code = code_template,
                            LV = LV
                            ))
    f.close()
    return

def saturation_pressure(Ref, ClassName):
    
    from CoolProp.CoolProp import Props
    
    Tc = Props(Ref,'Tcrit')
    pc = Props(Ref,'pcrit')
    rhoc = Props(Ref,'rhocrit')
    Tmin = Props(Ref,'Tmin')
    
    TT = np.linspace(Tmin+1e-6, Tc-0.00001, 300)
    p = [Props('P','T',T,'Q',0,Ref) for T in TT]
    rhoL = [Props('D','T',T,'Q',0,Ref) for T in TT]
    rhoV = [Props('D','T',T,'Q',1,Ref) for T in TT]
        
    Np = 60
    n = range(1,Np)
    max_abserror = 0
    while len(n) > 3:
    
        def f_p(B, x):
            # B is a vector of the parameters.
            # x is an array of the current x values.
            return sum([_B*x**(_n/2.0) for _B,_n in zip(B,n)])
        x = 1.0-TT/Tc
        y = (np.log(p)-np.log(pc))*TT/Tc
        
        linear = Model(f_p)
        mydata = Data(x, y)
        myodr = ODR(mydata, linear, beta0=[0]*len(n))
        myoutput = myodr.run()
        
        beta = myoutput.beta
        sd = myoutput.sd_beta
        
        p_fit = np.exp(f_p(myoutput.beta,x)*Tc/TT)*pc
        max_abserror = np.max(np.abs((p_fit/p)-1)*100)
        print max_abserror
        psat_error = max_abserror
        
        dropped_indices = [i for i in range(len(n)) if abs(sd[i])<1e-15 ]
        if dropped_indices:
            #for i in reversed(dropped_indices):
            # randomly drop one of them
            n.pop(random.choice(dropped_indices))
            continue
        
        if max_abserror < 0.1: #Max error is 0.1%
            Ncoeffs = str(list(myoutput.beta)).lstrip('[').rstrip(']')
            tcoeffs = str(n).lstrip('[').rstrip(']')
            maxerror = max_abserror
            N = len(n)
        else:
            break
        
        # Find the least significant entry (the one with the largest standard error)
        # and remove it
        n.pop(np.argmax(sd))
        
        #Remove elements that are not 
    import textwrap
    template = textwrap.dedent(
    """
    double {name:s}Class::psat(double T)
    {{
        // Maximum absolute error is {psat_error:f} % between {Tmin:f} K and {Tmax:f} K
        const double t[]={{0, {tcoeffs:s}}};
        const double N[]={{0, {Ncoeffs:s}}};
        double summer=0,theta;
        theta=1-T/reduce.T;
        for (int i=1;i<={N:d};i++)
        {{
            summer += N[i]*pow(theta,t[i]/2);
        }}
        return reduce.p*exp(reduce.T/T*summer);
    }}
    """)
    f = open('anc.txt','a')
    f.write(template.format(N = len(n),
                            tcoeffs = tcoeffs,
                            Ncoeffs = Ncoeffs,
                            name = ClassName,
                            Tmin = TT[0],
                            Tmax = TT[-1],
                            psat_error = maxerror,
                            ))
    f.close()
    return
                      
for RPFluid,Fluid in [('REFPROP-Methanol','Methanol')
                    ]:
#    saturation_pressure(RPFluid, Fluid)
#    saturation_density(RPFluid, Fluid, form='A', LV='L')
    saturation_density(RPFluid, Fluid, form='B', LV='V', perc_error_allowed = 0.4)
