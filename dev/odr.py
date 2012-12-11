import numpy as np
from scipy.odr import *

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
    
    
def example2(Ref, N = 3):
    
    from CoolProp.CoolProp import Props
    
    Tc = Props(Ref,'Tcrit')
    pc = Props(Ref,'pcrit')
    rhoc = Props(Ref,'rhocrit')
    Tmin = Props(Ref,'Tmin')
    
    TT = np.linspace(Tmin+1e-6, Tc-0.000001,10000)
    p = [Props('P','T',T,'Q',0,Ref) for T in TT]
    rhoL = [Props('D','T',T,'Q',0,Ref) for T in TT]
    rhoV = [Props('D','T',T,'Q',1,Ref) for T in TT]
    logppc = (np.log(p)-np.log(pc))*TT/Tc
    logrhoLrhoc = np.log(rhoL)-np.log(rhoc)
    logrhoVrhoc = np.log(rhoV)-np.log(rhoc)
        
    def f_p(B, x):
        # B is a vector of the parameters.
        # x is an array of the current x values.
        return B[0]*x + B[1]*x**1.5 + B[2]*x**2.3 + B[3]*x**3.6 + B[3]*x**5.2 + B[4]*x**7.3
    x = 1.0-TT/Tc
    y = logppc
    
    linear = Model(f_p)
    mydata = Data(x, y)
    myodr = ODR(mydata, linear, beta0=[0, 0, 0, 0, 0])
    myoutput = myodr.run()
    print Ref, 'p SSE =', myoutput.sum_square
    
    x = 1.0-TT/Tc
    y = np.log(np.array(rhoL)/rhoc)
    
    def f_coeffs(n,x):
        # The outer optimizer that will find the coefficients
        
        def f_rhoL(B, x):
            # B is a vector of the parameters.
            # x is an array of the current x values.
            return B[0]*x**n[0] + B[1]*x**n[1] + B[2]*x**n[2] #+ B[3]*x**n[3] #+  + B[4]*x**n[4]
        
        linear = Model(f_rhoL)
        mydata = Data(x, y)
        myodr = ODR(mydata, linear, beta0=[0]*N)
        myoutput = myodr.run()
    
        rho_fit = np.exp(f_rhoL(myoutput.beta,x))*rhoc
        max_abserror = np.max(np.abs((rho_fit/rhoL)-1)*100)
        
        print Ref, 'rhoL SSE =', myoutput.sum_square, myoutput.beta, max_abserror,'%'
        return myoutput.sum_square
    
    import scipy.optimize, random
    
#    n_library = [i/2.0 for i in range(20)]
#    Ngenes = 10
#    
#    genes = [random.sample(n_library,N) for gene in range(Ngenes)]
#        
#    while True:
#        for gene in genes:
#            f_coeffs(gene,x)
    n0 = [0.345,2.2,4.3]
    x = scipy.optimize.minimize(f_coeffs,n0, args=(x,), method = 'Powell')
    
#example1()

example2('REFPROP-R245fa')
