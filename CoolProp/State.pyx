#cython: embedsignature = True

cdef extern from "CoolProp.h":
    double _Props "Props" (char*,char,double,char,double,char*)
    void UseSinglePhaseLUT(bool)
    double DerivTerms(char *, double, double, char*)
    char * get_errstringc()
    void get_errstring(char*)
    int _set_1phase_LUT_params "set_1phase_LUT_params" (char*,int,int,double,double,double,double)
    void _debug "debug" (int)
    
from libc.math cimport pow, sin, cos, exp
from math import pow as pow_
cdef bint _LUT_Enabled
import CoolProp as CP

cpdef cmath_speed_test(float x, long N):
    """
    Run a few tests at the c++ level to test how fast math functions are
    
    Parameters
    ----------
    x : float
        Input value for trig functions
    N : int
        Number of calls to make
    """
    from time import clock
    cdef int i
    cdef double y
    
    t1=clock()
    for i in range(N):
        sin(x)
    t2=clock()
    print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,'sin',(t2-t1)/N*1e6)
    
    t1=clock()
    for i in range(N):
        cos(x)
    t2=clock()
    print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,'cos',(t2-t1)/N*1e6)
    
    t1=clock()
    for i in range(N):
        exp(x)
    t2=clock()
    print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,'exp',(t2-t1)/N*1e6)
    
    t1=clock()
    y=0
    for i in range(N):
        y += pow(<double>x,<int>4)
    t2=clock()
    print y
    print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,'pow(int)',(t2-t1)/N*1e6)
    
    t1=clock()
    y=0
    for i in range(N):
        y += pow(<double>x,<double>4)
    t2=clock()
    print y
    print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,'pow(double)',(t2-t1)/N*1e6)

cpdef int set_1phase_LUT_params(bytes Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax):
    """
    Set the 
    """
    #Set the LUT parameters in CoolProp copy of the LUT table
    CP.set_1phase_LUT_params(Ref, nT, np, Tmin, Tmax, pmin, pmax)
    
    #Set the LUT parameters in the local LUT copy build into the State module
    _set_1phase_LUT_params(Ref, nT, np, Tmin, Tmax, pmin, pmax)

    return 0

cpdef debug(int level):
    """
    Sets the debug level
    
    Parameters
    ----------
    level : int
        Flag indicating how verbose the debugging should be.
            0 : no debugging output
            ...
            ...
            10 : very annoying debugging output - every function call debugged
    """
    _debug(level)

cpdef LUT(bint LUTval):
    """
    
    LUTval : boolean
        If ``True``, turn on the use of lookup table.  Parameters must have 
        already been set through the use of set_1phase_LUT_params

    """
    if LUTval:
        _LUT_Enabled = True
        print 'Turning on singlephase LUT'
        UseSinglePhaseLUT(True)
    else:
        _LUT_Enabled = False
        UseSinglePhaseLUT(False)
        
cpdef double Props(bytes Parameter, bytes param1, float value1, bytes param2, float value2, bytes Fluid):
    """
    Expose the Props() function.  Uses the same call signature as the Props() function in CoolProp.CoolProp
    """
    cdef char _param1 = param1[0]
    cdef char _param2 = param2[0]  
    return _Props(Parameter, _param1, value1, _param2, value2, Fluid)

#from CoolProp import Props,UseSinglePhaseLUT,DerivTerms
cdef class State: 
    """
    A class that contains all the code that represents a thermodynamic state
    """
    
    def __init__(self,bytes Fluid, dict StateDict, double xL=-1.0, Liquid=''):
        self.Fluid=Fluid
        self.xL=xL
        self.Liquid=Liquid
        
        #Parse the inputs provided
        self.update(StateDict)
            
    def __reduce__(self):
        d={}
        d['xL']=self.xL
        d['Liquid']=self.Liquid
        d['Fluid']=self.Fluid
        d['T']=self.T_
        d['rho']=self.rho_
        return rebuildState,(d,)
          
    cpdef update(self,dict params, double xL=-1.0):
        """
        *params* is a list(or tuple) of strings that represent the parameters 
        that have been updated and will be used to fix the rest of the state. 
        ['T','P'] for temperature and pressure for instance
        """
            
        cdef double p
        cdef char* pchar
        cdef bytes errstr
        
        # If no value for xL is provided, it will have a value of -1 which is 
        # impossible, so don't update xL
        if xL > 0:
            #There is liquid
            self.xL=xL
            self.hasLiquid=True
        else:
            #There's no liquid
            self.xL=0.0
            self.hasLiquid=False
        
        #You passed in a dictionary, use the values to update the state
        if 'T' not in params:
            raise AttributeError('T must be provided in params dict in State.update')
            
        #Consume the 'T' key since it is required (TODO?)
        self.T_=float(params.pop('T'))
            
        #Given temperature and pressure, determine density of gas 
        # (or gas and oil if xL is provided)
        if abs(self.xL)<=1e-15:
            #Get the density if T,P provided, or pressure if T,rho provided
            if 'P' in params:
                self.p_=params['P']
                #Explicit type conversion
                pchar='D'
                rho = _Props(pchar,'T',self.T_,'P',self.p_,self.Fluid)
                if abs(rho)<1e90:
                    self.rho_=rho
                else:
                    errstr = get_errstringc()
                    raise ValueError(errstr)
            elif 'D' in params:
                self.rho_=params['D']
                #Explicit type conversion
                pchar='P'
                p = _Props(pchar,'T',self.T_,'D',self.rho_,self.Fluid)
                if abs(p)<1e90:
                    self.p_=p
                else:
                    errstr = get_errstringc()
                    raise ValueError(errstr+str(params))
            else:
                raise KeyError("Dictionary must contain the key 'T' and one of 'P' or 'D'")
            
        elif self.xL>0 and self.xL<=1:
            raise ValueError('Need more code here')
        else:
            raise ValueError('xL must be between 0 and 1')
        
    cpdef double get_MM(self):
        return _Props('M','T',0,'D',0,self.Fluid)
    
    property LUT:
        def __get__(self):
            return self.LUT
    
    cpdef double get_rho(self): 
        return self.rho_
    property rho:
        def __get__(self):
            return self.rho_
            
    cpdef double get_p(self): 
        return self.p_
    property p:
        def __get__(self):
            return self.p_
    
    cpdef double get_T(self): 
        return self.T_
    property T:
        def __get__(self):
            return self.T_
    
    cpdef double get_h(self): 
        return _Props("H",'T',self.T_,'D',self.rho_,self.Fluid)
    property h:
        def __get__(self):
            return self.get_h()
          
    cpdef double get_u(self): 
        return _Props("U",'T',self.T_,'D',self.rho_,self.Fluid)
    property u:
        def __get__(self):
            return self.get_u()
            
    cpdef double get_s(self): 
        return _Props("S",'T',self.T_,'D',self.rho_,self.Fluid)            
    property s:
        def __get__(self):
            return self.get_s()
    
    cpdef double get_cp0(self):
        return _Props("C0",'T',self.T_,'D',self.rho_,self.Fluid)
    
    cpdef double get_cp(self): 
        return _Props("C",'T',self.T_,'D',self.rho_,self.Fluid)
    property cp:
        def __get__(self):
            return self.get_cp()
            
    cpdef double get_cv(self): 
        return _Props("O",'T',self.T_,'D',self.rho_,self.Fluid)
    property cv:
        def __get__(self):
            return self.get_cv()
            
    cpdef double get_visc(self):
        return _Props('V','T',self.T_,'D',self.rho_,self.Fluid)
    property visc:
        def __get__(self):
            return self.get_visc()

    cpdef double get_cond(self):
        return _Props('L','T',self.T_,'D',self.rho_,self.Fluid)    
    property k:
        def __get__(self):
            return self.get_cond()
            
    property Prandtl:
        def __get__(self):
            return self.cp * self.visc / self.k
            
    cpdef double get_dpdT(self):
        return _Props("dpdT",'T',self.T_,'D',self.rho_,self.Fluid)
    property dpdT:
        def __get__(self):
            return self.get_dpdT()
        
    cpdef speed_test(self, int N):
        from time import clock
        cdef int i
        cdef char * k
        cdef char * Fluid = self.Fluid 
        print 'Direct c++ call to CoolProp without the Python call layer'
        print "'M' involves basically no computational effort and is a good measure of the function call overhead"
        keys = ['H','P','S','U','C','O','V','L','M','C0']
        for key in keys:
            t1=clock()
            for i in range(N):
                _Props(key,'T',self.T_,'D',self.rho_,Fluid)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,key,(t2-t1)/N*1e6)
    
    def __str__(self):
        units={'T': 'K', 
               'p': 'kPa', 
               'rho': 'kg/m^3',
               'h':'kJ/kg',
               'u':'kJ/kg',
               's':'kJ/kg/K',
               'visc':'Pa-s',
               'k':'kW/m/K',
               'cp':'kJ/kg/K',
               'cv':'kJ/kg/K',
               'dpdT':'kPa/K'}
        s=''
        for k in ['T','p','rho','h','u','s','visc','k','cp','cv','dpdT','Prandtl']:
            if k in units:
                s+=k+' = '+str(getattr(self,k))+' '+units[k]+'\n'
            else:
                s+=k+' = '+str(getattr(self,k))+' NO UNITS'+'\n'
        return s.rstrip()
        
    cpdef copy(self):
        cdef double T = self.T_*(1.0+1e-20)
        cdef double rho = self.rho_*(1.0+1e-20)
        ST=State(self.Fluid,{'T':T,'D':rho})
        return ST
    
def rebuildState(d):
    S=State(d['Fluid'],{'T':d['T'],'D':d['rho']})
    S.xL = d['xL']
    S.Liquid=d['Liquid']
    return S


