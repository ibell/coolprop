
#Functions defined at the module level
cpdef double Props(bytes Parameter, bytes Param1, float value1, bytes Param2, float value2, bytes Fluid)
cpdef LUT(bint LUTkey)
cpdef int set_1phase_LUT_params(bytes Ref, int nT,int np,double Tmin,double Tmax,double pmin,double pmax)
cpdef debug(int level)
cpdef cmath_speed_test(float x, long N)

cdef class State:
    cdef readonly bint hasLiquid
    cdef readonly bytes Liquid, Fluid
    cdef double T_, rho_, p_, xL
    
    cpdef speed_test(self, int N)
    cpdef update(self,dict params, double xL=?)
    cpdef copy(self)
    
    cpdef double get_T(self)
    cpdef double get_p(self)
    cpdef double get_h(self)
    cpdef double get_rho(self)
    cpdef double get_s(self)
    cpdef double get_u(self)
    cpdef double get_visc(self)
    cpdef double get_cond(self)
    cpdef double get_cp(self)
    cpdef double get_cp0(self)
    cpdef double get_cv(self)
    cpdef double get_MM(self)
    cpdef double get_dpdT(self)
    