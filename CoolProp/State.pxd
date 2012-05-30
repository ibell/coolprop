
cpdef LUT(bint LUTkey)
cpdef int set_1phase_LUT_params_(bytes Ref, int nT,int np,double Tmin,double Tmax,double pmin,double pmax)
cpdef debug_level(int level)

cdef class State:
    cdef readonly bint hasLiquid
    cdef readonly bytes Liquid, Fluid
    cdef double T_, rho_, p_, xL
    cpdef speed_test(self, int N)
    cpdef update(self,dict params, double xL=?)
    cpdef copy(self)
    