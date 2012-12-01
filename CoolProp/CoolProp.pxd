from libcpp.string cimport string

cdef class State:
    cdef readonly bint hasLiquid
    cdef readonly bytes Liquid, Fluid, phase
    cdef long iFluid,iParam1,iParam2,iOutput
    cdef double T_, rho_, p_, xL
    cdef readonly bint is_CPFluid
    
    cpdef speed_test(self, int N)
    cpdef update(self,dict params, double xL=*)
    cpdef update_Trho(self, double T, double rho)
    cpdef copy(self)
    cpdef double Props(self, long iOutput)
    
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
    
