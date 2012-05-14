
cdef class State:
    cdef readonly bint hasLiquid
    cdef readonly bytes Liquid, Fluid
    cdef double T_, rho_, p_, xL
    cpdef speed_test(self, int N)
    cpdef update(self,dict params, double xL=?)
    cpdef copy(self)
    