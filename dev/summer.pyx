
import numpy as np
cimport numpy as np

cpdef np.ndarray[np.float_t, ndim=1] sum_function(np.ndarray[np.float_t, ndim=1] B, np.ndarray[np.float_t, ndim=1] x, np.ndarray[np.float_t, ndim=1] n):

    cdef np.ndarray[np.float_t, ndim=1] out = np.zeros_like(x)
    cdef int i,j
    cdef double do
    cdef int Nx = len(x)
    cdef int Nn = len(n)
    
    for i in xrange(Nx):
        for j in range(Nn):
            out[i] += B[j]*x[i]**n[j]
            
    return out
        