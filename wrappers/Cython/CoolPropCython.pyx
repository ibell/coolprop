import cython
import numpy as np
cimport numpy as np

cdef extern from "CoolProp.h":
    double Props(char *,char,double,char,double,char *)
    void UseSinglePhaseLUT(int)
    
cpdef double nothing(void):
    UseSinglePhaseLUT(1)
    cdef double retval = 1.0
    return retval

cpdef PropsV(Output, bytes Name1, np.ndarray[np.float_t, ndim=1] Value1, bytes Name2, np.ndarray[np.float_t, ndim=1] Value2, bytes Ref):
    """
    Some documentation of this function
    """
    cdef int len1 = Value1.shape[0]
    cdef int len2 = Value2.shape[0]
    cdef int N = max(len1,len2)
    cpdef np.ndarray[np.float_t, ndim=1] OutArray=np.zeros((N),dtype=np.float)
    cdef char Name1_c=Name1[0]
    cdef char Name2_c=Name2[0]
    cdef char *Output_c=Output
    cdef char *Ref_c=Ref
    cdef int i
    cdef double Val1,Val2,OutVal
    
    if len1==len2:
        for i in range(N):
            Val1=Value1[i]
            Val2=Value2[i]
            OutArray[i]=Props(Output_c,Name1_c,Val1,Name2_c,Val2,Ref_c)
    elif len1==1:
        Val1=Value1[0]
        for i in range(N):
            Val2=Value2[i]
            OutArray[i]=Props(Output_c,Name1_c,Val1,Name2_c,Val2,Ref_c)
    elif len2==1:
        Val2=Value2[0]
        for i in range(N):
            Val1=Value1[i]
            OutVal=Props(Output_c,Name1_c,Val1,Name2_c,Val2,Ref_c)
            OutArray[i]=OutVal
    else:
        raise ValueError('Vectors must be the same length, or one must be a single value')
    return OutArray
    
cpdef np.ndarray[np.float_t, ndim=1] PropsVempty(bytes Output, bytes Name1, np.ndarray[np.float_t, ndim=1] Value1, bytes Name2, np.ndarray[np.float_t, ndim=1] Value2, bytes Ref):

    cdef int len1 = Value1.shape[0]
    cdef int len2 = Value2.shape[0]
    cdef int N = max(len1,len2)
    cpdef np.ndarray[np.float_t, ndim=1] OutArray=np.zeros((N),dtype=np.float)
    return OutArray