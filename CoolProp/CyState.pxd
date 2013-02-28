from libcpp cimport bool 
from libcpp.string cimport string

       
cdef class PureFluidClass:
    cdef CoolPropStateClass CPS     # hold a C++ instance which we're wrapping
    cpdef update(self, long iInput1, double Value1, long iInput2, double Value2)
    cpdef double rhoL(self)
    cpdef double rhoV(self)
    cpdef double pL(self)
    cpdef double pV(self)
    cpdef double TL(self)
    cpdef double TV(self)
    cpdef double sL(self)
    cpdef double sV(self)
    cpdef double hL(self)
    cpdef double hV(self)
    
    ## ---------------------------------------- 
    ##        Fluid property accessors
    ## ----------------------------------------
    
    cpdef double T(self)
    cpdef double rho(self)
    cpdef double p(self)
    cpdef double h(self)
    cpdef double s(self)
    cpdef double cp(self)
    cpdef double cv(self)
    cpdef double speed_sound(self)

    ## ---------------------------------------- 
    ##        TTSE LUT things
    ## ----------------------------------------

    
    cpdef enable_TTSE_LUT(self) # Enable the TTSE
    cpdef bool isenabled_TTSE_LUT(self) # Check if TTSE is enabled
    cpdef disable_TTSE_LUT(self) # Disable the TTSE

    cpdef double dTdp_along_sat(self)
    cpdef double d2Tdp2_along_sat(self)

    cpdef double dhdp_along_sat_vapor(self)
    cpdef double dhdp_along_sat_liquid(self)
    cpdef double d2hdp2_along_sat_vapor(self)
    cpdef double d2hdp2_along_sat_liquid(self)

    cpdef double dsdp_along_sat_vapor(self)
    cpdef double dsdp_along_sat_liquid(self)
    cpdef double d2sdp2_along_sat_vapor(self)
    cpdef double d2sdp2_along_sat_liquid(self)

    cpdef double drhodp_along_sat_vapor(self)
    cpdef double drhodp_along_sat_liquid(self)
    cpdef double d2rhodp2_along_sat_vapor(self)
    cpdef double d2rhodp2_along_sat_liquid(self)

    cpdef double drhodT_along_sat_vapor(self)
    cpdef double drhodT_along_sat_liquid(self)