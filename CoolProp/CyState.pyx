
cdef class PureFluidClass:
    def __cinit__(self, string name):
        self.CPS = CoolPropStateClass(name)
    cpdef update(self, long iInput1, double Value1, long iInput2, double Value2):
        self.CPS.update(iInput1,Value1,iInput2,Value2)
    
    cpdef double rhoL(self): 
        return self.CPS.rhoL()
    cpdef double rhoV(self):
        return self.CPS.rhoV()
    cpdef double pL(self):
        return self.CPS.pL()
    cpdef double pV(self):
        return self.CPS.pV()
    cpdef double TL(self):
        return self.CPS.TL()
    cpdef double TV(self):
        return self.CPS.TV()
    cpdef double sL(self):
        return self.CPS.sL()
    cpdef double sV(self):
        return self.CPS.sV()
    cpdef double hL(self):
        return self.CPS.hL()
    cpdef double hV(self):
        return self.CPS.hV()
    
    ## ----------------------------------------	
    ##        Fluid property accessors
    ## ----------------------------------------
    
    cpdef double T(self):
        return self.CPS.T()    
    cpdef double rho(self):
        return self.CPS.rho()
    cpdef double p(self):
        return self.CPS.p()
    cpdef double h(self):
        return self.CPS.h()
    cpdef double s(self):
        return self.CPS.s()
    cpdef double cp(self):
        return self.CPS.cp()
    cpdef double cv(self):
        return self.CPS.cv()
    cpdef double speed_sound(self):
        return self.CPS.speed_sound()
    
    cpdef double keyed_output(self, long iOutput):
        return self.CPS.keyed_output(iOutput)

    ## ----------------------------------------	
    ##        TTSE LUT things
    ## ----------------------------------------

    # Enable the TTSE
    cpdef enable_TTSE_LUT(self):
        self.CPS.enable_TTSE_LUT()
        
    # Check if TTSE is enabled
    cpdef bint isenabled_TTSE_LUT(self):
        return self.CPS.isenabled_TTSE_LUT()
        
    # Disable the TTSE
    cpdef disable_TTSE_LUT(self):
        self.CPS.disable_TTSE_LUT()
        

    cpdef double dTdp_along_sat(self):
        return self.CPS.dTdp_along_sat()
    cpdef double d2Tdp2_along_sat(self):
        return self.CPS.d2Tdp2_along_sat()

    cpdef double dhdp_along_sat_vapor(self):
        return self.CPS.dhdp_along_sat_vapor()
    cpdef double dhdp_along_sat_liquid(self):
        return self.CPS.dhdp_along_sat_liquid()
    cpdef double d2hdp2_along_sat_vapor(self):
        return self.CPS.d2hdp2_along_sat_vapor()
    cpdef double d2hdp2_along_sat_liquid(self):
        return self.CPS.d2hdp2_along_sat_liquid()

    cpdef double dsdp_along_sat_vapor(self):
        return self.CPS.dsdp_along_sat_vapor()
    cpdef double dsdp_along_sat_liquid(self):
        return self.CPS.dsdp_along_sat_liquid()
    cpdef double d2sdp2_along_sat_vapor(self):
        return self.CPS.d2sdp2_along_sat_vapor()
    cpdef double d2sdp2_along_sat_liquid(self):
        return self.CPS.d2sdp2_along_sat_liquid()

    cpdef double drhodp_along_sat_vapor(self):
        return self.CPS.drhodp_along_sat_vapor()
    cpdef double drhodp_along_sat_liquid(self):
        return self.CPS.drhodp_along_sat_liquid()
    cpdef double d2rhodp2_along_sat_vapor(self):
        return self.CPS.d2rhodp2_along_sat_vapor()
    cpdef double d2rhodp2_along_sat_liquid(self):
        return self.CPS.d2rhodp2_along_sat_liquid()

    cpdef double drhodT_along_sat_vapor(self):
        return self.CPS.drhodT_along_sat_vapor()
    cpdef double drhodT_along_sat_liquid(self):
        return self.CPS.drhodT_along_sat_liquid()