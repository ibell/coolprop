from libcpp.string cimport string

cdef extern from "CPState.h":
    cdef cppclass CoolPropStateClass:

        bint flag_SinglePhase, flag_TwoPhase

        ## Bulk values
        double _rho,_T,_p,_Q,_h,_s, tau, delta

        ## Phase flags
        bint TwoPhase, SinglePhase
        
        ## Nullary Constructor
        CoolPropStateClass() except +
        
        ## Constructor with fluid name
        CoolPropStateClass(string FluidName) except +

        ## Property updater
        ## Uses the indices in CoolProp for the input parameters
        void update(long iInput1, double Value1, long iInput2, double Value2) except +ValueError

        ## Property accessors for saturation parameters directly
        ## These all must be calculated every time if the state is saturated or two-phase
        double rhoL() except +
        double rhoV() except +
        double pL() except +
        double pV() except +
        double TL() except +
        double TV() except +
        ## Derived parameters for the saturation states
        double hL() except +
        double hV() except +
        double sL() except +
        double sV() except +

        ## Bulk properties accessors - temperature and density are directly calculated every time
        ## All other parameters are calculated on an as-needed basis
        ## If single-phase, just plug into the EOS, otherwise need to do two-phase analysis
        double T() except +
        double rho() except +
        double p() except +
        double h() except +
        double s() except +
        double cp() except +
        double cv() except +
        double speed_sound() except +
        double keyed_output(long iOutput) except +
        long phase() except +

        ## ---------------------------------------- 
        ## TTSE LUT things
        ## ----------------------------------------

        ## Enable the TTSE
        void enable_TTSE_LUT() except +
        ## Check if TTSE is enabled
        bint isenabled_TTSE_LUT() except +
        ## Disable the TTSE
        void disable_TTSE_LUT() except +
        ## Interpolate within the TTSE LUT
        double interpolate_in_TTSE_LUT(long iParam, long iInput1, double Input1, long iInput2, double Input2) except +

        ## ---------------------------------------- 
        ## Derivatives of properties
        ## ----------------------------------------

        double dvdp_constT() except +
        double dvdT_constp() except +

        double drhodT_constp() except +
        double drhodh_constp() except +
        double drhodp_consth() except +
        double drhodp_constT() except +
        double d2rhodp2_constT() except +
        double d2rhodTdp() except +
        double d2rhodhdQ() except +
        double d2rhodpdQ() except +
        double d2rhodhdp() except +
        double d2rhodh2_constp() except +
        double d2rhodT2_constp() except +
        
        double dpdrho_constT() except +
        double dpdrho_consth() except +
        double dpdT_constrho() except +
        double dpdT_consth() except +
        double d2pdrho2_constT() except +
        double d2pdrhodT() except +
        double d2pdT2_constrho() except +

        double dhdrho_constT() except +
        double dhdrho_constp() except +
        double dhdT_constrho() except +
        double dhdT_constp() except +
        double dhdp_constT() except +
        double d2hdrho2_constT() except +
        double d2hdrhodT() except +
        double d2hdT2_constrho() except +
        double d2hdT2_constp() except +
        double d2hdp2_constT() except +
        double d2hdTdp() except +

        double dsdrho_constT() except +
        double dsdT_constrho() except +
        double dsdrho_constp() except +
        double dsdT_constp() except +
        double dsdp_constT() except +
        double d2sdrho2_constT() except +
        double d2sdrhodT() except +
        double d2sdT2_constrho() except +
        double d2sdT2_constp() except +
        double d2sdp2_constT() except +
        double d2sdTdp() except +

        ## ---------------------------------------- 
        ## Derivatives along the saturation curve
        ## ----------------------------------------
        
        ## Derivative of temperature w.r.t. pressure along saturation curve
        double dTdp_along_sat() except +ValueError
        ## Second derivative of temperature w.r.t. pressure along saturation curve
        double d2Tdp2_along_sat() except +ValueError
        ## Partial derivative w.r.t. pressure of dTdp along saturation curve
        double ddp_dTdp_along_sat() except +ValueError
        ## Partial derivative w.r.t. temperature of dTdp along saturation curve
        double ddT_dTdp_along_sat() except +ValueError

        double dhdp_along_sat_vapor() except +ValueError
        double dhdp_along_sat_liquid() except +ValueError
        double d2hdp2_along_sat_vapor() except +ValueError
        double d2hdp2_along_sat_liquid() except +ValueError

        double dsdp_along_sat_vapor() except +ValueError
        double dsdp_along_sat_liquid() except +ValueError
        double d2sdp2_along_sat_vapor() except +ValueError
        double d2sdp2_along_sat_liquid() except +ValueError

        double drhodp_along_sat_vapor() except +ValueError
        double drhodp_along_sat_liquid() except +ValueError
        double d2rhodp2_along_sat_vapor() except +ValueError
        double d2rhodp2_along_sat_liquid() except +ValueError

        double drhodT_along_sat_vapor() except +ValueError
        double drhodT_along_sat_liquid() except +ValueError

        ## Clear out all the cached values
        void clear_cache()

        ## ---------------------------------------- 
        ## Helmholtz Energy Derivatives
        ## ----------------------------------------

        double phi0(double tau, double delta)
        double dphi0_dDelta(double tau, double delta)
        double dphi0_dTau(double tau, double delta)
        double d2phi0_dDelta2(double tau, double delta)
        double d2phi0_dDelta_dTau(double tau, double delta)
        double d2phi0_dTau2(double tau, double delta)
        double d3phi0_dDelta3(double tau, double delta)
        double d3phi0_dDelta2_dTau(double tau, double delta)
        double d3phi0_dDelta_dTau2(double tau, double delta)
        double d3phi0_dTau3(double tau, double delta)

        double phir(double tau, double delta)
        double dphir_dDelta(double tau, double delta)
        double dphir_dTau(double tau, double delta)
        double d2phir_dDelta2(double tau, double delta)
        double d2phir_dDelta_dTau(double tau, double delta)
        double d2phir_dTau2(double tau, double delta)
        double d3phir_dDelta3(double tau, double delta)
        double d3phir_dDelta2_dTau(double tau, double delta)
        double d3phir_dDelta_dTau2(double tau, double delta)
        double d3phir_dTau3(double tau, double delta)
        
cdef extern from "CoolPropDLL.h":
    int _set_reference_stateS "set_reference_stateS"(char *, char *)
    int _set_reference_stateD "set_reference_stateD"(char *, double T, double rho, double h0, double s0)
    
cdef extern from "Units.h":    
    double _fromSI "convert_from_SI_to_unit_system"(string input, double value, string new_system)
    double _toSI "convert_from_unit_system_to_SI"(string input, double value, string old_system)
    
cdef extern from "CoolProp.h":
    void _set_standard_unit_system "set_standard_unit_system"(int unit_system)
    double _Props "Props"(string Output, char Name1, double Prop1, char Name2, double Prop2, string Ref)
    double _PropsSI "PropsSI"(char* Output, char* Name1, double Prop1, char* Name2, double Prop2, char* Ref)
    double _IProps "IProps"(long Output, long Name1, double Prop1, long Name2, double Prop2, long Ref)
    double _Props1 "Props1"(string Ref, string Output)
    double _Props1SI "Props1SI"(char* Ref, char* Output)
    string _Phase "Phase"(char *Fluid, double T, double p)
    string _Phase_Tp "Phase_Tp"(string Fluid, double T, double p)
    string _Phase_Trho "Phase_Trho"(string Fluid, double T, double rho)
    double _DerivTerms "DerivTerms"(char* Term, double T, double rho, char* Ref)
    
    # Conversion functions
    double _F2K "F2K"(double T_F)
    double _K2F "K2F"(double T)
    
    string _get_global_param_string "get_global_param_string"(string ParamName)
    string _get_fluid_param_string "get_fluid_param_string"(string ParamName, string FluidName)
    
    long _set_phase "set_phase" (string phase)
    long _get_Fluid_index "get_Fluid_index" (string Fluid)
    long _get_param_index "get_param_index" (string param)
    string _get_index_units "get_index_units" (long index)
    
    string _get_TTSE_mode "get_TTSE_mode"(string Fluid)
    int _set_TTSE_mode "set_TTSE_mode"(char* Fluid, char* Value)
    
    #Ancillary equations
    double _psatL_anc "psatL_anc"(char*Fluid, double T)
    double _psatV_anc "psatV_anc"(char*Fluid, double T)
    double _rhosatL_anc "rhosatL_anc"(char*Fluid, double T)
    double _rhosatV_anc "rhosatV_anc"(char*Fluid, double T)
    
    double _viscosity_dilute "viscosity_dilute"(char* FluidName, double T)
    double _viscosity_residual "viscosity_residual"(char* FluidName, double T, double rho)
    double _conductivity_background "conductivity_background"(char* FluidName, double T, double rho)
    double _conductivity_critical "conductivity_critical"(char* FluidName, double T, double rho)
    double _conformal_Trho "conformal_Trho"(char* FluidName, char* ReferenceFluidName, double T, double rho, double *Tconform, double *rhoconform)
    
    string _add_REFPROP_fluid "add_REFPROP_fluid"(string FluidName)
    
    int _get_debug_level "get_debug_level"()
    void _set_debug_level "set_debug_level"(int level)
    
    string _get_BibTeXKey "get_BibTeXKey"(string Ref, string key)
    
    # Convenience functions
    int _IsFluidType "IsFluidType"(char* Ref, char* Type)
    
    # Enable the TTSE
    bint _enable_TTSE_LUT "enable_TTSE_LUT"(char *FluidName)
    # Check if TTSE is enabled
    bint _isenabled_TTSE_LUT "isenabled_TTSE_LUT"(char *FluidName)
    # Disable the TTSE
    bint _disable_TTSE_LUT "disable_TTSE_LUT"(char *FluidName)
    
    # Enable the writing of TTSE tables to file for this fluid
    bint _enable_TTSE_LUT_writing "enable_TTSE_LUT_writing"(char *FluidName)
    # Check if the writing of TTSE tables to file is enabled
    bint _isenabled_TTSE_LUT_writing "isenabled_TTSE_LUT_writing"(char *FluidName)
    # Disable the writing of TTSE tables to file for this fluid
    bint _disable_TTSE_LUT_writing "disable_TTSE_LUT_writing"(char *FluidName)
    
    # Over-ride the default size of both of the saturation LUT
    bint _set_TTSESat_LUT_size "set_TTSESat_LUT_size"(char *FluidName, int)
    # Over-ride the default size of the single-phase LUT
    bint _set_TTSESinglePhase_LUT_size "set_TTSESinglePhase_LUT_size"(char *FluidName, int Np, int Nh)
    # Over-ride the default range of the single-phase LUT
    bint _set_TTSESinglePhase_LUT_range "set_TTSESinglePhase_LUT_range"(char *FluidName, double hmin, double hmax, double pmin, double pmax)
    # Get the current range of the single-phase LUT
    bint _get_TTSESinglePhase_LUT_range "get_TTSESinglePhase_LUT_range"(char *FluidName, double *hmin, double *hmax, double *pmin, double *pmax)
        
cdef extern from "HumidAirProp.h":
    double _HAProps "HAProps"(char *OutputName, char *Input1Name, double Input1, char *Input2Name, double Input2, char *Input3Name, double Input3)
    double _HAProps_Aux "HAProps_Aux"(char* Name,double T, double p, double W, char *units)
    double _cair_sat "cair_sat"(double T)
       
cdef class State:
    cdef public PureFluidClass PFC      # hold a PureFluidClass which is the Cython wrapper of CoolPropStateClass
    cdef readonly bint hasLiquid
    cdef readonly bytes Liquid, Fluid, phase
    cdef long iFluid,iParam1,iParam2,iOutput
    cdef double T_, rho_, p_, xL
    cdef readonly bint is_CPFluid
    
    cpdef set_Fluid(self, str Fluid)
    cpdef speed_test(self, int N)
    cpdef update(self,dict params, double xL=*)
    cpdef update_ph(self, double p, double h)
    cpdef update_Trho(self, double T, double rho)
    cpdef State copy(self)
    cpdef double Props(self, long iOutput) except *
    cpdef long Phase(self) except *
    
    cpdef double get_Q(self) except *
    cpdef double get_T(self) except *
    cpdef double get_p(self) except *
    cpdef double get_h(self) except *
    cpdef double get_rho(self) except *
    cpdef double get_s(self) except *
    cpdef double get_u(self) except *
    cpdef double get_visc(self) except *
    cpdef double get_cond(self) except *
    cpdef double get_cp(self) except *
    cpdef double get_cp0(self) except *
    cpdef double get_cv(self) except *
    cpdef double get_MM(self) except *
    cpdef double get_dpdT(self) except *
    cpdef double get_speed_sound(self) except *
    cpdef double get_Tsat(self, double Q = *)
    cpdef double get_subcooling(self)
    cpdef double get_superheat(self)
    
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
    cpdef double keyed_output(self, long iOutput) except *
    cpdef long phase(self) except *
    
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
    cpdef bint isenabled_TTSE_LUT(self) # Check if TTSE is enabled
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
    
    cpdef double drhodT_constp(self)
    cpdef double drhodp_constT(self)
    cpdef double d2rhodp2_constT(self)
    cpdef double d2rhodTdp(self)
    cpdef double d2rhodT2_constp(self)
    cpdef double d2rhodhdQ(self)
    cpdef double d2rhodpdQ(self)
    cpdef double d2rhodhdp(self)
    cpdef double d2rhodh2_constp(self)

    cpdef double dpdrho_constT(self)
    cpdef double dpdrho_consth(self)
    cpdef double dpdT_constrho(self)
    cpdef double dpdT_consth(self)
    cpdef double d2pdrho2_constT(self)
    cpdef double d2pdrhodT(self)
    cpdef double d2pdT2_constrho(self)

    cpdef double dhdrho_constT(self)
    cpdef double dhdrho_constp(self)
    cpdef double dhdT_constrho(self)
    cpdef double dhdT_constp(self)
    cpdef double dhdp_constT(self)
    cpdef double d2hdrho2_constT(self)
    cpdef double d2hdrhodT(self)
    cpdef double d2hdT2_constrho(self)
    cpdef double d2hdT2_constp(self)
    cpdef double d2hdp2_constT(self)
    cpdef double d2hdTdp(self)

    cpdef double dsdrho_constT(self)
    cpdef double dsdT_constrho(self)
    cpdef double dsdrho_constp(self)
    cpdef double dsdT_constp(self)
    cpdef double dsdp_constT(self)
    cpdef double d2sdrho2_constT(self)
    cpdef double d2sdrhodT(self)
    cpdef double d2sdT2_constrho(self)
    cpdef double d2sdT2_constp(self)
    cpdef double d2sdp2_constT(self)
    cpdef double d2sdTdp(self)
