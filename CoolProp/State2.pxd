from libcpp cimport bool 
from libcpp.string cimport string

cdef extern from "CPState.h":    
    cdef cppclass CoolPropStateClass:

        bool flag_SinglePhase, flag_TwoPhase

        ## Bulk values
        double _rho,_T,_p,_Q,_h,_s, tau, delta

        ## Phase flags
        bool TwoPhase, SinglePhase
        
        ## Constructor with fluid name
        CoolPropStateClass()
        
        ## Constructor with fluid name
        CoolPropStateClass(string FluidName) except +

        ## Property updater
        ## Uses the indices in CoolProp for the input parameters
        void update(long iInput1, double Value1, long iInput2, double Value2)

        ## Property accessors for saturation parameters directly
        ## These all must be calculated every time if the state is saturated or two-phase
        double rhoL()
        double rhoV()
        double pL()
        double pV()
        double TL()
        double TV()
        ## Derived parameters for the saturation states
        double hL()
        double hV()
        double sL()
        double sV()

        ## Bulk properties accessors - temperature and density are directly calculated every time
        ## All other parameters are calculated on an as-needed basis
        ## If single-phase, just plug into the EOS, otherwise need to do two-phase analysis
        double T()
        double rho()
        double p()
        double h()
        double s()
        double cp()
        double cv()
        double speed_sound()

        ## ----------------------------------------	
        ## TTSE LUT things
        ## ----------------------------------------

        ## Enable the TTSE
        void enable_TTSE_LUT()
        ## Check if TTSE is enabled
        bool isenabled_TTSE_LUT()
        ## Disable the TTSE
        void disable_TTSE_LUT()
        ## Build the TTSE LUT
        bool build_TTSE_LUT()
        ## Interpolate within the TTSE LUT
        double interpolate_in_TTSE_LUT(long iParam, long iInput1, double Input1, long iInput2, double Input2)

        ## ----------------------------------------	
        ## Derivatives of properties
        ## ----------------------------------------

        double dvdp_constT()
        double dvdT_constp()

        double drhodT_constp()
        double drhodh_constp()
        double drhodp_consth()
        double drhodp_constT()
        double d2rhodp2_constT()
        double d2rhodTdp()
        double d2rhodT2_constp()
        
        double dpdrho_constT()
        double dpdrho_consth()
        double dpdT_constrho()
        double dpdT_consth()
        double d2pdrho2_constT()
        double d2pdrhodT()
        double d2pdT2_constrho()

        double dhdrho_constT()
        double dhdrho_constp()
        double dhdT_constrho()
        double dhdT_constp()
        double dhdp_constT()
        double d2hdrho2_constT()
        double d2hdrhodT()
        double d2hdT2_constrho()
        double d2hdT2_constp()
        double d2hdp2_constT()
        double d2hdTdp()

        double dsdrho_constT()
        double dsdT_constrho()
        double dsdrho_constp()
        double dsdT_constp()
        double dsdp_constT()
        double d2sdrho2_constT()
        double d2sdrhodT()
        double d2sdT2_constrho()
        double d2sdT2_constp()
        double d2sdp2_constT()
        double d2sdTdp()

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