from libcpp.string cimport string

cdef extern from "CoolProp.h":
    
    double _Props "Props"(string Output, char Name1, double Prop1, char Name2, double Prop2, string Ref)
    double _IProps "IProps"(long Output, long Name1, double Prop1, long Name2, double Prop2, long Ref)
    double _Props1 "Props"(string Ref, string Output)
    string _Phase "Phase"(char *Fluid, double T, double p)
    string _Phase_Tp "Phase_Tp"(string Fluid, double T, double p)
    string _Phase_Trho "Phase_Trho"(string Fluid, double T, double rho)
    double _DerivTerms "DerivTerms"(char *Term, double T, double rho, char * Ref)
    
    # Conversion functions
    double _F2K "F2K"(double T_F)
    double _K2F "K2F"(double T)
    
    string _FluidsList "FluidsList"()
    string _get_REFPROPname "get_REFPROPname"(string Ref)
    string _get_errstring "get_errstring"()
    string _get_aliases "get_aliases"(string Ref)
    long _set_phase "set_phase" (string phase)
    long _get_Fluid_index "get_Fluid_index" (string Fluid)
    long _get_param_index "get_param_index" (string param)
    string _get_index_units "get_index_units" (long index)
    char * get_errstringc()
    
    #Ancillary equations
    double _psatL_anc "psatL_anc"(char*Fluid, double T)
    double _psatV_anc "psatV_anc"(char*Fluid, double T)
    double _rhosatL_anc "rhosatL_anc"(char*Fluid, double T)
    double _rhosatV_anc "rhosatV_anc"(char*Fluid, double T)
    
    double _viscosity_dilute "viscosity_dilute"(char* FluidName, double T, double e_k, double sigma)
    
    
    int _get_debug "get_debug"()
    void _debug "debug"(int level)
    
    string _get_EOSReference "get_EOSReference"(string Ref)
    string _get_TransportReference "get_TransportReference"(string Ref)
    
    # Convenience functions
    int _IsFluidType "IsFluidType"(char *Ref, char *Type)
    
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

cdef extern from "CPState.h":
    cdef cppclass CoolPropStateClass:
        CoolPropStateClass() except +
        CoolPropStateClass(string FluidName) except +
        double T()
        double rho()
        double p()
        void update(long iInput1, double Value1, long iInput2, double Value2)
        double keyed_output(long iOutput)
        long phase()
        double Q()
        
cdef extern from "HumidAirProp.h":
    double _HAProps "HAProps"(char *OutputName, char *Input1Name, double Input1, char *Input2Name, double Input2, char *Input3Name, double Input3)
    double _HAProps_Aux "HAProps_Aux"(char* Name,double T, double p, double W, char *units)
       
cdef class State:
    cdef CoolPropStateClass *CPS      # hold a C++ instance which we're wrapping
    cdef readonly bint hasLiquid
    cdef readonly bytes Liquid, Fluid, phase
    cdef long iFluid,iParam1,iParam2,iOutput
    cdef double T_, rho_, p_, xL
    cdef readonly bint is_CPFluid
    
    cpdef speed_test(self, int N)
    cpdef update(self,dict params, double xL=*)
    cpdef update_ph(self, double p, double h)
    cpdef update_Trho(self, double T, double rho)
    cpdef copy(self)
    cpdef double Props(self, long iOutput)
    cpdef long Phase(self)
    
    cpdef double get_Q(self)
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
    
