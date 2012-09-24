from libcpp.string cimport string

cdef extern from "CoolProp.h":
    
    double _Props "Props"(string Output, char Name1, double Prop1, char Name2, double Prop2, string Ref)
    double _Props1 "Props"(string Ref, string Output)
    string _Phase "Phase"(char *Fluid, double T, double p)
    double _DerivTerms "DerivTerms"(char *Term, double T, double rho, char * Ref)
    
    # LUT functions
    void _UseSaturationLUT "UseSaturationLUT"(bint OnOff)
    bint _SaturationLUTStatus "SaturationLUTStatus" ()
    void _UseSinglePhaseLUT "UseSinglePhaseLUT"(bint OnOff)
    bint _SinglePhaseLUTStatus "SinglePhaseLUTStatus"()
    
    # Conversion functions
    double _F2K "F2K"(double T_F)
    double _K2F "K2F"(double T)
    
    string _FluidsList "FluidsList"()
    string _get_REFPROPname "get_REFPROPname"(string Ref)
    string _get_errstring "get_errstring"()
    long _get_param_index "get_param_index" (string param)
    string _get_index_units "get_index_units" (long index)
    
    int _get_debug "get_debug"()
    void _debug "debug"(int level)
    
    void _PrintSaturationTable "PrintSaturationTable"(char *FileName, char * Ref, double Tmin, double Tmax)
    
    string _get_EOSReference "get_EOSReference"(string Ref)
    string _get_TransportReference "get_TransportReference"(string Ref)

    int _set_1phase_LUT_params "set_1phase_LUT_params"(char *Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax, bint rebuild)
    void _get_1phase_LUT_params "get_1phase_LUT_params"(int *nT, int *np, double *Tmin, double *Tmax, double *pmin, double *pmax)
    
    # Convenience functions
    int _IsFluidType "IsFluidType"(char *Ref, char *Type)



    
    
