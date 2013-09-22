# A template the for the .h file
PPF_h_template = """
#ifndef {RefUpper:s}_H
#define {RefUpper:s}_H

    class {Ref:s}Class : public Fluid{{

    public:
        {Ref:s}Class();
        ~{Ref:s}Class(){{}};
        double psatL(double);
        double psatV(double);
        double rhosatL(double);
        double rhosatV(double);
    }};
#endif
"""

PPF_cpp_template = """
#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include "math.h"
#include "stdio.h"
#include <string.h>
#include "CoolProp.h"
#include "FluidClass.h"
#include "{Ref:s}.h"

{Ref:s}Class::{Ref:s}Class()
{{
    static double a[]={{{acoeffs:s}}};
    static double b[]={{{bcoeffs:s}}};

    static double N[]={{{Ncoeffs:s}}};
    static double t[]={{{tcoeffs:s}}};
    static double d[]={{{dcoeffs:s}}};
    static double l[]={{{Lcoeffs:s}}};

    phirlist.push_back(new phir_power(N,d,t,l,1,{N_phir:d}-1,{N_phir:d}));

    phi0list.push_back(new phi0_lead(0,0));
    phi0list.push_back(new phi0_logtau(-1.0));
    phi0list.push_back(new phi0_power(a[0],b[0]));
    phi0list.push_back(new phi0_cp0_Planck_Einstein(a,b,1,{N_cp0:d}-1,{N_cp0:d}));

    // Other fluid parameters
    params.molemass = {molemass:g}; //[kg/m^3]
    params.Ttriple = {Ttriple:g}; //[K]
    params.accentricfactor = {accentric:g}; //[-]
    params.R_u = 8.314472;
    isPure = false;
    
    // Critical parameters
    crit.rho = {rhocrit:g};
    crit.p = PressureUnit({pcrit:g},UNIT_KPA);
    crit.T = {Tcrit:g};
    crit.v = 1.0/crit.rho;

    // Limits of EOS
    limits.Tmin = params.Ttriple;

    name.assign("{Ref:s}");
}}
{pL:s}
{pV:s}
{rhoL:s}
{rhoV:s}
"""

