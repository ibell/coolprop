
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
#include "R407F.h"

R407FClass::R407FClass()
{
    // Constants for the ideal-gas contribution
    static double a[]={0, 1.440311, 2.000416, 5.360509, 3.471306};
    static double b[]={0, 0.256507331599, 696.9608/355.804, 1723.0555/355.804, 3875.0307/355.804};

    // Constants for the residual contribution
    static double N[]={0.5, 0.872372, -1.33046, -0.894762, 0.0589703, 0.000135205};
    static double t[]={0.50402, 0.838012, 0.032941, -0.41349, -0.0448794, -0.108068, -0.0232731};
    static double d[]={0, 1, 1, 1, 3, 7, 1, 2, 5, 1, 1, 4, 2};
    static double l[]={0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3};
    
    // Other fluid parameters
    params.molemass = 82.0583; //[kg/kmol]
    params.Ttriple = 200; //[K]
    params.accentricfactor = 0.7; //[-]
    params.R_u = 8.314472;
    isPure = false;
    
    // Critical parameters
    crit.rho = 477.285;
    crit.p = PressureUnit(4754.61,UNIT_KPA);
    crit.T = 355.804;
    crit.v = 1.0/crit.rho;
    
    phirlist.push_back(new phir_power(N,d,t,l,1,6-1,6));

    phi0list.push_back(new phi0_lead(0, 0));
    phi0list.push_back(new phi0_logtau(-1.0));
    phi0list.push_back(new phi0_cp0_poly(a[1],b[1],crit.T,298.15));
    phi0list.push_back(new phi0_Planck_Einstein(a,b,2,4,4+1));

    // Adjust to the IIR reference state (h=200 kJ/kg, s = 1 kJ/kg for sat. liq at 0C)
    params.HSReferenceState = "IIR";

    // Limits of EOS
    limits.Tmin = params.Ttriple;

    name.assign("R407F");
}

double R407FClass::psatL(double T)
{
    // Maximum absolute error is 0.480139 % between 149.496265 K and 355.803836 K
    const double t[]={0, 1, 2, 3, 8};
    const double N[]={0, 0.097471997278898881, -7.1139435312786263, 0.56195515936733975, -3.723502987238116};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1;i<=4;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}


double R407FClass::psatV(double T)
{
    // Maximum absolute error is 0.129074 % between 149.496265 K and 355.803836 K
    const double t[]={0, 1, 2, 3, 4, 7, 16};
    const double N[]={0, -0.24566974513211498, -6.8095186371477734, 0.010708848678671482, 0.58613267700126404, -4.665597037779353, -4.0021168321556981};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1;i<=6;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}


double R407FClass::rhosatL(double T)
{
    // Maximum absolute error is 0.279970 % between 149.496265 K and 354.803846 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.5};
    const double N[] = {0, -24.759370021283804, 274.45496486538042, -1287.8848034377916, 3346.6577012530388, -5076.8777933537731, 4355.7744436071553, -1738.6636645631254, 154.75307248529856};
    double summer=0,theta;
    theta=1-T/reduce.T;
    	
for (int i=1; i<=8; i++)
{
    summer += N[i]*pow(theta,t[i]);
}
return reduce.rho*(summer+1);

}


double R407FClass::rhosatV(double T)
{
    // Maximum absolute error is 0.081139 % between 149.496265 K and 354.803846 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.5, 1.8333333333333333, 2.1666666666666665};
    const double N[] = {0, 181.80471853012997, -2600.8741765435038, 15966.589090960564, -54677.117792695906, 112352.27942151741, -136064.66720260054, 81474.93533760219, -23039.880512878277, 7911.1923524425283, -1516.8853482173406};
    double summer=0,theta;
    theta=1-T/reduce.T;
    	
for (int i=1; i<=10; i++)
{
    summer += N[i]*pow(theta,t[i]);
}
return reduce.rho*exp(reduce.T/T*summer);

}

