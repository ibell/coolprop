
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
    static double a[]={0, 1.440000, 2.000227, 5.359371, 3.496132};
    static double b[]={0, 0.256551293851, 696.9434/355.804, 1723.0916/355.804, 3874.9951/355.804};

    // Constants for the residual contribution
    static double N[]={0.5, 0.919875, -1.82678, -0.352996, 0.0588362, 0.000129927, 0.259443, 0.708701, 0.0179878, -0.305208, -0.0510867, -0.0910294, -0.0300902};
    static double t[]={0, 0.25, 1.25, 1.5, 0.25, 0.875, 2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5};
    static double d[]={0, 1, 1, 1, 3, 7, 1, 2, 5, 1, 1, 4, 2};
    static double l[]={0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3};
    
    // Other fluid parameters
    params.molemass = 82.0583; //[kg/kmol]
    params.Ttriple = 200; //[K]
    params.accentricfactor = 0.360036864771; //[-]
    params.R_u = 8.314472;
    isPure = false;
    
    // Critical parameters
    crit.rho = 477.285;
    crit.p = PressureUnit(4754.61,UNIT_KPA);
    crit.T = 355.804;
    crit.v = 1.0/crit.rho;
    
    phirlist.push_back(new phir_power(N,d,t,l,1,13-1,13));

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
    // Maximum absolute error is 0.238386 % between 149.496265 K and 355.803836 K
    const double t[]={0, 1, 2, 3, 5, 9};
    const double N[]={0, 0.18082710436189095, -7.7983091066218071, 2.0533745861071391, -2.1140573561077987, -2.3513695979484255};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1;i<=5;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}


double R407FClass::psatV(double T)
{
    // Maximum absolute error is 0.465542 % between 149.496265 K and 355.803836 K
    const double t[]={0, 1, 2, 4, 8};
    const double N[]={0, -0.25845953186204634, -6.7270062079395379, 0.13564366673540645, -5.5138739927573255};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1;i<=4;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}


double R407FClass::rhosatL(double T)
{
    // Maximum absolute error is 0.110596 % between 149.496265 K and 354.803846 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.5};
    const double N[] = {0, -25.016336368491949, 275.16509661941035, -1282.9862444814557, 3317.3858644879506, -5014.0917140512256, 4290.9626091919372, -1709.960437659249, 151.99307236915118};
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
    // Maximum absolute error is 0.079941 % between 149.496265 K and 354.803846 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.5, 1.8333333333333333, 2.1666666666666665};
    const double N[] = {0, 181.82633126699102, -2600.9978883137055, 15966.622688274818, -54675.956461439942, 112349.32219323577, -136062.35313491867, 81475.183272830909, -23041.320486947832, 7912.2244219331069, -1517.1758761518065};
    double summer=0,theta;
    theta=1-T/reduce.T;
    	
for (int i=1; i<=10; i++)
{
    summer += N[i]*pow(theta,t[i]);
}
return reduce.rho*exp(reduce.T/T*summer);

}

