
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
    static double N[]={0, 0.998444, -1.02757, 1.17187, -0.13117, 0.00363706, -2.52408, -1.9824, -0.847253, 0.179562, -0.271224, -0.0785437, 0.670705, -0.636419, 0.0554872, -0.0631532, -0.000410728, 0.00745008, -0.00671065, -0.0172811, 0.0170546, -0.00525303};
    static double t[]={0, 0.425478, 1.20966, 2.9712, 2.94926, 0.203202, 1.91891, 1.78683, 3.00879, 0.192732, 0.74397, 2.9998, 2.10814, 4.29273, 0.250117, 7.00003, 4.69963, 13.0001, 15.9999, 25, 17, 7.39999};
    static double d[]={0, 1, 1, 1, 2, 5, 1, 2, 3, 5, 5, 5, 1, 1, 4, 4, 9, 2, 2, 4, 5, 6};
    static double l[]={0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3};
    
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
    
    phirlist.push_back(new phir_power(N,d,t,l,1,22-1,22));

    phi0list.push_back(new phi0_lead(0,0));
    phi0list.push_back(new phi0_logtau(-1.0));
    phi0list.push_back(new phi0_cp0_poly(a[1],b[1],crit.T,298.15));
    phi0list.push_back(new phi0_Planck_Einstein(a,b,2,4,4+1));

	params.HSReferenceState = "IIR";

    // Limits of EOS
    limits.Tmin = params.Ttriple;

    name.assign("R407F");

}



double R407FClass::psatL(double T)
{
    // Maximum absolute error is 0.320300 % between 149.496265 K and 355.803836 K
    const double t[]={0, 1, 2, 3, 7};
    const double N[]={0, 0.1526407322941514, -7.5051053732508395, 1.2368935727859358, -3.5752506378154121};
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
    // Maximum absolute error is 0.190135 % between 149.496265 K and 355.803836 K
    const double t[]={0, 1, 2, 3, 5, 9};
    const double N[]={0, -0.23308294258564521, -6.9866068352166923, 0.72568262854955323, -1.3377527378607106, -4.94392089320006};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1;i<=5;i++)
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
    // Maximum absolute error is 0.080552 % between 149.496265 K and 354.803846 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.5, 1.8333333333333333, 2.1666666666666665};
    const double N[] = {0, 181.83167536143111, -2601.0011679800691, 15966.232630027936, -54673.285181384716, 112341.17583197542, -136049.36963343105, 81465.618886849174, -23037.701739851385, 7910.7331971494532, -1516.858969623562};
    double summer=0,theta;
    theta=1-T/reduce.T;
    	
for (int i=1; i<=10; i++)
{
    summer += N[i]*pow(theta,t[i]);
}
return reduce.rho*exp(reduce.T/T*summer);

}

