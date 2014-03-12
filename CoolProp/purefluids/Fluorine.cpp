#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Fluorine.h"

FluorineClass::FluorineClass()
{
    double n[] = {0, 1.51144749736E+00, -2.98666288409E+00, 3.29644905098E+00, -2.98458624201E+00, -2.28688966459E+00, -1.09492193400E+00, 3.04775277572E+00, 1.15689564208E-01, -1.16100171627E+00, 2.95656394476E-01, 7.11482542928E-02, -1.71363832155E-03, 6.65317955515E-04, 5.06026676251E+00, -6.29268435440E+00, 6.17784808739E+00, -1.55366191788E+00, -2.87170687343E+00, 3.17214480494E+00, -2.67969025215E+00, 2.71865479252E+00, -1.07191065039E+00, 1.26597342291E+00, -7.06244695489E-01, 2.68707888826E-01, 5.27251190274E-02, 5.44411481926E-02, 2.28949994105E-04, -5.47908264304E-10, -9.64273224950E-02, 3.68084486225E-04};
	double d[] = {0, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 8, 9, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 6, 7, 8, 12, 4, 6, 6};
	double t[] = {0, 0, 0.5, 1.5, 2, 0.5, 1, 0.5, 2, 0.5, 1, 0, 0.5, 0, 1, 3, 4, 5, 1, 4, 5, 1, 3, 5, 4, 4, 1, 1, 5, 30, 20, 25};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
	double g[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.078102576, 1.078102576, 1.078102576, 1.078102576, 1.078102576, 1.078102576, 1.078102576, 1.078102576, 1.078102576, 1.078102576, 1.078102576, 1.078102576, 1.078102576, 1.078102576, 1.078102576, 2.156205153, 3.234307729, 3.234307729};

    // Critical parameters
    crit.rho = 15.603*37.99681;
    crit.p = PressureUnit(5172.4, UNIT_KPA);
    crit.T = 144.414;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 37.99681;
    params.Ttriple = 53.4811;
	params.ptriple = 0.23881;
    params.accentricfactor = 0.0449;
    params.R_u = 8.31448;

    // Limits of EOS
	limits.Tmin = params.Ttriple;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

	// Residual part
    phirlist.push_back(new phir_power(n,d,t,c,1,13,32));
	phirlist.push_back(new phir_exponential(n,d,t,c,g,14,31,32));

	// Ideal-gas part
	double f0[]={0.0, 3.0717001e-6, -5.2985762e-5, -16.372517, 3.6884682e-5, 2.5011231, 1.0127670, 8.9057501, 4.3887271};
	double g0[] = {0, -4, -3, 1, 2, 0, 0, 0, 0};
	phi0list.push_back(new phi0_power(f0, g0,1,4,9));
	phi0list.push_back(new phi0_logtau(f0[5]));
	phi0list.push_back(new phi0_Planck_Einstein2(f0[6],f0[7],-1));
	phi0list.push_back(new phi0_power(f0[8], 0));
	phi0list.push_back(new phi0_lead(0,0)); // This terms is needed, but not listed in Fluorine book

    name.assign("Fluorine");
    aliases.push_back(std::string("fluorine"));
    aliases.push_back(std::string("FLUORINE"));
    REFPROPname.assign("FLUORINE");

	BibTeXKeys.EOS = "deReuck-BOOK-1990";
	BibTeXKeys.CP0 = "deReuck-BOOK-1990";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double FluorineClass::psat(double T)
{
    const double t[]={0, 1, 1.5, 2, 3.5, 4.5};
    const double N[]={0, -6.1843496, 1.2095388, -6.1708665, -0.30963258, 0.068976007};
    double summer=0,theta;
    theta=reduce.T/T-1;
    for (int i=1;i<=5;i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(T/reduce.T*summer);
}

double FluorineClass::rhosatL(double T)
{
    // Maximum absolute error is 0.052236 % between 53.481101 K and 144.413990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.6666666666666667, 2.0, 2.3333333333333335};
    const double N[] = {0, 1.2896950515325272, -36.259697697303622, 389.20444135863613, -2176.3325373261714, 7366.9110799263972, -15496.555526995258, 19515.494134734785, -12219.632837082365, 3756.9802418741001, -1364.3462914968964, 265.8645475479338};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=11; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double FluorineClass::rhosatV(double T)
{
    // Maximum absolute error is 0.251225 % between 53.481101 K and 144.413990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.8333333333333333};
    const double N[] = {0, 5.9616454798694276, -183.74978318636624, 1940.0783642835002, -10472.357185702867, 33021.901889876041, -64463.557858790977, 78264.428052775431, -55910.789427023621, 19118.90631766606, -1327.7852195930525};
    double summer=0,theta;
    theta=1-T/reduce.T;    	
	for (int i=1; i<=10; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}
double FluorineClass::surface_tension_T(double T)
{
	// Mulero, JPCRD 2012
	return 0.03978*pow(1-T/reduce.T,1.218);
}
