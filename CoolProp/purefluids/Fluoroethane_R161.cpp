/* 
Properties of Fluoroethane (R161)
by Ian Bell
*/

#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include "math.h"
#include "stdio.h"
#include <string.h>
#include "CoolProp.h"
#include "FluidClass.h"
#include "Fluoroethane_R161.h"

FluoroethaneClass::FluoroethaneClass()
{
	// From Wu, IJT, 2012, 33:220-234

	double n[] = {0, 1.511, -2.3, -0.457, 0.1683, 0.04133, 0.62187, -0.0265, -1.03, -0.285, -0.476, 0.82, -0.3532, -0.116, -0.0220583, -1.63148};
	double d[] = {0, 1, 1, 2, 3, 4, 2, 7, 1, 2, 3, 1, 1, 3, 3, 3};
	double t[] = {0, 0.37, 0.97, 1.14, 0.744, 1, 1.26, 1, 1.8, 3, 2.25, 1, 1.2, 5.3, 1, 4};
	double l[] = {0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 0, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.96, 1.35, 1.26, 1.23, 16.8};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.7, 5.2, 3.9, 4.7, 413};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.9, 0.69, 0.67, 0.67, 1.15};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.683, 0.892, 0.785, 1.33, 0.86};

	phirlist.push_back(new phir_power(n,d,t,l,1,10,16));
	phirlist.push_back(new phir_gaussian(n,d,t,eta,epsilon,beta,gamma,11,15,16));

	phi0list.push_back(new phi0_lead(-6.9187,5.4788));
	phi0list.push_back(new phi0_logtau(3.0));
	
	// Critical parameters
	crit.rho = 6.28*48.0595;
	crit.p = PressureUnit(5010, UNIT_KPA);
	crit.T = 375.25;
	crit.v = 1.0/crit.rho;

	double u0[] = {0, 420/crit.T, 1548/crit.T, 3882/crit.T};
	double v0[] = {0, 2.059, 9.253, 6.088};
	phi0list.push_back(new phi0_Planck_Einstein(u0,v0,1,3,4));

	// Other fluid parameters
	params.molemass = 48.0595;
	params.Ttriple = 130;
	params.ptriple = 0.005512;
	params.accentricfactor = 0.21624284106618674;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = 130;
	limits.Tmax = 625.0;
	limits.pmax = 1000000.0;
	limits.rhomax = 20.6*params.molemass;

	name.assign("Fluoroethane");
	aliases.push_back("FLUOROETHANE");
	aliases.push_back("R161");
	REFPROPname.assign("R161");

	BibTeXKeys.EOS = "Wu-IJT-2012";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
//void FluoroethaneClass::ECSParams(double *e_k, double *sigma)
//{
//	*e_k = 263.88;  *sigma = 0.49748;
//}
double FluoroethaneClass::psat(double T)
{
	// Max error is  0.157883694009 % between 130.0 and 375.2499 K
	const double ti[]={0, 0.3525, 0.3685, 1.0, 2.5, 4.5, 11.666666666666666};
    const double Ni[]={0, 2.0009345767815701, -2.2094959866666275, -6.0595159878621461, -0.39826055495168905, -3.1969185812765044, -2.5214587112499913};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double FluoroethaneClass::rhosatL(double T)
{
	// Max error is  0.0929136336678 % between 130.0 and 375.2499 K
	const double ti[]={0, 0.376, 0.379, 0.383, 2.0, 3.1666666666666665, 13.5};
    const double Ni[]={0, 930.03044793228037, -1606.6478530343782, 678.07550448593327, -0.69545113594069996, 0.71935193663784658, -2.254513549911378};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double FluoroethaneClass::rhosatV(double T)
{
	// Max error is 0.0797344805483 % between 130.0 and 375.2499 K
	const double ti[]={0, 0.38149999999999995, 0.384, 0.38599999999999995, 2.0, 5.5, 7.0};
    const double Ni[]={0, -10358.6843563317, 23607.248931465219, -13253.107690006933, -1.6561138988731885, -6.748958860389263, 2.8895346478553039};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
	return reduce.rho*exp(reduce.T/T*summer);
}