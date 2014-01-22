#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "R161_Fluoroethane.h"

R161Class::R161Class()
{
	double n[] = {0.0, 1.511, -2.3, -0.457, 0.1683, 0.04133, 0.62187, -0.0265, -1.03, -0.285, -0.476, 0.82, -0.3532, -0.116, -0.0220583, -1.63148};
	double t[] = {0, 0.37, 0.97, 1.14, 0.744, 1, 1.26, 1, 1.8, 3, 2.25, 1, 1.2, 5.3, 1, 4};
	double d[] = {0, 1, 1, 2, 3, 4, 2, 7, 1, 2, 3, 1, 1, 3, 3, 3};
	double c[] = {0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 0, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.96, 1.35, 1.26, 1.23, 16.8};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.7, 5.2, 3.9, 4.7, 413};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.9, 0.69, 0.67, 0.67, 1.15};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.683, 0.892, 0.785, 1.33, 0.86};

	//Critical parameters
	crit.rho = 6.28*48.0595; //[kg/m^3]
	crit.p = PressureUnit(5010, UNIT_KPA); //[kPa]
	crit.T = 375.25; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 48.0595;
	params.Ttriple = 130;
	params.accentricfactor = 0.21624284106618674;
	params.R_u = 8.314472;
	params.ptriple = 0.005512;

	// Limits of EOS
	limits.Tmin = 130;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,10,16));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 11,15,16));

	const double a0 = -6.9187, a1 = 5.4788, n0 = 4;
	phi0list.push_back(new phi0_lead(a0, a1));
	phi0list.push_back(new phi0_logtau(n0-1));

	const double u0[] = {0, 420/crit.T, 1548/crit.T, 3882/crit.T};
	const double v0[] = {0, 2.059, 9.253, 6.088};

	phi0list.push_back(new phi0_Planck_Einstein(v0,u0,1,3,4));

	EOSReference.assign("Jiangtao Wu and Yong Zhou, \"An Equation of State for Fluoroethane (R161)\", Int J Thermophys (2012) 33:220-234");
	TransportReference.assign("Using ECS in fully predictive mode.");

	name.assign("R161");
	aliases.push_back(std::string("Fluoroethane"));
	aliases.push_back(std::string("FLUOROETHANE"));

	REFPROPname.assign("R161");

	BibTeXKeys.EOS = "Wu-IJT-2012";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R161Class::psat(double T)
{
	// Max error is  0.230661857744 % between 130.0 and 375.249999 K
    const double t[]={0, 0.3575, 0.3605, 0.8333333333333334, 1.1666666666666667, 3.6666666666666665, 5.5};
    const double N[]={0, 0.31348208716072901, -0.23919030666253932, -3.2448016166944127, -3.2492132423845081, -1.5585118618296889, -2.427095112126342};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double R161Class::rhosatL(double T)
{
    // Maximum absolute error is 0.240076 % between 130.000000 K and 375.250000 K
    const double t[] = {0, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.6666666666666667, 1.8333333333333333, 2.0};
    const double N[] = {0, -0.047548020147135293, 72.327097393829419, -819.54669834482479, 4946.726597164351, -18870.150565849312, 48463.342332873617, -84949.003918399481, 99927.349954656311, -75150.383145190368, 32522.693777539549, -6140.6005518003303};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=11; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double R161Class::rhosatV(double T)
{
    // Maximum absolute error is 0.294100 % between 130.000000 K and 375.250000 K
    const double t[] = {0, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.6666666666666667, 1.8333333333333333, 2.0, 2.1666666666666665};
    const double N[] = {0, 4.4900108859020449, -289.55507432056027, 4058.7708236096973, -29986.440551303382, 137204.7283140849, -417289.03201597766, 870535.22359268588, -1255033.5731927676, 1232137.9562410619, -787892.80993800913, 296355.72520564846, -49815.743244828758};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=12; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}
