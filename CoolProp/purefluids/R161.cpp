#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "R161.h"

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
	crit.p = 5010; //[kPa]
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
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,3));

	EOSReference.assign("Jiangtao Wu and Yong Zhou, \"An Equation of State for Fluoroethane (R161)\", Int J Thermophys (2012) 33:220–234");
	TransportReference.assign("Using ECS in fully predictive mode.");

	name.assign("R161");
	aliases.push_back(std::string("R161"));
	REFPROPname.assign("R161");

	BibTeXKeys.EOS = "Wu-IJT-2012";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R161Class::psat(double T)
{
    // Maximum absolute error is 0.070348 % between 130.000001 K and 375.299990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3,9};
    const double Ni[]={0,-7.5071365267951027, 2.8348327441979313, -2.8006268767049143, -0.41450590365286594, -0.88760138670279365, -3.2080318195360324, 2.7252970562347731 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=7;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double R161Class::rhosatL(double T)
{
    // Maximum absolute error is 0.664046 % between 130.000001 K and 375.299990 K
    const double ti[]={0,1.8100670231237792, 0.81784934496037065, 1.223141328833613, 1.949747072104139, 2.4490781003073709, 2.6934402021586124, 3.5463674040730853, 4.3189086783037105, 3.3175110382497359};
    const double Ni[]={0,11261.119384480522, 59.42203076754744, -468.8719880337257, -17445.615620624631, 22527.755455736948, -21001.562358183466, -6845.3992442215504, 301.91969935649598, 11612.799312825631};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=9;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double R161Class::rhosatV(double T)
{
    // Maximum absolute error is 0.551208 % between 130.000001 K and 375.299990 K
    const double ti[]={0,0.59815212632530579, 1.3112754501339567, 1.8174376922182702, 1.8980508024203997, 6.1976639514744205, 2.319522266066437};
    const double Ni[]={0,-7.9863505565397892, 37.174205426775309, -477.94789008409452, 499.45890087388472, -1.6853266358753007, -58.899195162678794};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}