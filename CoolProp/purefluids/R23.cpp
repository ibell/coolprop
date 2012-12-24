#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "R23.h"

static char EOSstr [] = "Steven G. Penoncello,Eric W. Lemmon,Richard T Jacobsen,Zhengjun Shan, \"A Fundamental Equation for Trifluoromethane (R-23)\", J. Phys. Chem. Ref. Data, Vol. 32, No. 4, 2003";

R23Class::R23Class()
{
	double n[] = {0.0, 7.041529, -8.259512, 0.008053040, -0.08617615, 0.006333410, -0.1863285, 0.3280510, 0.5191023, 0.06916144, -0.005045875, -0.01744221, -0.05003972, 0.04729813, -0.06164031, 0.01583585, -0.001795790, -0.001099007};
	double t[] = {0, 0.744, 0.94, 4.3, 1.46, 0.68, 4.8, 1.5, 2.07, 0.09, 9.6, 0.19, 11.2, 0.27, 1.6, 10.3, 14.0, 15.0};
	double d[] = {0, 1, 1, 1, 2, 5, 1, 2, 3, 5, 1, 2, 2, 4, 4, 4, 2, 2};
	double c[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 4};

	//Critical parameters
	crit.rho = 7.52*70.01385; //[kg/m^3]
	crit.p = 4832; //[kPa]
	crit.T = 299.293; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 70.01385;
	params.Ttriple = 118.02;
	params.accentricfactor = 0.262964892496154;
	params.R_u = 8.314472;
	params.ptriple = 0.058;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,17,18));

	const double a1 = 2.999, a2 = -8.31386064, a3 = 6.55087253;
	phi0list.push_back(new phi0_lead(a2,a3));
	phi0list.push_back(new phi0_logtau(a1));

	const double u0[] = {0, 744/crit.T, 1459/crit.T, 2135/crit.T, 4911/crit.T};
	const double v0[] = {0, 2.371, 3.237, 2.610, 0.8274};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("R23");
	REFPROPname.assign("R23");
}

double R23Class::psat(double T)
{
    const double ti[]={0,1.0,1.5,2.4,3.9};
    const double Ni[]={0,-7.2631,1.3140,-0.78507,-3.1991};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double R23Class::rhosatL(double T)
{
    const double ti[]={0,0.37,0.94,3.1};
    const double Ni[]={0,2.2636,0.47007,0.28660};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=3;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*(summer+1);
}
double R23Class::rhosatV(double T)
{
    // Maximum absolute error is 0.872847 % between 277.060001 K and 557.375990 K
    const double ti[]={0,0.43,1.4,3.7,8.0};
    const double Ni[]={0,-3.5136,-7.7491,-24.871,-65.637};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}