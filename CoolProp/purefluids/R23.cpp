#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "R23.h"

R23Class::R23Class()
{
	double n[] = {0.0, 7.041529, -8.259512, 0.008053040, -0.08617615, 0.006333410, -0.1863285, 0.3280510, 0.5191023, 0.06916144, -0.005045875, -0.01744221, -0.05003972, 0.04729813, -0.06164031, 0.01583585, -0.001795790, -0.001099007};
	double t[] = {0, 0.744, 0.94, 4.3, 1.46, 0.68, 4.8, 1.5, 2.07, 0.09, 9.6, 0.19, 11.2, 0.27, 1.6, 10.3, 14.0, 15.0};
	double d[] = {0, 1, 1, 1, 2, 5, 1, 2, 3, 5, 1, 2, 2, 4, 4, 4, 2, 2};
	double c[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 4};

	//Critical parameters
	crit.rho = 7.52*70.01385; //[kg/m^3]
	crit.p = PressureUnit(4832, UNIT_KPA); //[kPa]
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

	EOSReference.assign("Steven G. Penoncello,Eric W. Lemmon,Richard T Jacobsen,Zhengjun Shan, \"A Fundamental Equation for Trifluoromethane (R-23)\", J. Phys. Chem. Ref. Data, Vol. 32, No. 4, 2003");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("R23");
	REFPROPname.assign("R23");

	BibTeXKeys.EOS = "Penoncello-JPCRD-2003";
	BibTeXKeys.ECS_LENNARD_JONES = "Shan-ASHRAE-2000";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.CONDUCTIVITY = "Shan-ASHRAE-2000";
	BibTeXKeys.VISCOSITY = "Shan-ASHRAE-2000";
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
    return reduce.p.Pa*exp(reduce.T/T*summer);
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

double R23Class::viscosity_Trho(double T, double rho)
{
	double C1 = 1.3163, //
		   C2 = 0.1832,
		   DeltaGstar = 771.23,
		   rhoL = 32.174,
		   rhocbar = 7.5114,
		   DELTAeta_max = 3.967,
		   k =	1.380658e-23,
		   N_A = 6.022137e23,
		   pi = 3.141592654,
		   Ru = 8.31451;

	double a[] = {0.4425728, -0.5138403, 0.1547566, -0.02821844, 0.001578286};
	double e_k, sigma;
	this->ECSParams(&e_k,&sigma);
	double Tstar = T/e_k;
	double logTstar = log(Tstar);
	double Omega = exp(a[0]+a[1]*logTstar+a[2]*pow(logTstar,2)+a[3]*pow(logTstar,3)+a[4]*pow(logTstar,4));
	double eta_DG = 5.0/16.0*sqrt(params.molemass*k*T/(1000*pi*N_A))*1e24/(sigma*sigma*Omega); // uPa-s

	double rhobar = rho/params.molemass; // [mol/L]
	double eta_L = C2*(rhoL*rhoL)/(rhoL-rhobar)*sqrt(T)*exp(rhobar/(rhoL-rhobar)*DeltaGstar/(Ru*T));

	double chi = rhobar - rhocbar;
	double tau = T - reduce.T;

	double DELTAeta_c = 4*DELTAeta_max/((exp(chi)+exp(-chi))*(exp(tau)+exp(-tau)));

	return (pow((rhoL-rhobar)/rhoL,C1)*eta_DG+pow(rhobar/rhoL,C1)*eta_L+DELTAeta_c)/1e6;
}

double R23Class::conductivity_Trho(double T, double rho)
{
	double B1 = -2.5370, //
		   B2 = 0.05366,
		   C1 = 0.94215,
		   C2 = 0.14914,
		   DeltaGstar = 2508.58,
		   rhoL = 68.345,
		   rhocbar = 7.5114,
		   DELTAlambda_max = 25,
		   Ru = 8.31451;

	double lambda_DG = B1 + B2*T;

	double rhobar = rho/params.molemass; // [mol/L]
	double lambda_L = C2*(rhoL*rhoL)/(rhoL-rhobar)*sqrt(T)*exp(rhobar/(rhoL-rhobar)*DeltaGstar/(Ru*T));

	double chi = rhobar - rhocbar;
	double tau = T - reduce.T;

	double DELTAlambda_c = 4*DELTAlambda_max/((exp(chi)+exp(-chi))*(exp(tau)+exp(-tau)));

	return (pow((rhoL-rhobar)/rhoL,C1)*lambda_DG+pow(rhobar/rhoL,C1)*lambda_L+DELTAlambda_c)/1e3;
}
