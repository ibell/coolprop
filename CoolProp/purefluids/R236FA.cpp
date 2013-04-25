#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "R236FA.h"

R236FAClass::R236FAClass()
{
	double n[] = {0.0, 0.04453255, 1.777017, -2.230519, -0.6708606, 0.1587907, -1.425119, -0.6461628, 0.8469985, -0.5635356, -0.01535611, 1.156362, -0.4070310, -0.2172753, -1.007176, -0.00006902909};
	double t[] = {0, 1.07, 0.222, 0.66, 1.33, 0.227, 2.33, 1.94, 1.53, 2.65, 0.722, 1.11, 2.31, 3.68, 4.23, 0.614};
	double d[] = {0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3, 3, 2};
	double c[] = {0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.02, 1.336, 1.055, 5.84, 16.2};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.42, 2.31, 0.89, 80.0, 108.0};
	double gamma[] = {0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.13, 0.67, 0.46, 1.28, 1.2};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.712, 0.91, 0.677, 0.718, 1.64};

	//Critical parameters
	crit.rho = 3.626*152.0384; //[kg/m^3]
	crit.p = 3200; //[kPa]
	crit.T = 398.07; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 152.0384;
	params.Ttriple = 179.6;
	params.accentricfactor = 0.37691179763369909;
	params.R_u = 8.314472;
	params.ptriple = 0.16032679731173749;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,10,16));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 11,15,16));

	const double n5 = -17.5983849, n6 = 8.87150449, n0 = 10.175;
	phi0list.push_back(new phi0_lead(n5,n6));
	phi0list.push_back(new phi0_logtau(n0-1));

	const double u0[] = {0, 962/crit.T, 2394/crit.T, 5188/crit.T};
	const double v0[] = {0, 9.8782, 18.236, 49.934};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,3));

	EOSReference.assign("Jiang Pan and Xinfang Rui and Xiaodong Zhao and Liming Qiu, \"An equation of state for the thermodynamic properties of 1,1,1,3,3,3-hexafluoropropane (HFC-236fa)\", Fluid Phase Equilibria 321 (2012) 10--16\n\nNote: REFPROP 9.0 does not use this EOS, they use thermodynamic corresponding states");
	TransportReference.assign("Using ECS in fully predictive mode.");

	name.assign("R236FA");
	aliases.push_back(std::string("R236fa"));
	REFPROPname.assign("R236FA");

	BibTeXKeys.EOS = "Pan-FPE-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double R236FAClass::psat(double T)
{
    // Maximum absolute error is 0.030234 % between 179.520001 K and 398.069990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3,9};
    const double Ni[]={0,-7.8508917792245674, 1.4567640808160465, -0.062102900606892143, -9.3819597000175285, 14.550729344612401, -36.422971826769853, 37.644436183238703 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=7;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double R236FAClass::rhosatL(double T)
{
    // Maximum absolute error is 0.043993 % between 179.520001 K and 398.069990 K
    const double ti[]={0,0.67593967200386329, 1.0388728256250355, 1.0958022963770171, 1.2075666821666764, 1.4489255803047685, 1.534595334230606, 1.8422228222992407};
    const double Ni[]={0,33.66527737771591, -3018.4889875066406, 5655.5413015392987, -3680.3090670411439, 2448.5303944282114, -1498.7494751321274, 61.236298277462517};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=7;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double R236FAClass::rhosatV(double T)
{
    // Maximum absolute error is 0.099088 % between 179.520001 K and 398.069990 K
    const double ti[]={0,0.50195291798025177, 1.3220357032983603, 2.5731594638194166, 2.4474465250730582, 2.4731206226774058, 3.0334887073286887, 6.1686391096790096};
    const double Ni[]={0,-4.3707771378161109, -7.2992995710239832, 6024.9834521647063, 17995.745586913225, -23734.497572478816, -289.70356229231936, 4.6891629794522327};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=7;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}