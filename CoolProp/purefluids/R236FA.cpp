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
	double c[] = {0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.02, 1.336, 1.055, 5.84, 16.2};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.42, 2.31, 0.89, 80.0, 108.0};
	double gamma[] = {0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.13, 0.67, 0.46, 1.28, 1.2};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.712, 0.91, 0.677, 0.718, 1.64};

	//Critical parameters
	crit.rho = 3.626*152.0384; //[kg/m^3]
	crit.p = PressureUnit(3200, UNIT_KPA); //[kPa]
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
	// Max error is  0.17009505573 % between 179.6 and 398.069999 K
    const double t[]={0, 0.12200000000000001, 1.1666666666666667, 1.3333333333333333, 2.3333333333333335, 2.5, 3.5};
    const double N[]={0, -0.012979211638964282, -28.829423111516142, 28.838546594062379, -75.471830356030637, 78.305145880492105, -15.984987308278052};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double R236FAClass::rhosatL(double T)
{
    // Maximum absolute error is 0.239080 % between 179.600000 K and 398.070000 K
    const double t[] = {0, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.6666666666666667};
    const double N[] = {0, 110.99431711060721, -653.8278250306164, 1808.4970386386176, -2819.9588105709872, 2549.3846572841617, -1252.0268343097021, 259.90253247283982};
    double summer=0,theta;
    theta=1-T/reduce.T; 	
	for (int i=1; i<=7; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double R236FAClass::rhosatV(double T)
{
    // Maximum absolute error is 0.372237 % between 179.600000 K and 398.070000 K
    const double t[] = {0, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.6666666666666667, 1.8333333333333333};
    const double N[] = {0, -11.794124580366443, -645.63181864379703, 5244.661122543529, -18060.343097280878, 33869.097416810924, -36245.578523710719, 20887.554530748635, -5050.4457446082124};
    double summer=0,theta;
    theta=1-T/reduce.T; 	
	for (int i=1; i<=8; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);

}
