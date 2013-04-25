#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "R236EA.h"

R236EAClass::R236EAClass()
{
	double n[] = {0.0, 0.051074, 2.5584, -2.9180, -0.71485, 0.15534, -1.5894, -0.784, 0.85767, -0.67235, -0.017953, 1.3165, -0.42023, -0.28053, -1.4134, -0.0000062617};
	double t[] = {0, 1.0, 0.264, 0.5638, 1.306, 0.2062, 2.207, 2.283, 1.373, 2.33, 0.6376, 1.08, 1.67, 3.502, 4.357, 0.6945};
	double d[] = {0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3, 3, 2};
	double c[] = {0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.019, 1.341, 1.034, 5.264, 24.44};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.30, 2.479, 1.068, 79.850, 49.06};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.13, 0.6691, 0.465, 1.280, 0.8781};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.7119, 0.9102, 0.678, 0.7091, 1.727};

	//Critical parameters
	crit.rho = 3.716*152.0384; //[kg/m^3]
	crit.p = 3420; //[kPa]
	crit.T = 412.44; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 152.0384;
	params.Ttriple = 170;
	params.accentricfactor = 0.36878238039914035;
	params.R_u = 8.314472;
	params.ptriple = 17.525807103151166; // At Tmin of 243

	// Limits of EOS
	limits.Tmin = 243;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,10,16));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 11,15,16));

	const double n5 = -14.1214241350, n6 = 10.2355589225, n0 = 3.762;
	phi0list.push_back(new phi0_lead(n5,n6));
	phi0list.push_back(new phi0_logtau(n0-1));

	const double u0[] = {0, 144/crit.T, 385/crit.T, 1536/crit.T, 7121/crit.T};
	const double v0[] = {0, 0.7762, 10.41, 12.18, 3.332};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	EOSReference.assign("Xinfang Rui and Jiang Pan and Yugang Wang, \"An equation of state for the thermodynamic properties of 1,1,1,2,3,3-hexafluoropropane (R236ea)\", Fluid Phase Equilibria 341 (2013) 78-85\n\nNote: Erratum in paper: a1 should be -17.5983849 and a2 should be 8.87150449\n\nNote: REFPROP 9.0 does not use this EOS, they use thermodynamic corresponding states");
	TransportReference.assign("Using ECS in fully predictive mode.");

	name.assign("R236EA");
	aliases.push_back(std::string("R236ea"));
	REFPROPname.assign("R236EA");

	BibTeXKeys.EOS = "Rui-FPE-2013";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double R236EAClass::psat(double T)
{
    // Maximum absolute error is 3.866754 % between 242.000001 K and 412.439990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3,9};
    const double Ni[]={0,-12.435938949131765, 26.900144082811941, -88.321544596554602, 350.28330260511353, -1407.9376828993531, 5978.240328805251, -9529.2569566570437 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=7;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double R236EAClass::rhosatL(double T)
{
    // Maximum absolute error is 0.098519 % between 242.000001 K and 412.439990 K
    const double ti[]={0,-0.00047271337978119942, 1.0642837740271649, 1.9646287619396041, 2.099374228268986, 2.0213001573602218, 2.3007393125389775, 2.5046751761629613};
    const double Ni[]={0,0.35736193894626606, 14.98152652171032, -18598.008235627218, -27903.014272201213, 41667.430800419897, 5849.3107245345163, -1029.7512656495517};
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
double R236EAClass::rhosatV(double T)
{
    // Maximum absolute error is 0.059057 % between 242.000001 K and 412.439990 K
    const double ti[]={0,-0.00027253117723408836, 1.025130694967147, 2.4773331331060562, 2.1703827557176751, 1.6054739802753848, 7.621747384559848, 2.1631918578975449};
    const double Ni[]={0,-0.50154042927614517, -21.871785792621363, -220.5698284625384, 13839.09732823177, 99.188123538279541, -7.1564396559398347, -13706.010193186079};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=7;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}