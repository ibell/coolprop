#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Undecane.h"

UndecaneClass::UndecaneClass()
{
	double n[] = {0.0, -0.66172706, 1.3375396, -2.5608399, 0.10678910, 0.28873614e-3, 0.49587209e-1, 0.55407101e-7, 0.99754712, 1.5774025, 0.13108354e-2, -0.59326961, -0.93001876e-1, -0.17960228, -0.22560853e-1};
	double t[] = {0, 1.5,0.25,1.25,0.25,0.875,1.375,0,2.375,2.0,2.125,3.5,6.5,4.75,12.5};
	double d[] = {0, 1, 1, 1, 3, 7, 2, 1, 1, 2, 5, 1, 1, 4, 2};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3};

	//Critical parameters
	crit.rho = 1.5149*156.31; //[kg/m^3]
	crit.p = 1990.4; //[kPa]
	crit.T = 638.8; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 156.30826;
	params.Ttriple = 247.541;
	params.accentricfactor = 0.53903710137185668;
	params.R_u = 8.314472;
	params.ptriple = 0.00044605108038132317;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,14,15));

	phi0list.push_back(new phi0_lead(0, 0));
	phi0list.push_back(new phi0_logtau(-1));

	double alpha[] = {-1158848, 20321.8, -119.4274, 0.4284215, -4.157728e-4, 1.61828e-7};
	double beta[] = {-2, -1, 0, 1, 2, 3};

	std::vector<double> alpha_v(alpha, alpha+sizeof(alpha)/sizeof(double));
	std::vector<double> beta_v(beta, beta+sizeof(beta)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(alpha_v, beta_v, crit.T, 298.15, 0, 5));

	EOSReference.assign("Aleksandrov, I. S. and A. A. Gerasimov and B. A. Grigor’ev \"Using Fundamental Equations of State for Calculating the Thermodynamic Properties of Normal Undecane\" Thermal Engineering, 2011, Vol. 58, No. 8, pp. 691–698");
	TransportReference.assign("Using ECS in fully predictive mode.");

	name.assign("n-Undecane");
	aliases.push_back("Undecane");
	aliases.push_back("C11");
	REFPROPname.assign("C11");
	
	ECSReferenceFluid = "n-Dodecane";

	BibTeXKeys.EOS = "Aleksandrov-TE-2011";
}
double UndecaneClass::psat(double T)
{
    // Maximum absolute error is 0.027305 % between 247.541001 K and 638.799990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3,9};
    const double Ni[]={0,-9.3813452614877857, 4.2825084033495191, -5.9314504005650406, 1.6773321093108704, -11.636954318879175, 16.964446346543291, -14.058870766446166 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double UndecaneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.304330 % between 247.541001 K and 638.799990 K
    const double ti[]={0,0.63229621538870906, 1.0038056284523527, 0.86235851330609503, 0.87952338233778615, 0.90457520801426383, 8.8303918023981769};
    const double Ni[]={0,67.587954960804993, 622.19343233092457, -15756.38302838042, 28308.555888123774, -13240.532777759823, 0.15299119998547792};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double UndecaneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.198848 % between 247.541001 K and 638.799990 K
    const double ti[]={0,0.45209637851715323, 0.62014004645697651, 3.2948486604748419, 3.391538201868411, 3.5175920402216905, 3.4656311291203314};
    const double Ni[]={0,-3.2503273294249282, -2.1655950913673307, -3691.6518899679063, 14736.668843150595, 8908.4150123462787, -19962.521747567687};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}