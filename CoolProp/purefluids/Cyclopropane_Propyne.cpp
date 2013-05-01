#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Cyclopropane_Propyne.h"

CycloPropaneClass::CycloPropaneClass()
{
	double n[] = {0, -1.15634007067133000, 2.52574627508968000, -2.82266128097357000, 0.28357680160523500, -0.08427204963322530, 0.93108856598845400, -1.05296839887510000, 0.43202158160276800, -0.25110886434063600, 0.12772589248249800, 0.04836223357857030, -0.01164740783339940, 0.00033400656553506, -1.37016097592368000, 2.12444673007915000, -0.57890894273166200, 0.30494577050653800, -0.18427616517007700, -0.29211146040448100, 1.37016097592368000, -2.12444673007915000, 0.57890894273166200};
	double t[] = {0, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5};
	double d[] = {0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 0, 0, 0, 2, 2, 2, 0, 0, 0};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2};

	//Critical parameters
	crit.rho = 258.5; //[kg/m^3]
	crit.p = 5579.7; //[kPa]
	crit.T = 398.3; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 42.081;
	params.Ttriple = 145.7;
	params.accentricfactor = 0.13054956636875170;
	params.R_u = 8.3143;
	params.ptriple = 342.70692277003945;

	// Limits of EOS
	limits.Tmin = 273;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power(n,d,t,c,1,22,23));

	phi0list.push_back(new phi0_lead(0,0));
	phi0list.push_back(new phi0_logtau(-1));

	const double a0[] = {0, -9.05307e-3/0.197578, 5.05504e-5/0.197578, -7.72237e-8/0.197578, 4.0538e-11/0.197578};
	const double n0[] = {0, 1, 2, 3, 4};
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_constant(1.26016e0/0.197578,crit.T,298.15));
	phi0list.push_back(new phi0_cp0_poly(a0_v,n0_v,crit.T,298.15,1,4));

	name.assign("CycloPropane");
	aliases.push_back(std::string("cyclopropane"));
	aliases.push_back(std::string("Cyclopropane"));
	REFPROPname.assign("CYCLOPRO");

	BibTeXKeys.EOS = "Polt-CT-1992";
	BibTeXKeys.CP0 = "Polt-CT-1992";
}
double CycloPropaneClass::psat(double T)
{
    // Maximum absolute error is 0.187882 % between 273.000001 K and 398.299990 K
    const double t[]={0, 1, 2, 3, 7};
    const double N[]={0, -0.029915737385806201, -6.3990202086270314, 0.91387815294930119, -2.4469438384812348};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=3;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
	double r = reduce.p*exp(reduce.T/T*summer);
    return reduce.p*exp(reduce.T/T*summer);
}

double CycloPropaneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.127947 % between 273.000001 K and 398.299990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.1666666666666667, 1.3333333333333333, 1.6666666666666667};
    const double N[] = {0, 4.2130826166202509, -37.918092707580193, 149.7380622313288, -281.61370700759898, 244.50899198420026, -166.77358044074828, 101.52186468251415, -10.969640624299528};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
	for (i=1; i<=8; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	double r = reduce.rho*(summer+1);
	return reduce.rho*(summer+1);
}

double CycloPropaneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.171588 % between 273.000001 K and 398.299990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.8333333333333333};
    const double N[] = {0, -4.0162733322058601, 17.579157820549558, 129.28293085434478, -1511.7458291465116, 6396.6245746851037, -15200.674472181481, 21775.796916494346, -18139.92511444999, 7189.2193923624845, -661.11090761629998};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
	for (i=1; i<=10; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

//PropyneClass::PropyneClass()
//{
//	double n[] = {0, -1.15634007067133000, 2.52574627508968000, -2.82266128097357000, 0.28357680160523500, -0.08427204963322530, 0.93108856598845400, -1.05296839887510000, 0.43202158160276800, -0.25110886434063600, 0.12772589248249800, 0.04836223357857030, -0.01164740783339940, 0.00033400656553506, -1.37016097592368000, 2.12444673007915000, -0.57890894273166200, 0.30494577050653800, -0.18427616517007700, -0.29211146040448100, 1.37016097592368000, -2.12444673007915000, 0.57890894273166200};
//	double t[] = {0, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5};
//	double d[] = {0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 0, 0, 0, 2, 2, 2, 0, 0, 0};
//	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2};
//
//	//Critical parameters
//	crit.rho = 244.9; //[kg/m^3]
//	crit.p = 5626; //[kPa]
//	crit.T = 402.38; //[K]
//	crit.v = 1/crit.rho; 
//
//	// Other fluid parameters
//	params.molemass = 40.06;
//	params.Ttriple = _HUGE;
//	params.accentricfactor = _HUGE;
//	params.R_u = 8.3143;
//	params.ptriple = _HUGE;
//
//	// Limits of EOS
//	limits.Tmin = params.Ttriple;
//	limits.Tmax = 500.0;
//	limits.pmax = 100000.0;
//	limits.rhomax = 1000000.0*params.molemass;
//
//	phirlist.push_back(new phir_power(n,d,t,c,1,22,23));
//
//	const double n5 = 4.9916462, n6 = -0.1709449, n0 = 9.28421;
//	phi0list.push_back(new phi0_lead(n5,n6));
//	phi0list.push_back(new phi0_logtau(n0-1));
//
//	const double u0[] = {0, 21/crit.T, 1340/crit.T, 1672/crit.T, 7395/crit.T};
//	const double v0[] = {0, 1.48525, 0.822585, 16.2453, 1.15925};
//	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
//	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));
//
//	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));
//
//	EOSReference.assign(EOSstr);
//	TransportReference.assign("Using ECS in fully predictive mode");
//
//	name.assign("CycloPropane");
//	aliases.push_back(std::string("cyclopropane"));
//	aliases.push_back(std::string("Cyclopropane"));
//	REFPROPname.assign("CYCLOPRO");
//
//	BibTeXKeys.EOS = "";
//}
//double PropyneClass::psat(double T)
//{
//    // Maximum absolute error is 0.058428 % between 273.000001 K and 402.379990 K
//    const double t[]={0, 1, 2, 3, 7};
//    const double N[]={0, -0.0064125853696993685, -6.8243990527455569, 0.7209626970691978, -0.51952084490286898};
//    double summer=0,theta;
//    int i;
//    theta=1-T/reduce.T;
//    for (i=1;i<=3;i++)
//    {
//        summer += N[i]*pow(theta,t[i]/2);
//    }
//    return reduce.p*exp(reduce.T/T*summer);
//}
//
//double PropyneClass::rhosatL(double T)
//{
//    // Maximum absolute error is 0.126203 % between 273.000001 K and 402.379990 K
//    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.8333333333333333, 2.1666666666666665};
//    const double N[] = {0, 4.0434726582869249, -31.189422913936159, 68.555265694482571, 209.65666040673665, -1446.9272712891047, 3427.0000532850067, -4065.0391509043716, 2098.6288299939783, -400.20063454356421, 141.93833949828849};
//    double summer=0,theta;
//    int i;
//    theta=1-T/reduce.T;
//	for (i=1; i<=10; i++)
//	{
//		summer += N[i]*pow(theta,t[i]);
//	}
//	return reduce.rho*(summer+1);
//}
//
//double PropyneClass::rhosatV(double T)
//{
//    // Maximum absolute error is 0.058384 % between 273.000001 K and 402.379990 K
//    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.8333333333333333};
//    const double N[] = {0, -3.247356661315155, 21.032480035171591, -499.92815044636194, 2331.3544924983121, -5555.0311630324813, 7897.7172502284338, -6595.1414649864946, 2653.8699614585994, -257.78270300881718};
//    double summer=0,theta;
//    int i;
//    theta=1-T/reduce.T;
//	for (i=1; i<=9; i++)
//	{
//		summer += N[i]*pow(theta,t[i]);
//	}
//	return reduce.rho*exp(reduce.T/T*summer);
//
//}