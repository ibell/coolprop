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
	crit.rho = 42.081*6.1429149; //[kg/m^3]
	crit.p = PressureUnit(5579.7, UNIT_KPA); //[kPa]
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

	const double a0[] = {1.26016e0/0.197578, -9.05307e-3/0.197578, 5.05504e-5/0.197578, -7.72237e-8/0.197578, 4.0538e-11/0.197578};
	const double n0[] = {0, 1, 2, 3, 4};
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(a0_v,n0_v,crit.T,298.15,0,4));

	name.assign("CycloPropane");
	aliases.push_back(std::string("cyclopropane"));
	aliases.push_back(std::string("Cyclopropane"));
	aliases.push_back(std::string("CYCLOPROPANE"));
	REFPROPname.assign("CYCLOPRO");

	BibTeXKeys.EOS = "Polt-CT-1992";
	BibTeXKeys.CP0 = "Polt-CT-1992";
}
double CycloPropaneClass::psat(double T)
{
    // Max error is  0.0119532667406 % between 273.0 and 398.299999 K
    const double t[]={0, 0.060000000000000005, 0.373, 0.8333333333333334, 0.38699999999999996, 1.0, 6.166666666666667};
    const double N[]={0, -0.0066049611475060624, -0.033960224958691061, -1.0723103687158344, 0.09566816778003541, -4.8279317520422547, -21.663297696353155};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
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

PropyneClass::PropyneClass()
{
	double n[] = {0, -9.86950667682466E-01, 4.59528109357232E+00, -8.86063623531859E+00, 5.56346955560542E+00, -1.57450028544218E+00, -1.59068753573430E-01, 2.35738270184042E-01, 4.40755494598713E-01, 1.96126150614333E-01, -3.67759650330219E-01, 7.92931851007852E-03, 2.47509085735293E-03, 8.32903610194452E-03, 1.02590136933231E+00, -2.20786016506394E+00, 1.07889905203761E+00, -3.82188466985900E+00, 8.30345065618981E+00, -4.48323072602860E+00, -1.02590136933231E+00, 2.20786016506394E+00, -1.07889905203761E+00};
	double t[] = {0, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5};
	double d[] = {0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 0, 0, 0, 2, 2, 2, 0, 0, 0};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2};
	double g[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.65533788, 1.65533788, 1.65533788, 1.65533788, 1.65533788, 1.65533788};

	//Critical parameters
	crit.rho = 40.06*6.11333; //[kg/m^3]
	crit.p = PressureUnit(5626, UNIT_KPA); //[kPa]
	crit.T = 402.38; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 40.06;
	params.Ttriple = 170.5;
	params.accentricfactor = 0.204;
	params.R_u = 8.3143;
	params.ptriple = 263.566905989;

	// Limits of EOS
	limits.Tmin = 273;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power(n,d,t,c,1,16,23));
	phirlist.push_back(new phir_exponential(n,d,t,c,g,17,22,23));

	phi0list.push_back(new phi0_lead(0,0));
	phi0list.push_back(new phi0_logtau(-1));

	double R = params.R_u/params.molemass;
	const double a0[] = {3.42418e-1/R, 4.84403e-3/R, -3.47414e-6/R, 1.44887e-9/R, -2.6815e-13/R};
	const double n0[] = {0, 1, 2, 3, 4};
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(a0_v,n0_v,crit.T,298.15,0,4));

	EOSReference.assign("");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("Propyne");
	aliases.push_back(std::string("propyne"));
	aliases.push_back(std::string("PROPYNE"));
	REFPROPname.assign("PROPYNE");
  
	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Polt-CT-1992";
	BibTeXKeys.CP0 = "Polt-CT-1992";
}
double PropyneClass::psat(double T)
{
    // Max error is  0.0257131564429 % between 273.0 and 402.379999 K
    const double t[]={0, 0.082, 0.363, 1.0, 1.1666666666666667, 3.3333333333333335, 13.333333333333334};
    const double N[]={0, 0.00010080009642318073, 0.002891896060186386, -7.3925772455032153, 1.1171903341456655, 0.15103683905132995, -9781.2864014399875};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double PropyneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.135668 % between 273.000001 K and 402.379990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.8333333333333333, 2.1666666666666665};
    const double N[] = {0, 4.0417756773327778, -31.191719018579267, 68.399467163263139, 211.93208284987557, -1457.8363818716066, 3452.436253269464, -4095.4407531126485, 2114.2042391704863, -402.8056702552895, 142.73254644405606};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=10; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double PropyneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.071856 % between 273.000001 K and 402.379990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.8333333333333333};
    const double N[] = {0, -3.2611479340583425, 21.113350221071236, -500.60815996337084, 2331.103035511363, -5545.9781937510725, 7874.0174835532853, -6568.2932354931727, 2641.1681543521377, -256.41310994446786};
    double summer = 0, theta;
    theta = 1-T/reduce.T;
	for (int i=1; i<=9; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}
