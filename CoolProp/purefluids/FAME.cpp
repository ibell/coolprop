#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "FAME.h"

MethylPalmitateClass::MethylPalmitateClass()
{
	double n[] = {0.0, 0.4282821e-01, 0.2443162e01, -0.3757540e01, -0.1588526, 0.4055990e-01, -0.1524090e01, -0.7686167, 0.1799950e01, -0.1590967e01, -0.1267681e-01, 0.2198347e01, -0.7737211, -0.4314520};
	double d[] = {0.0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3};
	double t[] = {0.0, 1, 0.36, 1.22, 1.45, 0.7, 3, 3.9, 2.2, 2.9, 1.25, 2.6, 3.0, 3.2};
	double c[] = {0.0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 2, 2, 2};
	static double eta[] = {0,0,0,0,0,0,0,0,0,0,0,1.1,1.6,1.1};
	static double beta[] = {0,0,0,0,0,0,0,0,0,0,0,0.9,0.65,0.75};
	static double _gamma[] = {0,0,0,0,0,0,0,0,0,0,0,1.14,0.65,0.77};
	static double epsilon[] = {0,0,0,0,0,0,0,0,0,0,0,0.79,0.90,0.76};

	//Critical parameters
	crit.rho = 0.897*270.45066; //[kg/m^3]
	crit.p = PressureUnit(1350, UNIT_KPA); //[kPa]
	crit.T = 755; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 270.45066; // http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=8181#x27 (not provided in reference paper)
	params.Ttriple = 302.71; // From REFPROP, not provided in paper
	params.accentricfactor = 0.910316178;
	params.R_u = 8.314472;
	params.ptriple = 1.6401417670571351e-005;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,10,14));
	phirlist.push_back(new phir_gaussian(n,d,t,eta,epsilon,beta,_gamma,11,13,14));

	phi0list.push_back(new phi0_lead(-1,0));
	phi0list.push_back(new phi0_logtau(-1));

	const double u0[] = {0, 2952.37/crit.T, 734.653/crit.T, 1593.55/crit.T};
	const double v0[] = {0, 345.62/params.R_u, 289.038/params.R_u, 301.639/params.R_u};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(120.529/params.R_u,0.0801627,crit.T,298));
	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,3));

	EOSReference.assign("Marcia L. Huber, Eric W. Lemmon, Andrei Kazakov, Lisa S. Ott and Thomas J. Bruno, \"Model for the Thermodynamic Properties of a Biodiesel Fuel\", Energy & Fuels 2009, 23, 3790-3797");
	TransportReference.assign("Using ECS in fully predictive mode");

	ECSReferenceFluid = "Propane";

	name.assign("MethylPalmitate");
	aliases.push_back(std::string("METHYLPALMITATE"));
	REFPROPname.assign("MPALMITA");

	BibTeXKeys.EOS = "Huber-EF-2009";
}
double MethylPalmitateClass::psat(double T)
{
    // Maximum absolute error is 0.070151 % between 302.710001 K and 754.999990 K
    const double t[]={0, 1, 2, 3, 4, 6, 9, 11, 16};
    const double N[]={0, 0.083410990578241709, -14.613736697709468, 18.580645505755221, -22.712044922611589, 14.04235195679173, -39.379094664340776, 35.342507811099267, -17.247819090316987};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=8;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double MethylPalmitateClass::rhosatL(double T)
{
    // Maximum absolute error is 0.007307 % between 302.710001 K and 754.999990 K
    const double t[] = {0, 1.1666666666666667, 1.5, 1.8333333333333333, 2.6666666666666665, 3.5, 4.0, 4.833333333333333, 6.166666666666667, 6.833333333333333, 8.0, 9.166666666666666, 9.833333333333334, 0.060000000000000005, 0.07, 0.1, 0.16999999999999998, 0.4, 0.56, 0.7400000000000001};
    const double N[] = {0, -58108.800134320532, 94180.233245431518, -79572.801061490871, 75752.670383125995, -223566.33850559776, 312678.16645895044, -246748.24685768026, 388191.43607787171, -428069.82081930753, 237639.90465013348, -146725.07535923398, 56520.581117125519, 73771.341876351275, -128522.86289250721, 77497.56006890304, -30448.044056061943, 30564.730350469155, -54367.464275048667, 49463.823895339076};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
	for (i=1; i<=19; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double MethylPalmitateClass::rhosatV(double T)
{
    // Maximum absolute error is 0.001684 % between 302.710001 K and 754.999990 K
    const double t[] = {0, 1.1666666666666667, 1.5, 1.8333333333333333, 2.1666666666666665, 2.6666666666666665, 3.5, 4.0, 4.833333333333333, 5.333333333333333, 6.166666666666667, 6.833333333333333, 8.0, 9.166666666666666, 9.666666666666666, 9.833333333333334, 0.060000000000000005, 0.07, 0.1, 0.16999999999999998, 0.31, 0.4, 0.7400000000000001};
    const double N[] = {0, -1173849.2030583415, 3450276.4030197775, -6911104.3378707729, 8142481.7118763933, -7011858.4112871401, 14139914.657616941, -22973786.643250182, 43780568.168386273, -54186179.790569082, 46391211.20929148, -33312665.785844114, 16906307.380973198, -24042951.639010228, 40414367.531671472, -23844026.155786425, 2075549.6589411693, -3696932.4373428691, 2404278.827095673, -1220572.289794754, 1122552.655609034, -849206.49894605658, 395183.81437690801};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    	
	for (i=1; i<=22; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

MethylStearateClass::MethylStearateClass()
{
	double n[] = {0.0, 0.3959635e-01, 0.2466654e01, -0.3895950e01, -0.1167375, 0.4127229e-01, -0.1403734e01, -0.6465264, 0.1934675e01, -0.1608124e01, -0.1113813e-01, 0.2125325e01, -0.7772671, -0.4183684};
	double d[] = {0.0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3};
	double t[] = {0.0, 1, 0.3, 1.25, 1.65, 0.8, 3.1, 3.4, 2.3, 3.8, 1.2, 3.2, 3.8, 3.8};
	double c[] = {0.0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 2, 2, 2};
	static double eta[] = {0,0,0,0,0,0,0,0,0,0,0,1.1,1.6,1.1};
	static double beta[] = {0,0,0,0,0,0,0,0,0,0,0,0.9,0.65,0.75};
	static double _gamma[] = {0,0,0,0,0,0,0,0,0,0,0,1.14,0.65,0.77};
	static double epsilon[] = {0,0,0,0,0,0,0,0,0,0,0,0.79,0.90,0.76};

	//Critical parameters
	crit.rho = 0.7943*298.50382; //[kg/m^3]
	crit.p = PressureUnit(1239, UNIT_KPA); //[kPa]
	crit.T = 775; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 298.50382; // From REFPROP, not provided in paper (but should be!!)
	params.Ttriple = 311.84; // From REFPROP, not provided in paper
	params.accentricfactor = 1.0548242393764551;
	params.R_u = 8.314472;
	params.ptriple = 6.0109170319097108e-006;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,10,14));
	phirlist.push_back(new phir_gaussian(n,d,t,eta,epsilon,beta,_gamma,11,13,14));

	phi0list.push_back(new phi0_lead(-1,0));
	phi0list.push_back(new phi0_logtau(-1));

	const double u0[] = {0, 556.17/crit.T, 1311.85/crit.T, 2825.71/crit.T};
	const double v0[] = {0, 276.94/params.R_u, 408.997/params.R_u, 472.702/params.R_u};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(247.115/params.R_u,-0.0916606,crit.T,298));
	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,3));

	EOSReference.assign("Marcia L. Huber, Eric W. Lemmon, Andrei Kazakov, Lisa S. Ott and Thomas J. Bruno, \"Model for the Thermodynamic Properties of a Biodiesel Fuel\", Energy & Fuels 2009, 23, 3790-3797");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("MethylStearate");
	aliases.push_back(std::string("METHYLSTEARATE"));
	REFPROPname.assign("MSTEARAT");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Huber-EF-2009";
}
double MethylStearateClass::psat(double T)
{
	// Max error is  0.397566439203 % between 311.84 and 774.999999 K

    const double t[]={0, 0.14300000000000002, 0.3585, 0.38599999999999995, 1.3333333333333333, 3.0, 6.5};
    const double N[]={0, -1.9263224274547197, 52.202696741784585, -54.538101452238379, -5.5909231565565696, -7.4874462624058769, -6.5951175020983355};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1;i<=6;i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double MethylStearateClass::rhosatL(double T)
{
    // Maximum absolute error is 0.113259 % between 311.840000 K and 774.999990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.6666666666666667, 2.1666666666666665};
    const double N[] = {0, -2.6832450139941888, 82.902228478630065, -673.54375616477273, 2521.5321256247776, -4884.9974663474659, 4905.1903624150527, -2093.4654862717925, 164.05730980797304, -15.147352889107403};
    double summer=0,theta;
    theta=1-T/reduce.T;  	
	for (int i=1; i<=9; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double MethylStearateClass::rhosatV(double T)
{
    // Maximum absolute error is 0.374148 % between 311.840000 K and 774.999990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.6666666666666667, 2.1666666666666665, 2.5};
    const double N[] = {0, -37.619839972324314, 1301.622015374996, -13951.568553365443, 74366.125966100357, -228997.33885627409, 425858.84894835541, -463752.05805638887, 245590.08321458541, -47655.879185556769, 9709.4378268852961, -2453.2891902758365};
    double summer=0,theta;
    theta=1-T/reduce.T;	
	for (int i=1; i<=11; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

MethylOleateClass::MethylOleateClass()
{
	double n[] = {0.0, 0.4596121e-01, 0.2295400e01, -0.3554366e01, -0.2291674, 0.6854534e-01, -0.1535778e01, -0.7334697, 0.1712700e01, -0.1471394e01, -0.1724678e-01, 0.2115470e01, -0.7555374, -0.4134269};
	double d[] = {0.0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3};
	double t[] = {0.0, 1, 0.34, 1.14, 1.4, 0.6, 3.3, 4.1, 1.9, 3.8, 1.3, 3.4, 3.8, 4};
	double c[] = {0.0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 2, 2, 2};
	static double eta[] = {0,0,0,0,0,0,0,0,0,0,0,1.1,1.6,1.1};
	static double beta[] = {0,0,0,0,0,0,0,0,0,0,0,0.9,0.65,0.75};
	static double _gamma[] = {0,0,0,0,0,0,0,0,0,0,0,1.14,0.65,0.77};
	static double epsilon[] = {0,0,0,0,0,0,0,0,0,0,0,0.79,0.90,0.76};

	//Critical parameters
	crit.rho = 0.81285*296.48794; //[kg/m^3]
	crit.p = PressureUnit(1246, UNIT_KPA); //[kPa]
	crit.T = 782; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 296.48794; // From REFPROP, not provided in paper (but should be!!)
	params.Ttriple = 253.47; // From REFPROP, not provided in paper
	params.accentricfactor = 0.90584935998790117;
	params.R_u = 8.314472;
	params.ptriple = 7.283926368e-11;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,10,14));
	phirlist.push_back(new phir_gaussian(n,d,t,eta,epsilon,beta,_gamma,11,13,14));

	phi0list.push_back(new phi0_lead(-1,0));
	phi0list.push_back(new phi0_logtau(-1));

	const double u0[] = {0, 613.529/crit.T, 1405.31/crit.T, 2867.76/crit.T};
	const double v0[] = {0, 234.797/params.R_u, 335.768/params.R_u, 431.66/params.R_u};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(90.2385/params.R_u,0.146118,crit.T,298));
	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,3));

	EOSReference.assign("Marcia L. Huber, Eric W. Lemmon, Andrei Kazakov, Lisa S. Ott and Thomas J. Bruno, \"Model for the Thermodynamic Properties of a Biodiesel Fuel\", Energy & Fuels 2009, 23, 3790-3797");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("MethylOleate");
	aliases.push_back(std::string("METHYLOLEATE"));
	REFPROPname.assign("MOLEATE");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Huber-EF-2009";
}
double MethylOleateClass::rhosatL(double T)
{
    // Maximum absolute error is 0.006952 % between 253.470001 K and 781.999990 K
    const double t[] = {0, 1.0, 1.3333333333333333, 1.6666666666666667, 2.0, 2.8333333333333335, 3.5, 4.666666666666667, 5.666666666666667, 6.333333333333333, 7.333333333333333, 9.0, 9.833333333333334, 0.060000000000000005, 0.07, 0.11, 0.2, 0.33, 0.48, 0.7400000000000001};
    const double N[] = {0, -82480.411020539439, 107216.30136367604, -119464.42716646775, 74758.669815127287, -37136.377716522707, 32615.163775178549, -33283.720322499285, 56542.35765044095, -51430.873604313907, 19101.091932241863, -5559.0770770202598, 2011.9321975503494, 49510.737936082929, -80951.053122789395, 52033.249828680928, -42349.774841519036, 49448.088151712676, -50893.080522091506, 60315.428789212368};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    	
	for (i=1; i<=19; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

double MethylOleateClass::psat(double T)
{
    // Maximum absolute error is 0.057939 % between 253.470001 K and 781.999990 K
    const double N[]={0, 0.092460604354544326, -15.395996742115694, 24.659586025922085, -41.999693851027232, 30.836783373822726, -20.146144258482281, 2.9679745676418419, -14.048162885056424};
    const double t[]={0, 1, 2, 3, 4, 5, 7, 12, 20};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=8;i++)
    {
        summer=summer+N[i]*pow(theta,t[i]/2);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double MethylOleateClass::rhosatV(double T)
{
    // Maximum absolute error is 0.012897 % between 253.470001 K and 781.999990 K
    const double t[] = {0, 1.0, 1.3333333333333333, 1.6666666666666667, 2.0, 2.5, 2.8333333333333335, 4.0, 4.666666666666667, 5.666666666666667, 7.333333333333333, 8.166666666666666, 9.666666666666666, 9.833333333333334, 0.060000000000000005, 0.07, 0.11, 0.2, 0.33, 0.7400000000000001};
    const double N[] = {0, -145202.55426929996, 343645.79528647946, -653014.95731831074, 765417.52981076075, -820402.83962917572, 561642.7436295189, -262067.81151749901, 266048.23662481701, -150543.20978059477, 149870.54089901687, -143698.20280116031, 212379.19492188256, -167787.23282325375, -78861.475157157474, 125586.76790272589, -72083.387416849699, 43168.527566726785, -25955.484677461325, 51789.218413949959};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
	for (i=1; i<=19; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

MethylLinoleateClass::MethylLinoleateClass()
{
	// Typo in the paper, coefficients for n,d,t,c taken from REFPROP
	double n[] = {0.0, 0.3183187e-01, 0.1927286e01, -0.3685053e01, 0.8449312e-01, -0.9766643, -0.4323178, 0.2000470e01, -0.1752030e01, -0.1726895e-01, 0.2116515e01, -0.7884271, -0.3811699};
	double d[] = {0.0, 4, 1, 1, 3, 1, 3, 2, 2, 7, 1, 1, 3};
	double t[] = {0.0, 1, 0.2, 1.2, 1, 2.2, 2.5, 1.8, 1.92, 1.47, 1.7, 2.3, 2.1};
	double c[] = {0.0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 2, 2, 2};
	static double eta[] = {0,0,0,0,0,0,0,0,0,0,1.1,1.6,1.1};
	static double beta[] = {0,0,0,0,0,0,0,0,0,0,0.9,0.65,0.75};
	static double gamma[] = {0,0,0,0,0,0,0,0,0,0,1.14,0.65,0.77};
	static double epsilon[] = {0,0,0,0,0,0,0,0,0,0,0.79,0.90,0.76};

	//Critical parameters
	crit.rho = 0.8084*294.47206; //[kg/m^3]
	crit.p = PressureUnit(1341, UNIT_KPA); //[kPa]
	crit.T = 799; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 294.47206; // From REFPROP, not provided in paper (but should be!!)
	params.Ttriple = 260;
	params.accentricfactor = 0.80540638705564849;
	params.R_u = 8.314472;
	params.ptriple = 1.37170480867e-09;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,9,13));
	phirlist.push_back(new phir_gaussian(n,d,t,eta,epsilon,beta,gamma,10,12,13));

	phi0list.push_back(new phi0_lead(-1,0));
	phi0list.push_back(new phi0_logtau(-1));

	const double u0[] = {0, 3052.11/crit.T, 746.631/crit.T, 1624.33/crit.T};
	const double v0[] = {0, 437.371/params.R_u, 287.222/params.R_u, 321.956/params.R_u};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(190.986/params.R_u,0.020213,crit.T,298));
	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,3));

	EOSReference.assign("Marcia L. Huber, Eric W. Lemmon, Andrei Kazakov, Lisa S. Ott and Thomas J. Bruno, \"Model for the Thermodynamic Properties of a Biodiesel Fuel\", Energy & Fuels 2009, 23, 3790-3797");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("MethylLinoleate");
	aliases.push_back(std::string("METHYLLINOLEATE"));
	REFPROPname.assign("MLINOLEA");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Huber-EF-2009";
}
double MethylLinoleateClass::psat(double T)
{
	// Max error is  0.244108713916 % between 260.0 and 798.999999 K
    const double t[]={0, 0.3555, 0.3755, 1.0, 4.5, 5.0, 7.833333333333333};
    const double N[]={0, 3.112535413528013, -3.5179254519786887, -8.5946903762416351, -60.06552122814378, 60.4039958334563, -23.566970610754186};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double MethylLinoleateClass::rhosatL(double T)
{
    // Maximum absolute error is 0.226557 % between 238.100000 K and 798.999990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.6666666666666667};
    const double N[] = {0, 0.36771210270684251, -6.6511365876605586, 13.236184198363402, 130.22823748163637, -482.25108598635057, 581.39341839244594, -159.14828458251716, -123.51883688912748, 50.25077404563865};
    double summer=0,theta;
    theta=1-T/reduce.T; 	
	for (int i=1; i<=9; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);

}

double MethylLinoleateClass::rhosatV(double T)
{
	// Max error is  0.245037263793 % between 260.0 and 798.999999 K
    const double t[] = {0, 0.369, 0.39199999999999996, 1.5, 4.0, 8.5, 9.666666666666666};
    const double N[] = {0, 16.679966544797427, -21.905607838963768, -3.2363141775511135, -14.340084125414874, 28.437620296244379, -43.146148011693199};
    double summer=0,theta;
    theta=1-T/reduce.T; 	
	for (int i=1; i<=6; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

MethylLinolenateClass::MethylLinolenateClass()
{
	double n[] = {0.0, 0.4070829e-01, 0.2412375e01, -0.3756194e01, -0.1526466, 0.4682918e-01, -0.1470958e01, -0.7645500, 0.1908964e01, -0.1629366e01, -0.1242073e-01, 0.2180707e01, -0.7537264, -0.4347781};
	double d[] = {0.0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3};
	double t[] = {0.0, 1, 0.15, 1.24, 1.6, 1.28, 2.9, 3.15, 2.16, 2.8, 1.4, 2.5, 3, 3.1};
	double c[] = {0.0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 2, 2, 2};
	static double eta[] = {0,0,0,0,0,0,0,0,0,0,0,1.1,1.6,1.1};
	static double beta[] = {0,0,0,0,0,0,0,0,0,0,0,0.9,0.65,0.75};
	static double _gamma[] = {0,0,0,0,0,0,0,0,0,0,0,1.14,0.65,0.77};
	static double epsilon[] = {0,0,0,0,0,0,0,0,0,0,0,0.79,0.90,0.76};

	//Critical parameters
	crit.rho = 0.8473*292.45618; //[kg/m^3]
	crit.p = PressureUnit(1369, UNIT_KPA); //[kPa]
	crit.T = 772; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 292.45618; // From REFPROP, not provided in paper (but should be!!)
	params.Ttriple = 260;
	params.accentricfactor = 1.1426052586734956;
	params.R_u = 8.314472;
	params.ptriple = 4.25028035475e-10;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,10,14));
	phirlist.push_back(new phir_gaussian(n,d,t,eta,epsilon,beta,_gamma,11,13,14));

	phi0list.push_back(new phi0_lead(-1,0));
	phi0list.push_back(new phi0_logtau(-1));

	const double u0[] = {0, 1213.24/crit.T, 578.752/crit.T, 2799.79/crit.T};
	const double v0[] = {0, 290.379/params.R_u, 81.4323/params.R_u, 474.881/params.R_u};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(79.5913/params.R_u,0.214648,crit.T,298));
	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,3));

	EOSReference.assign("Marcia L. Huber, Eric W. Lemmon, Andrei Kazakov, Lisa S. Ott and Thomas J. Bruno, \"Model for the Thermodynamic Properties of a Biodiesel Fuel\", Energy & Fuels 2009, 23, 3790-3797");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("MethylLinolenate");
	aliases.push_back(std::string("METHYLLINOLENATE"));
	REFPROPname.assign("MLINOLEN");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Huber-EF-2009";
}
double MethylLinolenateClass::psat(double T)
{
    // Maximum absolute error is 0.026348 % between 218.650001 K and 771.999990 K
    const double t[]={0, 1, 2, 3, 4, 5, 6, 8, 11, 14, 21, 29};
    const double N[]={0, 0.031999337164269448, -14.019556524509508, -0.38748218869676004, 65.309462079102019, -205.16139295149546, 229.91621086423717, -164.96257816043214, 137.7433741144975, -95.337171228166142, 33.734161081404864, -28.882979361404388};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=11;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double MethylLinolenateClass::rhosatL(double T)
{
    // Maximum absolute error is 0.063421 % between 218.650001 K and 771.999990 K
    const double t[] = {0, 1.0, 1.3333333333333333, 1.6666666666666667, 2.0, 3.0, 3.8333333333333335, 4.666666666666667, 6.0, 7.0, 8.0, 8.833333333333334, 9.833333333333334, 0.060000000000000005, 0.07, 0.15000000000000002, 0.25, 0.44, 0.7400000000000001};
    const double N[] = {0, -15056.649583969527, 19767.391560897006, -19630.986590668665, 10050.751249613639, -2998.0445898068133, 2786.7789098153958, -1948.0170104170948, 1307.5061040707369, -1018.3717799977615, 412.03373595244182, -29.280234808601588, -28.607699681171002, -4583.5693514792756, 6339.8046170477346, -4366.0400072001057, 4824.0512586189707, -4809.9007379394625, 8985.3288323119759};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
	for (i=1; i<=18; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}
double MethylLinolenateClass::rhosatV(double T)
{
    // Maximum absolute error is 0.008934 % between 218.650001 K and 771.999990 K
    const double t[] = {0, 1.0, 1.3333333333333333, 1.6666666666666667, 2.5, 3.0, 3.8333333333333335, 4.666666666666667, 6.0, 7.0, 8.0, 8.833333333333334, 9.666666666666666, 9.833333333333334, 0.060000000000000005, 0.07, 0.15000000000000002, 0.25, 0.35, 0.44, 0.7400000000000001};
    const double N[] = {0, 238018.88922882653, -201749.38968114922, 122562.20339273071, -89540.488608563377, 99127.889705029229, -79694.990085788202, 67087.360360970793, -78456.900191036955, 118472.42330309846, -145973.75512104615, 133387.19835817578, -154615.40572718397, 99781.727480516856, -175521.64136248664, 250147.89781487608, -233658.40932227639, 504260.23874569894, -786103.91597105679, 547821.29534493643, -235391.86543613058};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;    	
	for (i=1; i<=20; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}
