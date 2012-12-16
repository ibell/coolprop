#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "FAME.h"

static char EOSstr [] = "Marcia L. Huber, Eric W. Lemmon, Andrei Kazakov, Lisa S. Ott and Thomas J. Bruno, \"Model for the Thermodynamic Properties of a Biodiesel Fuel\", Energy & Fuels 2009, 23, 3790–3797";

static double eta[] = {0,0,0,0,0,0,0,0,0,0,0,1.1,1.6,1.1};
static double beta[] = {0,0,0,0,0,0,0,0,0,0,0,0.9,0.65,0.75};
static double gamma[] = {0,0,0,0,0,0,0,0,0,0,0,1.14,0.65,0.77};
static double epsilon[] = {0,0,0,0,0,0,0,0,0,0,0,0.79,0.90,0.76};

MethylPalmitateClass::MethylPalmitateClass()
{
	double n[] = {0.0, 0.4282821e-01, 0.2443162e01, -0.3757540e01, -0.1588526, 0.4055990e-01, -0.1524090e01, -0.7686167, 0.1799950e01, -0.1590967e01, -0.1267681e-01, 0.2198347e01, -0.7737211, -0.4314520};
	double d[] = {0.0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3};
	double t[] = {0.0, 1, 0.36, 1.22, 1.45, 0.7, 3, 3.9, 2.2, 2.9, 1.25, 2.6, 3.0, 3.2};
	double c[] = {0.0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 2, 2, 2};

	//Critical parameters
	crit.rho = 0.897*270.45066; //[kg/m^3]
	crit.p = 1350; //[kPa]
	crit.T = 755; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 270.45066; // From REFPROP, not provided in paper (but should be!!)
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
	phirlist.push_back(new phir_gaussian(n,d,t,eta,epsilon,beta,gamma,11,13,14));

	phi0list.push_back(new phi0_lead(-1,0));
	phi0list.push_back(new phi0_logtau(-1));

	const double u0[] = {0, 2952.37/crit.T, 734.653/crit.T, 1593.55/crit.T};
	const double v0[] = {0, 345.62/params.R_u, 289.038/params.R_u, 301.639/params.R_u};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(120.529/params.R_u,0.0801627,crit.T,298));
	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,3));

	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("MethylPalmitate");
	REFPROPname.assign("MPALMITA");
}
double MethylPalmitateClass::psat(double T)
{
    // Maximum absolute error is 0.319997 % between 302.710001 K and 754.999990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-13.200136570613084, 10.232230408559996, -13.252398211520738, 3.8775597911663158, -10.836717437005664, 3.2428420037390784 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double MethylPalmitateClass::rhosatL(double T)
{
    // Maximum absolute error is 1.570443 % between 302.710001 K and 754.999990 K
    const double ti[]={0,0.94121321520207935, 1.179965834421876, 1.0335314735150354, 1.1331135023608645, 1.0731528828525154};
    const double Ni[]={0,3919.7677554060283, 13903.086414838357, -42231.757132307117, -44061.791844712265, 68472.286515812797};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double MethylPalmitateClass::rhosatV(double T)
{
    // Maximum absolute error is 1.817374 % between 302.710001 K and 754.999990 K
    const double ti[]={0,0.99777618191177053, 0.88325353722700706, 1.4297193748563828, 2.2337495332277939, 2.7022142965879965};
    const double Ni[]={0,187.10548874715568, -135.41079368224243, -92.128414912890577, 79.291535420423813, -58.151064811192136};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}



MethylStearateClass::MethylStearateClass()
{
	double n[] = {0.0, 0.3959635e-01, 0.2466654e01, -0.3895950e01, -0.1167375, 0.4127229e-01, -0.1403734e01, -0.6465264, 0.1934675e01, -0.1608124e01, -0.1113813e-01, 0.2125325e01, -0.7772671, -0.4183684};
	double d[] = {0.0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3};
	double t[] = {0.0, 1, 0.3, 1.25, 1.65, 0.8, 3.1, 3.4, 2.3, 3.8, 1.2, 3.2, 3.8, 3.8};
	double c[] = {0.0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 2, 2, 2};

	//Critical parameters
	crit.rho = 0.7943*298.50382; //[kg/m^3]
	crit.p = 1350; //[kPa]
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
	phirlist.push_back(new phir_gaussian(n,d,t,eta,epsilon,beta,gamma,11,13,14));

	phi0list.push_back(new phi0_lead(-1,0));
	phi0list.push_back(new phi0_logtau(-1));

	const double u0[] = {0, 556.17/crit.T, 1311.85/crit.T, 2825.71/crit.T};
	const double v0[] = {0, 276.94/params.R_u, 408.997/params.R_u, 472.702/params.R_u};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(247.115/params.R_u,-0.0916606,crit.T,298));
	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,3));

	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("MethylStearate");
	REFPROPname.assign("MSTEARAT");
}
double MethylStearateClass::psat(double T)
{
    // Maximum absolute error is 0.322583 % between 311.840001 K and 774.999990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-14.459239434241695, 12.37912791111332, -15.719518668396431, 4.217170718278541, -9.4079764656586296, 0.74036836337369816 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double MethylStearateClass::rhosatL(double T)
{
    // Maximum absolute error is 1.077838 % between 311.840001 K and 774.999990 K
    const double ti[]={0,0.99667371444881436, 1.0192006184798761, 1.1773853109891423, 1.1356383451447185, 1.2001188677161669};
    const double Ni[]={0,24351.830138590401, -37610.270404646813, -112006.8916967539, 68194.912021600103, 57072.048148294067};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double MethylStearateClass::rhosatV(double T)
{
    // Maximum absolute error is 2.073991 % between 311.840001 K and 774.999990 K
    const double ti[]={0,1.263349122460272, 1.0269387325970791, 1.2018903038634003, 3.457385121946936, 1.1000160007058244};
    const double Ni[]={0,2166.2353330100868, -1677.383335860987, -4629.1357796869388, -13.124058184257374, 4132.7636009779553};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}


MethylOleateClass::MethylOleateClass()
{
	double n[] = {0.0, 0.4596121e-01, 0.2295400e01, -0.3554366e01, -0.2291674, 0.6854534e-01, -0.1535778e01, -0.7334697, 0.1712700e01, -0.1471394e01, -0.1724678e-01, 0.2115470e01, -0.7555374, -0.4134269};
	double d[] = {0.0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3};
	double t[] = {0.0, 1, 0.34, 1.14, 1.4, 0.6, 3.3, 4.1, 1.9, 3.8, 1.3, 3.4, 3.8, 4};
	double c[] = {0.0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 2, 2, 2};

	//Critical parameters
	crit.rho = 0.81285*296.48794; //[kg/m^3]
	crit.p = 1246; //[kPa]
	crit.T = 782; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 296.48794; // From REFPROP, not provided in paper (but should be!!)
	params.Ttriple = 253.47; // From REFPROP, not provided in paper
	params.accentricfactor = 0.90584935998790117;
	params.R_u = 8.314472;
	params.ptriple = 3.7666888692665290e-010;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,10,14));
	phirlist.push_back(new phir_gaussian(n,d,t,eta,epsilon,beta,gamma,11,13,14));

	phi0list.push_back(new phi0_lead(-1,0));
	phi0list.push_back(new phi0_logtau(-1));

	const double u0[] = {0, 613.529/crit.T, 1405.31/crit.T, 2867.76/crit.T};
	const double v0[] = {0, 234.797/params.R_u, 335.768/params.R_u, 431.66/params.R_u};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(90.2385/params.R_u,0.146118,crit.T,298));
	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,3));

	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("MethylOleate");
	REFPROPname.assign("MOLEATE");
}
double MethylOleateClass::psat(double T)
{
    // Maximum absolute error is 1.390482 % between 253.470001 K and 781.999990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-13.442936432369059, 11.217243033108137, -14.142996397314308, 1.197878318149272, -5.227226412331718, -6.1399376753358634 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double MethylOleateClass::rhosatL(double T)
{
    // Maximum absolute error is 1.447773 % between 253.470001 K and 781.999990 K
    const double ti[]={0,0.94937796824402665, 1.1560674169934608, 1.107690161992271, 1.1772919406079165, 1.0735915727003718};
    const double Ni[]={0,2861.7320728922614, -164917.71707828148, 154829.05748253479, 83674.2097477434, -76445.69688436162};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double MethylOleateClass::rhosatV(double T)
{
    // Maximum absolute error is 1.734744 % between 253.470001 K and 781.999990 K
    const double ti[]={0,0.83103666935662368, 1.1245342629525954, 0.9247654651922157, 18.293052022059015, 4.1232486312225252};
    const double Ni[]={0,-138.23273376521226, -72.279725161050678, 200.20150501157326, -86.616459150126502, -12.030142644202229};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}


MethylLinoleateClass::MethylLinoleateClass()
{
	double n[] = {0.0, 0.3183187e-01, 0.1927286e01, -0.3685053e01, 0, 0.8449312e-01, -0.9766643, -0.4323178, 0.2000470e01, -0.1752030e01, -0.1726895e-01, 0.2116515e01, -0.7884271, -0.3811699};
	double d[] = {0.0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3};
	double t[] = {0.0, 1, 0.2, 1.2, 1, 1, 2.2, 2.5, 1.8, 1.92, 1.47, 1.7, 2.3, 2.1};
	double c[] = {0.0, 0, 0, 0, 2, 3, 1, 3, 2, 2, 7, 1, 2, 2};

	//Critical parameters
	crit.rho = 0.8084*294.47206; //[kg/m^3]
	crit.p = 1341; //[kPa]
	crit.T = 799; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 294.47206; // From REFPROP, not provided in paper (but should be!!)
	params.Ttriple = 238.1; // From REFPROP, not provided in paper
	params.accentricfactor = 0.80540638705564849;
	params.R_u = 8.314472;
	params.ptriple = 7.7198861810445706e-012;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,10,14));
	phirlist.push_back(new phir_gaussian(n,d,t,eta,epsilon,beta,gamma,11,13,14));

	phi0list.push_back(new phi0_lead(-1,0));
	phi0list.push_back(new phi0_logtau(-1));

	const double u0[] = {0, 3052.11/crit.T, 746.631/crit.T, 1624.33/crit.T};
	const double v0[] = {0, 437.371/params.R_u, 287.222/params.R_u, 321.956/params.R_u};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(190.986/params.R_u,0.020213,crit.T,298));
	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,3));

	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("MethylLinoleate");
	REFPROPname.assign("MLINOLEA");
}
double MethylLinoleateClass::psat(double T)
{
    // Maximum absolute error is 2.082260 % between 238.100001 K and 798.999990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-10.829759464094929, 4.139222605621188, -3.1010581570822278, -12.412418569076836, 9.9566411553044443, -17.618879595792713 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double MethylLinoleateClass::rhosatL(double T)
{
    // Maximum absolute error is 2.985828 % between 238.100001 K and 798.999990 K
    const double ti[]={0,0.67556883441082527, 1.0333449523440372, 0.99280162108812076, 0.98260403923198025, 1.1327290791382629};
    const double Ni[]={0,82.364542322461901, -20412.168602636255, 78938.606417978765, -60308.648358807543, 1701.4521743774508};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double MethylLinoleateClass::rhosatV(double T)
{
	// Maximum absolute error is 1.27997424387 % between 238.100001 K and 798.999990 K
    const double ti[]={0,0.115, 0.3595, 0.39799999999999996, 2.0, 5.333333333333333, 14.0};
    const double Ni[]={0,-0.15852177442423029, 14.437419845171069, -20.129739979921837, -19.606999329430558, -100.26985364196683, -452.37641056284605};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}


MethylLinolenateClass::MethylLinolenateClass()
{
	double n[] = {0.0, 0.4070829e-01, 0.2412375e01, -0.3756194e01, -0.1526466, 0.4682918e-01, -0.1470958e01, -0.7645500, 0.1908964e01, -0.1629366e01, -0.1242073e-01, 0.2180707e01, -0.7537264, -0.4347781};
	double d[] = {0.0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3};
	double t[] = {0.0, 1, 0.15, 1.24, 1.6, 1.28, 2.9, 3.15, 2.16, 2.8, 1.4, 2.5, 3, 3.1};
	double c[] = {0.0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 2, 2, 2};

	//Critical parameters
	crit.rho = 0.8473*292.45618; //[kg/m^3]
	crit.p = 1369; //[kPa]
	crit.T = 772; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 292.45618; // From REFPROP, not provided in paper (but should be!!)
	params.Ttriple = 218.65; // From REFPROP, not provided in paper
	params.accentricfactor = 1.1426052586734956;
	params.R_u = 8.314472;
	params.ptriple = 8.2813864418489102e-015;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,10,14));
	phirlist.push_back(new phir_gaussian(n,d,t,eta,epsilon,beta,gamma,11,13,14));

	phi0list.push_back(new phi0_lead(-1,0));
	phi0list.push_back(new phi0_logtau(-1));

	const double u0[] = {0, 1213.24/crit.T, 578.752/crit.T, 2799.79/crit.T};
	const double v0[] = {0, 290.379/params.R_u, 81.4323/params.R_u, 474.881/params.R_u};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(79.5913/params.R_u,0.214648,crit.T,298));
	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,3));

	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("MethylLinolenate");
	REFPROPname.assign("MLINOLEN");
}
double MethylLinolenateClass::psat(double T)
{
    // Maximum absolute error is 1.802313 % between 218.650001 K and 771.999990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-14.283400982043565, 9.2026490474933293, -9.3608488932207656, -8.3518807964920327, 9.8477856505972063, -16.403332986338672 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double MethylLinolenateClass::rhosatL(double T)
{
    // Maximum absolute error is 0.964872 % between 218.650001 K and 771.999990 K
    const double ti[]={0,0.86810315764419388, 0.81617509367367369, 0.8581917521048037, 1.6378897563909396, 1.7224731276916094};
    const double Ni[]={0,6266.7706628183305, 1329.8340171916366, -7580.8784420604525, -69.004729231325513, 54.873347153886208};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double MethylLinolenateClass::rhosatV(double T)
{
    // Maximum absolute error is 1.649883 % between 218.650001 K and 771.999990 K
    const double ti[]={0,0.97321524606468035, 0.95174096264607333, 0.92905617616399805, 11.439045948636716, 3.3562768628976345};
    const double Ni[]={0,-3249.8679572470774, 6317.5377524765099, -3077.573870048132, -16.257670459099018, -10.683025812276989};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
