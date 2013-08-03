#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "RC318_R21_R114_R13_R14.h"

RC318Class::RC318Class()
{
	double n[] = {0, 1.09415573664603E+00, -2.68265247887176E+00, 1.73403070801418E+00, -1.63611253452478E+00, 3.04834511239559E-01, 1.02771559991302E-01, -2.32367914087730E-02, 1.66151971110653E-01, -2.50103944487447E-02, 9.35681090149423E-02, 4.31929196624930E-02, -1.33439867188959E-01, 2.55416683254970E-02, -1.04729119771286E+00, 1.38034128674154E+00, -3.33625770182218E-01, -5.10485822124397E-01, 1.81840742529488E+00, -1.38530904967474E+00, 1.04729119771286E+00, -1.38034128674154E+00, 3.33625770182218E-01};
	double t[] = {0, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5};
	double d[] = {0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 0, 0, 0, 2, 2, 2, 0, 0, 0};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2};
	double g[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.99944, 0.99944, 0.99944, 0.99944, 0.99944, 0.99944};

	//Critical parameters
	crit.rho = 3.09938*200.0312; //[kg/m^3]
	crit.p = 3099.38; //[kPa]
	crit.T = 388.38; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 200.0312;
	params.Ttriple = 233.35;
	params.accentricfactor = 0.3553;
	params.R_u = 8.31451;
	params.ptriple = 19.461;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power(n,d,t,c,1,16,23));
	phirlist.push_back(new phir_exponential(n,d,t,c,g,17,22,23));

	phi0list.push_back(new phi0_lead(0,0));
	phi0list.push_back(new phi0_logtau(-1));

	double R = params.R_u/params.molemass;
	const double a0[] = {1.21e-1/R, 2.903e-3/R, -2.5327e-6/R, 7.7191e-10/R};
	const double n0[] = {0, 1, 2, 3, 4};
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(a0_v,n0_v,crit.T,298.15,0,3));

	name.assign("RC318");
	REFPROPname.assign("RC318");
  
	ECSReferenceFluid = "Propane";
	ECS_qd = 1/(0.356085e-9);

	BibTeXKeys.EOS = "Platzer-BOOK-1990";
	BibTeXKeys.CP0 = "Platzer-BOOK-1990";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double RC318Class::psat(double T)
{
    // Maximum absolute error is 0.014463 % between 233.350001 K and 388.379990 K
    const double t[]={0, 1, 2, 3, 4, 8};
    const double N[]={0, 0.0023402846842672542, -7.8996504796034381, 2.9361667746527802, -2.9819727078129805, -3.7079830133198022};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1;i<=4;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p*exp(reduce.T/T*summer);
}

double RC318Class::rhosatL(double T)
{
    // Maximum absolute error is 0.186636 % between 233.350001 K and 388.379990 K
    const double t[] = {0, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333};
    const double N[] = {0, 53.969051563565969, -197.00365240301335, 308.27001401289408, -229.95393573809568, 67.664593178058581};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=5; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double RC318Class::rhosatV(double T)
{
    // Maximum absolute error is 0.235209 % between 233.350001 K and 388.379990 K
    const double t[] = {0, 0.6666666666666666, 1.0, 1.1666666666666667, 1.3333333333333333, 1.6666666666666667};
    const double N[] = {0, -30.078583058035335, 296.09241082699941, -657.70383115910977, 480.88414718079628, -98.797638989440159};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=5; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

R21Class::R21Class()
{
	double n[] = {0, 5.04676623751950E-01, -7.32431416212257E-01, -8.68403860880684E-01, 1.46234705622829E-01, -2.80576335158724E-01, 8.64743657302055E-01, -2.70767234069002E+00, 3.30476391081085E+00, -2.10878239585469E-01, 4.49531450327333E-01, 1.20779813434928E-01, -2.77297954448155E-01, 3.05441292206290E-02, -4.43864848810647E+01, 9.26505601085111E+00, -5.51709104525115E-01, 1.21128809697908E+00, 1.67119476809363E-01, -5.04876793555323E-02, 4.43864848810647E+01, -9.26505601085111E+00, 5.51709104525115E-01};
	double t[] = {0, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5};
	double d[] = {0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 0, 0, 0, 2, 2, 2, 0, 0, 0};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2};
	double g[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.07470252, 0.07470252, 0.07470252, 0.07470252, 0.07470252, 0.07470252};

	//Critical parameters
	crit.rho = 5.1107656*102.9227; //[kg/m^3]
	crit.p = 5181.2; //[kPa]
	crit.T = 451.48; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 102.9227;
	params.Ttriple = 142.8;
	params.accentricfactor = 0.2061;
	params.R_u = 8.31451;
	params.ptriple = 0.872834974943;

	// Limits of EOS
	limits.Tmin = 200;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power(n,d,t,c,1,16,23));
	phirlist.push_back(new phir_exponential(n,d,t,c,g,17,22,23));

	phi0list.push_back(new phi0_lead(0,0));
	phi0list.push_back(new phi0_logtau(-1));

	double R = params.R_u/params.molemass;
	const double a0[] = {2.376576e-1/R, 1.271433e-3/R, 3.241352e-7/R, -2.492428e-9/R, 1.717208e-12/R};
	const double n0[] = {0, 1, 2, 3, 4};
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(a0_v,n0_v,crit.T,298.15,0,4));

	name.assign("R21");
	REFPROPname.assign("R21");
  
	ECSReferenceFluid = "R134a";

	BibTeXKeys.EOS = "Platzer-BOOK-1990";
	BibTeXKeys.CP0 = "Platzer-BOOK-1990";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R21Class::psat(double T)
{
    // Maximum absolute error is 0.099486 % between 200.000001 K and 451.479990 K
    const double t[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 18, 26};
    const double N[]={0, 0.094558463689140615, -10.027475351485936, 50.217905139105689, -502.56185292811892, 3415.7806161286767, -15504.258448977402, 47228.559805721052, -95048.131223023258, 118692.97399052662, -74947.50306001259, 22154.293393949894, -8414.6123934069092, 3440.3958662578339, -962.67543189295043};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1;i<=13;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p*exp(reduce.T/T*summer);
}

double R21Class::rhosatL(double T)
{
    // Maximum absolute error is 0.007797 % between 200.000001 K and 451.479990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.8333333333333333, 2.1666666666666665};
    const double N[] = {0, 8.1962989012252336, -126.69432982418643, 1052.8278007078461, -5339.195030700148, 17403.906139256076, -36983.75602624529, 50768.854983414814, -42575.104120024254, 17904.251888750605, -2429.6221450543758, 319.29162336471995};
    double summer=0,theta;
    theta=1-T/reduce.T;	
	for (int i=1; i<=11; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double R21Class::rhosatV(double T)
{
    // Maximum absolute error is 0.318686 % between 200.000001 K and 451.479990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.8333333333333333, 2.1666666666666665};
    const double N[] = {0, -11.965216242817384, 231.11620211923858, -2211.2434959327638, 12215.991323317601, -42096.616452521615, 93406.800295160952, -133237.72543166098, 115857.13120320938, -50458.426055336175, 7329.6177378866341, -1033.9002422581893};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=11; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

R114Class::R114Class()
{
	double n[] = {0, 1.07938940032879E+00, -1.99243731009857E+00, -1.55135220175231E-01, -1.21465858429101E-01, -1.65038674558161E-02, -1.86916017480622E-01, 3.08074956770949E-01, 1.15861545567346E-01, 2.76358779813890E-02, 1.08043424159349E-01, 4.60684822539207E-02, -1.74822007470687E-01, 3.17531741331376E-02, -3.40776521025164E-01, 3.23001398918284E-01, -4.24950543505197E-02, -1.66940287525002E+00, 4.08693538568874E+00, -2.41739233911370E+00, 3.40776521025164E-01, -3.23001398918284E-01, 4.24950543505197E-02};
	double t[] = {0, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5};
	double d[] = {0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 0, 0, 0, 2, 2, 2, 0, 0, 0};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2};
	double g[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.21104, 1.21104, 1.21104, 1.21104, 1.21104, 1.21104};

	//Critical parameters
	crit.rho = 3.3932*170.921; //[kg/m^3]
	crit.p = 3257.0; //[kPa]
	crit.T = 418.83; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 170.921;
	params.Ttriple = 180.63;
	params.accentricfactor = 0.2523;
	params.R_u = 8.31451;
	params.ptriple = 88.1623867625;

	// Limits of EOS
	limits.Tmin = 273.15;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power(n,d,t,c,1,16,23));
	phirlist.push_back(new phir_exponential(n,d,t,c,g,17,22,23));

	phi0list.push_back(new phi0_lead(0,0));
	phi0list.push_back(new phi0_logtau(-1));

	double R = params.R_u/params.molemass;
	const double a0[] = {9.765138e-2/R, 3.240861e-3/R, -5.895364e-6/R, 6.737929e-9/R, -3.546364e-12/R};
	const double n0[] = {0, 1, 2, 3, 4};
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(a0_v,n0_v,crit.T,298.15,0,4));

	name.assign("R114");
	REFPROPname.assign("R114");
  
	ECSReferenceFluid = "R134a";

	BibTeXKeys.EOS = "Platzer-BOOK-1990";
	BibTeXKeys.CP0 = "Platzer-BOOK-1990";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R114Class::psat(double T)
{
    // Maximum absolute error is 0.088440 % between 273.150001 K and 418.829990 K
    const double t[]={0, 1, 2, 3, 7};
    const double N[]={0, -0.015238961698740238, -7.0294363067308661, 0.9289618602692834, -3.6447987246533526};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1;i<=3;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p*exp(reduce.T/T*summer);
}

double R114Class::rhosatL(double T)
{
    // Maximum absolute error is 0.166287 % between 273.150001 K and 418.829990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5};
    const double N[] = {0, 8.0724210238996115, -96.761335399409077, 605.872851511985, -2291.2376321628576, 5509.0583210078776, -8366.137761880449, 7757.0256536890611, -4008.5271350615171, 885.72702898357159};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=9; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double R114Class::rhosatV(double T)
{
    // Maximum absolute error is 0.069288 % between 273.150001 K and 418.829990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.8333333333333333};
    const double N[] = {0, -9.6930811880906287, 117.96522621847706, -743.74723241355264, 2787.1929722048621, -6489.5985188809345, 9217.6641583192031, -7448.5109382036526, 2703.4677823598959, -143.67224624062317};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=9; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

R13Class::R13Class()
{
	double n[] = {0, 7.61142643481779E-01, -1.94465005029361E+00, 9.40938246527212E-01, -1.08106998126393E+00, 1.17501508367875E-01, 2.28304946999813E-01, -4.03338499859718E-01, 3.75585351263183E-01, -6.17542784117972E-02, 1.70325980529893E-01, 5.36611422401507E-02, -1.51602717890152E-01, 2.52032657506452E-02, -6.28346559897885E-01, 7.92797111317665E-01, -1.34038992719951E-01, -3.99863455306957E-02, 4.36410489713805E-01, -4.48724472286006E-01, 6.28346559897885E-01, -7.92797111317665E-01, 1.34038992719951E-01};
	double t[] = {0, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5};
	double d[] = {0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 0, 0, 0, 2, 2, 2, 0, 0, 0};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2};
	double g[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.9822996, 0.9822996, 0.9822996, 0.9822996, 0.9822996, 0.9822996};

	//Critical parameters
	crit.rho = 5.58*104.459; //[kg/m^3]
	crit.p = 3879.0; //[kPa]
	crit.T = 302.0; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 104.459;
	params.Ttriple = 92;
	params.accentricfactor = 0.1723;
	params.R_u = 8.31451;
	params.ptriple =  0.000896002851564;

	// Limits of EOS
	limits.Tmin = 98.15;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power(n,d,t,c,1,16,23));
	phirlist.push_back(new phir_exponential(n,d,t,c,g,17,22,23));

	phi0list.push_back(new phi0_lead(0,0));
	phi0list.push_back(new phi0_logtau(-1));

	double R = params.R_u/params.molemass;
	const double a0[] = {1.971309e-1/R, 1.438638e-3/R, 1.746775e-6/R, -6.830175e-9/R, 5.030396e-12/R};
	const double n0[] = {0, 1, 2, 3, 4};
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(a0_v,n0_v,crit.T,298.15,0,4));

	name.assign("R13");
	REFPROPname.assign("R13");
  
	ECSReferenceFluid = "Propane";
	ECS_qd = 1/(0.349636e-9);

	BibTeXKeys.EOS = "Platzer-BOOK-1990";
	BibTeXKeys.CP0 = "Platzer-BOOK-1990";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R13Class::psat(double T)
{
    // Maximum absolute error is 0.039440 % between 92.000001 K and 301.999990 K
    const double t[]={0, 1, 2, 3, 4, 6, 14, 23};
    const double N[]={0, -0.020600926775194043, -6.7055326313316268, 1.1165563585543907, 0.21670657964231119, -2.5524566919180347, -2.7436733707754968, 2.6945604974529349};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1;i<=6;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p*exp(reduce.T/T*summer);
}

double R13Class::rhosatL(double T)
{
    // Maximum absolute error is 0.024566 % between 92.000001 K and 301.999990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.8333333333333334, 1.0, 1.1666666666666667, 1.5, 1.8333333333333333};
    const double N[] = {0, 0.073472439822036884, -1.8276301609326788, 15.250990548179411, -114.39140927555329, 255.1979077969236, -197.87545054184721, 58.02985424357685, -11.786446645895108};
    double summer=0,theta;
    theta=1-T/reduce.T;   	
	for (int i=1; i<=8; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double R13Class::rhosatV(double T)
{
    // Maximum absolute error is 0.206366 % between 92.000001 K and 301.999990 K
    const double t[] = {0, 0.6666666666666666, 1.0, 1.1666666666666667, 1.5, 1.8333333333333333, 2.1666666666666665, 2.6666666666666665};
    const double N[] = {0, -36.979670238963841, 429.23940355104793, -835.52527231171302, 1016.3799655014985, -1028.9279725151914, 557.37521646657467, -110.75628417066763};
    double summer=0,theta;
    theta=1-T/reduce.T;	
	for (int i=1; i<=7; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

R14Class::R14Class()
{
	double n[] = {0, 1.03999039734947E+00, -2.45792025425590E+00, 7.99614558381524E-01, -7.49498955282594E-01, 1.52177772596398E-01, -2.93408332036298E-01, 7.17794503774445E-01, -4.26467444751902E-02, 2.26562749725952E-01, -3.91091694774020E-01, -2.57394805543082E-02, 5.54844886106659E-02, 6.10988262947185E-03, -3.34698748826510E-01, 5.86690904512625E-01, -1.47068929694523E-01, -1.90315426348019E-01, 7.16157134745809E-01, -7.03161905301754E-01, 3.34698748826510E-01, -5.86690904512625E-01, 1.47068929694523E-01};
	double t[] = {0, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5};
	double d[] = {0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 0, 0, 0, 2, 2, 2, 0, 0, 0};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2};
	double g[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.99832625, 0.99832625, 0.99832625, 0.99832625, 0.99832625, 0.99832625};

	//Critical parameters
	crit.rho = 7.1094194*88.0046; //[kg/m^3]
	crit.p = 3750.0; //[kPa]
	crit.T = 227.51; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 88.0046;
	params.Ttriple = 89.54;
	params.accentricfactor = 0.1785;
	params.R_u = 8.31451;
	params.ptriple =  11.2943290448;

	// Limits of EOS
	limits.Tmin = 120;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power(n,d,t,c,1,16,23));
	phirlist.push_back(new phir_exponential(n,d,t,c,g,17,22,23));

	phi0list.push_back(new phi0_lead(0,0));
	phi0list.push_back(new phi0_logtau(-1));

	double R = params.R_u/params.molemass;

	/// Refit in Excel based on REFPROP - seems the term in Platzer is way off
	const double a0[] = {3.728374e-01/R, -8.368995e-04/R, 1.316909e-05/R, -2.839480e-08/R, 1.937061e-11/R};
	const double n0[] = {0, 1, 2, 3, 4};
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

	phi0list.push_back(new phi0_cp0_poly(a0_v,n0_v,crit.T,298.15,0,4));

	name.assign("R14");
	REFPROPname.assign("R14");
  
	ECSReferenceFluid = "Nitrogen";
	ECS_qd = 1/(0.226566e-9);

	BibTeXKeys.EOS = "Platzer-BOOK-1990";
	BibTeXKeys.CP0 = "Platzer-BOOK-1990";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R14Class::psat(double T)
{
    // Maximum absolute error is 0.086626 % between 120.000001 K and 227.509990 K
    const double t[]={0, 1, 2, 3, 4, 5, 6, 9, 11, 19};
    const double N[]={0, 0.26304412336059479, -10.91771689667894, 29.08678750326963, -99.869651299377139, 188.80373892364335, -156.39025236765846, 105.33374357012718, -84.141452699003921, 51.949182533883437};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1;i<=8;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p*exp(reduce.T/T*summer);
}

double R14Class::rhosatL(double T)
{
    // Maximum absolute error is 0.128029 % between 120.000001 K and 227.509990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.8333333333333333, 2.1666666666666665, 2.5, 2.8333333333333335, 3.6666666666666665};
    const double N[] = {0, 199.3589449942568, -7047.335490112695, 92366.102607657609, -652077.6272373005, 2842186.9053680254, -8086044.4588733986, 15138109.773396848, -17804987.197905757, 10891965.511657136, -3675904.983491282, 1842460.5393635768, -741070.08139677474, 164387.52634921932, -4549.7773009256043};
    double summer=0,theta;
    theta=1-T/reduce.T;    	
	for (int i=1; i<=14; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double R14Class::rhosatV(double T)
{
    // Maximum absolute error is 0.212298 % between 120.000001 K and 227.509990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.8333333333333333, 2.1666666666666665, 2.5, 2.8333333333333335};
    const double N[] = {0, -100.7345626697176, 3444.3434008630493, -42830.254541518494, 284713.72386460466, -1162645.1073034331, 3084508.6432672697, -5359062.6069945302, 5819407.3052222747, -3268017.232780416, 910857.58529407543, -364141.81135523558, 111002.63233300112, -17148.507705086071};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=13; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}