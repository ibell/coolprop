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
	crit.p = PressureUnit(2.7775e+006, UNIT_PA); //[Pa]
	crit.T = 388.38; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 200.0312;
	params.Ttriple = 233.35;
	params.accentricfactor = 0.355345;
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
	// Max error is  0.00483960782838 % between 233.35 and 388.3799 K for RC318
    const double t[]={0, 0.353, 0.3585, 0.38899999999999996, 0.39499999999999996, 0.39649999999999996, 1.0, 2.6666666666666665, 4.166666666666667};
    const double N[]={0, -4692.4741793345729, 7715.2503820566762, -49464.460331933777, 188185.25751683454, -141744.53054929635, -5.7257304556280664, -1.8133107974678764, -3.2381748043228629};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=8; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
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
	crit.p = PressureUnit(5181.2, UNIT_KPA); //[kPa]
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
	// Max error is  0.11324761496 % between 200.0 and 451.479999 K
    const double t[]={0, 5.166666666666667, 0.5, 0.8333333333333334, 1.1666666666666667, 4.5, 15.5};
    const double N[]={0, 10.622994048827298, 0.24340767266291899, -3.6806156929627569, -2.8391976402596844, -11.694192853896972, -105.73416753738717};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
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
	crit.p = PressureUnit(3257.0, UNIT_KPA); //[kPa]
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
	// Max error is  0.0126361637689 % between 273.15 and 418.829999 K   
    const double t[]={0, 0.07400000000000001, 0.39849999999999997, 1.0, 1.3333333333333333, 1.6666666666666667, 3.6666666666666665};
    const double N[]={0, -0.0034645011523141735, 0.023157592718005465, -7.5688152387880105, 1.8469739464041124, -0.64855802507466331, -3.5293610438916647};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
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
	double n[] = {0, 7.61143010172E-01, -1.94465098795E+00, 9.40938700406E-01, -1.08107050239E+00, 1.17501564976E-01, 2.28305167217E-01, -4.03338888789E-01, 3.75585713420E-01, -6.17543677315E-02, 1.70326226881E-01, 5.36612457231E-02, -1.51603010301E-01, 2.52033265074E-02, -6.28346559920E-01, 7.92797111341E-01, -1.34038992692E-01, -3.99863840975E-02, 4.36410910529E-01, -4.48724904991E-01, 6.28346559920E-01, -7.92797111341E-01, 1.34038992692E-01};
	double t[] = {0, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5};
	double d[] = {0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 0, 0, 0, 2, 2, 2, 0, 0, 0};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2};
	double g[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.98230055, 0.98230055, 0.98230055, 0.98230055, 0.98230055, 0.98230055};

	//Critical parameters
	crit.rho = 5.58*104.459; //[kg/m^3]
	crit.p = PressureUnit(3879.0, UNIT_KPA); //[kPa]
	crit.T = 301.88; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 104.459;
	params.Ttriple = 92;
	params.accentricfactor = 0.174586327798;
	params.R_u = 8.31451;
	params.ptriple =  0.000907038339057;

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
	// Max error is 0.282358377629 % between 98.15 and 301.879999 K
    const double t[]={0, 0.383, 0.38999999999999996, 1.1666666666666667, 5.5, 0.3695, 23.333333333333332};
    const double N[]={0, 577.42725600666131, -400.21572594632761, -4.9874139876628041, -4.222492303172193, -178.69585671516839, -1720.3585941933911};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double R13Class::rhosatL(double T)
{
    // Max error is  0.090396245488 % between 98.15 and 301.879999 K
    const double t[] = {0, 0.113, 0.364, 0.3695, 0.3785, 1.3333333333333333, 6.666666666666667};
    const double N[] = {0, 2.4784634077003043, -4743.997952071094, 7615.934031115753, -2872.7859164099727, 1.0823469621368702, 0.12626188722376672};
    double summer=0,theta;
    theta=1-T/reduce.T;
    	
	for (int i=1; i<=6; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}
double R13Class::rhosatV(double T)
{
    // Max error is  0.0235817860212 % between 98.15 and 301.8799 K
    const double t[] = {0, 0.138, 0.14300000000000002, 0.38599999999999995, 0.38799999999999996, 0.38899999999999996, 0.39099999999999996, 0.39149999999999996, 4.0, 8.333333333333334, 21.5};
    const double N[] = {0, -1287.2952626479841, 1471.130388287577, -3382739168.0344748, 26487996709.339836, -37021553643.047638, 36901846787.213554, -22985550874.574944, -3.2136791021238045, -3.0923712297543609, -891.46368816511585};
    double summer=0,theta;
    theta=1-T/reduce.T;
    	
	for (int i=1; i<=10; i++)
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
	crit.p = PressureUnit(3750.0, UNIT_KPA); //[kPa]
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
	// Max error is  0.106188047838 % between 120.0 and 227.509999 K
    const double t[]={0, 0.056, 0.362, 0.379, 0.381, 1.0, 4.166666666666667};
    const double N[]={0, 0.022784802896396485, -47.710297004468607, 499.94408372159194, -452.62807116975046, -5.7182330867269302, -3.0735670244888449};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double R14Class::rhosatL(double T)
{
	// Max error is  0.442178557435 % between 120.0 and 227.5099 K

    const double t[] = {0, 0.069, 0.07600000000000001, 0.08, 0.082, 0.084, 0.136, 0.14800000000000002, 0.10900000000000001};
    const double N[] = {0, 1602566566.7780583, -26938218495.184853, 133716790849.44449, -173718665668.18839, 65516872229.093124, 26998449.766706958, -7249677.4848860716, -199094251.66107807};
    double summer=0,theta;
    theta=1-T/reduce.T;    	
	for (int i=1; i <= 8; i++)
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
