#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Butenes.h"

static double d_BUTENES[] =
{
0,
1.0, //[1]
1.0, //[2]
1.0, //[3]
2.0, //[4]
3.0, //[5]
7.0, //[6]
2.0, //[7]
5.0, //[8]
1.0, //[9]
4.0, //[10]
3.0, //[11]
4.0, //[12]
};

static double t_BUTENES[] =
{
0,
0.12,  //[1]
1.3, //[2]
1.74,   //[3]
2.1, //[4]
0.28,  //[5]
0.69, //[6]
0.75, //[7]
2.0,  //[8]
4.4, //[9]
4.7, //[10]
15.0,  //[11]
14.0,  //[12]
};

static double c_BUTENES[] =
{
0,
0.0, //[1]
0.0, //[2]
0.0, //[3]
0.0, //[4]
0.0, //[5]
0.0, //[6]
1.0, //[7]
1.0, //[8]
2.0, //[9]
2.0, //[10]
3.0, //[11]
3.0, //[12]
};

static char EOSstr_BUTENES [] = "Eric W. Lemmon, E. Christian Ihmels, \"Thermodynamic properties of the butenes Part II. Short fundamental equations of state\", Fluid Phase Equilibria v. 228-229 (2005) 173-187";

OneButeneClass::OneButeneClass()
{
	double n[] = {0.0, 0.78084, -2.8258, 0.99403, 0.017951, 0.088889, 0.00024673, 0.22846, -0.074009, -0.22913, -0.062334, -0.025385, 0.011040};

	//Critical parameters
	crit.rho = 4.24*56.10632; //[kg/m^3]
	crit.p = PressureUnit(4005.1, UNIT_KPA); //[kPa]
	crit.T = 419.29; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 56.10632;
	params.Ttriple = 87.8;
	params.accentricfactor = 0.191860647355;
	params.R_u = 8.314472;
	params.ptriple = 5.94529945955e-10;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d_BUTENES,t_BUTENES,c_BUTENES,1,12,13));

	const double a1 = -0.00101126, a2 = 2.3869174, c0 = 3.9197;
	phi0list.push_back(new phi0_lead(a1,a2));
	phi0list.push_back(new phi0_logtau(c0-1));

	const double u0[] = {0, 274/419.29, 951/419.29, 2127/419.29, 5752/419.29};
	const double v0[] = {0, 2.9406, 6.5395, 14.535, 5.8971};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	EOSReference.assign(EOSstr_BUTENES);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("1-Butene");
	aliases.push_back(std::string("1Butene"));
	aliases.push_back(std::string("1BUTENE"));
	aliases.push_back(std::string("1-BUTENE"));
	aliases.push_back(std::string("Butene"));
	REFPROPname.assign("1BUTENE");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Lemmon-FPE-2005";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double OneButeneClass::psat(double T)
{
    // Maximum absolute error is 0.162304 % between 87.800001 K and 419.289990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3,9};
    const double Ni[]={0,-7.1450182048398165, 2.2562618900495846, -2.36160856975917, 0.73673176213172731, -5.6018424878741815, 6.121828365893009, -4.2551550478203906 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=7;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double OneButeneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.644715 % between 87.800001 K and 419.289990 K
    const double ti[]={0,0.34250876414210879, 0.36355989633783242, 35.398692701297669, 3.5910129321917239, 7.3896726410407201};
    const double Ni[]={0,5.3874235581501173, -4.1100854559360691, -1.3372589031304434, 0.094501262405661535, -0.01350642745145569};
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
double OneButeneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.161887 % between 87.800001 K and 419.289990 K
    const double ti[]={0,0.38479607874880906, 0.82029186054153236, 23.635840599651146, 3.1770700106715823, 4.4973781676010844};
    const double Ni[]={0,-2.4416835351244481, -2.5245469947685257, -4.2144539186977754, -1.7867328713084489, -2.7501625926745277};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

IsoButeneClass::IsoButeneClass()
{
	double n[] = {0.0, 0.77111, -2.7971, 1.0118, 0.020730, 0.085086, 0.00021968, 0.20633, -0.078843, -0.23726, -0.080211, -0.027001, 0.013072};

	//Critical parameters
	crit.rho = 4.17*56.10632; //[kg/m^3]
	crit.p = PressureUnit(4009.8, UNIT_KPA); //[kPa]
	crit.T = 418.09; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 56.10632;
	params.Ttriple = 132.4;
	params.accentricfactor = 0.1925934521621;
	params.R_u = 8.314472;
	params.ptriple = 0.000676189903044909;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d_BUTENES,t_BUTENES,c_BUTENES,1,12,13));

	const double a1 = -0.12737888, a2 = 2.3125128, c0 = 4.0000;
	phi0list.push_back(new phi0_lead(a1,a2));
	phi0list.push_back(new phi0_logtau(c0-1));

	const double u0[] = {0,399/418.09,1270/418.09,2005/418.09,4017/418.09};
	const double v0[] = {0, 4.8924, 7.8320, 7.2867, 8.7293};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	EOSReference.assign(EOSstr_BUTENES);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("IsoButene");
	aliases.push_back(std::string("Isobutene"));
	aliases.push_back(std::string("ISOBUTENE"));
	REFPROPname.assign("IBUTENE");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Lemmon-FPE-2005";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double IsoButeneClass::psat(double T)
{
    // Maximum absolute error is 0.141011 % between 132.400001 K and 418.089990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-6.9116852201182075, 1.3539793916115146, -0.6053845232409808, -2.4488184963383248, -0.52147628515333433, -0.561826127052658 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double IsoButeneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.911818 % between 132.400001 K and 418.089990 K
    const double ti[]={0,0.30342768948442322, 0.68061509707049073, 1.4525834599224288, 2.3785976948607295, 2.8237458291611888};
    const double Ni[]={0,1.4523353540554691, -0.1572370124750222, -0.069052815100666839, 0.1081982691615444, 0.030812753121678148};
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
double IsoButeneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.218153 % between 132.400001 K and 418.089990 K
    const double ti[]={0,0.40535252257257232, 1.0656571941879163, 1.3231477111666543, 2.4463428843919233, 4.4302956092284793};
    const double Ni[]={0,-2.6824183183460883, -6.1489644495471669, 4.7663394841609668, -2.2153438450173208, -3.1991563607011462};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

Cis2ButeneClass::Cis2ButeneClass()
{
	double n[] = {0.0, 0.77827, -2.8064, 1.0030, 0.013762, 0.085514, 0.00021268, 0.22962, -0.072442, -0.23722, -0.074071, -0.026547, 0.012032};

	//Critical parameters
	crit.rho = 4.244*56.10632; //[kg/m^3]
	crit.p = PressureUnit(4225.5, UNIT_KPA); //[kPa]
	crit.T = 435.75; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 56.10632;
	params.Ttriple = 134.3;
	params.accentricfactor = 0.20235958587;
	params.R_u = 8.314472;
	params.ptriple = 0.0002636498688682175;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d_BUTENES,t_BUTENES,c_BUTENES,1,12,13));

	const double a1 = 0.2591542, a2 = 2.4189888, c0 = 3.9687;
	phi0list.push_back(new phi0_lead(a1,a2));
	phi0list.push_back(new phi0_logtau(c0-1));

	const double u0[] = {0, 248/435.75, 1183/435.75, 2092/435.75, 4397/435.75};
	const double v0[] = {0, 3.2375, 7.0437, 11.414, 7.3722};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	EOSReference.assign(EOSstr_BUTENES);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("cis-2-Butene");
	aliases.push_back(std::string("Cis-2-Butene"));
	aliases.push_back(std::string("CIS-2-BUTENE"));
	REFPROPname.assign("C2BUTENE");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Lemmon-FPE-2005";
}
double Cis2ButeneClass::psat(double T)
{
    // Maximum absolute error is 0.248652 % between 134.300001 K and 435.749990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-6.9258343900515573, 1.221853780238076, -0.31550093997536949, -3.1568049145274149, 0.23805014289458684, -1.279772722501435 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double Cis2ButeneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.450819 % between 134.300001 K and 435.749990 K
    const double ti[]={0,0.44413848705698888, 0.65980763681444066, 0.64276060618274145, 0.64094401703342718, 2.5283005524045095};
    const double Ni[]={0,11.567206951223904, -1109.7758269438145, 12460.578738205559, -11361.152072565803, 0.16407435859133529};
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
double Cis2ButeneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.382800 % between 134.300001 K and 435.749990 K
    const double ti[]={0,0.44290471958214006, 0.9919707337807756, 0.031337536771287318, -4.4942695307610876, 4.0064150025287013};
    const double Ni[]={0,-3.194141602461662, -1.9845364039394378, -0.025383887119914638, 6.1673175860603378e-37, -4.5162906962263722};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

Trans2ButeneClass::Trans2ButeneClass()
{
	double n[] = {0.0, 0.81107, -2.8846, 1.0265, 0.016591, 0.086511, 0.00023256, 0.22654, -0.072182, -0.24849, -0.071374, -0.024737, 0.011843};

	//Critical parameters
	crit.rho = 4.213*56.10632; //[kg/m^3]
	crit.p = PressureUnit(4027.3, UNIT_KPA); //[kPa]
	crit.T = 428.61; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 56.10632;
	params.Ttriple = 167.6;
	params.accentricfactor = 0.21007683443616;
	params.R_u = 8.314472;
	params.ptriple = 0.07481669961927020;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d_BUTENES,t_BUTENES,c_BUTENES,1,12,13));

	const double a1 = 0.5917816, a2 = 2.1427758, c0 = 3.9988;
	phi0list.push_back(new phi0_lead(a1,a2));
	phi0list.push_back(new phi0_logtau(c0-1));

	const double u0[] = {0, 362/428.61, 1603/428.61, 3729/428.61, 4527/428.61};
	const double v0[] = {0, 5.3276, 13.290, 9.6745, 0.40087};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	EOSReference.assign(EOSstr_BUTENES);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("trans-2-Butene");
	aliases.push_back(std::string("Trans-2-Butene"));
	aliases.push_back(std::string("TRANS-2-BUTENE"));
	REFPROPname.assign("T2BUTENE");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Lemmon-FPE-2005";
}
double Trans2ButeneClass::psat(double T)
{
    // Maximum absolute error is 0.205049 % between 167.600001 K and 428.609990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.417038031043373, 2.89238779867215, -3.4941759407935451, 2.3076971643105062, -6.477416610557051, 3.8758283450498365 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double Trans2ButeneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.702266 % between 167.600001 K and 428.609990 K
    const double ti[]={0,0.36279985966962308, 0.35461690571642823, 3.2508260848176342, 3.3292336646129717, 3.2254735273475723};
    const double Ni[]={0,-13.329703571695561, 14.617556102115302, 75.644376981618194, -18.667171640293407, -56.908406625616792};
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
double Trans2ButeneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.398092 % between 167.600001 K and 428.609990 K
    const double ti[]={0,0.3547482282481132, 0.73232141326808986, 2.8067184296149144, 15.282642347756052, 4.6162704278588027};
    const double Ni[]={0,-1.9653237149532947, -2.9314254366049646, -1.7095509890716998, 5.2783324320756684, -2.9407203414470202};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
