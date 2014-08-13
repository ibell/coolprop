/* Properties of SES36
by Ian Bell, 2012
*/

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include "math.h"
#include "stdio.h"
#include <string.h>
#include "CoolProp.h"
#include "FluidClass.h"
#include "R1234ze.h"

R1234zeClass::R1234zeClass()
{

	static const double n[]={0,4.434245E-02, 1.646369E+00, -2.437488E+00, -5.170560E-01, 1.815626E-01, -1.210104E+00, -5.944653E-01, 7.521992E-01, -6.747216E-01, -2.448183E-02, 1.379434E+00, -4.697024E-01, -2.036158E-01, -8.407447E-02, 5.109529E-04};

	// d used for consistency with CO2 correlation (corresponds to i from Span)
	static const double d[]={0,
	4, //[1]
	1, //[2]
	1, //[3]
	2, //[4]
	3, //[5]
	1, //[6]
	3, //[7]
	2, //[8]
	2, //[9]
	7, //[10]
	1, //[11]
	1, //[12]
	3, //[13]
	3, //[14]
	2, //[15
	};

	static const double t[]={0.00, 1.00,0.31,0.923,1.06,0.44,2.08,2.32,1.25,2.0,1.0,0.93,1.93,2.69,1.0,2.0};

	// c used for consistency with CO2 correlation (corresponds to l from Span)
	static const double c[]={0,
	0, //[1]
	0, //[2]
	0, //[3]
	0, //[4]
	0, //[5]
	2, //[6]
	2, //[7]
	1, //[8]
	2, //[9]
	1, //[10]
	0,0,0,0,0
	};

	static const double eta[]={
	0,0,0,0,0,0,0,0,0,0,0, // indices [0-10]
	1.0, //[11]
	1.4, //[12]
	1.134, //[13]
	7.68, //[14]
	24 //[15]
	};

	// epsilon is used here for consistency with the definitions in R744.c upon which Nitrogen.c is based
	// is the value unity in Span
	static const double beta[]={
	0,0,0,0,0,0,0,0,0,0,0, // indices [0-10]
	1.64, //[11]
	1.57, //[12]
	1.49, //[13]
	257.0, //[14]
	45, //[15]
	};

	static const double GAMMA[]={
	0,0,0,0,0,0,0,0,0,0,0, // indices [0-10]
	1.130, //[11]
	0.61, //[12]
	0.65, //[13]
	1.13, //[14]
	1.34 //[15]
	};

	// GAMMA is used here for consistency with the definitions in R744.c upon which Nitrogen.c is based
	static const double epsilon[]={
	0,0,0,0,0,0,0,0,0,0,0, // indices [0-10]
	0.711, //[11]
	0.856, //[12]
	0.753, //[13]
	0.772, //[14]
	1.88 //[15]
	};

	//Constants for ideal gas expression
	static const double a0[] = {0.0,6.259,7.303,8.597,2.333};
	static const double b0[] = {0.0,0.0,691/382.52,1705/382.520,4216/382.52}; // Terms divided by critical temp

	phirlist.push_back(new phir_power(n,d,t,c,1,10,16));
	phirlist.push_back(new phir_gaussian(n,d,t,eta,epsilon,beta,GAMMA,11,15,16));

	// phi0=log(delta)+a0[1]*log(tau)+a0[2]*log(1-exp(-b0*tau));
	phi0list.push_back(new phi0_lead(0,0)); // phi0_lead is like log(delta)+a1+a2*tau with a1=0, a2=0
	phi0list.push_back(new phi0_logtau(a0[1]));
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> b0_v(b0,b0+sizeof(b0)/sizeof(double));
	phi0list.push_back(new phi0_Planck_Einstein(a0_v,b0_v,2,4));

	// Critical parameters
	crit.rho = 4.29*114.0415928;
	crit.p = PressureUnit(3636.25, UNIT_KPA);
	crit.T = 382.52;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 114.0415928;
	params.Ttriple = 168.62;
	params.ptriple = 0.23;
	params.accentricfactor = 0.313;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 2000.0;
	limits.pmax = 2200000.0;
	limits.rhomax = 53.15*params.molemass;
	
	EOSReference.assign("Mark O. McLinden, Monika Thol, Eric W. Lemmon, \"Thermodynamic Properties of trans-1,3,3,3-tetrafluoropropene [R1234ze(E)]: Measurements of Density and Vapor Pressure and a Comprehensive Equation of State\", International Refrigeration and Air Conditioning Conference at Purdue, July 12-15, 2010");
	TransportReference.assign("Using ECS in fully predictive mode");

	
	name.assign("R1234ze(E)");
	aliases.push_back("R1234ZEE");
	aliases.push_back("R1234zeE");
	aliases.push_back(std::string("R1234ZE(E)"));
	REFPROPname.assign("R1234ZE");

	BibTeXKeys.EOS = "McLinden-PURDUE-2010";
	BibTeXKeys.CONDUCTIVITY = "Perkins-JCED-2011";
}
double R1234zeClass::psat(double T)
{
	// Max error is  0.164450626904 % between 168.62 and 382.519999 K
    const double ti[]={0, 0.3333333333333333, 0.39799999999999996, 0.38149999999999995, 1.1666666666666667, 2.5, 2.8333333333333335};
    const double Ni[]={0, -8.3421563693432379, -40.706881720081284, 47.898700907115519, -6.9428753841743065, 12.088882106622243, -14.335385281674474};
    double summer=0,theta;
    int i;
    theta=1 - T/reduce.T;
    for (i=1; i<=6; i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double R1234zeClass::rhosatL(double T)
{
    // Maximum absolute error is 0.630902 % between 168.620001 K and 382.519990 K
    const double ti[]={0,0.31509633757947603, 0.79111656912859174, 3.881097108987535, 7.7729566850460978, 2.6317441930191796};
    const double Ni[]={0,1.5412884996474399, -0.31992997858318134, -0.38423837781818643, 0.56607079429270613, 0.35255364714584947};
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
double R1234zeClass::rhosatV(double T)
{
    // Maximum absolute error is 0.449649 % between 168.620001 K and 382.519990 K
    const double ti[]={0,0.36034291391608658, 1.0029495942911635, 1.5238381484239822, 1.4251613173826694, 3.9498036504642666};
    const double Ni[]={0,-2.4332814409795982, -6.2539235529221004, -13.867477028090253, 16.426306950344401, -4.6886505464050012};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double R1234zeClass::conductivity_Trho(double T, double rho)
{
	double sumresid = 0;
	double A[] = {-0.0103589, 0.0308929, 0.000230348};
	double B1[] = {0, -0.0428296, 0.0927099, -0.0702107, 0.0249708, -0.00301838};
	double B2[] = {0, 0.0434288, -0.0605844, 0.0440187, -0.0155082, 0.00210190};

	double lambda_0 = A[0]*pow(T/reduce.T,0)+A[1]*pow(T/reduce.T,1)+A[2]*pow(T/reduce.T,2); //[W/m/K]

	for (int i = 1; i <= 5; i++)
	{
		sumresid += (B1[i]+B2[i]*T/reduce.T)*pow(rho/reduce.rho,i);
	}
	double lambda_r = sumresid; //[W/m/K]

	double lambda_c = this->conductivity_critical(T,rho,1/(0.585e-9)); //[W/m/K]

	return (lambda_0+lambda_r+lambda_c); //[W/m/K]
}












R1234zeZClass::R1234zeZClass()
{

	double n[] = {0, 0.77652368e01, -0.87025756e01, -0.28352251e00, 0.14534501e00, 0.92092105e-02, -0.24997382e00, 0.96674360e-01, 0.24685924e-01, -0.13255083e-01, -0.64231330e-01, 0.36638206e00, -0.25548847e00, -0.95592361e-01, 0.86271444e-01, 0.15997412e-01, -0.13127234e-01, 0.42293990e-02};
	double t[] = {0, 0.685, 0.8494, 1.87, 2., 0.142, 4.2, 0.08, 0.0, 1.1, 5.5, 6.6, 8.4, 7.2, 7.6, 8.5, 23.0, 18.};
	double d[] = {0, 1, 1, 1, 2, 5, 1, 3, 5, 7, 1, 2, 2, 3, 4, 2, 3, 5};
	double l[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3};
	
	phirlist.push_back(new phir_power(n,d,t,l,1,17,18));
	
	//Constants for ideal gas expression
	double t0[] = {0.0,0,1,2,3};
	double c0[] = {0.0,-1.6994, 24.527/423.27, -9.9249/423.27/423.27, 1.5158/423.27/423.27/423.27}; // Terms divided by critical temp
	std::vector<double> t0_v(t0,t0+sizeof(t0)/sizeof(double));
	std::vector<double> c0_v(c0,c0+sizeof(c0)/sizeof(double));
	phi0list.push_back(new phi0_lead(0,0)); // phi0_lead is like log(delta)+a1+a2*tau with a1=0, a2=0
	phi0list.push_back(new phi0_logtau(-1));
	phi0list.push_back(new phi0_cp0_poly(c0_v,t0_v,423.27,298,1,4));

	// Critical parameters
	crit.rho = 470;
	crit.p = PressureUnit(3533, UNIT_KPA);
	crit.T = 423.27;
	crit.v = 1.0/crit.rho;

	// Reducing parameters used in EOS
	reduce.p = PressureUnit(3533, UNIT_KPA);
	reduce.T = crit.T;
	reduce.rho = 470.615;
	reduce.v = 1.0/reduce.rho;

	preduce = &reduce;

	// Other fluid parameters
	params.molemass = 114.0415928;
	params.Ttriple = 273;
	params.ptriple = 67.8216;
	params.accentricfactor = 0.3274;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 2000.0;
	limits.pmax = 2200000.0;
	limits.rhomax = 53.15*params.molemass;

	name.assign("R1234ze(Z)");
	aliases.push_back("R1234ZE(Z)");
	REFPROPname = "N/A";

	BibTeXKeys.EOS = "Akasaka-DELFT-2013";
}
double R1234zeZClass::psat(double T)
{
    // Max error is  0.0362197220329 % between 273.0 and 423.269999 K
    const double t[]={0, 0.3675, 0.378, 0.8333333333333334, 1.1666666666666667, 3.5, 7.5};
    const double N[]={0, -3.4204004387395877, 3.7549450936800288, -4.2137431456192225, -3.0120279217374324, -2.6348957942753826, -6.8247581635654191};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double R1234zeZClass::rhosatL(double T)
{
    // Maximum absolute error is 0.109896 % between 273.000000 K and 423.207000 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.1666666666666667, 1.3333333333333333, 1.5};
    const double N[] = {0, 0.67225813646270216, -11.145228554195773, 77.35643204523889, -213.21240955222734, 262.97937232802644, -360.33055287979818, 362.6944623022427, -116.36985864325474};
    double summer=0,theta;
    theta=1-T/reduce.T;	
	for (int i=1; i<=8; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double R1234zeZClass::rhosatV(double T)
{
    // Maximum absolute error is 0.255509 % between 273.000000 K and 423.207000 K
    const double t[] = {0, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.6666666666666667, 1.8333333333333333};
    const double N[] = {0, 20.724125232917633, -334.7988571883364, 1630.4602570459422, -4085.9839380536523, 5626.8672277784344, -3717.4680060891073, 1458.3036712606865, -607.34248139716851};
    double summer=0,theta;
    theta=1-T/reduce.T;	
	for (int i=1; i<=8; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}
