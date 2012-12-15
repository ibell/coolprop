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

static const double n[]={0,0.0555630,1.66927,-2.53408,-0.475075,0.190055,-1.25154,-0.742195,0.537902,-0.741246,-0.0355064,1.58506,-0.502086,-0.191360,-0.975576};

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
};

static const double t[]={0.00, 1.00,0.34,0.91,1.23,0.46,2.26,2.50,2.00,2.24,0.90,1.06,1.79,3.75,0.92};

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
};


static const double eta[]={
0,0,0,0,0,0,0,0,0,0,0, // indices [0-10]
1.02, //[11]
1.34, //[12]
1.08, //[13]
6.41, //[14]
};

static const double GAMMA[]={
0,0,0,0,0,0,0,0,0,0,0, // indices [0-10]
1.140, //[11]
0.667, //[12]
0.505, //[13]
1.220, //[14]
};

// epsilon is used here for consistency with the definitions in R744.c upon which Nitrogen.c is based
// is the value unity in Span
static const double beta[]={
0,0,0,0,0,0,0,0,0,0,0, // indices [0-10]
1.19, //[11]
2.29, //[12]
1.15, //[13]
131.8, //[14]
};

// GAMMA is used here for consistency with the definitions in R744.c upon which Nitrogen.c is based
static const double epsilon[]={
0,0,0,0,0,0,0,0,0,0,0, // indices [0-10]
0.711, //[11]
0.914, //[12]
0.694, //[13]
0.731, //[14]
};

//Constants for ideal gas expression
static const double a0[] = {0.0,5.8887,7.0804,9.3371,2.5577};
static const double b0[] = {0.0,0.0,620/382.52,1570/382.520,3953/382.52}; // Terms divided by critical temp

R1234zeClass::R1234zeClass()
{
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> d_v(d,d+sizeof(d)/sizeof(double));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(c,c+sizeof(c)/sizeof(double));
	std::vector<double> eta_v(eta,eta+sizeof(eta)/sizeof(double));
	std::vector<double> epsilon_v(epsilon,epsilon+sizeof(epsilon)/sizeof(double));
	std::vector<double> beta_v(beta,beta+sizeof(beta)/sizeof(double));
	std::vector<double> gamma_v(GAMMA,GAMMA+sizeof(GAMMA)/sizeof(double));	

	phirlist.push_back(new phir_power(n_v,d_v,t_v,l_v,1,10));
	phirlist.push_back(new phir_gaussian(n_v,d_v,t_v,eta_v,epsilon_v,beta_v,gamma_v,11,14));

	// phi0=log(delta)+a0[1]*log(tau)+a0[2]*log(1-exp(-b0*tau));
	phi0list.push_back(new phi0_lead(0,0)); // phi0_lead is like log(delta)+a1+a2*tau with a1=0, a2=0
	phi0list.push_back(new phi0_logtau(a0[1]));
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> b0_v(b0,b0+sizeof(b0)/sizeof(double));
	phi0list.push_back(new phi0_Planck_Einstein(a0_v,b0_v,2,4));

	// Critical parameters
	crit.rho = 4.29*114.0415928;
	crit.p = 3636.25;
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
	
	EOSReference.assign("Mark O. McLinden, Monika Thol, Eric W. Lemmon, \"Thermodynamic Properties of trans-1,3,3,3-tetrafluoropropene [R1234ze(E)]: Measurements of Density and Vapor Pressure and a Comprehensive Equation of State\", International Refrigeration and Air Conditioning Conference at Purdue, July 12-15, 2010,  ");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("R1234ze");
}
double R1234zeClass::psat(double T)
{
    // Maximum absolute error is 0.028234 % between 168.620001 K and 382.519990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.6110008912993861, 1.9964392788056975, -2.3807295148484506, -0.34611211743257886, -4.3951234494260758, 3.0431424571889147 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
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
//
//double SES36Class::surface_tension_T(double T)
//{
//	double sigma0 = 0.04193, sigma1 = 1.188, sigma2 = -1.462;
//	return sigma0*pow(1-T/reduce.T,1.26)*(1+sigma1*pow(1-T/reduce.T,0.5)+sigma2*(1-T/reduce.T));
//}