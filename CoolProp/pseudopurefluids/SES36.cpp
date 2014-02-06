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
#include "SES36.h"

SES36Class::SES36Class()
{
	//Constants for ideal gas expression
	static const double a0[]={0.0,13.09,85.26};
	static const double b0 = 2197;

	static const double n[]={0,    
	0.0675748, //[1]
	1.76939, //[2]
	-2.7609, //[3]
	-0.566938, //[4]
	0.243576, //[5]
	-1.50937, //[6]
	-0.774081, //[7]
	0.953907, //[8]
	-1.43736, //[9]
	-0.0458514, //[10]
	2.46053, //[11]
	-0.903158, //[12]
	-0.288664, //[13]
	0.061038, //[14]
	};

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

	// t used for consistency with CO2 correlation (corresponds to j from Span)
	static const double t[]={0.00,
	1, //[1]
	0.3, //[2]
	0.947, //[3]
	1.08, //[4]
	0.44, //[5]
	1.7, //[6]
	1.5, //[7]
	1.35, //[8]
	2.1, //[9]
	0.97, //[10]
	0.8, //[11]
	2, //[12]
	2.5, //[13]
	4, //[14]
	};

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

	// alpha is used here for consistency with the definitions in R744.c upon which Nitrogen.c is based
	// is phi_k from Span
	static const double eta[]={
	0,0,0,0,0,0,0,0,0,0,0, // indices [0-10]
	1.023, //[11]
	1.383, //[12]
	1, //[13]
	7, //[14]
	};

	static const double GAMMA[]={
	0,0,0,0,0,0,0,0,0,0,0, // indices [0-10]
	1.1, //[11]
	0.64, //[12]
	0.5, //[13]
	1.26, //[14]
	};

	// epsilon is used here for consistency with the definitions in R744.c upon which Nitrogen.c is based
	// is the value unity in Span
	static const double beta[]={
	0,0,0,0,0,0,0,0,0,0,0, // indices [0-10]
	1.7, //[11]
	1.55, //[12]
	1.07, //[13]
	87, //[14]
	};

	// GAMMA is used here for consistency with the definitions in R744.c upon which Nitrogen.c is based
	static const double epsilon[]={
	0,0,0,0,0,0,0,0,0,0,0, // indices [0-10]
	0.713, //[11]
	0.917, //[12]
	0.69, //[13]
	0.748, //[14]
	};

	std::vector<double> eta_v(eta,eta+sizeof(eta)/sizeof(double));
	std::vector<double> epsilon_v(epsilon,epsilon+sizeof(epsilon)/sizeof(double));
	std::vector<double> beta_v(beta,beta+sizeof(beta)/sizeof(double));
	std::vector<double> gamma_v(GAMMA,GAMMA+sizeof(GAMMA)/sizeof(double));

	phirlist.push_back(new phir_power(n,d,t,c,1,10,11));
	phirlist.push_back(new phir_gaussian(n,d,t,eta,epsilon,beta,GAMMA,11,14,15));

	// phi0=log(delta)+a0[1]*log(tau)+a0[2]*log(1-exp(-b0*tau));
	phi0list.push_back(new phi0_lead(0,0)); // phi0_lead is like log(delta)+a1+a2*tau with a1=0, a2=0
	phi0list.push_back(new phi0_logtau(a0[1]-1)); // -1 value needed to yield the correct values for the ideal gas part
	phi0list.push_back(new phi0_Planck_Einstein(a0[2],b0/450.7));

	// Critical parameters
	crit.rho = 2.8*184.85;
	crit.p = PressureUnit(2849,UNIT_KPA);
	crit.T = 450.7;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 184.85;
	params.Ttriple = 200; // No value given, taken from PPF Refprop file
	params.ptriple = 0.573325755702; //Evaluated at 200 K
	params.accentricfactor = 0.352;
	params.R_u = 8.314472;
	isPure = false;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 2000.0;
	limits.pmax = 2200000.0;
	limits.rhomax = 53.15*params.molemass;
	
	EOSReference.assign("Unpublished report: Monika Thol, Eric W. Lemmon, Roland Span, \"Equation of State for a Refrigerant Mixture of R365mfc (1,1,1,3,3-Pentafluorobutane) and Galden® HT 55 (Perfluoropolyether)\",  ");
	TransportReference.assign("Using ECS in predictive mode\n\n"
		                      "Surface Tension: A. P. Fröba, H. Kremer, A. Leipertz, F. Flohr and C. Meurer, "
		                      "\"Thermophysical Properties of a Refrigerant Mixture of R365mfc (1,1,1,3,3-Pentafluorobutane) "
							  "and Galden® HT 55 (Perfluoropolyether)\", International Journal of Thermophysics, Volume 28, "
							  "Number 2 (2007), 449-480, DOI: 10.1007/s10765-007-0178-y");

	name.assign("SES36");
    
    BibTeXKeys.EOS = "Thol-2012";
}
double SES36Class::psatL(double T)
{
	const double ti[]={0,1.0,1.5,2.5,4.5};
    const double Ni[]={0,-7.9188,2.4297,-5.1434,11.301};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer += Ni[i]*pow(1-T/reduce.T,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double SES36Class::psatV(double T)
{
	return psatL(T);
}
double SES36Class::rhosatL(double T)
{
    const double Ni[]={0,-0.3968,17.9889,-52.375,67.9231,-31.7301};
    double summer=0;
    int i;
    for (i=1;i<=5;i++)
    {
		summer += Ni[i]*pow(1.0-T/reduce.T,((double)i)/3.0);
    }
    return 538*(1+summer);
}
double SES36Class::rhosatV(double T)
{
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 0.722490 %
	RHS = -0.155721*pow(theta,0.000000)-13.669499*pow(theta,2.452194)-6.186741*pow(theta,0.635276)-5.454902*pow(theta,2.573221)-3.512970*pow(theta,5.915373)-1.760205*pow(theta,7.207006)-0.547867*pow(theta,7.850921);
	rho = exp(RHS)*reduce.rho;
	return rho;
}
double SES36Class::surface_tension_T(double T)
{
	double sigma0 = 0.04193, sigma1 = 1.188, sigma2 = -1.462;
	return sigma0*pow(1-T/reduce.T,1.26)*(1+sigma1*pow(1-T/reduce.T,0.5)+sigma2*(1-T/reduce.T));
}
