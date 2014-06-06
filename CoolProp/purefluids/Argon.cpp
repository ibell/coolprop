/* Properties of Argon
by Ian Bell

Thermo properties from 
---------------------
"A New Equation of State for Argon Covering the Fluid Region
for Temperatures From the Melting Line to 700 K
at Pressures up to 1000 MPa"
Ch. Tegeler, R. Span, and W. Wagner
J. Phys. Chem. Ref. Data, Vol. 28, No. 3, 1999

Transport properties from
-------------------------
"Viscosity and Thermal Conductivity Equations for
Nitrogen, Oxygen, Argon, and Air"
E. W. Lemmon and R. T Jacobsen
International Journal of Thermophysics, Vol. 25, No. 1, January 2004

Note: Critical enhancement included

*/

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
// The most important line
//#define new new(_NORMAL_BLOCK, __FILE__, __LINE__)
#else
#include <stdlib.h>
#endif

#include "math.h"
#include "stdio.h"
#include <string.h>
#include "CoolProp.h"
#include <vector>
#include <iostream>
#include <list>
#include "Helmholtz.h"
#include "FluidClass.h"
#include "Argon.h"

ArgonClass::ArgonClass()
{
	static const double n[]={0,
	0.088722304990011,//[1]
	0.70514805167298,//[2]
	-1.682011565409,//[3]
	-0.14909014431486,//[4]
	-0.1202480460094,//[5]
	-0.12164978798599,//[6]
	0.40035933626752,//[7]
	-0.27136062699129,//[8]
	0.24211924579645,//[9]
	0.005788958318557,//[10]
	-0.041097335615341,//[11]
	0.024710761541614,//[12]
	-0.32181391750702,//[13]
	0.33230017695794,//[14]
	0.031019986287345,//[15]
	-0.030777086002437,//[16]
	0.093891137419581,//[17]
	-0.090643210682031,//[18]
	-0.00045778349276654,//[19]
	-0.000082659729025197,//[20]
	0.00013013415603147,//[21]
	-0.011397840001996,//[22]
	-0.024455169960535,//[23]
	-0.064324067175955,//[24]
	0.058889471093674,//[25]
	-0.00064933552112965,//[26]
	-0.013889862158435,//[27]
	0.4048983929691,//[28]
	-0.38612519594749,//[29]
	-0.18817142332233,//[30]
	0.15977647596482,//[31]
	0.053985518513856,//[32]
	-0.028953417958014,//[33]
	-0.013025413381384,//[34]
	0.0028948696775778,//[35]
	-0.0022647134304796,//[36]
	0.0017616456196368,//[37]
	0.0058552454482774,//[38]
	-0.69251908270028,//[39]
	1.5315490030516,//[40]
	-0.0027380447449783//[41]
	};

	static const double d[]={0,
	1,//[1]
	1,//[2]
	1,//[3]
	1,//[4]
	1,//[5]
	2,//[6]
	2,//[7]
	2,//[8]
	2,//[9]
	3,//[10]
	3,//[11]
	4,//[12]
	1,//[13]
	1,//[14]
	3,//[15]
	4,//[16]
	4,//[17]
	5,//[18]
	7,//[19]
	10,//[20]
	10,//[21]
	2,//[22]
	2,//[23]
	4,//[24]
	4,//[25]
	8,//[26]
	3,//[27]
	5,//[28]
	5,//[29]
	6,//[30]
	6,//[31]
	7,//[32]
	7,//[33]
	8,//[34]
	9,//[35]
	5,//[36]
	6,//[37]
	2,//[38]
	1,//[39]
	2,//[40]
	3,//[41]
	};

	static const double t[]={0.00,
	0,//[1]
	0.25,//[2]
	1,//[3]
	2.75,//[4]
	4,//[5]
	0,//[6]
	0.25,//[7]
	0.75,//[8]
	2.75,//[9]
	0,//[10]
	2,//[11]
	0.75,//[12]
	3,//[13]
	3.5,//[14]
	1,//[15]
	2,//[16]
	4,//[17]
	3,//[18]
	0,//[19]
	0.5,//[20]
	1,//[21]
	1,//[22]
	7,//[23]
	5,//[24]
	6,//[25]
	6,//[26]
	10,//[27]
	13,//[28]
	14,//[29]
	11,//[30]
	14,//[31]
	8,//[32]
	14,//[33]
	6,//[34]
	7,//[35]
	24,//[36]
	22,//[37]
	3,//[38]
	1,//[39]
	0,//[40]
	0//[41]
	};

	static const double c[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-12]
	1,//[13]
	1,//[14]
	1,//[15]
	1,//[16]
	1,//[17]
	1,//[18]
	1,//[19]
	1,//[20]
	1,//[21]
	2,//[22]
	2,//[23]
	2,//[24]
	2,//[25]
	2,//[26]
	3,//[27]
	3,//[28]
	3,//[29]
	3,//[30]
	3,//[31]
	3,//[32]
	3,//[33]
	3,//[34]
	3,//[35]
	4,//[36]
	4,//[37]
	0,0,0,0 // indices [38-41]
	};

	// alpha is used here for consistency with the definitions in R744.c upon which Argon.c is based
	static const double alpha[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-37]
	20,
	20,
	20,
	20
	};

	static const double beta[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-37]
	250,
	375,
	300,
	225
	};

	static const double GAMMA[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-37]
	1.11,
	1.14,
	1.17,
	1.11
	};

	static const double epsilon[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-37]
	1,
	1,
	1,
	1
	};

	//Constants for ideal gas expression
	static double a0[]={0.0,
		8.31666243,
		-4.94651164
	};

	phirlist.push_back(new phir_power(n,d,t,c,1,37,38));
	phirlist.push_back(new phir_gaussian(n,d,t,alpha,epsilon,beta,GAMMA,38,41,42));

	phi0list.push_back(new phi0_lead(a0[1],a0[2]));
	phi0list.push_back(new phi0_logtau(1.5));

	// Critical parameters
	crit.rho = 535.6;
	crit.p = PressureUnit(4863.0,UNIT_KPA);
	crit.T = 150.687;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 39.948;
	params.Ttriple = 83.806;
	params.ptriple = 68.9004210852;
	params.accentricfactor = -0.00219;
	params.R_u = 8.31451;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 2000.0;
	limits.pmax = 1000000.0;
	limits.rhomax = 50.65*params.molemass;
	
	EOSReference.assign("\"A New Equation of State for Argon Covering the Fluid Region"
						" for Temperatures From the Melting Line to 700 K"
						" at Pressures up to 1000 MPa\""
						" Ch. Tegeler, R. Span, and W. Wagner,"
						" J. Phys. Chem. Ref. Data, Vol. 28, No. 3, 1999");
	TransportReference.assign("\"Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air\" E. W. Lemmon and R. T Jacobsen International Journal of Thermophysics, Vol. 25, No. 1, January 2004 \nNote: Critical enhancement included");

	name.assign("Argon");
	aliases.push_back("argon");
	aliases.push_back("ARGON");

	BibTeXKeys.EOS = "Tegeler-JPCRD-1999";
	BibTeXKeys.VISCOSITY = "Lemmon-IJT-2004";
	BibTeXKeys.CONDUCTIVITY = "Lemmon-IJT-2004";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double ArgonClass::X_tilde(double T,double tau,double delta)
{
	// X_tilde is dimensionless
	// Equation 11 slightly rewritten
	double drho_dp,R_Argon;
	R_Argon=params.R_u/params.molemass;
	drho_dp=1.0/(R_Argon*T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta)));
	return reduce.p.Pa*delta/reduce.rho*drho_dp;
}
double ArgonClass::conductivity_Trho(double T, double rho)
{
	double e_k=143.2, //[K]
		   sigma=0.335, //[nm]
		   Tref=301.374, //[K]
		   zeta0=0.13, //[nm]
		   LAMBDA=0.055,
		   q_D=0.32; //[nm]
	double eta0,OMEGA,delta,tau,Tstar,lambda0,lambdar,num,
		cp,cv,OMEGA_tilde,OMEGA_tilde0,zeta,nu,gamma,R0,lambdac,k,
		pi=3.141592654,mu;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,0.8158,-0.4320,0.0,13.73,10.07,0.7375,-33.96,20.47,-2.274,-3.973};
	double t[]={0,0,-0.77,-1.0,0.0,0.0,0.0,0.8,1.2,0.8,0.5};
	double d[]={0,0,0,0,1,2,4,5,6,9,1};
	double l[]={0,0,0,0,0,0,0,2,2,2,4};
	double g[]={0,0,0,0,0,0,0,1,1,1,1};
	
	delta=rho/reduce.rho;
	tau=reduce.T/T;
	Tstar=T/(e_k);

	OMEGA=exp(b[0]*powInt(log(Tstar),0)
			 +b[1]*powInt(log(Tstar),1)
		     +b[2]*powInt(log(Tstar),2)
			 +b[3]*powInt(log(Tstar),3)
		     +b[4]*powInt(log(Tstar),4));

	eta0=0.0266958*sqrt(params.molemass*T)/(sigma*sigma*OMEGA);
	lambda0=N[1]*eta0+N[2]*pow(tau,t[2])+N[3]*pow(tau,t[3]);

	lambdar=N[4]*pow(tau,t[4])*pow(delta,d[4])*exp(-g[4]*pow(delta,l[4]))
		   +N[5]*pow(tau,t[5])*pow(delta,d[5])*exp(-g[5]*pow(delta,l[5]))
		   +N[6]*pow(tau,t[6])*pow(delta,d[6])*exp(-g[6]*pow(delta,l[6]))
		   +N[7]*pow(tau,t[7])*pow(delta,d[7])*exp(-g[7]*pow(delta,l[7]))
	 	   +N[8]*pow(tau,t[8])*pow(delta,d[8])*exp(-g[8]*pow(delta,l[8]))
		   +N[9]*pow(tau,t[9])*pow(delta,d[9])*exp(-g[9]*pow(delta,l[9]))
		   +N[10]*pow(tau,t[10])*pow(delta,d[10])*exp(-g[10]*pow(delta,l[10]));

	R0=1.01;
	nu=0.63;
	gamma=1.2415;
	k=1.380658e-23; //[J/K]

	num=X_tilde(T,reduce.T/T,delta)-X_tilde(Tref,reduce.T/Tref,delta)*Tref/T;

	// no critical enhancement if numerator of Eq. 10 is negative
	if (num<0)
		return (lambda0+lambdar)/1e3;

	cp = specific_heat_p_Trho(T,rho); //[J/kg/K]
	cv = specific_heat_v_Trho(T,rho); //[J/kg/K]
	mu = viscosity_Trho(T,rho)*1e6; //[uPa-s]

	zeta=zeta0*pow(num/LAMBDA,nu/gamma); //[nm]
	OMEGA_tilde=2.0/pi*((cp-cv)/cp*atan(zeta/q_D)+cv/cp*(zeta/q_D));
	OMEGA_tilde0=2.0/pi*(1.-exp(-1./(q_D/zeta+1.0/3.0*(zeta/q_D)*(zeta/q_D)/delta/delta)));
	lambdac=rho*cp*k*R0*T/(6*pi*zeta*mu)*(OMEGA_tilde-OMEGA_tilde0)*1e18; // 1e18 is conversion to W/m-K (not described in paper)

	return (lambda0+lambdar+lambdac)/1e3; //[W/m/K]
}
double ArgonClass::viscosity_Trho(double T, double rho)
{
	double e_k=143.2, //[K]
		   sigma=0.335; //[nm]
	double eta0,etar,OMEGA,delta,tau,Tstar;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,12.19,13.99,0.005027,-18.93,-6.698,-3.827};
	double t[]={0,0.42,0.0,0.95,0.5,0.9,0.8};
	double d[]={0,1,2,10,5,1,2};
	double l[]={0,0,0,0,2,4,4};
	double g[]={0,0,0,0,1,1,1};

	delta=rho/reduce.rho;
	tau=reduce.T/T;
	Tstar=T/(e_k);
	OMEGA=exp(b[0]*powInt(log(Tstar),0)
			 +b[1]*powInt(log(Tstar),1)
		     +b[2]*powInt(log(Tstar),2)
			 +b[3]*powInt(log(Tstar),3)
		     +b[4]*powInt(log(Tstar),4));

	eta0=0.0266958*sqrt(params.molemass*T)/(sigma*sigma*OMEGA);
	etar=N[1]*pow(tau,t[1])*pow(delta,d[1])*exp(-g[1]*pow(delta,l[1]))
		+N[2]*pow(tau,t[2])*pow(delta,d[2])*exp(-g[2]*pow(delta,l[2]))
		+N[3]*pow(tau,t[3])*pow(delta,d[3])*exp(-g[3]*pow(delta,l[3]))
		+N[4]*pow(tau,t[4])*pow(delta,d[4])*exp(-g[4]*pow(delta,l[4]))
		+N[5]*pow(tau,t[5])*pow(delta,d[5])*exp(-g[5]*pow(delta,l[5]))
		+N[6]*pow(tau,t[6])*pow(delta,d[6])*exp(-g[6]*pow(delta,l[6]));

	return (eta0+etar)/1e6; // uPa-s to Pa-s
}
double ArgonClass::psat(double T)
{
	const double ti[]={0,1.0,1.5,2.0,4.5};
    const double ai[]={0,-5.9409785,1.3553888,-0.46497607,-1.5399043};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1-T/reduce.T,ti[i]);
    }
	return reduce.p.Pa*exp(reduce.T/T*summer);
}
double ArgonClass::rhosatL(double T)
{
	const double ti[]={0,0.334,2.0/3.0,7.0/3.0,4.0};
    const double ai[]={0,1.5004262,-0.31381290,0.086461622,-0.041477525};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/reduce.T,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double ArgonClass::rhosatV(double T)
{
	const double ti[]={0,0.345,5.0/6.0,1.0,13.0/3.0};
    const double ai[]={0,-1.70695656,-4.02739448,1.55177558,-2.30683228};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/reduce.T,ti[i]);
    }
    return reduce.rho*exp(reduce.T/T*summer);
}
double ArgonClass::surface_tension_T(double T)
{
	// From Mulero, 2012, JPCRD
	return 0.037*pow(1-T/reduce.T,1.25);
}
