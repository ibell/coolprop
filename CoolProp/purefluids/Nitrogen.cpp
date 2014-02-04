/* Properties of Nitrogen
by Ian Bell

Thermo properties from 
---------------------
"A Reference Equation of State for the Thermodynamic Properties
of Nitrogen for Temperatures from 63.151 to 1000 K 
and Pressures to 2200 MPa", 
R. Span and E.W. Lemmon and R.T. Jacobsen and W. Wagner and A. Yokozeki, 
J. Phys. Chem. Ref. Data, v. 29, n. 6, 2000

Transport properties from
------------------------
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
#else
#include <stdlib.h>
#endif

#include "math.h"
#include "stdio.h"
#include <string.h>
#include "CoolProp.h"
#include "FluidClass.h"
#include "Nitrogen.h"




NitrogenClass::NitrogenClass()
{

	static const double n[]={0,    
	0.924803575275,//[1]
	-0.492448489428,//[2]
	0.661883336938,//[3]
	-1.92902649201,//[4]
	-0.0622469309629,//[5]
	0.349943957581,//[6]
	0.564857472498,//[7]
	-1.61720005987,//[8]
	-0.481395031883,//[9]
	0.421150636384,//[10]
	-0.0161962230825,//[11]
	0.172100994165,//[12]
	0.00735448924933,//[13]
	0.0168077305479,//[14]
	-0.00107626664179,//[15]
	-0.0137318088513,//[16]
	0.000635466899859,//[17]
	0.00304432279419,//[18]
	-0.0435762336045,//[19]
	-0.0723174889316,//[20]
	0.0389644315272,//[21]
	-0.021220136391,//[22]
	0.00408822981509,//[23]
	-0.0000551990017984,//[24]
	-0.0462016716479,//[25]
	-0.00300311716011,//[26]
	0.0368825891208,//[27]
	-0.0025585684622,//[28]
	0.00896915264558,//[29]
	-0.0044151337035,//[30]
	0.00133722924858,//[31]
	0.000264832491957,//[32]
	19.6688194015,//[33]
	-20.911560073,//[34]
	0.0167788306989,//[35]
	2627.67566274//[36]
	};

	// d used for consistency with CO2 correlation (corresponds to i from Span)
	static const double d[]={0,
	1,//[1]
	1,//[2]
	2,//[3]
	2,//[4]
	3,//[5]
	3,//[6]
	1,//[7]
	1,//[8]
	1,//[9]
	3,//[10]
	3,//[11]
	4,//[12]
	6,//[13]
	6,//[14]
	7,//[15]
	7,//[16]
	8,//[17]
	8,//[18]
	1,//[19]
	2,//[20]
	3,//[21]
	4,//[22]
	5,//[23]
	8,//[24]
	4,//[25]
	5,//[26]
	5,//[27]
	8,//[28]
	3,//[29]
	5,//[30]
	6,//[31]
	9,//[32]
	1,//[33]
	1,//[34]
	3,//[35]
	2//[36]
	};

	// t used for consistency with CO2 correlation (corresponds to j from Span)
	static const double t[]={0.00,
	0.25,//[1]
	0.875,//[2]
	0.5,//[3]
	0.875,//[4]
	0.375,//[5]
	0.75,//[6]
	0.5,//[7]
	0.75,//[8]
	2,//[9]
	1.25,//[10]
	3.5,//[11]
	1,//[12]
	0.5,//[13]
	3,//[14]
	0,//[15]
	2.75,//[16]
	0.75,//[17]
	2.5,//[18]
	4,//[19]
	6,//[20]
	6,//[21]
	3,//[22]
	3,//[23]
	6,//[24]
	16,//[25]
	11,//[26]
	15,//[27]
	12,//[28]
	12,//[29]
	7,//[30]
	4,//[31]
	16,//[32]
	0,//[33]
	1,//[34]
	2,//[35]
	3//[36]
	};

	// c used for consistency with CO2 correlation (corresponds to l from Span)
	static const double c[]={0,
	0,//[1]
	0,//[2]
	0,//[3]
	0,//[4]
	0,//[5]
	0,//[6]
	1,//[7]
	1,//[8]
	1,//[9]
	1,//[10]
	1,//[11]
	1,//[12]
	1,//[13]
	1,//[14]
	1,//[15]
	1,//[16]
	1,//[17]
	1,//[18]
	2,//[19]
	2,//[20]
	2,//[21]
	2,//[22]
	2,//[23]
	2,//[24]
	3,//[25]
	3,//[26]
	3,//[27]
	3,//[28]
	4,//[29]
	4,//[30]
	4,//[31]
	4,//[32]
	2,//[33]
	2,//[34]
	2,//[35]
	2//[36]
	};

	// alpha is used here for consistency with the definitions in R744.c upon which Nitrogen.c is based
	// is phi_k from Span
	static const double alpha[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-32]
	20,
	20,
	15,
	25
	};

	static const double beta[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-32]
	325,
	325,
	300,
	275
	};

	// epsilon is used here for consistency with the definitions in R744.c upon which Nitrogen.c is based
	// is the value unity in Span
	static const double epsilon[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-32]
	1,
	1,
	1,
	1
	};


	// GAMMA is used here for consistency with the definitions in R744.c upon which Nitrogen.c is based
	static const double GAMMA[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-32]
	1.16,
	1.16,
	1.13,
	1.25
	};

	//Constants for ideal gas expression
	static const double a0[]={0.0,
		2.5,
		-12.76952708,
		-0.00784163,
		-1.934819e-4,
		-1.247742e-5,
		6.678326e-8,
		1.012941,
		26.65788
	};
	

	phirlist.push_back(new phir_power(n,d,t,c,1,32,33));
	phirlist.push_back(new phir_gaussian(n,d,t,alpha,epsilon,beta,GAMMA,33,36,37));

	// phi0=log(delta)+a0[1]*log(tau)+a0[2]+a0[3]*tau+a0[4]/tau+a0[5]/tau/tau+a0[6]/tau/tau/tau+a0[7]*log(1-exp(-a0[8]*tau));
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> n0_v(sizeof(a0)/sizeof(double),0);
	n0_v[4]=-1.0;
	n0_v[5]=-2.0;
	n0_v[6]=-3.0;
	n0_v[7]=a0[8];
	phi_BC * phi0_lead_ = new phi0_lead(a0[2],a0[3]);
	phi_BC * phi0_logtau_ = new phi0_logtau(a0[1]);
	phi_BC * phi0_power_ = new phi0_power(a0_v,n0_v,4,6);
	phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(a0_v,n0_v,7,7);

	phi0list.push_back(phi0_lead_);
	phi0list.push_back(phi0_logtau_);
	phi0list.push_back(phi0_power_);
	phi0list.push_back(phi0_Planck_Einstein_);

	// Critical parameters
	crit.rho = 313.3;
	crit.p = PressureUnit(3395.8, UNIT_KPA);
	crit.T = 126.192;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 28.01348;
	params.Ttriple = 63.151;
	params.ptriple = 12.5220865181;
	params.accentricfactor = 0.0372 ;
	params.R_u = 8.31451;

	// Limits of EOS
	limits.Tmin = 63.151;
	limits.Tmax = 2000.0;
	limits.pmax = 2200000.0;
	limits.rhomax = 53.15*params.molemass;
	
	EOSReference.assign("\"A Reference Equation of State for the Thermodynamic Properties"
						"of Nitrogen for Temperatures from 63.151 to 1000 K "
						"and Pressures to 2200 MPa\", " 
						"R. Span and E.W. Lemmon and R.T. Jacobsen and W. Wagner and A. Yokozeki, "
						"J. Phys. Chem. Ref. Data, v. 29, n. 6, 2000");
	TransportReference.assign("\"Viscosity and Thermal Conductivity Equations for"
							  "Nitrogen, Oxygen, Argon, and Air\""
							  "E. W. Lemmon and R. T Jacobsen"
							  "International Journal of Thermophysics, Vol. 25, No. 1, January 2004");

	name.assign("Nitrogen");
	aliases.push_back("nitrogen");
	aliases.push_back(std::string("NITROGEN"));
	aliases.push_back("N2");

	BibTeXKeys.EOS = "Span-JPCRD-2000";
	BibTeXKeys.VISCOSITY = "Lemmon-IJT-2004";
	BibTeXKeys.CONDUCTIVITY = "Lemmon-IJT-2004";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";
}

double NitrogenClass::X_tilde(double T,double tau,double delta)
{
	// X_tilde is dimensionless
	// Equation 11 slightly rewritten
	double drho_dp,R_Nitrogen;
	R_Nitrogen=params.R_u/params.molemass;
	drho_dp=1.0/(R_Nitrogen*T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta)));
	return reduce.p.Pa*delta/reduce.rho*drho_dp;
}

double NitrogenClass::conductivity_dilute(double T)
{
	double e_k=98.94, //[K]
		   sigma=0.3656, //[nm]
		   Tstar, OMEGA, eta0, lambda0, tau;
	
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,1.511,2.117,-3.332,8.862,31.11,-73.13,20.03,-0.7096,0.2672};
	double t[]={0,0,-1.0,-0.7,0.0,0.03,0.2,0.8,0.6,1.9};
	
	tau=reduce.T/T;
	Tstar=T/(e_k);

	OMEGA=exp(b[0]*powInt(log(Tstar),0)
			 +b[1]*powInt(log(Tstar),1)
		     +b[2]*powInt(log(Tstar),2)
			 +b[3]*powInt(log(Tstar),3)
		     +b[4]*powInt(log(Tstar),4));

	eta0=0.0266958*sqrt(params.molemass*T)/(sigma*sigma*OMEGA);
	lambda0=N[1]*eta0+N[2]*pow(tau,t[2])+N[3]*pow(tau,t[3]);
	return lambda0/1e3; //[W/m/K]
}
double NitrogenClass::conductivity_background(double T, double rho)
{
	double tau, delta, lambdar;
	double N[]={0,1.511,2.117,-3.332,8.862,31.11,-73.13,20.03,-0.7096,0.2672};
	double t[]={0,0,-1.0,-0.7,0.0,0.03,0.2,0.8,0.6,1.9};
	double d[]={0,0,0,0,1,2,3,4,8,10};
	double l[]={0,0,0,0,0,0,1,2,2,2};
	double g[]={0,0,0,0,0,0,1,1,1,1};
	
	delta=rho/reduce.rho;
	tau=reduce.T/T;

	lambdar=N[4]*pow(tau,t[4])*pow(delta,d[4])*exp(-g[4]*pow(delta,l[4]))
		   +N[5]*pow(tau,t[5])*pow(delta,d[5])*exp(-g[5]*pow(delta,l[5]))
		   +N[6]*pow(tau,t[6])*pow(delta,d[6])*exp(-g[6]*pow(delta,l[6]))
		   +N[7]*pow(tau,t[7])*pow(delta,d[7])*exp(-g[7]*pow(delta,l[7]))
	 	   +N[8]*pow(tau,t[8])*pow(delta,d[8])*exp(-g[8]*pow(delta,l[8]))
		   +N[9]*pow(tau,t[9])*pow(delta,d[9])*exp(-g[9]*pow(delta,l[9]));
	return lambdar/1e3; // [W/m/K]
}

double NitrogenClass::conductivity_critical(double T, double rho)
{
	double num,delta,
		cp,cv,OMEGA_tilde,OMEGA_tilde0,zeta,nu,gamma,R0,lambdac,k,
		pi=3.141592654,mu,
		Tref=252.384, //[K]
	   zeta0=0.17, //[nm]
	   LAMBDA=0.055,
	   q_D=0.40; //[nm];
	
	R0=1.01;
	nu=0.63;
	gamma=1.2415;
	k=1.380658e-23; //[J/K]

	delta=rho/reduce.rho;

	num=X_tilde(T,reduce.T/T,delta)-X_tilde(Tref,reduce.T/Tref,delta)*Tref/T;

	// no critical enhancement if numerator of Eq. 10 is negative
	if (num<0)
		return 0;

	cp = specific_heat_p_Trho(T,rho); //[J/kg/K]
	cv = specific_heat_v_Trho(T,rho); //[J/kg/K]
	mu = viscosity_Trho(T,rho)*1e6; //[uPa-s]

	zeta = zeta0*pow(num/LAMBDA,nu/gamma); //[nm]
	OMEGA_tilde = 2.0/pi*((cp-cv)/cp*atan(zeta/q_D)+cv/cp*(zeta/q_D));
	OMEGA_tilde0 = 2.0/pi*(1.-exp(-1./(q_D/zeta+1.0/3.0*(zeta/q_D)*(zeta/q_D)/delta/delta)));
	lambdac = rho*cp*k*R0*T/(6*pi*zeta*mu)*(OMEGA_tilde-OMEGA_tilde0)*1e18; // 1e18 is conversion to mW/m-K (not described in paper)

	return lambdac/1e3; //[W/m/K]
}
double NitrogenClass::conductivity_Trho(double T, double rho)
{
	// Each term in kW/m/K
	return conductivity_dilute(T) + conductivity_background(T,rho) + conductivity_critical(T,rho);
}
double NitrogenClass::viscosity_dilute(double T)
{
	double e_k=98.94, //[K]
		   sigma=0.3656; //[nm]
	double eta0,OMEGA,Tstar;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};
	
	Tstar=T/(e_k);
	OMEGA=exp(b[0]*powInt(log(Tstar),0)
			 +b[1]*powInt(log(Tstar),1)
		     +b[2]*powInt(log(Tstar),2)
			 +b[3]*powInt(log(Tstar),3)
		     +b[4]*powInt(log(Tstar),4));

	eta0=0.0266958*sqrt(params.molemass*T)/(sigma*sigma*OMEGA); //[Pa-s]
	return eta0/1e6; //[Pa-s]
}
double NitrogenClass::viscosity_residual(double T, double rho)
{
	return viscosity_background(T,rho);
}
double NitrogenClass::viscosity_background(double T, double rho)
{
	double etar, tau, delta;
	delta=rho/reduce.rho;
	tau=reduce.T/T;
	double N[]={0,10.72,0.03989,0.001208,-7.402,4.620};
	double g[]={0,0,1,1,1,1};
	double t[]={0,0.1,0.25,3.2,0.9,0.3};
	double d[]={0,2,10,12,2,1};
	double l[]={0,0,1,1,2,3};

	etar=N[1]*pow(tau,t[1])*pow(delta,d[1])*exp(-g[1]*pow(delta,l[1]))
		+N[2]*pow(tau,t[2])*pow(delta,d[2])*exp(-g[2]*pow(delta,l[2]))
		+N[3]*pow(tau,t[3])*pow(delta,d[3])*exp(-g[3]*pow(delta,l[3]))
		+N[4]*pow(tau,t[4])*pow(delta,d[4])*exp(-g[4]*pow(delta,l[4]))
		+N[5]*pow(tau,t[5])*pow(delta,d[5])*exp(-g[5]*pow(delta,l[5]));

	return etar/1e6; // uPa-s to Pa-s
}
double NitrogenClass::viscosity_Trho(double T, double rho)
{
	return this->viscosity_dilute(T)+this->viscosity_background(T, rho);
}
double NitrogenClass::psat(double T)
{
	const double ti[]={0,1.0,1.5,2.5,5.0};
    const double Ni[]={0,-6.12445284,1.26327220,-0.765910082,-1.77570564};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+Ni[i]*pow(1-T/reduce.T,ti[i]);
    }
	return reduce.p.Pa*exp(reduce.T/T*summer);
}
double NitrogenClass::rhosatL(double T)
{
	const double ti[]={0,0.3294,2.0/3.0,8.0/3.0,35.0/6.0};
    const double Ni[]={0,1.48654237,-0.280476066,0.0894143085,-0.119879866};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+Ni[i]*pow(1.0-T/reduce.T,ti[i]);
    }
	return reduce.rho*exp(summer);
}
double NitrogenClass::rhosatV(double T)
{
	const double ti[]={0,0.34,5.0/6.0,7.0/6.0,13.0/6.0,14.0/3.0};
    const double Ni[]={0,-1.70127164,-3.70402649,1.29859383,-0.561424977,-2.68505381};
    double summer=0;
    int i;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(1.0-T/reduce.T,ti[i]);
    }
    return reduce.rho*exp(reduce.T/T*summer);
}
double NitrogenClass::surface_tension_T(double T)
{
	// From Mulero, 2012, JPCRD
	return 0.02898*pow(1-T/reduce.T,1.246);
}
