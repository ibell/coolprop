/* Properties of Oxygen
by Ian Bell

Thermo properties from 
---------------------
"Thermodynamic Properties of Oxygen from the Triple point to 300 K with pressures to 80 MPa", 
Richard B. Stewart, Richard T. Jacobsen, and W. Wagner,
J. Phys. Chem. Ref. Data, v. 20, n. 5, 1991

Transport properties from
------------------------
"Viscosity and Thermal Conductivity Equations for
Nitrogen, Oxygen, Argon, and Air"
E. W. Lemmon and R. T Jacobsen
International Journal of Thermophysics, Vol. 25, No. 1, January 2004

Surface Tension
---------------
Lemmon, E.W. and Penoncello, S.G.,
"The Surface Tension of Air and Air Component Mixtures,"
Adv. Cryo. Eng., 39:1927-1934, 1994.
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
#include "Oxygen.h"



OxygenClass::OxygenClass()
{
	static const double n[]={0,    
	 0.3983768749,    //[1]
	-1.846157454,     //[2]
	 0.4183473197,    //[3]
	 0.2370620711e-1, //[4]

	 0.9771730573e-1, //[5]
	 0.3017891294e-1, //[6]
	 0.2273353212e-1, //[7]
	 0.1357254086e-1, //[8]

	-0.4052698943e-1, //[9]
	 0.5454628515e-3, //[10]
	 0.5113182277e-3, //[11]
	 0.2953466883e-6, //[12]

	-0.8687645072e-4, //[13]
	-0.2127082589,    //[14
	 0.8735941958e-1, //[15]
	 0.1275509190,    //[16]

	-0.9067701064e-1, //[17]
	-0.3540084206e-1, //[18]
	-0.3623278059e-1, //[19]
	 0.1327699290e-1, //[20]

	-0.3254111865e-3, //[21]
	-0.8313582932e-2, //[22]
	 0.2124570559e-2, //[23]
	-0.8325206232e-3, //[24]

	-0.2626173276e-4, //[25]
	 0.2599581482e-2, //[26]
	 0.9984649663e-2, //[27]
	 0.2199923153e-2, //[28]

	-0.2591350486e-1, //[29]
	-0.1259630848,    //[30]
	 0.1478355637,    //[31]
	-0.1011251078e-1  //[32]

	};

	// d used for consistency (corresponds to i from Stewart)
	static const double d[]={0,
	1,//[1]
	1,//[2]
	1,//[3]
	2,//[4]

	2,//[5]
	2,//[6]
	3,//[7]
	3,//[8]

	3,//[9]
	6,//[10]
	7,//[11]
	7,//[12]

	8,//[13]
	1,//[14]
	1,//[15]
	2,//[16]

	2,//[17]
	3,//[18]
	3,//[19]
	5,//[20]

	6,//[21]
	7,//[22]
	8,//[23]
	10,//[24]

	2,//[25]
	3,//[26]
	3,//[27]
	4,//[28]

	4,//[29]
	5,//[30]
	5,//[31]
	5,//[32]
	};

	// t used for consistency (corresponds to j from Stewart)
	static const double t[]={0.00,
	 0.0,//[1]
	 1.5,//[2]
	 2.5,//[3]
	-0.5,//[4]

	 1.5,//[5]
	 2.0,//[6]
	 0.0,//[7]
	 1.0,//[8]

	 2.5,//[9]
	 0.0,//[10]
	 2.0,//[11]
	 5.0,//[12]

	 2.0,//[13]
	 5.0,//[14]
	 6.0,//[15]
	 3.5,//[16]

	 5.5,//[17]
	 3.0,//[18]
	 7.0,//[19]
	 6.0,//[20]

	 8.5,//[21]
	 4.0,//[22]
	 6.5,//[23]
	 5.5,//[24]

	 22.0,//[25]
	 11.0,//[26]
	 18.0,//[27]
	 11.0,//[28]

	 23.0,//[29]
	 17.0,//[30]
	 18.0,//[31]
	 23.0,//[32]
	};

	// c used for consistency (corresponds to l from Stewart)
	static const double cv[]={0,
	0,//[1]
	0,//[2]
	0,//[3]
	0,//[4]
	0,//[5]
	0,//[6]
	0,//[7]
	0,//[8]
	0,//[9]
	0,//[10]
	0,//[11]
	0,//[12]
	0,//[13]
	2,//[14]
	2,//[15]
	2,//[16]
	2,//[17]
	2,//[18]
	2,//[19]
	2,//[20]
	2,//[21]
	2,//[22]
	2,//[23]
	2,//[24]
	4,//[25]
	4,//[26]
	4,//[27]
	4,//[28]
	4,//[29]
	4,//[30]
	4,//[31]
	4 //[32]
	};

	//Constants for ideal gas expression
	static const double N0[]={0.0,
		1.06778,
		3.50042,
		0.166961e-7,
		1.01258,
		0.944365,
		2242.45,
		11580.4,
	};

	// The residual HE terms
	phirlist.push_back(new phir_power(n,d,t,cv,1,32,33));

	// Critical parameters
	crit.rho = 13.63*31.9988;
	crit.p = PressureUnit(5043.0, UNIT_KPA);
	crit.T = 154.581;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 31.9988;
	params.Ttriple = 54.361;
	params.ptriple = 0.146323903868;
	params.accentricfactor = 0.0222;
	params.R_u = 8.31434;

	double T0 = 298.15, 
		   p0 = 101.325, 
		   R_ = 8.31434/params.molemass,
		   rho0 = p0/(R_*T0),
		   m,
		   c,
		   R_u = 8.31434,
		   H0 = 8680.0, /// kJ/kmol
		   S0 = 205.043, /// kJ/kmol/K
		   Tc = crit.T,
		   tau0 = crit.T/T0, 
		   delta0 = rho0/crit.rho;
	
	// log(delta)+c+m*tau
	
	/// c is the constant term
	c=-S0/R_u-1+log(tau0/delta0);/*<< from the leading term*/

	/// m multiplies the tau term in the leading term (slope)
	m=H0/(R_u*Tc); /*<< from the leading term */

	std::vector<double> N0_v(N0,N0+sizeof(N0)/sizeof(double));

	phi_BC * phi0_lead_ = new phi0_lead(c,m);
	phi_BC * phi0_logtau_ = new phi0_logtau(-1);
	
	phi_BC * phi0_cp0_poly_1 = new phi0_cp0_poly(N0_v[1],-1.5,Tc,T0);/// checked - good
	phi_BC * phi0_cp0_constant_ = new phi0_cp0_constant(N0_v[2],Tc,T0);/// checked - good
	phi_BC * phi0_cp0_poly_2 = new phi0_cp0_poly(N0_v[3],2,Tc,T0);/// checked - good
	
	// The next term turns into one of the first form of the phi0_exponential
	// Term is of the form a_0*log(1-exp(-theta_0*tau))
	// theta_0 is N[6]/Tc
	phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(N0_v[4],N0_v[6]/Tc);/// checked - good

	// The last term turns into one of the second form of the phi0_exponential
	// Term is of the form a_0*log(c+exp(theta_0*tau))
	// c = 2/3
	// theta_0 is N[7]/Tc
	phi_BC * phi0_Planck_Einstein2_ = new phi0_Planck_Einstein2(N0_v[5],N0_v[7]/Tc,2.0/3.0);/// checked - good

	phi0list.push_back(phi0_lead_);
	phi0list.push_back(phi0_logtau_);
	phi0list.push_back(phi0_cp0_poly_1);
	phi0list.push_back(phi0_cp0_constant_);
	phi0list.push_back(phi0_cp0_poly_2);

	phi0list.push_back(phi0_Planck_Einstein_);
	phi0list.push_back(phi0_Planck_Einstein2_);

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 2000.0;
	limits.pmax = 2200000.0;
	limits.rhomax = 53.15*params.molemass;
	
	EOSReference.assign("\"Thermodynamic Properties of Oxygen from the Triple point to 300 K with pressures to 80 MPa\", "
						"Richard B. Stewart, Richard T. Jacobsen, and W. Wagner, "
						"J. Phys. Chem. Ref. Data, v. 20, n. 5, 1991");
	TransportReference.assign("Viscosity and Thermal Conductivity: \"Viscosity and Thermal Conductivity Equations for"
							  "Nitrogen, Oxygen, Argon, and Air\""
							  "E. W. Lemmon and R. T Jacobsen"
							  "International Journal of Thermophysics, Vol. 25, No. 1, January 2004");

	name.assign("Oxygen");
	aliases.push_back("oxygen");
	aliases.push_back(std::string("OXYGEN"));
	aliases.push_back("O2");

	BibTeXKeys.EOS = "Stewart-JPCRD-1991";
	BibTeXKeys.VISCOSITY = "Lemmon-IJT-2004";
	BibTeXKeys.CONDUCTIVITY = "Lemmon-IJT-2004";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double OxygenClass::X_tilde(double T,double tau,double delta)
{
	// X_tilde is dimensionless
	// Equation 11 slightly rewritten
	double drho_dp,R_Oxygen;
	R_Oxygen=params.R_u/params.molemass;
	drho_dp=1.0/(R_Oxygen*T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta)));
	return reduce.p.Pa*delta/reduce.rho*drho_dp;
}

double OxygenClass::conductivity_Trho(double T, double rho)
{
	double e_k=118.5, //[K]
		   sigma=0.3428, //[nm]
		   Tref=309.162, //[K]
		   zeta0=0.24, //[nm]
		   LAMBDA=0.055,
		   q_D=0.51; //[nm]
	double eta0,OMEGA,delta,tau,Tstar,lambda0,lambdar,num,
		cp,cv,OMEGA_tilde,OMEGA_tilde0,zeta,nu,gamma,R0,lambdac,k,
		pi=3.141592654,mu;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,1.036,6.283,-4.262,15.31,8.898,-0.7336,6.728,-4.374,-0.4747};
	double t[]={0,0,-0.9,-0.6,0.0,0.0,0.3,4.3,0.5,1.8};
	double d[]={0,0,0,0,1,3,4,5,7,10};
	double l[]={0,0,0,0,0,0,0,2,2,2};
	double g[]={0,0,0,0,0,0,0,1,1,1};
	
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
		   +N[9]*pow(tau,t[9])*pow(delta,d[9])*exp(-g[9]*pow(delta,l[9]));

	R0=1.01;
	nu=0.63;
	gamma=1.2415;
	k=1.380658e-23; //[J/K]

	num=X_tilde(T,reduce.T/T,delta)-X_tilde(Tref,reduce.T/Tref,delta)*Tref/T;

	// no critical enhancement if numerator of Eq. 10 is negative
	if (num<0)
		return (lambda0+lambdar)/1e6;

	cp = specific_heat_p_Trho(T,rho); //[J/kg/K]
	cv = specific_heat_v_Trho(T,rho); //[J/kg/K]
	mu = viscosity_Trho(T,rho)*1e6; //[uPa-s]

	zeta=zeta0*pow(num/LAMBDA,nu/gamma); //[nm]
	OMEGA_tilde=2.0/pi*((cp-cv)/cp*atan(zeta/q_D)+cv/cp*(zeta/q_D));
	OMEGA_tilde0=2.0/pi*(1.-exp(-1./(q_D/zeta+1.0/3.0*(zeta/q_D)*(zeta/q_D)/delta/delta)));
	lambdac=rho*(cp*1000.0)*k*R0*T/(6*pi*zeta*mu)*(OMEGA_tilde-OMEGA_tilde0)*1e18; // 1e18 is conversion to mW/m-K (not described in paper)

	return (lambda0+lambdar+lambdac)/1e3; //[W/m/K]
}
double OxygenClass::viscosity_Trho(double T, double rho)
{
	double e_k=118.5, //[K]
		   sigma=0.3428; //[nm]
	double eta0,etar,OMEGA,delta,tau,Tstar;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,17.67,0.4042,0.0001077,0.3510,-13.67};
	double t[]={0,0.05,0.0,2.10,0.0,0.5};
	double d[]={0,1,5,12,8,1};
	double l[]={0,0,0,0,1,2};
	double g[]={0,0,0,0,1,1};

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
		+N[5]*pow(tau,t[5])*pow(delta,d[5])*exp(-g[5]*pow(delta,l[5]));

	return (eta0+etar)/1e6; // uPa-s to Pa-s
}
double OxygenClass::psat(double T)
{
	const double ti[]={0,1.0,3.0/2.0,3.0,7.0,9.0};
    const double Ni[]={0,-6.043938,1.175627,-0.994086,-3.456781,3.361499};
    double summer=0;
    int i;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(1-T/reduce.T,ti[i]);
    }
	double p = reduce.p.Pa*exp(reduce.T/T*summer);
	return p;
}
double OxygenClass::rhosatV(double T)
{
	const double ti[]={0,1.0/3.0,2.0/3.0,1.0,5.0/3.0,4.0,9.0};
    const double Ni[]={0,-1.498431,-2.116826,-0.905713,-5.659990,-18.90964,-53.780774};
    double summer=0;
    int i;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(1.0-T/reduce.T,ti[i]);
    }
	double rho = reduce.rho*exp(summer);
	return rho;
}
double OxygenClass::rhosatL(double T)
{
	const double ti[]={0,1.0/3.0,2.0/3.0,3.0};
    const double Ni[]={0,1.507678,0.85810805,0.19035504};
    double summer=1;
    int i;
    for (i=1;i<=3;i++)
    {
        summer += Ni[i]*pow(1.0-T/reduce.T,ti[i]);
    }
    double rho = reduce.rho*summer;
	return rho;
}
double OxygenClass::surface_tension_T(double T)
{
	// From Mulero, 2012, JPCRD
	return 0.03843*pow(1-T/reduce.T,1.225);
}
