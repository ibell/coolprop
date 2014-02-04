/* Properties of Air
by Ian Bell
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
#include "Air.h"

AirClass::AirClass()
{
	static const double N[]={0,
	 0.118160747229,//[1]
	 0.713116392079,//[2]
	-0.161824192067e1,//[3]
	 0.714140178971e-1,//[4]
	-0.865421396646e-1,//[5]
	 0.134211176704,//[6]
	 0.112626704218e-1,//[7]
	-0.420533228842e-1,//[8]
	 0.349008431982e-1,//[9]
	 0.164957183186e-3,//[10]
	-0.101365037912,//[11]
	-0.173813690970,//[12]
	-0.472103183731e-1,//[13]
	-0.122523554253e-1,//[14]
	-0.146629609713,//[15]
	-0.316055879821e-1,//[16]
	 0.233594806142e-3,//[17]
	 0.148287891978e-1,//[18]
	-0.938782884667e-2//[19]
	};

	static const double d[]={0,
	1,//[1]
	1,//[2]
	1,//[3]
	2,//[4]
	3,//[5]
	3,//[6]
	4,//[7]
	4,//[8]
	4,//[9]
	6,//[10]
	1,//[11]
	3,//[12]
	5,//[13]
	6,//[14]
	1,//[15]
	3,//[16]
	11,//[17]
	1,//[18]
	3//[19]
	};

	static const double t[]={0.00,
	0,//[1]
	0.33,//[2]
	1.01,//[3]
	0,//[4]
	0,//[5]
	0.15,//[6]
	0,//[7]
	0.2,//[8]
	0.35,//[9]
	1.35,//[10]
	1.6,//[11]
	0.8,//[12]
	0.95,//[13]
	1.25,//[14]
	3.6,//[15]
	6,//[16]
	3.25,//[17]
	3.5,//[18]
	15//[19]
	};

	static const double l[]={
	0,0,0,0,0,0,0,0,0,0,0, // indices [0-10]
	1,//[11]
	1,//[12]
	1,//[13]
	1,//[14]
	2,//[15]
	2,//[16]
	2,//[17]
	3,//[18]
	3,//[19]
	};

	//Constants for ideal gas expression
	static const double N0[]={0.0,
	 0.605719400e-7,//[1]
	-0.210274769e-4,//[2]
	-0.158860716e-3,//[3]
	-13.841928076,//[4]
	 17.275266575,//[5]
	-0.195363420e-3,//[6]
	 2.490888032,//[7]
	 0.791309509,//[8]
	 0.212236768,//[9]
	-0.197938904,//[10]
	 25.36365,//[11]
	 16.90741,//[12]
	 87.31279//[13]
	};

	phirlist.push_back(new phir_power(N,d,t,l,1,19,20));
	
	/*
	phi0=log(delta);
    for (k=1;k<=5;k++)
    {
        phi0=phi0+N0[k]*powInt(tau,k-4);
    }
    phi0+=N0[6]*pow(tau,1.5)+N0[7]*log(tau);
    for (k=8;k<=9;k++)
    {
        phi0+=N0[k]*log(1.0-exp(-N0[k+3]*tau));
    }
    phi0+=N0[10]*log(2.0/3.0+exp(N0[13]*tau));
	*/
	std::vector<double> N0_v(N0,N0+sizeof(N0)/sizeof(double));
	std::vector<double> theta0_v(N0,N0+sizeof(N0)/sizeof(double));
	// Set some constants
	theta0_v[1]=1-4;
	theta0_v[2]=2-4;
	theta0_v[3]=3-4;
	theta0_v[4]=4-4;
	theta0_v[5]=5-4;
	theta0_v[6]=1.5;
	theta0_v[8]=N0_v[11];
	theta0_v[9]=N0_v[12];
	theta0_v[10]=N0_v[13];

	phi0list.push_back(new phi0_lead(0.0,0.0));
	phi0list.push_back(new phi0_power(N0_v,theta0_v,1,5));
	phi0list.push_back(new phi0_power(N0_v[6],1.5));
	phi0list.push_back(new phi0_logtau(N0_v[7]));
	phi0list.push_back(new phi0_Planck_Einstein(N0_v,theta0_v,8,9));
	phi0list.push_back(new phi0_Planck_Einstein2(N0_v[10],theta0_v[10],2.0/3.0));
	
	// Critical parameters (max condensing temperature)
	crit.rho = 11.8308*28.96546;
	crit.p = PressureUnit(3786.0,UNIT_KPA);
	crit.T = 132.5306;
	crit.v = 1.0/crit.rho;

	maxcondT.rho = 10.4477*28.96546;
	maxcondT.p = PressureUnit(3785.02,UNIT_KPA);
	maxcondT.T = 132.6312;
	maxcondT.v = 1.0/maxcondT.rho;

	// Other fluid parameters
	params.molemass = 28.96546;
	params.Ttriple = 59.75;
	params.ptriple = 2.44582352329;
	params.accentricfactor = 0.0816749632704;
	params.R_u = 8.31451;
	isPure = false;
	preduce = &maxcondT;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 2000.0;
	limits.pmax = 2000000.0;
	limits.rhomax = 14.21*params.molemass;
	
	EOSReference.assign("Lemmon, E.W., Jacobsen, R.T, Penoncello, S.G., and Friend, D.G.,"
						"\"Thermodynamic Properties of Air and Mixtures of Nitrogen, Argon, and"
						" Oxygen from 60 to 2000 K at Pressures to 2000 MPa,"
						" J. Phys. Chem. Ref. Data, 29(3):331-385, 2000.");
	TransportReference.assign("E.W. Lemmon and R. T. Jacobsen, \"Viscosity and Thermal Conductivity Equations for "
							  "Nitrogen, Oxygen, Argon, and Air\", "
							  "International Journal of Thermophysics, Vol. 25, No. 1, January 2004");

	name.assign("Air");
	aliases.push_back("air");
	aliases.push_back(std::string("AIR"));

	BibTeXKeys.EOS = "Lemmon-JPCRD-2000";
	BibTeXKeys.VISCOSITY = "Lemmon-IJT-2004";
	BibTeXKeys.CONDUCTIVITY = "Lemmon-IJT-2004";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double AirClass::rhosatL(double T)
{
	const double ti[]={0,0.65,0.85,0.95,1.1};
    const double Ni[]={0,44.3413,-240.073,285.139,-88.3366,-0.892181};
    double summer=0; int k;
    for (k=1;k<=4;k++)
    {
        summer=summer+Ni[k]*pow(1.0-T/reduce.T,ti[k]);
    }
    return reduce.rho*(1+summer+Ni[5]*log(T/reduce.T));
}

double AirClass::rhosatV(double T)
{
	const double ti[]={0,0.41,1.0,2.8,6.5};
    const double Ni[]={0,-2.0466,-4.7520,-13.259,-47.652};
    double summer=0; int k;
    for (k=1;k<=4;k++)
    {
        summer=summer+Ni[k]*pow(1.0-T/reduce.T,ti[k]);
    }
    return reduce.rho*exp(summer);
}

double AirClass::psatL(double T)
{
	const double Ni[]={0,0.2260724,-7.080499,5.700283,-12.44017,17.81926,-10.81364};
	double summer=0; int k;
    for (k=1;k<=6;k++)
    {
        summer=summer+Ni[k]*pow(1-T/reduce.T,((double)k)/2.0);
    }
	double p = reduce.p.Pa*exp(reduce.T/T*summer);
	return p;
}

double AirClass::psatV(double T)
{
	const double Ni[]={0,-0.1567266,-5.539635,0,0,0.7567212,0,0,-3.514322};
	double summer=0; int k;
    for (k=1;k<=8;k++)
    {
        summer=summer+Ni[k]*pow(1-T/reduce.T,((double)k)/2.0);
    }
	double p =  reduce.p.Pa*exp(reduce.T/T*summer);
	return p;
}

double AirClass::viscosity_Trho(double T, double rho)
{
	/*
	E.W. Lemmon and R.T Jacobsen, Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon and Air
	International Journal of Thermophysics, Vol. 25, No. 1, January 2004, p.28 
	*/

	double e_k=103.3, //[K]
		   sigma=0.360; //[nm]
	double eta0,etar,OMEGA,delta,tau,Tstar;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};
	
	// Ideal-gas part
	Tstar = T/(e_k);
	OMEGA=exp(b[0]*powInt(log(Tstar),0)
			 +b[1]*powInt(log(Tstar),1)
		     +b[2]*powInt(log(Tstar),2)
			 +b[3]*powInt(log(Tstar),3)
		     +b[4]*powInt(log(Tstar),4));

	eta0=0.0266958*sqrt(28.9586*T)/(sigma*sigma*OMEGA);

	// Residual part
	double N[]={0,10.72,1.122,0.002019,-8.876,-0.02916};
	double t[]={0,0.2,0.05,2.4,0.6,3.6};
	double d[]={0,1,4,9,1,8};
	double l[]={0,0,0,0,1,1};
	double g[]={0,0,0,0,1,1};
	delta = rho/(10.4477*28.9586);
	tau = 132.6312/T;
	etar=N[1]*pow(tau,t[1])*pow(delta,d[1])*exp(-g[1]*pow(delta,l[1]))
		+N[2]*pow(tau,t[2])*pow(delta,d[2])*exp(-g[2]*pow(delta,l[2]))
		+N[3]*pow(tau,t[3])*pow(delta,d[3])*exp(-g[3]*pow(delta,l[3]))
		+N[4]*pow(tau,t[4])*pow(delta,d[4])*exp(-g[4]*pow(delta,l[4]))
		+N[5]*pow(tau,t[5])*pow(delta,d[5])*exp(-g[5]*pow(delta,l[5]));

	return (eta0+etar)/1e6; // uPa-s to Pa-s
}

double AirClass::X_tilde(double T,double tau,double delta)
{
	// X_tilde is dimensionless
	// Equation 11 slightly rewritten
	double drho_dp,R_Air;
    R_Air = params.R_u/params.molemass;
	drho_dp=1.0/(R_Air*T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta)));
	return crit.p.Pa*delta/crit.rho*drho_dp;
}

double AirClass::conductivity_Trho(double T, double rho)
{
	/*
	E.W. Lemmon and R.T Jacobsen, Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon and Air
	International Journal of Thermophysics, Vol. 25, No. 1, January 2004, p.28 
	*/
	double e_k=103.3, //[K]
		   sigma=0.360, //[nm]
		   Tref=265.262, //[K]
		   zeta0=0.11, //[nm]
		   LAMBDA=0.055,
		   q_D=0.31; //[nm]
	double eta0,OMEGA,delta,tau,Tstar,lambda0,lambdar,num,
		cp,cv,OMEGA_tilde,OMEGA_tilde0,zeta,nu,gamma,R0,lambdac,k,
		pi=3.141592654,mu;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,1.308,1.405,-1.036,8.743,14.76,-16.62,3.793,-6.142,-0.3778};
	double t[]={0,0,-1.1,-0.3,0.1,0.0,0.5,2.7,0.3,1.3};
	double d[]={0,0,0,0,1,2,3,7,7,11};
	double l[]={0,0,0,0,0,0,2,2,2,2};
	double g[]={0,0,0,0,0,0,1,1,1,1};
	
	delta=rho/(10.4477*28.9586);
	tau=132.6312/T;
	Tstar=T/(e_k);

	OMEGA=exp(b[0]*powInt(log(Tstar),0)
			 +b[1]*powInt(log(Tstar),1)
		     +b[2]*powInt(log(Tstar),2)
			 +b[3]*powInt(log(Tstar),3)
		     +b[4]*powInt(log(Tstar),4));

	eta0=0.0266958*sqrt(28.9586*T)/(sigma*sigma*OMEGA);
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

	num=X_tilde(T,crit.T/T,delta)-X_tilde(Tref,crit.T/Tref,delta)*Tref/T;

	// no critical enhancement if numerator of Eq. 10 is negative
	if (num<0)
		return (lambda0+lambdar)/1e6;

	cp=Props('C','T',T,'D',rho,(char *)"Air");
	cv=Props('O','T',T,'D',rho,(char *)"Air");
	mu=Props('V','T',T,'D',rho,(char *)"Air")*1e6; //[uPa-s]

	zeta=zeta0*pow(num/LAMBDA,nu/gamma); //[nm]
	OMEGA_tilde=2.0/pi*((cp-cv)/cp*atan(zeta/q_D)+cv/cp*(zeta/q_D));
	OMEGA_tilde0=2.0/pi*(1.-exp(-1./(q_D/zeta+1.0/3.0*(zeta/q_D)*(zeta/q_D)/delta/delta)));
	lambdac=rho*(cp*1000.0)*k*R0*T/(6*pi*zeta*mu)*(OMEGA_tilde-OMEGA_tilde0)*1e18; // 1e18 is conversion to mW/m-K (not described in paper)

	return (lambda0+lambdar+lambdac)/1e3;
}

//
//double phi0_Air(double tau, double delta)
//{
//    double phi0=0;
//    int k;
//    
//    phi0=log(delta);
//    for (k=1;k<=5;k++)
//    {
//        phi0=phi0+N0[k]*powInt(tau,k-4);
//    }
//    phi0+=N0[6]*pow(tau,1.5)+N0[7]*log(tau);
//    for (k=8;k<=9;k++)
//    {
//        phi0+=N0[k]*log(1.0-exp(-N0[k+3]*tau));
//    }
//    phi0+=N0[10]*log(2.0/3.0+exp(N0[13]*tau));
//    return phi0;
//}
//
//double dphi0_dDelta_Air(double tau, double delta)
//{
//    return 1/delta;
//}
//
//double dphi02_dDelta2_Air(double tau, double delta)
//{
//    return -1.0/powInt(delta,2);
//}
//
//double dphi0_dTau_Air(double tau, double delta)
//{
//    double dphi0_dTau=0; int k;
//    for (k=1;k<=5;k++)
//    {
//        dphi0_dTau+=N0[k]*(k-4)*powInt(tau,k-5);
//    }
//    dphi0_dTau+=1.5*N0[6]*sqrt(tau)+N0[7]/tau;
//    dphi0_dTau+=N0[8]*N0[11]/(exp(N0[11]*tau)-1)+N0[9]*N0[12]/(exp(N0[12]*tau)-1);
//    dphi0_dTau+=N0[10]*N0[13]/(2.0/3.0*exp(-N0[13]*tau)+1);
//    
//    return dphi0_dTau;
//}
//
//double dphi02_dTau2_Air(double tau, double delta)
//{
//    double dphi02_dTau2=0;
//    int k;
//    
//    for (k=1;k<=5;k++)
//    {
//        dphi02_dTau2+=N0[k]*(k-4)*(k-5)*powInt(tau,k-6);
//    }
//    dphi02_dTau2+=0.75*N0[6]*pow(tau,-0.5)-N0[7]/(tau*tau);
//    dphi02_dTau2+=-N0[8]*N0[11]*N0[11]*exp(N0[11]*tau)/powInt(exp(N0[11]*tau)-1,2);
//    dphi02_dTau2+=-N0[9]*N0[12]*N0[12]*exp(N0[12]*tau)/powInt(exp(N0[12]*tau)-1,2);
//    dphi02_dTau2+=-(2.0/3.0)*N0[10]*N0[13]*N0[13]*exp(-N0[13]*tau)/powInt((2.0/3.0)*exp(-N0[13]*tau)+1,2);
//    return dphi02_dTau2;
//}
