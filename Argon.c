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

static const double Tc=150.687, rhoc=535.6, Pc=4863.0, M_Argon=39.948, _Ttriple=83.806;
             //           K          kg/m^3     kPa            kg/kmol          K
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

static const int d[]={0,
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

static const int c[]={
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
static const double a0[]={0.0,
	8.31666243,
	-4.94651164
};

int Load_Argon(struct fluidParamsVals *Fluid)
{
    // Function pointers
    Fluid->funcs.phir=phir_Argon;
    Fluid->funcs.dphir_dDelta=dphir_dDelta_Argon;
    Fluid->funcs.dphir2_dDelta2=dphir2_dDelta2_Argon;
    Fluid->funcs.dphir2_dDelta_dTau=dphir2_dDelta_dTau_Argon;
    Fluid->funcs.dphir_dTau=dphir_dTau_Argon;
    Fluid->funcs.dphir2_dTau2=dphir2_dTau2_Argon;
    Fluid->funcs.phi0=phi0_Argon;
    Fluid->funcs.dphi0_dDelta=dphi0_dDelta_Argon;
    Fluid->funcs.dphi02_dDelta2=dphi02_dDelta2_Argon;
    Fluid->funcs.dphi0_dTau=dphi0_dTau_Argon;
    Fluid->funcs.dphi02_dTau2=dphi02_dTau2_Argon;
    Fluid->funcs.rhosatL=rhosatL_Argon;
    Fluid->funcs.rhosatV=rhosatV_Argon;
    Fluid->funcs.psat=psat_Argon;

    Fluid->funcs.visc=Viscosity_Trho_Argon;
    Fluid->funcs.cond=Conductivity_Trho_Argon;

    //Lookup table parameters
    Fluid->LUT.Tmin=220.0;
    Fluid->LUT.Tmax=800.0;
    Fluid->LUT.pmin=50;
    Fluid->LUT.pmax=16000;

    //Fluid parameters
    Fluid->Type=FLUIDTYPE_REFRIGERANT_PURE;
    Fluid->Tc=Tc;
    Fluid->rhoc=rhoc;
    Fluid->MM=M_Argon;
    Fluid->pc=Pc;
    Fluid->Tt=_Ttriple;
    return 1;
}

double psat_Argon(double T)
{
    const double ti[]={0,1.0,1.5,2.0,4.5};
    const double ai[]={0,-5.9409785,1.3553888,-0.46497607,-1.5399043};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1-T/Tc,ti[i]);
    }
    return Pc*exp(Tc/T*summer);
}

double rhosatL_Argon(double T)
{
    const double ti[]={0,0.334,2.0/3.0,7.0/3.0,4.0};
    const double ai[]={0,1.5004262,-0.31381290,0.086461622,-0.041477525};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(summer);
}

double rhosatV_Argon(double T)
{
    const double ti[]={0,0.345,5.0/6.0,1.0,13.0/3.0};
    const double ai[]={0,-1.70695656,-4.02739448,1.55177558,-2.30683228};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(Tc/T*summer);
}


double Viscosity_Trho_Argon(double T, double rho)
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

	delta=rho/rhoc;
	tau=Tc/T;
	Tstar=T/(e_k);
	OMEGA=exp(b[0]*powInt(log(Tstar),0)
			 +b[1]*powInt(log(Tstar),1)
		     +b[2]*powInt(log(Tstar),2)
			 +b[3]*powInt(log(Tstar),3)
		     +b[4]*powInt(log(Tstar),4));

	eta0=0.0266958*sqrt(M_Argon*T)/(sigma*sigma*OMEGA);
	etar=N[1]*pow(tau,t[1])*pow(delta,d[1])*exp(-g[1]*pow(delta,l[1]))
		+N[2]*pow(tau,t[2])*pow(delta,d[2])*exp(-g[2]*pow(delta,l[2]))
		+N[3]*pow(tau,t[3])*pow(delta,d[3])*exp(-g[3]*pow(delta,l[3]))
		+N[4]*pow(tau,t[4])*pow(delta,d[4])*exp(-g[4]*pow(delta,l[4]))
		+N[5]*pow(tau,t[5])*pow(delta,d[5])*exp(-g[5]*pow(delta,l[5]))
		+N[6]*pow(tau,t[6])*pow(delta,d[6])*exp(-g[6]*pow(delta,l[6]));

	return (eta0+etar)/1e6; // uPa-s to Pa-s
}

static double X_tilde(double T,double tau,double delta)
{
	// X_tilde is dimensionless
	// Equation 11 slightly rewritten
	double drho_dp,R_Argon;
	R_Argon=8.31447215/M_Argon;
	drho_dp=1.0/(R_Argon*T*(1+2*delta*dphir_dDelta_Argon(tau,delta)+delta*delta*dphir2_dDelta2_Argon(tau,delta)));
	return Pc*delta/rhoc*drho_dp;
}

double Conductivity_Trho_Argon(double T, double rho)
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
	
	delta=rho/rhoc;
	tau=Tc/T;
	Tstar=T/(e_k);

	OMEGA=exp(b[0]*powInt(log(Tstar),0)
			 +b[1]*powInt(log(Tstar),1)
		     +b[2]*powInt(log(Tstar),2)
			 +b[3]*powInt(log(Tstar),3)
		     +b[4]*powInt(log(Tstar),4));

	eta0=0.0266958*sqrt(M_Argon*T)/(sigma*sigma*OMEGA);
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

	num=X_tilde(T,Tc/T,delta)-X_tilde(Tref,Tc/Tref,delta)*Tref/T;

	// no critical enhancement if numerator of Eq. 10 is negative
	if (num<0)
		return (lambda0+lambdar)/1e6;

	cp=Props('C','T',T,'D',rho,"Argon");
	cv=Props('O','T',T,'D',rho,"Argon");
	mu=Props('V','T',T,'D',rho,"Argon")*1e6; //[uPa-s]

	zeta=zeta0*pow(num/LAMBDA,nu/gamma); //[nm]
	OMEGA_tilde=2.0/pi*((cp-cv)/cp*atan(zeta/q_D)+cv/cp*(zeta/q_D));
	OMEGA_tilde0=2.0/pi*(1.-exp(-1./(q_D/zeta+1.0/3.0*(zeta/q_D)*(zeta/q_D)/delta/delta)));
	lambdac=rho*(cp*1000.0)*k*R0*T/(6*pi*zeta*mu)*(OMEGA_tilde-OMEGA_tilde0)*1e18; // 1e18 is conversion to mW/m-K (not described in paper)

	return (lambda0+lambdar+lambdac)/1e6;
}


double phir_Argon(double tau, double delta)
{ 
    
    int i;
    double phir=0,psi;
    
    for (i=1;i<=12;i++)
    {
        phir=phir+n[i]*powInt(delta,d[i])*pow(tau,t[i]);
    }
    
    for (i=13;i<=37;i++)
    {
        phir=phir+n[i]*powInt(delta,d[i])*pow(tau,t[i])*exp(-powInt(delta,c[i]));
    }
    
    for (i=38;i<=41;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        phir=phir+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi;
    }
    return phir;
}

double dphir_dDelta_Argon(double tau, double delta)
{ 
    int i;
    double dphir_dDelta=0,psi;
    double di, ci;
    for (i=1;i<=12;i++)
    {
        di=(double)d[i];
        dphir_dDelta=dphir_dDelta+n[i]*di*powInt(delta,d[i]-1)*pow(tau,t[i]);
    }
    for (i=13;i<=37;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir_dDelta=dphir_dDelta+n[i]*exp(-powInt(delta,c[i]))*(powInt(delta,d[i]-1)*pow(tau,t[i])*(di-ci*powInt(delta,c[i])));
    }
    for (i=38;i<=41;i++)
    {
        di=(double)d[i];        
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir_dDelta=dphir_dDelta+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(di/delta-2.0*alpha[i]*(delta-epsilon[i]));
    }
    return dphir_dDelta;
}

double dphir2_dDelta2_Argon(double tau, double delta)
{ 
    
    int i;
    double di,ci;
    double dphir2_dDelta2=0,psi;
    for (i=1;i<=12;i++)
    {
        di=(double)d[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*di*(di-1.0)*powInt(delta,d[i]-2)*pow(tau,t[i]);
    }
    for (i=13;i<=37;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*exp(-powInt(delta,c[i]))*(powInt(delta,d[i]-2)*pow(tau,t[i])*( (di-ci*powInt(delta,c[i]))*(di-1.0-ci*powInt(delta,c[i])) - ci*ci*powInt(delta,c[i])));
    }
    for (i=38;i<=41;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dDelta2=dphir2_dDelta2+n[i]*pow(tau,t[i])*psi*(-2.0*alpha[i]*powInt(delta,d[i])+4.0*powInt(alpha[i],2)*powInt(delta,d[i])*powInt(delta-epsilon[i],2)-4.0*di*alpha[i]*powInt(delta,d[i]-1)*(delta-epsilon[i])+di*(di-1.0)*powInt(delta,d[i]-2));
    }
    return dphir2_dDelta2;
}

    
double dphir2_dDelta_dTau_Argon(double tau, double delta)
{ 
    
    int i;
    double di, ci;
    double dphir2_dDelta_dTau=0,psi;

    for (i=1;i<=12;i++)
    {
        di=(double)d[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau + n[i]*di*t[i]*powInt(delta,d[i]-1)*pow(tau,t[i]-1.0);
    }
    for (i=13;i<=37;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau + n[i]*exp(-powInt(delta,c[i]))*powInt(delta,d[i]-1)*t[i]*pow(tau,t[i]-1.0)*(di-ci*powInt(delta,c[i]));
    }
    for (i=38;i<=41;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(di/delta-2.0*alpha[i]*(delta-epsilon[i]))*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    return dphir2_dDelta_dTau;
}

double dphir_dTau_Argon(double tau, double delta)
{ 
    
    int i;
    double dphir_dTau=0,psi;
    
    for (i=1;i<=12;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powInt(delta,d[i])*pow(tau,t[i]-1.0);
    }
    for (i=13;i<=37;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powInt(delta,d[i])*pow(tau,t[i]-1.0)*exp(-powInt(delta,c[i]));
    }
    for (i=38;i<=41;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir_dTau=dphir_dTau+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    return dphir_dTau;
}


double dphir2_dTau2_Argon(double tau, double delta)
{ 
    
    int i;
    double dphir2_dTau2=0,psi;
    
    for (i=1;i<=12;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powInt(delta,d[i])*pow(tau,t[i]-2.0);
    }
    for (i=13;i<=37;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powInt(delta,d[i])*pow(tau,t[i]-2.0)*exp(-powInt(delta,c[i]));
    }
    for (i=38;i<=41;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dTau2=dphir2_dTau2+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(powInt(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]),2)-t[i]/powInt(tau,2)-2.0*beta[i]);
    }
    return dphir2_dTau2;
}

double phi0_Argon(double tau, double delta)
{
    double phi0=0;
    
    phi0=log(delta)+a0[1]+a0[2]*tau+1.5*log(tau);
    return phi0;
}

double dphi0_dDelta_Argon(double tau, double delta)
{
    return 1/delta;
}

double dphi02_dDelta2_Argon(double tau, double delta)
{
    return -1.0/powInt(delta,2);
}

double dphi0_dTau_Argon(double tau, double delta)
{
    double dphi0_dTau=0;
    dphi0_dTau=a0[2]+1.5/tau;
    return dphi0_dTau;
}

double dphi02_dTau2_Argon(double tau, double delta)
{
    double dphi02_dTau2=0;
    dphi02_dTau2=-1.5/powInt(tau,2);
    return dphi02_dTau2;
}
