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

static const double Tc=126.192, rhoc=313.3, Pc=3395.8, M_Nitrogen=28.01348, _Ttriple=63.151;
             //          K           kg/m^3       kPa              g/mol              K
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
static const int d[]={0,
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
static const int c[]={0,
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

int Load_Nitrogen(struct fluidParamsVals *Fluid)
{
    // Function pointers
    Fluid->funcs.phir=phir_Nitrogen;
    Fluid->funcs.dphir_dDelta=dphir_dDelta_Nitrogen;
    Fluid->funcs.dphir2_dDelta2=dphir2_dDelta2_Nitrogen;
    Fluid->funcs.dphir2_dDelta_dTau=dphir2_dDelta_dTau_Nitrogen;
    Fluid->funcs.dphir_dTau=dphir_dTau_Nitrogen;
    Fluid->funcs.dphir2_dTau2=dphir2_dTau2_Nitrogen;
    Fluid->funcs.phi0=phi0_Nitrogen;
    Fluid->funcs.dphi0_dDelta=dphi0_dDelta_Nitrogen;
    Fluid->funcs.dphi02_dDelta2=dphi02_dDelta2_Nitrogen;
    Fluid->funcs.dphi0_dTau=dphi0_dTau_Nitrogen;
    Fluid->funcs.dphi02_dTau2=dphi02_dTau2_Nitrogen;
    Fluid->funcs.rhosatL=rhosatL_Nitrogen;
    Fluid->funcs.rhosatV=rhosatV_Nitrogen;
    Fluid->funcs.psat=psat_Nitrogen;

    Fluid->funcs.visc=Viscosity_Trho_Nitrogen;
    Fluid->funcs.cond=Conductivity_Trho_Nitrogen;

    //Lookup table parameters
    Fluid->LUT.Tmin=140.0;
    Fluid->LUT.Tmax=800.0;
    Fluid->LUT.pmin=50;
    Fluid->LUT.pmax=16000;

    //Fluid parameters
    Fluid->Type=FLUIDTYPE_REFRIGERANT_PURE;
    Fluid->Tc=Tc;
    Fluid->rhoc=rhoc;
    Fluid->MM=M_Nitrogen;
    Fluid->pc=Pc;
    Fluid->Tt=_Ttriple;
    return 1;
}

double rhosatL_Nitrogen(double T)
{
    const double ti[]={0,0.3294,2.0/3.0,8.0/3.0,35.0/6.0};
    const double Ni[]={0,1.48654237,-0.280476066,0.0894143085,-0.119879866};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+Ni[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(summer);
}

double rhosatV_Nitrogen(double T)
{
    const double ti[]={0,0.34,5.0/6.0,7.0/6.0,13.0/6.0,14.0/3.0};
    const double Ni[]={0,-1.70127164,-3.70402649,1.29859383,-0.561424977,-2.68505381};
    double summer=0;
    int i;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(Tc/T*summer);
}

double psat_Nitrogen(double T)
{
    const double ti[]={0,1.0,1.5,2.5,5.0};
    const double Ni[]={0,-6.12445284,1.26327220,-0.765910082,-1.77570564};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+Ni[i]*pow(1-T/Tc,ti[i]);
    }
    return Pc*exp(Tc/T*summer);
}

double Viscosity_Trho_Nitrogen(double T, double rho)
{
	double e_k=98.94, //[K]
		   sigma=0.3656; //[nm]
	double eta0,etar,OMEGA,delta,tau,Tstar;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,10.72,0.03989,0.001208,-7.402,4.620};
	double t[]={0,0.1,0.25,3.2,0.9,0.3};
	double d[]={0,2,10,12,2,1};
	double l[]={0,0,1,1,2,3};
	double g[]={0,0,1,1,1,1};

	delta=rho/rhoc;
	tau=Tc/T;
	Tstar=T/(e_k);
	OMEGA=exp(b[0]*powInt(log(Tstar),0)
			 +b[1]*powInt(log(Tstar),1)
		     +b[2]*powInt(log(Tstar),2)
			 +b[3]*powInt(log(Tstar),3)
		     +b[4]*powInt(log(Tstar),4));

	eta0=0.0266958*sqrt(M_Nitrogen*T)/(sigma*sigma*OMEGA);
	etar=N[1]*pow(tau,t[1])*pow(delta,d[1])*exp(-g[1]*pow(delta,l[1]))
		+N[2]*pow(tau,t[2])*pow(delta,d[2])*exp(-g[2]*pow(delta,l[2]))
		+N[3]*pow(tau,t[3])*pow(delta,d[3])*exp(-g[3]*pow(delta,l[3]))
		+N[4]*pow(tau,t[4])*pow(delta,d[4])*exp(-g[4]*pow(delta,l[4]))
		+N[5]*pow(tau,t[5])*pow(delta,d[5])*exp(-g[5]*pow(delta,l[5]));

	return (eta0+etar)/1e6; // uPa-s to Pa-s
}

double X_tilde(double T,double tau,double delta)
{
	// X_tilde is dimensionless
	// Equation 11 slightly rewritten
	double drho_dp,R_Nitrogen;
	R_Nitrogen=8.31447215/M_Nitrogen;
	drho_dp=1.0/(R_Nitrogen*T*(1+2*delta*dphir_dDelta_Nitrogen(tau,delta)+delta*delta*dphir2_dDelta2_Nitrogen(tau,delta)));
	return Pc*delta/rhoc*drho_dp;
}

double Conductivity_Trho_Nitrogen(double T, double rho)
{
	double e_k=98.94, //[K]
		   sigma=0.3656, //[nm]
		   Tref=252.384, //[K]
		   zeta0=0.17, //[nm]
		   LAMBDA=0.055,
		   q_D=0.40; //[nm]
	double eta0,OMEGA,delta,tau,Tstar,lambda0,lambdar,num,
		cp,cv,OMEGA_tilde,OMEGA_tilde0,zeta,nu,gamma,R0,lambdac,k,
		pi=3.141592654,mu;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,1.511,2.117,-3.332,8.862,31.11,-73.13,20.03,-0.7096,0.2672};
	double t[]={0,0,-1.0,-0.7,0.0,0.03,0.2,0.8,0.6,1.9};
	double d[]={0,0,0,0,1,2,3,4,8,10};
	double l[]={0,0,0,0,0,0,1,2,2,2};
	double g[]={0,0,0,0,0,0,1,1,1,1};
	
	delta=rho/rhoc;
	tau=Tc/T;
	Tstar=T/(e_k);

	OMEGA=exp(b[0]*powInt(log(Tstar),0)
			 +b[1]*powInt(log(Tstar),1)
		     +b[2]*powInt(log(Tstar),2)
			 +b[3]*powInt(log(Tstar),3)
		     +b[4]*powInt(log(Tstar),4));

	eta0=0.0266958*sqrt(M_Nitrogen*T)/(sigma*sigma*OMEGA);
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

	num=X_tilde(T,Tc/T,delta)-X_tilde(Tref,Tc/Tref,delta)*Tref/T;

	// no critical enhancement if numerator of Eq. 10 is negative
	if (num<0)
		return (lambda0+lambdar)/1e6;

	cp=Props('C','T',T,'D',rho,"Nitrogen");
	cv=Props('O','T',T,'D',rho,"Nitrogen");
	mu=Props('V','T',T,'D',rho,"Nitrogen")*1e6; //[uPa-s]

	zeta=zeta0*pow(num/LAMBDA,nu/gamma); //[nm]
	OMEGA_tilde=2.0/pi*((cp-cv)/cp*atan(zeta/q_D)+cv/cp*(zeta/q_D));
	OMEGA_tilde0=2.0/pi*(1.-exp(-1./(q_D/zeta+1.0/3.0*(zeta/q_D)*(zeta/q_D)/delta/delta)));
	lambdac=rho*(cp*1000.0)*k*R0*T/(6*pi*zeta*mu)*(OMEGA_tilde-OMEGA_tilde0)*1e18; // 1e18 is conversion to mW/m-K (not described in paper)

	return (lambda0+lambdar+lambdac)/1e6;
}

/**************************************************/
/*          Private Property Functions            */
/**************************************************/

double phir_Nitrogen(double tau, double delta)
{ 
    
    int i;
    double phir=0,psi;
    
    for (i=1;i<=6;i++)
    {
        phir=phir+n[i]*powInt(delta,d[i])*pow(tau,t[i]);
    }
    
    for (i=7;i<=32;i++)
    {
        phir=phir+n[i]*powInt(delta,d[i])*pow(tau,t[i])*exp(-powInt(delta,c[i]));
    }
    
    for (i=33;i<=36;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        phir=phir+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi;
    }
    return phir;
}

double dphir_dDelta_Nitrogen(double tau, double delta)
{ 
    int i;
    double dphir_dDelta=0,psi;
    double di, ci;
    for (i=1;i<=6;i++)
    {
        di=(double)d[i];
        dphir_dDelta=dphir_dDelta+n[i]*di*powInt(delta,d[i]-1)*pow(tau,t[i]);
    }
    for (i=7;i<=32;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir_dDelta=dphir_dDelta+n[i]*exp(-powInt(delta,c[i]))*(powInt(delta,d[i]-1)*pow(tau,t[i])*(di-ci*powInt(delta,c[i])));
    }
    for (i=33;i<=36;i++)
    {
        di=(double)d[i];        
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir_dDelta=dphir_dDelta+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(di/delta-2.0*alpha[i]*(delta-epsilon[i]));
    }
    return dphir_dDelta;
}

double dphir2_dDelta2_Nitrogen(double tau, double delta)
{ 
    
    int i;
    double di,ci;
    double dphir2_dDelta2=0,psi;
    for (i=1;i<=6;i++)
    {
        di=(double)d[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*di*(di-1.0)*powInt(delta,d[i]-2)*pow(tau,t[i]);
    }
    for (i=7;i<=32;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*exp(-powInt(delta,c[i]))*(powInt(delta,d[i]-2)*pow(tau,t[i])*( (di-ci*powInt(delta,c[i]))*(di-1.0-ci*powInt(delta,c[i])) - ci*ci*powInt(delta,c[i])));
    }
    for (i=33;i<=36;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dDelta2=dphir2_dDelta2+n[i]*pow(tau,t[i])*psi*(-2.0*alpha[i]*powInt(delta,d[i])+4.0*powInt(alpha[i],2)*powInt(delta,d[i])*powInt(delta-epsilon[i],2)-4.0*di*alpha[i]*powInt(delta,d[i]-1)*(delta-epsilon[i])+di*(di-1.0)*powInt(delta,d[i]-2));
    }
    return dphir2_dDelta2;
}

double dphir2_dDelta_dTau_Nitrogen(double tau, double delta)
{ 
    
    int i;
    double di, ci;
    double dphir2_dDelta_dTau=0,psi;

    for (i=1;i<=6;i++)
    {
        di=(double)d[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau + n[i]*di*t[i]*powInt(delta,d[i]-1)*pow(tau,t[i]-1.0);
    }
    for (i=7;i<=32;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau + n[i]*exp(-powInt(delta,c[i]))*powInt(delta,d[i]-1)*t[i]*pow(tau,t[i]-1.0)*(di-ci*powInt(delta,c[i]));
    }
    for (i=33;i<=36;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(di/delta-2.0*alpha[i]*(delta-epsilon[i]))*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    return dphir2_dDelta_dTau;
}

double dphir_dTau_Nitrogen(double tau, double delta)
{ 
    
    int i;
    double dphir_dTau=0,psi;
    
    for (i=1;i<=6;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powInt(delta,d[i])*pow(tau,t[i]-1.0);
    }
    for (i=7;i<=32;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powInt(delta,d[i])*pow(tau,t[i]-1.0)*exp(-powInt(delta,c[i]));
    }
    for (i=33;i<=36;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir_dTau=dphir_dTau+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    return dphir_dTau;
}


double dphir2_dTau2_Nitrogen(double tau, double delta)
{ 
    
    int i;
    double dphir2_dTau2=0,psi;
    
    for (i=1;i<=6;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powInt(delta,d[i])*pow(tau,t[i]-2.0);
    }
    for (i=7;i<=32;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powInt(delta,d[i])*pow(tau,t[i]-2.0)*exp(-powInt(delta,c[i]));
    }
    for (i=33;i<=36;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dTau2=dphir2_dTau2+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(powInt(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]),2)-t[i]/powInt(tau,2)-2.0*beta[i]);
    }
    return dphir2_dTau2;
}

/* Maxima code for derivatives of ideal gas residual Helmholtz:
A:log(delta)+a0[1]*log(tau)+a0[2]+a0[3]*tau+a0[4]/tau+a0[5]/tau/tau+a0[6]/tau/tau/tau+a0[7]*log(1-exp(-a0[8]*tau));
diff(A,delta);
diff(%,delta);
diff(A,tau);
diff(%,tau);
*/

double phi0_Nitrogen(double tau, double delta)
{
    double phi0=0;
    
    phi0=log(delta)+a0[1]*log(tau)+a0[2]+a0[3]*tau+a0[4]/tau+a0[5]/tau/tau+a0[6]/tau/tau/tau+a0[7]*log(1-exp(-a0[8]*tau));
    return phi0;
}

double dphi0_dDelta_Nitrogen(double tau, double delta)
{
    return 1/delta;
}

double dphi02_dDelta2_Nitrogen(double tau, double delta)
{
    return -1.0/powInt(delta,2);
}

double dphi0_dTau_Nitrogen(double tau, double delta)
{
    return a0[1]/tau+a0[3]-a0[4]/tau/tau-2.0*a0[5]/tau/tau/tau-3*a0[6]/tau/tau/tau/tau+(a0[7]*a0[8]*exp(-a0[8]*tau))/(1-exp(-a0[8]*tau));
}

double dphi02_dTau2_Nitrogen(double tau, double delta)
{
    return -(a0[7]*a0[8]*a0[8]*exp(-a0[8]*tau))/(1-exp(-a0[8]*tau))-(a0[7]*a0[8]*a0[8]*exp(-2*a0[8]*tau))/(1-exp(-a0[8]*tau))/(1-exp(-a0[8]*tau))-a0[1]/tau/tau+(2*a0[4])/powInt(tau,3)+(6*a0[5])/powInt(tau,4)+(12*a0[6])/powInt(tau,5);
}
