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
#include "PropErrorCodes.h"
#include "PropMacros.h"
#include "CoolProp.h"

static const double Tj=132.6312, R_Air=0.287117125828, rhoj=302.5507652, pj=3785.02, M_Air=28.9586, _Ttriple=59.75;
             //           K             kJ/kg-K       kg/m^3     kPa            kg/kmol          K
// Critical values (if needed)
// static const double Tc=132.5306,pc=3786.0;

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

static const int i[]={0,
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

static const double j[]={0.00,
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

static const int l[]={
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

int Load_Air(struct fluidParamsVals *Fluid)
{
    // Function pointers
    Fluid->funcs.phir=phir_Air;
    Fluid->funcs.dphir_dDelta=dphir_dDelta_Air;
    Fluid->funcs.dphir2_dDelta2=dphir2_dDelta2_Air;
    Fluid->funcs.dphir2_dDelta_dTau=dphir2_dDelta_dTau_Air;
    Fluid->funcs.dphir_dTau=dphir_dTau_Air;
    Fluid->funcs.dphir2_dTau2=dphir2_dTau2_Air;
    Fluid->funcs.phi0=phi0_Air;
    Fluid->funcs.dphi0_dDelta=dphi0_dDelta_Air;
    Fluid->funcs.dphi02_dDelta2=dphi02_dDelta2_Air;
    Fluid->funcs.dphi0_dTau=dphi0_dTau_Air;
    Fluid->funcs.dphi02_dTau2=dphi02_dTau2_Air;
    Fluid->funcs.rhosatL=rhosatL_Air;
    Fluid->funcs.rhosatV=rhosatV_Air;
    Fluid->funcs.p_dp=pdp_Air;
    Fluid->funcs.p_bp=pbp_Air;
    Fluid->funcs.visc=Viscosity_Trho_Air;
    Fluid->funcs.cond=Conductivity_Trho_Air;

    //Lookup table parameters
    Fluid->LUT.Tmin=220.0;
    Fluid->LUT.Tmax=800.0;
    Fluid->LUT.pmin=120.0;
    Fluid->LUT.pmax=16000.0;

    //Fluid parameters
    Fluid->Type=FLUIDTYPE_REFRIGERANT_PSEUDOPURE;
    Fluid->Tc=Tj;
    Fluid->rhoc=rhoj;
    Fluid->MM=M_Air;
    Fluid->pc=pj;
    Fluid->Tt=_Ttriple;
    return 1;
}

double rhosatL_Air(double T)
{
	const double ti[]={0,0.65,0.85,0.95,1.1};
    const double Ni[]={0,43.3413,-240.073,285.139,-88.3366,-0.892181};
    double summer=0; int k;
    for (k=1;k<=4;k++)
    {
        summer=summer+Ni[k]*pow(1.0-T/Tj,ti[k]);
    }
    return rhoj*(1+summer+Ni[5]*log(T/Tj));
}

double rhosatV_Air(double T)
{
	const double ti[]={0,0.41,1.0,2.8,6.5};
    const double Ni[]={0,-2.0466,-4.7520,-13.259,-47.652};
    double summer=0; int k;
    for (k=1;k<=4;k++)
    {
        summer=summer+Ni[k]*pow(1.0-T/Tj,ti[k]);
    }
    return rhoj*exp(summer);
}

double pbp_Air(double T)
{
	const double Ni[]={0,0.2260724,-7.080499,5.700283,-12.44017,17.81926,-10.81364};
	double summer=0; int k;
    for (k=1;k<=6;k++)
    {
        summer=summer+Ni[k]*pow(1-T/Tj,(double)k/2.0);
    }
	return pj*exp(Tj/T*summer);
}

double pdp_Air(double T)
{
	const double Ni[]={0,-0.1567266,-5.539635,0,0,0.7567212,0,0,-3.514322};
	double summer=0; int k;
    for (k=1;k<=8;k++)
    {
        summer=summer+Ni[k]*pow(1-T/Tj,(double)k/2.0);
    }
	return pj*exp(Tj/T*summer);
}

double Viscosity_Trho_Air(double T, double rho)
{
	/*
	E.W. Lemmon and R.T Jacobsen, Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon and Air
	International Journal of Thermophysics, Vol. 25, No. 1, January 2004, p.28 
	*/

	double e_k=103.3, //[K]
		   sigma=0.360; //[nm]
	double eta0,etar,OMEGA,delta,tau,Tstar;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,10.72,1.122,0.002019,-8.876,-0.02916};
	double t[]={0,0.2,0.05,2.4,0.6,3.6};
	double d[]={0,1,4,9,1,8};
	double l[]={0,0,0,0,1,1};
	double g[]={0,0,0,0,1,1};

	delta=rho/rhoj;
	tau=Tj/T;
	Tstar=T/(e_k);
	OMEGA=exp(b[0]*powInt(log(Tstar),0)
			 +b[1]*powInt(log(Tstar),1)
		     +b[2]*powInt(log(Tstar),2)
			 +b[3]*powInt(log(Tstar),3)
		     +b[4]*powInt(log(Tstar),4));

	eta0=0.0266958*sqrt(M_Air*T)/(sigma*sigma*OMEGA);
	etar=N[1]*pow(tau,t[1])*pow(delta,d[1])*exp(-g[1]*pow(delta,l[1]))
		+N[2]*pow(tau,t[2])*pow(delta,d[2])*exp(-g[2]*pow(delta,l[2]))
		+N[3]*pow(tau,t[3])*pow(delta,d[3])*exp(-g[3]*pow(delta,l[3]))
		+N[4]*pow(tau,t[4])*pow(delta,d[4])*exp(-g[4]*pow(delta,l[4]))
		+N[5]*pow(tau,t[5])*pow(delta,d[5])*exp(-g[5]*pow(delta,l[5]));

	return (eta0+etar)/1e6; // uPa-s to Pa-s
}

static double X_tilde(double T,double tau,double delta)
{
	// X_tilde is dimensionless
	// Equation 11 slightly rewritten
	double drho_dp;
	drho_dp=1.0/(R_Air*T*(1+2*delta*dphir_dDelta_Air(tau,delta)+delta*delta*dphir2_dDelta2_Air(tau,delta)));
	return pj*delta/rhoj*drho_dp;
}

double Conductivity_Trho_Air(double T, double rho)
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
	
	delta=rho/rhoj;
	tau=Tj/T;
	Tstar=T/(e_k);

	OMEGA=exp(b[0]*powInt(log(Tstar),0)
			 +b[1]*powInt(log(Tstar),1)
		     +b[2]*powInt(log(Tstar),2)
			 +b[3]*powInt(log(Tstar),3)
		     +b[4]*powInt(log(Tstar),4));

	eta0=0.0266958*sqrt(M_Air*T)/(sigma*sigma*OMEGA);
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

	num=X_tilde(T,Tj/T,delta)-X_tilde(Tref,Tj/Tref,delta)*Tref/T;

	// no critical enhancement if numerator of Eq. 10 is negative
	if (num<0)
		return (lambda0+lambdar)/1e6;

	cp=Props('C','T',T,'D',rho,"Air");
	cv=Props('O','T',T,'D',rho,"Air");
	mu=Props('V','T',T,'D',rho,"Air")*1e6; //[uPa-s]

	zeta=zeta0*pow(num/LAMBDA,nu/gamma); //[nm]
	OMEGA_tilde=2.0/pi*((cp-cv)/cp*atan(zeta/q_D)+cv/cp*(zeta/q_D));
	OMEGA_tilde0=2.0/pi*(1.-exp(-1./(q_D/zeta+1.0/3.0*(zeta/q_D)*(zeta/q_D)/delta/delta)));
	lambdac=rho*(cp*1000.0)*k*R0*T/(6*pi*zeta*mu)*(OMEGA_tilde-OMEGA_tilde0)*1e18; // 1e18 is conversion to mW/m-K (not described in paper)

	return (lambda0+lambdar+lambdac)/1e6;
}

double B_Air(double tau)
{
	// given by B*rhoc=lim(delta --> 0) [dphir_ddelta(tau)]
	return 1.0/rhoj*dphir_dDelta_Air(tau,1e-12);
}

double dBdT_Air(double tau)
{
	// given by B*rhoc^2=lim(delta --> 0) [dphir2_ddelta2(tau)]
	return -1.0/rhoj*tau*tau/Tj*dphir2_dDelta_dTau_Air(tau,1e-12);
}

double C_Air(double tau)
{
	// given by B*rhoc^2=lim(delta --> 0) [dphir2_ddelta2(tau)]
	return 1.0/(rhoj*rhoj)*dphir2_dDelta2_Air(tau,1e-12);
}

double dCdT_Air(double tau)
{
	// given by B*rhoc^2=lim(delta --> 0) [dphir2_ddelta2(tau)]
	return -1.0/(rhoj*rhoj)*tau*tau/Tj*dphir3_dDelta2_dTau_Air(tau,1e-12);
}

/**************************************************/
/*          Private Property Functions            */
/**************************************************/

double phir_Air(double tau, double delta)
{ 
    int k;
    double phir=0;
    
    for (k=1;k<=10;k++)
    {
        phir=phir+N[k]*powInt(delta,i[k])*pow(tau,j[k]);
    }
    
    for (k=11;k<=19;k++)
    {
        phir=phir+N[k]*powInt(delta,i[k])*pow(tau,j[k])*exp(-powInt(delta,l[k]));
    }
    return phir;
}

double dphir_dDelta_Air(double tau, double delta)
{ 
    int k;
    double dphir_dDelta=0;
    for (k=1;k<=10;k++)
    {
        dphir_dDelta+=i[k]*N[k]*powInt(delta,i[k]-1)*pow(tau,j[k]);
    }
    for (k=11;k<=19;k++)
    {
        dphir_dDelta+=N[k]*powInt(delta,i[k]-1)*pow(tau,j[k])*exp(-powInt(delta,l[k]))*(i[k]-l[k]*powInt(delta,l[k]));
    }
    return dphir_dDelta;
}

double dphir2_dDelta2_Air(double tau, double delta)
{ 
    int k;
    double di,ci;
    double dphir2_dDelta2=0;
    for (k=1;k<=10;k++)
    {
        di=(double)i[k];
        dphir2_dDelta2=dphir2_dDelta2+N[k]*di*(di-1.0)*powInt(delta,i[k]-2)*pow(tau,j[k]);
    }
    for (k=11;k<=19;k++)
    {
        di=(double)i[k];
        ci=(double)l[k];
        dphir2_dDelta2=dphir2_dDelta2+N[k]*exp(-powInt(delta,l[k]))*(powInt(delta,i[k]-2)*pow(tau,j[k])*( (di-ci*powInt(delta,l[k]))*(di-1.0-ci*powInt(delta,l[k])) - ci*ci*powInt(delta,l[k])));
    }
    return dphir2_dDelta2;
}

    
double dphir2_dDelta_dTau_Air(double tau, double delta)
{ 
    int k;
    double dphir2_dDelta_dTau=0;

    for (k=1;k<=10;k++)
    {
        dphir2_dDelta_dTau+=i[k]*j[k]*N[k]*powInt(delta,i[k]-1)*pow(tau,j[k]-1.0);
    }
    for (k=11;k<=19;k++)
    {
        dphir2_dDelta_dTau+=j[k]*N[k]*powInt(delta,i[k]-1)*pow(tau,j[k]-1.0)*exp(-powInt(delta,l[k]))*(i[k]-l[k]*powInt(delta,l[k]));
    }
    return dphir2_dDelta_dTau;
}

double dphir3_dDelta2_dTau_Air(double tau, double delta)
{ 
    int k;
    double dphir3_dDelta2_dTau=0;

    for (k=1;k<=10;k++)
    {
        dphir3_dDelta2_dTau+=i[k]*(i[k]-1)*j[k]*N[k]*powInt(delta,i[k]-2)*pow(tau,j[k]-1.0);
    }
    for (k=11;k<=19;k++)
    {
        dphir3_dDelta2_dTau+=j[k]*N[k]*powInt(delta,i[k]-2)*pow(tau,j[k]-1.0)*exp(-powInt(delta,l[k]))*((i[k]-l[k]*powInt(delta,l[k]))*(i[k]-1-l[k]*powInt(delta,l[k]))-l[k]*l[k]*powInt(delta,l[k]));
    }
    return dphir3_dDelta2_dTau;
}

double dphir_dTau_Air(double tau, double delta)
{ 
    int k;
    double dphir_dTau=0;
    
    for (k=1;k<=10;k++)
    {
        dphir_dTau=dphir_dTau+N[k]*j[k]*powInt(delta,i[k])*pow(tau,j[k]-1.0);
    }
    for (k=11;k<=19;k++)
    {
        dphir_dTau=dphir_dTau+N[k]*j[k]*powInt(delta,i[k])*pow(tau,j[k]-1.0)*exp(-powInt(delta,l[k]));
    }
    return dphir_dTau;
}

double dphir2_dTau2_Air(double tau, double delta)
{ 
    
    int k;
    double dphir2_dTau2=0;
    
    for (k=1;k<=10;k++)
    {
        dphir2_dTau2=dphir2_dTau2+N[k]*j[k]*(j[k]-1.0)*powInt(delta,i[k])*pow(tau,j[k]-2.0);
    }
    for (k=11;k<=19;k++)
    {
        dphir2_dTau2=dphir2_dTau2+N[k]*j[k]*(j[k]-1.0)*powInt(delta,i[k])*pow(tau,j[k]-2.0)*exp(-powInt(delta,l[k]));
    }
    return dphir2_dTau2;
}

double phi0_Air(double tau, double delta)
{
    double phi0=0;
    int k;
    
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
    return phi0;
}

double dphi0_dDelta_Air(double tau, double delta)
{
    return 1/delta;
}

double dphi02_dDelta2_Air(double tau, double delta)
{
    return -1.0/powInt(delta,2);
}

double dphi0_dTau_Air(double tau, double delta)
{
    double dphi0_dTau=0; int k;
    for (k=1;k<=5;k++)
    {
        dphi0_dTau+=N0[k]*(k-4)*powInt(tau,k-5);
    }
    dphi0_dTau+=1.5*N0[6]*sqrt(tau)+N0[7]/tau;
    dphi0_dTau+=N0[8]*N0[11]/(exp(N0[11]*tau)-1)+N0[9]*N0[12]/(exp(N0[12]*tau)-1);
    dphi0_dTau+=N0[10]*N0[13]/(2.0/3.0*exp(-N0[13]*tau)+1);
    
    return dphi0_dTau;
}

double dphi02_dTau2_Air(double tau, double delta)
{
    double dphi02_dTau2=0;
    int k;
    
    for (k=1;k<=5;k++)
    {
        dphi02_dTau2+=N0[k]*(k-4)*(k-5)*powInt(tau,k-6);
    }
    dphi02_dTau2+=0.75*N0[6]*pow(tau,-0.5)-N0[7]/(tau*tau);
    dphi02_dTau2+=-N0[8]*N0[11]*N0[11]*exp(N0[11]*tau)/powInt(exp(N0[11]*tau)-1,2);
    dphi02_dTau2+=-N0[9]*N0[12]*N0[12]*exp(N0[12]*tau)/powInt(exp(N0[12]*tau)-1,2);
    dphi02_dTau2+=-(2.0/3.0)*N0[10]*N0[13]*N0[13]*exp(-N0[13]*tau)/powInt((2.0/3.0)*exp(-N0[13]*tau)+1,2);
    return dphi02_dTau2;
}
