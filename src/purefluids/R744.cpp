/* Properties of Carbon Dioxide (R744)
by Ian Bell

Themo properties from 
"A New Equation of State for Carbon Dioxide Covering the Fluid Region from the 
Triple Point Temperature to 1100 K at Pressures up to 800 MPa", 
R. Span and W. Wagner, J. Phys. Chem. Ref. Data, v. 25, 1996

WARNING: Thermal conductivity not coded!!

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

static double alpha[40],beta[43],GAMMA[40],epsilon[40],a[43],b[43],A[43],B[43],C[43],D[43],a0[9],theta0[9];
static const double Tc=304.128, M_R744=44.01, rhoc=467.606, Pc=7377.3, _Ttriple=216.59;
             //    K               g/mol       kg/m^3          kPa              K
static const double n[]={0,    
 0.3885682320316100E+00,
 0.2938547594274000E+01,
-0.5586718853493400E+01,
-0.7675319959247700E+00,
 0.3172900558041600E+00,
 0.5480331589776700E+00,
 0.1227941122033500E+00,
 
 0.2165896154322000E+01,
 0.1584173510972400E+01,
-0.2313270540550300E+00,
 0.5811691643143600E-01,
-0.5536913720538200E-00,
 0.4894661590942200E-00,
-0.2427573984350100E-01,
 0.6249479050167800E-01,
-0.1217586022524600E+00,
-0.3705568527008600E+00,
-0.1677587970042600E-01,
-0.1196073663798700E+00,
-0.4561936250877800E-01,
 0.3561278927034600E-01, 
-0.7442772713205200E-02,
-0.1739570490243200E-02,
-0.2181012128952700E-01,
 0.2433216655923600E-01,
-0.3744013342346300E-01,
 0.1433871575687800E-00,
-0.1349196908328600E-00,
-0.2315122505348000E-01,
 0.1236312549290100E-01,
 0.2105832197294000E-02,
-0.3395851902636800E-03,
 0.5599365177159200E-02,
-0.3033511805564600E-03,

-0.2136548868832000E+03,
 0.2664156914927200E+05,
-0.2402721220455700E+05,
-0.2834160342399900E+03,
 0.2124728440017900E+03,
 
-0.6664227654075100E+00,
 0.7260863234989700E+00,
 0.5506866861284200E-01};

static const int d[]={0,
1,
1,
1,
1,
2,
2,
3,
1,
2,
4,
5,
5,
5,
6,
6,
6,
1,
1,
4,
4,
4,
7,
8,
2,
3,
3,
5,
5,
6,
7,
8,
10,
4,
8,
2,
2,
2,
3,
3};

static const double t[]={0.00,
0.00,
0.75,
1.00,
2.00,
0.75,
2.00,
0.75,
1.50,
1.50,
2.50,
0.00,
1.50,
2.00,
0.00,
1.00,
2.00,
3.00,
6.00,
3.00,
6.00,
8.00,
6.00,
0.00,
7.00,
12.00,
16.00,
22.00,
24.00,
16.00,
24.00,
8.00,
2.00,
28.00,
14.00,
1.00,
0.00,
1.00,
3.00,
3.00};

static const int c[]={0,0,0,0,0,0,0,0,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
2,
2,
2,
2,
2,
2,
3,
3,
3,
4,
4,
4,
4,
4,
4,
5,
6};

void setCoeffs(void)
{
alpha[35]=25.0;
alpha[36]=25.0;
alpha[37]=25.0;
alpha[38]=15.0;
alpha[39]=20.0;

beta[35]=325.0;
beta[36]=300.0;
beta[37]=300.0;
beta[38]=275.0;
beta[39]=275.0;

GAMMA[35]=1.16;
GAMMA[36]=1.19;
GAMMA[37]=1.19;
GAMMA[38]=1.25;
GAMMA[39]=1.22;

epsilon[35]=1.00;
epsilon[36]=1.00;
epsilon[37]=1.00;
epsilon[38]=1.00;
epsilon[39]=1.00;

a[40]=3.5;
a[41]=3.5;
a[42]=3.0;

b[40]=0.875;
b[41]=0.925;
b[42]=0.875;

beta[40]=0.300;
beta[41]=0.300;
beta[42]=0.300;

A[40]=0.700;
A[41]=0.700;
A[42]=0.700;

B[40]=0.3;
B[41]=0.3;
B[42]=1.0;

C[40]=10.0;
C[41]=10.0;
C[42]=12.5;

D[40]=275.0;
D[41]=275.0;
D[42]=275.0;

//Constants for ideal gas expression
a0[1]=8.37304456;
a0[2]=-3.70454304;
a0[3]=2.500000;
a0[4]=1.99427042;
a0[5]=0.62105248;
a0[6]=0.41195293;
a0[7]=1.04028922;
a0[8]=0.08327678;

theta0[4]=3.15163;
theta0[5]=6.11190;
theta0[6]=6.77708;
theta0[7]=11.32384;
theta0[8]=27.08792;
}

int Load_R744(struct fluidParamsVals *Fluid)
{
    // Function pointers
    Fluid->funcs.phir=phir_R744;
    Fluid->funcs.dphir_dDelta=dphir_dDelta_R744;
    Fluid->funcs.dphir2_dDelta2=dphir2_dDelta2_R744;
    Fluid->funcs.dphir2_dDelta_dTau=dphir2_dDelta_dTau_R744;
    Fluid->funcs.dphir_dTau=dphir_dTau_R744;
    Fluid->funcs.dphir2_dTau2=dphir2_dTau2_R744;
    Fluid->funcs.phi0=phi0_R744;
    Fluid->funcs.dphi0_dDelta=dphi0_dDelta_R744;
    Fluid->funcs.dphi02_dDelta2=dphi02_dDelta2_R744;
    Fluid->funcs.dphi0_dTau=dphi0_dTau_R744;
    Fluid->funcs.dphi02_dTau2=dphi02_dTau2_R744;
    Fluid->funcs.rhosatL=rhosatL_R744;
    Fluid->funcs.rhosatV=rhosatV_R744;
    Fluid->funcs.psat=psat_R744;

    Fluid->funcs.visc=Viscosity_Trho_R744;
    Fluid->funcs.cond=Conductivity_Trho_R744;

    //Lookup table parameters
    Fluid->LUT.Tmin=220.0;
    Fluid->LUT.Tmax=800.0;
    Fluid->LUT.pmin=500.0;
    Fluid->LUT.pmax=16000;

    //Fluid parameters
    Fluid->Type=FLUIDTYPE_REFRIGERANT_PURE;
    Fluid->Tc=Tc;
    Fluid->rhoc=rhoc;
    Fluid->MM=M_R744;
    Fluid->pc=Pc;
    Fluid->Tt=_Ttriple;
    return 1;
}

double rhosatV_R744(double T)
{
    const double ti[]={0,0.340,1.0/2.0,1.0,7.0/3.0,14.0/3.0};
    const double ai[]={0,-1.7074879,-0.82274670,-4.6008549,-10.111178,-29.742252};
    double summer=0;
    int i;
    for (i=1;i<=5;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(summer);
}

double rhosatL_R744(double T)
{
    const double ti[]={0,0.340,1.0/2.0,10.0/6.0,11.0/6.0};
    const double ai[]={0,1.9245108,-0.62385555,-0.32731127,0.39245142};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(summer);
}

double Viscosity_Trho_R744(double T,double rho)
{
	int i;
	double e_k=251.196,Tstar,sumGstar=0.0,Gstar,eta0,delta_eta;
	double a[]={0.235156,-0.491266,5.211155e-2,5.347906e-2,-1.537102e-2};
	double d11=0.4071119e-2,d21=0.7198037e-4,d64=0.2411697e-16,d81=0.2971072e-22,d82=-0.1627888e-22;

	Tstar=T/e_k;
	for (i=0;i<=4;i++)
	{
		sumGstar=sumGstar+a[i]*powInt(log(Tstar),i);
	}
	Gstar=exp(sumGstar);
	eta0=1.00697*sqrt(T)/Gstar;
	delta_eta=d11*rho+d21*rho*rho*d64*powInt(rho,6)/powInt(Tstar,3)+d81*powInt(rho,8)+d82*powInt(rho,8)/Tstar;

	return (eta0+delta_eta)/1e6;
}

double Conductivity_Trho_R744(double T,double rho)
{
	fprintf(stderr,"Thermal conductivity not coded for R744 (CO2).  Sorry.\n");
	return _HUGE;
}

/**************************************************/
/*          Private Property Functions            */
/**************************************************/

double phir_R744(double tau, double delta)
{ 
    
    int i;
    double phir=0,theta,DELTA,PSI,psi;
    
    for (i=1;i<=7;i++)
    {
        phir=phir+n[i]*powInt(delta,d[i])*pow(tau,t[i]);
    }
    
    for (i=8;i<=34;i++)
    {
        phir=phir+n[i]*powInt(delta,d[i])*pow(tau,t[i])*exp(-powInt(delta,c[i]));
    }
    
    for (i=35;i<=39;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        phir=phir+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi;
    }
    
    for (i=40;i<=42;i++)
    {
        theta=(1.0-tau)+A[i]*pow(powInt(delta-1.0,2),1/(2*beta[i]));
        DELTA=powInt(theta,2)+B[i]*pow(powInt(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powInt(delta-1.0,2)-D[i]*powInt(tau-1.0,2));
        phir=phir+n[i]*pow(DELTA,b[i])*delta*PSI;
    }
    
    return phir;
}

double dphir_dDelta_R744(double tau, double delta)
{ 
    int i;
    double dphir_dDelta=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta,psi;
    double di, ci;
    for (i=1;i<=7;i++)
    {
        di=(double)d[i];

        dphir_dDelta=dphir_dDelta+n[i]*di*powInt(delta,d[i]-1)*pow(tau,t[i]);
    }
    for (i=8;i<=34;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir_dDelta=dphir_dDelta+n[i]*exp(-powInt(delta,c[i]))*(powInt(delta,d[i]-1)*pow(tau,t[i])*(di-ci*powInt(delta,c[i])));
    }
    for (i=35;i<=39;i++)
    {
        di=(double)d[i];        
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir_dDelta=dphir_dDelta+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(di/delta-2.0*alpha[i]*(delta-epsilon[i]));
    }
    for (i=40;i<=42;i++)
    {
        theta=(1.0-tau)+A[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=powInt(theta,2)+B[i]*pow(powInt(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powInt(delta-1.0,2)-D[i]*powInt(tau-1.0,2));
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(powInt(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        dphir_dDelta=dphir_dDelta+n[i]*(pow(DELTA,b[i])*(PSI+delta*dPSI_dDelta)+dDELTAbi_dDelta*delta*PSI);
    }
    return dphir_dDelta;
}

double dphir2_dDelta2_R744(double tau, double delta)
{ 
    
    int i;
    double di,ci;
    
    double dphir2_dDelta2=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta,psi,dPSI2_dDelta2,dDELTAbi2_dDelta2,dDELTA2_dDelta2;
    for (i=1;i<=7;i++)
    {
        di=(double)d[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*di*(di-1.0)*powInt(delta,d[i]-2)*pow(tau,t[i]);
    }
    for (i=8;i<=34;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*exp(-powInt(delta,c[i]))*(powInt(delta,d[i]-2)*pow(tau,t[i])*( (di-ci*powInt(delta,c[i]))*(di-1.0-ci*powInt(delta,c[i])) - ci*ci*powInt(delta,c[i])));
    }
    for (i=35;i<=39;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dDelta2=dphir2_dDelta2+n[i]*pow(tau,t[i])*psi*(-2.0*alpha[i]*powInt(delta,d[i])+4.0*powInt(alpha[i],2)*powInt(delta,d[i])*powInt(delta-epsilon[i],2)-4.0*di*alpha[i]*powInt(delta,d[i]-1)*(delta-epsilon[i])+di*(di-1.0)*powInt(delta,d[i]-2));
    }
    for (i=40;i<=42;i++)
    {
               
        theta=(1.0-tau)+A[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=powInt(theta,2)+B[i]*pow(powInt(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powInt(delta-1.0,2)-D[i]*powInt(tau-1.0,2));
        
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(powInt(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
        dPSI2_dDelta2=(2.0*C[i]*powInt(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
        dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+powInt(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(powInt(delta-1.0,2),a[i]-2.0)+2.0*powInt(A[i]/beta[i],2)*powInt(pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
        dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*powInt(dDELTA_dDelta,2));
        
        dphir2_dDelta2=dphir2_dDelta2+n[i]*(pow(DELTA,b[i])*(2.0*dPSI_dDelta+delta*dPSI2_dDelta2)+2.0*dDELTAbi_dDelta*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta2*delta*PSI);
    }
    return dphir2_dDelta2;
}

    
double dphir2_dDelta_dTau_R744(double tau, double delta)
{ 
    
    int i;
    double di, ci;
    double dphir2_dDelta_dTau=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta,psi;
    double dPSI2_dDelta_dTau, dDELTAbi2_dDelta_dTau, dPSI_dTau, dDELTAbi_dTau;
    for (i=1;i<=7;i++)
    {
        di=(double)d[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*di*t[i]*powInt(delta,d[i]-1)*pow(tau,t[i]-1.0);
    }
    for (i=8;i<=34;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*exp(-powInt(delta,c[i]))*powInt(delta,d[i]-1)*t[i]*pow(tau,t[i]-1.0)*(di-ci*powInt(delta,c[i]));
    }
    for (i=35;i<=39;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(di/delta-2.0*alpha[i]*(delta-epsilon[i]))*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    for (i=40;i<=42;i++)
    {
        
        theta=(1.0-tau)+A[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=powInt(theta,2)+B[i]*pow(powInt(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powInt(delta-1.0,2)-D[i]*powInt(tau-1.0,2));
        
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(powInt(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        
        dPSI2_dDelta_dTau=4.0*C[i]*D[i]*(delta-1.0)*(tau-1.0)*PSI;
        dDELTAbi2_dDelta_dTau=-A[i]*b[i]*2.0/beta[i]*pow(DELTA,b[i]-1.0)*(delta-1.0)*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0)-2.0*theta*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0)*dDELTA_dDelta;
        
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*(pow(DELTA,b[i])*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+delta*dDELTAbi_dDelta*dPSI_dTau+ dDELTAbi_dTau*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta_dTau*delta*PSI);
    }
    return dphir2_dDelta_dTau;
}

double dphir_dTau_R744(double tau, double delta)
{ 
    
    int i;
    double dphir_dTau=0,theta,DELTA,PSI,dPSI_dTau,dDELTAbi_dTau,psi;
    
    for (i=1;i<=7;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powInt(delta,d[i])*pow(tau,t[i]-1.0);
    }
    
    for (i=8;i<=34;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powInt(delta,d[i])*pow(tau,t[i]-1.0)*exp(-powInt(delta,c[i]));
    }
    
    for (i=35;i<=39;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir_dTau=dphir_dTau+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    
    for (i=40;i<=42;i++)
    {
        theta=(1.0-tau)+A[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=powInt(theta,2)+B[i]*pow(powInt(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powInt(delta-1.0,2)-D[i]*powInt(tau-1.0,2));
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        dphir_dTau=dphir_dTau+n[i]*delta*(dDELTAbi_dTau*PSI+pow(DELTA,b[i])*dPSI_dTau);
    }
    
    return dphir_dTau;
}


double dphir2_dTau2_R744(double tau, double delta)
{ 
    
    int i;
    double dphir2_dTau2=0,theta,DELTA,PSI,dPSI_dTau,dDELTAbi_dTau,psi,dPSI2_dTau2,dDELTAbi2_dTau2;
    
    for (i=1;i<=7;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powInt(delta,d[i])*pow(tau,t[i]-2.0);
    }
    
    for (i=8;i<=34;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powInt(delta,d[i])*pow(tau,t[i]-2.0)*exp(-powInt(delta,c[i]));
    }
    
    for (i=35;i<=39;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dTau2=dphir2_dTau2+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(powInt(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]),2)-t[i]/powInt(tau,2)-2.0*beta[i]);
    }
    
    for (i=40;i<=42;i++)
    {
        theta=(1.0-tau)+A[i]*pow(powInt(delta-1.0,2),1/(2*beta[i]));
        DELTA=powInt(theta,2)+B[i]*pow(powInt(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powInt(delta-1.0,2)-D[i]*powInt(tau-1.0,2));
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        dPSI2_dTau2=(2.0*D[i]*powInt(tau-1.0,2)-1.0)*2.0*D[i]*PSI;
        dDELTAbi2_dTau2=2.0*b[i]*pow(DELTA,b[i]-1.0)+4.0*powInt(theta,2)*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0);
        dphir2_dTau2=dphir2_dTau2+n[i]*delta*(dDELTAbi2_dTau2*PSI+2.0*dDELTAbi_dTau*dPSI_dTau+pow(DELTA,b[i])*dPSI2_dTau2);
    }
    
    return dphir2_dTau2;
}

double dphi0_dDelta_R744(double tau, double delta)
{
    return 1/delta;
}

double dphi02_dDelta2_R744(double tau, double delta)
{
    return -1.0/powInt(delta,2);
}

double phi0_R744(double tau, double delta)
{
    double phi0=0;
    int i;
    
    phi0=log(delta)+a0[1]+a0[2]*tau+a0[3]*log(tau);
    for (i=4;i<=8;i++)
    {
        phi0=phi0+a0[i]*log(1.0-exp(-theta0[i]*tau));
    }
    return phi0;
}


double dphi0_dTau_R744(double tau, double delta)
{
    double dphi0_dTau=0;
    int i;
    dphi0_dTau=a0[2]+a0[3]/tau;

    for (i=4;i<=8;i++)
    {
        dphi0_dTau=dphi0_dTau+a0[i]*theta0[i]*(1.0/(1.0-exp(-theta0[i]*tau))-1.0);
    }
    return dphi0_dTau;
}

double dphi02_dTau2_R744(double tau, double delta)
{
    double dphi02_dTau2=0;
    int i;
    
    dphi02_dTau2=-a0[3]/powInt(tau,2);
    for (i=4;i<=8;i++)
    {
        dphi02_dTau2=dphi02_dTau2-a0[i]*powInt(theta0[i],2)*exp(-theta0[i]*tau)/powInt(1.0-exp(-theta0[i]*tau),2);
    }
    return dphi02_dTau2;
}

double psat_R744(double T)
{
    const double ti[]={0,1.0,1.5,2.0,4.0};
    const double ai[]={0,-7.0602087,1.9391218,-1.6463597,-3.2995634};
    double summer=0;
    int i;
    setCoeffs();
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1-T/Tc,ti[i]);
    }
    return Pc*exp(Tc/T*summer);
}
