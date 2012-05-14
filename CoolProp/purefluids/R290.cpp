/* Properties of Propane (R290)
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

static const double Tc=369.89, rhoc=220.4781, Pc=4251.2, M_R290=44.09562, _Ttriple=85.525;
             //           K          kg/m^3       kPa            kg/kmol          K
static const double n[]={0,
0.042910051,
1.7313671,
-2.4516524,
0.34157466,
-0.46047898,
-0.66847295,
0.20889705,
0.19421381,
-0.22917851,
-0.60405866,
0.066680654,
0.017534618,
0.33874242,
0.22228777,
-0.23219062,
-0.09220694,
-0.47575718,
-0.017486824};

static const int d[]={0,
4, //[ 1]
1, //[ 2]
1, //[ 3]
2, //[ 4]
2, //[ 5]
1, //[ 6]
3, //[ 7]
6, //[ 8]
6, //[ 9]
2, //[10]
3, //[11]
1, //[12]
1, //[13]
1, //[14]
2, //[15]
2, //[16]
4, //[17]
1  //[18]
};

static const double t[]={0.00, //offset for natural indices
1.00,
0.33,
0.80,
0.43,
0.90,
2.46,
2.09,
0.88,
1.09,
3.25,
4.62,
0.76,
2.50,
2.75,
3.05,
2.55,
8.40,
6.75};

static const int c[]={
0,0,0,0,0,0, // indices [0-5]
1,
1,
1,
1,
2,
2,
0,0,0,0,0,0,0 // indices [12-18]
};

// alpha instead of eta is used here for consistency with the definitions in R744.c upon which R290.c is based
static const double alpha[]={
0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-11]
0.963,
1.977,
1.917,
2.307,
2.546,
3.28,
14.6};

static const double beta[]={
0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-11]
2.33,
3.47,
3.15,
3.19,
0.92,
18.8,
547.8};

static const double GAMMA[]={
0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-11]
0.684,
0.829,
1.419,
0.817,
1.5,
1.426,
1.093};

static const double epsilon[]={
0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-11]
1.283,
0.6936,
0.788,
0.473,
0.8577,
0.271,
0.948};

//Constants for ideal gas expression
static const double a0[]={0.0,
    -4.970583,
    4.29352,
    3.043,
    5.874,
    9.337,
    7.922
};

static const double b0[]={0.0,
    0,0, //[1 and 2 are not used]
    1.062478,
    3.344237,
    5.363757,
    11.762957
};

//For the viscosity
static const double tv[]={
    0.0,			//[0]
    0.0,			//[1]
    0.0,			//[2]
    0.0,			//[3]
    0.0,			//[4]
    1.0,			//[5]
    1.0,			//[6]
    2.0,			//[7]
    2.0,			//[8]
    2.0,			//[9]
    3.0,			//[10]
    4.0,			//[11]
    1.0,			//[12]
    2.0				//[13]
};

static const double dv[]={
    0.0,			//[0]
    1.0,			//[1]
    2.0,			//[2]
    3.0,			//[3]
    13.0,			//[4]
    12.0,			//[5]
    16.0,			//[6]
    0.0,			//[7]
    18.0,			//[8]
    20.0,			//[9]
    13.0,			//[10]
    4.0,			//[11]
    0.0,			//[12]
    1.0				//[13]
};

static const double nv[]={
    0.0,			//[0]
    -0.7548580e-1,	//[1]
    0.7607150,		//[2]
    -0.1665680,		//[3]
    0.1627612e-5,	//[4]
    0.1443764e-4,	//[5]
    -0.2759086e-6,	//[6]
    -0.1032756,		//[7]
    -0.2498159e-7,	//[8]
    0.4069891e-8,	//[9]
    -0.1513434e-5,	//[10]
    0.2591327e-2,	//[11]
    0.5650076,		//[12]
    0.1207253		//[13]
};

int Load_R290(struct fluidParamsVals *Fluid)
{
    // Function pointers
    Fluid->funcs.phir=phir_R290;
    Fluid->funcs.dphir_dDelta=dphir_dDelta_R290;
    Fluid->funcs.dphir2_dDelta2=dphir2_dDelta2_R290;
    Fluid->funcs.dphir2_dDelta_dTau=dphir2_dDelta_dTau_R290;
    Fluid->funcs.dphir_dTau=dphir_dTau_R290;
    Fluid->funcs.dphir2_dTau2=dphir2_dTau2_R290;
    Fluid->funcs.phi0=phi0_R290;
    Fluid->funcs.dphi0_dDelta=dphi0_dDelta_R290;
    Fluid->funcs.dphi02_dDelta2=dphi02_dDelta2_R290;
    Fluid->funcs.dphi0_dTau=dphi0_dTau_R290;
    Fluid->funcs.dphi02_dTau2=dphi02_dTau2_R290;
    Fluid->funcs.rhosatL=rhosatL_R290;
    Fluid->funcs.rhosatV=rhosatV_R290;
    Fluid->funcs.psat=psat_R290;

    Fluid->funcs.visc=Viscosity_Trho_R290;
    Fluid->funcs.cond=Conductivity_Trho_R290;

    //Lookup table parameters
    Fluid->LUT.Tmin=200.0;
    Fluid->LUT.Tmax=800.0;
    Fluid->LUT.pmin=500;
    Fluid->LUT.pmax=16000;

    //Fluid parameters
    Fluid->Type=FLUIDTYPE_REFRIGERANT_PURE;
    Fluid->Tc=Tc;
    Fluid->rhoc=rhoc;
    Fluid->MM=M_R290;
    Fluid->pc=Pc;
    Fluid->Tt=_Ttriple;
    return 1;
}

double rhosatL_R290(double T)
{
    const double ti[]={0,0.345,0.74,2.6,7.2};
    const double Ni[]={0,1.82205,0.65802,0.21109,0.083973};
    double summer=1;
    int i;
    double theta;
    theta=1-T/Tc;
    for (i=1;i<=4;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return rhoc*summer;
}

double rhosatV_R290(double T)
{
    const double ti[]={0,0.3785,1.07,2.7,5.5,10,20};
    const double Ni[]={0,-2.4887,-5.1069,-12.174,-30.495,-52.192,-134.89};
    double summer=0,theta;
    int i;
    theta=1.0-T/Tc;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return rhoc*exp(summer);
}

double psat_R290(double T)
{
    const double ti[]={0,1.0,1.5,2.2,4.8,6.2};
    const double Ni[]={0,-6.7722,1.6938,-1.3341,-3.1876,0.94937};
    double summer=0,theta;
    int i;
    theta=1-T/Tc;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return Pc*exp(Tc/T*summer);
}

double Viscosity_Trho_R290(double T, double rho)
{
    /* 
    Properties taken from "A Reference Multiparameter Viscosity Equation 
    for Propane with an Optimized Functional Form" 
    by G. Scalabrin and P. Marchi and R. Span
    J. Phys. Chem. Ref. Data, Vol. 35, No. 3, 2006, 1415-1442
    */

    // inputs in T [K], and p [kPa]
    // output in Pa-s

    int i;
    double Tr,rhor,Hc=17.1045,etar, sum=0;

    //Reduced Temperature
    Tr=T/Tc;
    rhor=rho/rhoc;
    for(i=1;i<=11;i++)
    {
        sum+=nv[i]*pow(Tr,tv[i])*pow(rhor,dv[i]);
    }
    for(i=12;i<=13;i++)
    {
        sum+=exp(-rhor*rhor/2.0)*nv[i]*pow(Tr,tv[i])*pow(rhor,dv[i]);
    }
    etar=sum;

    return (exp(etar)-1)*Hc/1e6;
}

double Conductivity_Trho_R290(double T, double rho)
{
    /*Properties taken from "Measurement and Correlation of the Thermal Conductivity of 
    Propane from 86 K to 600 K at Pressures to 70 MPa" 
    by Kenneth N. Marsh, Richard A. Perkins, and Maria L. V. Ramires
    J. Chem. Eng. Data 2002, 47, 932-940

    The empirical critical enhancement is implemented
    */
    
    // output in kW/m-K

    double delta,lambda0,lambdar,lambdac,sum=0,tau;
    double DELTAT_c,DELTArho_c;
    int i;

    //Set constants required
    double B1[]={
        0.0,			//[0]
        -3.51153e-2,	//[1]
         1.70890e-1,	//[2]
        -1.47688e-1,	//[3]
         5.19283e-2,	//[4]
        -6.18662e-3		//[5]
    };
    double B2[]={
        0.0,			//[0]
         4.69195e-2,	//[1]
        -1.48616e-1,	//[2]
         1.32457e-1,	//[3]
        -4.85636e-2,	//[4]
         6.60414e-3		//[5]
    };
    double C[]={
        0.0,			//[0]
         3.66486e-4,	//[1]
        -2.21696e-3,	//[2]
         2.64213e+0		//[3]
    };
    double A[]={
         0.0,			//[0]
        -1.24778e-3,	//[1]
         8.16371e-3,	//[2]
         1.99374e-2,	//[3]
    };
    delta=rho/rhoc;
    tau=Tc/T;
    lambda0=A[1]+A[2]/tau+A[3]/(tau*tau);
    for(i=1;i<=5;i++)
    {
        sum+=(B1[i]+B2[i]/tau)*pow(delta,(double)i);
    }
    lambdar=sum;
    DELTAT_c=(1.0/tau-1.0);
    DELTArho_c=delta-1.0;
    lambdac=C[1]/(C[2]+fabs(DELTAT_c))*exp(-(C[3]*DELTArho_c)*(C[3]*DELTArho_c));

    return (lambda0+lambdar+lambdac)/1000.0;
}

/**************************************************/
/*          Private Property Functions            */
/**************************************************/

double phir_R290(double tau, double delta)
{ 
    
    int i;
    double phir=0,psi;
    
    for (i=1;i<=5;i++)
    {
        phir=phir+n[i]*powInt(delta,d[i])*pow(tau,t[i]);
    }
    
    for (i=6;i<=11;i++)
    {
        phir=phir+n[i]*powInt(delta,d[i])*pow(tau,t[i])*exp(-powInt(delta,c[i]));
    }
    
    for (i=12;i<=18;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        phir=phir+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi;
    }
    return phir;
}

double dphir_dDelta_R290(double tau, double delta)
{ 
    int i;
    double dphir_dDelta=0,psi;
    for (i=1;i<=5;i++)
    {
        dphir_dDelta+=n[i]*powInt(delta,d[i]-1)*pow(tau,t[i])*d[i];
    }
    for (i=6;i<=11;i++)
    {
        dphir_dDelta+=n[i]*powInt(delta,d[i]-1)*pow(tau,t[i])*exp(-powInt(delta,c[i]))*(d[i]-c[i]*powInt(delta,c[i]));
    }
    for (i=12;i<=18;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir_dDelta+=n[i]*powInt(delta,d[i]-1)*pow(tau,t[i])*psi*(d[i]-2.0*alpha[i]*delta*(delta-epsilon[i]));
    }
    return dphir_dDelta;
}

double dphir2_dDelta2_R290(double tau, double delta)
{ 
    
    int i;
    double di,ci;
    double dphir2_dDelta2=0,psi;
    for (i=1;i<=5;i++)
    {
        di=(double)d[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*di*(di-1.0)*powInt(delta,d[i]-2)*pow(tau,t[i]);
    }
    for (i=6;i<=11;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*exp(-powInt(delta,c[i]))*(powInt(delta,d[i]-2)*pow(tau,t[i])*( (di-ci*powInt(delta,c[i]))*(di-1.0-ci*powInt(delta,c[i])) - ci*ci*powInt(delta,c[i])));
    }
    for (i=12;i<=18;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dDelta2=dphir2_dDelta2+n[i]*pow(tau,t[i])*psi*(-2.0*alpha[i]*powInt(delta,d[i])+4.0*powInt(alpha[i],2)*powInt(delta,d[i])*powInt(delta-epsilon[i],2)-4.0*di*alpha[i]*powInt(delta,d[i]-1)*(delta-epsilon[i])+di*(di-1.0)*powInt(delta,d[i]-2));
    }
    return dphir2_dDelta2;
}

    
double dphir2_dDelta_dTau_R290(double tau, double delta)
{ 
    
    int i;
    double di, ci;
    double dphir2_dDelta_dTau=0,psi;

    for (i=1;i<=5;i++)
    {
        di=(double)d[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau + n[i]*di*t[i]*powInt(delta,d[i]-1)*pow(tau,t[i]-1.0);
    }
    for (i=6;i<=11;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau + n[i]*exp(-powInt(delta,c[i]))*powInt(delta,d[i]-1)*t[i]*pow(tau,t[i]-1.0)*(di-ci*powInt(delta,c[i]));
    }
    for (i=12;i<=18;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(di/delta-2.0*alpha[i]*(delta-epsilon[i]))*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    return dphir2_dDelta_dTau;
}

double dphir_dTau_R290(double tau, double delta)
{ 
    
    int i;
    double dphir_dTau=0,psi;
    
    for (i=1;i<=5;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powInt(delta,d[i])*pow(tau,t[i]-1.0);
    }
    for (i=6;i<=11;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powInt(delta,d[i])*pow(tau,t[i]-1.0)*exp(-powInt(delta,c[i]));
    }
    for (i=12;i<=18;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir_dTau=dphir_dTau+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    return dphir_dTau;
}

double dphir2_dTau2_R290(double tau, double delta)
{ 
    
    int i;
    double dphir2_dTau2=0,psi;
    
    for (i=1;i<=5;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powInt(delta,d[i])*pow(tau,t[i]-2.0);
    }
    for (i=6;i<=11;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powInt(delta,d[i])*pow(tau,t[i]-2.0)*exp(-powInt(delta,c[i]));
    }
    for (i=12;i<=18;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dTau2=dphir2_dTau2+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(powInt(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]),2)-t[i]/powInt(tau,2)-2.0*beta[i]);
    }
    return dphir2_dTau2;
}

double phi0_R290(double tau, double delta)
{
    double phi0=0;
    phi0=log(delta)+3*log(tau)+a0[1]+a0[2]*tau
        +a0[3]*log(1-exp(-b0[3]*tau))
        +a0[4]*log(1-exp(-b0[4]*tau))
        +a0[5]*log(1-exp(-b0[5]*tau))
        +a0[6]*log(1-exp(-b0[6]*tau));
    return phi0;
}

double dphi0_dDelta_R290(double tau, double delta)
{
    return 1/delta;
}

double dphi02_dDelta2_R290(double tau, double delta)
{
    return -1.0/powInt(delta,2);
}

double dphi0_dTau_R290(double tau, double delta)
{
    double dphi0_dTau=0;
    dphi0_dTau=3.0/tau+a0[2]
        +a0[3]*b0[3]*(1/(exp(b0[3]*tau)-1))
        +a0[4]*b0[4]*(1/(exp(b0[4]*tau)-1))
        +a0[5]*b0[5]*(1/(exp(b0[5]*tau)-1))
        +a0[6]*b0[6]*(1/(exp(b0[6]*tau)-1));
    return dphi0_dTau;
}

double dphi02_dTau2_R290(double tau, double delta)
{
    double dphi02_dTau2=0;
    dphi02_dTau2=-3.0/powInt(tau,2)
        -a0[3]*b0[3]*b0[3]*exp(b0[3]*tau)/powInt(exp(b0[3]*tau)-1.0,2)
        -a0[4]*b0[4]*b0[4]*exp(b0[4]*tau)/powInt(exp(b0[4]*tau)-1.0,2)
        -a0[5]*b0[5]*b0[5]*exp(b0[5]*tau)/powInt(exp(b0[5]*tau)-1.0,2)
        -a0[6]*b0[6]*b0[6]*exp(b0[6]*tau)/powInt(exp(b0[6]*tau)-1.0,2);
    return dphi02_dTau2;
}
