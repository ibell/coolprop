/*
Properties for R410A.  
by Ian Bell

Pseudo-pure fluid thermo props from 
"Pseudo-pure fluid Equations of State for the Refrigerant Blends R410A, R404A, R507C and R407C" 
by E.W. Lemmon, Int. J. Thermophys. v. 24, n4, 2003

In order to call the exposed functions, rho_, h_, s_, cp_,...... 
there are three different ways the inputs can be passed, and this is expressed by the Types integer flag.  
These macros are defined in the PropMacros.h header file:
1) First parameter temperature, second parameter pressure ex: h_R410A(260,354.7,1)=274
    In this case, the lookup tables are built if needed and then interpolated
2) First parameter temperature, second parameter density ex: h_R410A(260,13.03,2)=274
    Density and temp plugged directly into EOS
3) First parameter temperature, second parameter pressure ex: h_R410A(260,354.7,3)=274
    Density solved for, then plugged into EOS (can be quite slow)
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

static const double a[]={
    36.8871,		//[0]
    7.15807,		//[1]
    -46.87575,		//[2]
    2.0623,			//[3]
    5.9751,			//[4]
    1.5612			//[5]
};

static const double b[]={
    0,				//[0]
    0,				//[1]
    -0.1,			//[2]
    2.02326,		//[3]
    5.00154,		//[4]
    11.2484			//[5]
};

static const double N[]={
     0.0,			//[0]
     0.987252,		//[1]
    -1.03017,		//[2]
     1.17666,		//[3]
    -0.138991,		//[4]
     0.00302373,	//[5]
    -2.53639,		//[6]
    -1.96680,		//[7]
    -0.830480,		//[8]
     0.172477,		//[9]
    -0.261116,		//[10]
    -0.0745473,		//[11]
     0.679757,		//[12]
    -0.652431,		//[13]
     0.0553849,		//[14]
    -0.0710970,		//[15]
    -0.000875332,	//[16]
     0.0200760,		//[17]
    -0.0139761,		//[18]
    -0.0185110,		//[19]
     0.0171939,		//[20]
    -0.00482049		//[21]
};

static const double j[]={
    0.0,			//[0]
    0.44,			//[1]
    1.2,			//[2]
    2.97,			//[3]
    2.95,			//[4]
    0.2,			//[5]
    1.93,			//[6]
    1.78,			//[7]
    3.0,			//[8]
    0.2,			//[9]
    0.74,			//[10]
    3.0,			//[11]
    2.1,			//[12]
    4.3,			//[13]
    0.25,			//[14]
    7.0,			//[15]
    4.7,			//[16]
    13.0,			//[17]
    16.0,			//[18]
    25.0,			//[19]
    17.0,			//[20]
    7.4				//[21]
};

static const int i[]={
    0,				//[0]
    1,				//[1]
    1,				//[2]
    1,				//[3]
    2,				//[4]
    5,				//[5]
    1,				//[6]
    2,				//[7]
    3,				//[8]
    5,				//[9]
    5,				//[10]
    5,				//[11]
    1,				//[12]
    1,				//[13]
    4,				//[14]
    4,				//[15]
    9,				//[16]
    2,				//[17]
    2,				//[18]
    4,				//[19]
    5,				//[20]
    6				//[21]
};

static const int L[]={
    0,				//[0]
    0,				//[1]
    0,				//[2]
    0,				//[3]
    0,				//[4]
    0,				//[5]
    1,				//[6]
    1,				//[7]
    1,				//[8]
    1,				//[9]
    1,				//[10]
    1,				//[11]
    2,				//[12]
    2,				//[13]
    2,				//[14]
    2,				//[15]
    2,				//[16]
    3,				//[17]
    3,				//[18]
    3,				//[19]
    3,				//[20]
    3				//[21]
};

static const int g[]={
    0,				//[0]
    0,				//[1]
    0,				//[2]
    0,				//[3]
    0,				//[4]
    0,				//[5]
    1,				//[6]
    1,				//[7]
    1,				//[8]
    1,				//[9]
    1,				//[10]
    1,				//[11]
    1,				//[12]
    1,				//[13]
    1,				//[14]
    1,				//[15]
    1,				//[16]
    1,				//[17]
    1,				//[18]
    1,				//[19]
    1,				//[20]
    1,				//[21]
};

static const double Nbp[]={
    0.0,			//[0]
    -7.2818,		//[1]
    2.5093,			//[2]
    -3.2695,		//[3]
    -2.8022			//[4]
};
    
static const double tbp[]={
    0.0,			//[0]
    1.0,			//[1]
    1.8,			//[2]
    2.4,			//[3]
    4.9				//[4]
};

static const double Ndp[]={
    0.0,			//[0]
    -7.4411,		//[1]
    1.9883,			//[2]
    -2.4925,		//[3]
    -3.2633			//[4]
};

static const double tdp[]={
    0.0,			//[0]
    1.0,			//[1]
    1.6,			//[2]
    2.4,			//[3]
    5.0			    //[4]
};

static const double M=72.5824; //[g/mol]
static const double Tm=344.494; //[K]
static const double pm=4901.2; //[MPa--> kPa]
static const double pc=4810; //[MPa--> kPa] From (Calm 2007 HPAC Engineering)
static const double rhom=459.0300696; //6.324*M; //[mol/dm^3--> kg/m^3]
static const double _Ttriple=200.0; //[K]

int Load_R410A(struct fluidParamsVals *Fluid)
{
    // Function pointers
    Fluid->funcs.phir=ar_R410A;
    Fluid->funcs.dphir_dDelta=dar_ddelta_R410A;
    Fluid->funcs.dphir2_dDelta2=d2ar_ddelta2_R410A;
    Fluid->funcs.dphir2_dDelta_dTau=d2ar_ddelta_dtau_R410A;
    Fluid->funcs.dphir_dTau=dar_dtau_R410A;
    Fluid->funcs.dphir2_dTau2=d2ar_dtau2_R410A;
    Fluid->funcs.phi0=a0_R410A;
    Fluid->funcs.dphi0_dDelta=da0_ddelta_R410A;
    Fluid->funcs.dphi02_dDelta2=d2a0_ddelta2_R410A;
    Fluid->funcs.dphi0_dTau=da0_dtau_R410A;
    Fluid->funcs.dphi02_dTau2=d2a0_dtau2_R410A;
    Fluid->funcs.rhosatL=rhosatL_R410A;
    Fluid->funcs.rhosatV=rhosatV_R410A;
    Fluid->funcs.p_dp=p_dp_R410A;
    Fluid->funcs.p_bp=p_bp_R410A;
    Fluid->funcs.visc=Viscosity_Trho_R410A;
    Fluid->funcs.cond=Conductivity_Trho_R410A;
    
    //Lookup table parameters
    Fluid->LUT.Tmin=150.0;
    Fluid->LUT.Tmax=550.0;
    Fluid->LUT.pmin=70.03;
    Fluid->LUT.pmax=6000.0;
    
    //Fluid parameters
    Fluid->Type=FLUIDTYPE_REFRIGERANT_PSEUDOPURE;
    Fluid->Tc=Tm;
    Fluid->rhoc=rhom;
    Fluid->MM=M;
    Fluid->pc=pc;
    Fluid->Tt=_Ttriple;
    return 1;
}

double p_bp_R410A(double T)
{
    //Bubble point of R410A
    double sum=0,theta;
    int k;
    theta=1-T/Tm;
    
    for (k=1;k<=4;k++)
    {
        sum+=Nbp[k]*pow(theta,tbp[k]);
    }
    return pm*exp(Tm/T*sum);
}
    
double p_dp_R410A(double T)
{
    //Dew point of R410A
    double sum=0,theta;
    int k;
    theta=1-T/Tm;

    for (k=1;k<=4;k++)
    {
        sum+=Ndp[k]*pow(theta,tdp[k]);
    }
    return pm*exp(Tm/T*sum);
}

double rhosatV_R410A(double T)
{
    double THETA,a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8;
    THETA=1-T/344.5;

    a1 =       -4.02;
    a2 =       -18.8;
    a3 =      -18.16; 
    a4 =      -17.37;
    a5 =      -16.15;
    a6 =       -7.87;
    a7 =      -20.46;
    a8 =      -18.59;
    b1 =      0.4803; 
    b2 =       9.217;  
    b3 =        7.52;  
    b4 =       3.471; 
    b5 =       5.124; 
    b6 =       1.498; 
    b7 =       9.831; 
    b8 =        8.99; 

    return exp(a1*pow(THETA,b1)+a2*pow(THETA,b2)+a3*pow(THETA,b3)+a4*pow(THETA,b4)+a5*pow(THETA,b5)+a6*pow(THETA,b6)+a7*pow(THETA,b7)+a8*pow(THETA,b8))*459.53;
}

double rhosatL_R410A(double T)
{
    double THETA,a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6;
    THETA=1-T/344.5;

    a1 =       2.876;
    a2 =      -2.465;
    a3 =     -0.9235;
    a4 =     -0.4798;
    a5 =      0.7394;
    a6 =       1.762; 
    b1 =       0.365;
    b2 =      0.6299;
    b3 =       1.864; 
    b4 =       1.932; 
    b5 =       2.559; 
    b6 =       1.149; 

    return exp(a1*pow(THETA,b1)+a2*pow(THETA,b2)+a3*pow(THETA,b3)+a4*pow(THETA,b4)+a5*pow(THETA,b5)+a6*pow(THETA,b6))*459.53;
}

double Viscosity_Trho_R410A(double T, double rho)
{
    // Properties taken from "Viscosity of Mixed 
    // Refrigerants R404A,R407C,R410A, and R507A" 
    // by Vladimir Geller, 
    // 2000 Purdue Refrigeration conferences

    // inputs in T [K], and p [kPa]
    // output in Pa-s

   double eta_microPa_s;

   //Set constants required
   double a_0=-2.695e0,a_1=5.850e-2,a_2=-2.129e-5,b_1=9.047e-3,b_2=5.784e-5,
       b_3=1.309e-7,b_4=-2.422e-10,b_5=9.424e-14,b_6=3.933e-17;

   eta_microPa_s=a_0+a_1*T+a_2*T*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho+b_5*rho*rho*rho*rho*rho+b_6*rho*rho*rho*rho*rho*rho;
   return eta_microPa_s/1e6;
}

double Conductivity_Trho_R410A(double T, double rho)
{
    // Properties taken from "Thermal Conductivity 
    // of the Refrigerant mixtures R404A,R407C,R410A, and R507A" 
    // by V.Z. Geller, B.Z. Nemzer, and U.V. Cheremnykh 
    // Int. J. Thermophysics, v. 22, n 4 2001

    // inputs in T [K], and p [kPa] or rho [kg/m^3]
    // output in W/m-K

    //Set constants required
    double a_0=-8.872e0,a_1=7.410e-2,b_1=3.576e-2,b_2=-9.045e-6,b_3=4.343e-8,b_4=-3.705e-12;

    return (a_0+a_1*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho)/1.e6; // from mW/m-K to kW/m-K
}

// ***************************************************
//                HELMHOLTZ DERIVATIVES
// ***************************************************

double a0_R410A(double tau, double delta)
{
    double sum;
    int k;
    sum=log(delta)-log(tau)+a[0]+a[1]*tau+a[2]*pow(tau,b[2]);
    for(k=3;k<=5;k++)
    {
        sum+=a[k]*log(1.0-exp(-b[k]*tau));
    }
    return sum;
}
double da0_dtau_R410A(double tau, double delta)
{
    double sum;
    int k;
    sum=-1/tau+a[1]+a[2]*b[2]*pow(tau,b[2]-1);
    for(k=3;k<=5;k++)
    {
        sum+=a[k]*b[k]*exp(-b[k]*tau)/(1-exp(-b[k]*tau));
    }
    return sum;
}

double d2a0_dtau2_R410A(double tau, double delta)
{
    double sum;
    int k;
    sum=1.0/(tau*tau)+a[2]*b[2]*(b[2]-1)*pow(tau,b[2]-2);
    for(k=3;k<=5;k++)
    {
        sum+=-a[k]*b[k]*b[k]/(4.0*powInt(sinh(b[k]*tau/2.0),2));
    }
    return sum;
}

double da0_ddelta_R410A(double tau, double delta)
{
    return 1/delta;
}

double d2a0_ddelta2_R410A(double tau, double delta)
{
    return -1/(delta*delta);
}

double ar_R410A(double tau, double delta)
{
    double sum=0;
    int k;
    for  (k=1;k<=21;k++)
    {
        sum+=N[k]*pow(delta,i[k])*pow(tau,j[k])*exp(-g[k]*pow(delta,L[k]));
    }
    return sum;
}
double dar_dtau_R410A(double tau,double delta)
{
    double sum=0;
    int k;

    for  (k=1;k<=21;k++)
    {
        sum+=j[k]*N[k]*pow(delta,i[k])*pow(tau,j[k]-1)*exp(-g[k]*pow(delta,L[k]));
    }
    return sum;
}

double d2ar_dtau2_R410A(double tau, double delta)
{
    double sum=0;
    int k;
    for  (k=1;k<=21;k++)
    {
        sum+=N[k]*pow(delta,i[k])*j[k]*(j[k]-1)*pow(tau,j[k]-2)*exp(-g[k]*pow(delta,L[k]));
    }
    return sum;
}

double d2ar_ddelta2_R410A(double tau,double delta)
{
    double dar2_dDelta2=0;
    int k;
    for  (k=1;k<=21;k++)
    {
        dar2_dDelta2=dar2_dDelta2+N[k]*pow(tau,j[k])*exp(-g[k]*pow(delta,L[k]))*(pow(delta,i[k]-2)*pow(i[k],2)-pow(delta,i[k]-2)*i[k]-2*pow(delta,i[k]-2+L[k])*i[k]*g[k]*L[k]-pow(delta,i[k]-2+L[k])*g[k]*pow(L[k],2)+pow(delta,i[k]-2+L[k])*g[k]*L[k]+pow(delta,i[k]+2*L[k]-2)*pow(g[k],2)*pow(L[k],2));
    }
    return dar2_dDelta2;
}

double dar_ddelta_R410A(double tau,double delta)
{
    double sum=0,gk,ik,Lk;
    int k;
    for  (k=1;k<=21;k++)
    {
        gk=(double)g[k];
        ik=(double)i[k];
        Lk=(double)L[k];
        sum+=N[k]*powInt(delta,i[k]-1)*pow(tau,j[k])*exp(-gk*powInt(delta,L[k]))*(ik-powInt(delta,L[k])*gk*Lk);
    }
    return sum;
}

double d2ar_ddelta_dtau_R410A(double tau,double delta)
{
    double dar_dDelta_dTau=0;
    int k;
    for  (k=1;k<=21;k++)
    {
        dar_dDelta_dTau=dar_dDelta_dTau+N[k]*j[k]*pow(tau,j[k]-1)*(i[k]*pow(delta,i[k]-1)*exp(-g[k]*pow(delta,L[k]))+pow(delta,i[k])*exp(-g[k]*pow(delta,L[k]))*(-g[k]*L[k]*pow(delta,L[k]-1)));
    }
    return dar_dDelta_dTau;
}
