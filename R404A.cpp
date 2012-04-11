/*
Properties for R404A.  
by Ian Bell

Pseudo-pure fluid thermo props from 
"Pseudo-pure fluid Equations of State for the Refrigerant Blends R410A, R404A, R507C and R407C" 
by E.W. Lemmon, Int. J. Thermophys. v. 24, n4, 2003

In order to call the exposed functions, rho_, h_, s_, cp_,...... 
there are three different ways the inputs can be passed, and this is expressed by the Types integer flag.  
These macros are defined in the PropMacros.h header file:
1) First parameter temperature, second parameter pressure ex: h_R404A(260,354.7,1)=274
	In this case, the lookup tables are built if needed and then interpolated
2) First parameter temperature, second parameter density ex: h_R404A(260,13.03,2)=274
	Density and temp plugged directly into EOS
3) First parameter temperature, second parameter pressure ex: h_R404A(260,354.7,3)=274
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
    7.00407,		//[0]
    7.98695,		//[1]
    -18.8664,		//[2]
    0.63078,		//[3]
    3.5979,			//[4]
    5.0335			//[5]
};

static const double b[]={
    0,				//[0]
    0,				//[1]
    -0.3,			//[2]
    1.19617,		//[3]
    2.32861,		//[4]
    5.00188			//[5]
};

static const double N[]={
	 0.0,			//[0]
     6.10984,		//[1]
    -7.79453,		//[2]
     0.0183377,		//[3]
	 0.262270,		//[4]
    -0.00351688,	//[5]
     0.0116181,		//[6]
     0.00105992,	//[7]
     0.850922,		//[8]
    -0.520084,		//[9]
    -0.0464225,		//[10]
     0.621190,		//[11]
    -0.195505,		//[12]
     0.336159,		//[13]
    -0.0376062,		//[14]
    -0.00636579,	//[15]
    -0.0758262,		//[16]
    -0.0221041,		//[17]
	 0.0310441,		//[18]
     0.0132798,		//[19]
     0.0689437,		//[20]
    -0.0507525,		//[21]
	 0.0161382,		//[22]
};

static const double j[]={
	0.0,			//[0]
    0.67,			//[1]
    0.91,			//[2]
    5.96,			//[3]
    0.7,			//[4]
    6.0,			//[5]
    0.3,			//[6]
    0.7,			//[7]
    1.7,			//[8]
    3.3,			//[9]
    7.0,			//[10]
    2.05,			//[11]
    4.3,			//[12]
    2.7,			//[13]
    1.8,			//[14]
    1.25,			//[15]
    12.0,			//[16]
    6.0,			//[17]
    8.7,			//[18]
    11.6,			//[19]
    13.0,			//[20]
    17.0,			//[21]
	16.0			//[22]
};

static const int i[]={
	0,				//[0]
    1,				//[1]
    1,				//[2]
    1,				//[3]
    2,				//[4]
    2,				//[5]
    4,				//[6]
    6,				//[7]
    1,				//[8]
    1,				//[9]
    1,				//[10]
    2,				//[11]
    2,				//[12]
    3,				//[13]
    4,				//[14]
    7,				//[15]
    2,				//[16]
    3,				//[17]
    4,				//[18]
    4,				//[19]
    2,				//[20]
    3,				//[21]
	5				//[22]
};

static const int L[]={
    0,				//[0]
	0,				//[1]
    0,				//[2]
    0,				//[3]
    0,				//[4]
    0,				//[5]
    0,				//[6]
    0,				//[7]
    1,				//[8]
    1,				//[9]
    1,				//[10]
    1,				//[11]
    1,				//[12]
    1,				//[13]
    1,				//[14]
    1,				//[15]
    2,				//[16]
    2,				//[17]
    2,				//[18]
    2,				//[19]
    3,				//[20]
    3,				//[21]
	3				//[22]
};

static const int g[]={
	0,				//[0]
	0,				//[1]
	0,				//[2]
	0,				//[3]
	0,				//[4]
	0,				//[5]
	0,				//[6]
	0,				//[7]
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
	1				//[22]
};

static const double Nbp[]={
	0.0,			//[0]
	0.061067,		//[1]
    -6.5646,		//[2]
    -3.6162,		//[3]
    -3.9771			//[4]
};
    
static const double tbp[]={
	0.0,			//[0]
    0.54,			//[1]
    0.965,			//[2]
    3.7,			//[3]
    9.0				//[4]
};

static const double Ndp[]={
    0.0,			//[0]
	-0.00026863,	//[1]
    -6.5757,		//[2]
    -4.1802,		//[3]
    -7.9102			//[4]
};

static const double tdp[]={
	0.0,			//[0]
    0.1,			//[1]
    0.972,			//[2]
    3.8,			//[3]
    9.0				//[4]
};

static const double M=97.6038; //[g/mol]
static const double Tm=345.27; //[K]
static const double pm=3734.8; //[MPa--> kPa]
static const double pc=3720; //[MPa--> kPa] From (Calm 2007 HPAC Engineering)
static const double rhom=482.163; //4.940*M; //[mol/dm^3--> kg/m^3]
static const double _Ttriple=200.0; //[K]

int Load_R404A(struct fluidParamsVals *Fluid)
{
    // Function pointers
    Fluid->funcs.phir=ar_R404A;
    Fluid->funcs.dphir_dDelta=dar_ddelta_R404A;
    Fluid->funcs.dphir2_dDelta2=d2ar_ddelta2_R404A;
    Fluid->funcs.dphir2_dDelta_dTau=d2ar_ddelta_dtau_R404A;
    Fluid->funcs.dphir_dTau=dar_dtau_R404A;
    Fluid->funcs.dphir2_dTau2=d2ar_dtau2_R404A;
    Fluid->funcs.phi0=a0_R404A;
    Fluid->funcs.dphi0_dDelta=da0_ddelta_R404A;
    Fluid->funcs.dphi02_dDelta2=d2a0_ddelta2_R404A;
    Fluid->funcs.dphi0_dTau=da0_dtau_R404A;
    Fluid->funcs.dphi02_dTau2=d2a0_dtau2_R404A;
    Fluid->funcs.rhosatL=rhosatL_R404A;
    Fluid->funcs.rhosatV=rhosatV_R404A;
    Fluid->funcs.p_dp=p_dp_R404A;
    Fluid->funcs.p_bp=p_bp_R404A;
    Fluid->funcs.visc=Viscosity_Trho_R404A;
    Fluid->funcs.cond=Conductivity_Trho_R404A;

    //Lookup table parameters
    Fluid->LUT.Tmin=220.0;
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

double p_bp_R404A(double T)
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
    
double p_dp_R404A(double T)
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

double rhosatV_R404A(double T)
{
	//int counter;
	//double rho1,rho2,rho3,r1,r2,r3,change,eps=1e-8;
	///* Solve for the density which gives the 
	//pressure equal to the dew pressure */
	//counter=1;
	//change=999;
 //   rho1=10;
 //   rho2=10+.001;
 //   r1=p_dp_R404A(T)-Pressure_Trho(T,rho1);
 //   r2=p_dp_R404A(T)-Pressure_Trho(T,rho2);
 //   while(counter==1 || fabs(change)>eps)
 //   {
 //       rho3=rho2-r2/(r2-r1)*(rho2-rho1);
 //       r3=p_dp_R404A(T)-Pressure_Trho(T,rho3);
 //       change=r2/(r2-r1)*(rho2-rho1);
 //       rho1=rho2; rho2=rho3;
 //       r1=r2; r2=r3;
 //       counter=counter+1;
 //   }
 //   return rho3;

	double THETA,a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,b6,b7;
	THETA=1-T/345.27;

	a1 =      -17.61; 
    a2 =      -16.61; 
    a3 =      -6.059; 
    a4 =      -16.43; 
    a5 =      -16.85; 
    a6 =      -16.26; 
    a7 =      -2.926; 
    b1 =       2.939; 
    b2 =       6.538; 
    b3 =       1.144; 
    b4 =       6.536; 
    b5 =        6.54;  
    b6 =       6.533; 
    b7 =      0.3967;

	return exp(a1*pow(THETA,b1)+a2*pow(THETA,b2)+a3*pow(THETA,b3)+a4*pow(THETA,b4)+a5*pow(THETA,b5)+a6*pow(THETA,b6)+a7*pow(THETA,b7))*482.13;
}

double rhosatL_R404A(double T)
{
	//int counter;

	//double rho1,rho2,rho3,r1,r2,r3,change,eps=1e-10;
	///* Solve for the density which gives the 
	//pressure equal to the bubble pressure */
	//counter=1;
	//change=999;
 //   rho1=-2.57878924896E-10*powInt(T,6) + 3.58401823351E-07*powInt(T,5) - 2.04856639109E-04*powInt(T,4) + 6.15631400205E-02*powInt(T,3) - 1.02527033660E+01*powInt(T,2) + 8.94030144564E+02*powInt(T,1) - 3.02025468937E+04;
 //   rho2=rho1+0.1;
 //   r1=p_bp_R404A(T)-Pressure_Trho(T,rho1);
 //   r2=p_bp_R404A(T)-Pressure_Trho(T,rho2);
 //   while(counter==1 || fabs(r1)>eps)
 //   {
 //       rho3=rho2-r2/(r2-r1)*(rho2-rho1);
 //       r3=p_bp_R404A(T)-Pressure_Trho(T,rho3);
 //       change=r2/(r2-r1)*(rho2-rho1);
 //       rho1=rho2; rho2=rho3;
 //       r1=r2; r2=r3;
 //       counter=counter+1;
 //   }
 //   return rho3;

	double THETA,a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,b6,b7;
	THETA=1-T/345.27;

	a1 =      0.3897;
	a2 =      0.2768;
	a3 =       0.546; 
	a4 =      0.3473;
	a5 =     -0.1371;
	a6 =      0.2254;
	a7 =     -0.1426;
	b1 =      0.2581;
	b2 =       2.985;
	b3 =       0.255;
	b4 =      0.2596;
	b5 =       1.526;
	b6 =      0.2661;
	b7 =      0.0742;

	return exp(a1*pow(THETA,b1)+a2*pow(THETA,b2)+a3*pow(THETA,b3)+a4*pow(THETA,b4)+a5*pow(THETA,b5)+a6*pow(THETA,b6)+a7*pow(THETA,b7))*482.163;

}

double Viscosity_Trho_R404A(double T, double rho)
{
    // Properties taken from "Viscosity of Mixed 
	// Refrigerants R404A,R407C,R410A, and R507A" 
	// by Vladimir Geller, 
	// 2000 Purdue Refrigeration conferences

    // inputs in T [K], and p [kPa]
    // output in Pa-s

   double 
	eta_microPa_s;

   //Set constants required
   double a_0=9.766e-1,a_1=3.676e-2,a_2=2.938e-6,b_1=2.260e-3,b_2=1.786e-4,
	   b_3=-4.202e-7,b_4=8.489e-10,b_5=-8.670e-13,b_6=3.566e-16;

   eta_microPa_s=a_0+a_1*T+a_2*T*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho+b_5*rho*rho*rho*rho*rho+b_6*rho*rho*rho*rho*rho*rho;
   return eta_microPa_s/1e6;
}

double Conductivity_Trho_R404A(double T, double rho)
{
	// Properties taken from "Thermal Conductivity 
	// of the Refrigerant mixtures R404A,R407C,R410A, and R507A" 
	// by V.Z. Geller, B.Z. Nemzer, and U.V. Cheremnykh 
	// Int. J. Thermophysics, v. 22, n 4 2001

	// inputs in T [K], and p [kPa] or rho [kg/m^3]
	// output in W/m-K

	//Set constants required
	double a_0=-8.624e0,a_1=7.360e-2,b_1=3.222e-2,b_2=2.569e-5,b_3=-2.693e-8,b_4=2.007e-11;

	return (a_0+a_1*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho)/1.e6; // from mW/m-K to kW/m-K
}

// ***************************************************
//                HELMHOLTZ DERIVATIVES
// ***************************************************

double a0_R404A(double tau, double delta)
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
double da0_dtau_R404A(double tau, double delta)
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

double d2a0_dtau2_R404A(double tau, double delta)
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

double da0_ddelta_R404A(double tau, double delta)
{
	return 1/delta;
}

double d2a0_ddelta2_R404A(double tau, double delta)
{
	return -1/(delta*delta);
}

double ar_R404A(double tau, double delta)
{
	double sum=0;
	int k;
	for  (k=1;k<=22;k++)
    {
        sum+=N[k]*pow(delta,i[k])*pow(tau,j[k])*exp(-g[k]*pow(delta,L[k]));
    }
	return sum;
}
double dar_dtau_R404A(double tau,double delta)
{
	double sum=0;
	int k;

	for  (k=1;k<=22;k++)
    {
        sum+=j[k]*N[k]*pow(delta,i[k])*pow(tau,j[k]-1)*exp(-g[k]*pow(delta,L[k]));
    }
	return sum;
}

double d2ar_dtau2_R404A(double tau, double delta)
{
	double sum=0;
	int k;
	for  (k=1;k<=22;k++)
    {
        sum+=N[k]*pow(delta,i[k])*j[k]*(j[k]-1)*pow(tau,j[k]-2)*exp(-g[k]*pow(delta,L[k]));
    }
	return sum;
}

double d2ar_ddelta2_R404A(double tau,double delta)
{
    double dar2_dDelta2=0;
    int k;
    for  (k=1;k<=22;k++)
    {
        dar2_dDelta2=dar2_dDelta2+N[k]*pow(tau,j[k])*exp(-g[k]*pow(delta,L[k]))*(pow(delta,i[k]-2)*pow((double)i[k],2)-pow(delta,i[k]-2)*i[k]-2*pow(delta,i[k]-2+L[k])*i[k]*g[k]*L[k]-pow(delta,i[k]-2+L[k])*g[k]*pow((double)L[k],2)+pow(delta,i[k]-2+L[k])*g[k]*L[k]+pow(delta,i[k]+2*L[k]-2)*pow((double)g[k],2)*pow((double)L[k],2));
    }
    return dar2_dDelta2;
}

double dar_ddelta_R404A(double tau,double delta)
{
    double sum=0,gk,ik,Lk;
    int k;
    for  (k=1;k<=22;k++)
    {
		gk=(double)g[k];
		ik=(double)i[k];
		Lk=(double)L[k];
		sum+=N[k]*powInt(delta,i[k]-1)*pow(tau,j[k])*exp(-gk*powInt(delta,L[k]))*(ik-powInt(delta,L[k])*gk*Lk);
    }
    return sum;
}

double d2ar_ddelta_dtau_R404A(double tau,double delta)
{
    double dar_dDelta_dTau=0;
    int k;
    for  (k=1;k<=22;k++)
    {
        dar_dDelta_dTau=dar_dDelta_dTau+N[k]*j[k]*pow(tau,j[k]-1)*(i[k]*pow(delta,i[k]-1)*exp(-g[k]*pow(delta,L[k]))+pow(delta,i[k])*exp(-g[k]*pow(delta,L[k]))*(-g[k]*L[k]*pow(delta,L[k]-1)));
    }
    return dar_dDelta_dTau;
}
