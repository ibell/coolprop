/* 
Properties for R407C.  
by Ian Bell (ihb2@cornell.edu)

Pseudo-pure fluid thermo props from 
"Pseudo-pure fluid Equations of State for the Refrigerant Blends R410A, R404A, R507C and R407C" 
by E.W. Lemmon, Int. J. Thermophys. v. 24, n4, 2003

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
    2.13194,		//[0]
    8.05008,		//[1]
    -14.3914,		//[2]
    1.4245,			//[3]
    3.9419,			//[4]
    3.1209			//[5]
};

static const double b[]={
    0,				//[0]
    0,				//[1]
    -0.4,			//[2]
    2.40437,		//[3]
    5.25122,		//[4]
    13.3632			//[5]
};

static const double N[]={
	 0.0,			//[0]
	 1.0588,		//[1]
	-1.12018,		//[2]
	 0.629064,		//[3]
	-0.351953,		//[4]
	 0.00455978,	//[5]
	-1.75725,		//[6]
	-1.12009,		//[7]
	 0.0277353,		//[8]
	 0.898881,		//[9]
	-1.17591,		//[10]
	 0.0818591,		//[11]
	-0.0794097,		//[12]
	-0.0000104047,	//[13]
	 0.233779,		//[14]
	-0.291790,		//[15]
	 0.0154776,		//[16]
	-0.0314579,		//[17]
	-0.00442552,	//[18]
	-0.0101254,		//[19]
	 0.00915953,	//[20]
	-0.003615		//[21]
};

static const double j[]={
	0.0,			//[0]
    0.241,			//[1]
    0.69,			//[2]
    2.58,			//[3]
    1.15,			//[4]
    0.248,			//[5]
    2.15,			//[6]
    2.43,			//[7]
    5.3,			//[8]
    0.76,			//[9]
    1.48,			//[10]
    0.24,			//[11]
    2.86,			//[12]
    8.0,			//[13]
    3.3,			//[14]
    4.7,			//[15]
    0.45,			//[16]
    8.4,			//[17]
    16.2,			//[18]
    26.0,			//[19]
    16.0,			//[20]
    8.7				//[21]
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
    2,				//[8]
    3,				//[9]
    3,				//[10]
    5,				//[11]
    5,				//[12]
    5,				//[13]
    1,				//[14]
    1,				//[15]
    4,				//[16]
    4,				//[17]
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
    1,				//[12]
    1,				//[13]
    2,				//[14]
    2,				//[15]
    2,				//[16]
    2,				//[17]
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
	0.48722,		//[1]
    -6.6959,		//[2]
    -1.4165,		//[3]
    -2.5109			//[4]
};
    
static const double tbp[]={
	0.0,			//[0]
    0.54,			//[1]
    0.925,			//[2]
    2.7,			//[3]
    4.7				//[4]
};

static const double Ndp[]={
    0.0,			//[0]
	-0.086077,		//[1]
    -6.6364,		//[2]
    -2.4648,		//[3]
    -3.4776			//[4]
};

static const double tdp[]={
	0.0,			//[0]
    0.4,			//[1]
    0.965,			//[2]
    3.1,			//[3]
    5.0				//[4]
};

static const double M=86.2036; //[g/mol]
static const double Tm=359.345; //[K]
static const double pm=4631.7; //[MPa--> kPa]
static const double pc=4600; //[MPa--> kPa] From (Calm 2007 HPAC Engineering)
static const double rhom=453.43094; //5.260*M; //[mol/dm^3--> kg/m^3]
static const double _Ttriple=200.0; //[K]

int Load_R407C(struct fluidParamsVals *Fluid)
{
    // Function pointers
    Fluid->funcs.phir=ar_R407C;
    Fluid->funcs.dphir_dDelta=dar_ddelta_R407C;
    Fluid->funcs.dphir2_dDelta2=d2ar_ddelta2_R407C;
    Fluid->funcs.dphir2_dDelta_dTau=d2ar_ddelta_dtau_R407C;
    Fluid->funcs.dphir_dTau=dar_dtau_R407C;
    Fluid->funcs.dphir2_dTau2=d2ar_dtau2_R407C;
    Fluid->funcs.phi0=a0_R407C;
    Fluid->funcs.dphi0_dDelta=da0_ddelta_R407C;
    Fluid->funcs.dphi02_dDelta2=d2a0_ddelta2_R407C;
    Fluid->funcs.dphi0_dTau=da0_dtau_R407C;
    Fluid->funcs.dphi02_dTau2=d2a0_dtau2_R407C;
    Fluid->funcs.rhosatL=rhosatL_R407C;
    Fluid->funcs.rhosatV=rhosatV_R407C;
    Fluid->funcs.p_dp=p_dp_R407C;
    Fluid->funcs.p_bp=p_bp_R407C;
    Fluid->funcs.visc=Viscosity_Trho_R407C;
    Fluid->funcs.cond=Conductivity_Trho_R407C;

    //Lookup table parameters
    Fluid->LUT.Tmin=220.0;
    Fluid->LUT.Tmax=550.0;
    Fluid->LUT.pmin=70.03;
    Fluid->LUT.pmax=5000.0;

    //Fluid parameters
    Fluid->Type=FLUIDTYPE_REFRIGERANT_PSEUDOPURE;
    Fluid->Tc=Tm;
    Fluid->rhoc=rhom;
    Fluid->MM=M;
    Fluid->pc=pc;
    Fluid->Tt=_Ttriple;
    return 1;
}

double p_bp_R407C(double T)
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
    
double p_dp_R407C(double T)
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

double rhosatV_R407C(double T)
{
	double theta;
	theta=1-T/359.345;
	return exp(+2.481666724+2.974063949110e+00*pow(theta,-0.023725)-14.334177672*theta+60.7010706462*theta*theta-376.514657192*powInt(theta,3)+1178.57631747*powInt(theta,4)-1999.37929072*powInt(theta,5)+1307.74729667*powInt(theta,6));
}

double rhosatL_R407C(double T)
{
	double theta;
	theta=1-T/359.345;
	return exp(5.544589249+1.764403419125e+00*pow(theta,0.12337)+0.544950396285*theta-0.784102758738*theta*theta+0.741332715649*theta*theta*theta);
}

double Viscosity_Trho_R407C(double T, double rho)
{
    // Properties taken from "Viscosity of Mixed 
	// Refrigerants R404A,R407C,R410A, and R507A" 
	// by Vladimir Geller, 
	// 2000 Purdue Refrigeration conferences

    // inputs in T [K], and p [kPa]
    // output in Pa-s

   double eta_microPa_s;

   //Set constants required
   double a_0=-1.507e0,a_1=4.894e-2,a_2=-9.305e-6,b_1=-3.038e-3,b_2=2.927e-4,
	   b_3=-9.559e-7,b_4=1.739e-9,b_5=-1.455e-12,b_6=4.756e-16;

   eta_microPa_s=a_0+a_1*T+a_2*T*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho+b_5*rho*rho*rho*rho*rho+b_6*rho*rho*rho*rho*rho*rho;
   return eta_microPa_s/1e6;
}

double Conductivity_Trho_R407C(double T, double rho)
{
	// Properties taken from "Thermal Conductivity 
	// of the Refrigerant mixtures R404A,R407C,R410A, and R507A" 
	// by V.Z. Geller, B.Z. Nemzer, and U.V. Cheremnykh 
	// Int. J. Thermophysics, v. 22, n 4 2001

	// inputs in T [K], and p [kPa] or rho [kg/m^3]
	// output in W/m-K

	//Set constants required
	double a_0=-9.628e0,a_1=7.638e-2,b_1=2.715e-2,b_2=4.963e-5,b_3=-4.912e-8,b_4=2.884e-11;

	return (a_0+a_1*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho)/1.e6; // from mW/m-K to kW/m-K
}

// ***************************************************
//                HELMHOLTZ DERIVATIVES
// ***************************************************

double a0_R407C(double tau, double delta)
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
double da0_dtau_R407C(double tau, double delta)
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

double d2a0_dtau2_R407C(double tau, double delta)
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

double da0_ddelta_R407C(double tau, double delta)
{
	return 1/delta;
}

double d2a0_ddelta2_R407C(double tau, double delta)
{
	return -1/(delta*delta);
}

double ar_R407C(double tau, double delta)
{
	double sum=0;
	int k;
	for  (k=1;k<=21;k++)
    {
        sum+=N[k]*pow(delta,i[k])*pow(tau,j[k])*exp(-g[k]*pow(delta,L[k]));
    }
	return sum;
}
double dar_dtau_R407C(double tau,double delta)
{
	double sum=0;
	int k;

	for  (k=1;k<=21;k++)
    {
        sum+=j[k]*N[k]*pow(delta,i[k])*pow(tau,j[k]-1)*exp(-g[k]*pow(delta,L[k]));
    }
	return sum;
}

double d2ar_dtau2_R407C(double tau, double delta)
{
	double sum=0;
	int k;
	for  (k=1;k<=21;k++)
    {
        sum+=N[k]*pow(delta,i[k])*j[k]*(j[k]-1)*pow(tau,j[k]-2)*exp(-g[k]*pow(delta,L[k]));
    }
	return sum;
}

double d2ar_ddelta2_R407C(double tau,double delta)
{
    double dar2_dDelta2=0;
    int k;
    for  (k=1;k<=21;k++)
    {
        dar2_dDelta2=dar2_dDelta2+N[k]*pow(tau,j[k])*exp(-g[k]*pow(delta,L[k]))*(pow(delta,i[k]-2)*pow((double)i[k],2)-pow(delta,i[k]-2)*i[k]-2*pow(delta,i[k]-2+L[k])*i[k]*g[k]*L[k]-pow(delta,i[k]-2+L[k])*g[k]*pow((double)L[k],2)+pow(delta,i[k]-2+L[k])*g[k]*L[k]+pow(delta,i[k]+2*L[k]-2)*pow((double)g[k],2)*pow((double)L[k],2));
    }
    return dar2_dDelta2;
}

double dar_ddelta_R407C(double tau,double delta)
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

double d2ar_ddelta_dtau_R407C(double tau,double delta)
{
    double dar_dDelta_dTau=0;
    int k;
    for  (k=1;k<=21;k++)
    {
        dar_dDelta_dTau=dar_dDelta_dTau+N[k]*j[k]*pow(tau,j[k]-1)*(i[k]*pow(delta,i[k]-1)*exp(-g[k]*pow(delta,L[k]))+pow(delta,i[k])*exp(-g[k]*pow(delta,L[k]))*(-g[k]*L[k]*pow(delta,L[k]-1)));
    }
    return dar_dDelta_dTau;
}
