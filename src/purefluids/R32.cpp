/*
Properties for R32.  
by Ian Bell

???
Thermo props from
"A International Standard Formulation for the Thermodynamic Properties of 1,1,1,2-Tetrafluoroethane 
(HFC-134a) for Temperatures from 170 K to 455 K and Pressures up to 70 MPa"
by Reiner Tillner-Roth and Hans Dieter Baehr, J. Phys. Chem. Ref. Data, v. 23, 1994, pp 657-729

*/

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include <math.h>
#include "string.h"
#include "stdio.h"
#include <stdlib.h>
#include "CoolProp.h"

static const double a[]={
	 0.0,			//[0]
	 0.1046634e+1, 	//[1]
	-0.5451165,		//[2]
	-0.2448595e-2,	//[3]
	-0.4877002e-1,	//[4]
	 0.3520158e-1,	//[5]
	 0.1622750e-2,	//[6]
	 0.2377225e-4,	//[7]
	 0.29149e-1,	//[8]
	 0.3386203e-2,	//[9]
	-0.4202444e-2,	//[10]
	 0.4782025e-3,	//[11]
	-0.5504323e-2,	//[12]
	-0.2418396e-1,	//[13]
	 0.4209034,		//[14]
 	-0.4616537,		//[15]
	-0.1200513e+1,	//[16]
	-0.2591550e+1,	//[17]
	-0.1400145e+1,	//[18]
	 0.8263017		//[19]
};

static const int d[]={
	0,			//[0]
	1, 			//[1]
	2, 			//[2]
	5, 			//[3]
	1, 			//[4]
	1, 			//[5]
	3, 			//[6]
	8, 			//[7]
	4, 			//[8]
	4, 			//[9]
	4, 			//[10]
	8, 			//[11]
	3, 			//[12]
	5, 			//[13]
	1, 			//[14]
	1, 			//[15]
	3, 			//[16]
	1, 			//[17]
	2, 			//[18]
	3 			//[19]
};

static const double t[]={
	0.0,		//[0]
	1.0/4.0,	//[1]
	1.0,		//[2]
	-1.0/4.0,	//[3]
	-1.0,		//[4]
	2.0,		//[5]
	2.0,		//[6]
	3.0/4.0,	//[7]
	1.0/4.0,	//[8]
	18.0,		//[9]
	26.0,		//[10]
	-1.0,		//[11]
	25.0, 		//[12]
	7.0/4.0,	//[13]
	4.0,		//[14]
	5.0,		//[15]
	1.0,		//[16]
	3.0/2.0,	//[17]
	1.0,		//[18]
	1.0/2.0		//[19]
};

static const int e[]={
	0,			//[0]
	0,			//[1]
	0,			//[2]
	0,			//[3]
	0,			//[4]
	0,			//[5]
	0,			//[6]
	0,			//[7]
	0,			//[8]
	4,			//[9]
	3,			//[10]
	1,			//[11]
	4,			//[12]
	1,			//[13]
	2,			//[14]
	2,			//[15]
	1,			//[16]
	1,			//[17]
	1,			//[18]
	1			//[19]
};

static const double a0[]={
	-8.258096,	//[0]
	6.353098,	//[1]
	3.004486,	//[2]
	1.160761,	//[3]
	2.645151,	//[4]
	5.794987,	//[5]
	1.129475	//[6]
};
static const double n0[]={
	0.0,		//[0]
	0.0,		//[1]
	0.0,		//[2]
	2.2718538,	//[3]
	11.9144210,	//[4]
	5.1415638,	//[5]
	32.7682170	//[6]
};

static const double M=52.024; //[kg/kmol]
static const double Tc=351.255; //[K]
static const double rhoc=424; //[kg/m^3]
static const double pc=5784; //[kPa]
static const double _Ttriple=136.34; //[K]

int Load_R32(struct fluidParamsVals *Fluid)
{
    // Function pointers
    Fluid->funcs.phir=phir_R32;
    Fluid->funcs.dphir_dDelta=dphir_ddelta_R32;
    Fluid->funcs.dphir2_dDelta2=d2phir_ddelta2_R32;
    Fluid->funcs.dphir2_dDelta_dTau=d2phir_ddelta_dtau_R32;
    Fluid->funcs.dphir_dTau=dphir_dtau_R32;
    Fluid->funcs.dphir2_dTau2=d2phir_dtau2_R32;
    Fluid->funcs.phi0=phi0_R32;
    Fluid->funcs.dphi0_dDelta=dphi0_ddelta_R32;
    Fluid->funcs.dphi02_dDelta2=d2phi0_ddelta2_R32;
    Fluid->funcs.dphi0_dTau=dphi0_dtau_R32;
    Fluid->funcs.dphi02_dTau2=d2phi0_dtau2_R32;
    Fluid->funcs.rhosatL=rhosatL_R32;
    Fluid->funcs.rhosatV=rhosatV_R32;
    Fluid->funcs.psat=psat_R32;

    Fluid->funcs.visc=Viscosity_Trho_R32;
    Fluid->funcs.cond=Conductivity_Trho_R32;

    //Lookup table parameters
    Fluid->LUT.Tmin=220.0;
    Fluid->LUT.Tmax=470.0;
    Fluid->LUT.pmin=24.46;
    Fluid->LUT.pmax=1973;

    //Fluid parameters
    Fluid->Type=FLUIDTYPE_REFRIGERANT_PURE;
    Fluid->Tc=Tc;
    Fluid->rhoc=rhoc;
    Fluid->MM=M;
    Fluid->pc=pc;
    Fluid->Tt=_Ttriple;
    return 1;
}

double psat_R32(double T)
{
	double theta,phi,p0=5781.16;

	phi=T/Tc;
	theta=1-phi;

	return p0*exp((-7.44892*theta+1.6886*pow(theta,3.0/2.0)-1.908*pow(theta,5.0/2.0)-2.810*theta*theta*theta*theta*theta)/phi);
}

double rhosatL_R32(double T)
{
	double theta, phi;
	phi=T/Tc;
	theta=1-phi;

	return 424.0+434.55*pow(theta,1.0/4.0)+1296.53*pow(theta,2.0/3.0)-777.49*theta+366.84*pow(theta,5.0/3.0);
}

double rhosatV_R32(double T)
{
	double theta, phi;
	phi=T/Tc;
	theta=1-phi;

	return rhoc*exp(-1.969*pow(theta,1.0/3.0)-2.0222*pow(theta,2.0/3.0)-6.7409*pow(theta,4.0/3.0)-27.479*pow(theta,11.0/3.0));
}

double Viscosity_Trho_R32(double T, double rho)
{
	//ERROR
	fprintf(stderr,"Viscosity_Trho for R32 not coded");
	return _HUGE;
}
double Conductivity_Trho_R32(double T, double rho)
{
	//ERROR
	fprintf(stderr,"Conductivity_Trho for R32 not coded");
	return _HUGE;
}
	
//**********************************************
//                 Derivatives
//**********************************************

double phi0_R32(double tau,double delta)
{
	return log(delta)+a0[0]+a0[1]*tau+a0[2]*log(tau)+a0[3]*log(1-exp(-n0[3]*tau))+a0[4]*log(1-exp(-n0[4]*tau))+a0[5]*log(1-exp(-n0[5]*tau))+a0[6]*log(1-exp(-n0[6]*tau));
}

double dphi0_ddelta_R32(double tau,double delta)
{
	return 1.0/delta;
}

double dphi0_dtau_R32(double tau,double delta)
{
	int j;
	double sum; 
	sum=a0[1]+a0[2]/tau;
	for (j=3;j<=6;j++)
	{
		sum+=a0[j]*n0[j]/(exp(n0[j]*tau)-1.0);
	}
	return sum;
}
double d2phi0_ddelta2_R32(double tau,double delta)
{
	return -1.0/(delta*delta);
}
double d2phi0_dtau2_R32(double tau,double delta)
{
	int j;
	double sum;
	sum=-a0[2]/(tau*tau);
	for (j=3;j<=6;j++)
	{
		sum -= a0[j]*n0[j]*n0[j]*exp(-n0[j]*tau)/powInt(1.0-exp(-n0[j]*tau),2);
	}
	return sum;
}
double d2phi0_ddelta_dtau_R32(double tau, double delta)
{
	return 0.0;
}

double phir_R32(double tau, double delta)
{
	double sum=0;
	int i;

	for (i=1;i<=8;i++)
	{
		sum += a[i]*pow(tau,t[i])*powInt(delta,d[i]);
	}
	for (i=9;i<=19;i++)
	{
		sum += exp(-powInt(delta,e[i]))*a[i]*pow(tau,t[i])*powInt(delta,d[i]);
	}
	return sum;
}

double dphir_ddelta_R32(double tau, double delta)
{
	double sum=0;
	int i;

	for (i=1;i<=8;i++)
	{
		sum += a[i]*((double)d[i])*powInt(delta,d[i]-1)*pow(tau,t[i]);
	}
	for (i=9;i<=19;i++)
	{
		sum += a[i]*exp(-powInt(delta,e[i]))*((double)d[i]-e[i]*powInt(delta,e[i]))*powInt(delta,d[i]-1)*pow(tau,t[i]);
	}
	return sum;
}
double dphir_dtau_R32(double tau, double delta)
{
	double sum=0;
	int i;

	for (i=1;i<=8;i++)
	{
		sum += a[i]*t[i]*powInt(delta,d[i])*pow(tau,t[i]-1.0);
	}
	for (i=9;i<=19;i++)
	{
		sum += a[i]*t[i]*exp(-powInt(delta,e[i]))*powInt(delta,d[i])*pow(tau,t[i]-1.0);
	}
	return sum;
}
double d2phir_ddelta2_R32(double tau, double delta)
{
	double sum=0,di,ei; 
	int i;
	for (i=1;i<=8;i++)
	{
		di=(double)d[i];
		sum += a[i]*di*(di-1.0)*powInt(delta,d[i]-2)*pow(tau,t[i]);
	}
	for (i=9;i<=19;i++)
	{
		di=(double)d[i];
		ei=(double)e[i];
		sum += a[i]* exp(-powInt(delta,e[i])) *(di*di-di-ei*powInt(delta,e[i])*(2.0*di+ei-1.0-ei*powInt(delta,e[i])))*powInt(delta,d[i]-2)*pow(tau,t[i]);
	}
	return sum;
}

double d2phir_dtau2_R32(double tau, double delta)
{
	double sum=0;
	int i;

	for (i=1;i<=8;i++)
		sum += a[i]*t[i]*(t[i]-1.0)*pow(tau,t[i]-2.0)*powInt(delta,d[i]);
	for (i=9;i<=19;i++)
	{
		sum += a[i]*t[i]*(t[i]-1.0)*exp(-powInt(delta,e[i]))*powInt(delta,d[i])*pow(tau,t[i]-2.0);
	}
	return sum;
}

double d2phir_ddelta_dtau_R32(double tau, double delta)
{
	double sum=0,di,ei;
	int i;

	for (i=1;i<=8;i++)
	{
		di=(double)d[i];
		sum += a[i]*di*t[i]*powInt(delta,d[i]-1)*pow(tau,t[i]-1.0);
	}
	for (i=9;i<=19;i++)
	{
		di=(double)d[i];
		ei=(double)e[i];
		sum += a[i]*t[i]*(di-ei*powInt(delta,e[i]))*exp(-powInt(delta,e[i]))*powInt(delta,d[i]-1)*pow(tau,t[i]-1.0);
	}
	return sum;
}
