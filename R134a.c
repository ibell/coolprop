/*
Properties for R134a.  
by Ian Bell

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
#include "CoolProp.h"

static const double a[]={
	 0.0,			//[0]
	 0.5586817e-1, 	//[1]
	 0.4982230e0,	//[2]
	 0.2458698e-1,	//[3]
	 0.8570145e-3,	//[4]
	 0.4788584e-3,	//[5]
	-0.1800808e1,	//[6]
	 0.2671641e0,	//[7]
	-0.4781652e-1,	//[8]
	 0.1423987e-1,	//[9]
	 0.3324062e0,	//[10]
	-0.7485907e-2,	//[11]
	 0.1017263e-3,	//[12]
	-0.5184567e+0,	//[13]
	-0.8692288e-1,	//[14]
 	 0.2057144e+0,	//[15]
	-0.5000457e-2,	//[16]
	 0.4603262e-3,	//[17]
	-0.3497836e-2,	//[18]
	 0.6995038e-2,	//[19]
	-0.1452184e-1,	//[20]
	-0.1285458e-3,	//[21]
};

static const int d[]={
	0,			//[0]
	2, 			//[1]
	1, 			//[2]
	3, 			//[3]
	6, 			//[4]
	6, 			//[5]
	1, 			//[6]
	1, 			//[7]
	2, 			//[8]
	5, 			//[9]
	2, 			//[10]
	2, 			//[11]
	4, 			//[12]
	1, 			//[13]
	4, 			//[14]
	1, 			//[15]
	2, 			//[16]
	4, 			//[17]
	1, 			//[18]
	5, 			//[19]
	3, 			//[20]
	10 			//[21]
};

static const double t[]={
	0.0,		//[0]
	-1.0/2.0,	//[1]
	0.0,		//[2]
	0.0,		//[3]
	0.0,		//[4]
	3.0/2.0,	//[5]
	3.0/2.0,	//[6]
	2.0,		//[7]
	2.0,		//[8]
	1.0,		//[9]
	3.0,		//[10]
	5.0,		//[11]
	1.0,		//[12]
	5.0, 		//[13]
	5.0,		//[14]
	6.0,		//[15]
	10.0,		//[16]
	10.0,		//[17]
	10.0,		//[18]
	18.0,		//[19]
	22.0,		//[20]
	50.0		//[21]
};

static const int N[]={
	8,			//[0]
	11,			//[1]
	17,			//[2]
	20,			//[3]
	21			//[4]
};

static const double a0[]={
	0.0,		//[0]
	-1.019535,	//[1]
	9.047135,	//[2]
	-1.629789,	//[3]
	-9.723916,	//[4]
	-3.927170	//[5]
};
static const double t0[]={
	0.0,		//[0]
	0.0,		//[1]
	0.0,		//[2]
	0.0,		//[3]
	-1.0/2.0,	//[4]
	-3.0/4.0	//[5]
};

static const double M=102.032; //[kg/kmol]
static const double Tc=374.18; //[K]
static const double rhoc=508; //[kg/m^3]
static const double pc=4056.29; //[kPa]
static const double Tt=169.85; //[K]

//Microsoft version of math.h doesn't include acosh.h
#if defined(_MSC_VER)
static double acosh(double x)
{
 	return log(x + sqrt(x*x - 1.0) );
}
#endif

int Load_R134a(struct fluidParamsVals *Fluid)
{
    // Function pointers
    Fluid->funcs.phir=phir_R134a;
    Fluid->funcs.dphir_dDelta=dphir_dDelta_R134a;
    Fluid->funcs.dphir2_dDelta2=dphir2_dDelta2_R134a;
    Fluid->funcs.dphir2_dDelta_dTau=dphir2_dDelta_dTau_R134a;
    Fluid->funcs.dphir_dTau=dphir_dTau_R134a;
    Fluid->funcs.dphir2_dTau2=dphir2_dTau2_R134a;
    Fluid->funcs.phi0=phi0_R134a;
    Fluid->funcs.dphi0_dDelta=dphi0_dDelta_R134a;
    Fluid->funcs.dphi02_dDelta2=dphi02_dDelta2_R134a;
    Fluid->funcs.dphi0_dTau=dphi0_dTau_R134a;
    Fluid->funcs.dphi02_dTau2=dphi02_dTau2_R134a;
    Fluid->funcs.rhosatL=rhosatL_R134a;
    Fluid->funcs.rhosatV=rhosatV_R134a;
    Fluid->funcs.psat=psat_R134a;
    
    Fluid->funcs.visc=Viscosity_Trho_R134a;
    Fluid->funcs.cond=Conductivity_Trho_R134a;

    //Lookup table parameters
    Fluid->LUT.Tmin=220.0;
    Fluid->LUT.Tmax=470.0;
    Fluid->LUT.pmin=24;
    Fluid->LUT.pmax=1973;
    
    //Fluid parameters
    Fluid->Type=FLUIDTYPE_REFRIGERANT_PURE;
    Fluid->Tc=Tc;
    Fluid->rhoc=rhoc;
    Fluid->MM=M;
    Fluid->pc=pc;
    Fluid->Tt=Tt;
    return 1;
}

double psat_R134a(double T)
{
	double theta,phi;
	
	phi=T/374.18;
	theta=1-phi;

	return pc*exp((-7.686556*theta+2.311791*pow(theta,1.5)-2.039554*theta*theta-3.583758*theta*theta*theta*theta)/phi);
}

double rhosatL_R134a(double T)
{
	double theta, THETA;
	theta=T/374.15;
	THETA=1-theta;

	return 518.20+884.13*pow(THETA,1.0/3.0)+485.84*pow(THETA,2.0/3.0)+193.29*pow(THETA,10.0/3.0);
}

double rhosatV_R134a(double T)
{
	double theta, THETA,rho0=516.86;
	theta=T/374.15;
	THETA=1-theta;

	return rho0*exp(-2.837294*pow(THETA,1.0/3.0)-7.875988*pow(THETA,2.0/3.0)+4.478586*pow(THETA,1.0/2.0)-14.140125*pow(THETA,9.0/4.0)-52.361297*pow(THETA,11.0/2.0));
}

double Viscosity_Trho_R134a(double T, double rho)
{
	/* 
	From "A Reference Multiparameter Viscosity Equation for R134a
	with an Optimized Functional Form"
	by G. Scalabrin and P. Marchi, R. Span
	J. Phys. Chem. Ref. Data, Vol. 35, No. 2, 2006 
	*/
	double sum=0, Tr,rhor;
	int i;
	double g[]={0.0,0.0,0.0,1.0,1.0,2.0,2.0,4.0,0.0,2.0,5.0};
	double h[]={0.0,2.0,20.0,0.0,3.0,0.0,4.0,14.0,1.0,1.0,3.0};
	double n[]={0.0,0.6564868,0.6882417e-10,0.5668165,-0.2989820,-0.1061795,
		0.6245080e-1,0.2758366e-6,-0.1621088,0.1675102,-0.9224693e-1};

	Tr=T/Tc;
	rhor=rho/rhoc;

	for (i=1;i<=7;i++)
	{
		sum += n[i]*pow(Tr,g[i])*pow(rhor,h[i]);
	}
	for (i=8;i<=10;i++)
	{
		sum += exp(-2*rhor*rhor)*n[i]*pow(Tr,g[i])*pow(rhor,h[i]);
	}
	return (exp(sum)-1.0)*25.17975/1e6;
}
double Conductivity_Trho_R134a(double T, double rho)
{
	/* 
	From "A multiparameter thermal conductivity equation
	for R134a with an optimized functional form"
	by G. Scalabrin, P. Marchi, F. Finezzo, 
	Fluid Phase Equilibria 245 (2006) 37�51 
	*/
	int i;
	double sum=0, Tr,rhor,alpha,lambda_r_ce,lambda_r,num,den;
	double g[]={0.0,0.0,0.5,1.0,1.5,4.0,5.5,6.0,0.0};
	double h[]={0.0,1.0,1.0,6.0,0.0,3.0,0.0,0.0,1.0};
	double n[]={0.0,23.504800,-15.261689,0.064403724,7.9735850,0.28675949,
		8.1819842,-6.4852285,-4.8298888};
	double nc=1.2478242;
	double a[]={0.0,1.0,0.0,0.0,0.30,0.30,0.36525,
		0.61221,0.94930,0.92162,0.15,0.08,0.14};

	Tr=T/Tc;
	rhor=rho/rhoc;
	alpha=1.0-a[10]*acosh(1+a[11]*pow((1-Tr)*(1-Tr),a[12]));
	num=rhor*exp(-pow(rhor,a[1])/a[1]-powInt(a[2]*(Tr-1.0),2)-powInt(a[3]*(rhor-1.0),2));
	den=pow(pow(powInt((1.0-1.0/Tr)+a[4]*pow((rhor-1.0)*(rhor-1.0),1.0/(2.0*a[5])),2),a[6])+pow(a[7]*a[7]*(rhor-alpha)*(rhor-alpha),a[8]),a[9]);
	lambda_r_ce=num/den;
	for(i=1;i<=7;i++)
	{
		sum+=n[i]*pow(Tr,g[i])*pow(rhor,h[i]);
	}
	lambda_r=sum+n[8]*exp(-5.0*rhor*rhor)*pow(Tr,g[8])*pow(rhor,h[8])+nc*lambda_r_ce;
	return 2.0547*lambda_r/1e6;
}


//**********************************************
//                 Derivatives
//**********************************************

double phi0_R134a(double tau,double delta)
{
	return a0[1]+a0[2]*tau+a0[3]*log(tau)+log(delta)+a0[4]*pow(tau,-1.0/2.0)+a0[5]*pow(tau,-3.0/4.0);
}

double phir_R134a(double tau, double delta)
{
	double sum=0;
	int i;

	for (i=1;i<=8;i++)
	{
		sum += a[i]*pow(tau,t[i])*powInt(delta,d[i]);
	}
	for (i=9;i<=11;i++)
	{
		sum += exp(-delta)*a[i]*pow(tau,t[i])*powInt(delta,d[i]);
	}
	for (i=12;i<=17;i++)
	{
		sum += exp(-delta*delta)*a[i]*pow(tau,t[i])*powInt(delta,d[i]);
	}
	for (i=18;i<=20;i++)
	{
		sum += exp(-delta*delta*delta)*a[i]*pow(tau,t[i])*powInt(delta,d[i]);
	}
	sum += exp(-delta*delta*delta*delta)*a[21]*pow(tau,t[21])*powInt(delta,d[21]);
	return sum;
}

double dphi0_dDelta_R134a(double tau,double delta)
{
	return 1.0/delta;
}

double dphi0_dTau_R134a(double tau,double delta)
{
	int j;
	double sum; 
	sum=a0[2]+a0[3]/tau;
	for (j=4;j<=5;j++)
	{
		sum+=a0[j]*t0[j]*pow(tau,t0[j]-1.0);
	}
	return sum;
}
double dphi02_dDelta2_R134a(double tau,double delta)
{
	return -1.0/(delta*delta);
}
double dphi02_dTau2_R134a(double tau,double delta)
{
	int j;
	double sum;
	sum=-a0[3]/(tau*tau);
	for (j=4;j<=5;j++)
	{
		sum+=a0[j]*t0[j]*(t0[j]-1.0)*pow(tau,t0[j]-2.0);
	}
	return sum;
}
double dphi02_dDelta_dtau_R134a(double tau, double delta)
{
	return 0.0;
}

double dphir_dDelta_R134a(double tau, double delta)
{
	double sum=0;
	int i,k;

	for (i=1;i<=N[0];i++)
	{
		sum += a[i]*((double)d[i])*pow(tau,t[i])*powInt(delta,d[i]-1);
	}
	for (k=1;k<=4;k++)
	{
		for (i=N[k-1]+1;i<=N[k];i++)
		{
			sum += exp(-powInt(delta,k))*a[i]*((double)d[i]-k*powInt(delta,k))*pow(tau,t[i])*powInt(delta,d[i]-1);
		}
	}
	return sum;
}
double dphir_dTau_R134a(double tau, double delta)
{
	double sum=0;
	int i,k;

	for (i=1;i<=N[0];i++)
	{
		sum += a[i]*t[i]*pow(tau,t[i]-1.0)*powInt(delta,d[i]);
	}

	for (k=1;k<=4;k++)
	{
		for (i=N[k-1]+1;i<=N[k];i++)
		{
			sum += exp(-powInt(delta,k))*a[i]*t[i]*pow(tau,t[i]-1.0)*powInt(delta,d[i]);
		}
	}
	return sum;
}
double dphir2_dDelta2_R134a(double tau, double delta)
{
	double sum=0,d_,k_;
	int i,k;

	for (i=1;i<=N[0];i++)
	{
		d_=(double)d[i];
		sum += a[i]*d_*(d_-1.0)*pow(tau,t[i])*powInt(delta,d[i]-2);
	}

	for (k=1;k<=4;k++)
	{
		for (i=N[k-1]+1;i<=N[k];i++)
		{
			d_=(double)d[i];
			k_=(double)k;
			sum += exp(-powInt(delta,k))*a[i]*(d_*d_-d_-k_*powInt(delta,k)*(2.0*d_+k_-1.0-k_*powInt(delta,k)))*pow(tau,t[i])*powInt(delta,d[i]-2);
		}
	}
	return sum;
}

double dphir2_dTau2_R134a(double tau, double delta)
{
	double sum=0;
	int i,k;

	for (i=1;i<=N[0];i++)
	{
		sum += a[i]*t[i]*(t[i]-1.0)*pow(tau,t[i]-2.0)*powInt(delta,d[i]);
	}

	for (k=1;k<=4;k++)
	{
		for (i=N[k-1]+1;i<=N[k];i++)
		{
			sum += exp(-powInt(delta,k))*a[i]*t[i]*(t[i]-1.0)*pow(tau,t[i]-2.0)*powInt(delta,d[i]);
		}
	}
	return sum;
}

double dphir2_dDelta_dTau_R134a(double tau, double delta)
{
	double sum=0,d_,k_;
	int i,k;

	for (i=1;i<=N[0];i++)
	{
		d_=(double)d[i];
		sum += a[i]*t[i]*d_*pow(tau,t[i]-1.0)*powInt(delta,d[i]-1);
	}

	for (k=1;k<=4;k++)
	{
		for (i=N[k-1]+1;i<=N[k];i++)
		{
			d_=(double)d[i];
			k_=(double)k;
			sum += exp(-powInt(delta,k))*a[i]*t[i]*(d_-k_*powInt(delta,k))*pow(tau,t[i]-1.0)*powInt(delta,d[i]-1);
		}
	}
	return sum;
}
