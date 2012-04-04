/*
Properties for R717 (ammonia).  
by Ian Bell

Thermo props from
"Eine neue Fundamentalgleichung fur Ammoniak (A new Equation of State for Ammonia)"
by R. Tillner-Roth and F. Harms-Watzenberg and H.D. Baehr, Deutscher Kaelte- und Klimatechnischer Verein Tagung 1993

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

/*
From REFPROP documentation:
The original paper has a typographical error that shows a positive coefficient
instead of negative.  The correct value should be -0.3497111e-01.
*/

static const double a[]={
	 0.0,			//[0]
	 0.4554431e-1, 	//[1]
	 0.7238548e+0,	//[2]
	 0.1229470e-1,	//[3]
	-0.1858814e+1,	//[4]
	 0.2141882e-10,	//[5]
	-0.1430020e-1,	//[6]
	 0.3441324e+0,	//[7]
	-0.2873571e+0,	//[8]
	 0.2352589e-4,	//[9]
	-0.3497111e-1,	//[10]
	 0.2397852e-1,	//[11]
	 0.1831117e-2,	//[12]
	-0.4085375e-1,	//[13]
	 0.2379275e+0,	//[14]
 	-0.3548972e-1,	//[15]
	-0.1823729e+0,	//[16]
	 0.2281556e-1,	//[17]
	-0.6663444e-2,	//[18]
	-0.8847486e-2,	//[19]
	 0.2272635e-2,	//[20]
	-0.5588655e-3,	//[21]
};

static const int d[]={
	0,			//[0]
	2, 			//[1]
	1, 			//[2]
	4, 			//[3]
	1, 			//[4]
	15, 		//[5]
	3, 			//[6]
	3, 			//[7]
	1, 			//[8]
	8, 			//[9]
	2, 			//[10]
	1, 			//[11]
	8, 			//[12]
	1, 			//[13]
	2, 			//[14]
	3, 			//[15]
	2, 			//[16]
	4, 			//[17]
	3, 			//[18]
	1, 			//[19]
	2, 			//[20]
	4 			//[21]
};

static const double t[]={
	0.0,		//[0]
	-1.0/2.0,	//[1]
	1.0/2.0,	//[2]
	1.0,		//[3]
	3.0/2.0,	//[4]
	3.0,		//[5]
	0.0,		//[6]
	3.0,		//[7]
	4.0,		//[8]
	4.0,		//[9]
	5.0,		//[10]
	3.0,		//[11]
	5.0,		//[12]
	6.0, 		//[13]
	8.0,		//[14]
	8.0,		//[15]
	10.0,		//[16]
	10.0,		//[17]
	5.0,		//[18]
	15.0/2.0,	//[19]
	15.0,		//[20]
	30.0		//[21]
};

static const int N[]={
	5,			//[0]
	10,			//[1]
	17,			//[2]
	21,			//[3]
	21			//[4]
};

static const double a0[]={
	0.0,		//[0]
	-15.815020,	//[1]
	4.255726,	//[2]
	11.474340,	//[3]
	-1.296211,	//[4]
	0.5706757	//[5]
};
static const double t0[]={
	0.0,		//[0]
	0.0,		//[1]
	0.0,		//[2]
	1.0/3.0,	//[3]
	-3.0/2.0,	//[4]
	-7.0/4.0	//[5]
};

static const double M=17.03; //[kg/kmol]
static const double Tc=405.40; //[K]
static const double rhoc=225; //[kg/m^3]
static const double pc=11333; //[kPa]s
static const double _Ttriple=195.5; //[K]

int Load_R717(struct fluidParamsVals *Fluid)
{
    // Function pointers
    Fluid->funcs.phir=phir_R717;
    Fluid->funcs.dphir_dDelta=dphir_dDelta_R717;
    Fluid->funcs.dphir2_dDelta2=dphir2_dDelta2_R717;
    Fluid->funcs.dphir2_dDelta_dTau=dphir2_dDelta_dTau_R717;
    Fluid->funcs.dphir_dTau=dphir_dTau_R717;
    Fluid->funcs.dphir2_dTau2=dphir2_dTau2_R717;
    Fluid->funcs.phi0=phi0_R717;
    Fluid->funcs.dphi0_dDelta=dphi0_dDelta_R717;
    Fluid->funcs.dphi02_dDelta2=dphi02_dDelta2_R717;
    Fluid->funcs.dphi0_dTau=dphi0_dTau_R717;
    Fluid->funcs.dphi02_dTau2=dphi02_dTau2_R717;
    Fluid->funcs.rhosatL=rhosatL_R717;
    Fluid->funcs.rhosatV=rhosatV_R717;
    Fluid->funcs.psat=psat_R717;

    Fluid->funcs.visc=Viscosity_Trho_R717;
    Fluid->funcs.cond=Conductivity_Trho_R717;

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

double psat_R717(double T)
{
	double theta,phi;

	phi=T/Tc;
	theta=1-phi;

	return pc*exp(5.64633073E-04 - 7.13654961E+00*theta-2.46962841E+00*theta*theta-2.56218797E+01*powInt(theta,3)+4.02270547E+01*powInt(theta,4)-6.72626052E+01*powInt(theta,5));
}

double rhosatL_R717(double T)
{
	double theta, THETA;
	theta=T/Tc;
	THETA=1-theta;
	// |Error| is < 0.1% between 197.6 and 396.2 K
	return exp(1.485*pow(THETA,0.2524)-0.0753)*rhoc;
}

double rhosatV_R717(double T)
{
	double theta, THETA;
	double a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,rhog_hat;
	theta=T/Tc;
	THETA=1-theta;

	// |Error| is < 0.1% between 225 and 375 K (approximately)

	a1 = -0.1912;
	a2 = -5.934;
	a3 = -17.75;  
	a4 = -18.07;  
	a5 = -18.08;  
	a6 = -13.59;  
	b1 = 0.02903; 
	b2 = 0.6403;  
	b3 = 4.829;  
	b4 = 6.766;  
	b5 = 6.767;  
	b6 = 2.333;  
	
	rhog_hat = a1*pow(THETA,b1)+a2*pow(THETA,b2)+a3*pow(THETA,b3)+a4*pow(THETA,b4)+a5*pow(THETA,b5)+a6*pow(THETA,b6);
	return rhoc*exp(rhog_hat);
}

double Viscosity_Trho_R717(double T, double rho)
{
	/* 
	From "The Viscosity of Ammonia"
	by A. Fenghour and W.A. Wakeham and V. Vesovic and J.T.R. Watson and J. Millat and E. Vogel
	J. Phys. Chem. Ref. Data, Vol. 24, No. 5, 1995 
	*/
	double sum=0, e_k=386.0,sigma=0.2957,M=17.03026,sum2=0.0;
	int i,j;
	double T_star,G_eta_star,eta0,B_eta_star,B_eta,b_1,deltaeta_h;
	double a[]={4.99318220,-0.61122364,0.0,0.18535124,-0.11160946};
	double c[]={-0.17999496e1,0.46692621e2,-0.53460794e3,0.33604074e4,-0.13019164e5,
		0.33414230e5,-0.58711743e5,0.71426686e5,-0.59834012e5,0.33652741e5,
		-0.12027350e5,0.24348205e4,-0.20807957e3};
	// indices are backwards from paper
	double d[5][3]={{0,0.17366936e-2,0.0},
	                {0.0,-0.64250359e-2,0.0},
	                {2.19664285e-1,0.0,1.67668649e-4},
	                {0.0,0.0,-1.49710093e-4},
	                {-0.83651107e-1,0.0,0.77012274e-4}};

	// rho is units of mol/L, so convert the density from kg/m^3 to mol/L (poorly documented in paper)
	rho=rho/M;

	sum=0;
	T_star=T/e_k;
	for (i=0;i<=4;i++)
	{
		sum+=a[i]*powInt(log(T_star),i);
	}
	G_eta_star=exp(sum);

	// From REFPROP fluid file: !=0.021357*SQRT(MW)*(unknown factor of 100)  [Chapman-Enskog term]
	// Seems like there is a typo in Fenghour - or am I missing something?
	eta0=0.021357*sqrt(T*M)*100/(sigma*sigma*G_eta_star);
	
	sum=0;
	for (i=0;i<=12;i++)
	{
		sum+=c[i]*powInt(sqrt(T_star),-i);
	}
	B_eta_star=sum;
	B_eta=B_eta_star*(0.6022137*sigma*sigma*sigma);
	b_1=B_eta*eta0;

	sum=0;
	for (i=2;i<=4;i++)
	{
		sum2=0.0;
		for (j=0;j<=4;j++)
		{
			// indices of d are backwards from paper
			sum2+=d[j][i-2]/powInt(T_star,j);
		}
		sum+=sum2*powInt(rho,i);
	}
	deltaeta_h=sum;
	return (eta0+b_1*rho+deltaeta_h)/1e6;
}
double Conductivity_Trho_R717(double T, double rho)
{
	/* 
	From "Thermal Conductivity of Ammonia in a Large 
	Temperature and Pressure Range Including the Critical Region"
	by R. Tufeu, D.Y. Ivanov, Y. Garrabos, B. Le Neindre, 
	Bereicht der Bunsengesellschaft Phys. Chem. 88 (1984) 422-427
	*/

	/* 
	Does not include the critical enhancement.  Comparison of EES (without enhancement) and Refprop (with enhancement) give 
	errors in saturated liquid and saturated vapor conductivities of 

	T < 325K, error < 0.1%
	325K < T < 355 K, error <1%

	Nearly all practical conditions will be in the <325K range
	*/


	double lambda_0,lambda_tilde;

	double a= 0.3589e-1;
	double b=-0.1750e-3;
	double c= 0.4551e-6;
	double d= 0.1685e-9;
	double e=-0.4828e-12;
	double lambda_1= 0.16207e-3;
	double lambda_2= 0.12038e-5;
	double lambda_3=-0.23139e-8;
	double lambda_4= 0.32749e-11;

	// This variable appears to be unused in the code:
	//double a_zeta_plus=0.7;

	double LAMBDA=1.2, nu=0.63, gamma =1.24, DELTA=0.50, rhoc_visc=235,t,zeta_0_plus=1.34e-10,a_zeta=1,GAMMA_0_plus=0.423e-8;
	double pi=3.141592654,a_chi,k_B=1.3806504e-23,X_T,DELTA_lambda,dPdT,eta_B,DELTA_lambda_id,DELTA_lambda_i;

	lambda_0=a+b*T+c*T*T+d*T*T*T+e*T*T*T*T;
	lambda_tilde=lambda_1*rho+lambda_2*rho*rho+lambda_3*rho*rho*rho+lambda_4*rho*rho*rho*rho;
	
	if (0)//(T>Tc)
	{
		t=fabs((T-Tc)/Tc);
		a_chi=a_zeta/0.7;
		eta_B=(2.60*1.6*t)*1e-5;
		dPdT=(2.18-0.12/exp(17.8*t))*1e5; //[Pa-K]
		X_T=0.61*rhoc_visc+16.5*log(t);
		// Along the critical isochore (only a function of temperature) (Eq. 9)
		DELTA_lambda_i=LAMBDA*(k_B*T*T)/(6*pi*eta_B*(zeta_0_plus*pow(t,-nu)*(1+a_zeta*pow(t,DELTA))))*dPdT*dPdT*GAMMA_0_plus*pow(t,-gamma)*(1+a_chi*pow(t,DELTA));
		DELTA_lambda_id=DELTA_lambda_i*exp(-36*t*t);
		if (rho<0.6*rhoc)
		{
			DELTA_lambda=DELTA_lambda_id*(X_T*X_T)/(X_T*X_T+powInt(0.6*rhoc_visc-0.96*rhoc_visc,2))*powInt(rho,2)/powInt(0.6*rhoc_visc,2);
		}
		else
		{
			DELTA_lambda=DELTA_lambda_id*(X_T*X_T)/(X_T*X_T+powInt(rho-0.96*rhoc_visc,2));
		}
	}
	else
	{
		DELTA_lambda=0.0;
	}

	return (lambda_0+lambda_tilde+DELTA_lambda)/1000.0;
}
	
//**********************************************
//                 Derivatives
//**********************************************

double phi0_R717(double tau,double delta)
{
	return log(delta)+a0[1]+a0[2]*tau-log(tau)+a0[3]*pow(tau,1.0/3.0)+a0[4]*pow(tau,-3.0/2.0)+a0[5]*pow(tau,-7.0/4.0);
}

double phir_R717(double tau, double delta)
{
	double sum=0;
	int i;

	for (i=1;i<=5;i++)
	{
		sum += a[i]*pow(tau,t[i])*powInt(delta,d[i]);
	}
	for (i=6;i<=10;i++)
	{
		sum += exp(-delta)*a[i]*pow(tau,t[i])*powInt(delta,d[i]);
	}
	for (i=11;i<=17;i++)
	{
		sum += exp(-delta*delta)*a[i]*pow(tau,t[i])*powInt(delta,d[i]);
	}
	for (i=18;i<=21;i++)
	{
		sum += exp(-delta*delta*delta)*a[i]*pow(tau,t[i])*powInt(delta,d[i]);
	}
	return sum;
}

double dphi0_dDelta_R717(double tau,double delta)
{
	return 1.0/delta;
}

double dphi0_dTau_R717(double tau,double delta)
{
	int j;
	double sum; 
	sum=a0[2]-1/tau;
	for (j=3;j<=5;j++)
	{
		sum+=a0[j]*t0[j]*pow(tau,t0[j]-1.0);
	}
	return sum;
}
double dphi02_dDelta2_R717(double tau,double delta)
{
	return -1.0/(delta*delta);
}
double dphi02_dTau2_R717(double tau,double delta)
{
	int j;
	double sum;
	sum=+1.0/(tau*tau);
	for (j=3;j<=5;j++)
	{
		sum+=a0[j]*t0[j]*(t0[j]-1.0)*pow(tau,t0[j]-2.0);
	}
	return sum;
}
double dphi02_dDelta_dTau_R717(double tau, double delta)
{
	return 0.0;
}

double dphir_dDelta_R717(double tau, double delta)
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
double dphir_dTau_R717(double tau, double delta)
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
double dphir2_dDelta2_R717(double tau, double delta)
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

double dphir2_dTau2_R717(double tau, double delta)
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

double dphir2_dDelta_dTau_R717(double tau, double delta)
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
