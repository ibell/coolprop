/*
Properties for R134a.  
by Ian Bell

Thermo props from
"A International Standard Formulation for the Thermodynamic Properties of 1,1,1,2-Tetrafluoroethane 
(HFC-134a) for Temperatures from 170 K to 455 K and Pressures up to 70 MPa"
by Reiner Tillner-Roth and Hans Dieter Baehr, J. Phys. Chem. Ref. Data, v. 23, 1994, pp 657-729

In order to call the exposed functions, rho_, h_, s_, cp_,...... there are three different 
ways the inputs can be passed, and this is expressed by the Types integer flag.  
These macros are defined in the PropMacros.h header file:
1) First parameter temperature, second parameter pressure ex: h_R410A(260,1785,1)=-67.53
	In this case, the lookup tables are built if needed and then interpolated
2) First parameter temperature, second parameter density ex: h_R410A(260,43.29,2)=-67.53
	Density and temp plugged directly into EOS
3) First parameter temperature, second parameter pressure ex: h_R410A(260,1785,3)=-67.53
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

#include <math.h>
#include "string.h"
#include "stdio.h"
#include <stdlib.h>
#include "PropErrorCodes.h"
#include "PropMacros.h"
#include "R134a.h"

static int errCode;
static char errStr[ERRSTRLENGTH];

#define nP 100
#define nT 100

static double Tmin=220,Tmax=470, Pmin=24.46,Pmax=1973;
static double hmat[nT][nP];
static double rhomat[nT][nP];
static double cpmat[nT][nP];
static double smat[nT][nP];
static double cvmat[nT][nP];
static double umat[nT][nP];
static double viscmat[nT][nP];
static double kmat[nT][nP];
static double Tvec[nT];
static double pvec[nP];

static int TablesBuilt;

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
static const double R=0.08148885644; //[kJ/kg-K]
// R found from Ru/M, or 8.314471/0.102032/1000

// Function prototypes
static double Pressure_Trho(double T, double rho);
static double IntEnergy_Trho(double T, double rho);
static double Enthalpy_Trho(double T, double rho);
static double Entropy_Trho(double T, double rho);
static double SpecHeatV_Trho(double T, double rho);
static double SpecHeatP_Trho(double T, double rho);
static double Viscosity_Trho(double T, double rho);
static double Conductivity_Trho(double T, double rho);

static double get_Delta(double T, double p);
static double LookupValue(const char *Prop, double T, double p);

static double phi0(double tau,double delta);
static double dphi0_ddelta(double tau,double delta);
static double d2phi0_ddelta2(double tau,double delta);
static double dphi0_dtau(double tau,double delta);
static double d2phi0_dtau2(double tau,double delta);
static double d2phi0_ddelta_dtau(double tau, double delta);

static double phir(double tau,double delta);
static double dphir_ddelta(double tau, double delta);
static double d2phir_ddelta2(double tau, double delta);
static double dphir_dtau(double tau, double delta);
static double d2phir_dtau2(double tau, double delta);
static double d2phir_ddelta_dtau(double tau, double delta);
static double QuadInterp(double x0, double x1, double x2, double f0, double f1, double f2, double x);
static double powInt(double x, int y);

static double dhdT(double tau, double delta);
static double dhdrho(double tau, double delta);
static double dpdT(double tau, double delta);
static double dpdrho(double tau, double delta);

//Microsoft version of math.h doesn't include acosh.h
#if defined(_MSC_VER)
static double acosh(double x)
{
 	return log(x + sqrt(x*x - 1.0) );
}
#endif

static void WriteLookup(void)
{
	int i,j;
	FILE *fp_h,*fp_s,*fp_rho,*fp_u,*fp_cp,*fp_cv,*fp_visc,*fp_k;
	fp_h=fopen("h.csv","w");
	fp_s=fopen("s.csv","w");
	fp_u=fopen("u.csv","w");
	fp_cp=fopen("cp.csv","w");
	fp_cv=fopen("cv.csv","w");
	fp_rho=fopen("rho.csv","w");
	fp_visc=fopen("visc.csv","w");
	fp_k=fopen("k.csv","w");

	// Write the pressure header row
	for (j=0;j<nP;j++)
	{
		fprintf(fp_h,",%0.12f",pvec[j]);
		fprintf(fp_s,",%0.12f",pvec[j]);
		fprintf(fp_rho,",%0.12f",pvec[j]);
		fprintf(fp_u,",%0.12f",pvec[j]);
		fprintf(fp_cp,",%0.12f",pvec[j]);
		fprintf(fp_cv,",%0.12f",pvec[j]);
		fprintf(fp_visc,",%0.12f",pvec[j]);
		fprintf(fp_k,",%0.12f",pvec[j]);
	}
	fprintf(fp_h,"\n");
	fprintf(fp_s,"\n");
	fprintf(fp_rho,"\n");
	fprintf(fp_u,"\n");
	fprintf(fp_cp,"\n");
	fprintf(fp_cv,"\n");
	fprintf(fp_visc,"\n");
	fprintf(fp_k,"\n");
	
	for (i=1;i<nT;i++)
	{
		fprintf(fp_h,"%0.12f",Tvec[i]);
		fprintf(fp_s,"%0.12f",Tvec[i]);
		fprintf(fp_rho,"%0.12f",Tvec[i]);
		fprintf(fp_u,"%0.12f",Tvec[i]);
		fprintf(fp_cp,"%0.12f",Tvec[i]);
		fprintf(fp_cv,"%0.12f",Tvec[i]);
		fprintf(fp_visc,"%0.12f",Tvec[i]);
		fprintf(fp_k,"%0.12f",Tvec[i]);
		for (j=0;j<nP;j++)
		{
			fprintf(fp_h,",%0.12f",hmat[i][j]);
			fprintf(fp_s,",%0.12f",smat[i][j]);
			fprintf(fp_rho,",%0.12f",rhomat[i][j]);
			fprintf(fp_u,",%0.12f",umat[i][j]);
			fprintf(fp_cp,",%0.12f",cpmat[i][j]);
			fprintf(fp_cv,",%0.12f",cvmat[i][j]);
			fprintf(fp_visc,",%0.12f",viscmat[i][j]);
			fprintf(fp_k,",%0.12f",kmat[i][j]);
		}
		fprintf(fp_h,"\n");
		fprintf(fp_s,"\n");
		fprintf(fp_rho,"\n");
		fprintf(fp_u,"\n");
		fprintf(fp_cp,"\n");
		fprintf(fp_cv,"\n");
		fprintf(fp_visc,"\n");
		fprintf(fp_k,"\n");
	}
	fclose(fp_h);
	fclose(fp_s);
	fclose(fp_rho);
	fclose(fp_u);
	fclose(fp_cp);
	fclose(fp_cv);
	fclose(fp_visc);
	fclose(fp_k);
}
static void BuildLookup(void)
{
	int i,j;

	/*      Supercritical
	||X	X X X X	X X X X X X X X X X X X X X X X X
	||X	X X X X	X X X X X X X X X X X X X X X X X
	||			--------  X X X X X X X X X X X X
	||		  /			 \  X X X X X X X X X X X
p	||		 /			  | X X X X Superheated Gas
	||		/	Two  	  / X X X X X X X X X X X
	||	   /	Phase	 /X X X X X X X X X X X X 
	||    /				/ X X X X X X X X X X X X 
	||=	= = = = = = = = = = = = = = = = = = = = =
						Enthalpy										
	*/
	if (!TablesBuilt)
	{
		printf("Building Lookup Tables... Please wait...");

		for (i=0;i<nT;i++)
		{
			Tvec[i]=Tmin+i*(Tmax-Tmin)/(nT-1);
		}
		for (j=0;j<nP;j++)
		{
			pvec[j]=Pmin+j*(Pmax-Pmin)/(nP-1);
		}
		for (i=0;i<nT;i++)
		{
			for (j=0;j<nP;j++)
			{
				if (Tvec[i]>Tc || pvec[j]<psat_R134a(Tvec[i]))
				{					
					rhomat[i][j]=get_Delta(Tvec[i],pvec[j])*rhoc;
					hmat[i][j]=h_R134a(Tvec[i],rhomat[i][j],TYPE_Trho);
					smat[i][j]=s_R134a(Tvec[i],rhomat[i][j],TYPE_Trho);
					umat[i][j]=u_R134a(Tvec[i],rhomat[i][j],TYPE_Trho);
					cpmat[i][j]=cp_R134a(Tvec[i],rhomat[i][j],TYPE_Trho);
					cvmat[i][j]=cv_R134a(Tvec[i],rhomat[i][j],TYPE_Trho);
					viscmat[i][j]=visc_R134a(Tvec[i],rhomat[i][j],TYPE_Trho);
					kmat[i][j]=k_R134a(Tvec[i],rhomat[i][j],TYPE_Trho);
				}
				else
				{
					hmat[i][j]=_HUGE;
					smat[i][j]=_HUGE;
					umat[i][j]=_HUGE;
					rhomat[i][j]=_HUGE;
					umat[i][j]=_HUGE;
					cpmat[i][j]=_HUGE;
					cvmat[i][j]=_HUGE;
					viscmat[i][j]=_HUGE;
					kmat[i][j]=_HUGE;
				}
			}
		}
		TablesBuilt=1;
		printf("Tables Built\n");
		//WriteLookup();
	}
}

int errCode_R134a(void)
{
	return errCode;
}

double rho_R134a(double T, double p, int Types)
{
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_TPNoLookup:
			return get_Delta(T,p)*rhoc;
		case TYPE_TP:
			BuildLookup();
			return LookupValue("rho",T,p);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double p_R134a(double T, double rho)
{
	errCode=0; // Reset error code
	return Pressure_Trho(T,rho);
}
double h_R134a(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return Enthalpy_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoc;
			return Enthalpy_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("h",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double s_R134a(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return Entropy_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoc;
			return Entropy_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("s",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double u_R134a(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return IntEnergy_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoc;
			return IntEnergy_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("u",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double cp_R134a(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return SpecHeatP_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoc;
			return SpecHeatP_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("cp",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double cv_R134a(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return SpecHeatV_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoc;
			return SpecHeatV_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("cv",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double visc_R134a(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return Viscosity_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoc;
			return Viscosity_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("visc",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double k_R134a(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return Conductivity_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoc;
			return Conductivity_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("k",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}

double w_R134a(double T, double p_rho, int Types)
{
	double rho;
	double delta,tau,c1,c2;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho; break;
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoc; break;
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
	
	delta=rho/rhoc;
	tau=Tc/T;

	c1=-SpecHeatV_Trho(T,rho)/R;
	c2=(1.0+2.0*delta*dphir_ddelta(tau,delta)+delta*delta*d2phir_ddelta2(tau,delta));
    return sqrt(-c2*T*SpecHeatP_Trho(T,rho)*1000/c1);
}

double pcrit_R134a(void)
{
	return pc;
}
double Tcrit_R134a(void)
{
	return Tc;
}
double MM_R134a(void)
{
	return M;
}
double psat_R134a(double T)
{
	double theta,phi;
	errCode=0; // Reset error code

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

static double Pressure_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return R*T*rho*(1.0+delta*dphir_ddelta(tau,delta));
}
static double IntEnergy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return R*T*tau*(dphi0_dtau(tau,delta)+dphir_dtau(tau,delta));
}
static double Enthalpy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return R*T*(1.0+tau*(dphi0_dtau(tau,delta)+dphir_dtau(tau,delta))+delta*dphir_ddelta(tau,delta));
}
static double Entropy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return R*(tau*(dphi0_dtau(tau,delta)+dphir_dtau(tau,delta))-phi0(tau,delta)-phir(tau,delta));
}
static double SpecHeatV_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return -R*tau*tau*(d2phi0_dtau2(tau,delta)+d2phir_dtau2(tau,delta));
}
static double SpecHeatP_Trho(double T, double rho)
{
	double delta,tau,c1,c2;
	delta=rho/rhoc;
	tau=Tc/T;

	c1=1.0+delta*dphir_ddelta(tau,delta)-delta*tau*d2phir_ddelta_dtau(tau,delta);
	c2=1.0+2.0*delta*dphir_ddelta(tau,delta)+delta*delta*d2phir_ddelta2(tau,delta);
	
	return R*(-tau*tau*(d2phi0_dtau2(tau,delta)+d2phir_dtau2(tau,delta))+c1*c1/c2);
}
static double Viscosity_Trho(double T, double rho)
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
static double Conductivity_Trho(double T, double rho)
{
	/* 
	From "A multiparameter thermal conductivity equation
	for R134a with an optimized functional form"
	by G. Scalabrin, P. Marchi, F. Finezzo, 
	Fluid Phase Equilibria 245 (2006) 37–51 
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

double dhdT_R134a(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho;
			return dhdT(Tc/T,rho/rhoc);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoc;
			return dhdT(Tc/T,rho/rhoc);
		case 99:
			rho=get_Delta(T,p_rho)*rhoc;
			return (Enthalpy_Trho(T+0.001,rho)-Enthalpy_Trho(T,rho))/0.001;
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}

double dhdrho_R134a(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho;
			return dhdrho(Tc/T,rho/rhoc);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoc;
			return dhdrho(Tc/T,rho/rhoc);
		case 99:
			rho=get_Delta(T,p_rho)*rhoc;
			return (Enthalpy_Trho(T,rho+0.001)-Enthalpy_Trho(T,rho))/0.001;
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}

double dpdT_R134a(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho;
			return dpdT(Tc/T,rho/rhoc);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoc;
			return dpdT(Tc/T,rho/rhoc);
		case 99:
			rho=get_Delta(T,p_rho)*rhoc;
			return (Pressure_Trho(T+0.01,rho)-Pressure_Trho(T,rho))/0.01;
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}

double dpdrho_R134a(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho;
			return dpdrho(Tc/T,rho/rhoc);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoc;
			return dpdrho(Tc/T,rho/rhoc);
		case 99:
			rho=get_Delta(T,p_rho)*rhoc;
			return (Pressure_Trho(T,rho+0.0001)-Pressure_Trho(T,rho))/0.0001;
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}


/**************************************************/
/*           Property Derivatives                 */
/**************************************************/
// See Lemmon, 2000 for more information
static double dhdrho(double tau, double delta)
{
	double T;
	T=Tc/tau;
	//Note: dphi02_dDelta_dTau(tau,delta) is equal to zero
	return R*T/rhoc*(tau*(d2phir_ddelta_dtau(tau,delta))+dphir_ddelta(tau,delta)+delta*d2phir_ddelta2(tau,delta));
}
static double dhdT(double tau, double delta)
{
	double dhdT_rho,T,dhdtau;
	T=Tc/tau;   
	dhdT_rho=R*tau*(dphi0_dtau(tau,delta)+dphir_dtau(tau,delta))+R*delta*dphir_ddelta(tau,delta)+R;
	dhdtau=R*T*(dphi0_dtau(tau,delta)+ dphir_dtau(tau,delta))+R*T*tau*(d2phi0_dtau2(tau,delta)+d2phir_dtau2(tau,delta))+R*T*delta*d2phir_ddelta_dtau(tau,delta);
	return dhdT_rho+dhdtau*(-Tc/T/T);
}
static double dpdT(double tau, double delta)
{
	double T,rho;
	T=Tc/tau;   rho=delta*rhoc;
	return rho*R*(1+delta*dphir_ddelta(tau,delta)-delta*tau*d2phir_ddelta_dtau(tau,delta));
}
static double dpdrho(double tau, double delta)
{
	double T,rho;
	T=Tc/tau;   rho=delta*rhoc;
	return R*T*(1+2*delta*dphir_ddelta(tau,delta)+delta*delta*d2phir_ddelta2(tau,delta));
}


//**********************************************
//                 Derivatives
//**********************************************



static double phi0(double tau,double delta)
{
	return a0[1]+a0[2]*tau+a0[3]*log(tau)+log(delta)+a0[4]*pow(tau,-1.0/2.0)+a0[5]*pow(tau,-3.0/4.0);
}

static double phir(double tau, double delta)
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

static double dphi0_ddelta(double tau,double delta)
{
	return 1.0/delta;
}

static double dphi0_dtau(double tau,double delta)
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
static double d2phi0_ddelta2(double tau,double delta)
{
	return -1.0/(delta*delta);
}
static double d2phi0_dtau2(double tau,double delta)
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
static double d2phi0_ddelta_dtau(double tau, double delta)
{
	return 0.0;
}

static double dphir_ddelta(double tau, double delta)
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
static double dphir_dtau(double tau, double delta)
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
static double d2phir_ddelta2(double tau, double delta)
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

static double d2phir_dtau2(double tau, double delta)
{
	double sum=0,d_,k_;
	int i,k;

	for (i=1;i<=N[0];i++)
	{
		sum += a[i]*t[i]*(t[i]-1.0)*pow(tau,t[i]-2.0)*powInt(delta,d[i]);
	}

	for (k=1;k<=4;k++)
	{
		for (i=N[k-1]+1;i<=N[k];i++)
		{
			d_=(double)d[i];
			k_=(double)k;
			sum += exp(-powInt(delta,k))*a[i]*t[i]*(t[i]-1.0)*pow(tau,t[i]-2.0)*powInt(delta,d[i]);
		}
	}
	return sum;
}

static double d2phir_ddelta_dtau(double tau, double delta)
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

/* 
Utility functions 
*/

static double get_Delta(double T, double p)
{
    double change,eps=1e-8, tau,delta_guess;
    int counter=1;
    double r1,r2,r3,delta1,delta2,delta3;
    
	
	if (T>Tc || p<=psat_R134a(T)) 
	// short-circuits, so avoids problems with 
	// psat not being defined for super-critical temps
    {
		// Superheated vapor and/or supercritical
        delta_guess=p/(R*T)/rhoc;
    }
    else
    {
		// Subcooled liquid
        delta_guess=10;
    }
	change=999;
	tau=Tc/T;
    delta1=delta_guess;
    delta2=delta_guess+.001;
    r1=p/(delta1*rhoc*R*T)-1.0-delta1*dphir_ddelta(tau,delta1);
    r2=p/(delta2*rhoc*R*T)-1.0-delta2*dphir_ddelta(tau,delta2);
    while(counter==1 || fabs(change)>eps)
    {
        delta3=delta2-r2/(r2-r1)*(delta2-delta1);
        r3=p/(delta3*rhoc*R*T)-1.0-delta3*dphir_ddelta(tau,delta3);
        change=r2/(r2-r1)*(delta2-delta1);
        delta1=delta2;
        delta2=delta3;
        r1=r2;
        r2=r3;
        counter=counter+1;
        /*mexPrintf("%g \t %g \t %g \t %g \t %g \n",delta1,r1,delta2,r2,change);*/
    }
	//printf("T: %g p: %g rho: %g\n",T,p,delta3*rhoc);
    return delta3;
}

static double LookupValue(const char *Prop, double T, double p)
{
	int iPlow, iPhigh, iTlow, iThigh,L,R,M;
	double T1, T2, T3, P1, P2, P3, y1, y2, y3, a1, a2, a3;
	double (*mat)[nT][nP];

	if (T>Tmax || T<Tmin)
	{
		errCode=OUT_RANGE_T;
	}
	if (p>Pmax || p<Pmin)
	{
		errCode=OUT_RANGE_P;
	}
	if (T>Tmax || T<Tmin || p>Pmax ||p<Pmin)
	{
		printf("Input to LookupValue() for %s is out of bounds [T:%g p:%g]\n",Prop,T,p);
		return -1e6;
	}
	

	L=0;
	R=nP-1;
	M=(L+R)/2;
	// Use interval halving to find the indices which bracket the temperature of interest
	while (R-L>1)
	{
		if (T>=Tvec[M])
		{ L=M; M=(L+R)/2; }
		if (T<Tvec[M])
		{ R=M; M=(L+R)/2; }
	}
	iTlow=L; iThigh=R;

	L=0;
	R=nP-1;
	M=(L+R)/2;
	// Use interval halving to find the indices which bracket the pressure of interest
	while (R-L>1)
	{
		if (p>=pvec[M])
		{ L=M; M=(L+R)/2; }
		if (p<pvec[M])
		{ R=M; M=(L+R)/2; }
	}
	iPlow=L; iPhigh=R;

	/* Depending on which property is desired, 
	make the matrix "mat" a pointer to the 
	desired property matrix */
	if (!strcmp(Prop,"rho"))
		mat=&rhomat;
	if (!strcmp(Prop,"cp"))
		mat=&cpmat;
	if (!strcmp(Prop,"cv"))
		mat=&cvmat;
	if (!strcmp(Prop,"h"))
		mat=&hmat;
	if (!strcmp(Prop,"s"))
		mat=&smat;
	if (!strcmp(Prop,"u"))
		mat=&umat;
	if (!strcmp(Prop,"visc"))
		mat=&viscmat;
	if (!strcmp(Prop,"k"))
		mat=&kmat;
	
	//At Low Temperature Index
	y1=(*mat)[iTlow][iPlow];
	y2=(*mat)[iTlow][iPhigh];
	y3=(*mat)[iTlow][iPhigh+1];
	P1=pvec[iPlow];
	P2=pvec[iPhigh];
	P3=pvec[iPhigh+1];
	a1=QuadInterp(P1,P2,P3,y1,y2,y3,p);

	//At High Temperature Index
	y1=(*mat)[iThigh][iPlow];
	y2=(*mat)[iThigh][iPhigh];
	y3=(*mat)[iThigh][iPhigh+1];
	a2=QuadInterp(P1,P2,P3,y1,y2,y3,p);

	//At High Temperature Index+1 (for QuadInterp() )
	y1=(*mat)[iThigh+1][iPlow];
	y2=(*mat)[iThigh+1][iPhigh];
	y3=(*mat)[iThigh+1][iPhigh+1];
	a3=QuadInterp(P1,P2,P3,y1,y2,y3,p);

	//At Final Interpolation
	T1=Tvec[iTlow];
	T2=Tvec[iThigh];
	T3=Tvec[iThigh+1];
	return QuadInterp(T1,T2,T3,a1,a2,a3,T);
	
}

/* Helper functions */



static double QuadInterp(double x0, double x1, double x2, double f0, double f1, double f2, double x)
{
    double L0, L1, L2;
    L0=((x-x1)*(x-x2))/((x0-x1)*(x0-x2));
    L1=((x-x0)*(x-x2))/((x1-x0)*(x1-x2));
    L2=((x-x0)*(x-x1))/((x2-x0)*(x2-x1));
    return L0*f0+L1*f1+L2*f2;
}

static double powInt(double x, int y)
{
    int i;
    double product=1.0;
    double x_in;
    int y_in;
    
    if (y==0)
    {
        return 1.0;
    }
    
    if (y<0)
    {
        x_in=1/x;
        y_in=-y;
    }
	else
	{
		x_in=x;
		y_in=y;
	}

    if (y_in==1)
    {
        return x_in;
    }    
    
    product=x_in;
    for (i=1;i<y_in;i++)
    {
        product=product*x_in;
    }
    
    return product;
}
