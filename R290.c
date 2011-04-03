/* Properties of Propane (R290)
by Ian Bell

Themo properties from 
Miyamoto, H., and Watanabe, K., "A thermodynamic property model for fluid-phase propane,"
Int. J. Thermophys., 21(5):1045-1072, 2000.  

In order to call the exposed functions, rho_, h_, s_, cp_,...... there are three 
different ways the inputs can be passed, and this is expressed by the Types integer flag.  
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

#include "math.h"
#include "stdio.h"
#include <string.h>
#include "R290.h"
#include "PropMacros.h"
#include "PropErrorCodes.h"

static int errCode;
static char errStr[ERRSTRLENGTH];

#define nP 200
#define nT 200

static double hmat[nT][nP];
static double rhomat[nT][nP];
static double cpmat[nT][nP];
static double smat[nT][nP];
static double cvmat[nT][nP];
static double umat[nT][nP];
static double viscmat[nT][nP];
static double Tvec[nT];
static double pvec[nP];

static int TablesBuilt;

static const double c[]={
	2.8749,			//[0]
    2.0623,			//[1]
    5.9751,			//[2]
	1.5612			//[3]
};

static const double e[]={
    0.1,			//[0]
    697,			//[1]
    1723.0,			//[2]
    3875.0			//[3]
};

static const double a0[]={
    0.0,			//[0]
	-4.992402,		//[1]
    4.291476,		//[2]
    3.021394,		//[3]
    2.889980,		//[4]
    4.474243,		//[5]
    8.139803,		//[6]
	10.48251		//[7]
};

static const double n[]={
    0.0,			//[0]
    0.0,			//[1]
    0.0,			//[2]
    0.0,			//[3]
    1.048309,		//[4]    
	3.053170,		//[5]
	11.42280,		//[6]
	5.042815		//[7]
};

static const double a[]={
	 0.0,			//[0]
     2.698378e-1,	//[1]
    -1.339252e+0,	//[2]
    -2.273858e-2,	//[3]
     2.414973e-1,	//[4]
    -3.321461e-2,	//[5]
     2.203323e-3,	//[6]
     5.935588e-5,	//[7]
    -1.137457e-6,	//[8]
    -2.379299e+0,	//[9]
     2.337373e+0,	//[10]
     1.242344e-3,	//[11]
    -7.352787e-3,	//[12]
     1.965751e-3,	//[13]
    -1.402666e-1,	//[14]
    -2.093360e-2,	//[15]
    -2.475221e-4,	//[16]
    -1.482723e-2,	//[17]
    -1.303038e-2,	//[18]
     3.634670e-5	//[19]
};

static const double t[]={
	0.0,			//[0]
    -0.25,			//[1]
    1.5,			//[2]
    -0.75,			//[3]
    0.0,			//[4]
    1.25,			//[5]
    1.5,			//[6]
    0.5,			//[7]
    2.5,			//[8]
    1.5,			//[9]
    1.75,			//[10]
    -0.25,			//[11]
    3.0,			//[12]
    3.0,			//[13]
    4.0,			//[14]
    2.0,			//[15]
    -1.0,			//[16]
    2.0,			//[17]
    19.0,			//[18]
    5.0				//[19]
};

static const int d[]={
	0,				//[0]
    1,				//[1]
    1,				//[2]
    2,				//[3]
    2,				//[4]
    3,				//[5]
    5,				//[6]
    8,				//[7]
    8,				//[8]
    3,				//[9]
    3,				//[10]
    8,				//[11]
    5,				//[12]
    6,				//[13]
    1,				//[14]
    5,				//[15]
    7,				//[16]
    2,				//[17]
    3,				//[18]
    15				//[19]
};

static const int N[]={
	8,
	13,
	16,
	19
};

static const double Asat[]={
	 0.0,			//[0]
	-6.741653,		//[1]
     1.455497,		//[2]
    -1.312986,		//[3]
    -2.111039		//[4]
};

static const double Bsat[]={
	 0.0,			//[0]
	 0.5312985,		//[1]
    -1.702073,		//[2]
    -4.998449,		//[3]
   -12.18881,		//[4]
   -42.75035,		//[5]
  -107.8777			//[6]
};

static const double Csat[]={
	 0.0,			//[0]
	 0.2758388,		//[1]
	 1.810924,		//[2]
    -0.8907309,		//[3]
     0.1273854		//[4]    
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




static const double R=0.188555507; // 8.314472 / 44.09562;
static const double M=44.09562; //[g/mol]
static const double Tc=369.825; //[K]
static const double pc=4247.09; //[MPa--> kPa]
static const double rhoc=218.5; //[kg/m^3]


// Local function prototypes
static double Pressure_Trho(double T, double rho);
static double IntEnergy_Trho(double T, double rho);
static double Enthalpy_Trho(double T, double rho);
static double Entropy_Trho(double T, double rho);
static double SpecHeatV_Trho(double T, double rho);
static double SpecHeatP_Trho(double T, double rho);

static double get_Delta(double T, double p);
static double LookupValue(char *Prop,double T, double p);

static double phi0(double tau, double delta);
static double dphi0_dtau(double tau, double delta);
static double d2phi0_dtau2(double tau, double delta);
static double dphi0_ddelta(double tau, double delta);

static double phir(double tau, double delta);
static double dphir_dtau(double tau,double delta);
static double d2phir_dtau2(double tau, double delta);
static double dphir_ddelta(double tau,double delta);
static double d2phir_ddelta2(double tau,double delta);
static double d2phir_ddelta_dtau(double tau,double delta);

static double powInt(double x, int y);
static double QuadInterp(double x0, double x1, double x2, double f0, double f1, double f2, double x);


static void WriteLookup(void)
{
	int i,j;
	FILE *fp_h,*fp_s,*fp_rho,*fp_u,*fp_cp,*fp_cv,*fp_visc;
	fp_h=fopen("h.csv","w");
	fp_s=fopen("s.csv","w");
	fp_u=fopen("u.csv","w");
	fp_cp=fopen("cp.csv","w");
	fp_cv=fopen("cv.csv","w");
	fp_rho=fopen("rho.csv","w");
	fp_visc=fopen("visc.csv","w");

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
	}
	fprintf(fp_h,"\n");
	fprintf(fp_s,"\n");
	fprintf(fp_rho,"\n");
	fprintf(fp_u,"\n");
	fprintf(fp_cp,"\n");
	fprintf(fp_cv,"\n");
	fprintf(fp_visc,"\n");
	
	for (i=1;i<nT;i++)
	{
		fprintf(fp_h,"%0.12f",Tvec[i]);
		fprintf(fp_s,"%0.12f",Tvec[i]);
		fprintf(fp_rho,"%0.12f",Tvec[i]);
		fprintf(fp_u,"%0.12f",Tvec[i]);
		fprintf(fp_cp,"%0.12f",Tvec[i]);
		fprintf(fp_cv,"%0.12f",Tvec[i]);
		fprintf(fp_visc,"%0.12f",Tvec[i]);
		for (j=0;j<nP;j++)
		{
			fprintf(fp_h,",%0.12f",hmat[i][j]);
			fprintf(fp_s,",%0.12f",smat[i][j]);
			fprintf(fp_rho,",%0.12f",rhomat[i][j]);
			fprintf(fp_u,",%0.12f",umat[i][j]);
			fprintf(fp_cp,",%0.12f",cpmat[i][j]);
			fprintf(fp_cv,",%0.12f",cvmat[i][j]);
			fprintf(fp_visc,",%0.12f",viscmat[i][j]);
		}
		fprintf(fp_h,"\n");
		fprintf(fp_s,"\n");
		fprintf(fp_rho,"\n");
		fprintf(fp_u,"\n");
		fprintf(fp_cp,"\n");
		fprintf(fp_cv,"\n");
		fprintf(fp_visc,"\n");
	}
	fclose(fp_h);
	fclose(fp_s);
	fclose(fp_u);
	fclose(fp_cp);
	fclose(fp_cv);
	fclose(fp_visc);
	fclose(fp_rho);
}
static void BuildLookup(void)
{
	int i,j;
	double Tmin=220,Tmax=340;
	double Pmin=70.03,Pmax=3327;

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
				if (pvec[j]<psat_R290(Tvec[i]) || pvec[j]>psat_R290(Tvec[i]))
				{					
					rhomat[i][j]=get_Delta(Tvec[i],pvec[j])*rhoc;
					hmat[i][j]=h_R290(Tvec[i],rhomat[i][j],TYPE_Trho);
					smat[i][j]=s_R290(Tvec[i],rhomat[i][j],TYPE_Trho);
					umat[i][j]=u_R290(Tvec[i],rhomat[i][j],TYPE_Trho);
					cpmat[i][j]=cp_R290(Tvec[i],rhomat[i][j],TYPE_Trho);
					cvmat[i][j]=cv_R290(Tvec[i],rhomat[i][j],TYPE_Trho);
					viscmat[i][j]=visc_R290(Tvec[i],rhomat[i][j],TYPE_Trho);
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
				}
			}
		}
		TablesBuilt=1;
		printf("Tables Built!\n");
		//WriteLookup();
	}
}

int errCode_R290(void)
{
	return errCode;
}

double rho_R290(double T, double p_rho, int Types)
{
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_TPNoLookup:
			return get_Delta(T,p_rho)*rhoc;
		case TYPE_TP:
			BuildLookup();
			return LookupValue("rho",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double p_R290(double T, double rho)
{
	return Pressure_Trho(T,rho);
}
double h_R290(double T, double p_rho, int Types)
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
double s_R290(double T, double p_rho, int Types)
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
double u_R290(double T, double p_rho, int Types)
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
double cp_R290(double T, double p_rho, int Types)
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
double cv_R290(double T, double p_rho, int Types)
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

double w_R290(double T, double p_rho, int Types)
{
	double rho;
	double delta,tau,c1,c2;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho;
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoc;
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


double pcrit_R290(void)
{
	return pc;
}

double Tcrit_R290(void)
{
	return Tc;
}


double MM_R290(void)
{
	return M;
}
double psat_R290(double T)
{
	double sum=0,x,Tr;
	Tr=T/Tc;
    x=1-Tr;
    
    sum=Asat[1]*x+Asat[2]*pow(x,1.5)+Asat[3]*pow(x,2.5)+Asat[4]*pow(x,4.5);
    return pc*exp(1.0/(1.0-x)*sum);
}

double rhosatV_R290(double T)
{
	double sum=0,x,Tr;
	Tr=T/Tc;
    x=1-Tr;
    
    sum=Bsat[1]*pow(x,0.1)+Bsat[2]*pow(x,0.2)+Bsat[3]*pow(x,0.8)+Bsat[4]*pow(x,2.4)+Bsat[5]*pow(x,5.8)+Bsat[6]*pow(x,13.9);
    return rhoc*exp(sum);
}


double rhosatL_R290(double T)
{
	double sum=0,x,Tr;
	Tr=T/Tc;
    x=1-Tr;
    
    sum=Csat[1]*pow(x,0.2)+Csat[2]*pow(x,0.4)+Csat[3]*pow(x,0.6)+Csat[4]*pow(x,1.8);
    return rhoc*exp(sum);
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


static double phi0(double tau, double delta)
{
	double sum;
	int k;
	sum=log(delta)+a0[1]+a0[2]*tau+a0[3]*log(tau);
	for(k=4;k<=7;k++)
	{
		sum+=a0[k]*log(1.0-exp(-n[k]*tau));
	}
	return sum;
}
static double dphi0_dtau(double tau, double delta)
{
	double sum;
	int k;
	sum=a0[2]+a0[3]/tau;
	for(k=4;k<=7;k++)
	{
		sum+=a0[k]*n[k]*exp(-n[k]*tau)/(1-exp(-n[k]*tau));
	}
	return sum;
}

static double d2phi0_dtau2(double tau, double delta)
{
	double sum;
	int k;
	sum=-a0[3]/(tau*tau);
	for(k=4;k<=7;k++)
	{
		sum+=-a0[k]*n[k]*n[k]/(4.0*powInt(sinh(n[k]*tau/2.0),2));
	}
	return sum;
}

static double dphi0_ddelta(double tau, double delta)
{
	return 1/delta;
}

static double phir(double tau, double delta)
{
	double sum=0;
	int i;

	for (i=1;i<=8;i++)
	{
		sum += a[i]*pow(tau,t[i])*powInt(delta,d[i]);
	}
	for (i=9;i<=13;i++)
	{
		sum += exp(-delta)*a[i]*pow(tau,t[i])*powInt(delta,d[i]);
	}
	for (i=14;i<=16;i++)
	{
		sum += exp(-delta*delta)*a[i]*pow(tau,t[i])*powInt(delta,d[i]);
	}
	for (i=17;i<=19;i++)
	{
		sum += exp(-delta*delta*delta)*a[i]*pow(tau,t[i])*powInt(delta,d[i]);
	}
	return sum;
}

static double dphir_ddelta(double tau, double delta)
{
	double sum=0;
	int i,k;

	for (i=1;i<=N[0];i++)
	{
		sum += a[i]*((double)d[i])*pow(tau,t[i])*powInt(delta,d[i]-1);
	}
	for (k=1;k<=3;k++)
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

	for (k=1;k<=3;k++)
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
	/* wxMaxima code:
	diff(diff(a[i]*tau^t[i]*delta^d[i]*exp(-delta^k),delta),delta);
	factor(%);
	*/
	for (k=1;k<=3;k++)
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

	for (k=1;k<=3;k++)
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

	/* 
	wxMaxima code:
	diff(diff(a[i]*tau^t[i]*delta^d[i]*exp(-delta^k),delta),tau);
	factor(%);
	*/
	for (k=1;k<=3;k++)
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

double visc_R290(double T, double p_rho, int Types)
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

	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rhor=p_rho/rhoc; break;
		case TYPE_TPNoLookup:
			rhor=get_Delta(T,p_rho); break;
		case TYPE_TP:
			rhor=get_Delta(T,p_rho); break;
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}

	//Reduced Temperature
	Tr=T/Tc;

	for(i=1;i<=11;i++)
	{
		sum+=                    nv[i]*pow(Tr,tv[i])*pow(rhor,dv[i]);
	}
	for(i=12;i<=13;i++)
	{
		sum+=exp(-rhor*rhor/2.0)*nv[i]*pow(Tr,tv[i])*pow(rhor,dv[i]);
	}
	etar=sum;

	return (exp(etar)-1)*Hc/1e6;
}

double k_R290(double T, double T_rho, int Types)
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
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho: 
			delta=T_rho/rhoc; break;
		case (TYPE_TP | TYPE_TPNoLookup):
			delta=get_Delta(T,T_rho); break;
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
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

static double get_Delta(double T, double p)
{
    double change,eps=1e-8, tau,delta_guess,ps;
    int counter=1;
    double r1,r2,r3,delta1,delta2,delta3;
    
	if (T>Tc)
	{
	delta_guess=p/(R*T)/rhoc;
	}
	else
	{
		if (p<=psat_R290(T))
	    {
			// Superheated vapor
	        delta_guess=p/(R*T)/rhoc;
	    }
	    else
	    {
			// Subcooled liquid
	        delta_guess=10;
	    }
	}
	ps=psat_R290(T);
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
    return delta3;
}

static double LookupValue(char *Prop, double T, double p)
{
	int iPlow, iPhigh, iTlow, iThigh,L,R,M;
	double T1, T2, T3, P1, P2, P3, y1, y2, y3, a1, a2, a3;
	double (*mat)[nT][nP];

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
	make the matrix mat a pointer to the 
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

static double QuadInterp(double x0, double x1, double x2, double f0, double f1, double f2, double x)
{
    double L0, L1, L2;
    L0=((x-x1)*(x-x2))/((x0-x1)*(x0-x2));
    L1=((x-x0)*(x-x2))/((x1-x0)*(x1-x2));
    L2=((x-x0)*(x-x1))/((x2-x0)*(x2-x1));
    return L0*f0+L1*f1+L2*f2;
}
