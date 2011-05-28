/* Properties of Argon
by Ian Bell

Thermo properties from 
---------------------
"A New Equation of State for Argon Covering the Fluid Region
for Temperatures From the Melting Line to 700 K
at Pressures up to 1000 MPa"
Ch. Tegeler, R. Span, and W. Wagner
J. Phys. Chem. Ref. Data, Vol. 28, No. 3, 1999

Transport properties from
------------------------
"Viscosity and Thermal Conductivity Equations for
Nitrogen, Oxygen, Argon, and Air"
E. W. Lemmon and R. T Jacobsen
International Journal of Thermophysics, Vol. 25, No. 1, January 2004

Note: Critical enhancement included


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

//You can include any C libraries that you normally use
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
#include "PropErrorCodes.h"
#include "PropMacros.h"
#include "Argon.h"
#include "time.h"

static int errCode;
static char errStr[ERRSTRLENGTH];

#define nP 200
#define nT 200
const static double Tmin=220,Tmax=800, Pmin=500.03,Pmax=16000;

static double hmat[nT][nP];
static double rhomat[nT][nP];
static double cpmat[nT][nP];
static double smat[nT][nP];
static double cvmat[nT][nP];
static double cmat[nT][nP];
static double umat[nT][nP];
static double viscmat[nT][nP];
static double Tvec[nT];
static double pvec[nP];

static int TablesBuilt;

static const double Tc=150.687, R_Argon=0.2081333, rhoc=535.6, Pc=4863.0, M_Argon=39.948, Ttriple=83.806;
             //           K             kJ/kg-K         kg/m^3     kPa            kg/kmol          K
static const double n[]={0,
0.088722304990011,//[1]
0.70514805167298,//[2]
-1.682011565409,//[3]
-0.14909014431486,//[4]
-0.1202480460094,//[5]
-0.12164978798599,//[6]
0.40035933626752,//[7]
-0.27136062699129,//[8]
0.24211924579645,//[9]
0.005788958318557,//[10]
-0.041097335615341,//[11]
0.024710761541614,//[12]
-0.32181391750702,//[13]
0.33230017695794,//[14]
0.031019986287345,//[15]
-0.030777086002437,//[16]
0.093891137419581,//[17]
-0.090643210682031,//[18]
-0.00045778349276654,//[19]
-0.000082659729025197,//[20]
0.00013013415603147,//[21]
-0.011397840001996,//[22]
-0.024455169960535,//[23]
-0.064324067175955,//[24]
0.058889471093674,//[25]
-0.00064933552112965,//[26]
-0.013889862158435,//[27]
0.4048983929691,//[28]
-0.38612519594749,//[29]
-0.18817142332233,//[30]
0.15977647596482,//[31]
0.053985518513856,//[32]
-0.028953417958014,//[33]
-0.013025413381384,//[34]
0.0028948696775778,//[35]
-0.0022647134304796,//[36]
0.0017616456196368,//[37]
0.0058552454482774,//[38]
-0.69251908270028,//[39]
1.5315490030516,//[40]
-0.0027380447449783//[41]
};

static const int d[]={0,
1,//[1]
1,//[2]
1,//[3]
1,//[4]
1,//[5]
2,//[6]
2,//[7]
2,//[8]
2,//[9]
3,//[10]
3,//[11]
4,//[12]
1,//[13]
1,//[14]
3,//[15]
4,//[16]
4,//[17]
5,//[18]
7,//[19]
10,//[20]
10,//[21]
2,//[22]
2,//[23]
4,//[24]
4,//[25]
8,//[26]
3,//[27]
5,//[28]
5,//[29]
6,//[30]
6,//[31]
7,//[32]
7,//[33]
8,//[34]
9,//[35]
5,//[36]
6,//[37]
2,//[38]
1,//[39]
2,//[40]
3,//[41]
};

static const double t[]={0.00,
0,//[1]
0.25,//[2]
1,//[3]
2.75,//[4]
4,//[5]
0,//[6]
0.25,//[7]
0.75,//[8]
2.75,//[9]
0,//[10]
2,//[11]
0.75,//[12]
3,//[13]
3.5,//[14]
1,//[15]
2,//[16]
4,//[17]
3,//[18]
0,//[19]
0.5,//[20]
1,//[21]
1,//[22]
7,//[23]
5,//[24]
6,//[25]
6,//[26]
10,//[27]
13,//[28]
14,//[29]
11,//[30]
14,//[31]
8,//[32]
14,//[33]
6,//[34]
7,//[35]
24,//[36]
22,//[37]
3,//[38]
1,//[39]
0,//[40]
0//[41]
};

static const int c[]={
0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-12]
1,//[13]
1,//[14]
1,//[15]
1,//[16]
1,//[17]
1,//[18]
1,//[19]
1,//[20]
1,//[21]
2,//[22]
2,//[23]
2,//[24]
2,//[25]
2,//[26]
3,//[27]
3,//[28]
3,//[29]
3,//[30]
3,//[31]
3,//[32]
3,//[33]
3,//[34]
3,//[35]
4,//[36]
4,//[37]
0,0,0,0 // indices [38-41]
};

// alpha is used here for consistency with the definitions in R744.c upon which Argon.c is based
static const double alpha[]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-37]
20,
20,
20,
20
};

static const double beta[]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-37]
250,
375,
300,
225
};

static const double GAMMA[]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-37]
1.11,
1.14,
1.17,
1.11
};

static const double epsilon[]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-37]
1,
1,
1,
1
};

//Constants for ideal gas expression
static const double a0[]={0.0,
	8.31666243,
	-4.94651164
};

static double powI(double x, int y);
static double QuadInterpolate(double x0, double x1, double x2, double f0, double f1, double f2, double x);

static double get_Delta(double T, double P);
static double Pressure_Trho(double T, double rho);
static double IntEnergy_Trho(double T, double rho);
static double Enthalpy_Trho(double T, double rho);
static double Entropy_Trho(double T, double rho);
static double SpecHeatV_Trho(double T, double rho);
static double SpecHeatP_Trho(double T, double rho);
static double SpeedSound_Trho(double T, double rho);

static double dhdT(double tau, double delta);
static double dhdrho(double tau, double delta);
static double dpdT(double tau, double delta);
static double dpdrho(double tau, double delta);

static double LookupValue(char *Prop,double T, double p);
static int isNAN(double x);
static int isINFINITY(double x);

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
	fclose(fp_rho);
	fclose(fp_u);
	fclose(fp_cp);
	fclose(fp_cv);
	fclose(fp_visc);
	
}
static void BuildLookup(void)
{
	int i,j;

	if (!TablesBuilt)
	{
		printf("Building Lookup Tables... Please wait...\n");

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
				if (Tvec[i]>Tc || pvec[j]<psat_Argon(Tvec[i]))
				{					
					rhomat[i][j]=get_Delta(Tvec[i],pvec[j])*rhoc;
					hmat[i][j]=h_Argon(Tvec[i],rhomat[i][j],TYPE_Trho);
					smat[i][j]=s_Argon(Tvec[i],rhomat[i][j],TYPE_Trho);
					umat[i][j]=u_Argon(Tvec[i],rhomat[i][j],TYPE_Trho);
					cpmat[i][j]=cp_Argon(Tvec[i],rhomat[i][j],TYPE_Trho);
					cvmat[i][j]=cv_Argon(Tvec[i],rhomat[i][j],TYPE_Trho);
					cmat[i][j]=c_Argon(Tvec[i],rhomat[i][j],TYPE_Trho);
					viscmat[i][j]=_HUGE;
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
					cmat[i][j]=_HUGE;
					viscmat[i][j]=_HUGE;
				}
			}
		}
		TablesBuilt=1;
		//WriteLookup();
	}
}


/**************************************************/
/*          Public Property Functions            */
/**************************************************/

double rho_Argon(double T, double p, int Types)
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
double p_Argon(double T, double rho)
{
	return Pressure_Trho(T,rho);
}
double h_Argon(double T, double p_rho, int Types)
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
double s_Argon(double T, double p_rho, int Types)
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
double u_Argon(double T, double p_rho, int Types)
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
double cp_Argon(double T, double p_rho, int Types)
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
double cv_Argon(double T, double p_rho, int Types)
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
double w_Argon(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return SpeedSound_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoc;
			return SpeedSound_Trho(T,rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}

double rhosatL_Argon(double T)
{
    const double ti[]={0,0.334,2.0/3.0,7.0/3.0,4.0};
    const double ai[]={0,1.5004262,-0.31381290,0.086461622,-0.041477525};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(summer);
}

double rhosatV_Argon(double T)
{
    const double ti[]={0,0.345,5.0/6.0,1.0,13.0/3.0};
    const double ai[]={0,-1.70695656,-4.02739448,1.55177558,-2.30683228};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(Tc/T*summer);
}

double c_Argon(double T, double p_rho, int Types)
{
	double rho;
	if (Types==TYPE_Trho)
		return SpeedSound_Trho(T,p_rho);
	else
	{
		if (Types==TYPE_TP)
		{
			BuildLookup();
			return LookupValue("c",T,p_rho);
		}
		else //Types==TYPE_TPNoTable
		{
			rho=get_Delta(T,p_rho)*rhoc;
			return SpeedSound_Trho(T,rho);
		}
	}
}

double MM_Argon(void)
{
	return M_Argon;
}

double pcrit_Argon(void)
{
	return Pc;
}

double Tcrit_Argon(void)
{
	return Tc;
}
double Ttriple_Argon(void)
{
	return Ttriple;
}
double rhocrit_Argon(void)
{
	return rhoc;
}
int errCode_Argon(void)
{
	return errCode;
}

double psat_Argon(double T)
{
    const double ti[]={0,1.0,1.5,2.0,4.5};
    const double ai[]={0,-5.9409785,1.3553888,-0.46497607,-1.5399043};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1-T/Tc,ti[i]);
    }
    return Pc*exp(Tc/T*summer);
}

double Tsat_Argon(double P)
{
    double change,eps=.00005;
    int counter=1;
    double r1,r2,r3,T1,T2,T3;
 
    T1=275;
    T2=275+.01;
    r1=psat_Argon(T1)-P;
    r2=psat_Argon(T2)-P;
    
    // End at change less than 0.5%
    while(counter==1 || (fabs(change)/fabs(T2)>eps && counter<40))
    {
        T3=T2-0.5*r2/(r2-r1)*(T2-T1);
        r3=psat_Argon(T3)-P;
        change=0.5*r2/(r2-r1)*(T2-T1);
        T1=T2;
        T2=T3;
        r1=r2;
        r2=r3;
        counter=counter+1;
    }
    return T3;
}   

double rhosat_Argon(double T, double x)
{
    
    if (x>0.5)
        return rhosatV_Argon(T);
    else
        return rhosatL_Argon(T);
}


double visc_Argon(double T, double p_rho, int Types)
{
	double e_k=143.2, //[K]
		   sigma=0.335; //[nm]
	double rho,eta0,etar,OMEGA,delta,tau,Tstar;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,12.19,13.99,0.005027,-18.93,-6.698,-3.827};
	double t[]={0,0.42,0.0,0.95,0.5,0.9,0.8};
	double d[]={0,1,2,10,5,1,2};
	double l[]={0,0,0,0,2,4,4};
	double g[]={0,0,0,0,1,1,1};


	if (Types==TYPE_Trho)
		rho=p_rho;
	else if (Types==TYPE_TPNoLookup)
	{
		rho=get_Delta(T,p_rho)*rhoc;
	}
	else if (Types==TYPE_TP)
	{
		BuildLookup();
		return LookupValue("d",T,p_rho);
	}
	delta=rho/rhoc;
	tau=Tc/T;
	Tstar=T/(e_k);
	OMEGA=exp(b[0]*powI(log(Tstar),0)
			 +b[1]*powI(log(Tstar),1)
		     +b[2]*powI(log(Tstar),2)
			 +b[3]*powI(log(Tstar),3)
		     +b[4]*powI(log(Tstar),4));

	eta0=0.0266958*sqrt(M_Argon*T)/(sigma*sigma*OMEGA);
	etar=N[1]*pow(tau,t[1])*pow(delta,d[1])*exp(-g[1]*pow(delta,l[1]))
		+N[2]*pow(tau,t[2])*pow(delta,d[2])*exp(-g[2]*pow(delta,l[2]))
		+N[3]*pow(tau,t[3])*pow(delta,d[3])*exp(-g[3]*pow(delta,l[3]))
		+N[4]*pow(tau,t[4])*pow(delta,d[4])*exp(-g[4]*pow(delta,l[4]))
		+N[5]*pow(tau,t[5])*pow(delta,d[5])*exp(-g[5]*pow(delta,l[5]))
		+N[6]*pow(tau,t[6])*pow(delta,d[6])*exp(-g[6]*pow(delta,l[6]));

	return (eta0+etar)/1e6; // uPa-s to Pa-s
}

static double X_tilde(double T,double tau,double delta)
{
	// X_tilde is dimensionless
	// Equation 11 slightly rewritten
	double drho_dp;
	drho_dp=1.0/(R_Argon*T*(1+2*delta*dphir_dDelta_Argon(tau,delta)+delta*delta*dphir2_dDelta2_Argon(tau,delta)));
	return Pc*delta/rhoc*drho_dp;
}

double k_Argon(double T, double p_rho, int Types)
{
	double e_k=143.2, //[K]
		   sigma=0.335, //[nm]
		   Tref=301.374, //[K]
		   zeta0=0.13, //[nm]
		   LAMBDA=0.055,
		   q_D=0.32; //[nm]
	double rho,eta0,OMEGA,delta,tau,Tstar,lambda0,lambdar,num,
		cp,cv,OMEGA_tilde,OMEGA_tilde0,zeta,nu,gamma,R0,lambdac,k,
		pi=3.141592654,mu;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,0.8158,-0.4320,0.0,13.73,10.07,0.7375,-33.96,20.47,-2.274,-3.973};
	double t[]={0,0,-0.77,-1.0,0.0,0.0,0.0,0.8,1.2,0.8,0.5};
	double d[]={0,0,0,0,1,2,4,5,6,9,1};
	double l[]={0,0,0,0,0,0,0,2,2,2,4};
	double g[]={0,0,0,0,0,0,0,1,1,1,1};
	
	if (Types==TYPE_Trho)
		rho=p_rho;
	else if (Types==TYPE_TPNoLookup)
	{
		rho=get_Delta(T,p_rho)*rhoc;
	}
	else if (Types==TYPE_TP)
	{
		BuildLookup();
		return LookupValue("D",T,p_rho);
	}
	delta=rho/rhoc;
	tau=Tc/T;
	Tstar=T/(e_k);

	OMEGA=exp(b[0]*powI(log(Tstar),0)
			 +b[1]*powI(log(Tstar),1)
		     +b[2]*powI(log(Tstar),2)
			 +b[3]*powI(log(Tstar),3)
		     +b[4]*powI(log(Tstar),4));

	eta0=0.0266958*sqrt(M_Argon*T)/(sigma*sigma*OMEGA);
	lambda0=N[1]*eta0+N[2]*pow(tau,t[2])+N[3]*pow(tau,t[3]);

	lambdar=N[4]*pow(tau,t[4])*pow(delta,d[4])*exp(-g[4]*pow(delta,l[4]))
		   +N[5]*pow(tau,t[5])*pow(delta,d[5])*exp(-g[5]*pow(delta,l[5]))
		   +N[6]*pow(tau,t[6])*pow(delta,d[6])*exp(-g[6]*pow(delta,l[6]))
		   +N[7]*pow(tau,t[7])*pow(delta,d[7])*exp(-g[7]*pow(delta,l[7]))
	 	   +N[8]*pow(tau,t[8])*pow(delta,d[8])*exp(-g[8]*pow(delta,l[8]))
		   +N[9]*pow(tau,t[9])*pow(delta,d[9])*exp(-g[9]*pow(delta,l[9]))
		   +N[10]*pow(tau,t[10])*pow(delta,d[10])*exp(-g[10]*pow(delta,l[10]));

	R0=1.01;
	nu=0.63;
	gamma=1.2415;
	k=1.380658e-23; //[J/K]

	num=X_tilde(T,Tc/T,delta)-X_tilde(Tref,Tc/Tref,delta)*Tref/T;

	// no critical enhancement if numerator of Eq. 10 is negative
	if (num<0)
		return (lambda0+lambdar)/1e6;

	cp=cp_Argon(T,rho,TYPE_Trho);
	cv=cv_Argon(T,rho,TYPE_Trho);
	mu=visc_Argon(T,rho,TYPE_Trho)*1e6; //[uPa-s]

	zeta=zeta0*pow(num/LAMBDA,nu/gamma); //[nm]
	OMEGA_tilde=2.0/pi*((cp-cv)/cp*atan(zeta/q_D)+cv/cp*(zeta/q_D));
	OMEGA_tilde0=2.0/pi*(1.-exp(-1./(q_D/zeta+1.0/3.0*(zeta/q_D)*(zeta/q_D)/delta/delta)));
	lambdac=rho*(cp*1000.0)*k*R0*T/(6*pi*zeta*mu)*(OMEGA_tilde-OMEGA_tilde0)*1e18; // 1e18 is conversion to mW/m-K (not described in paper)

	return (lambda0+lambdar+lambdac)/1e6;
}

double dhdT_Argon(double T, double p_rho, int Types)
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

double dhdrho_Argon(double T, double p_rho, int Types)
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

double dpdT_Argon(double T, double p_rho, int Types)
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

double dpdrho_Argon(double T, double p_rho, int Types)
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
/*          Private Property Functions            */
/**************************************************/

static double Pressure_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	return R_Argon*T*rho*(1.0+delta*dphir_dDelta_Argon(tau,delta));
}
static double IntEnergy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return R_Argon*T*tau*(dphi0_dTau_Argon(tau,delta)+dphir_dTau_Argon(tau,delta));
}
static double Enthalpy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return R_Argon*T*(1+tau*(dphi0_dTau_Argon(tau,delta)+dphir_dTau_Argon(tau,delta))+delta*dphir_dDelta_Argon(tau,delta));
}
static double Entropy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return R_Argon*(tau*(dphi0_dTau_Argon(tau,delta)+dphir_dTau_Argon(tau,delta))-phi0_Argon(tau,delta)-phir_Argon(tau,delta));
}
static double SpecHeatV_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return -R_Argon*powI(tau,2)*(dphi02_dTau2_Argon(tau,delta)+dphir2_dTau2_Argon(tau,delta));
}
static double SpecHeatP_Trho(double T, double rho)
{
	double delta,tau,c1,c2;
	delta=rho/rhoc;
	tau=Tc/T;

	c1=powI(1.0+delta*dphir_dDelta_Argon(tau,delta)-delta*tau*dphir2_dDelta_dTau_Argon(tau,delta),2);
    c2=(1.0+2.0*delta*dphir_dDelta_Argon(tau,delta)+powI(delta,2)*dphir2_dDelta2_Argon(tau,delta));
    return R_Argon*(-powI(tau,2)*(dphi02_dTau2_Argon(tau,delta)+dphir2_dTau2_Argon(tau,delta))+c1/c2);
}

static double SpeedSound_Trho(double T, double rho)
{
	double delta,tau,c1,c2;
	delta=rho/rhoc;
	tau=Tc/T;

	c1=-SpecHeatV_Trho(T,rho)/R_Argon;
	c2=(1.0+2.0*delta*dphir_dDelta_Argon(tau,delta)+powI(delta,2)*dphir2_dDelta2_Argon(tau,delta));
    return sqrt(-c2*T*SpecHeatP_Trho(T,rho)*1000/c1);
}

/**************************************************/
/*           Property Derivatives                 */
/**************************************************/
// See Lemmon, 2000 for more information
static double dhdrho(double tau, double delta)
{
	double T,R;
	T=Tc/tau;   R=R_Argon;
	//Note: dphi02_dDelta_dTau(tau,delta) is equal to zero
	return R*T/rhoc*(tau*(dphir2_dDelta_dTau_Argon(tau,delta))+dphir_dDelta_Argon(tau,delta)+delta*dphir2_dDelta2_Argon(tau,delta));
}
static double dhdT(double tau, double delta)
{
	double dhdT_rho,T,R,dhdtau;
	T=Tc/tau;   R=R_Argon;
	dhdT_rho=R*tau*(dphi0_dTau_Argon(tau,delta)+dphir_dTau_Argon(tau,delta))+R*delta*dphir_dDelta_Argon(tau,delta)+R;
	dhdtau=R*T*(dphi0_dTau_Argon(tau,delta)+ dphir_dTau_Argon(tau,delta))+R*T*tau*(dphi02_dTau2_Argon(tau,delta)+dphir2_dTau2_Argon(tau,delta))+R*T*delta*dphir2_dDelta_dTau_Argon(tau,delta);
	return dhdT_rho+dhdtau*(-Tc/T/T);
}
static double dpdT(double tau, double delta)
{
	double T,R,rho;
	T=Tc/tau;   R=R_Argon;  rho=delta*rhoc;
	return rho*R*(1+delta*dphir_dDelta_Argon(tau,delta)-delta*tau*dphir2_dDelta_dTau_Argon(tau,delta));
}
static double dpdrho(double tau, double delta)
{
	double T,R,rho;
	T=Tc/tau;   R=R_Argon;  rho=delta*rhoc;
	return R*T*(1+2*delta*dphir_dDelta_Argon(tau,delta)+delta*delta*dphir2_dDelta2_Argon(tau,delta));
}

/**************************************************/
/*          Private Property Functions            */
/**************************************************/

double phir_Argon(double tau, double delta)
{ 
    
    int i;
    double phir=0,psi;
    
    for (i=1;i<=12;i++)
    {
        phir=phir+n[i]*powI(delta,d[i])*pow(tau,t[i]);
    }
    
    for (i=13;i<=37;i++)
    {
        phir=phir+n[i]*powI(delta,d[i])*pow(tau,t[i])*exp(-powI(delta,c[i]));
    }
    
    for (i=38;i<=41;i++)
    {
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        phir=phir+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi;
    }
    return phir;
}

double dphir_dDelta_Argon(double tau, double delta)
{ 
    int i;
    double dphir_dDelta=0,psi;
    double di, ci;
    for (i=1;i<=12;i++)
    {
        di=(double)d[i];
        dphir_dDelta=dphir_dDelta+n[i]*di*powI(delta,d[i]-1)*pow(tau,t[i]);
    }
    for (i=13;i<=37;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir_dDelta=dphir_dDelta+n[i]*exp(-powI(delta,c[i]))*(powI(delta,d[i]-1)*pow(tau,t[i])*(di-ci*powI(delta,c[i])));
    }
    for (i=38;i<=41;i++)
    {
        di=(double)d[i];        
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir_dDelta=dphir_dDelta+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi*(di/delta-2.0*alpha[i]*(delta-epsilon[i]));
    }
    return dphir_dDelta;
}

double dphir2_dDelta2_Argon(double tau, double delta)
{ 
    
    int i;
    double di,ci;
    double dphir2_dDelta2=0,psi;
    for (i=1;i<=12;i++)
    {
        di=(double)d[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*di*(di-1.0)*powI(delta,d[i]-2)*pow(tau,t[i]);
    }
    for (i=13;i<=37;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*exp(-powI(delta,c[i]))*(powI(delta,d[i]-2)*pow(tau,t[i])*( (di-ci*powI(delta,c[i]))*(di-1.0-ci*powI(delta,c[i])) - ci*ci*powI(delta,c[i])));
    }
    for (i=38;i<=41;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir2_dDelta2=dphir2_dDelta2+n[i]*pow(tau,t[i])*psi*(-2.0*alpha[i]*powI(delta,d[i])+4.0*powI(alpha[i],2)*powI(delta,d[i])*powI(delta-epsilon[i],2)-4.0*di*alpha[i]*powI(delta,d[i]-1)*(delta-epsilon[i])+di*(di-1.0)*powI(delta,d[i]-2));
    }
    return dphir2_dDelta2;
}

    
double dphir2_dDelta_dTau_Argon(double tau, double delta)
{ 
    
    int i;
    double di, ci;
    double dphir2_dDelta_dTau=0,psi;

    for (i=1;i<=12;i++)
    {
        di=(double)d[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau + n[i]*di*t[i]*powI(delta,d[i]-1)*pow(tau,t[i]-1.0);
    }
    for (i=13;i<=37;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau + n[i]*exp(-powI(delta,c[i]))*powI(delta,d[i]-1)*t[i]*pow(tau,t[i]-1.0)*(di-ci*powI(delta,c[i]));
    }
    for (i=38;i<=41;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi*(di/delta-2.0*alpha[i]*(delta-epsilon[i]))*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    return dphir2_dDelta_dTau;
}

double dphir_dTau_Argon(double tau, double delta)
{ 
    
    int i;
    double dphir_dTau=0,psi;
    
    for (i=1;i<=12;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powI(delta,d[i])*pow(tau,t[i]-1.0);
    }
    for (i=13;i<=37;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powI(delta,d[i])*pow(tau,t[i]-1.0)*exp(-powI(delta,c[i]));
    }
    for (i=38;i<=41;i++)
    {
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir_dTau=dphir_dTau+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    return dphir_dTau;
}


double dphir2_dTau2_Argon(double tau, double delta)
{ 
    
    int i;
    double dphir2_dTau2=0,psi;
    
    for (i=1;i<=12;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powI(delta,d[i])*pow(tau,t[i]-2.0);
    }
    for (i=13;i<=37;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powI(delta,d[i])*pow(tau,t[i]-2.0)*exp(-powI(delta,c[i]));
    }
    for (i=38;i<=41;i++)
    {
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir2_dTau2=dphir2_dTau2+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi*(powI(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]),2)-t[i]/powI(tau,2)-2.0*beta[i]);
    }
    return dphir2_dTau2;
}

double phi0_Argon(double tau, double delta)
{
    double phi0=0;
    
    phi0=log(delta)+a0[1]+a0[2]*tau+1.5*log(tau);
    return phi0;
}

double dphi0_dDelta_Argon(double tau, double delta)
{
    return 1/delta;
}

double dphi02_dDelta2_Argon(double tau, double delta)
{
    return -1.0/powI(delta,2);
}

double dphi0_dTau_Argon(double tau, double delta)
{
    double dphi0_dTau=0;
    dphi0_dTau=a0[2]+1.5/tau;
    return dphi0_dTau;
}

double dphi02_dTau2_Argon(double tau, double delta)
{
    double dphi02_dTau2=0;
    dphi02_dTau2=-1.5/powI(tau,2);
    return dphi02_dTau2;
}

static double get_Delta(double T, double P)
{
    double change,eps=.00005;
    int counter=1;
    double r1,r2,r3,delta1,delta2,delta3;
    double tau;
    double delta_guess;
 
    if (P>Pc)
    {
        if (T>Tc)
        {
            delta_guess=P/(R_Argon*T)/rhoc;
        }
        else
        {
            delta_guess=1000/rhoc;
        }
    }
    else
    {
        if (T>Tc)
	{
		// Supercritical (try ideal gas)
		delta_guess=P/(R_Argon*T)/rhoc;
	}
	else if(P<psat_Argon(T))
        {
		//Superheated vapor
		delta_guess=(rhosatV_Argon(T)*P/psat_Argon(T))/rhoc;
        }
	else //(T<Tsat_Argon(P))
	{
		delta_guess=(rhosatL_Argon(T))/rhoc;
	}
    }
    tau=Tc/T;
    delta1=delta_guess;
    delta2=delta_guess+.00001;
    r1=P/(delta1*rhoc*R_Argon*T)-1.0-delta1*dphir_dDelta_Argon(tau,delta1);
    r2=P/(delta2*rhoc*R_Argon*T)-1.0-delta2*dphir_dDelta_Argon(tau,delta2);
    
    // End at change less than 0.05%
    while(counter==1 || (fabs(r2)/delta2>eps && counter<40))
    {
        delta3=delta2-r2/(r2-r1)*(delta2-delta1);
        r3=P/(delta3*rhoc*R_Argon*T)-1.0-delta3*dphir_dDelta_Argon(tau,delta3);
        change=r2/(r2-r1)*(delta2-delta1);
        delta1=delta2;
        delta2=delta3;
        r1=r2;
        r2=r3;
        counter=counter+1;
    }
	
    return delta3;
}



static double LookupValue(char *Prop, double T, double p)
{
	int iPlow, iPhigh, iTlow, iThigh,L,R,M,iter;
	double T1, T2, T3, P1, P2, P3, y1, y2, y3, a1, a2, a3;
	double (*mat)[nT][nP];

	//Input checking
	if (T>Tmax || T<Tmin)
	{
		errCode=OUT_RANGE_T;
	}
	if (p>Pmax || p<Pmin)
	{
		errCode=OUT_RANGE_P;
	}
	if (T<Tmin || T > Tmax || p<Pmin || p>Pmax || isNAN(T) || isINFINITY(T) ||isNAN(p) || isINFINITY(p))
    {
		printf("Inputs to LookupValue(%s) of Argon out of bounds:  T=%g K \t P=%g kPa \n",Prop,T,p);
    }

	L=0;
	R=nT-1;
	M=(L+R)/2;
	iter=0;
	// Use interval halving to find the indices which bracket the temperature of interest
	while (R-L>1)
	{
		if (T>=Tvec[M])
		{ L=M; M=(L+R)/2; continue;}
		if (T<Tvec[M])
		{ R=M; M=(L+R)/2; continue;}
		iter++;
		if (iter>100)
			printf("Problem with T(%g K)\n",T);
	}
	iTlow=L; iThigh=R;

	L=0;
	R=nP-1;
	M=(L+R)/2;
	iter=0;
	// Use interval halving to find the indices which bracket the pressure of interest
	while (R-L>1)
	{
		if (p>=pvec[M])
		{ L=M; M=(L+R)/2; continue;}
		if (p<pvec[M])
		{ R=M; M=(L+R)/2; continue;}
		iter++;
		if (iter>100)
			printf("Problem with p(%g kPa)\n",p);
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
	a1=QuadInterpolate(P1,P2,P3,y1,y2,y3,p);

	//At High Temperature Index
	y1=(*mat)[iThigh][iPlow];
	y2=(*mat)[iThigh][iPhigh];
	y3=(*mat)[iThigh][iPhigh+1];
	a2=QuadInterpolate(P1,P2,P3,y1,y2,y3,p);

	//At High Temperature Index+1 (for QuadInterpolate() )
	y1=(*mat)[iThigh+1][iPlow];
	y2=(*mat)[iThigh+1][iPhigh];
	y3=(*mat)[iThigh+1][iPhigh+1];
	a3=QuadInterpolate(P1,P2,P3,y1,y2,y3,p);

	//At Final Interpolation
	T1=Tvec[iTlow];
	T2=Tvec[iThigh];
	T3=Tvec[iThigh+1];
	return QuadInterpolate(T1,T2,T3,a1,a2,a3,T);
	
}

static double powI(double x, int y)
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

static double QuadInterpolate(double x0, double x1, double x2, double f0, double f1, double f2, double x)
{
    double L0, L1, L2;
    L0=((x-x1)*(x-x2))/((x0-x1)*(x0-x2));
    L1=((x-x0)*(x-x2))/((x1-x0)*(x1-x2));
    L2=((x-x0)*(x-x1))/((x2-x0)*(x2-x1));
    return L0*f0+L1*f1+L2*f2;
}

static int isNAN(double x)
{
	// recommendation from http://www.devx.com/tips/Tip/42853
	return x != x;
}

static int isINFINITY(double x)
{
	// recommendation from http://www.devx.com/tips/Tip/42853
	if ((x == x) && ((x - x) != 0.0)) 
		return 1;//return (x < 0.0 ? -1 : 1); // This will tell you whether positive or negative infinity
	else 
		return 0;
}





