/* Properties of Carbon Dioxide (R744)
by Ian Bell

Themo properties from 
"A New Equation of State for Carbon Dioxide Covering the Fluid Region from the 
Triple Point Temperature to 1100 K at Pressures up to 800 MPa", 
R. Span and W. Wagner, J. Phys. Chem. Ref. Data, v. 25, 1996

WARNING: TRANSPORT PROPERTIES NOT CODED!!

In order to call the exposed functions, rho_, h_, s_, cp_,...... there are three 
different ways the inputs can be passed, and this is expressed by the Types integer flag.  
These macros are defined in the PropMacros.h header file:
1) First parameter temperature, second parameter pressure ex: h_R410A(260,1785,1)=-67.53
	In this case, the lookup tables are built if needed and then interpolated
2) First parameter temperature, second parameter density ex: h_R410A(260,43.29,2)=-67.53
	Density and temp plugged directly into EOS
3) First parameter temperature, second parameter pressure ex: h_R410A(260,1785,3)=-67.53
	Density solved for, then plugged into EOS (can be quite slow)


Sanity Checks based on data from EES

P[i]	T[i]	Cp[i]		Cv[i]		D[i]		H[i]	S[i]		V[i]		k[i]	mu[i]
[kPa]	[K]		[kJ/kg-K]	[kJ/kg-K]	[kg/m^3]	[kJ/kg]	[kJ/kg-K]	[m^3/kg]	[W/m-K]	[Pa-s]
1785	240		2.043		0.9458		1091		-379.9	-2.017		0.000917	0.1508	0.0001788
1785	250		2.131		0.9364		1046		-359.1	-1.932		0.000956	0.1377	0.0001523
1785	250		1.219		0.7458		46.64		-69.74	-0.7749		0.02144		0.01421	0.00001288
1785	260		1.188		0.7334		43.29		-67.53	-0.7571		0.0231		0.01483	0.00001335

P[i]	T[i]	Cp[i]		Cv[i]		D[i]		H[i]	S[i]		V[i]		k[i]	mu[i]
[kPa]	[K]		[kJ/kg-K]	[kJ/kg-K]	[kg/m^3]	[kJ/kg]	[kJ/kg-K]	[m^3/kg]	[W/m-K]	[Pa-s]
3659	265		2.298		0.9325		978.9		-326.5	-1.813		0.001022	0.1205	0.000122
3659	275		2.599		0.9484		916.2		-302.2	-1.723		0.001092	0.1063	0.0001011
3659	275		1.868		0.8847		103.6		-77.02	-0.9042		0.009655	0.01944	0.00001494
3659	285		1.701		0.8447		92.33		-70.85	-0.8726		0.01083		0.01896	0.00001523

P[i]	T[i]	Cp[i]		Cv[i]		D[i]		H[i]	S[i]		V[i]		k[i]	mu[i]
[kPa]	[K]		[kJ/kg-K]	[kJ/kg-K]	[kg/m^3]	[kJ/kg]	[kJ/kg-K]	[m^3/kg]	[W/m-K]	[Pa-s]
6713	290		3.044		0.9574		834.4		-265.9	-1.607		0.001199	0.09088	0.00007991
6713	300		7.332		1.111		681.9		-223.9	-1.465		0.001466	0.07834	0.00005326
6713	300		6.703		1.246		268.2		-119.6	-1.117		0.003728	0.04769	0.00002118
6713	310		3.33		0.9873		191.9		-88.7	-1.005		0.005211	0.02957	0.00001885

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
#include "R744.h"
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

static double alpha[40],beta[43],GAMMA[40],epsilon[40],a[43],b[43],A[43],B[43],C[43],D[43],a0[9],theta0[9];
static const double Tc=304.128, R_R744=0.1889241, rhoc=467.6, Pc=7377.3;
             //    K               kJ/kg-K       kg/m^3          kPa
static const double n[]={0,    
 0.3885682320316100E+00,
 0.2938547594274000E+01,
-0.5586718853493400E+01,
-0.7675319959247700E+00,
 0.3172900558041600E+00,
 0.5480331589776700E+00,
 0.1227941122033500E+00,
 
 0.2165896154322000E+01,
 0.1584173510972400E+01,
-0.2313270540550300E+00,
 0.5811691643143600E-01,
-0.5536913720538200E-00,
 0.4894661590942200E-00,
-0.2427573984350100E-01,
 0.6249479050167800E-01,
-0.1217586022524600E+00,
-0.3705568527008600E+00,
-0.1677587970042600E-01,
-0.1196073663798700E+00,
-0.4561936250877800E-01,
 0.3561278927034600E-01, 
-0.7442772713205200E-02,
-0.1739570490243200E-02,
-0.2181012128952700E-01,
 0.2433216655923600E-01,
-0.3744013342346300E-01,
 0.1433871575687800E-00,
-0.1349196908328600E-00,
-0.2315122505348000E-01,
 0.1236312549290100E-01,
 0.2105832197294000E-02,
-0.3395851902636800E-03,
 0.5599365177159200E-02,
-0.3033511805564600E-03,

-0.2136548868832000E+03,
 0.2664156914927200E+05,
-0.2402721220455700E+05,
-0.2834160342399900E+03,
 0.2124728440017900E+03,
 
-0.6664227654075100E+00,
 0.7260863234989700E+00,
 0.5506866861284200E-01};

static const int d[]={0,
1,
1,
1,
1,
2,
2,
3,
1,
2,
4,
5,
5,
5,
6,
6,
6,
1,
1,
4,
4,
4,
7,
8,
2,
3,
3,
5,
5,
6,
7,
8,
10,
4,
8,
2,
2,
2,
3,
3};

static const double t[]={0.00,
0.00,
0.75,
1.00,
2.00,
0.75,
2.00,
0.75,
1.50,
1.50,
2.50,
0.00,
1.50,
2.00,
0.00,
1.00,
2.00,
3.00,
6.00,
3.00,
6.00,
8.00,
6.00,
0.00,
7.00,
12.00,
16.00,
22.00,
24.00,
16.00,
24.00,
8.00,
2.00,
28.00,
14.00,
1.00,
0.00,
1.00,
3.00,
3.00};

static const int c[]={0,0,0,0,0,0,0,0,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
2,
2,
2,
2,
2,
2,
3,
3,
3,
4,
4,
4,
4,
4,
4,
5,
6};

void setCoeffs(void)
{
alpha[35]=25.0;
alpha[36]=25.0;
alpha[37]=25.0;
alpha[38]=15.0;
alpha[39]=20.0;

beta[35]=325.0;
beta[36]=300.0;
beta[37]=300.0;
beta[38]=275.0;
beta[39]=275.0;

GAMMA[35]=1.16;
GAMMA[36]=1.19;
GAMMA[37]=1.19;
GAMMA[38]=1.25;
GAMMA[39]=1.22;

epsilon[35]=1.00;
epsilon[36]=1.00;
epsilon[37]=1.00;
epsilon[38]=1.00;
epsilon[39]=1.00;

a[40]=3.5;
a[41]=3.5;
a[42]=3.0;

b[40]=0.875;
b[41]=0.925;
b[42]=0.875;

beta[40]=0.300;
beta[41]=0.300;
beta[42]=0.300;

A[40]=0.700;
A[41]=0.700;
A[42]=0.700;

B[40]=0.3;
B[41]=0.3;
B[42]=1.0;

C[40]=10.0;
C[41]=10.0;
C[42]=12.5;

D[40]=275.0;
D[41]=275.0;
D[42]=275.0;

//Constants for ideal gas expression
a0[1]=8.37304456;
a0[2]=-3.70454304;
a0[3]=2.500000;
a0[4]=1.99427042;
a0[5]=0.62105248;
a0[6]=0.41195293;
a0[7]=1.04028922;
a0[8]=0.08327678;

theta0[4]=3.15163;
theta0[5]=6.11190;
theta0[6]=6.77708;
theta0[7]=11.32384;
theta0[8]=27.08792;
}

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

static double phir(double tau, double delta);
static double phi0(double tau, double delta);
static double dphir_dDelta(double tau, double delta);
static double dphir2_dDelta2(double tau, double delta);
static double dphir_dTau(double tau, double delta);
static double dphi0_dDelta(double tau, double delta);
static double dphi02_dDelta2(double tau, double delta);
static double dphi0_dTau(double tau, double delta);
static double dphi02_dTau2(double tau, double delta);
static double dphir2_dTau2(double tau, double delta);
static double dphir2_dDelta_dTau(double tau, double delta);

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
				if (Tvec[i]>Tc || pvec[j]<Psat_R744(Tvec[i]))
				{					
					rhomat[i][j]=get_Delta(Tvec[i],pvec[j])*rhoc;
					hmat[i][j]=h_R744(Tvec[i],rhomat[i][j],TYPE_Trho);
					smat[i][j]=s_R744(Tvec[i],rhomat[i][j],TYPE_Trho);
					umat[i][j]=u_R744(Tvec[i],rhomat[i][j],TYPE_Trho);
					cpmat[i][j]=cp_R744(Tvec[i],rhomat[i][j],TYPE_Trho);
					cvmat[i][j]=cv_R744(Tvec[i],rhomat[i][j],TYPE_Trho);
					cmat[i][j]=c_R744(Tvec[i],rhomat[i][j],TYPE_Trho);
					viscmat[i][j]=visc_R744(Tvec[i],rhomat[i][j],TYPE_Trho);
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

double rho_R744(double T, double p, int Types)
{
	setCoeffs();
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
double p_R744(double T, double rho)
{
	setCoeffs();
	return Pressure_Trho(T,rho);
}
double h_R744(double T, double p_rho, int Types)
{
	double rho;
	setCoeffs();
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
double s_R744(double T, double p_rho, int Types)
{
	double rho;
	setCoeffs();
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
double u_R744(double T, double p_rho, int Types)
{
	double rho;
	setCoeffs();
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
double cp_R744(double T, double p_rho, int Types)
{
	double rho;
	setCoeffs();
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
double cv_R744(double T, double p_rho, int Types)
{
	double rho;
	setCoeffs();
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

double rhosatV_R744(double T)
{
    const double ti[]={0,0.340,1.0/2.0,1.0,7.0/3.0,14.0/3.0};
    const double ai[]={0,-1.7074879,-0.82274670,-4.6008549,-10.111178,-29.742252};
    double summer=0;
    int i;
    for (i=1;i<=5;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(summer);
    
}

double rhosatL_R744(double T)
{
    const double ti[]={0,0.340,1.0/2.0,10.0/6.0,11.0/6.0};
    const double ai[]={0,1.9245108,-0.62385555,-0.32731127,0.39245142};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(summer);
}

double w_R744(double T, double p_rho, int Types)
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

double c_R744(double T, double p_rho, int Types)
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

double visc_R744(double T,double p_rho, int Types)
{
	int i;
	double e_k=251.196,Tstar,sumGstar=0.0,Gstar,eta0,delta_eta,rho;
	double a[]={0.235156,-0.491266,5.211155e-2,5.347906e-2,-1.537102e-2};
	double d11=0.4071119e-2,d21=0.7198037e-4,d64=0.2411697e-16,d81=0.2971072e-22,d82=-0.1627888e-22;

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

	Tstar=T/e_k;
	for (i=0;i<=4;i++)
	{
		sumGstar=sumGstar+a[i]*powI(log(Tstar),i);
	}
	Gstar=exp(sumGstar);
	eta0=1.00697*sqrt(T)/Gstar;
	delta_eta=d11*rho+d21*rho*rho*d64*powI(rho,6)/powI(Tstar,3)+d81*powI(rho,8)+d82*powI(rho,8)/Tstar;

	return (eta0+delta_eta)/1e6;


}

double k_R744(double T,double p_rho, int Types)
{
	
	fprintf(stderr,"Thermal conductivity not coded for R744 (CO2).  Sorry.\n");
	return _HUGE;
}

double MM_R744(void)
{
	return 44.01;
}

double pcrit_R744(void)
{
	return Pc;
}

double Tcrit_R744(void)
{
	return Tc;
}

int errCode_R744(void)
{
	return errCode;
}


double dhdT_R744(double T, double p_rho, int Types)
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

double dhdrho_R744(double T, double p_rho, int Types)
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

double dpdT_R744(double T, double p_rho, int Types)
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

double dpdrho_R744(double T, double p_rho, int Types)
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
	return R_R744*T*rho*(1.0+delta*dphir_dDelta(tau,delta));
}
static double IntEnergy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return R_R744*T*tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta));
}
static double Enthalpy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return R_R744*T*(1+tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+delta*dphir_dDelta(tau,delta));
}
static double Entropy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return R_R744*(tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))-phi0(tau,delta)-phir(tau,delta));
}
static double SpecHeatV_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return -R_R744*powI(tau,2)*(dphi02_dTau2(tau,delta)+dphir2_dTau2(tau,delta));
}
static double SpecHeatP_Trho(double T, double rho)
{
	double delta,tau,c1,c2;
	delta=rho/rhoc;
	tau=Tc/T;

	c1=powI(1.0+delta*dphir_dDelta(tau,delta)-delta*tau*dphir2_dDelta_dTau(tau,delta),2);
    c2=(1.0+2.0*delta*dphir_dDelta(tau,delta)+powI(delta,2)*dphir2_dDelta2(tau,delta));
    return R_R744*(-powI(tau,2)*(dphi02_dTau2(tau,delta)+dphir2_dTau2(tau,delta))+c1/c2);
}

static double SpeedSound_Trho(double T, double rho)
{
	double delta,tau,c1,c2;
	delta=rho/rhoc;
	tau=Tc/T;

	c1=-SpecHeatV_Trho(T,rho)/R_R744;
	c2=(1.0+2.0*delta*dphir_dDelta(tau,delta)+powI(delta,2)*dphir2_dDelta2(tau,delta));
    return sqrt(-c2*T*SpecHeatP_Trho(T,rho)*1000/c1);
}

/**************************************************/
/*           Property Derivatives                 */
/**************************************************/
// See Lemmon, 2000 for more information
static double dhdrho(double tau, double delta)
{
	double T,R;
	T=Tc/tau;   R=R_R744;
	//Note: dphi02_dDelta_dTau(tau,delta) is equal to zero
	return R*T/rhoc*(tau*(dphir2_dDelta_dTau(tau,delta))+dphir_dDelta(tau,delta)+delta*dphir2_dDelta2(tau,delta));
}
static double dhdT(double tau, double delta)
{
	double dhdT_rho,T,R,dhdtau;
	T=Tc/tau;   R=R_R744;
	dhdT_rho=R*tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+R*delta*dphir_dDelta(tau,delta)+R;
	dhdtau=R*T*(dphi0_dTau(tau,delta)+ dphir_dTau(tau,delta))+R*T*tau*(dphi02_dTau2(tau,delta)+dphir2_dTau2(tau,delta))+R*T*delta*dphir2_dDelta_dTau(tau,delta);
	return dhdT_rho+dhdtau*(-Tc/T/T);
}
static double dpdT(double tau, double delta)
{
	double T,R,rho;
	T=Tc/tau;   R=R_R744;  rho=delta*rhoc;
	return rho*R*(1+delta*dphir_dDelta(tau,delta)-delta*tau*dphir2_dDelta_dTau(tau,delta));
}
static double dpdrho(double tau, double delta)
{
	double T,R,rho;
	T=Tc/tau;   R=R_R744;  rho=delta*rhoc;
	return R*T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*dphir2_dDelta2(tau,delta));

}

/**************************************************/
/*          Private Property Functions            */
/**************************************************/

static double phir(double tau, double delta)
{ 
    
    int i;
    double phir=0,theta,DELTA,PSI,psi;
    
    for (i=1;i<=7;i++)
    {
        phir=phir+n[i]*powI(delta,d[i])*pow(tau,t[i]);
    }
    
    for (i=8;i<=34;i++)
    {
        phir=phir+n[i]*powI(delta,d[i])*pow(tau,t[i])*exp(-powI(delta,c[i]));
    }
    
    for (i=35;i<=39;i++)
    {
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        phir=phir+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi;
    }
    
    for (i=40;i<=42;i++)
    {
        theta=(1.0-tau)+A[i]*pow(powI(delta-1.0,2),1/(2*beta[i]));
        DELTA=powI(theta,2)+B[i]*pow(powI(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powI(delta-1.0,2)-D[i]*powI(tau-1.0,2));
        phir=phir+n[i]*pow(DELTA,b[i])*delta*PSI;
    }
    
    return phir;
}



static double dphir_dDelta(double tau, double delta)
{ 
    
    int i;
    double dphir_dDelta=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta,psi;
    double di, ci;
//     mexPrintf("dphir_dDelta\t");
    for (i=1;i<=7;i++)
    {
        di=(double)d[i];

        dphir_dDelta=dphir_dDelta+n[i]*di*powI(delta,d[i]-1)*pow(tau,t[i]);
    }
//    mexPrintf("%10.10g\t",dphir_dDelta);
    for (i=8;i<=34;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir_dDelta=dphir_dDelta+n[i]*exp(-powI(delta,c[i]))*(powI(delta,d[i]-1)*pow(tau,t[i])*(di-ci*powI(delta,c[i])));
    }
//     mexPrintf("%10.10g\t",dphir_dDelta);
    for (i=35;i<=39;i++)
    {
        di=(double)d[i];        
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir_dDelta=dphir_dDelta+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi*(di/delta-2.0*alpha[i]*(delta-epsilon[i]));
    }
// mexPrintf("%10.10g\t",dphir_dDelta);
    for (i=40;i<=42;i++)
    {
        theta=(1.0-tau)+A[i]*pow(powI(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=powI(theta,2)+B[i]*pow(powI(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powI(delta-1.0,2)-D[i]*powI(tau-1.0,2));
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(powI(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(powI(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        dphir_dDelta=dphir_dDelta+n[i]*(pow(DELTA,b[i])*(PSI+delta*dPSI_dDelta)+dDELTAbi_dDelta*delta*PSI);
    }
//     mexPrintf("%10.10g\n",dphir_dDelta);
    return dphir_dDelta;
}

static double dphir2_dDelta2(double tau, double delta)
{ 
    
    int i;
    double di,ci;
    
    double dphir2_dDelta2=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta,psi,dPSI2_dDelta2,dDELTAbi2_dDelta2,dDELTA2_dDelta2;
//     mexPrintf("dphir2_dDelta2\t");
    for (i=1;i<=7;i++)
    {
        di=(double)d[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*di*(di-1.0)*powI(delta,d[i]-2)*pow(tau,t[i]);
    }
//     mexPrintf("%10.10g\t",dphir2_dDelta2);
    for (i=8;i<=34;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*exp(-powI(delta,c[i]))*(powI(delta,d[i]-2)*pow(tau,t[i])*( (di-ci*powI(delta,c[i]))*(di-1.0-ci*powI(delta,c[i])) - ci*ci*powI(delta,c[i])));
    }
//     mexPrintf("%10.10g\t",dphir2_dDelta2);
    for (i=35;i<=39;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir2_dDelta2=dphir2_dDelta2+n[i]*pow(tau,t[i])*psi*(-2.0*alpha[i]*powI(delta,d[i])+4.0*powI(alpha[i],2)*powI(delta,d[i])*powI(delta-epsilon[i],2)-4.0*di*alpha[i]*powI(delta,d[i]-1)*(delta-epsilon[i])+di*(di-1.0)*powI(delta,d[i]-2));
    }
//     mexPrintf("%10.10g\t",dphir2_dDelta2);
    for (i=40;i<=42;i++)
    {
               
        theta=(1.0-tau)+A[i]*pow(powI(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=powI(theta,2)+B[i]*pow(powI(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powI(delta-1.0,2)-D[i]*powI(tau-1.0,2));
        
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(powI(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(powI(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
        dPSI2_dDelta2=(2.0*C[i]*powI(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
        dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+powI(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(powI(delta-1.0,2),a[i]-2.0)+2.0*powI(A[i]/beta[i],2)*powI(pow(powI(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(powI(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
        dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*powI(dDELTA_dDelta,2));
        
        dphir2_dDelta2=dphir2_dDelta2+n[i]*(pow(DELTA,b[i])*(2.0*dPSI_dDelta+delta*dPSI2_dDelta2)+2.0*dDELTAbi_dDelta*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta2*delta*PSI);
    }
//     mexPrintf("%10.10g\n",dphir2_dDelta2);
    return dphir2_dDelta2;
}

    
static double dphir2_dDelta_dTau(double tau, double delta)
{ 
    
    int i;
    double di, ci;
    double dphir2_dDelta_dTau=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta,psi,dPSI2_dDelta2,dDELTAbi2_dDelta2,dDELTA2_dDelta2;
    double dPSI2_dDelta_dTau, dDELTAbi2_dDelta_dTau, dPSI_dTau, dDELTAbi_dTau, dPSI2_dTau2, dDELTAbi2_dTau2;
//     mexPrintf("dphir2_dDelta_dTau\t");
    for (i=1;i<=7;i++)
    {
        di=(double)d[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*di*t[i]*powI(delta,d[i]-1)*pow(tau,t[i]-1.0);
    }
//     mexPrintf("%10.10g\t",dphir2_dDelta_dTau);
    for (i=8;i<=34;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*exp(-powI(delta,c[i]))*powI(delta,d[i]-1)*t[i]*pow(tau,t[i]-1.0)*(di-ci*powI(delta,c[i]));
    }
//     mexPrintf("%10.10g\t",dphir2_dDelta_dTau);
    for (i=35;i<=39;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi*(di/delta-2.0*alpha[i]*(delta-epsilon[i]))*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
//     mexPrintf("%10.10g\t",dphir2_dDelta_dTau);
    for (i=40;i<=42;i++)
    {
        
        theta=(1.0-tau)+A[i]*pow(powI(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=powI(theta,2)+B[i]*pow(powI(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powI(delta-1.0,2)-D[i]*powI(tau-1.0,2));
        
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(powI(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(powI(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
        dPSI2_dDelta2=(2.0*C[i]*powI(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
        dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+powI(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(powI(delta-1.0,2),a[i]-2.0)+2.0*powI(A[i]/beta[i],2)*powI(pow(powI(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(powI(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
        dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*powI(dDELTA_dDelta,2));
        
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        dPSI2_dTau2=(2.0*D[i]*powI(tau-1.0,2)-1.0)*2.0*D[i]*PSI;
        dDELTAbi2_dTau2=2.0*b[i]*pow(DELTA,b[i]-1.0)+4.0*powI(theta,2)*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0);
        
        dPSI2_dDelta_dTau=4.0*C[i]*D[i]*(delta-1.0)*(tau-1.0)*PSI;
        dDELTAbi2_dDelta_dTau=-A[i]*b[i]*2.0/beta[i]*pow(DELTA,b[i]-1.0)*(delta-1.0)*pow(powI(delta-1.0,2),1.0/(2.0*beta[i])-1.0)-2.0*theta*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0)*dDELTA_dDelta;
        
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*(pow(DELTA,b[i])*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+delta*dDELTAbi_dDelta*dPSI_dTau+ dDELTAbi_dTau*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta_dTau*delta*PSI);
    }
//     mexPrintf("%10.10g\n",dphir2_dDelta_dTau);
    return dphir2_dDelta_dTau;
}

static double dphir_dTau(double tau, double delta)
{ 
    
    int i;
    double dphir_dTau=0,theta,DELTA,PSI,dPSI_dTau,dDELTAbi_dTau,psi;
    
    for (i=1;i<=7;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powI(delta,d[i])*pow(tau,t[i]-1.0);
    }
    
    for (i=8;i<=34;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powI(delta,d[i])*pow(tau,t[i]-1.0)*exp(-powI(delta,c[i]));
    }
    
    for (i=35;i<=39;i++)
    {
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir_dTau=dphir_dTau+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    
    for (i=40;i<=42;i++)
    {
        theta=(1.0-tau)+A[i]*pow(powI(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=powI(theta,2)+B[i]*pow(powI(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powI(delta-1.0,2)-D[i]*powI(tau-1.0,2));
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        dphir_dTau=dphir_dTau+n[i]*delta*(dDELTAbi_dTau*PSI+pow(DELTA,b[i])*dPSI_dTau);
    }
    
    return dphir_dTau;
}


static double dphir2_dTau2(double tau, double delta)
{ 
    
    int i;
    double dphir2_dTau2=0,theta,DELTA,PSI,dPSI_dTau,dDELTAbi_dTau,psi,dPSI2_dTau2,dDELTAbi2_dTau2;
    
    for (i=1;i<=7;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powI(delta,d[i])*pow(tau,t[i]-2.0);
    }
    
    for (i=8;i<=34;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powI(delta,d[i])*pow(tau,t[i]-2.0)*exp(-powI(delta,c[i]));
    }
    
    for (i=35;i<=39;i++)
    {
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir2_dTau2=dphir2_dTau2+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi*(powI(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]),2)-t[i]/powI(tau,2)-2.0*beta[i]);
    }
    
    for (i=40;i<=42;i++)
    {
        theta=(1.0-tau)+A[i]*pow(powI(delta-1.0,2),1/(2*beta[i]));
        DELTA=powI(theta,2)+B[i]*pow(powI(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powI(delta-1.0,2)-D[i]*powI(tau-1.0,2));
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        dPSI2_dTau2=(2.0*D[i]*powI(tau-1.0,2)-1.0)*2.0*D[i]*PSI;
        dDELTAbi2_dTau2=2.0*b[i]*pow(DELTA,b[i]-1.0)+4.0*powI(theta,2)*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0);
        dphir2_dTau2=dphir2_dTau2+n[i]*delta*(dDELTAbi2_dTau2*PSI+2.0*dDELTAbi_dTau*dPSI_dTau+pow(DELTA,b[i])*dPSI2_dTau2);
    }
    
    return dphir2_dTau2;
}

static double dphi0_dDelta(double tau, double delta)
{
    return 1/delta;
}

static double dphi02_dDelta2(double tau, double delta)
{
    return -1.0/powI(delta,2);
}

static double phi0(double tau, double delta)
{
    double phi0=0;
    int i;
    
    phi0=log(delta)+a0[1]+a0[2]*tau+a0[3]*log(tau);
    for (i=4;i<=8;i++)
    {
        phi0=phi0+a0[i]*log(1.0-exp(-theta0[i]*tau));
    }
    return phi0;
}


static double dphi0_dTau(double tau, double delta)
{
    double dphi0_dTau=0;
    int i;
    dphi0_dTau=a0[2]+a0[3]/tau;

    for (i=4;i<=8;i++)
    {
        dphi0_dTau=dphi0_dTau+a0[i]*theta0[i]*(1.0/(1.0-exp(-theta0[i]*tau))-1.0);
    }
    return dphi0_dTau;
}

static double dphi02_dTau2(double tau, double delta)
{
    double dphi02_dTau2=0;
    int i;
    
    dphi02_dTau2=-a0[3]/powI(tau,2);
    for (i=4;i<=8;i++)
    {
        dphi02_dTau2=dphi02_dTau2-a0[i]*powI(theta0[i],2)*exp(-theta0[i]*tau)/powI(1.0-exp(-theta0[i]*tau),2);
    }
    return dphi02_dTau2;
}

static double get_Delta(double T, double P)
{
    
    double change,eps=.0005;
    int counter=1;
    double r1,r2,r3,delta1,delta2,delta3;
    double tau;
    double delta_guess;
     setCoeffs();
 
    if (P>Pc)
    {
        if (T>Tc)
        {
            delta_guess=P/(R_R744*T)/rhoc;
        }
        else
        {
            delta_guess=1000/rhoc;
        }
//         return R744_rho_fit_supercrit(T,P)/rhoc;
    }
    else
    {
//         return R744_rho_fit_subcrit(T,P)/rhoc;
        if (T>Tsat_R744(P))
        {
            delta_guess=P/(R_R744*T)/rhoc;
        }
        else
        {
            delta_guess=1000/rhoc;
        }
    }
    
   
    tau=Tc/T;
    delta1=delta_guess;
    delta2=delta_guess+.00001;
    r1=P/(delta1*rhoc*R_R744*T)-1.0-delta1*dphir_dDelta(tau,delta1);
    r2=P/(delta2*rhoc*R_R744*T)-1.0-delta2*dphir_dDelta(tau,delta2);
    
    // End at change less than 0.05%
    while(counter==1 || (fabs(change)/fabs(delta2)>eps && counter<40))
    {
        delta3=delta2-r2/(r2-r1)*(delta2-delta1);
        r3=P/(delta3*rhoc*R_R744*T)-1.0-delta3*dphir_dDelta(tau,delta3);
        change=r2/(r2-r1)*(delta2-delta1);
        delta1=delta2;
        delta2=delta3;
        r1=r2;
        r2=r3;
        counter=counter+1;

//         mexPrintf("Iteration: %i \n",counter);    
//         mexPrintf("%g \t %g \t %g \t %g \t %g\n",delta2,r2,change,T,P);    
        
    }
   //printf("Iteration: %d \n",counter);
    return delta3;
}

double Psat_R744(double T)
{
    const double ti[]={0,1.0,1.5,2.0,4.0};
    const double ai[]={0,-7.0602087,1.9391218,-1.6463597,-3.2995634};
    double summer=0;
    int i;
    setCoeffs();
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1-T/Tc,ti[i]);
    }
    return Pc*exp(Tc/T*summer);
}

double Tsat_R744(double P)
{
    double change,eps=.00005;
    int counter=1;
    double r1,r2,r3,T1,T2,T3;
 
    T1=275;
    T2=275+.01;
    r1=Psat_R744(T1)-P;
    r2=Psat_R744(T2)-P;
    
    // End at change less than 0.5%
    while(counter==1 || (fabs(change)/fabs(T2)>eps && counter<40))
    {
        T3=T2-0.5*r2/(r2-r1)*(T2-T1);
        r3=Psat_R744(T3)-P;
        change=0.5*r2/(r2-r1)*(T2-T1);
        T1=T2;
        T2=T3;
        r1=r2;
        r2=r3;
        counter=counter+1;
    }
    return T3;
}   

double hsat_R744(double T, double x)
{
    double delta,tau;
    
    if (x>0.5)
    {
        delta=rhosatV_R744(T)/rhoc;
        tau=Tc/T;
        return R_R744*T*(1+tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+delta*dphir_dDelta(tau,delta));
    }
    else
    {
        delta=rhosatL_R744(T)/rhoc;
        tau=Tc/T;
        return R_R744*T*(1+tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+delta*dphir_dDelta(tau,delta));
    }   
}

double ssat_R744(double T, double x)
{
    double delta,tau;
    
    if (x>0.5)
    {
        delta=rhosatV_R744(T)/rhoc;
        tau=Tc/T;
        return R_R744*(tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))-phi0(tau,delta)-phir(tau,delta));
    }
    else
    {
        delta=rhosatL_R744(T)/rhoc;
        tau=Tc/T;
        return R_R744*(tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))-phi0(tau,delta)-phir(tau,delta));
    }   
}

double rhosat_R744(double T, double x)
{
    
    if (x>0.5)
    {
        return rhosatV_R744(T);
    }
    else
    {
        return rhosatL_R744(T);
    }   
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
		printf("Inputs to LookupValue(%s) of R744 out of bounds:  T=%g K \t P=%g kPa \n",Prop,T,p);
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





