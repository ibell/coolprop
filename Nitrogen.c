/* Properties of Nitrogen
by Ian Bell

Thermo properties from 
---------------------
"A Reference Equation of State for the Thermodynamic Properties
of Nitrogen for Temperatures from 63.151 to 1000 K 
and Pressures to 2200 MPa", 
R. Span and E.W. Lemmon and R.T. Jacobsen and W. Wagner and A. Yokozeki, 
J. Phys. Chem. Ref. Data, v. 29, n. 6, 2000

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
1) First parameter temperature, second parameter pressure ex: h_Nitrogen(260,1785,1)
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
#include "PropErrorCodes.h"
#include "PropMacros.h"
#include "Nitrogen.h"
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

static const double Tc=126.192, R_Nitrogen=0.296803, rhoc=313.3, Pc=3395.8, M_Nitrogen=28.01348,Ttriple=63.151;
             //    K               kJ/kg-K             kg/m^3          kPa
static const double n[]={0,    
0.924803575275,//[1]
-0.492448489428,//[2]
0.661883336938,//[3]
-1.92902649201,//[4]
-0.0622469309629,//[5]
0.349943957581,//[6]
0.564857472498,//[7]
-1.61720005987,//[8]
-0.481395031883,//[9]
0.421150636384,//[10]
-0.0161962230825,//[11]
0.172100994165,//[12]
0.00735448924933,//[13]
0.0168077305479,//[14]
-0.00107626664179,//[15]
-0.0137318088513,//[16]
0.000635466899859,//[17]
0.00304432279419,//[18]
-0.0435762336045,//[19]
-0.0723174889316,//[20]
0.0389644315272,//[21]
-0.021220136391,//[22]
0.00408822981509,//[23]
-0.0000551990017984,//[24]
-0.0462016716479,//[25]
-0.00300311716011,//[26]
0.0368825891208,//[27]
-0.0025585684622,//[28]
0.00896915264558,//[29]
-0.0044151337035,//[30]
0.00133722924858,//[31]
0.000264832491957,//[32]
19.6688194015,//[33]
-20.911560073,//[34]
0.0167788306989,//[35]
2627.67566274//[36]
};

// d used for consistency with CO2 correlation (corresponds to i from Span)
static const int d[]={0,
1,//[1]
1,//[2]
2,//[3]
2,//[4]
3,//[5]
3,//[6]
1,//[7]
1,//[8]
1,//[9]
3,//[10]
3,//[11]
4,//[12]
6,//[13]
6,//[14]
7,//[15]
7,//[16]
8,//[17]
8,//[18]
1,//[19]
2,//[20]
3,//[21]
4,//[22]
5,//[23]
8,//[24]
4,//[25]
5,//[26]
5,//[27]
8,//[28]
3,//[29]
5,//[30]
6,//[31]
9,//[32]
1,//[33]
1,//[34]
3,//[35]
2//[36]
};

// t used for consistency with CO2 correlation (corresponds to j from Span)
static const double t[]={0.00,
0.25,//[1]
0.875,//[2]
0.5,//[3]
0.875,//[4]
0.375,//[5]
0.75,//[6]
0.5,//[7]
0.75,//[8]
2,//[9]
1.25,//[10]
3.5,//[11]
1,//[12]
0.5,//[13]
3,//[14]
0,//[15]
2.75,//[16]
0.75,//[17]
2.5,//[18]
4,//[19]
6,//[20]
6,//[21]
3,//[22]
3,//[23]
6,//[24]
16,//[25]
11,//[26]
15,//[27]
12,//[28]
12,//[29]
7,//[30]
4,//[31]
16,//[32]
0,//[33]
1,//[34]
2,//[35]
3//[36]
};

// c used for consistency with CO2 correlation (corresponds to l from Span)
static const int c[]={0,
0,//[1]
0,//[2]
0,//[3]
0,//[4]
0,//[5]
0,//[6]
1,//[7]
1,//[8]
1,//[9]
1,//[10]
1,//[11]
1,//[12]
1,//[13]
1,//[14]
1,//[15]
1,//[16]
1,//[17]
1,//[18]
2,//[19]
2,//[20]
2,//[21]
2,//[22]
2,//[23]
2,//[24]
3,//[25]
3,//[26]
3,//[27]
3,//[28]
4,//[29]
4,//[30]
4,//[31]
4,//[32]
2,//[33]
2,//[34]
2,//[35]
2//[36]
};

// alpha is used here for consistency with the definitions in R744.c upon which Nitrogen.c is based
// is phi_k from Span
static const double alpha[]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-32]
20,
20,
15,
25
};

static const double beta[]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-32]
325,
325,
300,
275
};

// epsilon is used here for consistency with the definitions in R744.c upon which Nitrogen.c is based
// is the value unity in Span
static const double epsilon[]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-32]
1,
1,
1,
1
};


// GAMMA is used here for consistency with the definitions in R744.c upon which Nitrogen.c is based
static const double GAMMA[]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-32]
1.16,
1.16,
1.13,
1.25
};



//Constants for ideal gas expression
static const double a0[]={0.0,
	2.5,
	-12.76952708,
	-0.00784163,
	-1.934819e-4,
	-1.247742e-5,
	6.678326e-8,
	1.012941,
	26.65788
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
				if (Tvec[i]>Tc || pvec[j]<psat_Nitrogen(Tvec[i]))
				{					
					rhomat[i][j]=get_Delta(Tvec[i],pvec[j])*rhoc;
					hmat[i][j]=h_Nitrogen(Tvec[i],rhomat[i][j],TYPE_Trho);
					smat[i][j]=s_Nitrogen(Tvec[i],rhomat[i][j],TYPE_Trho);
					umat[i][j]=u_Nitrogen(Tvec[i],rhomat[i][j],TYPE_Trho);
					cpmat[i][j]=cp_Nitrogen(Tvec[i],rhomat[i][j],TYPE_Trho);
					cvmat[i][j]=cv_Nitrogen(Tvec[i],rhomat[i][j],TYPE_Trho);
					cmat[i][j]=c_Nitrogen(Tvec[i],rhomat[i][j],TYPE_Trho);
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

double rho_Nitrogen(double T, double p, int Types)
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
double p_Nitrogen(double T, double rho)
{
	return Pressure_Trho(T,rho);
}
double h_Nitrogen(double T, double p_rho, int Types)
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
double s_Nitrogen(double T, double p_rho, int Types)
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
double u_Nitrogen(double T, double p_rho, int Types)
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
double cp_Nitrogen(double T, double p_rho, int Types)
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
double cv_Nitrogen(double T, double p_rho, int Types)
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
double w_Nitrogen(double T, double p_rho, int Types)
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
double rhosatL_Nitrogen(double T)
{
    const double ti[]={0,0.3294,2.0/3.0,8.0/3.0,35.0/6.0};
    const double Ni[]={0,1.48654237,-0.280476066,0.0894143085,-0.119879866};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+Ni[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(summer);
}

double rhosatV_Nitrogen(double T)
{
    const double ti[]={0,0.34,5.0/6.0,7.0/6.0,13.0/6.0,14.0/3.0};
    const double Ni[]={0,-1.70127164,-3.70402649,1.29859383,-0.561424977,-2.68505381};
    double summer=0;
    int i;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(Tc/T*summer);
}

double c_Nitrogen(double T, double p_rho, int Types)
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

double MM_Nitrogen(void)
{
	return M_Nitrogen;
}

double pcrit_Nitrogen(void)
{
	return Pc;
}

double Tcrit_Nitrogen(void)
{
	return Tc;
}

double Ttriple_Nitrogen(void)
{
	return Ttriple;
}

int errCode_Nitrogen(void)
{
	return errCode;
}

double psat_Nitrogen(double T)
{
    const double ti[]={0,1.0,1.5,2.5,5.0};
    const double Ni[]={0,-6.12445284,1.26327220,-0.765910082,-1.77570564};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+Ni[i]*pow(1-T/Tc,ti[i]);
    }
    return Pc*exp(Tc/T*summer);
}

double Tsat_Nitrogen(double P)
{
    double change,eps=.00005;
    int counter=1;
    double r1,r2,r3,T1,T2,T3;
 
    T1=275;
    T2=275+.01;
    r1=psat_Nitrogen(T1)-P;
    r2=psat_Nitrogen(T2)-P;
    
    // End at change less than 0.5%
    while(counter==1 || (fabs(change)/fabs(T2)>eps && counter<40))
    {
        T3=T2-0.5*r2/(r2-r1)*(T2-T1);
        r3=psat_Nitrogen(T3)-P;
        change=0.5*r2/(r2-r1)*(T2-T1);
        T1=T2;
        T2=T3;
        r1=r2;
        r2=r3;
        counter=counter+1;
    }
    return T3;
}   

double hsat_Nitrogen(double T, double x)
{
    double delta,tau;
    
    if (x>0.5)
    {
        delta=rhosatV_Nitrogen(T)/rhoc;
        tau=Tc/T;
        return R_Nitrogen*T*(1+tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+delta*dphir_dDelta(tau,delta));
    }
    else
    {
        delta=rhosatL_Nitrogen(T)/rhoc;
        tau=Tc/T;
        return R_Nitrogen*T*(1+tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+delta*dphir_dDelta(tau,delta));
    }   
}

double ssat_Nitrogen(double T, double x)
{
    double delta,tau;
    
    if (x>0.5)
    {
        delta=rhosatV_Nitrogen(T)/rhoc;
        tau=Tc/T;
        return R_Nitrogen*(tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))-phi0(tau,delta)-phir(tau,delta));
    }
    else
    {
        delta=rhosatL_Nitrogen(T)/rhoc;
        tau=Tc/T;
        return R_Nitrogen*(tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))-phi0(tau,delta)-phir(tau,delta));
    }   
}

double rhosat_Nitrogen(double T, double x)
{
    
    if (x>0.5)
        return rhosatV_Nitrogen(T);
    else
        return rhosatL_Nitrogen(T);
}


double visc_Nitrogen(double T, double p_rho, int Types)
{
	double e_k=98.94, //[K]
		   sigma=0.3656; //[nm]
	double rho,eta0,etar,OMEGA,delta,tau,Tstar;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,10.72,0.03989,0.001208,-7.402,4.620};
	double t[]={0,0.1,0.25,3.2,0.9,0.3};
	double d[]={0,2,10,12,2,1};
	double l[]={0,0,1,1,2,3};
	double g[]={0,0,1,1,1,1};


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

	eta0=0.0266958*sqrt(M_Nitrogen*T)/(sigma*sigma*OMEGA);
	etar=N[1]*pow(tau,t[1])*pow(delta,d[1])*exp(-g[1]*pow(delta,l[1]))
		+N[2]*pow(tau,t[2])*pow(delta,d[2])*exp(-g[2]*pow(delta,l[2]))
		+N[3]*pow(tau,t[3])*pow(delta,d[3])*exp(-g[3]*pow(delta,l[3]))
		+N[4]*pow(tau,t[4])*pow(delta,d[4])*exp(-g[4]*pow(delta,l[4]))
		+N[5]*pow(tau,t[5])*pow(delta,d[5])*exp(-g[5]*pow(delta,l[5]));

	return (eta0+etar)/1e6; // uPa-s to Pa-s
}

double X_tilde(double T,double tau,double delta)
{
	// X_tilde is dimensionless
	// Equation 11 slightly rewritten
	double drho_dp;
	drho_dp=1.0/(R_Nitrogen*T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*dphir2_dDelta2(tau,delta)));
	return Pc*delta/rhoc*drho_dp;
}

double k_Nitrogen(double T, double p_rho, int Types)
{
	double e_k=98.94, //[K]
		   sigma=0.3656, //[nm]
		   Tref=252.384, //[K]
		   zeta0=0.17, //[nm]
		   LAMBDA=0.055,
		   q_D=0.40; //[nm]
	double rho,eta0,OMEGA,delta,tau,Tstar,lambda0,lambdar,num,
		cp,cv,OMEGA_tilde,OMEGA_tilde0,zeta,nu,gamma,R0,lambdac,k,
		pi=3.141592654,mu;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,1.511,2.117,-3.332,8.862,31.11,-73.13,20.03,-0.7096,0.2672};
	double t[]={0,0,-1.0,-0.7,0.0,0.03,0.2,0.8,0.6,1.9};
	double d[]={0,0,0,0,1,2,3,4,8,10};
	double l[]={0,0,0,0,0,0,1,2,2,2};
	double g[]={0,0,0,0,0,0,1,1,1,1};
	
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

	eta0=0.0266958*sqrt(M_Nitrogen*T)/(sigma*sigma*OMEGA);
	lambda0=N[1]*eta0+N[2]*pow(tau,t[2])+N[3]*pow(tau,t[3]);

	lambdar=N[4]*pow(tau,t[4])*pow(delta,d[4])*exp(-g[4]*pow(delta,l[4]))
		   +N[5]*pow(tau,t[5])*pow(delta,d[5])*exp(-g[5]*pow(delta,l[5]))
		   +N[6]*pow(tau,t[6])*pow(delta,d[6])*exp(-g[6]*pow(delta,l[6]))
		   +N[7]*pow(tau,t[7])*pow(delta,d[7])*exp(-g[7]*pow(delta,l[7]))
	 	   +N[8]*pow(tau,t[8])*pow(delta,d[8])*exp(-g[8]*pow(delta,l[8]))
		   +N[9]*pow(tau,t[9])*pow(delta,d[9])*exp(-g[9]*pow(delta,l[9]));

	R0=1.01;
	nu=0.63;
	gamma=1.2415;
	k=1.380658e-23; //[J/K]

	num=X_tilde(T,Tc/T,delta)-X_tilde(Tref,Tc/Tref,delta)*Tref/T;

	// no critical enhancement if numerator of Eq. 10 is negative
	if (num<0)
		return (lambda0+lambdar)/1e6;

	cp=cp_Nitrogen(T,rho,TYPE_Trho);
	cv=cv_Nitrogen(T,rho,TYPE_Trho);
	mu=visc_Nitrogen(T,rho,TYPE_Trho)*1e6; //[uPa-s]

	zeta=zeta0*pow(num/LAMBDA,nu/gamma); //[nm]
	OMEGA_tilde=2.0/pi*((cp-cv)/cp*atan(zeta/q_D)+cv/cp*(zeta/q_D));
	OMEGA_tilde0=2.0/pi*(1.-exp(-1./(q_D/zeta+1.0/3.0*(zeta/q_D)*(zeta/q_D)/delta/delta)));
	lambdac=rho*(cp*1000.0)*k*R0*T/(6*pi*zeta*mu)*(OMEGA_tilde-OMEGA_tilde0)*1e18; // 1e18 is conversion to mW/m-K (not described in paper)

	return (lambda0+lambdar+lambdac)/1e6;
}

double dhdT_Nitrogen(double T, double p_rho, int Types)
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

double dhdrho_Nitrogen(double T, double p_rho, int Types)
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

double dpdT_Nitrogen(double T, double p_rho, int Types)
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

double dpdrho_Nitrogen(double T, double p_rho, int Types)
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
	return R_Nitrogen*T*rho*(1.0+delta*dphir_dDelta(tau,delta));
}
static double IntEnergy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return R_Nitrogen*T*tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta));
}
static double Enthalpy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return R_Nitrogen*T*(1+tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+delta*dphir_dDelta(tau,delta));
}
static double Entropy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return R_Nitrogen*(tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))-phi0(tau,delta)-phir(tau,delta));
}
static double SpecHeatV_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoc;
	tau=Tc/T;
	
	return -R_Nitrogen*powI(tau,2)*(dphi02_dTau2(tau,delta)+dphir2_dTau2(tau,delta));
}
static double SpecHeatP_Trho(double T, double rho)
{
	double delta,tau,c1,c2;
	delta=rho/rhoc;
	tau=Tc/T;

	c1=powI(1.0+delta*dphir_dDelta(tau,delta)-delta*tau*dphir2_dDelta_dTau(tau,delta),2);
    c2=(1.0+2.0*delta*dphir_dDelta(tau,delta)+powI(delta,2)*dphir2_dDelta2(tau,delta));
    return R_Nitrogen*(-powI(tau,2)*(dphi02_dTau2(tau,delta)+dphir2_dTau2(tau,delta))+c1/c2);
}

static double SpeedSound_Trho(double T, double rho)
{
	double delta,tau,c1,c2;
	delta=rho/rhoc;
	tau=Tc/T;

	c1=-SpecHeatV_Trho(T,rho)/R_Nitrogen;
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
	T=Tc/tau;   R=R_Nitrogen;
	//Note: dphi02_dDelta_dTau(tau,delta) is equal to zero
	return R*T/rhoc*(tau*(dphir2_dDelta_dTau(tau,delta))+dphir_dDelta(tau,delta)+delta*dphir2_dDelta2(tau,delta));
}
static double dhdT(double tau, double delta)
{
	double dhdT_rho,T,R,dhdtau;
	T=Tc/tau;   R=R_Nitrogen;
	dhdT_rho=R*tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+R*delta*dphir_dDelta(tau,delta)+R;
	dhdtau=R*T*(dphi0_dTau(tau,delta)+ dphir_dTau(tau,delta))+R*T*tau*(dphi02_dTau2(tau,delta)+dphir2_dTau2(tau,delta))+R*T*delta*dphir2_dDelta_dTau(tau,delta);
	return dhdT_rho+dhdtau*(-Tc/T/T);
}
static double dpdT(double tau, double delta)
{
	double T,R,rho;
	T=Tc/tau;   R=R_Nitrogen;  rho=delta*rhoc;
	return rho*R*(1+delta*dphir_dDelta(tau,delta)-delta*tau*dphir2_dDelta_dTau(tau,delta));
}
static double dpdrho(double tau, double delta)
{
	double T,R,rho;
	T=Tc/tau;   R=R_Nitrogen;  rho=delta*rhoc;
	return R*T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*dphir2_dDelta2(tau,delta));

}
/**************************************************/
/*          Private Property Functions            */
/**************************************************/

static double phir(double tau, double delta)
{ 
    
    int i;
    double phir=0,psi;
    
    for (i=1;i<=6;i++)
    {
        phir=phir+n[i]*powI(delta,d[i])*pow(tau,t[i]);
    }
    
    for (i=7;i<=32;i++)
    {
        phir=phir+n[i]*powI(delta,d[i])*pow(tau,t[i])*exp(-powI(delta,c[i]));
    }
    
    for (i=33;i<=36;i++)
    {
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        phir=phir+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi;
    }
    return phir;
}

static double dphir_dDelta(double tau, double delta)
{ 
    int i;
    double dphir_dDelta=0,psi;
    double di, ci;
    for (i=1;i<=6;i++)
    {
        di=(double)d[i];
        dphir_dDelta=dphir_dDelta+n[i]*di*powI(delta,d[i]-1)*pow(tau,t[i]);
    }
    for (i=7;i<=32;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir_dDelta=dphir_dDelta+n[i]*exp(-powI(delta,c[i]))*(powI(delta,d[i]-1)*pow(tau,t[i])*(di-ci*powI(delta,c[i])));
    }
    for (i=33;i<=36;i++)
    {
        di=(double)d[i];        
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir_dDelta=dphir_dDelta+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi*(di/delta-2.0*alpha[i]*(delta-epsilon[i]));
    }
    return dphir_dDelta;
}

static double dphir2_dDelta2(double tau, double delta)
{ 
    
    int i;
    double di,ci;
    double dphir2_dDelta2=0,psi;
    for (i=1;i<=6;i++)
    {
        di=(double)d[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*di*(di-1.0)*powI(delta,d[i]-2)*pow(tau,t[i]);
    }
    for (i=7;i<=32;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*exp(-powI(delta,c[i]))*(powI(delta,d[i]-2)*pow(tau,t[i])*( (di-ci*powI(delta,c[i]))*(di-1.0-ci*powI(delta,c[i])) - ci*ci*powI(delta,c[i])));
    }
    for (i=33;i<=36;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir2_dDelta2=dphir2_dDelta2+n[i]*pow(tau,t[i])*psi*(-2.0*alpha[i]*powI(delta,d[i])+4.0*powI(alpha[i],2)*powI(delta,d[i])*powI(delta-epsilon[i],2)-4.0*di*alpha[i]*powI(delta,d[i]-1)*(delta-epsilon[i])+di*(di-1.0)*powI(delta,d[i]-2));
    }
    return dphir2_dDelta2;
}

    
static double dphir2_dDelta_dTau(double tau, double delta)
{ 
    
    int i;
    double di, ci;
    double dphir2_dDelta_dTau=0,psi;

    for (i=1;i<=6;i++)
    {
        di=(double)d[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau + n[i]*di*t[i]*powI(delta,d[i]-1)*pow(tau,t[i]-1.0);
    }
    for (i=7;i<=32;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau + n[i]*exp(-powI(delta,c[i]))*powI(delta,d[i]-1)*t[i]*pow(tau,t[i]-1.0)*(di-ci*powI(delta,c[i]));
    }
    for (i=33;i<=36;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi*(di/delta-2.0*alpha[i]*(delta-epsilon[i]))*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    return dphir2_dDelta_dTau;
}

static double dphir_dTau(double tau, double delta)
{ 
    
    int i;
    double dphir_dTau=0,psi;
    
    for (i=1;i<=6;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powI(delta,d[i])*pow(tau,t[i]-1.0);
    }
    for (i=7;i<=32;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powI(delta,d[i])*pow(tau,t[i]-1.0)*exp(-powI(delta,c[i]));
    }
    for (i=33;i<=36;i++)
    {
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir_dTau=dphir_dTau+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    return dphir_dTau;
}


static double dphir2_dTau2(double tau, double delta)
{ 
    
    int i;
    double dphir2_dTau2=0,psi;
    
    for (i=1;i<=6;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powI(delta,d[i])*pow(tau,t[i]-2.0);
    }
    for (i=7;i<=32;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powI(delta,d[i])*pow(tau,t[i]-2.0)*exp(-powI(delta,c[i]));
    }
    for (i=33;i<=36;i++)
    {
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir2_dTau2=dphir2_dTau2+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi*(powI(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]),2)-t[i]/powI(tau,2)-2.0*beta[i]);
    }
    return dphir2_dTau2;
}

/* Maxima code for derivatives of ideal gas residual Helmholtz:
A:log(delta)+a0[1]*log(tau)+a0[2]+a0[3]*tau+a0[4]/tau+a0[5]/tau/tau+a0[6]/tau/tau/tau+a0[7]*log(1-exp(-a0[8]*tau));
diff(A,delta);
diff(%,delta);
diff(A,tau);
diff(%,tau);
*/

static double phi0(double tau, double delta)
{
    double phi0=0;
    
    phi0=log(delta)+a0[1]*log(tau)+a0[2]+a0[3]*tau+a0[4]/tau+a0[5]/tau/tau+a0[6]/tau/tau/tau+a0[7]*log(1-exp(-a0[8]*tau));
    return phi0;
}

static double dphi0_dDelta(double tau, double delta)
{
    return 1/delta;
}

static double dphi02_dDelta2(double tau, double delta)
{
    return -1.0/powI(delta,2);
}

static double dphi0_dTau(double tau, double delta)
{
    return a0[1]/tau+a0[3]-a0[4]/tau/tau-2.0*a0[5]/tau/tau/tau-3*a0[6]/tau/tau/tau/tau+(a0[7]*a0[8]*exp(-a0[8]*tau))/(1-exp(-a0[8]*tau));
}

static double dphi02_dTau2(double tau, double delta)
{
    return -(a0[7]*a0[8]*a0[8]*exp(-a0[8]*tau))/(1-exp(-a0[8]*tau))-(a0[7]*a0[8]*a0[8]*exp(-2*a0[8]*tau))/(1-exp(-a0[8]*tau))/(1-exp(-a0[8]*tau))-a0[1]/tau/tau+(2*a0[4])/powI(tau,3)+(6*a0[5])/powI(tau,4)+(12*a0[6])/powI(tau,5);
}

static double get_Delta(double T, double P)
{
    double change,eps=1e-8;
    int counter=1;
    double r1,r2,r3,delta1,delta2,delta3;
    double tau;
    double delta_guess;
 
    if (P>Pc)
    {
        if (T>Tc)
        {
            delta_guess=P/(R_Nitrogen*T)/rhoc;
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
		delta_guess=P/(R_Nitrogen*T)/rhoc;
	}
	else if(P<psat_Nitrogen(T))
        {
		//Superheated vapor
		delta_guess=(rhosatV_Nitrogen(T)*P/psat_Nitrogen(T))/rhoc;
        }
	else 
	{
		delta_guess=(rhosatL_Nitrogen(T))/rhoc;
	}
    }
    
   
    tau=Tc/T;
    delta1=delta_guess;
    delta2=delta_guess+.00001;
    r1=P/(delta1*rhoc*R_Nitrogen*T)-1.0-delta1*dphir_dDelta(tau,delta1);
    r2=P/(delta2*rhoc*R_Nitrogen*T)-1.0-delta2*dphir_dDelta(tau,delta2);
    
    // End at change less than 0.05%
    while(counter==1 || (fabs(r2)/delta2>eps && counter<40))
    {
        delta3=delta2-r2/(r2-r1)*(delta2-delta1);
        r3=P/(delta3*rhoc*R_Nitrogen*T)-1.0-delta3*dphir_dDelta(tau,delta3);
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





