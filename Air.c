/* Properties of Air
by Ian Bell


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
#include "PropErrorCodes.h"
#include "PropMacros.h"
#include "Air.h"
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

static const double Tj=132.6312, R_Air=0.287117125828, rhoj=302.5507652, pj=3785.02, M_Air=28.9586, Ttriple=59.75;
             //           K             kJ/kg-K       kg/m^3     kPa            kg/kmol          K
static const double Tc=132.5306,pc=3786.0;

static const double N[]={0,
 0.118160747229,//[1]
 0.713116392079,//[2]
-0.161824192067e1,//[3]
 0.714140178971e-1,//[4]
-0.865421396646e-1,//[5]
 0.134211176704,//[6]
 0.112626704218e-1,//[7]
-0.420533228842e-1,//[8]
 0.349008431982e-1,//[9]
 0.164957183186e-3,//[10]
-0.101365037912,//[11]
-0.173813690970,//[12]
-0.472103183731e-1,//[13]
-0.122523554253e-1,//[14]
-0.146629609713,//[15]
-0.316055879821e-1,//[16]
 0.233594806142e-3,//[17]
 0.148287891978e-1,//[18]
-0.938782884667e-2//[19]
};

static const int i[]={0,
1,//[1]
1,//[2]
1,//[3]
2,//[4]
3,//[5]
3,//[6]
4,//[7]
4,//[8]
4,//[9]
6,//[10]
1,//[11]
3,//[12]
5,//[13]
6,//[14]
1,//[15]
3,//[16]
11,//[17]
1,//[18]
3//[19]
};

static const double j[]={0.00,
0,//[1]
0.33,//[2]
1.01,//[3]
0,//[4]
0,//[5]
0.15,//[6]
0,//[7]
0.2,//[8]
0.35,//[9]
1.35,//[10]
1.6,//[11]
0.8,//[12]
0.95,//[13]
1.25,//[14]
3.6,//[15]
6,//[16]
3.25,//[17]
3.5,//[18]
15//[19]
};

static const int l[]={
0,0,0,0,0,0,0,0,0,0,0, // indices [0-10]
1,//[11]
1,//[12]
1,//[13]
1,//[14]
2,//[15]
2,//[16]
2,//[17]
3,//[18]
3,//[19]
};

//Constants for ideal gas expression
static const double N0[]={0.0,
 0.605719400e-7,//[1]
-0.210274769e-4,//[2]
-0.158860716e-3,//[3]
-13.841928076,//[4]
 17.275266575,//[5]
-0.195363420e-3,//[6]
 2.490888032,//[7]
 0.791309509,//[8]
 0.212236768,//[9]
-0.197938904,//[10]
 25.36365,//[11]
 16.90741,//[12]
 87.31279//[13]
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
				if (Tvec[i]>Tj || pvec[j]<pbp_Air(Tvec[i]))
				{					
					rhomat[i][j]=get_Delta(Tvec[i],pvec[j])*rhoj;
					hmat[i][j]=h_Air(Tvec[i],rhomat[i][j],TYPE_Trho);
					smat[i][j]=s_Air(Tvec[i],rhomat[i][j],TYPE_Trho);
					umat[i][j]=u_Air(Tvec[i],rhomat[i][j],TYPE_Trho);
					cpmat[i][j]=cp_Air(Tvec[i],rhomat[i][j],TYPE_Trho);
					cvmat[i][j]=cv_Air(Tvec[i],rhomat[i][j],TYPE_Trho);
					cmat[i][j]=c_Air(Tvec[i],rhomat[i][j],TYPE_Trho);
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

double rho_Air(double T, double p, int Types)
{
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_TPNoLookup:
			return get_Delta(T,p)*rhoj;
		case TYPE_TP:
			BuildLookup();
			return LookupValue("rho",T,p);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double p_Air(double T, double rho)
{
	return Pressure_Trho(T,rho);
}
double h_Air(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return Enthalpy_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoj;
			return Enthalpy_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("h",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double s_Air(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return Entropy_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoj;
			return Entropy_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("s",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double u_Air(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return IntEnergy_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoj;
			return IntEnergy_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("u",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double cp_Air(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return SpecHeatP_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoj;
			return SpecHeatP_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("cp",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double cv_Air(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return SpecHeatV_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoj;
			return SpecHeatV_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("cv",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double w_Air(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return SpeedSound_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoj;
			return SpeedSound_Trho(T,rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}

double rhosatL_Air(double T)
{
	const double ti[]={0,0.65,0.85,0.95,1.1};
    const double Ni[]={0,43.3413,-240.073,285.139,-88.3366,-0.892181};
    double summer=0; int k;
    for (k=1;k<=4;k++)
    {
        summer=summer+Ni[k]*pow(1.0-T/Tj,ti[k]);
    }
    return rhoj*(1+summer+Ni[5]*log(T/Tj));
}

double rhosatV_Air(double T)
{
	const double ti[]={0,0.41,1.0,2.8,6.5};
    const double Ni[]={0,-2.0466,-4.7520,-13.259,-47.652};
    double summer=0; int k;
    for (k=1;k<=4;k++)
    {
        summer=summer+Ni[k]*pow(1.0-T/Tj,ti[k]);
    }
    return rhoj*exp(summer);
}

double c_Air(double T, double p_rho, int Types)
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
			rho=get_Delta(T,p_rho)*rhoj;
			return SpeedSound_Trho(T,rho);
		}
	}
}

double MM_Air(void)
{
	return M_Air;
}

double pcrit_Air(void)
{
	return pc;
}

double Tcrit_Air(void)
{
	return Tc;
}
double Ttriple_Air(void)
{
	return Ttriple;
}
double rhojrit_Air(void)
{
	return rhoj;
}
int errCode_Air(void)
{
	return errCode;
}

double pbp_Air(double T)
{
	const double Ni[]={0,0.2260724,-7.080499,5.700283,-12.44017,17.81926,-10.81364};
	double summer=0; int k;
    for (k=1;k<=6;k++)
    {
        summer=summer+Ni[k]*pow(1-T/Tj,(double)k/2.0);
    }
	return pj*exp(Tj/T*summer);
}

double pdp_Air(double T)
{
	const double Ni[]={0,-0.1567266,-5.539635,0,0,0.7567212,0,0,-3.514322};
	double summer=0; int k;
    for (k=1;k<=8;k++)
    {
        summer=summer+Ni[k]*pow(1-T/Tj,(double)k/2.0);
    }
	return pj*exp(Tj/T*summer);
}

double Tsat_Air(double P)
{
    double change,eps=.00005;
    int counter=1;
    double r1,r2,r3,T1,T2,T3;
 
    printf("Warning!!! Not properly implemented for dew and bubble\n");
    T1=275;
    T2=275+.01;
    r1=pdp_Air(T1)-P;
    r2=pdp_Air(T2)-P;
    
    // End at change less than 0.5%
    while(counter==1 || (fabs(change)/fabs(T2)>eps && counter<40))
    {
        T3=T2-0.5*r2/(r2-r1)*(T2-T1);
        r3=pdp_Air(T3)-P;
        change=0.5*r2/(r2-r1)*(T2-T1);
        T1=T2;
        T2=T3;
        r1=r2;
        r2=r3;
        counter=counter+1;
    }
    return T3;
}   

double rhosat_Air(double T, double x)
{
    
    if (x>0.5)
        return rhosatV_Air(T);
    else
        return rhosatL_Air(T);
}


double visc_Air(double T, double p_rho, int Types)
{
	double e_k=103.3, //[K]
		   sigma=0.360; //[nm]
	double rho,eta0,etar,OMEGA,delta,tau,Tstar;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,10.72,1.122,0.002019,-8.876,-0.02916};
	double t[]={0,0.2,0.05,2.4,0.6,3.6};
	double d[]={0,1,4,9,1,8};
	double l[]={0,0,0,0,1,1};
	double g[]={0,0,0,0,1,1};


	if (Types==TYPE_Trho)
		rho=p_rho;
	else if (Types==TYPE_TPNoLookup)
	{
		rho=get_Delta(T,p_rho)*rhoj;
	}
	else if (Types==TYPE_TP)
	{
		BuildLookup();
		return LookupValue("d",T,p_rho);
	}
	delta=rho/rhoj;
	tau=Tj/T;
	Tstar=T/(e_k);
	OMEGA=exp(b[0]*powI(log(Tstar),0)
			 +b[1]*powI(log(Tstar),1)
		     +b[2]*powI(log(Tstar),2)
			 +b[3]*powI(log(Tstar),3)
		     +b[4]*powI(log(Tstar),4));

	eta0=0.0266958*sqrt(M_Air*T)/(sigma*sigma*OMEGA);
	etar=N[1]*pow(tau,t[1])*pow(delta,d[1])*exp(-g[1]*pow(delta,l[1]))
		+N[2]*pow(tau,t[2])*pow(delta,d[2])*exp(-g[2]*pow(delta,l[2]))
		+N[3]*pow(tau,t[3])*pow(delta,d[3])*exp(-g[3]*pow(delta,l[3]))
		+N[4]*pow(tau,t[4])*pow(delta,d[4])*exp(-g[4]*pow(delta,l[4]))
		+N[5]*pow(tau,t[5])*pow(delta,d[5])*exp(-g[5]*pow(delta,l[5]));

	return (eta0+etar)/1e6; // uPa-s to Pa-s
}

static double X_tilde(double T,double tau,double delta)
{
	// X_tilde is dimensionless
	// Equation 11 slightly rewritten
	double drho_dp;
	drho_dp=1.0/(R_Air*T*(1+2*delta*dphir_dDelta_Air(tau,delta)+delta*delta*dphir2_dDelta2_Air(tau,delta)));
	return pj*delta/rhoj*drho_dp;
}

double k_Air(double T, double p_rho, int Types)
{
	double e_k=103.3, //[K]
		   sigma=0.360, //[nm]
		   Tref=265.262, //[K]
		   zeta0=0.11, //[nm]
		   LAMBDA=0.055,
		   q_D=0.31; //[nm]
	double rho,eta0,OMEGA,delta,tau,Tstar,lambda0,lambdar,num,
		cp,cv,OMEGA_tilde,OMEGA_tilde0,zeta,nu,gamma,R0,lambdac,k,
		pi=3.141592654,mu;
	double b[]={0.431,-0.4623,0.08406,0.005341,-0.00331};

	double N[]={0,1.308,1.405,-1.036,8.743,14.76,-16.62,3.793,-6.142,-0.3778};
	double t[]={0,0,-1.1,-0.3,0.1,0.0,0.5,2.7,0.3,1.3};
	double d[]={0,0,0,0,1,2,3,7,7,11};
	double l[]={0,0,0,0,0,0,2,2,2,2};
	double g[]={0,0,0,0,0,0,0,1,1,1};
	
	if (Types==TYPE_Trho)
		rho=p_rho;
	else if (Types==TYPE_TPNoLookup)
	{
		rho=get_Delta(T,p_rho)*rhoj;
	}
	else if (Types==TYPE_TP)
	{
		BuildLookup();
		return LookupValue("D",T,p_rho);
	}
	delta=rho/rhoj;
	tau=Tj/T;
	Tstar=T/(e_k);

	OMEGA=exp(b[0]*powI(log(Tstar),0)
			 +b[1]*powI(log(Tstar),1)
		     +b[2]*powI(log(Tstar),2)
			 +b[3]*powI(log(Tstar),3)
		     +b[4]*powI(log(Tstar),4));

	eta0=0.0266958*sqrt(M_Air*T)/(sigma*sigma*OMEGA);
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

	num=X_tilde(T,Tj/T,delta)-X_tilde(Tref,Tj/Tref,delta)*Tref/T;

	// no critical enhancement if numerator of Eq. 10 is negative
	if (num<0)
		return (lambda0+lambdar)/1e6;

	cp=cp_Air(T,rho,TYPE_Trho);
	cv=cv_Air(T,rho,TYPE_Trho);
	mu=visc_Air(T,rho,TYPE_Trho)*1e6; //[uPa-s]

	zeta=zeta0*pow(num/LAMBDA,nu/gamma); //[nm]
	OMEGA_tilde=2.0/pi*((cp-cv)/cp*atan(zeta/q_D)+cv/cp*(zeta/q_D));
	OMEGA_tilde0=2.0/pi*(1.-exp(-1./(q_D/zeta+1.0/3.0*(zeta/q_D)*(zeta/q_D)/delta/delta)));
	lambdac=rho*(cp*1000.0)*k*R0*T/(6*pi*zeta*mu)*(OMEGA_tilde-OMEGA_tilde0)*1e18; // 1e18 is conversion to mW/m-K (not described in paper)

	return (lambda0+lambdar+lambdac)/1e6;
}

double B_Air(double tau)
{
	// given by B*rhoc=lim(delta --> 0) [dphir_ddelta(tau)]
	return 1.0/rhoj*dphir_dDelta_Air(tau,1e-12);
}

double dBdT_Air(double tau)
{
	// given by B*rhoc^2=lim(delta --> 0) [dphir2_ddelta2(tau)]
	return -1.0/rhoj*tau*tau/Tj*dphir2_dDelta_dTau_Air(tau,1e-12);
}

double C_Air(double tau)
{
	// given by B*rhoc^2=lim(delta --> 0) [dphir2_ddelta2(tau)]
	return 1.0/(rhoj*rhoj)*dphir2_dDelta2_Air(tau,1e-12);
}

double dCdT_Air(double tau)
{
	// given by B*rhoc^2=lim(delta --> 0) [dphir2_ddelta2(tau)]
	return -1.0/(rhoj*rhoj)*tau*tau/Tj*dphir3_dDelta2_dTau_Air(tau,1e-12);
}

double dhdT_Air(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho;
			return dhdT(Tj/T,rho/rhoj);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoj;
			return dhdT(Tj/T,rho/rhoj);
		case 99:
			rho=get_Delta(T,p_rho)*rhoj;
			return (Enthalpy_Trho(T+0.001,rho)-Enthalpy_Trho(T,rho))/0.001;
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}

double dhdrho_Air(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho;
			return dhdrho(Tj/T,rho/rhoj);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoj;
			return dhdrho(Tj/T,rho/rhoj);
		case 99:
			rho=get_Delta(T,p_rho)*rhoj;
			return (Enthalpy_Trho(T,rho+0.001)-Enthalpy_Trho(T,rho))/0.001;
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}

double dpdT_Air(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho;
			return dpdT(Tj/T,rho/rhoj);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoj;
			return dpdT(Tj/T,rho/rhoj);
		case 99:
			rho=get_Delta(T,p_rho)*rhoj;
			return (Pressure_Trho(T+0.01,rho)-Pressure_Trho(T,rho))/0.01;
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}

double dpdrho_Air(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho;
			return dpdrho(Tj/T,rho/rhoj);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhoj;
			return dpdrho(Tj/T,rho/rhoj);
		case 99:
			rho=get_Delta(T,p_rho)*rhoj;
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
	delta=rho/rhoj;
	tau=Tj/T;
	return R_Air*T*rho*(1.0+delta*dphir_dDelta_Air(tau,delta));
}
static double IntEnergy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoj;
	tau=Tj/T;
	
	return R_Air*T*tau*(dphi0_dTau_Air(tau,delta)+dphir_dTau_Air(tau,delta));
}
static double Enthalpy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoj;
	tau=Tj/T;
	
	return R_Air*T*(1+tau*(dphi0_dTau_Air(tau,delta)+dphir_dTau_Air(tau,delta))+delta*dphir_dDelta_Air(tau,delta));
}
static double Entropy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoj;
	tau=Tj/T;
	
	return R_Air*(tau*(dphi0_dTau_Air(tau,delta)+dphir_dTau_Air(tau,delta))-phi0_Air(tau,delta)-phir_Air(tau,delta));
}
static double SpecHeatV_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhoj;
	tau=Tj/T;
	
	return -R_Air*powI(tau,2)*(dphi02_dTau2_Air(tau,delta)+dphir2_dTau2_Air(tau,delta));
}
static double SpecHeatP_Trho(double T, double rho)
{
	double delta,tau,c1,c2;
	delta=rho/rhoj;
	tau=Tj/T;

	c1=powI(1.0+delta*dphir_dDelta_Air(tau,delta)-delta*tau*dphir2_dDelta_dTau_Air(tau,delta),2);
    c2=(1.0+2.0*delta*dphir_dDelta_Air(tau,delta)+powI(delta,2)*dphir2_dDelta2_Air(tau,delta));
    return R_Air*(-powI(tau,2)*(dphi02_dTau2_Air(tau,delta)+dphir2_dTau2_Air(tau,delta))+c1/c2);
}

static double SpeedSound_Trho(double T, double rho)
{
	double delta,tau,c1,c2;
	delta=rho/rhoj;
	tau=Tj/T;

	c1=-SpecHeatV_Trho(T,rho)/R_Air;
	c2=(1.0+2.0*delta*dphir_dDelta_Air(tau,delta)+powI(delta,2)*dphir2_dDelta2_Air(tau,delta));
    return sqrt(-c2*T*SpecHeatP_Trho(T,rho)*1000/c1);
}

/**************************************************/
/*           Property Derivatives                 */
/**************************************************/
// See Lemmon, 2000 for more information
static double dhdrho(double tau, double delta)
{
	double T,R;
	T=Tj/tau;   R=R_Air;
	//Note: dphi02_dDelta_dTau(tau,delta) is equal to zero
	return R*T/rhoj*(tau*(dphir2_dDelta_dTau_Air(tau,delta))+dphir_dDelta_Air(tau,delta)+delta*dphir2_dDelta2_Air(tau,delta));
}
static double dhdT(double tau, double delta)
{
	double dhdT_rho,T,R,dhdtau;
	T=Tj/tau;   R=R_Air;
	dhdT_rho=R*tau*(dphi0_dTau_Air(tau,delta)+dphir_dTau_Air(tau,delta))+R*delta*dphir_dDelta_Air(tau,delta)+R;
	dhdtau=R*T*(dphi0_dTau_Air(tau,delta)+ dphir_dTau_Air(tau,delta))+R*T*tau*(dphi02_dTau2_Air(tau,delta)+dphir2_dTau2_Air(tau,delta))+R*T*delta*dphir2_dDelta_dTau_Air(tau,delta);
	return dhdT_rho+dhdtau*(-Tj/T/T);
}
static double dpdT(double tau, double delta)
{
	double T,R,rho;
	T=Tj/tau;   R=R_Air;  rho=delta*rhoj;
	return rho*R*(1+delta*dphir_dDelta_Air(tau,delta)-delta*tau*dphir2_dDelta_dTau_Air(tau,delta));
}
static double dpdrho(double tau, double delta)
{
	double T,R,rho;
	T=Tj/tau;   R=R_Air;  rho=delta*rhoj;
	return R*T*(1+2*delta*dphir_dDelta_Air(tau,delta)+delta*delta*dphir2_dDelta2_Air(tau,delta));
}

/**************************************************/
/*          Private Property Functions            */
/**************************************************/

double phir_Air(double tau, double delta)
{ 
    int k;
    double phir=0;
    
    for (k=1;k<=10;k++)
    {
        phir=phir+N[k]*powI(delta,i[k])*pow(tau,j[k]);
    }
    
    for (k=11;k<=19;k++)
    {
        phir=phir+N[k]*powI(delta,i[k])*pow(tau,j[k])*exp(-powI(delta,l[k]));
    }
    return phir;
}

double dphir_dDelta_Air(double tau, double delta)
{ 
    int k;
    double dphir_dDelta=0;
    for (k=1;k<=10;k++)
    {
        dphir_dDelta+=i[k]*N[k]*powI(delta,i[k]-1)*pow(tau,j[k]);
    }
    for (k=11;k<=19;k++)
    {
        dphir_dDelta+=N[k]*powI(delta,i[k]-1)*pow(tau,j[k])*exp(-powI(delta,l[k]))*(i[k]-l[k]*powI(delta,l[k]));
    }
    return dphir_dDelta;
}

double dphir2_dDelta2_Air(double tau, double delta)
{ 
    int k;
    double di,ci;
    double dphir2_dDelta2=0;
    for (k=1;k<=10;k++)
    {
        di=(double)i[k];
        dphir2_dDelta2=dphir2_dDelta2+N[k]*di*(di-1.0)*powI(delta,i[k]-2)*pow(tau,j[k]);
    }
    for (k=11;k<=19;k++)
    {
        di=(double)i[k];
        ci=(double)l[k];
        dphir2_dDelta2=dphir2_dDelta2+N[k]*exp(-powI(delta,l[k]))*(powI(delta,i[k]-2)*pow(tau,j[k])*( (di-ci*powI(delta,l[k]))*(di-1.0-ci*powI(delta,l[k])) - ci*ci*powI(delta,l[k])));
    }
    return dphir2_dDelta2;
}

    
double dphir2_dDelta_dTau_Air(double tau, double delta)
{ 
    int k;
    double dphir2_dDelta_dTau=0;

    for (k=1;k<=10;k++)
    {
        dphir2_dDelta_dTau+=i[k]*j[k]*N[k]*powI(delta,i[k]-1)*pow(tau,j[k]-1.0);
    }
    for (k=11;k<=19;k++)
    {
        dphir2_dDelta_dTau+=j[k]*N[k]*powI(delta,i[k]-1)*pow(tau,j[k]-1.0)*exp(-powI(delta,l[k]))*(i[k]-l[k]*powI(delta,l[k]));
    }
    return dphir2_dDelta_dTau;
}

double dphir3_dDelta2_dTau_Air(double tau, double delta)
{ 
    int k;
    double dphir3_dDelta2_dTau=0;

    for (k=1;k<=10;k++)
    {
        dphir3_dDelta2_dTau+=i[k]*(i[k]-1)*j[k]*N[k]*powI(delta,i[k]-2)*pow(tau,j[k]-1.0);
    }
    for (k=11;k<=19;k++)
    {
        dphir3_dDelta2_dTau+=j[k]*N[k]*powI(delta,i[k]-2)*pow(tau,j[k]-1.0)*exp(-powI(delta,l[k]))*((i[k]-l[k]*powI(delta,l[k]))*(i[k]-1-l[k]*powI(delta,l[k]))-l[k]*l[k]*powI(delta,l[k]));
    }
    return dphir3_dDelta2_dTau;
}

double dphir_dTau_Air(double tau, double delta)
{ 
    
    int k;
    double dphir_dTau=0;
    
    for (k=1;k<=10;k++)
    {
        dphir_dTau=dphir_dTau+N[k]*j[k]*powI(delta,i[k])*pow(tau,j[k]-1.0);
    }
    for (k=11;k<=19;k++)
    {
        dphir_dTau=dphir_dTau+N[k]*j[k]*powI(delta,i[k])*pow(tau,j[k]-1.0)*exp(-powI(delta,l[k]));
    }
    return dphir_dTau;
}


double dphir2_dTau2_Air(double tau, double delta)
{ 
    
    int k;
    double dphir2_dTau2=0;
    
    for (k=1;k<=10;k++)
    {
        dphir2_dTau2=dphir2_dTau2+N[k]*j[k]*(j[k]-1.0)*powI(delta,i[k])*pow(tau,j[k]-2.0);
    }
    for (k=11;k<=19;k++)
    {
        dphir2_dTau2=dphir2_dTau2+N[k]*j[k]*(j[k]-1.0)*powI(delta,i[k])*pow(tau,j[k]-2.0)*exp(-powI(delta,l[k]));
    }
    return dphir2_dTau2;
}

double phi0_Air(double tau, double delta)
{
    double phi0=0;
    int k;
    
    phi0=log(delta);
    for (k=1;k<=5;k++)
    {
        phi0=phi0+N0[k]*powI(tau,k-4);
    }
    phi0+=N0[6]*pow(tau,1.5)+N0[7]*log(tau);
    for (k=8;k<=9;k++)
    {
        phi0+=N0[k]*log(1.0-exp(-N0[k+3]*tau));
    }
    phi0+=N0[10]*log(2.0/3.0+exp(N0[13]*tau));
    return phi0;
}

double dphi0_dDelta_Air(double tau, double delta)
{
    return 1/delta;
}

double dphi02_dDelta2_Air(double tau, double delta)
{
    return -1.0/powI(delta,2);
}

double dphi0_dTau_Air(double tau, double delta)
{
    double dphi0_dTau=0; int k;
    for (k=1;k<=5;k++)
    {
        dphi0_dTau+=N0[k]*(k-4)*powI(tau,k-5);
    }
    dphi0_dTau+=1.5*N0[6]*sqrt(tau)+N0[7]/tau;
    dphi0_dTau+=N0[8]*N0[11]/(exp(N0[11]*tau)-1)+N0[9]*N0[12]/(exp(N0[12]*tau)-1);
    dphi0_dTau+=N0[10]*N0[13]/(2.0/3.0*exp(-N0[13]*tau)+1);
    
    return dphi0_dTau;
}

double dphi02_dTau2_Air(double tau, double delta)
{
    double dphi02_dTau2=0;
    int k;
    
    for (k=1;k<=5;k++)
    {
        dphi02_dTau2+=N0[k]*(k-4)*(k-5)*powI(tau,k-6);
    }
    dphi02_dTau2+=0.75*N0[6]*pow(tau,-0.5)-N0[7]/(tau*tau);
    dphi02_dTau2+=-N0[8]*N0[11]*N0[11]*exp(N0[11]*tau)/powI(exp(N0[11]*tau)-1,2);
    dphi02_dTau2+=-N0[9]*N0[12]*N0[12]*exp(N0[12]*tau)/powI(exp(N0[12]*tau)-1,2);
    dphi02_dTau2+=-(2.0/3.0)*N0[10]*N0[13]*N0[13]*exp(-N0[13]*tau)/powI((2.0/3.0)*exp(-N0[13]*tau)+1,2);
    return dphi02_dTau2;
}

static double get_Delta(double T, double P)
{
    double change,eps=.00005;
    int counter=1;
    double r1,r2,r3,delta1,delta2,delta3;
    double tau;
    double delta_guess;
 
    if (P>pc)
    {
        if (T>Tj)
        {
            delta_guess=P/(R_Air*T)/rhoj;
        }
        else
        {
            delta_guess=1000/rhoj;
        }
    }
    else
    {
        if (T>Tj)
	{
		// Supercritical (try ideal gas)
		delta_guess=P/(R_Air*T)/rhoj;
	}
	else if(P<pdp_Air(T))
        {
		//Superheated vapor
		delta_guess=(rhosatV_Air(T)*P/pdp_Air(T))/rhoj;
        }
	else //(T<Tsat_Air(P))
	{
		delta_guess=(rhosatL_Air(T))/rhoj;
	}
    }
    tau=Tj/T;
    delta1=delta_guess;
    delta2=delta_guess+.00001;
    r1=P/(delta1*rhoj*R_Air*T)-1.0-delta1*dphir_dDelta_Air(tau,delta1);
    r2=P/(delta2*rhoj*R_Air*T)-1.0-delta2*dphir_dDelta_Air(tau,delta2);
    
    // End at change less than 0.05%
    while(counter==1 || (fabs(r2)/delta2>eps && counter<40))
    {
        delta3=delta2-r2/(r2-r1)*(delta2-delta1);
        r3=P/(delta3*rhoj*R_Air*T)-1.0-delta3*dphir_dDelta_Air(tau,delta3);
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
    int k;
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
    for (k=1;k<y_in;k++)
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





