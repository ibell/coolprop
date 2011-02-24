/*
Properties for R410A.  
by Ian Bell

Pseudo-pure fluid thermo props from 
"Pseudo-pure fluid Equations of State for the Refrigerant Blends R410A, R404A, R507C and R407C" 
by E.W. Lemmon, Int. J. Thermophys. v. 24, n4, 2003

In order to call the exposed functions, rho_, h_, s_, cp_,...... 
there are three different ways the inputs can be passed, and this is expressed by the Types integer flag.  
These macros are defined in the PropMacros.h header file:
1) First parameter temperature, second parameter pressure ex: h_R410A(260,354.7,1)=274
	In this case, the lookup tables are built if needed and then interpolated
2) First parameter temperature, second parameter density ex: h_R410A(260,13.03,2)=274
	Density and temp plugged directly into EOS
3) First parameter temperature, second parameter pressure ex: h_R410A(260,354.7,3)=274
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
#include "R410A.h"

static int errCode;
static char errStr[ERRSTRLENGTH];

#define nP 200
#define nT 200

static double Tmin=150, Tmax=550, Pmin=70.03, Pmax=6000;
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

static const double a[]={
    36.8871,		//[0]
    7.15807,		//[1]
    -46.87575,		//[2]
    2.0623,			//[3]
    5.9751,			//[4]
    1.5612			//[5]
};

static const double b[]={
    0,				//[0]
    0,				//[1]
    -0.1,			//[2]
    2.02326,		//[3]
    5.00154,		//[4]
    11.2484			//[5]
};

static const double N[]={
	 0.0,			//[0]
     0.987252,		//[1]
    -1.03017,		//[2]
     1.17666,		//[3]
    -0.138991,		//[4]
     0.00302373,	//[5]
    -2.53639,		//[6]
    -1.96680,		//[7]
    -0.830480,		//[8]
     0.172477,		//[9]
    -0.261116,		//[10]
    -0.0745473,		//[11]
     0.679757,		//[12]
    -0.652431,		//[13]
     0.0553849,		//[14]
    -0.0710970,		//[15]
    -0.000875332,	//[16]
     0.0200760,		//[17]
    -0.0139761,		//[18]
    -0.0185110,		//[19]
     0.0171939,		//[20]
    -0.00482049		//[21]
};

static const double j[]={
	0.0,			//[0]
    0.44,			//[1]
    1.2,			//[2]
    2.97,			//[3]
    2.95,			//[4]
    0.2,			//[5]
    1.93,			//[6]
    1.78,			//[7]
    3.0,			//[8]
    0.2,			//[9]
    0.74,			//[10]
    3.0,			//[11]
    2.1,			//[12]
    4.3,			//[13]
    0.25,			//[14]
    7.0,			//[15]
    4.7,			//[16]
    13.0,			//[17]
    16.0,			//[18]
    25.0,			//[19]
    17.0,			//[20]
    7.4				//[21]
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
    3,				//[8]
    5,				//[9]
    5,				//[10]
    5,				//[11]
    1,				//[12]
    1,				//[13]
    4,				//[14]
    4,				//[15]
    9,				//[16]
    2,				//[17]
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
    2,				//[12]
    2,				//[13]
    2,				//[14]
    2,				//[15]
    2,				//[16]
    3,				//[17]
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
	-7.2818,		//[1]
    2.5093,			//[2]
    -3.2695,		//[3]
    -2.8022			//[4]
};
    
static const double tbp[]={
	0.0,			//[0]
    1.0,			//[1]
    1.8,			//[2]
    2.4,			//[3]
    4.9				//[4]
};

static const double Ndp[]={
    0.0,			//[0]
	-7.4411,		//[1]
    1.9883,			//[2]
    -2.4925,		//[3]
    -3.2633			//[4]
};

static const double tdp[]={
	0.0,			//[0]
    1.0,			//[1]
    1.6,			//[2]
    2.4,			//[3]
    5.0				//[4]
};

static const double R=0.114547443; //8.314472/72.5854;
static const double M=72.5824; //[g/mol]
static const double Tm=344.494; //[K]
static const double pm=4901.2; //[MPa--> kPa]
static const double pc=4810; //[MPa--> kPa] From (Calm 2007 HPAC Engineering)
static const double rhom=459.0300696; //6.324*M; //[mol/dm^3--> kg/m^3]


// Local function prototypes
static double Pressure_Trho(double T, double rho);
static double IntEnergy_Trho(double T, double rho);
static double Enthalpy_Trho(double T, double rho);
static double Entropy_Trho(double T, double rho);
static double SpecHeatV_Trho(double T, double rho);
static double SpecHeatP_Trho(double T, double rho);
static double Viscosity_Trho(double T, double rho);
static double Conductivity_Trho(double T, double rho);


static double get_Delta(double T, double p);
static double LookupValue(char *Prop,double T, double p);
static double powInt(double x, int y);
static double QuadInterp(double x0, double x1, double x2, double f0, double f1, double f2, double x);

static double a0(double tau, double delta);
static double da0_dtau(double tau, double delta);
static double d2a0_dtau2(double tau, double delta);
static double da0_ddelta(double tau, double delta);

static double dhdT(double tau, double delta);
static double dhdrho(double tau, double delta);
static double dpdT(double tau, double delta);
static double dpdrho(double tau, double delta);

static double ar(double tau, double delta);
static double dar_dtau(double tau,double delta);
static double d2ar_dtau2(double tau, double delta);
static double dar_ddelta(double tau,double delta);
static double d2ar_ddelta2(double tau,double delta);
static double d2ar_ddelta_dtau(double tau,double delta);


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
}
static void BuildLookup(void)
{
	int i,j;
	
	// Properties evaluated at all points with X in the 
	// following p-h plot:
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
				if (Tvec[i]>Tm || pvec[j]>p_dp_R410A(Tvec[i]) || pvec[j]<p_bp_R410A(Tvec[i]))
				{					
					rhomat[i][j]=get_Delta(Tvec[i],pvec[j])*rhom;
					hmat[i][j]=h_R410A(Tvec[i],rhomat[i][j],TYPE_Trho);
					smat[i][j]=s_R410A(Tvec[i],rhomat[i][j],TYPE_Trho);
					umat[i][j]=u_R410A(Tvec[i],rhomat[i][j],TYPE_Trho);
					cpmat[i][j]=cp_R410A(Tvec[i],rhomat[i][j],TYPE_Trho);
					cvmat[i][j]=cv_R410A(Tvec[i],rhomat[i][j],TYPE_Trho);
					viscmat[i][j]=visc_R410A(Tvec[i],rhomat[i][j],TYPE_Trho);
					kmat[i][j]=k_R410A(Tvec[i],rhomat[i][j],TYPE_Trho);
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
		//WriteLookup();
	}
}
double rho_R410A(double T, double p, int Types)
{
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_TPNoLookup:
			return get_Delta(T,p)*rhom;
		case TYPE_TP:
			BuildLookup();
			return LookupValue("rho",T,p);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double p_R410A(double T, double rho)
{
	return Pressure_Trho(T,rho);
}
double h_R410A(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return Enthalpy_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhom;
			return Enthalpy_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("h",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double s_R410A(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return Entropy_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhom;
			return Entropy_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("s",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double u_R410A(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return IntEnergy_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhom;
			return IntEnergy_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("u",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double cp_R410A(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return SpecHeatP_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhom;
			return SpecHeatP_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("cp",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double cv_R410A(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return SpecHeatV_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhom;
			return SpecHeatV_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("cv",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double visc_R410A(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return Viscosity_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhom;
			return Viscosity_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("visc",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}
double k_R410A(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			return Conductivity_Trho(T,p_rho);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhom;
			return Conductivity_Trho(T,rho);
		case TYPE_TP:
			BuildLookup();
			return LookupValue("k",T,p_rho);
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}

double w_R410A(double T, double p_rho, int Types)
{
	double rho;
	double delta,tau,c1,c2;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho; break;
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhom; break;
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
	
	delta=rho/rhom;
	tau=Tm/T;

	c1=-SpecHeatV_Trho(T,rho)/R;
	c2=(1.0+2.0*delta*dar_ddelta(tau,delta)+delta*delta*d2ar_ddelta2(tau,delta));
    return sqrt(-c2*T*SpecHeatP_Trho(T,rho)*1000/c1);
}

double pcrit_R410A(void)
{
	return pm;
}

double Tcrit_R410A(void)
{
	return Tm;
}

double MM_R410A(void)
{
	return M;
}

int errCode_R410A(void)
{
	return errCode;
}

double p_bp_R410A(double T)
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
    
double p_dp_R410A(double T)
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

double rhosatV_R410A(double T)
{
	//int counter;
	//double rho1,rho2,rho3,r1,r2,r3,change,eps=1e-8;
	///* Solve for the density which gives the 
	//pressure equal to the dew pressure */
	//counter=1;
	//change=999;
 //   rho1=10;
 //   rho2=10+.001;
 //   r1=p_dp_R410A(T)-Pressure_Trho(T,rho1);
 //   r2=p_dp_R410A(T)-Pressure_Trho(T,rho2);
 //   while(counter==1 || fabs(change)>eps)
 //   {
 //       rho3=rho2-r2/(r2-r1)*(rho2-rho1);
 //       r3=p_dp_R410A(T)-Pressure_Trho(T,rho3);
 //       change=r2/(r2-r1)*(rho2-rho1);
 //       rho1=rho2; rho2=rho3;
 //       r1=r2; r2=r3;
 //       counter=counter+1;
 //   }
 //   return rho3;

	double THETA,a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8;
	THETA=1-T/344.5;

	a1 =       -4.02;
	a2 =       -18.8;
	a3 =      -18.16; 
	a4 =      -17.37;
	a5 =      -16.15;
	a6 =       -7.87;
	a7 =      -20.46;
	a8 =      -18.59;
	b1 =      0.4803; 
	b2 =       9.217;  
	b3 =        7.52;  
	b4 =       3.471; 
	b5 =       5.124; 
	b6 =       1.498; 
	b7 =       9.831; 
	b8 =        8.99; 

	return exp(a1*pow(THETA,b1)+a2*pow(THETA,b2)+a3*pow(THETA,b3)+a4*pow(THETA,b4)+a5*pow(THETA,b5)+a6*pow(THETA,b6)+a7*pow(THETA,b7)+a8*pow(THETA,b8))*459.53;
}

double rhosatL_R410A(double T)
{
	//int counter;

	//double rho1,rho2,rho3,r1,r2,r3,change,eps=1e-10;
	///* Solve for the density which gives the 
	//pressure equal to the bubble pressure */
	//counter=1;
	//change=999;
 //   rho1=-2.57878924896E-10*powInt(T,6) + 3.58401823351E-07*powInt(T,5) - 2.04856639109E-04*powInt(T,4) + 6.15631400205E-02*powInt(T,3) - 1.02527033660E+01*powInt(T,2) + 8.94030144564E+02*powInt(T,1) - 3.02025468937E+04;
 //   rho2=rho1+0.1;
 //   r1=p_bp_R410A(T)-Pressure_Trho(T,rho1);
 //   r2=p_bp_R410A(T)-Pressure_Trho(T,rho2);
 //   while(counter==1 || fabs(r1)>eps)
 //   {
 //       rho3=rho2-r2/(r2-r1)*(rho2-rho1);
 //       r3=p_bp_R410A(T)-Pressure_Trho(T,rho3);
 //       change=r2/(r2-r1)*(rho2-rho1);
 //       rho1=rho2; rho2=rho3;
 //       r1=r2; r2=r3;
 //       counter=counter+1;
 //   }
 //   return rho3;

	double THETA,a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6;
	THETA=1-T/344.5;

	a1 =       2.876;
	a2 =      -2.465;
	a3 =     -0.9235;
	a4 =     -0.4798;
	a5 =      0.7394;
	a6 =       1.762; 
	b1 =       0.365;
	b2 =      0.6299;
	b3 =       1.864; 
	b4 =       1.932; 
	b5 =       2.559; 
	b6 =       1.149; 

	return exp(a1*pow(THETA,b1)+a2*pow(THETA,b2)+a3*pow(THETA,b3)+a4*pow(THETA,b4)+a5*pow(THETA,b5)+a6*pow(THETA,b6))*459.53;

}

double Viscosity_Trho(double T, double rho)
{
    // Properties taken from "Viscosity of Mixed 
	// Refrigerants R404A,R407C,R410A, and R507A" 
	// by Vladimir Geller, 
	// 2000 Purdue Refrigeration conferences

    // inputs in T [K], and p [kPa]
    // output in Pa-s

   double eta_microPa_s;

   //Set constants required
   double a_0=-2.695e0,a_1=5.850e-2,a_2=-2.129e-5,b_1=9.047e-3,b_2=5.784e-5,
	   b_3=1.309e-7,b_4=-2.422e-10,b_5=9.424e-14,b_6=3.933e-17;

   eta_microPa_s=a_0+a_1*T+a_2*T*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho+b_5*rho*rho*rho*rho*rho+b_6*rho*rho*rho*rho*rho*rho;
   return eta_microPa_s/1e6;
}

double Conductivity_Trho(double T, double rho)
{
	// Properties taken from "Thermal Conductivity 
	// of the Refrigerant mixtures R404A,R407C,R410A, and R507A" 
	// by V.Z. Geller, B.Z. Nemzer, and U.V. Cheremnykh 
	// Int. J. Thermophysics, v. 22, n 4 2001

	// inputs in T [K], and p [kPa] or rho [kg/m^3]
	// output in W/m-K

	//Set constants required
	double a_0=-8.872e0,a_1=7.410e-2,b_1=3.576e-2,b_2=-9.045e-6,b_3=4.343e-8,b_4=-3.705e-12;

	return (a_0+a_1*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho)/1.e6; // from mW/m-K to kW/m-K
}


double dhdT_R410A(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho;
			return dhdT(Tm/T,rho/rhom);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhom;
			return dhdT(Tm/T,rho/rhom);
		case 99:
			rho=get_Delta(T,p_rho)*rhom;
			return (Enthalpy_Trho(T+0.001,rho)-Enthalpy_Trho(T,rho))/0.001;
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}

double dhdrho_R410A(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho;
			return dhdrho(Tm/T,rho/rhom);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhom;
			return dhdrho(Tm/T,rho/rhom);
		case 99:
			rho=get_Delta(T,p_rho)*rhom;
			return (Enthalpy_Trho(T,rho+0.001)-Enthalpy_Trho(T,rho))/0.001;
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}

double dpdT_R410A(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho;
			return dpdT(Tm/T,rho/rhom);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhom;
			return dpdT(Tm/T,rho/rhom);
		case 99:
			rho=get_Delta(T,p_rho)*rhom;
			return (Pressure_Trho(T+0.01,rho)-Pressure_Trho(T,rho))/0.01;
		default:
			errCode=BAD_PROPCODE;
			return _HUGE;
	}
}

double dpdrho_R410A(double T, double p_rho, int Types)
{
	double rho;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho;
			return dpdrho(Tm/T,rho/rhom);
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhom;
			return dpdrho(Tm/T,rho/rhom);
		case 99:
			rho=get_Delta(T,p_rho)*rhom;
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
	T=Tm/tau;  
	//Note: dphi02_dDelta_dTau(tau,delta) is equal to zero
	return R*T/rhom*(tau*(d2ar_ddelta_dtau(tau,delta))+dar_ddelta(tau,delta)+delta*d2ar_ddelta2(tau,delta));
}
static double dhdT(double tau, double delta)
{
	double dhdT_rho,T,dhdtau;
	T=Tm/tau;   
	dhdT_rho=R*tau*(da0_dtau(tau,delta)+dar_dtau(tau,delta))+R*delta*dar_ddelta(tau,delta)+R;
	dhdtau=R*T*(da0_dtau(tau,delta)+ dar_dtau(tau,delta))+R*T*tau*(d2a0_dtau2(tau,delta)+d2ar_dtau2(tau,delta))+R*T*delta*d2ar_ddelta_dtau(tau,delta);
	return dhdT_rho+dhdtau*(-Tm/T/T);
}
static double dpdT(double tau, double delta)
{
	double T,rho;
	T=Tm/tau;     rho=delta*rhom;
	return rho*R*(1+delta*dar_ddelta(tau,delta)-delta*tau*d2ar_ddelta_dtau(tau,delta));
}
static double dpdrho(double tau, double delta)
{
	double T,rho;
	T=Tm/tau;   rho=delta*rhom;
	return R*T*(1+2*delta*dar_ddelta(tau,delta)+delta*delta*d2ar_ddelta2(tau,delta));
}




static double get_Delta(double T, double p)
{
    double change,eps=1e-8, tau,delta_guess;
    int counter=1;
    double r1,r2,r3,delta1,delta2,delta3;
    
	if (T>Tm)
	{
		 delta_guess=p/(8.314/M*T)/rhom/0.7; //0.7 for compressibility factor
	}
	else
	{
		if (p<=(0.5*p_dp_R410A(T)+0.5*p_bp_R410A(T)))
		{
			// Superheated vapor
			delta_guess=p/(8.314/M*T)/rhom;
		}
		else
		{
			// Subcooled liquid
			delta_guess=10;
		}
	}


	tau=Tm/T;
    delta1=delta_guess;
    delta2=delta_guess+.001;
    r1=p/(delta1*rhom*R*T)-1.0-delta1*dar_ddelta(tau,delta1);
    r2=p/(delta2*rhom*R*T)-1.0-delta2*dar_ddelta(tau,delta2);
    while(counter==1 || fabs(change)>eps)
    {
        delta3=delta2-r2/(r2-r1)*(delta2-delta1);
        r3=p/(delta3*rhom*R*T)-1.0-delta3*dar_ddelta(tau,delta3);
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

// ******************************************
//          THERMODYNAMIC FUNCTIONS
// ******************************************

static double Pressure_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhom;
	tau=Tm/T;

	return R*T*rho*(1.0+delta*dar_ddelta(tau,delta));
}
static double IntEnergy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhom;
	tau=Tm/T;
	
	return R*T*tau*(da0_dtau(tau,delta)+dar_dtau(tau,delta));
}
static double Enthalpy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhom;
	tau=Tm/T;
	
	return R*T*(1.0+tau*(da0_dtau(tau,delta)+dar_dtau(tau,delta))+delta*dar_ddelta(tau,delta));
}
static double Entropy_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhom;
	tau=Tm/T;
	
	return R*(tau*(da0_dtau(tau,delta)+dar_dtau(tau,delta))-a0(tau,delta)-ar(tau,delta));
}
static double SpecHeatV_Trho(double T, double rho)
{
	double delta,tau;
	delta=rho/rhom;
	tau=Tm/T;
	
	return -R*tau*tau*(d2a0_dtau2(tau,delta)+d2ar_dtau2(tau,delta));
}
static double SpecHeatP_Trho(double T, double rho)
{
	double delta,tau,c1,c2;
	delta=rho/rhom;
	tau=Tm/T;

	c1=1.0+delta*dar_ddelta(tau,delta)-delta*tau*d2ar_ddelta_dtau(tau,delta);
	c2=1.0+2.0*delta*dar_ddelta(tau,delta)+delta*delta*d2ar_ddelta2(tau,delta);
	
	return R*(-tau*tau*(d2a0_dtau2(tau,delta)+d2ar_dtau2(tau,delta))+c1*c1/c2);
}

// ***************************************************
//                HELMHOLTZ DERIVATIVES
// ***************************************************

static double a0(double tau, double delta)
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
static double da0_dtau(double tau, double delta)
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

static double d2a0_dtau2(double tau, double delta)
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

static double da0_ddelta(double tau, double delta)
{
	return 1/delta;
}

static double ar(double tau, double delta)
{
	double sum=0;
	int k;
	for  (k=1;k<=21;k++)
    {
        sum+=N[k]*pow(delta,i[k])*pow(tau,j[k])*exp(-g[k]*pow(delta,L[k]));
    }
	return sum;
}
static double dar_dtau(double tau,double delta)
{
	double sum=0;
	int k;

	for  (k=1;k<=21;k++)
    {
        sum+=j[k]*N[k]*pow(delta,i[k])*pow(tau,j[k]-1)*exp(-g[k]*pow(delta,L[k]));
    }
	return sum;
}

static double d2ar_dtau2(double tau, double delta)
{
	double sum=0;
	int k;
	for  (k=1;k<=21;k++)
    {
        sum+=N[k]*pow(delta,i[k])*j[k]*(j[k]-1)*pow(tau,j[k]-2)*exp(-g[k]*pow(delta,L[k]));
    }
	return sum;
}

static double d2ar_ddelta2(double tau,double delta)
{
    double dar2_dDelta2=0;
    int k;
    for  (k=1;k<=21;k++)
    {
        dar2_dDelta2=dar2_dDelta2+N[k]*pow(tau,j[k])*exp(-g[k]*pow(delta,L[k]))*(pow(delta,i[k]-2)*pow(i[k],2)-pow(delta,i[k]-2)*i[k]-2*pow(delta,i[k]-2+L[k])*i[k]*g[k]*L[k]-pow(delta,i[k]-2+L[k])*g[k]*pow(L[k],2)+pow(delta,i[k]-2+L[k])*g[k]*L[k]+pow(delta,i[k]+2*L[k]-2)*pow(g[k],2)*pow(L[k],2));
    }
    return dar2_dDelta2;
}

static double dar_ddelta(double tau,double delta)
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

static double d2ar_ddelta_dtau(double tau,double delta)
{
    double dar_dDelta_dTau=0;
    int k;
    for  (k=1;k<=21;k++)
    {
        dar_dDelta_dTau=dar_dDelta_dTau+N[k]*j[k]*pow(tau,j[k]-1)*(i[k]*pow(delta,i[k]-1)*exp(-g[k]*pow(delta,L[k]))+pow(delta,i[k])*exp(-g[k]*pow(delta,L[k]))*(-g[k]*L[k]*pow(delta,L[k]-1)));
    }
    return dar_dDelta_dTau;
}

// ********************************************************
//                  MAINTENANCE FUNCTIONS
// ********************************************************

static double LookupValue(char *Prop, double T, double p)
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

static double powInt(double x, int y)
{
	// Raise a double to an integer power
	// Overload not provided in math.h
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
	/* Quadratic interpolation.  
	Based on method from Kreyszig, 
	Advanced Engineering Mathematics, 9th Edition 
	*/
    double L0, L1, L2;
    L0=((x-x1)*(x-x2))/((x0-x1)*(x0-x2));
    L1=((x-x0)*(x-x2))/((x1-x0)*(x1-x2));
    L2=((x-x0)*(x-x1))/((x2-x0)*(x2-x1));
    return L0*f0+L1*f1+L2*f2;
}
