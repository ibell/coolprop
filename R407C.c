/*
Properties for R407C.  
by Ian Bell (ihb2@cornell.edu)

Pseudo-pure fluid thermo props from 
"Pseudo-pure fluid Equations of State for the Refrigerant Blends R410A, R404A, R507C and R407C" 
by E.W. Lemmon, Int. J. Thermophys. v. 24, n4, 2003

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
#include "R407C.h"

static int errCode;
static char errStr[ERRSTRLENGTH];

#define nP 200
#define nT 200

static double Tmin=220, Tmax=450, Pmin=70.03, Pmax=5000;
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
	0.76575,			//[0]
    1.4245,			//[1]
    3.9419,			//[2]
	3.1209			//[3]
};

static const double e[]={
    0.4,			//[0]
    864,			//[1]
    1887,			//[2]
    4802.0			//[3]
};

static const double a[]={
    2.13194,		//[0]
    8.05008,		//[1]
    -14.3914,		//[2]
    1.4245,			//[3]
    3.9419,			//[4]
    3.1209			//[5]
};

static const double b[]={
    0,				//[0]
    0,				//[1]
    -0.4,			//[2]
    2.40437,		//[3]
    5.25122,		//[4]
    13.3632			//[5]
};

static const double N[]={
	 0.0,			//[0]
	 1.0588,		//[1]
	-1.12018,		//[2]
	 0.629064,		//[3]
	-0.351953,		//[4]
	 0.00455978,	//[5]
	-1.75725,		//[6]
	-1.12009,		//[7]
	 0.0277353,		//[8]
	 0.898881,		//[9]
	-1.17591,		//[10]
	 0.0818591,		//[11]
	-0.0794097,		//[12]
	-0.0000104047,	//[13]
	 0.233779,		//[14]
	-0.291790,		//[15]
	 0.0154776,		//[16]
	-0.0314579,		//[17]
	-0.00442552,	//[18]
	-0.0101254,		//[19]
	 0.00915953,	//[20]
	-0.003615		//[21]
};

static const double j[]={
	0.0,			//[0]
    0.241,			//[1]
    0.69,			//[2]
    2.58,			//[3]
    1.15,			//[4]
    0.248,			//[5]
    2.15,			//[6]
    2.43,			//[7]
    5.3,			//[8]
    0.76,			//[9]
    1.48,			//[10]
    0.24,			//[11]
    2.86,			//[12]
    8.0,			//[13]
    3.3,			//[14]
    4.7,			//[15]
    0.45,			//[16]
    8.4,			//[17]
    16.2,			//[18]
    26.0,			//[19]
    16.0,			//[20]
    8.7				//[21]
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
    2,				//[8]
    3,				//[9]
    3,				//[10]
    5,				//[11]
    5,				//[12]
    5,				//[13]
    1,				//[14]
    1,				//[15]
    4,				//[16]
    4,				//[17]
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
    1,				//[12]
    1,				//[13]
    2,				//[14]
    2,				//[15]
    2,				//[16]
    2,				//[17]
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
	0.48722,		//[1]
    -6.6959,		//[2]
    -1.4165,		//[3]
    -2.5109			//[4]
};
    
static const double tbp[]={
	0.0,			//[0]
    0.54,			//[1]
    0.925,			//[2]
    2.7,			//[3]
    4.7				//[4]
};

static const double Ndp[]={
    0.0,			//[0]
	-0.086077,		//[1]
    -6.6364,		//[2]
    -2.4648,		//[3]
    -3.4776			//[4]
};

static const double tdp[]={
	0.0,			//[0]
    0.4,			//[1]
    0.965,			//[2]
    3.1,			//[3]
    5.0				//[4]
};

static const double R=0.096451; //8.314472/86.2036;
static const double M=86.2036; //[g/mol]
static const double Tm=359.345; //[K]
static const double pm=4631.7; //[MPa--> kPa]
static const double pc=4600; //[MPa--> kPa] From (Calm 2007 HPAC Engineering)
static const double rhom=453.43094; //5.260*M; //[mol/dm^3--> kg/m^3]


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
				if (Tvec[i]>Tm || pvec[j]>p_dp_R407C(Tvec[i]) || pvec[j]<p_bp_R407C(Tvec[i]))
				{					
					rhomat[i][j]=get_Delta(Tvec[i],pvec[j])*rhom;
					hmat[i][j]=h_R407C(Tvec[i],rhomat[i][j],TYPE_Trho);
					smat[i][j]=s_R407C(Tvec[i],rhomat[i][j],TYPE_Trho);
					umat[i][j]=u_R407C(Tvec[i],rhomat[i][j],TYPE_Trho);
					cpmat[i][j]=cp_R407C(Tvec[i],rhomat[i][j],TYPE_Trho);
					cvmat[i][j]=cv_R407C(Tvec[i],rhomat[i][j],TYPE_Trho);
					viscmat[i][j]=visc_R407C(Tvec[i],rhomat[i][j],TYPE_Trho);
					kmat[i][j]=k_R407C(Tvec[i],rhomat[i][j],TYPE_Trho);
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
		WriteLookup();
	}
}
double rho_R407C(double T, double p, int Types)
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
double p_R407C(double T, double rho)
{
	return Pressure_Trho(T,rho);
}
double h_R407C(double T, double p_rho, int Types)
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
double s_R407C(double T, double p_rho, int Types)
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
double u_R407C(double T, double p_rho, int Types)
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
double cp_R407C(double T, double p_rho, int Types)
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
double cv_R407C(double T, double p_rho, int Types)
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
double visc_R407C(double T, double p_rho, int Types)
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
double k_R407C(double T, double p_rho, int Types)
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

double w_R407C(double T, double p_rho, int Types)
{
	double rho;
	double delta,tau,c1,c2;
	errCode=0; // Reset error code
	switch(Types)
	{
		case TYPE_Trho:
			rho=p_rho;
		case TYPE_TPNoLookup:
			rho=get_Delta(T,p_rho)*rhom;
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

double pcrit_R407C(void)
{
	return pm;
}

double Tcrit_R407C(void)
{
	return Tm;
}

double MM_R407C(void)
{
	return M;
}

int errCode_R407C(void)
{
	return errCode;
}

double p_bp_R407C(double T)
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
    
double p_dp_R407C(double T)
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

double rhosatV_R407C(double T)
{
	double theta;
	theta=1-T/359.345;
	return exp(+2.481666724+2.974063949110e+00*pow(theta,-0.023725)-14.334177672*theta+60.7010706462*theta*theta-376.514657192*powInt(theta,3)+1178.57631747*powInt(theta,4)-1999.37929072*powInt(theta,5)+1307.74729667*powInt(theta,6));
}

double rhosatL_R407C(double T)
{
	double theta;
	theta=1-T/359.345;
	return exp(5.544589249+1.764403419125e+00*pow(theta,0.12337)+0.544950396285*theta-0.784102758738*theta*theta+0.741332715649*theta*theta*theta);
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
   double a_0=-1.507e0,a_1=4.894e-2,a_2=-9.305e-6,b_1=-3.038e-3,b_2=2.927e-4,
	   b_3=-9.559e-7,b_4=1.739e-9,b_5=-1.455e-12,b_6=4.756e-16;

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
	double a_0=-9.628e0,a_1=7.638e-2,b_1=2.715e-2,b_2=4.963e-5,b_3=-4.912e-8,b_4=2.884e-11;

	return (a_0+a_1*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho)/1.e6; // from mW/m-K to kW/m-K
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
		if (p<=(0.5*p_dp_R407C(T)+0.5*p_bp_R407C(T)))
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
