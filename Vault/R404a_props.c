//You can include any C libraries that you normally use
#if defined(__WIN32__) || defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#endif
#include "math.h"
#include "stdio.h"
#include "R404a_props.h"

#if !defined(_HUGE)
#define _HUGE HUGE
#endif

#define nP 200
#define nT 200

double hmat[nT][nP];
double rhomat[nT][nP];
double cpmat[nT][nP];
double smat[nT][nP];
double cvmat[nT][nP];
double umat[nT][nP];
double viscmat[nT][nP];
double T[nT];
double p[nP];

int TablesBuilt;

static int const1=47;
static double c[4], e[4], a[6], b[6], N[23], j[23], g[23], Nbp[5], Ndp[5], tbp[5], tdp[5];
static int i[23],L[23];
static double Tm, pm, rhom, M, rho,R;
static double cv,cp,h,s;
    
// Private function prototypes
double dar_dDelta_R404a(double tau,double delta);
double dar_dDelta_dTau_R404a(double tau,double delta);
double dar2_dDelta2_R404a(double tau,double delta);
double ar_R404a(double tau,double delta);
double dar_dTau_R404a(double tau,double delta);
double dar2_dTau2_R404a(double tau,double delta);
double a0_R404a(double tau,double delta);
double da0_dTau_R404a(double tau,double delta);
double da0_dDelta_R404a(double tau,double delta);
double da02_dTau2_R404a(double tau,double delta);
double p_bp(double T);
double p_dp(double T);

void WriteLookup()
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
		fprintf(fp_h,",%0.12f",p[j]);
		fprintf(fp_s,",%0.12f",p[j]);
		fprintf(fp_rho,",%0.12f",p[j]);
		fprintf(fp_u,",%0.12f",p[j]);
		fprintf(fp_cp,",%0.12f",p[j]);
		fprintf(fp_cv,",%0.12f",p[j]);
		fprintf(fp_visc,",%0.12f",p[j]);
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
		fprintf(fp_h,"%0.12f",T[i]);
		fprintf(fp_s,"%0.12f",T[i]);
		fprintf(fp_rho,"%0.12f",T[i]);
		fprintf(fp_u,"%0.12f",T[i]);
		fprintf(fp_cp,"%0.12f",T[i]);
		fprintf(fp_cv,"%0.12f",T[i]);
		fprintf(fp_visc,"%0.12f",T[i]);
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
void BuildLookup()
{
	int i,j;
	double rho;
	double Tmin=220,Tmax=340;
	double Pmin=70.03,Pmax=3327;

	if (!TablesBuilt)
	{
		printf_plus("Building Lookup Tables... Please wait...\n");

		for (i=0;i<nT;i++)
		{
			T[i]=Tmin+i*(Tmax-Tmin)/(nT-1);
		}
		for (j=0;j<nP;j++)
		{
			p[j]=Pmin+j*(Pmax-Pmin)/(nP-1);
		}
		for (i=0;i<nT;i++)
		{
			for (j=0;j<nP;j++)
			{
				if (p[j]>p_dp(T[i]))
				{
					rho=get_Delta_R404a(T[i],p[j])*rhom;
					rhomat[i][j]=rho;
					hmat[i][j]=h_R404a(T[i],rhomat[i][j],2);
					smat[i][j]=s_R404a(T[i],rhomat[i][j],2);
					umat[i][j]=u_R404a(T[i],rhomat[i][j],2);
					cpmat[i][j]=cp_R404a(T[i],rhomat[i][j],2);
					cvmat[i][j]=cv_R404a(T[i],rhomat[i][j],2);
					viscmat[i][j]=visc_R404a(T[i],rhomat[i][j],2);
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
	}
}

void setCoeffs_R404a()
{
    int k;
    M=97.6038; //[g/mol]
    Tm=345.270; //[K]
    pm=3.7348*1000; //[MPa--> kPa]
    rhom=4.940*M; //[mol/dm^3--> kg/m^3]
    
    R=8.314472/M;

    c[0]=1.2744;
    c[1]=0.63078;
    c[2]=3.5979;
    c[3]=5.0335;

    e[0]=0.3;
    e[1]=413.0;
    e[2]=804.0;
    e[3]=1727.0;

    a[0]=7.00407;
    a[1]=7.98695;
    a[2]=-18.8664;
    a[3]=0.63078;
    a[4]=3.6979;
    a[5]=5.0335;

    b[0]=0;
    b[1]=0;
    b[2]=-0.3;
    b[3]=1.19617;
    b[4]=2.32861;
    b[5]=5.00188;

    N[1]=6.10984;
    N[2]=-7.79453;
    N[3]=0.0183377;
    N[4]=0.262270;
    N[5]=-0.00351688;
    N[6]=0.0116181;
    N[7]=0.00105992;
    N[8]=0.850922;
    N[9]=-0.520084;
    N[10]=-0.0464225;
    N[11]=0.621190;
    N[12]=-0.195505;
    N[13]=0.336159;
    N[14]=-0.0376062;
    N[15]=-0.00636579;
    N[16]=-0.0758262;
    N[17]=-0.0221041;
    N[18]=0.0310441;
    N[19]=0.0132798;
    N[20]=0.0689437;
    N[21]=-0.0507525;
    N[22]=0.0161382;

    j[1]=0.67;
    j[2]=0.91;
    j[3]=5.96;
    j[4]=0.7;
    j[5]=6.0;
    j[6]=0.3;
    j[7]=0.7;
    j[8]=1.7;
    j[9]=3.3;
    j[10]=7.0;
    j[11]=2.05;
    j[12]=4.3;
    j[13]=2.7;
    j[14]=1.8;
    j[15]=1.25;
    j[16]=12.0;
    j[17]=6.0;
    j[18]=8.7;
    j[19]=11.6;
    j[20]=13.0;
    j[21]=17.0;
    j[22]=16.0;
    

    i[1]=1;
    i[2]=1;
    i[3]=1;
    i[4]=2;
    i[5]=2;
    i[6]=4;
    i[7]=6;
    i[8]=1;
    i[9]=1;
    i[10]=1;
    i[11]=2;
    i[12]=2;
    i[13]=3;
    i[14]=4;
    i[15]=7;
    i[16]=2;
    i[17]=3;
    i[18]=4;
    i[19]=4;
    i[20]=2;
    i[21]=3;
    i[22]=5;

    L[1]=0;
    L[2]=0;
    L[3]=0;
    L[4]=0;
    L[5]=0;
    L[6]=0;
    L[7]=0;
    L[8]=1;
    L[9]=1;
    L[10]=1;
    L[11]=1;
    L[12]=1;
    L[13]=1;
    L[14]=1;
    L[15]=1;
    L[16]=2;
    L[17]=2;
    L[18]=2;
    L[19]=2;
    L[20]=3;
    L[21]=3;
    L[22]=3;
    
    Nbp[1]=0.061067;
    Nbp[2]=-6.5646;
    Nbp[3]=-3.6162;
    Nbp[4]=-3.9771;
    
    tbp[1]=0.54;
    tbp[2]=0.965;
    tbp[3]=3.7;
    tbp[4]=9.0;
    
    Ndp[1]=-0.00026863;
    Ndp[2]=-6.5757;
    Ndp[3]=-4.1802;
    Ndp[4]=-7.9102;

    tdp[1]=0.1;
    tdp[2]=0.972;
    tdp[3]=3.8;
    tdp[4]=9.0;

    for (k=1;k<=7;k++)
    {
    g[k]=0;
    }

    for (k=8;k<=22;k++)
    {
    g[k]=1;
    }
/*
    mexPrintf("set Coeffs\n");
*/
}

double p_bp(double T)
{
    double tau,theta,summer;
    int k;
    setCoeffs_R404a();
    tau=Tm/T;
    theta=1-T/Tm;
    if (T>Tm)
    {
		printf("Temperature above critical temperature (in p_bp)");
        //mexErrMsgTxt("Temperature above critical temperature (in p_bp)");
    }
    // Calculation of bubble point for given temperature
    summer=0;
    for (k=1;k<=4;k++)
    {
        summer=summer+Nbp[k]*pow(theta,tbp[k]);
    }
    return pm*exp(tau*summer);
}
    
double p_dp(double T)
{
    double tau,theta,summer;
    int k;
    setCoeffs_R404a();

    if (T>Tm)
    {
		printf("Temperature above critical temperature (in p_dp)");
    }
    tau=Tm/T;
    theta=1-T/Tm;
    // Calculation of dew point for given temperature
    summer=0;
    for (k=1;k<=4;k++)
    {
        summer=summer+Ndp[k]*pow(theta,tdp[k]);
    }
    return pm*exp(tau*summer);
}


double visc_R404a(double T, double param2, int inputsType)
{
    // Properties taken from "Viscosity of Mixed Refrigerants R404A,R407C,R410A, and R507A" 
	// by Vladimir Geller, 2000 Purdue Refrigeration conferences
    // inputs in T [K], and p [kPa]
    // output in [Pa-s]

   double delta,eta_microPa_s, rho;
   double a0=9.766e-1,a1=3.676e-2,a2=2.938e-6,b1=2.260e-3,b2=1.786e-4,b3=-4.202e-7,b4=8.489e-10,b5=-8.670e-13,b6=3.566e-16;
   setCoeffs_R404a();
   if (inputsType==TYPE_TP)
    {
        delta=get_Delta_R404a(T,param2);
    }
    else
    {
        delta=param2/rhom;
    }
   rho=delta*rhom;

   eta_microPa_s=a0+a1*T+a2*T*T+b1*rho+b2*rho*rho+b3*rho*rho*rho+b4*rho*rho*rho*rho+b5*rho*rho*rho*rho*rho+b6*rho*rho*rho*rho*rho*rho;
   return eta_microPa_s/1e6;
}

double rho_R404a(double T, double p)
{
    double delta;
    setCoeffs_R404a();
    delta=get_Delta_R404a(T,p);
	BuildLookup();
	WriteLookup();
    return delta*rhom;
}

double h_R404a(double T, double param2, int inputsType)
{
    double delta, tau;
    setCoeffs_R404a();
    if (inputsType==1)
    {
        delta=get_Delta_R404a(T,param2);
    }
    else
    {
        delta=param2/rhom;
    }
    tau=Tm/T;
    return R*T*(tau*(da0_dTau_R404a(tau,delta)+dar_dTau_R404a(tau,delta))+delta*dar_dDelta_R404a(tau,delta)+1); 
}

double u_R404a(double T, double param2, int inputsType)
{
    double delta, tau;
    setCoeffs_R404a();
    if (inputsType==1)
    {
        delta=get_Delta_R404a(T,param2);
    }
    else
    {
        delta=param2/rhom;
    }
    tau=Tm/T;
    return R*T*(tau*(da0_dTau_R404a(tau,delta)+dar_dTau_R404a(tau,delta))+delta*dar_dDelta_R404a(tau,delta)+1)-1-delta*dar_dDelta_R404a(tau,delta); 
}

double s_R404a(double T, double param2, int inputsType)
{
    double delta, tau;    
    setCoeffs_R404a();
    if (inputsType==1)
    {
        delta=get_Delta_R404a(T,param2);
    }
    else
    {
        delta=param2/rhom;
    }
    tau=Tm/T;
    return R*(tau*(da0_dTau_R404a(tau,delta)+dar_dTau_R404a(tau,delta))-a0_R404a(tau,delta)-ar_R404a(tau,delta));
}

// double s_R404a_Tr(double T, double rho)
// {
//     double delta, tau;    
//     setCoeffs_R404a();
//     delta=rho/rhom;
//     tau=Tm/T;
//     return R*(tau*(da0_dTau_R404a(tau,delta)+dar_dTau_R404a(tau,delta))-a0_R404a(tau,delta)-ar_R404a(tau,delta));
// }

double cp_R404a(double T, double param2, int inputsType)
{
    double delta, tau,cv;  
    setCoeffs_R404a();
    if (inputsType==1)
    {
        delta=get_Delta_R404a(T,param2);
    }
    else
    {
        delta=param2/rhom;
    }
    tau=Tm/T;
    cv=-R*pow(tau,2.0)*(dar2_dTau2_R404a(tau,delta)+da02_dTau2_R404a(tau,delta));
    return cv+R*pow(1.0+delta*dar_dDelta_R404a(tau,delta)-delta*tau*dar_dDelta_dTau_R404a(tau,delta),2.0)/(1.0+2.0*delta*dar_dDelta_R404a(tau,delta)+pow(delta,2.0)*dar2_dDelta2_R404a(tau,delta));
}

double cv_R404a(double T, double param2, int inputsType)
{
    double delta, tau;    
    setCoeffs_R404a();
    if (inputsType==1)
    {
        delta=get_Delta_R404a(T,param2);
    }
    else
    {
        delta=param2/rhom;
    }
    tau=Tm/T;
    return -R*pow(tau,2.0)*(dar2_dTau2_R404a(tau,delta)+da02_dTau2_R404a(tau,delta));
}

double MM_R404a()
{
	return M;
}

double Psat_R404a(double T, double x)
{
    if (x==0)
    {
        return p_bp(T);
    }
    else
    {
        return p_dp(T);
    }
}



// ***********************************************************************
// ***********************************************************************
// ********************** Residual Derivatives ***************************
// ***********************************************************************
// ***********************************************************************


double a0_R404a(double tau,double delta)
{
     return log(delta)-log(tau)+a[0]+a[1]*tau+a[2]*pow(tau,b[2])+a[3]*log(1.0-exp(-b[3]*tau))+a[4]*log(1.0-exp(-b[4]*tau))+a[5]*log(1.0-exp(-b[5]*tau));
}

double da0_dDelta_R404a(double tau,double delta)
{
    return 1.0/delta;
}

double da0_dTau_R404a(double tau,double delta)
{
    return -1.0/tau+a[1]+a[2]*b[2]*pow(tau,b[2]-1.0)+a[3]/(1.0-exp(-b[3]*tau))*(-exp(-b[3]*tau))*(-b[3])+a[4]/(1.0-exp(-b[4]*tau))*(-exp(-b[4]*tau))*(-b[4])+a[5]/(1.0-exp(-b[5]*tau))*(-exp(-b[5]*tau))*(-b[5]);
}

double da02_dTau2_R404a(double tau,double delta)
{
 
    return pow(tau,-2.0)+a[2]*b[2]*(b[2]-1.0)*pow(tau,b[2]-2.0)-a[3]*pow(b[3],2.0)*exp(-b[3]*tau)/pow(-1.0+exp(-b[3]*tau),2.0)-a[4]*pow(b[4],2.0)*exp(-b[4]*tau)/pow(-1.0+exp(-b[4]*tau),2.0)-a[5]*pow(b[5],2.0)*exp(-b[5]*tau)/pow(-1.0+exp(-b[5]*tau),2.0);
}
    
double ar_R404a(double tau,double delta)
{
    double summer=0; 
    int k;
     for  (k=1;k<=22;k++)
    {
        summer=summer+N[k]*powInt(delta,i[k])*pow(tau,j[k])*exp(-g[k]*powInt(delta,L[k]));
    }
    return summer;
}

double dar_dTau_R404a(double tau,double delta)
{
    double summer=0; 
    int k;
     for  (k=1;k<=22;k++)
    {
        summer=summer+N[k]*powInt(delta,i[k])*j[k]*pow(tau,j[k]-1)*exp(-g[k]*powInt(delta,L[k]));
    }
    return summer;
}

double dar2_dTau2_R404a(double tau,double delta)
{
    double summer=0; 
    int k;
    for (k=1;k<=22;k++)
    {
        summer=summer+N[k]*powInt(delta,i[k])*j[k]*(j[k]-1)*pow(tau,j[k]-2)*exp(-g[k]*powInt(delta,L[k]));
    }
    return summer;
}

double dar2_dDelta2_R404a(double tau,double delta)
{
    double dar2_dDelta2=0;
    int k;
    for  (k=1;k<=22;k++)
    {
        dar2_dDelta2=dar2_dDelta2+N[k]*powInt(delta,i[k])*pow(tau,j[k])*exp(-g[k]*powInt(delta,L[k]))*(powInt(i[k],2)-i[k]-2*i[k]*g[k]*powInt(delta,L[k])*L[k]-g[k]*powInt(delta,L[k])*powInt(L[k],2)+g[k]*powInt(delta,L[k])*L[k]+pow(g[k],2.0)*powInt(powInt(delta,L[k]),2)*powInt((double)L[k],2))/pow(delta,2.0);
    }
    return dar2_dDelta2;
}

double dar_dDelta_R404a(double tau,double delta)
{
    double dar_dDelta=0;
    int k;
    for  (k=1;k<=22;k++)
    {
        dar_dDelta=dar_dDelta+N[k]*powInt(delta,i[k])*pow(tau,j[k])*exp(-g[k]*powInt(delta,L[k]))*((double)i[k]-g[k]*powInt(delta,L[k])*(double)L[k])/delta;
    }
    return dar_dDelta;
}

double dar_dDelta_dTau_R404a(double tau,double delta)
{
    double dar_dDelta_dTau=0;
    int k;
    for  (k=1;k<=22;k++)
    {
        dar_dDelta_dTau=dar_dDelta_dTau+N[k]*j[k]*pow(tau,j[k]-1)*(i[k]*pow(delta,i[k]-1)*exp(-g[k]*pow(delta,L[k]))+pow(delta,i[k])*exp(-g[k]*pow(delta,L[k]))*(-g[k]*L[k]*pow(delta,L[k]-1)));
    }
    return dar_dDelta_dTau;
}

double get_Delta_R404a(double T, double p)
{
    double delta_guess,tau;
    double change,eps=1e-5;
    int counter=1;
    double r1,r2,r3,delta1,delta2,delta3;

    if (p<=(0.5*p_dp(T)+0.5*p_bp(T)))
    {
        delta_guess=p/(8.314/M*T)/rhom;
    }
    else
    {
        delta_guess=10;
    }

    tau=Tm/T;
    delta1=delta_guess;
    delta2=delta_guess+.001;
    r1=p/(delta1*rhom*R*T)-1-delta1*dar_dDelta_R404a(tau,delta1);
    r2=p/(delta2*rhom*R*T)-1-delta2*dar_dDelta_R404a(tau,delta2);
    while(counter==1 || fabs(change)>eps)
    {
        delta3=delta2-r2/(r2-r1)*(delta2-delta1);
        r3=p/(delta3*rhom*R*T)-1-delta3*dar_dDelta_R404a(tau,delta3);
        change=r2/(r2-r1)*(delta2-delta1);
        delta1=delta2;
        delta2=delta3;
        r1=r2;
        r2=r3;
        counter=counter+1;
//         mexPrintf("%g \t %g \t %g \t %g \t %g \t %g \t %g \n",delta1,r1,delta2,r2,change,T,p);
    }

	//printf("delta: %g delta_guess: %g, ideal: %g ---- %d\n",delta3,delta_guess,p/(8.314/97.6038*T)/rhom,counter);
    return delta3;

}

