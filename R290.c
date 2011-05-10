/* Properties of Propane (R290)
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
#include "R290.h"
#include "time.h"

static int errCode;
static char errStr[ERRSTRLENGTH];

#define nP 200
#define nT 200
const static double Tmin=200,Tmax=800, Pmin=500.03,Pmax=16000;

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

static const double Tc=369.89, R_R290=0.188555507, rhoc=220.476, Pc=4251.2, M_R290=44.09562, Ttriple=85.525;
             //           K             kJ/kg-K         kg/m^3     kPa            kg/kmol          K
static const double n[]={0,
0.042910051,
1.7313671,
-2.4516524,
0.34157466,
-0.46047898,
-0.66847295,
0.20889705,
0.19421381,
-0.22917851,
-0.60405866,
0.066680654,
0.017534618,
0.33874242,
0.22228777,
-0.23219062,
-0.09220694,
-0.47575718,
-0.017486824};

static const int d[]={0,
4, //[ 1]
1, //[ 2]
1, //[ 3]
2, //[ 4]
2, //[ 5]
1, //[ 6]
3, //[ 7]
6, //[ 8]
6, //[ 9]
2, //[10]
3, //[11]
1, //[12]
1, //[13]
1, //[14]
2, //[15]
2, //[16]
4, //[17]
1  //[18]
};

static const double t[]={0.00, //offset for natural indices
1.0,
0.33,
0.8,
0.43,
0.9,
2.46,
2.09,
0.88,
1.09,
3.25,
4.62,
0.76,
2.5,
2.75,
3.05,
2.55,
8.4,
6.75};

static const int c[]={
0,0,0,0,0,0, // indices [0-5]
1,
1,
1,
1,
2,
2,
0,0,0,0,0,0,0 // indices [12-18]
};

// alpha instead of eta is used here for consistency with the definitions in R744.c upon which R290.c is based
static const double alpha[]={
0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-11]
0.963,
1.977,
1.917,
2.307,
2.546,
3.28,
14.6};

static const double beta[]={
0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-11]
2.33,
3.47,
3.15,
3.19,
0.92,
18.8,
547.8};

static const double GAMMA[]={
0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-11]
0.684,
0.829,
1.419,
0.817,
1.5,
1.426,
1.093};

static const double epsilon[]={
0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-11]
1.283,
0.6936,
0.788,
0.473,
0.8577,
0.271,
0.948};

//Constants for ideal gas expression
static const double a0[]={0.0,
    -4.970583,
    4.29352,
    3.043,
    5.874,
    9.337,
    7.922
};

static const double b0[]={0.0,
    0,0, //[1 and 2 are not used]
    1.062478,
    3.344237,
    5.363757,
    11.762957
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
                if (Tvec[i]>Tc || pvec[j]<psat_R290(Tvec[i]))
                {					
                    rhomat[i][j]=get_Delta(Tvec[i],pvec[j])*rhoc;
                    hmat[i][j]=h_R290(Tvec[i],rhomat[i][j],TYPE_Trho);
                    smat[i][j]=s_R290(Tvec[i],rhomat[i][j],TYPE_Trho);
                    umat[i][j]=u_R290(Tvec[i],rhomat[i][j],TYPE_Trho);
                    cpmat[i][j]=cp_R290(Tvec[i],rhomat[i][j],TYPE_Trho);
                    cvmat[i][j]=cv_R290(Tvec[i],rhomat[i][j],TYPE_Trho);
                    cmat[i][j]=c_R290(Tvec[i],rhomat[i][j],TYPE_Trho);
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

double rho_R290(double T, double p, int Types)
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

double rhosatL_R290(double T)
{
    const double ti[]={0,0.345,0.74,2.6,7.2};
    const double Ni[]={0,1.82205,0.65802,0.21109,0.083973};
    double summer=1;
    int i;
    double theta;
    theta=1-T/Tc;
    for (i=1;i<=4;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return rhoc*summer;
}

double rhosatV_R290(double T)
{
    const double ti[]={0,0.3785,1.07,2.7,5.5,10,20};
    const double Ni[]={0,-2.4887,-5.1069,-12.174,-30.495,-52.192,-134.89};
    double summer=0,theta;
    int i;
    theta=1.0-T/Tc;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return rhoc*exp(summer);
}

double c_R290(double T, double p_rho, int Types)
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

double MM_R290(void)
{
    return M_R290;
}

double pcrit_R290(void)
{
    return Pc;
}

double Tcrit_R290(void)
{
    return Tc;
}
double Ttriple_R290(void)
{
    return Ttriple;
}
int errCode_R290(void)
{
    return errCode;
}

double psat_R290(double T)
{
    const double ti[]={0,1.0,1.5,2.2,4.8,6.2};
    const double Ni[]={0,-6.7722,1.6938,-1.3341,-3.1876,0.94937};
    double summer=0,theta;
    int i;
    theta=1-T/Tc;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return Pc*exp(Tc/T*summer);
}

double Tsat_R290(double P)
{
    double change,eps=.00005;
    int counter=1;
    double r1,r2,r3,T1,T2,T3;
 
    T1=275;
    T2=275+.01;
    r1=psat_R290(T1)-P;
    r2=psat_R290(T2)-P;
    
    // End at change less than 0.5%
    while(counter==1 || (fabs(change)/fabs(T2)>eps && counter<40))
    {
        T3=T2-0.5*r2/(r2-r1)*(T2-T1);
        r3=psat_R290(T3)-P;
        change=0.5*r2/(r2-r1)*(T2-T1);
        T1=T2;
        T2=T3;
        r1=r2;
        r2=r3;
        counter=counter+1;
    }
    return T3;
}   

double hsat_R290(double T, double x)
{
    double delta,tau;
    
    if (x>0.5)
    {
        delta=rhosatV_R290(T)/rhoc;
        tau=Tc/T;
        return R_R290*T*(1+tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+delta*dphir_dDelta(tau,delta));
    }
    else
    {
        delta=rhosatL_R290(T)/rhoc;
        tau=Tc/T;
        return R_R290*T*(1+tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+delta*dphir_dDelta(tau,delta));
    }   
}

double ssat_R290(double T, double x)
{
    double delta,tau;
    
    if (x>0.5)
    {
        delta=rhosatV_R290(T)/rhoc;
        tau=Tc/T;
        return R_R290*(tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))-phi0(tau,delta)-phir(tau,delta));
    }
    else
    {
        delta=rhosatL_R290(T)/rhoc;
        tau=Tc/T;
        return R_R290*(tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))-phi0(tau,delta)-phir(tau,delta));
    }   
}

double rhosat_R290(double T, double x)
{
    
    if (x>0.5)
        return rhosatV_R290(T);
    else
        return rhosatL_R290(T);
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
        sum+=nv[i]*pow(Tr,tv[i])*pow(rhor,dv[i]);
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

double dhdT_R290(double T, double p_rho, int Types)
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

double dhdrho_R290(double T, double p_rho, int Types)
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

double dpdT_R290(double T, double p_rho, int Types)
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

double dpdrho_R290(double T, double p_rho, int Types)
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
    return R_R290*T*rho*(1.0+delta*dphir_dDelta(tau,delta));
}
static double IntEnergy_Trho(double T, double rho)
{
    double delta,tau;
    delta=rho/rhoc;
    tau=Tc/T;
    
    return R_R290*T*tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta));
}
static double Enthalpy_Trho(double T, double rho)
{
    double delta,tau;
    delta=rho/rhoc;
    tau=Tc/T;
    
    return R_R290*T*(1+tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+delta*dphir_dDelta(tau,delta));
}
static double Entropy_Trho(double T, double rho)
{
    double delta,tau;
    delta=rho/rhoc;
    tau=Tc/T;
    
    return R_R290*(tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))-phi0(tau,delta)-phir(tau,delta));
}
static double SpecHeatV_Trho(double T, double rho)
{
    double delta,tau;
    delta=rho/rhoc;
    tau=Tc/T;
    
    return -R_R290*powI(tau,2)*(dphi02_dTau2(tau,delta)+dphir2_dTau2(tau,delta));
}
static double SpecHeatP_Trho(double T, double rho)
{
    double delta,tau,c1,c2;
    delta=rho/rhoc;
    tau=Tc/T;

    c1=powI(1.0+delta*dphir_dDelta(tau,delta)-delta*tau*dphir2_dDelta_dTau(tau,delta),2);
    c2=(1.0+2.0*delta*dphir_dDelta(tau,delta)+powI(delta,2)*dphir2_dDelta2(tau,delta));
    return R_R290*(-powI(tau,2)*(dphi02_dTau2(tau,delta)+dphir2_dTau2(tau,delta))+c1/c2);
}

static double SpeedSound_Trho(double T, double rho)
{
    double delta,tau,c1,c2;
    delta=rho/rhoc;
    tau=Tc/T;

    c1=-SpecHeatV_Trho(T,rho)/R_R290;
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
    T=Tc/tau;   R=R_R290;
    //Note: dphi02_dDelta_dTau(tau,delta) is equal to zero
    return R*T/rhoc*(tau*(dphir2_dDelta_dTau(tau,delta))+dphir_dDelta(tau,delta)+delta*dphir2_dDelta2(tau,delta));
}
static double dhdT(double tau, double delta)
{
    double dhdT_rho,T,R,dhdtau;
    T=Tc/tau;   R=R_R290;
    dhdT_rho=R*tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+R*delta*dphir_dDelta(tau,delta)+R;
    dhdtau=R*T*(dphi0_dTau(tau,delta)+ dphir_dTau(tau,delta))+R*T*tau*(dphi02_dTau2(tau,delta)+dphir2_dTau2(tau,delta))+R*T*delta*dphir2_dDelta_dTau(tau,delta);
    return dhdT_rho+dhdtau*(-Tc/T/T);
}
static double dpdT(double tau, double delta)
{
    double T,R,rho;
    T=Tc/tau;   R=R_R290;  rho=delta*rhoc;
    return rho*R*(1+delta*dphir_dDelta(tau,delta)-delta*tau*dphir2_dDelta_dTau(tau,delta));
}
static double dpdrho(double tau, double delta)
{
    double T,R,rho;
    T=Tc/tau;   R=R_R290;  rho=delta*rhoc;
    return R*T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*dphir2_dDelta2(tau,delta));
}

/**************************************************/
/*          Private Property Functions            */
/**************************************************/

static double phir(double tau, double delta)
{ 
    
    int i;
    double phir=0,psi;
    
    for (i=1;i<=5;i++)
    {
        phir=phir+n[i]*powI(delta,d[i])*pow(tau,t[i]);
    }
    
    for (i=6;i<=11;i++)
    {
        phir=phir+n[i]*powI(delta,d[i])*pow(tau,t[i])*exp(-powI(delta,c[i]));
    }
    
    for (i=12;i<=18;i++)
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
    for (i=1;i<=5;i++)
    {
        di=(double)d[i];
        dphir_dDelta=dphir_dDelta+n[i]*di*powI(delta,d[i]-1)*pow(tau,t[i]);
    }
    for (i=6;i<=11;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir_dDelta=dphir_dDelta+n[i]*exp(-powI(delta,c[i]))*(powI(delta,d[i]-1)*pow(tau,t[i])*(di-ci*powI(delta,c[i])));
    }
    for (i=12;i<=18;i++)
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
    for (i=1;i<=5;i++)
    {
        di=(double)d[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*di*(di-1.0)*powI(delta,d[i]-2)*pow(tau,t[i]);
    }
    for (i=6;i<=11;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*exp(-powI(delta,c[i]))*(powI(delta,d[i]-2)*pow(tau,t[i])*( (di-ci*powI(delta,c[i]))*(di-1.0-ci*powI(delta,c[i])) - ci*ci*powI(delta,c[i])));
    }
    for (i=12;i<=18;i++)
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

    for (i=1;i<=5;i++)
    {
        di=(double)d[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau + n[i]*di*t[i]*powI(delta,d[i]-1)*pow(tau,t[i]-1.0);
    }
    for (i=6;i<=11;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau + n[i]*exp(-powI(delta,c[i]))*powI(delta,d[i]-1)*t[i]*pow(tau,t[i]-1.0)*(di-ci*powI(delta,c[i]));
    }
    for (i=12;i<=18;i++)
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
    
    for (i=1;i<=5;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powI(delta,d[i])*pow(tau,t[i]-1.0);
    }
    for (i=6;i<=11;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powI(delta,d[i])*pow(tau,t[i]-1.0)*exp(-powI(delta,c[i]));
    }
    for (i=12;i<=18;i++)
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
    
    for (i=1;i<=5;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powI(delta,d[i])*pow(tau,t[i]-2.0);
    }
    for (i=6;i<=11;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powI(delta,d[i])*pow(tau,t[i]-2.0)*exp(-powI(delta,c[i]));
    }
    for (i=12;i<=18;i++)
    {
        psi=exp(-alpha[i]*powI(delta-epsilon[i],2)-beta[i]*powI(tau-GAMMA[i],2));
        dphir2_dTau2=dphir2_dTau2+n[i]*powI(delta,d[i])*pow(tau,t[i])*psi*(powI(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]),2)-t[i]/powI(tau,2)-2.0*beta[i]);
    }
    return dphir2_dTau2;
}

static double phi0(double tau, double delta)
{
    double phi0=0;
    phi0=log(delta)+3*log(tau)+a0[1]+a0[2]*tau
        +a0[3]*log(1-exp(-b0[3]*tau))
        +a0[4]*log(1-exp(-b0[4]*tau))
        +a0[5]*log(1-exp(-b0[5]*tau))
        +a0[6]*log(1-exp(-b0[6]*tau));
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
    double dphi0_dTau=0;
    dphi0_dTau=3.0/tau+a0[2]
        +a0[3]*b0[3]*(1/(exp(b0[3]*tau)-1))
        +a0[4]*b0[4]*(1/(exp(b0[4]*tau)-1))
        +a0[5]*b0[5]*(1/(exp(b0[5]*tau)-1))
        +a0[6]*b0[6]*(1/(exp(b0[6]*tau)-1));
    return dphi0_dTau;
}

static double dphi02_dTau2(double tau, double delta)
{
    double dphi02_dTau2=0;
    dphi02_dTau2=-3.0/powI(tau,2)
        -a0[3]*b0[3]*b0[3]*exp(b0[3]*tau)/powI(exp(b0[3]*tau)-1.0,2)
        -a0[4]*b0[4]*b0[4]*exp(b0[4]*tau)/powI(exp(b0[4]*tau)-1.0,2)
        -a0[5]*b0[5]*b0[5]*exp(b0[5]*tau)/powI(exp(b0[5]*tau)-1.0,2)
        -a0[6]*b0[6]*b0[6]*exp(b0[6]*tau)/powI(exp(b0[6]*tau)-1.0,2);
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
            delta_guess=P/(R_R290*T)/rhoc;
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
        delta_guess=P/(R_R290*T)/rhoc;
    }
    else if(P<psat_R290(T))
        {
        //Superheated vapor
        delta_guess=(rhosatV_R290(T)*P/psat_R290(T))/rhoc;
        }
    else //(T<Tsat_R290(P))
    {
        delta_guess=(rhosatL_R290(T))/rhoc;
    }
    }
    tau=Tc/T;
    delta1=delta_guess;
    delta2=delta_guess+.00001;
    r1=P/(delta1*rhoc*R_R290*T)-1.0-delta1*dphir_dDelta(tau,delta1);
    r2=P/(delta2*rhoc*R_R290*T)-1.0-delta2*dphir_dDelta(tau,delta2);
    
    // End at change less than 0.05%
    while(counter==1 || (fabs(r2)/delta2>eps && counter<40))
    {
        delta3=delta2-r2/(r2-r1)*(delta2-delta1);
        r3=P/(delta3*rhoc*R_R290*T)-1.0-delta3*dphir_dDelta(tau,delta3);
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





