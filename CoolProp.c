#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#endif

#if defined(_WIN32) || defined(__WIN32__) || defined(_WIN64) || defined(__WIN64__)
#define __ISWINDOWS__
#endif

#if defined(__ISWINDOWS__)
#include <windows.h>
#include "REFPROP.h"
#endif

#include <stdlib.h>
#include "string.h"
#include <stdio.h>
#include "CoolPropTools.h"
#include "CoolProp.h"

// Some constants for REFPROP... defined by macros for ease of use 
#define refpropcharlength 255
#define filepathlength 255
#define lengthofreference 3
#define errormessagelength 255
#define ncmax 20		// Note: ncmax is the max number of components
#define numparams 72 
#define maxcoefs 50

char LoadedFluid[255];
char LoadedREFPROPRef[255];
int FluidType;

double R_u=8.314472; //From Lemmon et al 2009 (propane)

// Function prototypes
double _Dekker_Tsat(double x_min, double x_max, double eps, double p, double Q,char *Ref);
int LoadFluid(char *Ref);

// Global residual Helmholtz derivatives function pointers
double (*p_func)(double,double);
double (*h_func)(double,double,int);
double (*s_func)(double,double,int);
double (*u_func)(double,double,int);
double (*rho_func)(double,double,int);
double (*cp_func)(double,double,int);
double (*cv_func)(double,double,int);
double (*visc_func)(double,double,int);
double (*k_func)(double,double,int);
double (*w_func)(double,double,int);

double (*MM_func)(void);
double (*rhocrit_func)(void);
double (*Tcrit_func)(void);
double (*pcrit_func)(void);
double (*Ttriple_func)(void);

//For the pure fluids
double (*psat_func)(double); 
// For the psedo-pure fluids
double (*pdp_func)(double);
double (*pbp_func)(double);

// Global residual Helmholtz derivatives function pointers
double (*phir_func)(double,double);
double (*dphir_dDelta_func)(double,double);
double (*dphir2_dDelta2_func)(double,double);
double (*dphir2_dDelta_dTau_func)(double,double);
double (*dphir_dTau_func)(double,double);
double (*dphir2_dTau2_func)(double,double);
double (*phi0_func)(double,double);
double (*dphi0_dDelta_func)(double,double);
double (*dphi02_dDelta2_func)(double,double);
double (*dphi0_dTau_func)(double,double);
double (*dphi02_dTau2_func)(double,double);

double (*rhosatV_func)(double);
double (*rhosatL_func)(double);

#define NLUTFLUIDS 30 // How many saturation curve LUTs you can have in one instance of CoolProp
#define NLUT 300 // How many values to use for each lookup table

double Tsat_LUT[NLUTFLUIDS][NLUT],rhosatL_LUT[NLUTFLUIDS][NLUT],rhosatV_LUT[NLUTFLUIDS][NLUT],psat_LUT[NLUTFLUIDS][NLUT];
char RefLUT[NLUTFLUIDS][255]; // The names of the fluids that are loaded
int FlagUseSaturationLUT=0; //Default to not use LUT
int FlagUseSinglePhaseLUT=0; //Default to not use LUT

static struct fluidParamsVals Fluid;
    
static double get_Delta(double T, double p)
{
    double eps=1e-8, tau,delta_guess,M,Tc,rhoc,R;
    int counter=1;
    double r1,r2,r3,delta1,delta2,delta3;
    
    //Molecular mass of the fluid
    M=Fluid.MM;
    if (M<0.001 || M>1000) {}
    rhoc=Fluid.rhoc;
    Tc=Fluid.Tc;
    R=R_u/M;
        
    if (T>(Fluid.Tc) )
    {
        delta_guess=p/(8.314/M*T)/Fluid.rhoc/0.7; //0.7 for compressibility factor
    }
    else
    {
        if (Fluid.Type == FLUIDTYPE_REFRIGERANT_PURE && fabs(p-Fluid.funcs.psat(T))/p>0.0001)
        {
            // Superheated vapor
            delta_guess=p/(8.314/M*T)/Fluid.rhoc;
        }
        else if (Fluid.Type == FLUIDTYPE_REFRIGERANT_PSEUDOPURE &&(p<=0.5*Fluid.funcs.p_dp(T)+0.5*Fluid.funcs.p_bp(T)))
        {
            // Superheated vapor
            delta_guess=p/(8.314/M*T)/Fluid.rhoc;
        }
        else
        {
            // Subcooled liquid
            delta_guess=10;
        }
    }

    tau=Tc/T;
    delta1=delta_guess;
    delta2=delta_guess+.001;
    r1=p/(delta1*rhoc*R*T)-1.0-delta1*dar_ddelta_R410A(tau,delta1);
    r2=p/(delta2*rhoc*R*T)-1.0-delta2*dar_ddelta_R410A(tau,delta2);
    while( counter==1 || fabs(r2)>eps)
    {
        delta3=delta2-r2/(r2-r1)*(delta2-delta1);
        r3=p/(delta3*rhoc*R*T)-1.0-delta3*dar_ddelta_R410A(tau,delta3);
        delta1=delta2;
        delta2=delta3;
        r1=r2;
        r2=r3;
        counter=counter+1;
    }
    return delta3;
}

int BuildSaturationLUT(char *Ref)
{
    // returns the index of the LUT
    int i,k;
    double T_t,T_c,dT;
    
    // Check if the LUT is already in memory
    for (k=0;k<NLUTFLUIDS;k++)
    {
        // If it is found, break out of loop and stop; no LUT needed
        if (!strcmp(RefLUT[k],Ref)) return k;
        // You made it to an empty string, need a LUT, jump out of for loop
        if (!strcmp(RefLUT[k],"")) break;
    }
    if (k==NLUTFLUIDS)
    {
        // Ran out of fluids
        printf("Sorry, ran out of saturation LUT slots, increase NLUTFLUIDS and recompile\n");
    }
    
    // Therefore it hasn't been found yet, assign the refrigerant name
    strcpy(RefLUT[k],Ref);
    
    // Linearly space the values from the triple point of the fluid to the just shy of the critical point
    T_t=Props('R','T',0,'P',0,Ref);
    T_c=Props('B','T',0,'P',0,Ref)-0.000001;
    
    dT=(T_c-T_t)/(NLUT-1);
    for (i=0;i<NLUT;i++)
    {
        // Linearly space the elements
        Tsat_LUT[k][i]=T_t+i*dT;
        // Calculate the saturation properties
        rhosatPure(Ref, Tsat_LUT[k][i], &(rhosatL_LUT[k][i]), &(rhosatV_LUT[k][i]), &(psat_LUT[k][i]));
    }
    return k;
}   

void UseSaturationLUT(int OnOff)
{
    if (OnOff==1 || OnOff==0)
    {
        FlagUseSaturationLUT=OnOff;
    }
    else
    {
        printf("Sorry, UseSaturationLUT() takes an integer input, either 0 (no) or 1 (yes)\n");
    }
}

void UseSinglePhaseLUT(int OnOff)
{
    if (OnOff==1 || OnOff==0)
    {
        FlagUseSinglePhaseLUT=OnOff;
    }
    else
    {
        printf("Sorry, UseSinglePhaseLUT() takes an integer input, either 0 (no) or 1 (yes)\n");
    }
}
int SinglePhaseLUTStatus(void)
{
    return FlagUseSinglePhaseLUT;
}

static double QuadInterpolate(double x0, double x1, double x2, double f0, double f1, double f2, double x)
{
    double L0, L1, L2;
    L0=((x-x1)*(x-x2))/((x0-x1)*(x0-x2));
    L1=((x-x0)*(x-x2))/((x1-x0)*(x1-x2));
    L2=((x-x0)*(x-x1))/((x2-x0)*(x2-x1));
    return L0*f0+L1*f1+L2*f2;
}

static double ApplySaturationLUT(char *OutPropName,char *InPropName,double T_or_p, char *Ref)
{
    int i,k;
    double x0,x1,x2,y0,y1,y2;
    // pointers to the x and y vectors for later
    double (*y)[NLUT],(*x)[NLUT];
    
    // First try to build the LUT, get index of LUT
    k=BuildSaturationLUT(Ref);
    
    // Then find the index to the left of the temp
    for(i=0;i<NLUT;i++)
    {
        if (Tsat_LUT[k][i]>=T_or_p) break;
    }
    
    if (!strcmp(InPropName,"T") || !strcmp(InPropName,"Tsat"))
        x=&(Tsat_LUT[k]);
    
    if (!strcmp(OutPropName,"rhoL"))
        y=&(rhosatL_LUT[k]);
    else if (!strcmp(OutPropName,"rhoV"))
        y=&(rhosatV_LUT[k]);
    else if (!strcmp(OutPropName,"p"))
        y=&(psat_LUT[k]);
    
    // Need a three-point set to interpolate using a quadratic.
    if (i<NLUT-2)
    {
        x0=(*x)[i];
        x1=(*x)[i+1];
        x2=(*x)[i+2];
        y0=(*y)[i];
        y1=(*y)[i+1];
        y2=(*y)[i+2];
    }
    else
    {
        // Go "backwards" with the interpolation range
        x0=(*x)[i+1];
        x1=(*x)[i];
        x2=(*x)[i-1];
        y0=(*y)[i+1];
        y1=(*y)[i];
        y2=(*y)[i-1];
    }
    return QuadInterpolate(x0,x1,x2,y0,y1,y2,T_or_p);
}

void MatInv_2(double A[2][2] , double B[2][2])
{
    double Det;
    //Using Cramer's Rule to solve

    Det=A[0][0]*A[1][1]-A[1][0]*A[0][1];
    B[0][0]=1.0/Det*A[1][1];
    B[1][1]=1.0/Det*A[0][0];
    B[1][0]=-1.0/Det*A[1][0];
    B[0][1]=-1.0/Det*A[0][1];
}

double Pressure_Trho(double T, double rho)
{
    double delta,tau,R;
    R=R_u/Fluid.MM;
    tau=Fluid.Tc/T;
    delta=rho/Fluid.rhoc;
    return R*T*rho*(1.0+delta*dphir_dDelta_func(tau,delta));
}
double IntEnergy_Trho(double T, double rho)
{
    double delta,tau,R;
    R=R_u/Fluid.MM;
    tau=Fluid.Tc/T;
    delta=rho/Fluid.rhoc;
    return R*T*tau*(dphi0_dTau_func(tau,delta)+dphir_dTau_func(tau,delta));
}
double Enthalpy_Trho(double T, double rho)
{
    double delta,tau,R;
    R=R_u/Fluid.MM;
    tau=Fluid.Tc/T;
    delta=rho/Fluid.rhoc;
    return R*T*(1.0+tau*(dphi0_dTau_func(tau,delta)+dphir_dTau_func(tau,delta))+delta*dphir_dDelta_func(tau,delta));
}
double Entropy_Trho(double T, double rho)
{
    double delta,tau,R;
    R=R_u/Fluid.MM;
    tau=Fluid.Tc/T;
    delta=rho/Fluid.rhoc;
    return R*(tau*(dphi0_dTau_func(tau,delta)+dphir_dTau_func(tau,delta))-phi0_func(tau,delta)-phir_func(tau,delta));
}
double SpecHeatV_Trho(double T, double rho)
{
    double delta,tau,R;
    R=R_u/Fluid.MM;
    tau=Fluid.Tc/T;
    delta=rho/Fluid.rhoc;
    
    return -R*powInt(tau,2)*(dphi02_dTau2_func(tau,delta)+dphir2_dTau2_func(tau,delta));
}
double SpecHeatP_Trho(double T, double rho)
{
    double delta,tau,c1,c2,R;
    R=R_u/Fluid.MM;
    tau=Fluid.Tc/T;
    delta=rho/Fluid.rhoc;

    c1=powInt(1.0+delta*dphir_dDelta_func(tau,delta)-delta*tau*dphir2_dDelta_dTau_func(tau,delta),2);
    c2=(1.0+2.0*delta*dphir_dDelta_func(tau,delta)+powInt(delta,2)*dphir2_dDelta2_func(tau,delta));
    return R*(-powInt(tau,2)*(dphi02_dTau2_func(tau,delta)+dphir2_dTau2_func(tau,delta))+c1/c2);
}

double SpeedSound_Trho(double T, double rho)
{
    double delta,tau,c1,c2,R;
    R=R_u/Fluid.MM;
    tau=Fluid.Tc/T;
    delta=rho/Fluid.rhoc;

    c1=-SpecHeatV_Trho(T,rho)/R;
    c2=(1.0+2.0*delta*dphir_dDelta_func(tau,delta)+powInt(delta,2)*dphir2_dDelta2_func(tau,delta));
    return sqrt(-c2*T*SpecHeatP_Trho(T,rho)*1000/c1);
}

double Gibbs_Trho(double T,double rho)
{
    return Enthalpy_Trho(T,rho)-T*Entropy_Trho(T,rho);
}

double Density_Tp(double T, double p, double rho)
{
    double delta,tau,dpdrho,error=999,R;
    R=R_u/Fluid.MM;
    tau=Fluid.Tc/T;
    delta=rho/Fluid.rhoc;
    while (fabs(error)>1e-7)
    {
        delta=rho/Fluid.rhoc;
        // Use Newton's method to find the saturation density since the derivative of pressure w.r.t. density is known from EOS
        dpdrho=R*T*(1+2*delta*dphir_dDelta_func(tau,delta)+delta*delta*dphir2_dDelta2_func(tau,delta));
        // Update the step using Newton's method
        rho=rho-(Pressure_Trho(T,rho)-p)/dpdrho;
        // Residual
        error=Pressure_Trho(T,rho)-p;
    }		
    return rho;
}

void rhosatPure(char *Ref, double T, double *rhoLout, double *rhoVout, double *pout)
{
    // Only works for pure fluids (no blends)
    // At equilibrium, saturated vapor and saturated liquid are at the same pressure and the same Gibbs energy
    double rhoL,rhoV,p,error=999,x1,x2,x3,y1,y2,f,p_guess;
    int iter;

    LoadFluid(Ref);

    // Use the density ancillary function as the starting point for the secant solver
    rhoL=rhosatL_func(T);
    rhoV=rhosatV_func(T);
    p_guess=Pressure_Trho(T,rhoV);
    
    //printf("L %g V %G p %g\n",rhoL,rhoV,p_guess);

    iter=1;
    // Use a secant method to obtain pressure
    while ((iter<=3 || fabs(error)>1e-7) && iter<100)
    {
        if (iter==1){x1=p_guess; p=x1;}
        if (iter==2){x2=1.0001*p_guess; p=x2;}
        if (iter>2) {p=x2;}
            //Recalculate the densities based on the current pressure
            rhoL=Density_Tp(T,p,rhoL);
            rhoV=Density_Tp(T,p,rhoV);
            // Residual between saturated liquid and saturated vapor gibbs function
            f=Gibbs_Trho(T,rhoL)-Gibbs_Trho(T,rhoV);
            if (iter>10)
            {
            	printf("iter>100:: L %g V %G p %g\n",rhoL,rhoV,p_guess);
            }
        if (iter==1){y1=f;}
        if (iter>1)
        {
            y2=f;
            x3=x2-y2/(y2-y1)*(x2-x1);
            error=f;
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
        if (iter>90)
        {
            printf("rhosatPure failed, current values of rhoL and rhoV are %g,%g\n",rhoL,rhoV);
            return;
        }
    }
    *rhoLout=rhoL;
    *rhoVout=rhoV;
    *pout=p;
    return;
}

int LoadFluid(char *Ref)
{
    if (!strcmp(LoadedFluid,Ref))
    {
        // Already Loaded, don't do anything else
        return 0;
    }
    else
    {
        strcpy(LoadedFluid,Ref);
        printf("Loading fluid...");
        // Wire up the function pointers for the given refrigerant
        if (!strcmp(Ref,"Argon"))
        {
            FluidType=FLUIDTYPE_REFRIGERANT_PURE;
            p_func=p_Argon;
            h_func=h_Argon;
            s_func=s_Argon;
            u_func=u_Argon;
            rho_func=rho_Argon;
            cp_func=cp_Argon;
            cv_func=cv_Argon;
            visc_func=visc_Argon;
            k_func=k_Argon;
            w_func=w_Argon;
            Ttriple_func=Ttriple_Argon;
            Tcrit_func=Tcrit_Argon;
            pcrit_func=pcrit_Argon;
            rhocrit_func=rhocrit_Argon;
            MM_func=MM_Argon;
            rhosatV_func=rhosatV_Argon;
            rhosatL_func=rhosatL_Argon;
            psat_func=psat_Argon;

            phir_func=phir_Argon;
            dphir_dDelta_func=dphir_dDelta_Argon;
            dphir2_dDelta2_func=dphir2_dDelta2_Argon;
            dphir2_dDelta_dTau_func=dphir2_dDelta_dTau_Argon;
            dphir_dTau_func=dphir_dTau_Argon;
            dphir2_dTau2_func=dphir2_dTau2_Argon;
            phi0_func=phi0_Argon;
            dphi0_dDelta_func=dphi0_dDelta_Argon;
            dphi02_dDelta2_func=dphi02_dDelta2_Argon;
            dphi0_dTau_func=dphi0_dTau_Argon;
            dphi02_dTau2_func=dphi02_dTau2_Argon;

        }
        else if (!strcmp(Ref,"Nitrogen") || !strcmp(Ref,"N2"))
        {
            FluidType=FLUIDTYPE_REFRIGERANT_PURE;
            p_func=p_Nitrogen;
            h_func=h_Nitrogen;
            s_func=s_Nitrogen;
            u_func=u_Nitrogen;
            rho_func=rho_Nitrogen;
            cp_func=cp_Nitrogen;
            cv_func=cv_Nitrogen;
            visc_func=visc_Nitrogen;
            k_func=k_Nitrogen;
            w_func=w_Nitrogen;
            Ttriple_func=Ttriple_Nitrogen;
            Tcrit_func=Tcrit_Nitrogen;
            pcrit_func=pcrit_Nitrogen;
            rhocrit_func=rhocrit_Nitrogen;
            MM_func=MM_Nitrogen;
            rhosatV_func=rhosatV_Nitrogen;
            rhosatL_func=rhosatL_Nitrogen;
            psat_func=psat_Nitrogen;

            phir_func=phir_Nitrogen;
            dphir_dDelta_func=dphir_dDelta_Nitrogen;
            dphir2_dDelta2_func=dphir2_dDelta2_Nitrogen;
            dphir2_dDelta_dTau_func=dphir2_dDelta_dTau_Nitrogen;
            dphir_dTau_func=dphir_dTau_Nitrogen;
            dphir2_dTau2_func=dphir2_dTau2_Nitrogen;
            phi0_func=phi0_Nitrogen;
            dphi0_dDelta_func=dphi0_dDelta_Nitrogen;
            dphi02_dDelta2_func=dphi02_dDelta2_Nitrogen;
            dphi0_dTau_func=dphi0_dTau_Nitrogen;
            dphi02_dTau2_func=dphi02_dTau2_Nitrogen;
        }
        else if (!strcmp(Ref,"R744") || !strcmp(Ref,"CO2"))
        {
            FluidType=FLUIDTYPE_REFRIGERANT_PURE;
            p_func=p_R744;
            h_func=h_R744;
            s_func=s_R744;
            u_func=u_R744;
            rho_func=rho_R744;
            cp_func=cp_R744;
            cv_func=cv_R744;
            visc_func=visc_R744;
            k_func=k_R744;
            w_func=w_R744;
            Ttriple_func=Ttriple_R744;
            Tcrit_func=Tcrit_R744;
            pcrit_func=pcrit_R744;
            rhocrit_func=rhocrit_R744;
            MM_func=MM_R744;
            rhosatV_func=rhosatV_R744;
            rhosatL_func=rhosatL_R744;
            psat_func=psat_R744;

            phir_func=phir_R744;
            dphir_dDelta_func=dphir_dDelta_R744;
            dphir2_dDelta2_func=dphir2_dDelta2_R744;
            dphir2_dDelta_dTau_func=dphir2_dDelta_dTau_R744;
            dphir_dTau_func=dphir_dTau_R744;
            dphir2_dTau2_func=dphir2_dTau2_R744;
            phi0_func=phi0_R744;
            dphi0_dDelta_func=dphi0_dDelta_R744;
            dphi02_dDelta2_func=dphi02_dDelta2_R744;
            dphi0_dTau_func=dphi0_dTau_R744;
            dphi02_dTau2_func=dphi02_dTau2_R744;
        }
        else if (!strcmp(Ref,"R718") || !strcmp(Ref,"Water") || !strcmp(Ref,"H2O"))
        {
            FluidType=FLUIDTYPE_REFRIGERANT_PURE;
            p_func=p_Water;
            h_func=h_Water;
            s_func=s_Water;
            u_func=u_Water;
            rho_func=rho_Water;
            cp_func=cp_Water;
            cv_func=cv_Water;
            visc_func=visc_Water;
            k_func=k_Water;
            w_func=w_Water;
            Ttriple_func=Ttriple_Water;
            Tcrit_func=Tcrit_Water;
            pcrit_func=pcrit_Water;
            rhocrit_func=rhocrit_Water;
            MM_func=MM_Water;
            rhosatV_func=rhosatV_Water;
            rhosatL_func=rhosatL_Water;
            psat_func=psat_Water;

            phir_func=phir_Water;
            dphir_dDelta_func=dphir_dDelta_Water;
            dphir2_dDelta2_func=dphir2_dDelta2_Water;
            dphir2_dDelta_dTau_func=dphir2_dDelta_dTau_Water;
            dphir_dTau_func=dphir_dTau_Water;
            dphir2_dTau2_func=dphir2_dTau2_Water;
            phi0_func=phi0_Water;
            dphi0_dDelta_func=dphi0_dDelta_Water;
            dphi02_dDelta2_func=dphi02_dDelta2_Water;
            dphi0_dTau_func=dphi0_dTau_Water;
            dphi02_dTau2_func=dphi02_dTau2_Water;
        }
        else if (!strcmp(Ref,"R134a"))
        {
            FluidType=FLUIDTYPE_REFRIGERANT_PURE;
            p_func=p_R134a;
            h_func=h_R134a;
            s_func=s_R134a;
            u_func=u_R134a;
            rho_func=rho_R134a;
            cp_func=cp_R134a;
            cv_func=cv_R134a;
            visc_func=visc_R134a;
            k_func=k_R134a;
            w_func=w_R134a;
            Ttriple_func=Ttriple_R134a;
            Tcrit_func=Tcrit_R134a;
            pcrit_func=pcrit_R134a;
            rhocrit_func=rhocrit_R134a;
            MM_func=MM_R134a;
            rhosatV_func=rhosatV_R134a;
            rhosatL_func=rhosatL_R134a;
            psat_func=psat_R134a;
            
            Load_R134a(&Fluid);
            phir_func=Fluid.funcs.phir;
            dphir_dDelta_func=Fluid.funcs.dphir_dDelta;
            dphir2_dDelta2_func=Fluid.funcs.dphir2_dDelta2;
            dphir2_dDelta_dTau_func=Fluid.funcs.dphir2_dDelta_dTau;
            dphir_dTau_func=Fluid.funcs.dphir_dTau;
            dphir2_dTau2_func=Fluid.funcs.dphir2_dTau2;
            phi0_func=Fluid.funcs.phi0;
            dphi0_dDelta_func=Fluid.funcs.dphi0_dDelta;
            dphi02_dDelta2_func=Fluid.funcs.dphi02_dDelta2;
            dphi0_dTau_func=Fluid.funcs.dphi0_dTau;
            dphi02_dTau2_func=Fluid.funcs.dphi02_dTau2;
        }
        else if (!strcmp(Ref,"R290"))
        {
            FluidType=FLUIDTYPE_REFRIGERANT_PURE;
            p_func=p_R290;
            h_func=h_R290;
            s_func=s_R290;
            u_func=u_R290;
            rho_func=rho_R290;
            cp_func=cp_R290;
            cv_func=cv_R290;
            visc_func=visc_R290;
            k_func=k_R290;
            w_func=w_R290;
            Ttriple_func=Ttriple_R290;
            Tcrit_func=Tcrit_R290;
            pcrit_func=pcrit_R290;
            rhocrit_func=rhocrit_R290;
            MM_func=MM_R290;
            rhosatV_func=rhosatV_R290;
            rhosatL_func=rhosatL_R290;
            psat_func=psat_R290;

            phir_func=phir_R290;
            dphir_dDelta_func=dphir_dDelta_R290;
            dphir2_dDelta2_func=dphir2_dDelta2_R290;
            dphir2_dDelta_dTau_func=dphir2_dDelta_dTau_R290;
            dphir_dTau_func=dphir_dTau_R290;
            dphir2_dTau2_func=dphir2_dTau2_R290;
            phi0_func=phi0_R290;
            dphi0_dDelta_func=dphi0_dDelta_R290;
            dphi02_dDelta2_func=dphi02_dDelta2_R290;
            dphi0_dTau_func=dphi0_dTau_R290;
            dphi02_dTau2_func=dphi02_dTau2_R290;
        }
        else if (!strcmp(Ref,"R717") || !strcmp(Ref,"NH3") || !strcmp(Ref,"Ammonia") || !strcmp(Ref,"ammonia"))
        {
            FluidType=FLUIDTYPE_REFRIGERANT_PURE;
            p_func=p_R717;
            h_func=h_R717;
            s_func=s_R717;
            u_func=u_R717;
            rho_func=rho_R717;
            cp_func=cp_R717;
            cv_func=cv_R717;
            visc_func=visc_R717;
            k_func=k_R717;
            w_func=w_R717;
            Ttriple_func=Ttriple_R717;
            Tcrit_func=Tcrit_R717;
            pcrit_func=pcrit_R717;
            rhocrit_func=rhocrit_R717;
            MM_func=MM_R717;
            rhosatV_func=rhosatV_R717;
            rhosatL_func=rhosatL_R717;
            psat_func=psat_R717;

            phir_func=phir_R717;
            dphir_dDelta_func=dphir_dDelta_R717;
            dphir2_dDelta2_func=dphir2_dDelta2_R717;
            dphir2_dDelta_dTau_func=dphir2_dDelta_dTau_R717;
            dphir_dTau_func=dphir_dTau_R717;
            dphir2_dTau2_func=dphir2_dTau2_R717;
            phi0_func=phi0_R717;
            dphi0_dDelta_func=dphi0_dDelta_R717;
            dphi02_dDelta2_func=dphi02_dDelta2_R717;
            dphi0_dTau_func=dphi0_dTau_R717;
            dphi02_dTau2_func=dphi02_dTau2_R717;
        }
        else if (!strcmp(Ref,"Air"))
        {
            FluidType=FLUIDTYPE_REFRIGERANT_PSEUDOPURE;
            p_func=p_Air;
            h_func=h_Air;
            s_func=s_Air;
            u_func=u_Air;
            rho_func=rho_Air;
            cp_func=cp_Air;
            cv_func=cv_Air;
            visc_func=visc_Air;
            k_func=k_Air;
            w_func=w_Air;
            Ttriple_func=Ttriple_Air;
            Tcrit_func=Tcrit_Air;
            pcrit_func=pcrit_Air;
            MM_func=MM_Air;
            rhosatV_func=rhosatV_Air;
            rhosatL_func=rhosatL_Air;
            pdp_func=pdp_Air;
            pbp_func=pbp_Air;
        }
        else if (!strcmp(Ref,"R410A"))
        {
            FluidType=FLUIDTYPE_REFRIGERANT_PSEUDOPURE;
            p_func=p_R410A;
            h_func=h_R410A;
            s_func=s_R410A;
            u_func=u_R410A;
            rho_func=rho_R410A;
            cp_func=cp_R410A;
            cv_func=cv_R410A;
            visc_func=visc_R410A;
            k_func=k_R410A;
            w_func=w_R410A;
            Ttriple_func=Ttriple_R410A;
            Tcrit_func=Tcrit_R410A;
            pcrit_func=pcrit_R410A;
            rhocrit_func=rhocrit_R410A;
            MM_func=MM_R410A;
            rhosatV_func=rhosatV_R410A;
            rhosatL_func=rhosatL_R410A;
            pdp_func=p_dp_R410A;
            pbp_func=p_bp_R410A;
            
            Load_R410A(&Fluid);
            phir_func=Fluid.funcs.phir;
            dphir_dDelta_func=Fluid.funcs.dphir_dDelta;
            dphir2_dDelta2_func=Fluid.funcs.dphir2_dDelta2;
            dphir2_dDelta_dTau_func=Fluid.funcs.dphir2_dDelta_dTau;
            dphir_dTau_func=Fluid.funcs.dphir_dTau;
            dphir2_dTau2_func=Fluid.funcs.dphir2_dTau2;
            phi0_func=Fluid.funcs.phi0;
            dphi0_dDelta_func=Fluid.funcs.dphi0_dDelta;
            dphi02_dDelta2_func=Fluid.funcs.dphi02_dDelta2;
            dphi0_dTau_func=Fluid.funcs.dphi0_dTau;
            dphi02_dTau2_func=Fluid.funcs.dphi02_dTau2;
        }
        else if (!strcmp(Ref,"R404A"))
        {
            FluidType=FLUIDTYPE_REFRIGERANT_PSEUDOPURE;
            p_func=p_R404A;
            h_func=h_R404A;
            s_func=s_R404A;
            u_func=u_R404A;
            rho_func=rho_R404A;
            cp_func=cp_R404A;
            cv_func=cv_R404A;
            visc_func=visc_R404A;
            k_func=k_R404A;
            w_func=w_R404A;
            Ttriple_func=Ttriple_R404A;
            Tcrit_func=Tcrit_R404A;
            pcrit_func=pcrit_R404A;
            MM_func=MM_R404A;
            rhosatV_func=rhosatV_R404A;
            rhosatL_func=rhosatL_R404A;
            pdp_func=p_dp_R404A;
            pbp_func=p_bp_R404A;
        }
        else if (!strcmp(Ref,"R407C"))
        {
            FluidType=FLUIDTYPE_REFRIGERANT_PSEUDOPURE;
            p_func=p_R407C;
            h_func=h_R407C;
            s_func=s_R407C;
            u_func=u_R407C;
            rho_func=rho_R407C;
            cp_func=cp_R407C;
            cv_func=cv_R407C;
            visc_func=visc_R407C;
            k_func=k_R407C;
            w_func=w_R407C;
            Ttriple_func=Ttriple_R407C;
            Tcrit_func=Tcrit_R407C;
            pcrit_func=pcrit_R407C;
            MM_func=MM_R407C;
            rhosatV_func=rhosatV_R407C;
            rhosatL_func=rhosatL_R407C;
            pdp_func=p_dp_R407C;
            pbp_func=p_bp_R407C;
        }
        else if (!strcmp(Ref,"R507A"))
        {
            FluidType=FLUIDTYPE_REFRIGERANT_PSEUDOPURE;
            p_func=p_R507A;
            h_func=h_R507A;
            s_func=s_R507A;
            u_func=u_R507A;
            rho_func=rho_R507A;
            cp_func=cp_R507A;
            cv_func=cv_R507A;
            visc_func=visc_R507A;
            k_func=k_R507A;
            w_func=w_R507A;
            Ttriple_func=Ttriple_R507A;
            Tcrit_func=Tcrit_R507A;
            pcrit_func=pcrit_R507A;
            MM_func=MM_R507A;
            rhosatV_func=rhosatV_R507A;
            rhosatL_func=rhosatL_R507A;
            pdp_func=p_dp_R507A;
            pbp_func=p_bp_R507A;
        }
        else
        {
            sprintf(CP_errString,"Refrigerant %s not allowed\n",Ref);
            fprintf(stderr,"%s\n",CP_errString);
            return -1;
        }
        printf("Fluid loaded\n");
    }
    return 0;
}

void Help()
{
    
    //~ printf("CoolProp Help\n");
    //~ printf("CoolProp is written by Ian Bell (ihb2@cornell.edu)\n");
    //~ printf("\n");
    //~ printf("Following the naming conventions of MATLAB linked with REFPROP,\n");
    //~ printf("each output property is represented by one character:\n");
    //~ printf("\n");
    //~ printf("P   Pressure [kPa]\n");
    //~ printf("T   Temperature [K]\n");
    //~ printf("D   Density [kg/m3]\n");
    //~ printf("H   Enthalpy [kJ/kg]\n");
    //~ printf("S   Entropy [kJ/(kg/K)]\n");
    //~ printf("U   Internal energy [kJ/kg]\n");
    //~ printf("C   Cp [kJ/(kg K)]\n");
    //~ printf("O   Cv [kJ/(kg K)]\n");
    //~ printf("K   Ratio of specific heats (Cp/Cv) [-]\n");
    //~ printf("A   Speed of sound [m/s]\n");
    //~ printf("X   liquid phase and gas phase composition (mass fractions)\n");
    //~ printf("V   Dynamic viscosity [Pa*s]\n");
    //~ printf("L   Thermal conductivity [kW/(m K)]\n");
    //~ printf("Q   Quality (vapor fraction) (kg/kg)\n");
    //~ printf("I   Surface tension [N/m]\n");
    //~ printf("F   Freezing point of secondary fluid [K] **NOT IN MATLAB-REFPROP **\n");
    //~ printf("M   Maximum temperature for secondary fluid [K] **NOT IN MATLAB-REFPROP **\n");
    //~ printf("M   Molar mass for non-secondary fluid [g/mol] **NOT IN MATLAB-REFPROP **\n");
    //~ printf("B   Critical Temperature [K] **NOT IN MATLAB-REFPROP **\n");
    //~ printf("E   Critical Pressure [kPa] **NOT IN MATLAB-REFPROP **\n");
    //~ printf("R   Triple point temperature [K] **NOT IN MATLAB-REFPROP **\n");
    //~ printf("\n");
    //~ printf("******** To call **************\n");
    //~ printf("To call the function Props, for instance for R410A at 300K, 400 kPa, you would do:\n");
    //~ printf("Props(\"H\",\"T\",300,\"P\",400,\"R410A\")\n");
    //~ printf("\n");
    //~ printf("Or to call a pure fluid from REFPROP (for instance Propane).  \n");
    //~ printf("The name of the refrigerant is \"REPFROP-\" plus the REFPROP defined name of the fluid, for instance\n");
    //~ printf("\"Propane\" for propane (R290)\n");
    //~ printf("\n");
    //~ printf("See the folder C:\\Program Files\\REFPROP\\fluids for the names of the fluids\n");
    //~ printf("\n");
    //~ printf("To call Propane from REFPROP:\n");
    //~ printf("Props(\"H\",\"T\",300,\"P\",400,\"REFPROP-Propane\")\n");
    //~ printf("\n");
    //~ printf("**************** Inputs ***************\n");
    //~ printf("The limited list of inputs that are allowed are:\n");
    //~ printf("\n");
    //~ printf("Prop1    ||    Prop2\n");
    //~ printf("--------------------\n");
    //~ printf("  T      ||      P\n");
    //~ printf("  T      ||      Q\n");
    //~ printf("  T      ||      D\n");
    
}

int IsFluidType(char *Ref, char *Type)
{
    // Call a function to set the fluid-specific properties
    Props('M','T',0,'P',0,Ref);
    
    if (((FluidType==FLUIDTYPE_REFRIGERANT_PURE) || (FluidType==FLUIDTYPE_REFPROP)) && !strcmp(Type,"PureFluid"))
    {
        return 1;
    }
    else if (FluidType==FLUIDTYPE_BRINE && !strcmp(Type,"Brine"))
    {
        return 1;
    }
    else if (FluidType==FLUIDTYPE_REFRIGERANT_PSEUDOPURE && !strcmp(Type,"PseudoPure"))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
int Phase(double T, double rho, char * Ref)
{
    double rhosatL,rhosatV,p,Tbubble,Tdew;
    if (fabs(SecFluids('D',0,0,Ref))<1e10)
    {
        //It's a secondary fluid, always subcooled
        return PHASE_SUBCOOLED;
    }
    p=Props('P','T',T,'D',rho,Ref);
    if (p>pcrit(Ref))
    {
        return PHASE_SUPERCRITICAL;
    }
    else
    {
        Tbubble=Tsat(Ref, p, 0.0,T);
        Tdew=Tsat(Ref, p, 1.0,T);
        rhosatV=Props('D','T',Tdew,'Q',1,Ref);
        rhosatL=Props('D','T',Tbubble,'Q',0,Ref);

        if (rho<rhosatV)
            return PHASE_SUPERHEATED;
        else if (rho>rhosatL)
            return PHASE_SUBCOOLED;
        else
            return PHASE_TWOPHASE;
    }
}

double Props(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
    double T,p,Q,rhoV,rhoL,Value,rho,pdp,pbp;
    int isTwoPhase;
    char errString[ERRSTRLENGTH];
    
    /*
    Following the naming conventions of MATLAB linked with REFPROP,
    each output property is represented by one character:

    P   Pressure [kPa]
    T   Temperature [K]
    D   Density [kg/m3]
    H   Enthalpy [kJ/kg]
    S   Entropy [kJ/(kg/K)]
    U   Internal energy [kJ/kg]
    C   Cp [kJ/(kg K)]
    O   Cv [kJ/(kg K)]
    K   Ratio of specific heats (Cp/Cv) [-]
    A   Speed of sound [m/s]
    X   liquid phase and gas phase composition (mass fractions)
    V   Dynamic viscosity [Pa*s]
    L   Thermal conductivity [kW/(m K)]
    Q   Quality (vapor fraction) (kg/kg)
    I   Surface tension [N/m]
    F	Freezing point of secondary fluid [K] **NOT IN MATLAB-REFPROP **
    M	Maximum temperature for secondary fluid [K] **NOT IN MATLAB-REFPROP **
    B	Critical Temperature [K] **NOT IN MATLAB-REFPROP **
    E	Critical Pressure [K] **NOT IN MATLAB-REFPROP **
    R   Triple point temperature [K]
    */

    // **********************************************************************************
    // **********************************************************************************
    //                                   REFPROP
    // **********************************************************************************
    // **********************************************************************************
    
    /* 
    If the fluid name is not actually a refrigerant name, but a string beginning with "REFPROP-",
    then REFPROP is used to calculate the desired property.
    */
    #if defined(__ISWINDOWS__)
    if (strncmp(Ref,"REFPROP-",8)==0)  // First eight characters match "REFPROP-"
    {
    #else
    if (0) // Automatically skip it because REFPROP is not supported on this platform
    {
    #endif
        FluidType=FLUIDTYPE_REFPROP;
        #if defined(__ISWINDOWS__)
        return REFPROP(Output,Name1,Prop1,Name2,Prop2,Ref);
        #endif
    }

    // **********************************************************************************
    // **********************************************************************************
    //                                Normal Property evaluation
    // **********************************************************************************
    // **********************************************************************************

    // It's a brine, call the brine routine
    else if (!strcmp(Ref,"HC-10") || (Ref[0]=='E' && Ref[1]=='G') || 
        (Ref[0]=='P' && Ref[1]=='G') || strncmp(Ref,"Methanol",8)==0 || strncmp(Ref,"NH3/H2O",7)==0)
    {
        if (Name1!='T' || Name2!='P')
        {
            printf("Warning: For brine, Name1 must be 'T' and Name2 must be 'P'\n");
        }
        FluidType=FLUIDTYPE_BRINE;
        return SecFluids(Output,Prop1,Prop2,Ref);
    }
    else 
    {
        T=Prop1; 
        // Load the fluid-specific parameters and function pointers
        LoadFluid(Ref);
        
        // Check if it is an output that doesn't require a state input
        // Deal with it and quit
        switch (Output)
        {
            case 'M':
                Value=MM_func(); break;
            case 'E':
                Value=pcrit_func(); break;
            case 'B':
                Value=Tcrit_func(); break;
            case 'R':
                Value=Ttriple_func(); break;
        }

        if (Name1=='T' && (Name2=='P' || Name2=='D' || Name2=='Q'))
        {
            // Temperature and Pressure are the inputs
            if (Name2=='P')
            {
                if (FlagUseSinglePhaseLUT==1)
                {
                    BuildLookupTable(Ref,&Fluid);
                    Value= LookupValue(Output, T, p, Ref, &Fluid);
                    return Value;
                }
                else
                {
                    // Collect the value for pressure
                    p=Prop2;
                    //First find density as a function of temp and pressure (all parameters need it)
                    rho=get_Delta(T,p)*Fluid.rhoc;
                    if (Output=='D')
                    {
                        return rho;
                    }
                    else
                    {
                        //Recursively call Props with the calculated temperature and density
                        return Props(Output,'T',T,'D',rho,Ref);
                    }
                }
                return Value;
            }
            // Temperature and density are the inputs
            else if (Name2=='D')
            {
                rho=Prop2;
                
                // Determine if it is somewhere in the two-phase region.

                // First check if the condition is well away from saturation 
                // by using the ancillary equations because this is much faster
                //printf("T: %g L: %g V: %g \n",T,rhosatL_func(T),rhosatV_func(T));
                if (rho>1.002*rhosatL_func(T) || rho<0.998*rhosatV_func(T))
                {
                    isTwoPhase=0;
                }
                else
                {
                    if (!strcmp(Ref,"R404A") || !strcmp(Ref,"R410A") || !strcmp(Ref,"R407C") 
                        || !strcmp(Ref,"R507A") || !strcmp(Ref,"Air"))
                    {
                        rhoV=rhosatV_func(T);
                        rhoL=rhosatL_func(T);
                    }
                    else
                    {
                        rhosatPure(Ref,T,&rhoL,&rhoV,&p);
                    }
                    if (rho<=rhoL && rho>=rhoV)
                        isTwoPhase=1;
                    else
                        isTwoPhase=0;
                    printf("Is two phase, rhol: %g rhoV: %g\n",rhoL,rhoV);
                }
                if (isTwoPhase==0)
                {
                    // It is not two-phase, and use EOS or transport relations
                    switch (Output)
                    {
                        case 'P':
                            Value=Pressure_Trho(T,rho); break;
                        case 'H':
                            Value=Enthalpy_Trho(T,rho); break;
                        case 'S':
                            Value=Entropy_Trho(T,rho); break;
                        case 'U':
                            Value=IntEnergy_Trho(T,rho); break;
                        case 'C':
                            Value=SpecHeatP_Trho(T,rho); break;
                        case 'O':
                            Value=SpecHeatV_Trho(T,rho); break;
                        case 'A':
							Value=w_func(T,rho,TYPE_Trho); break;
                        case 'V':
                            Value=visc_func(T,rho,TYPE_Trho); break;
                        case 'L':
                            Value=k_func(T,rho,TYPE_Trho); break;
                        default:
                            strcpy(errString,"Invalid Output Name");
                            return -100;
                    }
                    return Value;
                }
                else
                {
                    // It is two-phase. Find the quality and call Props again with the quality
                    Q=(1/rho-1/rhoL)/(1/rhoV-1/rhoL);
                    return Props(Output,'T',T,'Q',Q,Ref);
                }
            }

            // Temperature and quality are the inputs
            else if (Name2=='Q')
            {
                Q=Prop2;
                if (!strcmp(Ref,"R290") || !strcmp(Ref,"Argon") ||!strcmp(Ref,"Nitrogen")  ||!strcmp(Ref,"R134a") ||!strcmp(Ref,"R717") || !strcmp(Ref,"CO2") || !strcmp(Ref,"R744") || !strcmp(Ref,"Water") )
                {
                    if (FlagUseSaturationLUT)
                    {
                        // Use the saturation Lookup Table;
                        rhoL=ApplySaturationLUT("rhoL","T",T,Ref);
                        rhoV=ApplySaturationLUT("rhoV","T",T,Ref);
                        p=ApplySaturationLUT("p","T",T,Ref);
                    }
                    else
                    {
                        rhosatPure(Ref,T,&rhoL,&rhoV,&p);
                    }
                }
                else
                {
                    pbp=pbp_func(T);
                    pdp=pdp_func(T);
                    rhoV=rhosatV_func(T);
                    rhoL=rhosatL_func(T);
                    p=Q*pdp+(1-Q)*pbp;
                }
                rho=1/(Q/rhoV+(1-Q)/rhoL);

                switch (Output)
                {
                    case 'P':
                        Value=p; 
                        break;
                    case 'H':
                        Value=Q*h_func(T,rhoV,TYPE_Trho)+(1-Q)*h_func(T,rhoL,TYPE_Trho); 
                        break;
                    case 'D':
                        Value=rho; break;
                    case 'S':
                        Value=Q*s_func(T,rhoV,TYPE_Trho)+(1-Q)*s_func(T,rhoL,TYPE_Trho); 
                        break;
                    case 'U':
                        Value=Q*u_func(T,rhoV,TYPE_Trho)+(1-Q)*u_func(T,rhoL,TYPE_Trho); 
                        break;
                    case 'C':
                        Value=Q*cp_func(T,rhoV,TYPE_Trho)+(1-Q)*cp_func(T,rhoL,TYPE_Trho); 
                        break;
                    case 'O':
                        Value=Q*cv_func(T,rhoV,TYPE_Trho)+(1-Q)*cv_func(T,rhoL,TYPE_Trho); 
                        break;
                    case 'V':
                        Value=Q*visc_func(T,rhoV,TYPE_Trho)+(1-Q)*visc_func(T,rhoL,TYPE_Trho); 
                        break;
                    case 'L':
                        Value=Q*k_func(T,rhoV,TYPE_Trho)+(1-Q)*k_func(T,rhoL,TYPE_Trho); 
                        break;
                    case 'A':
                        Value=Q*w_func(T,rhoV,TYPE_Trho)+(1-Q)*w_func(T,rhoL,TYPE_Trho); 
                        break;
                    case 'M':
                        Value=MM_func(); break;
                    case 'E':
                        Value=pcrit_func(); break;
                    case 'B':
                        Value=Tcrit_func(); break;
                    case 'R':
                        Value=Ttriple_func(); break;
                    default:
                        strcpy(errString,"Invalid Output Name");
                        return -100;
                }
                return Value;
            }
        }
        else if (Output=='T' && Name1=='P' && Name2=='Q')
        {
            return Tsat(Ref,Prop1,Prop2,0);
        }
        else if (Name1=='P' && Name2=='Q')
        {
            T=Tsat(Ref,Prop1,Prop2,0);
            return Props(Output,'T',T,'Q',Prop2,Ref);
        }
        else
        {
            fprintf(stderr,"Names of input properties invalid (%c,%c) with refrigerant %s.  Valid choices are T,P or T,Q or T,D or P,Q",Name1,Name2,Ref);
            return _HUGE;
        }
    }
    return 0;
}

double pcrit(char *Ref)
{
    // Call a function to set the global constants
    Props('M','T',0,'P',0,Ref);
    
    // Brines do not have critical pressures, set it to a big number
    if (IsFluidType(Ref,"Brine"))
    {
        //ERROR!
    	return 100000;
    }
    #if defined(__ISWINDOWS__) 
    if (strncmp(Ref,"REFPROP-",8)==0)  // First eight characters match "REFPROP-"
    {
        char *REFPROPRef=NULL,*RefCopy=NULL;
        double pcrit;
        RefCopy=malloc(strlen(Ref)+1);
        strcpy(RefCopy,Ref);
        REFPROPRef = strtok(RefCopy,"-");
        REFPROPRef = strtok(NULL,"-");
        // 'E' is the code for the critical pressure (ran out of sensible characters).  
        pcrit=REFPROP('E','T',0.0,'Q',0.0,REFPROPRef);
        free(RefCopy);
        return pcrit;
    }
    #else
    if (0){} // Automatically skip it because REFPROP is not supported on this platform
    #endif
    else{
        return Props('E','T',0.0,'Q',0.0,Ref);
    }
    //ERROR!
    return 0;
}

double Tcrit(char *Ref)
{	
    // Call a function to set the global constants
    Props('M','T',0,'P',0,Ref);
    
    // Brines do not have critical temperatures, set it to a big number
    if (IsFluidType(Ref,"Brine"))
    {
    	//ERROR!
    	return 100000;
    }
    #if defined(__ISWINDOWS__) 
    if (strncmp(Ref,"REFPROP-",8)==0)  // First eight characters match "REFPROP-"
    {
        char *REFPROPRef=NULL,*RefCopy=NULL;
        double Tc;
        RefCopy=malloc(strlen(Ref)+1);
        strcpy(RefCopy,Ref);
        REFPROPRef = strtok(RefCopy,"-");
        REFPROPRef = strtok(NULL,"-");
        // 'B' is the code for the critical pressure (ran out of sensible characters).
        Tc=REFPROP('B','T',0.0,'Q',0.0,REFPROPRef);
        free(RefCopy);
        return Tc;
    }
    #else
    if (0){} // Automatically skip it because REFPROP is not supported on this platform
    #endif
    else
    {
        return Props('B','T',0.0,'P',0.0,Ref);
    }
    //ERROR!
    return 0;
}
double Ttriple(char *Ref)
{
    // Call a function to set the global constants
    Props('M','T',0,'P',0,Ref);
    
    // Brines do not have triple point temperatures, set it to a big number
    if (IsFluidType(Ref,"Brine"))
    {
        return 100000;
    }
    #if defined(__ISWINDOWS__) 
    if (strncmp(Ref,"REFPROP-",8)==0)  // First eight characters match "REFPROP-"
    {
        char *REFPROPRef=NULL,*RefCopy=NULL;
        double Ttriple;
        RefCopy=malloc(strlen(Ref)+1);
        strcpy(RefCopy,Ref);
        REFPROPRef = strtok(RefCopy,"-");
        REFPROPRef = strtok(NULL,"-");
        // 'R' is the code for the triple point temperature (ran out of sensible characters).
        Ttriple=REFPROP('R','T',0.0,'Q',0.0,REFPROPRef);
        free(RefCopy);
        return Ttriple;
    }
    #else
    if (0){} // Automatically skip it because REFPROP is not supported on this platform
    #endif
    else
    {
        return Props('R','T',0.0,'Q',0.0,Ref);
    }
    //ERROR!
    return 0;
}

int errCode(char * Ref)
{
    if (!strcmp(Ref,"R290"))
        return errCode_R290();
    if (!strcmp(Ref,"R744"))
        return errCode_R744();
    if (!strcmp(Ref,"R717"))
        return errCode_R717();
    if (!strcmp(Ref,"R32"))
        return errCode_R32();
    if (!strcmp(Ref,"R404A"))
        return errCode_R404A();
    return -1;
}


double T_hp(char *Ref, double h, double p, double T_guess)
{
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f,T=300;
    int iter=1;
    while ((iter<=3 || change>eps) && iter<100)
    {
        if (iter==1){x1=T_guess; T=x1;}
        if (iter==2){x2=T_guess+0.1; T=x2;}
        if (iter>2) {T=x2;}
            f=Props('H','T',T,'P',p,Ref)-h;
        if (iter==1){y1=f;}
        if (iter>1)
        {
            y2=f;
            x3=x2-y2/(y2-y1)*(x2-x1);
            change=fabs(y2/(y2-y1)*(x2-x1));
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
        if (iter>60)
        {
            //ERROR
        	//printf("%d: T_hp not converging with inputs(%s,%g,%g,%g) value: %0.12g\n",iter,Ref,h,p,T_guess,f);
        }
    }
    return T;
}

double h_sp(char *Ref, double s, double p, double T_guess)
{
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f=999,T=300;
    int iter=1;

    
    while ((iter<=3 || change>eps) && iter<100)
    {
        if (iter==1){x1=T_guess; T=x1;}
        if (iter==2){x2=T_guess+1.0; T=x2;}
        if (iter>2) {T=x2;}

            // Find the temperature which gives the same entropy
            f=Props('S','T',T,'P',p,Ref)-s;

        if (iter==1){y1=f;}
        if (iter>1)
        {
            y2=f;
            x3=x2-y2/(y2-y1)*(x2-x1);
            change=fabs(y2/(y2-y1)*(x2-x1));
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
        if (iter>50)
        {
        	//ERROR
            printf("h_sp not converging with inputs(%s,%g,%g,%g)\n",Ref,s,p,T_guess);
        }
    }
    return Props('H','T',T,'P',p,Ref);
}

double Tsat(char *Ref, double p, double Q, double T_guess)
{
    double x1=0,x2=0,y1=0,y2=0,tau1,tau3,tau2,logp1,logp2,logp3,T1,T2,T3;
    int i,k;
    double Tc,Tmax,Tmin;
    // pointers to the x and y vectors for later
    double (*y)[NLUT],(*x)[NLUT],x0,y0;

    // Brines do not have saturation temperatures, set it to a big number
    if (IsFluidType(Ref,"Brine"))
    {
        return 100000;
    }
    #if defined(__ISWINDOWS__) 
        // It's a REFPROP fluid - use REFPROP to do all the calculations
        if (!strncmp(Ref,"REFPROP-",8))
        {
            return REFPROP('T','P',p,'Q',Q,Ref);
        }
    #endif
    
    // Do reverse interpolation in the Saturation Lookup Table
    if (FlagUseSaturationLUT && FluidType==FLUIDTYPE_REFRIGERANT_PURE)
    {
        // First try to build the LUT, get index of LUT
        k=BuildSaturationLUT(Ref);
        
        // Then find the index to the left of the pressure
        for(i=0;i<NLUT;i++)
        {
            if (psat_LUT[k][i]>=p) break;
        }
        
        x=&(psat_LUT[k]);
        y=&(Tsat_LUT[k]);
        
        // Need a three-point set to interpolate using a quadratic.
        if (i<NLUT-2)
        {
            x0=(*x)[i];
            x1=(*x)[i+1];
            x2=(*x)[i+2];
            y0=(*y)[i];
            y1=(*y)[i+1];
            y2=(*y)[i+2];
        }
        else
        {
            // Go "backwards" with the interpolation range
            x0=(*x)[i+1];
            x1=(*x)[i];
            x2=(*x)[i-1];
            y0=(*y)[i+1];
            y1=(*y)[i];
            y2=(*y)[i-1];
        }
        return QuadInterpolate(x0,x1,x2,y0,y1,y2,p);
        
    }
    Tc=Tcrit(Ref);
    Tmax=Tc-0.000001;
    Tmin=Ttriple(Ref)+1;
    
    // Plotting Tc/T versus log(p) tends to give very close to straight line
    // Use this fact to figure out a reasonable guess temperature
    
    if (FluidType==FLUIDTYPE_REFRIGERANT_PURE)
    {
        T1=Ttriple(Ref)+1;
        T3=Tcrit(Ref)-1;
        T2=(T1+T3)/2;
        tau1=Tc/T1;
        tau2=Tc/T2;
        tau3=Tc/T3;
        logp1=log(psat_func(T1));
        logp2=log(psat_func(T2));
        logp3=log(psat_func(T3));

        T_guess=Tc/QuadInterpolate(logp1,logp2,logp3,tau1,tau2,tau3,log(p));
        if (T_guess+5<Tmax)
            Tmax=T_guess+5;
        if (T_guess-5>Tmin)
            Tmin=T_guess-5;
    }

    return _Dekker_Tsat(Tmin,Tmax,0.001,p,Q,Ref);
    
    //~ *** This is old code, kept as a reference
    //~ while ((iter<=3 || exp(f)-1>eps) && iter<100)
    //~ {
        //~ if (iter==1){x1=Tc/Tmin; tau=x1;}
        //~ if (iter==2){x2=Tc/Tmax; tau=x2;}
        //~ if (iter>2) {tau=x2;}

            //~ T=Tc/tau;
            //~ f=log(Props('P','T',T,'Q',Q,Ref))-log(p);

            //~ //printf("T: %g f %g\n",T,f);

        //~ if (iter==1){y1=f;}
        //~ if (iter>1)
        //~ {
            //~ y2=f;
            //~ x3=x2-y2/(y2-y1)*(x2-x1);
            //~ change=fabs(y2/(y2-y1)*(x2-x1));
            //~ y1=y2; x1=x2; x2=x3;
        //~ }
        //~ iter=iter+1;
        //~ if (iter>50)
        //~ {
            //~ printf("Tsat not converging with inputs(%s,%g,%g,%g)\n",Ref,p,Q,T_guess);
        //~ }
    //~ }
    //~ return Tc/x2;
}

double K2F(double T)
{
    return T * 9 / 5 - 459.67;
}

double F2K(double T_F)
{
    return (T_F + 459.67) * 5 / 9;
}

void PrintSaturationTable(char *FileName, char * Ref,double Tmin, double Tmax)
{

    double T;
    FILE *f;
    f=fopen(FileName,"w");
    fprintf(f,"T,pL,pV,rhoL,rhoV,uL,uV,hL,hV,sL,sV,cvL,cvV,cpL,cpV,kL,kV,muL,muV\n");
    fprintf(f,"[K],[kPa],[kPa],[kg/m^3],[kg/m^3],[kJ/kg],[kJ/kg],[kJ/kg],[kJ/kg],[kJ/kg-K],[kJ/kg-K],[kJ/kg-K],[kJ/kg-K],[kJ/kg-K],[kJ/kg-K],[kW/m-K],[kW/m-K],[Pa-s],[Pa-s]\n");
    fprintf(f,"--------------------------------------------------------------------------\n");

    for (T=Tmin;T<Tmax;T+=1.0)
    {
    fprintf(f,"%0.3f,",T);
    fprintf(f,"%0.6f,",Props('P','T',T,'Q',0,Ref));
    fprintf(f,"%0.6f,",Props('P','T',T,'Q',1,Ref));
    fprintf(f,"%0.6f,",Props('D','T',T,'Q',0,Ref));
    fprintf(f,"%0.6f,",Props('D','T',T,'Q',1,Ref));
    fprintf(f,"%0.6f,",Props('U','T',T,'Q',0,Ref));
    fprintf(f,"%0.6f,",Props('U','T',T,'Q',1,Ref));
    fprintf(f,"%0.6f,",Props('H','T',T,'Q',0,Ref));
    fprintf(f,"%0.6f,",Props('H','T',T,'Q',1,Ref));
    fprintf(f,"%0.6f,",Props('S','T',T,'Q',0,Ref));
    fprintf(f,"%0.6f,",Props('S','T',T,'Q',1,Ref));
    fprintf(f,"%0.6f,",Props('O','T',T,'Q',0,Ref));
    fprintf(f,"%0.6f,",Props('O','T',T,'Q',1,Ref));
    fprintf(f,"%0.6f,",Props('C','T',T,'Q',0,Ref));
    fprintf(f,"%0.6f,",Props('C','T',T,'Q',1,Ref));
    fprintf(f,"%0.9f,",Props('L','T',T,'Q',0,Ref));
    fprintf(f,"%0.9f,",Props('L','T',T,'Q',1,Ref));
    fprintf(f,"%0.6g,",Props('V','T',T,'Q',0,Ref));
    fprintf(f,"%0.6g,",Props('V','T',T,'Q',1,Ref));
    fprintf(f,"\n");
    }
    fclose(f);
}

double DerivTerms(char *Term,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
    // Pointers to the functions
    double (*dhdT)(double,double,int);
    double (*dpdT)(double,double,int);
    double (*dhdrho)(double,double,int);
    double (*dpdrho)(double,double,int);
    
    // Wire up the pointers for the given refrigerant
    if (!strcmp(Ref,"Argon"))
    {
        dhdT=dhdT_Argon;
        dpdT=dpdT_Argon;
        dhdrho=dhdrho_Argon;
        dpdrho=dpdrho_Argon;
    }
    else if (!strcmp(Ref,"Nitrogen") || !strcmp(Ref,"N2"))
    {
        dhdT=dhdT_Nitrogen;
        dpdT=dpdT_Nitrogen;
        dhdrho=dhdrho_Nitrogen;
        dpdrho=dpdrho_Nitrogen;
    }
    else if (!strcmp(Ref,"R744") || !strcmp(Ref,"CO2"))
    {
        dhdT=dhdT_R744;
        dpdT=dpdT_R744;
        dhdrho=dhdrho_R744;
        dpdrho=dpdrho_R744;
    }
    else if (!strcmp(Ref,"R410A"))
    {
        dhdT=dhdT_R410A;
        dpdT=dpdT_R410A;
        dhdrho=dhdrho_R410A;
        dpdrho=dpdrho_R410A;
    }
    else
    {
        printf("Bad Refrigerant Name in DerivTerms [%s]\n",Ref);
    }
    
    if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdT"))
        return dhdT(Prop1,Prop2,TYPE_TPNoLookup);
    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdrho"))
        return dhdrho(Prop1,Prop2,TYPE_TPNoLookup);
    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdT"))
        return dpdT(Prop1,Prop2,TYPE_TPNoLookup);
    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdrho"))
        return dpdrho(Prop1,Prop2,TYPE_TPNoLookup);
    
    else if (Name1=='T' && Name2=='D' && !strcmp(Term,"dhdT"))
        return dhdT(Prop1,Prop2,TYPE_Trho);
    else if (Name1=='T' && Name2=='D' && !strcmp(Term,"dhdrho"))
        return dhdrho(Prop1,Prop2,TYPE_Trho);
    else if (Name1=='T' && Name2=='D' && !strcmp(Term,"dpdT"))
        return dpdT(Prop1,Prop2,TYPE_Trho);
    else if (Name1=='T' && Name2=='D' && !strcmp(Term,"dpdrho"))
        return dpdrho(Prop1,Prop2,TYPE_Trho);
    
    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdTnum"))
        return dhdT(Prop1,Prop2,99);
    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdrhonum"))
        return dhdrho(Prop1,Prop2,99);
    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdTnum"))
        return dpdT(Prop1,Prop2,99);
    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdrhonum"))
        return dpdrho(Prop1,Prop2,99);

    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdTnum"))
        return dhdT(Prop1,Prop2,99);
    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdrhonum"))
        return dhdrho(Prop1,Prop2,99);
    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdTnum"))
        return dpdT(Prop1,Prop2,99);
    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdrhonum"))
        return dpdrho(Prop1,Prop2,99);
    else
    {
        printf("Bad pair of properties[%c,%c] and derivative name [%s]\n",Name1,Name2,Term);
        return _HUGE;
    }	
}

static void swap(double *x, double *y)
{
    double temp;
    temp = *x;
    *x = *y;
    *y = temp;
}

double _Dekker_Tsat(double x_min, double x_max, double eps, double p, double Q,char *Ref)
{
    double a_k,b_k,f_ak,f_bk,f_bkn1,error, 
        b_kn1,b_kp1,s,m,a_kp1,f_akp1,f_bkp1,x,Tc;

    int iter=1;

    Tc=Tcrit_func();
    x_min=Tc/x_min;
    x_max=Tc/x_max;

    // Loop for the solver method
    while ((iter <= 1 || fabs(error) > eps) && iter < 100)
    {
        // Start with the maximum value
        if (iter == 1)
        {
            a_k = x_max;
            x = a_k;
        }
        // End with the minimum value
        if (iter == 2)
        {
            b_k = x_min;
            x = b_k;
        }
        if (iter > 2)
            x = b_k;

        // Evaluate residual
        error = log(Props('P','T',Tc/x,'Q',Q,Ref))-log(p);

        // First time through, store the outputs
        if(iter == 1)
        {
            f_ak = error;
            b_kn1 = a_k;
            f_bkn1 = error;
        }
        if (iter > 1)
        {
            f_bk = error;
            //Secant solution
            s = b_k - (b_k - b_kn1) / (f_bk - f_bkn1) * f_bk;
            //Midpoint solution
            m = (a_k + b_k) / 2.0;

            if (s > b_k && s < m)
            {
                //Use the secant solution
                b_kp1 = s;
            }
            else
            {
                //Use the midpoint solution
                b_kp1 = m;
            }

            //See if the signs of iterate and contrapoint are the same
            f_bkp1 = log(Props('P','T',Tc/b_kp1,'Q',Q,Ref))-log(p);

            if (f_ak / fabs(f_ak) != f_bkp1 / fabs(f_bkp1))
            {
                // If a and b have opposite signs, 
                //  keep the same contrapoint
                a_kp1 = a_k;
                f_akp1 = f_ak;
            }
            else
            {
                //Otherwise, keep the iterate
                a_kp1 = b_k;
                f_akp1 = f_bk;
            }

            if (fabs(f_akp1) < fabs(f_bkp1))
            {
                //a_k+1 is a better guess than b_k+1, so swap a and b values
                swap(&a_kp1, &b_kp1);
                swap(&f_akp1, &f_bkp1);
            }

            //Update variables
            //Old values
            b_kn1 = b_k;
            f_bkn1 = f_bk;
            //values at this iterate
            b_k = b_kp1;
            a_k = a_kp1;
            f_ak = f_akp1;
            f_bk = f_bkp1;
        }
        iter++;
        if (iter>90 && fabs(error)>eps)
        {
            printf("Dekker for Tsat has failed with inputs (%g,%g,%g,%g,%g,%s); value: %g\n",x_min,x_max,eps,p,Q,Ref,error);
        }
    }
    return Tc/b_k;
}
