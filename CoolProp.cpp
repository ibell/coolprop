#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#endif

#include "CoolProp.h"

#if defined(__ISWINDOWS__)
#include <windows.h>
#include "REFPROP.h"
#endif

#include <iostream>
#include <stdlib.h>
#include "string.h"
#include <stdio.h>
#include "CoolPropTools.h"

using namespace std;
char LoadedFluid[255];

double R_u=8.314472; //From Lemmon et al 2009 (propane)

// Function prototypes
double _Dekker_Tsat(double x_min, double x_max, double eps, double p, double Q,char *Ref);
int LoadFluid(char *Ref);
static void rhosatPure(char *Ref, double T, double *rhoLout, double *rhoVout, double *pout);

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

double (*Viscosity_Trho)(double,double);
double (*Conductivity_Trho)(double,double);
double (*rhosatV_func)(double);
double (*rhosatL_func)(double);

#define NLUTFLUIDS 30 // How many saturation curve LUTs you can have in one instance of CoolProp
#define NLUT 300 // How many values to use for each lookup table

double Tsat_LUT[NLUTFLUIDS][NLUT],rhosatL_LUT[NLUTFLUIDS][NLUT],rhosatV_LUT[NLUTFLUIDS][NLUT],psat_LUT[NLUTFLUIDS][NLUT];
char RefSatLUT[NLUTFLUIDS][255]; // The names of the fluids that are loaded
int FlagUseSaturationLUT=0; //Default to not use LUT
int FlagUseSinglePhaseLUT=0; //Default to not use LUT
int FlagDebug=1;

char CP_errString[5000];
int ErrorFlag;

// The structure that contains all the fluid-specific parameters and pointers to functions
static struct fluidParamsVals Fluid;
    
static double get_Delta(double T, double p)
{
    double eps=1e-10, tau,delta_guess,M,Tc,rhoc,R;
    int counter=1;
    double r1,r2,r3,delta1,delta2,delta3;
    
    //Fluid Properties
    M=Fluid.MM;
    rhoc=Fluid.rhoc;
    Tc=Fluid.Tc;
    R=R_u/M;
        
    if (T>(Fluid.Tc) )
    {
    	//0.7 for compressibility factor since supercritical fluid is quite dense
        delta_guess=p/(8.314/M*T)/Fluid.rhoc/0.7;
    }
    else
    {
        if (Fluid.Type == FLUIDTYPE_REFRIGERANT_PURE && p<Fluid.funcs.psat(T))
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
    r1=p/(delta1*rhoc*R*T)-1.0-delta1*dphir_dDelta_func(tau,delta1);
    r2=p/(delta2*rhoc*R*T)-1.0-delta2*dphir_dDelta_func(tau,delta2);
    while( counter==1 || fabs(r2)>eps)
    {
        delta3=delta2-r2/(r2-r1)*(delta2-delta1);
        r3=p/(delta3*rhoc*R*T)-1.0-delta3*dphir_dDelta_func(tau,delta3);
        delta1=delta2;
        delta2=delta3;
        r1=r2;
        r2=r3;
        counter=counter+1;
        if (counter>40)
        {
        	//ERROR
        	if (fabs(r3)/10<eps)
        		return _HUGE;
        	else
            {
				printf("counter for get_Delta > 40, fcn is %g and delta is %g w/ inputs of T: %g p: %g\n",r3,delta3,T,p);
        		return _HUGE;
            }
        }
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
        if (!strcmp(RefSatLUT[k],Ref)) return k;
        // You made it to an empty string, need a LUT, jump out of for loop
        if (!strcmp(RefSatLUT[k],"")) break;
    }
    if (k==NLUTFLUIDS)
    {
        // Ran out of fluids
        printf("Sorry, ran out of saturation LUT slots, increase NLUTFLUIDS and recompile\n");
    }
    
    // Therefore it hasn't been found yet, assign the refrigerant name
    strcpy(RefSatLUT[k],Ref);
    
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
    while (fabs(error)>1e-10)
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

static void rhosatPure(char *Ref, double T, double *rhoLout, double *rhoVout, double *pout)
{
    // Only works for pure fluids (no blends)
    // At equilibrium, saturated vapor and saturated liquid are at the same pressure and the same Gibbs energy
    double rhoL,rhoV,p,error=999,x1,x2,x3,y1,y2,f,p_guess;
    int iter;
    char Local_errString[300];

    if (T>Fluid.Tc || T<Fluid.Tt)
    {
    	sprintf(Local_errString,"rhosatPure: Temperature [%g] is out of two-phase range [%g,%g]",T,Fluid.Tt,Fluid.Tc);
		Append2ErrorString(Local_errString);
		return;
    }

    // Use the density ancillary function as the starting point for the secant solver
    rhoL=rhosatL_func(T);
    rhoV=rhosatV_func(T);
    p_guess=Pressure_Trho(T,rhoV);

    if (!ValidNumber(rhoL) || !ValidNumber(rhoV))
	{
		//ERROR
		sprintf(Local_errString,"rhosatPure: rhoL [%g] or rhoV [%g] is invalid number\n",rhoL,rhoV);
		Append2ErrorString(Local_errString);
		return;
	}
    iter=1;
    // Use a secant method to obtain pressure
    while ((iter<=3 || fabs(error)>1e-10) && iter<100)
    {
        if (iter==1){x1=p_guess; p=x1;}
        if (iter==2){x2=1.0001*p_guess; p=x2;}
        if (iter>2) {p=x2;}
            //Recalculate the densities based on the current pressure
            rhoL=Density_Tp(T,p,rhoL);
            rhoV=Density_Tp(T,p,rhoV);
            // Residual between saturated liquid and saturated vapor gibbs function
            f=Gibbs_Trho(T,rhoL)-Gibbs_Trho(T,rhoV);
            if (!ValidNumber(rhoL) || !ValidNumber(rhoV) || !ValidNumber(f))
            {
            	//ERROR
            	sprintf(Local_errString,"rhosatPure: rhoL [%g] rhoV [%g] or f [%g] is invalid number\n",rhoL,rhoV,f);
				Append2ErrorString(Local_errString);
				return;
            }
            if (iter>100)
            {
            	//ERROR
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
        if (iter>100)
        {
        	//ERROR
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
	char Local_errString[100];

    //Always load the fluid if possible
    if (0)//(!strcmp(LoadedFluid,Ref))
    {
        // Already Loaded, don't do anything else
        return OK;
    }
    else
    {
    	//Copy the refrigerant name
        strcpy(LoadedFluid,Ref);
        // Wire up the function pointers for the given refrigerant
        if (!strcmp(Ref,"Argon"))
        {
        	Load_Argon(&Fluid);
        }
        else if (!strcmp(Ref,"Nitrogen") || !strcmp(Ref,"N2"))
        {
        	Load_Nitrogen(&Fluid);
        }
        else if (!strcmp(Ref,"R744") || !strcmp(Ref,"CO2"))
        {
        	Load_R744(&Fluid);
        }
        else if (!strcmp(Ref,"R718") || !strcmp(Ref,"Water") || !strcmp(Ref,"H2O"))
        {
        	Load_Water(&Fluid);
        }
        else if (!strcmp(Ref,"R134a"))
        {
            Load_R134a(&Fluid);
        }
        else if (!strcmp(Ref,"R290"))
        {
        	Load_R290(&Fluid);
        }
        else if (!strcmp(Ref,"R717") || !strcmp(Ref,"NH3") || !strcmp(Ref,"Ammonia") || !strcmp(Ref,"ammonia"))
        {
        	Load_R717(&Fluid);
        }
        else if (!strcmp(Ref,"Air"))
        {
        	Load_Air(&Fluid);
        }
        else if (!strcmp(Ref,"R410A"))
        {
            Load_R410A(&Fluid);
        }
        else if (!strcmp(Ref,"R404A"))
        {
        	Load_R404A(&Fluid);
        }
        else if (!strcmp(Ref,"R407C"))
        {
        	Load_R407C(&Fluid);
        }
        else if (!strcmp(Ref,"R507A"))
        {
        	Load_R507A(&Fluid);
        }
        else
        {
            sprintf(Local_errString,"Refrigerant %s not allowed",Ref);
            Append2ErrorString(Local_errString);
            return FAIL;
        }
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
		Viscosity_Trho=Fluid.funcs.visc;
		Conductivity_Trho=Fluid.funcs.cond;
		rhosatV_func=Fluid.funcs.rhosatV;
		rhosatL_func=Fluid.funcs.rhosatL;
		psat_func=Fluid.funcs.psat;
		pbp_func=Fluid.funcs.p_bp;
		pdp_func=Fluid.funcs.p_dp;
        //printf("Loaded Fluid %s\n",Ref);
        return OK;
    }
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

static int IsBrine(char *Ref)
{
    if (
        strcmp(Ref,"HC-10")==0 || 
        strncmp(Ref,"EG",2)==0 || 
        strncmp(Ref,"PG",2)==0 || 
        strncmp(Ref,"Methanol",8)==0 || 
        strncmp(Ref,"NH3/H2O",7)==0
       )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
static int IsREFPROP(char *Ref)
{
    if (strncmp(Ref,"REFPROP-",8)==0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
static int IsPseudoPure(char *Ref)
{
    if (strncmp(Ref,"Air",3)==0 ||
        strncmp(Ref,"R410A",4)==0 ||
        strncmp(Ref,"R407C",4)==0 ||
        strncmp(Ref,"R507A",4)==0 ||
        strncmp(Ref,"R404A",4)==0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
int IsFluidType(char *Ref,char *Type)
{
    if (IsBrine(Ref) && !strcmp(Type,"Brine"))
    {
        return 1;
    }
    else if (IsPseudoPure(Ref) && !strcmp(Type,"PseudoPure"))
    {
        return 1;
    }
    else if (((Fluid.Type==FLUIDTYPE_REFRIGERANT_PURE) || IsREFPROP(Ref)) && !strcmp(Type,"PureFluid"))
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

void PropsV(char *Output,char Name1, double *Prop1, int len1, char Name2, double *Prop2, int len2, char * Ref, double *OutVec, int n)
{
    // This is a wrapper function that allows for vectors to be passed through 
    // the SWIG-c++ interface, all the calculations are done at the c++ level,
    // then all the values are passed back to the Python level.  For large arrays,
    // the savings in function overhead should be quite large
    
    for (int i=0; i<n; i++)
    {
        OutVec[i]=Props(Output,Name1,Prop1[i],Name2,Prop2[i],Ref);
    }
}

double Props(char *Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
	// This is a wrapper function to allow for overloading, and the passing of a string rather than a single character
    if (strlen(Output)==1)
    {
        return Props(Output[0],Name1,Prop1,Name2,Prop2,Ref);
    }
}
double Props(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
    double T,p,Q,rhoV,rhoL,Value,rho,pdp,pbp;
    int isTwoPhase,success;
    char Local_errString[300];

    //Flush out any errors from the CoolProp error bubbling stack
	strcpy(CP_errString,"");
	ErrorFlag=OK;

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
    if (IsREFPROP(Ref))  // First eight characters match "REFPROP-"
    {
        Fluid.Type=FLUIDTYPE_REFPROP;
        #if defined(__ISWINDOWS__)
        return REFPROP(Output,Name1,Prop1,Name2,Prop2,Ref);
		#else
        sprintf(Local_errString,"Your refrigerant [%s] is from REFPROP, but REFPROP not supported on this platform",Ref);
        Append2ErrorString(Local_errString);
        if (FlagDebug==1)
        	PrintError();
        return -_HUGE;
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
        	sprintf(Local_errString,"For brine, Name1 must be 'T' and Name2 must be 'P'");
			Append2ErrorString(Local_errString);
			if (FlagDebug==1)
				PrintError();
			return -_HUGE;
        }
        Fluid.Type=FLUIDTYPE_BRINE;
        return SecFluids(Output,Prop1,Prop2,Ref);
    }
    else // It is something based on CoolProp routines
    {
        T=Prop1; 
        
        // Load the fluid-specific parameters and function pointers
        success=LoadFluid(Ref);
        //Error if bad fluid name
        if (ErrorFlag!=OK)
        {
        	if (FlagDebug==1)
				PrintError();
        	return -_HUGE;
        }
        
        // Check if it is an output that doesn't require a state input
        // Deal with it and return
        switch (Output)
        {
            case 'M':
                return Fluid.MM;
            case 'E':
                return Fluid.pc;
            case 'B':
                return Fluid.Tc;
            case 'R':
                return Fluid.Tt;
        }

        if (Name1=='T' && (Name2=='P' || Name2=='D' || Name2=='Q'))
        {
            // Temperature and Pressure are the inputs
            if (Name2=='P')
            {
            	// Collect the value for pressure
            	p=Prop2;
                if (FlagUseSinglePhaseLUT==1)
                {
                    return LookupValue(Output, T, p, Ref, &Fluid);
                }
                else
                {
                    //First find density as a function of temp and pressure (all parameters need it)
                    rho=get_Delta(T,p)*Fluid.rhoc;
					if (!ValidNumber(rho))
					{
						//ERROR
						printf("rho is definitely not a valid number [%g] with inputs of T: %g p: %g\n",rho,T,p);
						return _HUGE;
					}
                    else if (fabs(rho)>10000)
					{
						//ERROR
						printf("rho is probably not a valid value [%g] with inputs of T: %g p: %g\n",rho,T,p);
						return _HUGE;
					}
                    else if (Output=='D')
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

                if (T>Fluid.Tc || rho>1.002*rhosatL_func(T) || rho<0.998*rhosatV_func(T))
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
                    if (!ValidNumber(rhoL) || !ValidNumber(rhoV))
					{
						//ERROR
						sprintf(Local_errString,"rhoL [%g] or rhoV [%g] is invalid number\n",rhoL,rhoV);
						Append2ErrorString(Local_errString);
						if (FlagDebug==1)
							PrintError();
						return -_HUGE;
					}
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
							Value=SpeedSound_Trho(T,rho); break;
                        case 'V':
                            Value=Viscosity_Trho(T,rho); break;
                        case 'L':
                            Value=Conductivity_Trho(T,rho); break;
                        default:
                        	//ERROR
                            sprintf(Local_errString,"Invalid Output Name: %c",Output);
							Append2ErrorString(Local_errString);
							if (FlagDebug==1)
								PrintError();
							return -_HUGE;
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
                        if (ErrorFlag==FAIL)
                        {
                        if (FlagDebug==1)
							PrintError();
                        	return -_HUGE;
                        }
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
                        Value=Q*Enthalpy_Trho(T,rhoV)+(1-Q)*Enthalpy_Trho(T,rhoL);
                        break;
                    case 'D':
                        Value=rho;
                        break;
                    case 'S':
                        Value=Q*Entropy_Trho(T,rhoV)+(1-Q)*Entropy_Trho(T,rhoL);
                        break;
                    case 'U':
                        Value=Q*IntEnergy_Trho(T,rhoV)+(1-Q)*IntEnergy_Trho(T,rhoL);
                        break;
                    case 'C':
                        Value=Q*SpecHeatP_Trho(T,rhoV)+(1-Q)*SpecHeatP_Trho(T,rhoL);
                        break;
                    case 'O':
                        Value=Q*SpecHeatV_Trho(T,rhoV)+(1-Q)*SpecHeatV_Trho(T,rhoL);
                        break;
                    case 'V':
                        Value=Q*Viscosity_Trho(T,rhoV)+(1-Q)*Viscosity_Trho(T,rhoL);
                        break;
                    case 'L':
                        Value=Q*Conductivity_Trho(T,rhoV)+(1-Q)*Conductivity_Trho(T,rhoL);
                        break;
                    case 'A':
                        Value=Q*SpeedSound_Trho(T,rhoV)+(1-Q)*SpeedSound_Trho(T,rhoL);
                        break;
                    default:
                        sprintf(Local_errString,"Invalid Output Name: %c",Output);
						Append2ErrorString(Local_errString);
						if (FlagDebug==1)
							PrintError();
						return -_HUGE;
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
        	sprintf(Local_errString,"Names of input properties invalid (%c,%c) with refrigerant %s.  Valid choices are T,P or T,Q or T,D or P,Q",Name1,Name2,Ref);
			Append2ErrorString(Local_errString);
			if (FlagDebug==1)
				PrintError();
			return -_HUGE;
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
    else
    {
        return Props('B','T',0.0,'P',0.0,Ref);
    }
    //ERROR!
    return 0;
}
double Ttriple(char *Ref)
{
    // Set the global constants
    LoadFluid(Ref);
    
    // Brines do not have triple point temperatures, set it to a big number
    if (IsFluidType(Ref,"Brine"))
    {
        return 100000;
    }
    else
    {
        return Props('R','T',0.0,'Q',0.0,Ref);
    }
    //ERROR!
    return 0;
}

double T_hp(char *Ref, double h, double p, double T_guess)
{
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f,T=300;
    int iter=1;
    char Local_errString[300];
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
        if (iter>100)
        {
            //ERROR
        	sprintf(Local_errString,"iter %d: T_hp not converging with inputs(%s,%g,%g,%g) value: %0.12g\n",iter,Ref,h,p,T_guess,f);
			Append2ErrorString(Local_errString);
            printf("%s\n",Local_errString);
			return -_HUGE;
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
        return _HUGE;
    }
    #if defined(__ISWINDOWS__) 
        // It's a REFPROP fluid - use REFPROP to do all the calculations
        if (!strncmp(Ref,"REFPROP-",8))
        {
            return REFPROP('T','P',p,'Q',Q,Ref);
        }
    #endif
    
    // Do reverse interpolation in the Saturation Lookup Table
    if (FlagUseSaturationLUT && Fluid.Type==FLUIDTYPE_REFRIGERANT_PURE)
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
    
    if (Fluid.Type==FLUIDTYPE_REFRIGERANT_PURE)
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

    return _Dekker_Tsat(Tmin,Tmax,0.0001,p,Q,Ref);
    
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

double DerivTerms(char *Term,double T, double rho,char * Ref)
{
	LoadFluid(Ref);

	if (!strcmp(Term,"dpdT"))
	{
		double delta=rho/Fluid.rhoc;
		double tau=Fluid.Tc/T;
		double R=R_u/Fluid.MM;
		return rho*R*(1+delta*dphir_dDelta_func(tau,delta)-delta*tau*dphir2_dDelta_dTau_func(tau,delta));
	}
    else if (!strcmp(Term,"dvdp"))
	{
		double delta=rho/Fluid.rhoc;
		double tau=Fluid.Tc/T;
		double R=R_u/Fluid.MM;
        double dpdrho=R*T*(1+2*delta*dphir_dDelta_func(tau,delta)+delta*delta*dphir2_dDelta2_func(tau,delta));
		return -1/dpdrho/(rho*rho);
	}
	else
	{
		printf("Sorry DerivTerms is a work in progress, your derivative term [%s] is not available!!",Term);
		return _HUGE;
	}
    
//    // Wire up the pointers for the given refrigerant
//    if (!strcmp(Ref,"AA"))
//    {
//    }
//    else
//    {
//        printf("Bad Refrigerant Name in DerivTerms [%s]\n",Ref);
//        return _HUGE;
//    }
//
//    if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdT"))
//        return dhdT(Prop1,Prop2,TYPE_TPNoLookup);
//    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdrho"))
//        return dhdrho(Prop1,Prop2,TYPE_TPNoLookup);
//    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdT"))
//        return dpdT(Prop1,Prop2,TYPE_TPNoLookup);
//    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdrho"))
//        return dpdrho(Prop1,Prop2,TYPE_TPNoLookup);
//
//    else if (Name1=='T' && Name2=='D' && !strcmp(Term,"dhdT"))
//        return dhdT(Prop1,Prop2,TYPE_Trho);
//    else if (Name1=='T' && Name2=='D' && !strcmp(Term,"dhdrho"))
//        return dhdrho(Prop1,Prop2,TYPE_Trho);
//    else if (Name1=='T' && Name2=='D' && !strcmp(Term,"dpdT"))
//        return dpdT(Prop1,Prop2,TYPE_Trho);
//    else if (Name1=='T' && Name2=='D' && !strcmp(Term,"dpdrho"))
//        return dpdrho(Prop1,Prop2,TYPE_Trho);
//
//    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdTnum"))
//        return dhdT(Prop1,Prop2,99);
//    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdrhonum"))
//        return dhdrho(Prop1,Prop2,99);
//    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdTnum"))
//        return dpdT(Prop1,Prop2,99);
//    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdrhonum"))
//        return dpdrho(Prop1,Prop2,99);
//
//    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdTnum"))
//        return dhdT(Prop1,Prop2,99);
//    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdrhonum"))
//        return dhdrho(Prop1,Prop2,99);
//    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdTnum"))
//        return dpdT(Prop1,Prop2,99);
//    else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdrhonum"))
//        return dpdrho(Prop1,Prop2,99);
//    else
//    {
//        printf("Bad pair of properties[%c,%c] and derivative name [%s]\n",Name1,Name2,Term);
//        return _HUGE;
//    }
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

    Tc=Fluid.Tc;
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
