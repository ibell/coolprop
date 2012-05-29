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
#include <list>
#include <exception>
#include <stdio.h>
#include "string.h"
#include "FluidClass.h"
#include "CoolPropTools.h"
#include "CPExceptions.h"
#include "Brine.h"


// Function prototypes
void _T_hp(char *Ref, double h, double p, double *T, double *rho);
double rho_TP(double T, double p);
double _Props(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref);

bool FlagUseSaturationLUT=false; //Default to use LUT since they are used so much and don't take long to build
bool FlagUseSinglePhaseLUT=false; //Default to not use LUT

std::string err_string;

FluidsContainer Fluids = FluidsContainer();
Fluid * pFluid;

std::string get_errstring(void){return err_string;}
char * get_errstringc(void){return (char*)err_string.c_str();}

int set_1phase_LUT_params(std::string Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax)
{ return set_1phase_LUT_params(Ref, nT, np, Tmin, Tmax, pmin, pmax, false); }

int set_1phase_LUT_params(char *Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax)
{ return set_1phase_LUT_params(Ref, nT, np, Tmin, Tmax, pmin, pmax, false); }

int set_1phase_LUT_params(std::string Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax, bool rebuild)
{
	try{
		// Get a pointer to the fluid (if possible)
		Fluid * LUTFluid = Fluids.get_fluid(Ref);
		// Set the LUT parameters
		LUTFluid->set_1phase_LUT_params(nT,np,Tmin,Tmax,pmin,pmax);

		if (rebuild==true)
			LUTFluid->BuildLookupTable();
	}
	catch (NotImplementedError &e)
	{
		err_string = std::string("CoolProp Error: ").append(e.what());
		return -1;
	}
	return 0;
}

int set_1phase_LUT_params(char *Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax, bool rebuild)
{
	// Overload to take "normal" c-string and convert to std::string
	return set_1phase_LUT_params(std::string(Ref), nT, np, Tmin, Tmax, pmin, pmax,rebuild);
}

void UseSaturationLUT(bool OnOff)
{
    if (OnOff==true || OnOff==false)
    {
        FlagUseSaturationLUT=(bool)OnOff;
    }
    else
    {
        printf("Sorry, UseSaturationLUT() takes an integer input, either 0 (no) or 1 (yes)\n");
    }
}
bool SaturationLUTStatus()
{
	return FlagUseSaturationLUT;
}

void UseSinglePhaseLUT(bool OnOff)
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
bool SinglePhaseLUTStatus(void)
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
	else if ((pFluid->pure() || IsREFPROP(Ref)) && !strcmp(Type,"PureFluidStruct"))
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
    if (p>Props(Ref,"pcrit") && T>Props(Ref,"Tcrit"))
    {
        return PHASE_SUPERCRITICAL;
    }
    else
    {
        Tbubble=Props('T', 'P',p, 'Q',0.0,Ref);
        Tdew=Props('T','P', p, 'Q',1.0,Ref);
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

double Props(std::string Ref,std::string Output)
{
	return Props((char*)Ref.c_str(),(char*)Output.c_str());
}

double Props(char *Fluid, char *Output)
{
	// Fluid was loaded successfully
	try{	
		try{
			// Try to load the CoolProp Fluid
			pFluid = Fluids.get_fluid(Fluid);
		}
		catch(NotImplementedError &e){
			// It didn't load properly.  Perhaps it is a REFPROP fluid.
			try{
				if (!strcmp(Output,"Ttriple"))
					return _Props('R','T',0,'P',0,Fluid);
				else if (!strcmp(Output,"Tcrit"))
					return _Props('B','T',0,'P',0,Fluid);
				else if (!strcmp(Output,"pcrit"))
					return _Props('E','T',0,'P',0,Fluid);
				else if (!strcmp(Output,"molemass"))
					return _Props('M','T',0,'P',0,Fluid);
				else 
					throw ValueError(format("Output parameter \"%s\" is invalid",Output));
			}
			// Catch any error that subclasses the std::exception
			catch(std::exception &e){
				err_string = std::string("CoolProp error: ").append(e.what());
				return _HUGE;
			}
		}

		if (!strcmp(Output,"Ttriple"))
			return pFluid->params.Ttriple;
		else if (!strcmp(Output,"Tcrit"))
			return pFluid->crit.T;
		else if (!strcmp(Output,"pcrit"))
			return pFluid->crit.p;
		else if (!strcmp(Output,"rhocrit"))
			return pFluid->crit.rho;
		else if (!strcmp(Output,"molemass"))
			return pFluid->params.molemass;
		else
		{
			throw ValueError(format("Output parameter \"%s\" is invalid",Output));
			return _HUGE;
		}
	}
	// Catch any error that subclasses the std::exception
	catch(std::exception &e){
		err_string = std::string("CoolProp error: ").append(e.what());
		return _HUGE;
	}
	catch(...)
	{
		std::cout << "Indeterminate error" << std::endl;
		return _HUGE;
	}
}

double Props(char *Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
	// This is a wrapper function to allow for overloading, and the passing of a string rather than a single character
    if (strlen(Output)==1)
    {
		try
		{
			return _Props(Output[0],Name1,Prop1,Name2,Prop2,Ref);
		}
		// Catch any error that subclasses the std::exception
		catch(std::exception &e){
			err_string = std::string("CoolProp error: ").append(e.what());
			return _HUGE;
		}
		catch(...){
			std::cout << "Indeterminate error:" << std::endl;
			return _HUGE;
		}
    }
	else
	{
		printf("Currently Props with a string first input is not supported");
		return _HUGE;
	}
}
double Props(char Output,char Name1, double Prop1, char Name2, double Prop2, std::string Ref)
{
	return Props(Output,Name1,Prop1,Name2,Prop2,(char*)Ref.c_str());
}

double Props(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
	try{
		return _Props(Output,Name1,Prop1,Name2,Prop2,Ref);
	}
	catch(std::exception &e){
		err_string=std::string("CoolProp error: ").append(e.what());
		return _HUGE;
	}
	catch(...){
		std::cout << "Indeterminate error" << std::endl;
		return _HUGE;
	}
}
// Make this a wrapped function so that error bubbling can be done properly
double _Props(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
    double T,Q,rhoV,rhoL,Value,rho,pL,pV;

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
			throw ValueError("For brine, Name1 must be 'T' and Name2 must be 'P'");
		}
        return SecFluids(Output,Prop1,Prop2,Ref);
    }
    else // It is something based on CoolProp routines
    {

		//Load the fluid - throws a NotImplementedError if not matched
		pFluid=Fluids.get_fluid(Ref);

        T=Prop1; 
        
        // Check if it is an output that doesn't require a state input
        // Deal with it and return
        switch (Output)
        {
            case 'M':
				return pFluid->params.molemass;
            case 'E':
				return pFluid->crit.p;
            case 'B':
				return pFluid->crit.T;
            case 'R':
				return pFluid->params.Ttriple;
        }

		// In any case, you want to get a (temperature,density) pair
		if (Name1=='T' && Name2=='P')
		{
			if (Output=='P') return Prop2;
			else if (Output=='T') return Prop1;
			// Get density as a function of T&p
			if (FlagUseSinglePhaseLUT==true){
                return pFluid->LookupValue_TP(Output, T, Prop2);
            }
			else{
				rho = rho_TP(Prop1,Prop2);
			}
			
			if (Output=='D')
				return rho;
			else
				return Props(Output,Name1,Prop1,'D',rho,Ref);
		}
		else if (Name1=='T' && Name2=='Q')
		{
			Q=Prop2;
			// Get the saturation properties
			pFluid->saturation(Prop1,FlagUseSaturationLUT,&pL,&pV,&rhoL,&rhoV);
			// Find the effective density to use
			rho=1/(Q/rhoV+(1-Q)/rhoL);
			if (Output=='P') 
				return Q*pV+(1-Q)*pL;

			// Recurse and call Props again with the calculated density
			if (Output=='D')
				return 1/(Q/rhoV+(1-Q)/rhoL);
			else if (Output=='C' || Output=='O')
				return Props(Output,'T',Prop1,'D',rho,Ref);
			else
				return Q*Props(Output,'T',Prop1,'D',rhoV,Ref)+(1-Q)*Props(Output,'T',Prop1,'D',rhoL,Ref);
		}
		else if (Name1=='T' && Name2=='D')
		{
			T=Prop1;
			rho=Prop2;
			if (Output=='D')
				return Prop2;
			// If you are using LUT, use it
			if (FlagUseSinglePhaseLUT==1){
                Value = pFluid->LookupValue_Trho(Output, T, rho);
                return Value;
            }
			rho = Prop2;
            switch (Output)
            {
				case 'D':
					Value=rho; break;
                case 'P':
                    Value=pFluid->pressure_Trho(T,rho); break;
                case 'H':
                    Value=pFluid->enthalpy_Trho(T,rho); break;
                case 'S':
                    Value=pFluid->entropy_Trho(T,rho); break;
                case 'U':
                    Value=pFluid->internal_energy_Trho(T,rho); break;
                case 'C':
                    Value=pFluid->specific_heat_p_Trho(T,rho); break;
                case 'O':
                    Value=pFluid->specific_heat_v_Trho(T,rho); break;
                case 'A':
					Value=pFluid->speed_sound_Trho(T,rho); break;
				case 'G':
					Value=pFluid->gibbs_Trho(T,rho); break;
                case 'V':
                    Value=pFluid->viscosity_Trho(T,rho); break;
                case 'L':
                    Value=pFluid->conductivity_Trho(T,rho); break;
                default:
					throw ValueError(format("Invalid Output Name: %c",Output));
					return _HUGE;
            }
            return Value;
		}
        else if (Output=='T' && Name1=='P' && Name2=='Q')
        {
			return pFluid->Tsat(Prop1,Prop2,0);
        }
        else if (Name1=='P' && Name2=='Q')
        {
            T=pFluid->Tsat(Prop1,Prop2,0);
            return Props(Output,'T',T,'Q',Prop2,Ref);
        }
        else if (Output=='T' && Name1=='H' && Name2=='P')
        {
        	_T_hp(Ref,Prop1,Prop2,&T, &rho);
            return T;
        }
        else
        {
			throw ValueError(format("Not a valid pair of input keys %c,%c and output key %c",Name1,Name2,Output));
        }
    }
    return 0;
}
double rho_TP(double T, double p)
{
	// Calculate the density as a function of T&p, either using EOS or LUT
	if (FlagUseSinglePhaseLUT==true)
    {
		return pFluid->LookupValue_TP('D', T, p);
    }
    else
    {
        //Find density as a function of temp and pressure (all parameters need it)
		return pFluid->density_Tp(T,p);
    }
}
void _T_hp(char *Ref, double h, double p, double *Tout, double *rhoout)
{
	int iter;
	double A[2][2], B[2][2],T_guess,R;
	double dar_ddelta,da0_dtau,d2a0_dtau2,dar_dtau,d2ar_ddelta_dtau,d2ar_ddelta2,d2ar_dtau2,d2a0_ddelta_dtau;
	double f1,f2,df1_dtau,df1_ddelta,df2_ddelta,df2_dtau,h_hot;
    double hsatL,hsatV,TsatL,TsatV,tau,delta,worst_error;
	//First figure out where you are
	
	R=pFluid->R();
	if (p > Props(Ref,"pcrit"))
	{
        //Supercritical pressure
		T_guess = Props(Ref,"Tcrit")+30.0;
		delta = p/(R*T_guess)/pFluid->reduce.rho/0.7;
	}
	else
	{
		hsatL=Props('H','P',p,'Q',0.0,Ref);
		hsatV=Props('H','P',p,'Q',1.0,Ref);
		TsatL=Props('T','P',p,'Q',0.0,Ref);
		TsatV=Props('T','P',p,'Q',1.0,Ref);

		if (h>hsatV)
		{
			//Superheated vapor
			// Very superheated
			h_hot = Props('H','T',TsatV+40.0,'P',p,Ref);
			T_guess = TsatV+(h-hsatV)/(h_hot-hsatV)*40.0;
			delta = Props('D','T',T_guess,'P',p,Ref)/ (pFluid->reduce.rho);
		}
		else if (h<hsatL)
		{
			// Subcooled liquid
			T_guess = TsatL+(h-hsatL)/Props('C','P',p,'Q',0.0,Ref);
			delta = 10;
		}
		else
		{

		}
	}

	double Tc = pFluid->reduce.T;
	tau=Tc/T_guess;
	double rhoc = pFluid->reduce.rho;

    worst_error=999;
    iter=0;
    while (worst_error>1e-6)
    {
    	// All the required partial derivatives

    	da0_dtau = pFluid->dphi0_dTau(tau,delta);
    	d2a0_dtau2 = pFluid->d2phi0_dTau2(tau,delta);
    	d2a0_ddelta_dtau = 0.0;
    	dar_dtau = pFluid->dphir_dTau(tau,delta);
		dar_ddelta = pFluid->dphir_dDelta(tau,delta);
		d2ar_ddelta_dtau = pFluid->d2phir_dDelta_dTau(tau,delta);
		d2ar_ddelta2 = pFluid->d2phir_dDelta2(tau,delta);
		d2ar_dtau2 = pFluid->d2phir_dTau2(tau,delta);

		f1 = delta/tau*(1+delta*dar_ddelta)-p/(rhoc*R*Tc);
		f2 = (1+delta*dar_ddelta)+tau*(da0_dtau+dar_dtau)-tau*h/(R*Tc);
		df1_dtau = (1+delta*dar_ddelta)*(-delta/tau/tau)+delta/tau*(delta*d2ar_ddelta_dtau);
		df1_ddelta = (1.0/tau)*(1+2.0*delta*dar_ddelta+delta*delta*d2ar_ddelta2);
		df2_dtau = delta*d2ar_ddelta_dtau+da0_dtau+dar_dtau+tau*(d2a0_dtau2+d2ar_dtau2)-h/(R*Tc);
		df2_ddelta = (dar_ddelta+delta*d2ar_ddelta2)+tau*(d2a0_ddelta_dtau+d2ar_ddelta_dtau);

		//First index is the row, second index is the column
		A[0][0]=df1_dtau;
		A[0][1]=df1_ddelta;
		A[1][0]=df2_dtau;
		A[1][1]=df2_ddelta;

		MatInv_2(A,B);
		tau -= B[0][0]*f1+B[0][1]*f2;
		delta -= B[1][0]*f1+B[1][1]*f2;

		worst_error = max(fabs(f1),fabs(f2));
		iter+=1;
		if (iter>100)
		{
			printf("_Thp did not converge\n");
			*Tout = _HUGE;
            *rhoout = _HUGE;
		}
    }
    *Tout = Tc/tau;
    *rhoout = delta*rhoc;
}

double T_hp(char *Ref, double h, double p, double T_guess)
{
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f=999,T=300;
    int iter=1;

    while ((iter<=3 || fabs(f)>eps) && iter<100)
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
			throw SolutionError(format("iter %d: T_hp not converging with inputs(%s,%g,%g,%g) value: %0.12g\n",iter,Ref,h,p,T_guess,f));
			return _HUGE;
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

double K2F(double T)
{ return T * 9 / 5 - 459.67; }

double F2K(double T_F)
{ return (T_F + 459.67) * 5 / 9;}

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

std::string FluidsList()
{
	return Fluids.FluidList();
}

double DerivTerms(char *Term,double T, double rho, char * Ref)
{
	pFluid=Fluids.get_fluid(Ref);

    double rhoc =pFluid->reduce.rho;
	double delta=rho/rhoc;
	double tau=pFluid->reduce.T/T;
	double R=pFluid->R();
	double dtau_dT = -pFluid->reduce.T/T/T;

	if (!strcmp(Term,"dpdT"))
	{
		return rho*R*(1+delta*pFluid->dphir_dDelta(tau,delta)-delta*tau*pFluid->d2phir_dDelta_dTau(tau,delta));
	}
    else if (!strcmp(Term,"dpdrho"))
	{
        double dpdrho=R*T*(1+2*delta*pFluid->dphir_dDelta(tau,delta)+delta*delta*pFluid->d2phir_dDelta2(tau,delta));
		return dpdrho;
	}
    else if (!strcmp(Term,"dvdp"))
	{
        double dpdrho=R*T*(1+2*delta*pFluid->dphir_dDelta(tau,delta)+delta*delta*pFluid->d2phir_dDelta2(tau,delta));
		return -1/dpdrho/(rho*rho);
	}
	else if (!strcmp(Term,"Z"))
	{
		return 1+delta*pFluid->dphir_dDelta(tau,delta);
	}
	else if (!strcmp(Term,"dZ_dDelta"))
    {
        return delta*pFluid->d2phir_dDelta2(tau,delta)+pFluid->dphir_dDelta(tau,delta);
    }
	else if (!strcmp(Term,"dZ_dTau"))
    {
        return delta*pFluid->d2phir_dDelta_dTau(tau,delta);
    }
	else if (!strcmp(Term,"B"))
	{
		// given by B*rhoc=lim(delta --> 0) [dphir_ddelta(tau)]
		return 1.0/rhoc*pFluid->dphir_dDelta(tau,1e-12);
	}
	else if (!strcmp(Term,"dBdT"))
	{
		return 1.0/rhoc*pFluid->d2phir_dDelta_dTau(tau,1e-12)*dtau_dT;
	}
	else if (!strcmp(Term,"C"))
	{
		// given by C*rhoc^2=lim(delta --> 0) [d2phir_dDelta2(tau)]
		return 1.0/(rhoc*rhoc)*pFluid->d2phir_dDelta2(tau,1e-12);
	}
	else if (!strcmp(Term,"dCdT"))
	{
		return 1.0/(rhoc*rhoc)*pFluid->d3phir_dDelta2_dTau(tau,1e-12)*dtau_dT;
    }
	else if (!strcmp(Term,"phir"))
    {
        return pFluid->phir(tau,delta);
    }
	else if (!strcmp(Term,"phi0"))
    {
        return pFluid->phi0(tau,delta);
    }
    else if (!strcmp(Term,"dphi0_dTau"))
    {
        return pFluid->dphi0_dTau(tau,delta);
    }
	else if (!strcmp(Term,"dphir_dTau"))
    {
        return pFluid->dphir_dTau(tau,delta);
    }
	else if (!strcmp(Term,"dphir_dDelta"))
    {
        return pFluid->dphir_dDelta(tau,delta);
    }
	else if (!strcmp(Term,"d2phir_dTau2"))
    {
        return pFluid->d2phir_dTau2(tau,delta);
    }
	else if (!strcmp(Term,"d2phir_dDelta2"))
    {
        return pFluid->d2phir_dDelta2(tau,delta);
    }
	
	else if (!strcmp(Term,"IsothermalCompressibility"))
	{
		double dpdrho=R*T*(1+2*delta*pFluid->dphir_dDelta(tau,delta)+delta*delta*pFluid->d2phir_dDelta2(tau,delta));
		return 1/(rho*dpdrho);
	}
	else
	{
		printf("Sorry DerivTerms is a work in progress, your derivative term [%s] is not available!!",Term);
		return _HUGE;
	}
   
}


