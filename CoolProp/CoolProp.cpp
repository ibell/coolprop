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
void _T_hp(std::string Ref, double h, double p, double *T, double *rho);
void _T_sp(std::string Ref, double s, double p, double *T, double *rho);
double rho_TP(double T, double p);
double _Props(std::string Output,char Name1, double Prop1, char Name2, double Prop2, std::string Ref);

bool FlagUseSaturationLUT=false; //Default to not use LUT
bool FlagUseSinglePhaseLUT=false; //Default to not use LUT

static std::string err_string;
int debug_level=0;
FluidsContainer Fluids = FluidsContainer();
Fluid * pFluid;

void swap (double *x, double *y)
{
	double tmp;
	tmp = *y;
	*y = *x;
	*x = tmp;
}
void swap (char *x, char *y)
{
	char tmp;
	tmp = *y;
	*y = *x;
	*x = tmp;
}

int get_debug(){return debug_level;}
int debug(){return debug_level;}
void debug(int level){debug_level=level;}

std::string get_errstring(void){return err_string;}
void get_errstring(char* str){str=(char*) get_errstring().c_str();};
char * get_errstringc(void){return (char*)err_string.c_str();}

int set_1phase_LUT_params(std::string Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax)
{ return set_1phase_LUT_params(Ref, nT, np, Tmin, Tmax, pmin, pmax, false); }

int set_1phase_LUT_params(char *Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax)
{ return set_1phase_LUT_params(Ref, nT, np, Tmin, Tmax, pmin, pmax, false); }

int set_1phase_LUT_params(char *Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax, bool rebuild)
{
	// Overload to take "normal" c-string and convert to std::string
	return set_1phase_LUT_params(std::string(Ref), nT, np, Tmin, Tmax, pmin, pmax,rebuild);
}
int set_1phase_LUT_params(std::string Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax, bool rebuild)
{
	try{
		// Get a pointer to the fluid (if possible)
		Fluid * LUTFluid = Fluids.get_fluid(Ref);
		// Set the LUT parameters
		LUTFluid->set_1phase_LUT_params(nT,np,Tmin,Tmax,pmin,pmax);
		// If you call this function, it will build the LUT on the next function call
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
void get_1phase_LUT_params(int *nT, int *np, double *Tmin, double *Tmax, double *pmin, double *pmax){
	pFluid->get_1phase_LUT_params(nT,np,Tmin,Tmax,pmin,pmax);
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
static int IsREFPROP(std::string Ref)
{
    if (!Ref.compare(0,8,"REFPROP-"))
        return 1;
    else
        return 0;
}
int IsFluidType(char *Ref, char *Type)
{
    if (IsBrine(Ref) && !strcmp(Type,"Brine"))
    {
        return 1;
    }
    else if (!pFluid->pure() && (!strcmp(Type,"PseudoPure") || !strcmp(Type,"PseudoPureFluid")))
    {
        return 1;
    }
	else if ((pFluid->pure() || IsREFPROP(Ref)) && !strcmp(Type,"PureFluid"))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
void Phase(char *Fluid,double T, double p, char *Phase_str)
{
	strcpy(Phase_str,(char*)Phase(std::string(Fluid),T,p).c_str());
}
std::string Phase(std::string Fluid, double T, double p)
{
	try{
		// Try to load the CoolProp Fluid
		pFluid = Fluids.get_fluid(Fluid);
		return pFluid->phase_Tp(T,p);
	}
	catch(NotImplementedError){
		return std::string("");
	}
}

// All the function interfaces that point to the single-input Props function
double Props1(char* Ref, char * Output)
{
	// Redirect to the Props function - should have called it Props1 from the outset
	return Props(Ref, Output);
}
double Props1(std::string Ref, std::string Output)
{
	// Redirect to the Props function - should have called it Props1 from the outset
	return Props((char*)Ref.c_str(), (char*)Output.c_str());
}
double Props(std::string Ref,std::string Output)
{
	return Props((char*)Ref.c_str(), (char*)Output.c_str());
}

double Props(char *Fluid, char *Output)
{
	// Fluid was loaded successfully
	try{	
		try{
			// Try to load the CoolProp Fluid
			pFluid = Fluids.get_fluid(Fluid);
		}
		catch(NotImplementedError){
			// It didn't load properly.  Perhaps it is a REFPROP fluid.
			try{
				if (!strcmp(Output,"Ttriple"))
					return Props('R','T',0,'P',0,Fluid);
				else if (!strcmp(Output,"Tcrit"))
					return Props('B','T',0,'P',0,Fluid);
				else if (!strcmp(Output,"pcrit"))
					return Props('E','T',0,'P',0,Fluid);
				else if (!strcmp(Output,"molemass"))
					return Props('M','T',0,'P',0,Fluid);
				else 
					throw ValueError(format("Output parameter \"%s\" is invalid",Output));
			}
			// Catch any error that subclasses the std::exception
			catch(std::exception &e){
				err_string = std::string("CoolProp error: ").append(e.what());
				std::cout << err_string <<std::endl;
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
		std::cout << err_string <<std::endl;
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
	// Go to the std::string, std::string version
	return Props(std::string(Output),Name1,Prop1,Name2,Prop2,std::string(Ref));
}
double Props(char Output,char Name1, double Prop1, char Name2, double Prop2, char* Ref)
{
	// Go to the std::string, std::string version
	return Props(std::string(1,Output),Name1,Prop1,Name2,Prop2,std::string(Ref));
}
double Props(std::string Output,char Name1, double Prop1, char Name2, double Prop2, std::string Ref)
{
	// In this function the error catching happens;
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
double _Props(std::string Output,char Name1, double Prop1, char Name2, double Prop2, std::string Ref)
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

	if (debug()>5){
		std::cout<<__FILE__<<": "<<Output<<","<<Name1<<","<<Prop1<<","<<Name2<<","<<Prop2<<","<<Ref<<std::endl;
		std::cout<<__FILE__<<": Using SinglePhase LUT is "<<FlagUseSinglePhaseLUT<<std::endl;
	}

    /* 
    If the fluid name is not actually a refrigerant name, but a string beginning with "REFPROP-",
    then REFPROP is used to calculate the desired property.
    */
    if (IsREFPROP(Ref))  // First eight characters match "REFPROP-"
    {
        #if defined(__ISWINDOWS__)
        return REFPROP(Output.c_str()[0],Name1,Prop1,Name2,Prop2,(char*)Ref.c_str());
		#else
        throw AttributeError(format("Your refrigerant [%s] is from REFPROP, but REFPROP not supported on this platform",Ref.c_str()));
        return -_HUGE;
        #endif
    }

    // **********************************************************************************
    // **********************************************************************************
    //                                Normal Property evaluation
    // **********************************************************************************
    // **********************************************************************************

    // It's a brine, call the brine routine
	else if (!Ref.compare("HC-10") || 
		!Ref.compare(0,2,"EG") || 
		!Ref.compare(0,2,"PG") || 
		!Ref.compare(0,8,"Methanol") || 
		!Ref.compare(0,7,"NH3/H2O"))
    {
        if (Name1!='T' || Name2!='P')
        {
			throw ValueError("For brine, Name1 must be 'T' and Name2 must be 'P'");
		}
        return SecFluids(Output[0],Prop1,Prop2,(char*)Ref.c_str());
    }
    else // It is something based on CoolProp routines
    {

		//Load the fluid - throws a NotImplementedError if not matched
		pFluid=Fluids.get_fluid(Ref);
        
        // Check if it is an output that doesn't require a state input
        // Deal with it and return
		
        if (Output[0]=='M')
			return pFluid->params.molemass;
		else if (Output[0]=='E')
			return pFluid->crit.p;
		else if (Output[0]=='B')
			return pFluid->crit.T;
		else if (Output[0]=='R')
			return pFluid->params.Ttriple;

		//Surface tension is only a function of temperature
		if (!Output.compare("I") || !Output.compare("SurfaceTension")){
			if (Name1=='T')
				return pFluid->surface_tension_T(Prop1);
			else
				throw ValueError(format("If output is surface tension ['I' or 'SurfaceTension'], Param1 must be temperature"));
		}

		// In any case, you want to get a (temperature, density) pair
		if ((Name1=='T' && Name2=='P') || (Name1=='P' && Name2 == 'T'))
		{
			//Swap values and keys
			if (Name1=='P' && Name2 == 'T')
			{
				swap(&Prop1,&Prop2);
				swap(&Name1,&Name2);
			}

			// Handle trivial outputs
			if (!Output.compare("P")) 
				return Prop2;
			else if (!Output.compare("T")) 
				return Prop1;
			
			// Get density as a function of T&p, then call again
			if (FlagUseSinglePhaseLUT==true){
				return pFluid->LookupValue_TP(std::string(Output), Prop1, Prop2);
            }
			else{
				rho = rho_TP(Prop1,Prop2);
			}
			
			if (!Output.compare("D")){
				if (debug()>5){
					std::cout<<__FILE__<<": "<<Output<<","<<Name1<<","<<Prop1<<","<<Name2<<","<<Prop2<<","<<Ref<<"="<<rho<<std::endl;
				}
				return rho;
			}
			else
				return Props(Output,Name1,Prop1,'D',rho,Ref);
		}
		else if ((Name1=='T' && Name2=='Q') || (Name1=='Q' && Name2=='T'))
		{
			if (Name1=='Q' && Name2=='T'){
				//Swap values and keys to get T,Q
				swap(&Prop1,&Prop2);
				swap(&Name1,&Name2);
			}
			Q=Prop2;
			// Get the saturation properties
			pFluid->saturation(Prop1,FlagUseSaturationLUT,&pL,&pV,&rhoL,&rhoV);
			// Find the effective density to use
			rho=1/(Q/rhoV+(1-Q)/rhoL);
			if (!Output.compare("P")) 
				return Q*pV+(1-Q)*pL;

			// Recurse and call Props again with the calculated density
			if (!Output.compare("D"))
				return 1/(Q/rhoV+(1-Q)/rhoL);
			else if (!Output.compare("C") || !Output.compare("O"))
				return Props(Output,'T',Prop1,'D',rho,Ref);
			else
				if (fabs(Q)<1e-12)
					return Props(Output,'T',Prop1,'D',rhoL,Ref);
				else if (fabs(Q-1)<1e-12)
					return Props(Output,'T',Prop1,'D',rhoV,Ref);
				else
					return Q*Props(Output,'T',Prop1,'D',rhoV,Ref)+(1-Q)*Props(Output,'T',Prop1,'D',rhoL,Ref);
		}
		else if ((Name1=='T' && Name2=='D') || (Name1=='D' && Name2=='T'))
		{
			if (Name1=='D' && Name2=='T')
			{
				//Swap values and keys to get T,D
				swap(&Prop1,&Prop2);
				swap(&Name1,&Name2);
			}
			T=Prop1;
			rho=Prop2;
			if (!Output.compare("D"))
				return Prop2;
			// If you are using LUT, use it
			if (FlagUseSinglePhaseLUT==1){
				// Try to use the LUT, if the parameter is not included in the LUT,
				// allow it to fall back to the conventional analysis
				try{
					Value = pFluid->LookupValue_Trho(std::string(Output), T, rho);
					return Value;
				}
                catch(ValueError){

                }
            }
			rho = Prop2; 
			if (!Output.compare("D"))
				Value=rho;
			else if (!Output.compare("T"))
				Value=T;
			else if (!Output.compare("P"))
				Value=pFluid->pressure_Trho(T,rho);
			else if (!Output.compare("H"))
				Value=pFluid->enthalpy_Trho(T,rho);
			else if (!Output.compare("S"))
				Value=pFluid->entropy_Trho(T,rho);
			else if (!Output.compare("U"))
				Value=pFluid->internal_energy_Trho(T,rho);
			else if (!Output.compare("C"))
				Value=pFluid->specific_heat_p_Trho(T,rho);
			else if (!Output.compare("C0"))
				Value=pFluid->specific_heat_p_ideal_Trho(T);
			else if (!Output.compare("O"))
				Value=pFluid->specific_heat_v_Trho(T,rho);
			else if (!Output.compare("A"))
				Value=pFluid->speed_sound_Trho(T,rho);
			else if (!Output.compare("G"))
				Value=pFluid->gibbs_Trho(T,rho);
			else if (!Output.compare("V"))
				Value=pFluid->viscosity_Trho(T,rho);
			else if (!Output.compare("L"))
				Value=pFluid->conductivity_Trho(T,rho);
			else if (!Output.compare("dpdT"))
				Value=pFluid->dpdT_Trho(T,rho);
			else{
				throw ValueError(format("Invalid Output Name: %s ",Output.c_str()));
				return _HUGE;
            }
			if (debug()>5){
				std::cout<<__FILE__<<": "<<Output<<","<<Name1<<","<<Prop1<<","<<Name2<<","<<Prop2<<","<<Ref<<"="<<Value<<std::endl;
			}
            return Value;
		}
        else if (Name1=='P' && Name2=='Q')
        {
            T=pFluid->Tsat(Prop1,Prop2,0);
            return Props(Output,'T',T,'Q',Prop2,Ref);
        }
		else if (Name1=='Q' && Name2=='P')
        {
            T=pFluid->Tsat(Prop2,Prop1,0);
            return Props(Output,'T',T,'Q',Prop1,Ref);
        }
        else if (Name1=='H' && Name2=='P')
        {
        	_T_hp(Ref,Prop1,Prop2,&T, &rho);
			return Props(Output,'T',T,'D',rho,Ref);
        }
		else if (Name1=='P' && Name2=='H')
        {
        	_T_hp(Ref,Prop2,Prop1,&T, &rho);
			return Props(Output,'T',T,'D',rho,Ref);
        }
		else if (Name1=='S' && Name2=='P')
        {
        	_T_sp(Ref,Prop1,Prop2,&T, &rho);
			return Props(Output,'T',T,'D',rho,Ref);
        }
		else if (Name1=='P' && Name2=='S')
        {
        	_T_sp(Ref,Prop2,Prop1,&T, &rho);
			return Props(Output,'T',T,'D',rho,Ref);
        }
        else
        {
			throw ValueError(format("Not a valid pair of input keys %c,%c and output key %s",Name1,Name2,Output.c_str()));
        }
    }
    return 0;
}
double rho_TP(double T, double p)
{
	// Calculate the density as a function of T&p, either using EOS or LUT
	if (FlagUseSinglePhaseLUT==true)
    {
		return pFluid->LookupValue_TP(std::string("D"), T, p);
    }
    else
    {
        //Find density as a function of temp and pressure (all parameters need it)
		return pFluid->density_Tp(T,p);
    }
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

void _T_sp(std::string Ref, double s, double p, double *Tout, double *rhoout)
{
	int iter;
	double A[2][2], B[2][2],T_guess,R;
	double dar_ddelta,da0_dtau,d2a0_dtau2,dar_dtau,d2ar_ddelta_dtau,d2ar_ddelta2,d2ar_dtau2,d2a0_ddelta_dtau,ar,a0,da0_ddelta;
	double f1,f2,df1_dtau,df1_ddelta,df2_ddelta,df2_dtau,s_hot;
    double rhosatL,ssatL,ssatV,TsatL,TsatV,tau,delta,worst_error;
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
		rhosatL=Props(std::string("D"),'P',p,'Q',0.0,Ref);
		ssatL=Props(std::string("S"),'P',p,'Q',0.0,Ref);
		ssatV=Props(std::string("S"),'P',p,'Q',1.0,Ref);
		TsatL=Props(std::string("T"),'P',p,'Q',0.0,Ref);
		TsatV=Props(std::string("T"),'P',p,'Q',1.0,Ref);

		if (s>ssatV)
		{
			// Superheated vapor
			// Very superheated
			s_hot = Props(std::string("S"),'T',TsatV+40.0,'P',p,Ref);
			T_guess = TsatV+(s-ssatV)/(s_hot-ssatV)*40.0;
			delta = Props(std::string("D"),'T',T_guess,'P',p,Ref)/ (pFluid->reduce.rho);
		}
		else if (s<ssatL)
		{
			// Subcooled liquid
			T_guess = TsatL*log((s-ssatL)/Props(std::string("C"),'P',p,'Q',0.0,Ref));
			delta = rhosatL/pFluid->reduce.rho;
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
		a0 = pFluid->phi0(tau,delta);
    	da0_dtau = pFluid->dphi0_dTau(tau,delta);
    	d2a0_dtau2 = pFluid->d2phi0_dTau2(tau,delta);
    	da0_ddelta = pFluid->dphi0_dDelta(tau,delta);
		d2a0_ddelta_dtau = 0.0;

    	ar = pFluid->phir(tau,delta);
		dar_dtau = pFluid->dphir_dTau(tau,delta);
		d2ar_dtau2 = pFluid->d2phir_dTau2(tau,delta);
		dar_ddelta = pFluid->dphir_dDelta(tau,delta);
		d2ar_ddelta2 = pFluid->d2phir_dDelta2(tau,delta);
		d2ar_ddelta_dtau = pFluid->d2phir_dDelta_dTau(tau,delta);
		
		// Residual and derivatives thereof for entropy
		f1 = tau*(da0_dtau+dar_dtau)-ar-a0-s/R;
		df1_dtau = tau*(d2a0_dtau2 + d2ar_dtau2)+(da0_dtau+dar_dtau)-dar_dtau-da0_dtau;
		df1_ddelta = tau*(d2a0_ddelta_dtau+d2ar_ddelta_dtau)-dar_ddelta-da0_ddelta;

		// Residual and derivatives thereof for pressure
		f2 = delta/tau*(1+delta*dar_ddelta)-p/(rhoc*R*Tc);
		df2_dtau = (1+delta*dar_ddelta)*(-delta/tau/tau)+delta/tau*(delta*d2ar_ddelta_dtau);
		df2_ddelta = (1.0/tau)*(1+2.0*delta*dar_ddelta+delta*delta*d2ar_ddelta2);

		//First index is the row, second index is the column
		A[0][0]=df1_dtau;
		A[0][1]=df1_ddelta;
		A[1][0]=df2_dtau;
		A[1][1]=df2_ddelta;

		MatInv_2(A,B);
		tau -= B[0][0]*f1+B[0][1]*f2;
		delta -= B[1][0]*f1+B[1][1]*f2;

        if (fabs(f1)>fabs(f2))
            worst_error=fabs(f1);
        else
            worst_error=fabs(f2);

		iter+=1;
		if (iter>100)
		{
			printf("_Tsp did not converge\n");
			*Tout = _HUGE;
            *rhoout = _HUGE;
		}
    }
    *Tout = Tc/tau;
    *rhoout = delta*rhoc;
}

void _T_hp(std::string Ref, double h, double p, double *Tout, double *rhoout)
{
	int iter;
	double A[2][2], B[2][2],T_guess,R;
	double dar_ddelta,da0_dtau,d2a0_dtau2,dar_dtau,d2ar_ddelta_dtau,d2ar_ddelta2,d2ar_dtau2,d2a0_ddelta_dtau;
	double f1,f2,df1_dtau,df1_ddelta,df2_ddelta,df2_dtau,h_hot;
    double rhosatL,hsatL,hsatV,TsatL,TsatV,tau,delta,worst_error;
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
		rhosatL=Props(std::string("D"),'P',p,'Q',0.0,Ref);
		hsatL=Props(std::string("H"),'P',p,'Q',0.0,Ref);
		hsatV=Props(std::string("H"),'P',p,'Q',1.0,Ref);
		TsatL=Props(std::string("T"),'P',p,'Q',0.0,Ref);
		TsatV=Props(std::string("T"),'P',p,'Q',1.0,Ref);

		if (h>hsatV)
		{
			//Superheated vapor
			// Very superheated
			h_hot = Props(std::string("H"),'T',TsatV+40.0,'P',p,Ref);
			T_guess = TsatV+(h-hsatV)/(h_hot-hsatV)*40.0;
			delta = Props(std::string("D"),'T',T_guess,'P',p,Ref)/ (pFluid->reduce.rho);
		}
		else if (h<hsatL)
		{
			// Subcooled liquid
			T_guess = TsatL+(h-hsatL)/Props(std::string("C"),'P',p,'Q',0.0,Ref);
			delta = rhosatL/pFluid->reduce.rho;
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

        if (fabs(f1)>fabs(f2))
            worst_error=fabs(f1);
        else
            worst_error=fabs(f2);

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
			//throw SolutionError(format("iter %d: T_hp not converging with inputs(%s,%g,%g,%g) value: %0.12g\n",iter,Ref,h,p,T_guess,f));
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

void FluidsList(char* str)
{
	str=(char*)FluidsList().c_str();
	return;
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
		// not sure this is working properly, so not documented
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
	else if (!strcmp(Term,"d2phir_dDelta_dTau"))
    {
        return pFluid->d2phir_dDelta_dTau(tau,delta);
    }
	else if (!strcmp(Term,"d3phir_dDelta2_dTau"))
    {
        return pFluid->d3phir_dDelta2_dTau(tau,delta);
    }
	else if (!strcmp(Term,"phi0"))
    {
        return pFluid->phi0(tau,delta);
    }
    else if (!strcmp(Term,"dphi0_dTau"))
    {
        return pFluid->dphi0_dTau(tau,delta);
    }
	else if (!strcmp(Term,"d2phi0_dTau2"))
    {
        return pFluid->d2phi0_dTau2(tau,delta);
    }
	else if (!strcmp(Term,"dphi0_dDelta"))
    {
        return pFluid->dphi0_dDelta(tau,delta);
    }
	else if (!strcmp(Term,"d2phi0_dDelta2"))
    {
        return pFluid->d2phi0_dDelta2(tau,delta);
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
std::string get_EOSReference(std::string Ref)
{
	pFluid=Fluids.get_fluid(Ref);
	if (pFluid!=NULL)
	{
		return pFluid->get_EOSReference();
	}
	else
		return std::string("");
}
std::string get_TransportReference(std::string Ref)
{
	pFluid=Fluids.get_fluid(Ref);
	if (pFluid!=NULL)
	{
		return pFluid->get_TransportReference();
	}
	else
		return std::string("");
}
void get_REFPROPname(char* Ref, char * str)
{
	str= (char*)get_REFPROPname(std::string(Ref)).c_str();
}
std::string get_REFPROPname(std::string Ref)
{
	pFluid=Fluids.get_fluid(Ref);

	if ( (pFluid->get_REFPROPname()).size()!=0 ){
		return pFluid->get_REFPROPname();
	}
	else{
		return pFluid->get_name();
	}
}
