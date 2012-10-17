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
double _Props(std::string Output,std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref);
double _CoolProp_Fluid_Props(long iOutput, long iName1, double Value1, long iName2, double Value2, Fluid *pFluid);

bool FlagUseSaturationLUT=false; //Default to not use LUT
bool FlagUseSinglePhaseLUT=false; //Default to not use LUT

static std::string err_string;
int debug_level=0;
FluidsContainer Fluids = FluidsContainer();
Fluid * pFluid;

// Define some constants that will be used throughout
enum params {iB,iT,iP,iD,iC,iC0,iO,iU,iH,iS,iA,iG,iQ,iV,iL,iI,iMM,iTcrit,iTtriple,iPcrit,iRhocrit,iAccentric,iDpdT};

// This is a map of all possible strings to a unique identifier
std::pair<std::string, long> map_data[] = {
	std::make_pair(std::string("E"),iPcrit),
	std::make_pair(std::string("M"),iMM),
	std::make_pair(std::string("w"),iAccentric),
	std::make_pair(std::string("R"),iTtriple),
	std::make_pair(std::string("N"),iRhocrit),
	std::make_pair(std::string("B"),iTcrit),

	std::make_pair(std::string("pcrit"),iPcrit),
	std::make_pair(std::string("molemass"),iMM),
	std::make_pair(std::string("accentric"),iAccentric),
	std::make_pair(std::string("Ttriple"),iTtriple),
	std::make_pair(std::string("rhocrit"),iRhocrit),
	std::make_pair(std::string("Tcrit"),iTcrit),

	std::make_pair(std::string("Q"),iQ),
	std::make_pair(std::string("T"),iT),
    std::make_pair(std::string("P"),iP),
	std::make_pair(std::string("D"),iD),
	std::make_pair(std::string("C"),iC),
	std::make_pair(std::string("C0"),iC0),
	std::make_pair(std::string("O"),iO),
	std::make_pair(std::string("U"),iU),
	std::make_pair(std::string("H"),iH),
	std::make_pair(std::string("S"),iS),
	std::make_pair(std::string("A"),iA),
	std::make_pair(std::string("G"),iG),
	std::make_pair(std::string("V"),iV),
	std::make_pair(std::string("L"),iL),
	std::make_pair(std::string("I"),iI),
	std::make_pair(std::string("SurfaceTension"),iI),
	std::make_pair(std::string("dpdT"),iDpdT),
	std::make_pair(std::string("M"),iMM),
};
//Now actually construct the map
std::map<std::string, long> param_map(map_data,
    map_data + sizeof map_data / sizeof map_data[0]);

// This is a map of all unique identifiers to std::strings with the default units
// that are used internally
std::pair<long, std::string> units_data[] = {
	std::make_pair(iPcrit, std::string("kPa")),
	std::make_pair(iMM, std::string("kg/kmol")),
	std::make_pair(iAccentric, std::string("-")),
	std::make_pair(iTtriple, std::string("K")),
	std::make_pair(iRhocrit, std::string("kg/m^3")),
	std::make_pair(iTcrit, std::string("K")),

	std::make_pair(iQ, std::string("-")),
	std::make_pair(iT, std::string("K")),
    std::make_pair(iP, std::string("kPa")),
	std::make_pair(iD, std::string("kg/m^3")),
	std::make_pair(iC, std::string("kJ/kg/K")),
	std::make_pair(iC0, std::string("kJ/kg/K")),
	std::make_pair(iO, std::string("kJ/kg/K")),
	std::make_pair(iU, std::string("kJ/kg")),
	std::make_pair(iH, std::string("kJ/kg")),
	std::make_pair(iS, std::string("kJ/kg/K")),
	std::make_pair(iA, std::string("m/s")),
	std::make_pair(iG, std::string("kJ/kg")),
	std::make_pair(iV, std::string("Pa*s")),
	std::make_pair(iL, std::string("kW/m/K")),
	std::make_pair(iI, std::string("N/m")),
	std::make_pair(iDpdT, std::string("kPa/K"))
};

//Now actually construct the map
std::map<long, std::string> units_map(units_data,
		units_data + sizeof units_data / sizeof units_data[0]);

double _T_hp_secant(std::string Ref, double h, double p, double T_guess)
{
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f=999,T=300;
    int iter=1;

    while ((iter<=3 || fabs(f)>eps) && iter<100)
    {
        if (iter==1){x1=T_guess; T=x1;}
        if (iter==2){x2=T_guess+0.1; T=x2;}
        if (iter>2) {T=x2;}
            f=Props("H",'T',T,'P',p,Ref)-h;
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
        }
    }
    return T;
}

void swap (double *x, double *y)
{
	double tmp;
	tmp = *y;
	*y = *x;
	*x = tmp;
}
void swap (std::string *x, std::string *y)
{
	std::string tmp;
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
void swap (long *x, long *y)
{
	long tmp;
	tmp = *y;
	*y = *x;
	*x = tmp;
}


EXPORT_CODE int CONVENTION get_debug(){return debug_level;}
int  debug(){return debug_level;}
EXPORT_CODE void CONVENTION debug(int level){debug_level=level;}

std::string get_errstring(void){return err_string;}
EXPORT_CODE void CONVENTION get_errstring(char* str){str=(char*) get_errstring().c_str();};
EXPORT_CODE char * CONVENTION get_errstringc(void){return (char*)err_string.c_str();}

int set_1phase_LUT_params(std::string Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax)
{ return set_1phase_LUT_params(Ref, nT, np, Tmin, Tmax, pmin, pmax, false); }

int set_1phase_LUT_params(char *Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax)
{ return set_1phase_LUT_params(Ref, nT, np, Tmin, Tmax, pmin, pmax, false); }

EXPORT_CODE int CONVENTION set_1phase_LUT_params(char *Ref, int nT, int np, double Tmin, double Tmax, double pmin, double pmax, bool rebuild)
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
EXPORT_CODE void CONVENTION get_1phase_LUT_params(int *nT, int *np, double *Tmin, double *Tmax, double *pmin, double *pmax){
	pFluid->get_1phase_LUT_params(nT,np,Tmin,Tmax,pmin,pmax);
}

EXPORT_CODE void CONVENTION UseSaturationLUT(bool OnOff)
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
EXPORT_CODE bool CONVENTION SaturationLUTStatus()
{
	return FlagUseSaturationLUT;
}

EXPORT_CODE void CONVENTION UseSinglePhaseLUT(bool OnOff)
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


EXPORT_CODE bool CONVENTION SinglePhaseLUTStatus(void)
{
    return FlagUseSinglePhaseLUT;
}

long get_Fluid_index(std::string FluidName)
{
	// Try to get the fluid from Fluids by name
	pFluid = Fluids.get_fluid(FluidName);
	// If NULL, didn't find it (or its alias)
	if (pFluid!=NULL)
	{
		// Find the fluid index
		return Fluids.get_fluid_index(pFluid);
	}
	else
		return -1;
}
EXPORT_CODE long CONVENTION get_Fluid_index(char * param)
{
	return get_Fluid_index(std::string(param));
}

std::string get_index_units(long index)
{
	std::map<long, std::string>::iterator it;
	// Try to find using the map
	it = units_map.find(index);
	// If it is found the iterator will not be equal to end
	if (it != units_map.end() )
	{
		// Return the index of the parameter
		return (*it).second;
	}
	else
	{
		std::cout << "Didn't match parameter: " << index << std::endl;
		return std::string("Didn't match parameter");
	}
}
EXPORT_CODE void CONVENTION get_index_units(long param, char * units)
{
	strcpy(units, (char*)get_index_units(param).c_str());
	return;
}

long get_param_index(std::string param)
{
	std::map<std::string,long>::iterator it;
	// Try to find using the map
	it = param_map.find(param);
	// If it is found the iterator will not be equal to end
	if (it != param_map.end() )
	{
		// Return the index of the parameter
		return (*it).second;
	}
	else
	{
		std::cout << "Didn't match parameter: " << param << std::endl;
		return -1;
	}
}
EXPORT_CODE long CONVENTION get_param_index(char * param)
{
	return get_param_index(std::string(param));
}
static int IsCoolPropFluid(std::string FluidName)
{
	// Try to get the fluid from Fluids by name
	try
	{
		pFluid = Fluids.get_fluid(FluidName);
	}
	catch (NotImplementedError)
	{
		return NULL;
	}
	// If NULL, didn't find it (or its alias)
	if (pFluid!=NULL)
	{
		return true;
	}
	else
		return false;
}

static int IsCoolPropFluid(char* Fluid)
{
	return IsCoolPropFluid(std::string(Fluid));
}

static int IsBrine(char* Ref)
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
EXPORT_CODE int CONVENTION IsFluidType(char *Ref, char *Type)
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
EXPORT_CODE void CONVENTION Phase(char *Fluid,double T, double p, char *Phase_str)
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
EXPORT_CODE double CONVENTION Props1(char* Ref, char * Output)
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
				else if (!strcmp(Output,"rhocrit"))
					return Props('N','T',0,'P',0,Fluid);
				else if (!strcmp(Output,"accentric"))
					return Props('w','T',0,'P',0,Fluid);
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
		else if (!strcmp(Output,"accentric"))
			return pFluid->params.accentricfactor;
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

EXPORT_CODE double CONVENTION Props(char *Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
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
		return _Props(Output,std::string(1,Name1),Prop1,std::string(1,Name2),Prop2,Ref);
	}
	catch(std::exception &e){
		err_string=std::string("CoolProp error: ").append(e.what());
		return _HUGE;
	}
	catch(...){
		err_string=std::string("CoolProp error: Indeterminate error");
		return _HUGE;
	}
}


// Make this a wrapped function so that error bubbling can be done properly
double _Props(std::string Output,std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref)
{
    
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

	if (IsCoolPropFluid(Ref))
	{
		pFluid = Fluids.get_fluid(Ref);
		// Convert all the parameters to integers
		long iOutput = get_param_index(Output);
		long iName1 = get_param_index(std::string(Name1));
		long iName2 = get_param_index(std::string(Name2));
		// Call the internal method that uses the parameters converted to longs
		return _CoolProp_Fluid_Props(iOutput,iName1,Prop1,iName2,Prop2,pFluid);
	}

    /* 
    If the fluid name is not actually a refrigerant name, but a string beginning with "REFPROP-",
    then REFPROP is used to calculate the desired property.
    */
    else if (IsREFPROP(Ref))  // First eight characters match "REFPROP-"
    {
        #if defined(__ISWINDOWS__)
		return REFPROP(Output.c_str()[0],Name1.c_str()[0],Prop1,Name2.c_str()[0],Prop2,(char*)Ref.c_str());
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
	else if (IsBrine((char*)Ref.c_str()))
    {
		//Enthalpy and pressure are the inputs
		if ((Name1.c_str()[0]=='H' && Name2.c_str()[0]=='P') || (Name2.c_str()[0]=='H' && Name1.c_str()[0]=='P'))
        {
			if (Name2.c_str()[0]=='H' && Name1.c_str()[0]=='P')
			{
				swap(&Prop1,&Prop2);
				swap(&Name1,&Name2);
			}
			// Start with a guess of 10 K below max temp of fluid
			double Tguess = SecFluids('M',Prop1,Prop2,(char*)Ref.c_str())-10;
			// Get the temperature
			double T =_T_hp_secant(Ref,Prop1,Prop2,Tguess);
			// Return whatever property is desired
			return SecFluids(Output[0],T,Prop2,(char*)Ref.c_str());
		}
		else if (Name1.c_str()[0]!='T' || Name2.c_str()[0]!='P')
        {
			throw ValueError("For brine, Name1 must be 'T' and Name2 must be 'P' or Name1 must be 'H' and Name2 must be 'P'");
		}
		else
		{
			return SecFluids(Output[0],Prop1,Prop2,(char*)Ref.c_str());
		}
    }
    return 0;
}
double _CoolProp_Fluid_Props(long iOutput, long iName1, double Prop1, long iName2, double Prop2, Fluid *pFluid)
{
	double T,Q,rhoV,rhoL,Value,rho,pL,pV;

	// This private method uses the indices directly for speed

	// Check if it is an output that doesn't require a state input
    // Deal with it and return

	switch (iOutput)
	{
		case iMM:
			return pFluid->params.molemass;
			break;
		case iPcrit:
			return pFluid->crit.p;
			break;
		case iTcrit:
			return pFluid->crit.T;
			break;
		case iTtriple:
			return pFluid->params.Ttriple;
			break;
		case iRhocrit:
			return pFluid->reduce.rho;
			break;
		case iAccentric: 
			return pFluid->params.accentricfactor;
			break;
	}

	if (iName1 == iT && Prop1 < pFluid->limits.Tmin) 
		throw ValueError(format("Input temperature to Props function [%f K] is below the fluid minimum temp [%f K]",Prop1,pFluid->limits.Tmin));
	if (iName2 == iT && Prop2 < pFluid->limits.Tmin) 
		throw ValueError(format("Input temperature to Props function [%f K] is below the fluid minimum temp [%f K]",Prop2,pFluid->limits.Tmin));

	//Surface tension is only a function of temperature
	if (iOutput == iI){
		if (iName1 == iT)
			return pFluid->surface_tension_T(Prop1)/1000;
		else if (iName2 == iT)
			return pFluid->surface_tension_T(Prop2)/1000;
		else
			throw ValueError(format("If output is surface tension ['I' or 'SurfaceTension'], Param1 must be temperature"));
	}

	// In any case, you want to get a (temperature, density) pair
	if ((iName1 == iT && iName2 == iP) || (iName1 == iP && iName2 == iT))
	{
		//Swap values and keys
		if (iName1 == iP && iName2 == iT)
		{
			swap(&Prop1,&Prop2);
			swap(&iName1,&iName2);
		}

		// Handle trivial outputs
		if (iOutput == iP) 
			return Prop2;
		else if (iOutput == iT) 
			return Prop1;
		
		// Get density as a function of T&p, then call again
		if (FlagUseSinglePhaseLUT==true){
			return pFluid->LookupValue_TP(std::string((char*)iOutput), Prop1, Prop2);
        }
		else{
			rho = rho_TP(Prop1,Prop2);
		}
		
		if (iOutput == iD){
			if (debug()>5){
				std::cout<<__FILE__<<": "<<iOutput<<","<<iName1<<","<<Prop1<<","<<iName2<<","<<Prop2<<","<<pFluid->get_name()<<"="<<rho<<std::endl;
			}
			return rho;
		}
		else
			return _CoolProp_Fluid_Props(iOutput,iName1,Prop1,iD,rho,pFluid);
	}
	else if ((iName1 == iT && iName2 == iQ) || (iName1 == iQ && iName2 == iT))
	{
		if (iName1 == iQ && iName2 == iT){
			//Swap values and keys to get order of T, Q
			swap(&Prop1,&Prop2);
			swap(&iName1,&iName2);
		}
		
		T = Prop1;
		Q = Prop2;
		if (T <= pFluid->params.Ttriple || T >= pFluid->reduce.T){
			throw ValueError(format("Your saturation temperature [%f K] is out of range [%f K, %f K]",T,pFluid->params.Ttriple,pFluid->reduce.T ));
		}
		if (Q>1+1e-13 || Q<-1e-13){
			throw ValueError(format("Your quality [%f] is out of range (0, 1)",Q ));
		}
		// Get the saturation properties
		pFluid->saturation(Prop1,FlagUseSaturationLUT,&pL,&pV,&rhoL,&rhoV);
		// Find the effective density to use
		rho=1/(Q/rhoV+(1-Q)/rhoL);

		// Trivial output
		if (iOutput == iP)
			return Q*pV+(1-Q)*pL;

		// Recurse and call Props again with the calculated density
		if (iOutput == iD)
			return 1/(Q/rhoV+(1-Q)/rhoL);
		else if (iOutput == iC || iOutput == iO)
			return _CoolProp_Fluid_Props(iOutput,iT,Prop1,iD,rho,pFluid);
		else
			if (fabs(Q)<1e-12)
				return _CoolProp_Fluid_Props(iOutput,iT,Prop1,iD,rhoL,pFluid);
			else if (fabs(Q-1)<1e-12)
				return _CoolProp_Fluid_Props(iOutput,iT,Prop1,iD,rhoV,pFluid);
			else
				return Q*_CoolProp_Fluid_Props(iOutput,iT,Prop1,iD,rhoV,pFluid)+(1-Q)*_CoolProp_Fluid_Props(iOutput,iT,Prop1,iD,rhoL,pFluid);
	}
	else if ((iName1 == iT && iName2 == iD) || (iName1 == iD && iName2 == iT))
	{
		if (iName1 == iD && iName2 == iT)
		{
			//Swap values and keys to get T,D
			swap(&Prop1,&Prop2);
			swap(&iName1,&iName2);
		}
		T=Prop1;
		rho=Prop2;
		// Trivial output
		if (iOutput == iD)
			return rho;
		else if (iOutput == iT)
			return T;

		// If you are using LUT, use it
		if (FlagUseSinglePhaseLUT==1){
			// Try to use the LUT, if the parameter is not included in the LUT,
			// allow it to fall back to the conventional analysis
			try{
				Value = pFluid->LookupValue_Trho(std::string((char*)iOutput), T, rho);
				return Value;
			}
            catch(ValueError){

            }
        }
		rho = Prop2; 
		switch (iOutput)
		{
			case iP:
				Value=pFluid->pressure_Trho(T,rho);
				break;
			case iH:
				Value=pFluid->enthalpy_Trho(T,rho);
				break;
			case iS:
				Value=pFluid->entropy_Trho(T,rho);
				break;
			case iU:
				Value=pFluid->internal_energy_Trho(T,rho);
				break;
			case iC:
				Value=pFluid->specific_heat_p_Trho(T,rho);
				break;
			case iC0:
				Value=pFluid->specific_heat_p_ideal_Trho(T);
				break;
			case iO:
				Value=pFluid->specific_heat_v_Trho(T,rho);
				break;
			case iA:
				Value=pFluid->speed_sound_Trho(T,rho);
				break;
			case iG:
				Value=pFluid->gibbs_Trho(T,rho);
				break;
			case iV:
				Value=pFluid->viscosity_Trho(T,rho);
				break;
			case iL:
				Value=pFluid->conductivity_Trho(T,rho);
				break;
			case iDpdT:
				Value=pFluid->dpdT_Trho(T,rho);
				break;
			default:
				throw ValueError(format("Invalid Output index: %d ",iOutput));
				return _HUGE;
        }
		if (debug()>5){
			std::cout<<__FILE__<<"More output" <<std::endl;
			std::cout<<__FILE__<<__LINE__<<": "<<iOutput<<","<<iName1<<","<<Prop1<<","<<iName2<<","<<Prop2<<","<<pFluid->get_name()<<"="<<Value<<std::endl;
			std::cout<<__FILE__<<__LINE__<<": "<<iOutput<<","<<iName1<<","<<Prop1<<","<<iName2<<","<<Prop2<<","<<pFluid->get_name()<<"="<<Value<<std::endl;
		}
        return Value;
	}
    else if (iName1 == iP && iName2 == iQ)
    {
        T=pFluid->Tsat(Prop1,Prop2,0);
        return _CoolProp_Fluid_Props(iOutput,iT,T,iQ,Prop2,pFluid);
    }
	else if (iName1 == iQ && iName2 == iP)
    {
        T=pFluid->Tsat(Prop2,Prop1,0);
        return _CoolProp_Fluid_Props(iOutput,iT,T,iQ,Prop1,pFluid);
    }
    else if (iName1 == iH && iName2 == iP)
    {
    	_T_hp(pFluid->get_name(),Prop1,Prop2,&T, &rho);
		return _CoolProp_Fluid_Props(iOutput,iT,T,iD,rho,pFluid);
    }
	else if (iName1 == iP && iName2 == iH)
    {
    	_T_hp(pFluid->get_name(),Prop2,Prop1,&T, &rho);
		return _CoolProp_Fluid_Props(iOutput,iT,T,iD,rho,pFluid);
    }
	else if (iName1 == iS && iName2 == iP)
    {
    	_T_sp(pFluid->get_name(),Prop1,Prop2,&T, &rho);
		return _CoolProp_Fluid_Props(iOutput,iT,T,iD,rho,pFluid);
    }
	else if (iName1 == iP && iName2 == iS)
    {
		_T_sp(pFluid->get_name(),Prop2,Prop1,&T, &rho);
		return _CoolProp_Fluid_Props(iOutput,iT,T,iD,rho,pFluid);
    }
    else
    {
		throw ValueError(format("Not a valid pair of input keys %d,%d and output key %d",iName1,iName2,iOutput));
    }
}
EXPORT_CODE double CONVENTION IProps(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid)
{
	pFluid = Fluids.get_fluid(iFluid);
	// Didn't work
	if (pFluid == NULL)
		return _HUGE;
	else
		return _CoolProp_Fluid_Props(iOutput,iName1,Prop1,iName2,Prop2,pFluid);
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
	double f1,f2,df1_dtau,df1_ddelta,df2_ddelta,df2_dtau,s_hot,cp,s1,T1,p1;
    double rhosatL,rhosatV,ssatL,ssatV,TsatL,TsatV,tau,delta,worst_error;
	
	//First figure out where you are

	double Tc = pFluid->reduce.T;
	double rhoc = pFluid->reduce.rho;
	double pc = pFluid->reduce.p;

	R=pFluid->R();
	if (p > Props(Ref,"pcrit"))
	{
        // Supercritical; use critical point as anchor to determine guess for T
		// by assuming it is ideal gas
		// For isentropic, T2/T1 = (p2/p1)^(R/cp); T1,p1, are critical state
		T1 = Tc + 5.0;
		p1 = pc + 5.0;
		s1 = Props(std::string("S"),'T',T1,'P',p1,Ref);
		cp = Props(std::string("C"),'T',T1,'P',p1,Ref);
		T_guess = T1*exp((s-s1+R*log(p/p1))/cp);
		delta = p/(R*T_guess)/pFluid->reduce.rho;
	}
	else
	{
		rhosatL=Props(std::string("D"),'P',p,'Q',0.0,Ref);
		rhosatV=Props(std::string("D"),'P',p,'Q',1.0,Ref);
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
			// It is two-phase
			// Return the quality weighted values
			double quality = (s-ssatL)/(ssatV-ssatL);
			*Tout = quality*TsatV+(1-quality)*TsatL;
			double v = quality*(1/rhosatV)+(1-quality)*1/rhosatL;
			*rhoout = 1/v;

		}
	}

	tau=Tc/T_guess;

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
    double rhosatL,rhosatV,hsatL,hsatV,TsatL,TsatV,tau,delta,worst_error;
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
		rhosatV=Props(std::string("D"),'P',p,'Q',1.0,Ref);
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
			// It is two-phase
			// Return the quality weighted values
			double quality = (h-hsatL)/(hsatV-hsatL);
			*Tout = quality*TsatV+(1-quality)*TsatL;
			double v = quality*(1/rhosatV)+(1-quality)*1/rhosatL;
			*rhoout = 1/v;
			return;
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


//
//double h_sp(char *Ref, double s, double p, double T_guess)
//{
//    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f=999,T=300;
//    int iter=1;
//    
//    while ((iter<=3 || change>eps) && iter<100)
//    {
//        if (iter==1){x1=T_guess; T=x1;}
//        if (iter==2){x2=T_guess+1.0; T=x2;}
//        if (iter>2) {T=x2;}
//
//            // Find the temperature which gives the same entropy
//            f=Props('S','T',T,'P',p,Ref)-s;
//
//        if (iter==1){y1=f;}
//        if (iter>1)
//        {
//            y2=f;
//            x3=x2-y2/(y2-y1)*(x2-x1);
//            change=fabs(y2/(y2-y1)*(x2-x1));
//            y1=y2; x1=x2; x2=x3;
//        }
//        iter=iter+1;
//        if (iter>50)
//        {
//        	//ERROR
//            printf("h_sp not converging with inputs(%s,%g,%g,%g)\n",Ref,s,p,T_guess);
//        }
//    }
//    return Props('H','T',T,'P',p,Ref);
//}

EXPORT_CODE double CONVENTION K2F(double T)
{ return T * 9 / 5 - 459.67; }

EXPORT_CODE double CONVENTION F2K(double T_F)
{ return (T_F + 459.67) * 5 / 9;}

EXPORT_CODE void CONVENTION PrintSaturationTable(char *FileName, char * Ref,double Tmin, double Tmax)
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

EXPORT_CODE void CONVENTION FluidsList(char* str)
{
	str=(char*)FluidsList().c_str();
	return;
}
std::string FluidsList()
{
	return Fluids.FluidList();
}

EXPORT_CODE double CONVENTION DerivTerms(char *Term,double T, double rho, char * Ref)
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
EXPORT_CODE void CONVENTION get_REFPROPname(char* Ref, char * str)
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
