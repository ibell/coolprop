#include "coolpropsolver.h"
#include "CoolPropTools.h"
#include "CoolProp.h"
#include "CPState.h"
#include <iostream>
#include <string>
#include <stdlib.h>

long iFluid;
CoolPropSolver::CoolPropSolver(const string &mediumName, const string &libraryName, const string &substanceName)
	: BaseSolver(mediumName, libraryName, substanceName){
	setFluidConstants();

	iFluid = get_Fluid_index(substanceName);

	// Fluid name can be used to pass in other parameters.  
	// The string can be composed like "Propane|enable_TTSE=1|calc_transport=0"

	std::vector<string> name_options = strsplit(substanceName,'|');
	
	// Set the defaults
	enable_TTSE = false;
	debug_level = 0;
	calc_transport = false;
	twophase_derivsmoothing_xend = 0;

	if (name_options.size()>1)
	{
		for (unsigned int i = 1; i<name_options.size(); i++)
		{
			// Split around the equals sign
			std::vector<string> param_val = strsplit(name_options[i],'=');
			if (param_val.size() != 2)
			{
				errorMessage((char*)format("Could not parse the option [%s], must be in the form param=value",name_options[i].c_str()).c_str());
			}

			// Check each of the options in turn
			if (!param_val[0].compare("enable_TTSE"))
			{
				if (!param_val[1].compare("1") || !param_val[1].compare("true"))
				{
					std::cout << "TTSE is on";
					enable_TTSE = true;
				}
				else if (!param_val[1].compare("0") || !param_val[1].compare("false"))
				{
					std::cout << "TTSE is off";
					enable_TTSE = false;
				}
				else
					errorMessage((char*)format("I don't know how to handle this option [%s]",name_options[i].c_str()).c_str());
			}
			else if (!param_val[0].compare("calc_transport"))
			{
				if (!param_val[1].compare("1") || !param_val[1].compare("true"))
					calc_transport = true;
				else if (!param_val[1].compare("0") || !param_val[1].compare("false"))
					calc_transport = false;
				else
					errorMessage((char*)format("I don't know how to handle this option [%s]",name_options[i].c_str()).c_str());
			}
			else if (!param_val[0].compare("twophase_derivsmoothing_xend"))
			{
				twophase_derivsmoothing_xend = strtod(param_val[1].c_str(),NULL);
				if (twophase_derivsmoothing_xend<0 || twophase_derivsmoothing_xend > 1)
					errorMessage((char*)format("I don't know how to handle this twophase_derivsmoothing_xend value [%d]",param_val[0].c_str()).c_str());
			}
			else if (!param_val[0].compare("debug"))
			{
				debug_level = (int)strtol(param_val[1].c_str(),NULL,NULL);
				if (debug_level<0 || debug_level > 1000)
					errorMessage((char*)format("I don't know how to handle this debug level [%s]",param_val[0].c_str()).c_str());
			}
			else
			{
				errorMessage((char*)format("This option [%s] was not understood",name_options[i].c_str()).c_str());
			}
			
			// Some options were passed in, lets see what we have
			std::cout << param_val[0] << " has the value of " << param_val[1] << std::endl;
		}
	}

	state = new CoolPropStateClass(name_options[0]);
}

void CoolPropSolver::setFluidConstants(){
	_fluidConstants.pc = Props(substanceName,"pcrit");
	_fluidConstants.Tc = Props(substanceName,"Tcrit");
	_fluidConstants.MM = Props(substanceName,"molemass");
	_fluidConstants.dc = Props(substanceName,"rhocrit");
}

void CoolPropSolver::setSat_p(double &p, ExternalSaturationProperties *const properties){

	// Convert to CoolProp input units
	p /= 1000.0;

	if (debug_level > 5)
		std::cout << format("setSat_p(%0.16e)\n",p);

	if (enable_TTSE)
		state->enable_TTSE_LUT();
	else
		state->disable_TTSE_LUT();
	try
	{
		state->update(iP,p,iQ,0); // quality only matters for pseudo-pure fluids
		properties->psat = p*1000; //[Pa --> kPa]
		properties->Tsat = state->TL(); // Not correct for pseudo-pure fluids
		properties->dl = state->rhoL();
		properties->dv = state->rhoV();
		properties->hl = state->hL()*1000;
		properties->hv = state->hV()*1000;
		properties->dTp = state->dTdp_along_sat()/1000; //[1/kPa --> 1/Pa]
		
		properties->ddldp = state->drhodp_along_sat_liquid()/1000; //[1/kPa -- > 1/Pa]
		properties->ddvdp = state->drhodp_along_sat_vapor()/1000; //[1/kPa -- > 1/Pa]
		properties->dhldp = state->dhdp_along_sat_liquid(); // [kJ/kg/kPa --> J/kg/Pa]
		properties->dhvdp = state->dhdp_along_sat_vapor(); // [kJ/kg/kPa --> J/kg/Pa]
	}
	catch(std::exception &e)
	{
		errorMessage((char*)e.what());
	}
}

void CoolPropSolver::setSat_T(double &T, ExternalSaturationProperties *const properties){

	if (debug_level > 5)
		std::cout << format("setSat_T(%0.16e)\n",T);

	if (enable_TTSE)
		state->enable_TTSE_LUT();
	else
		state->disable_TTSE_LUT();
	try
	{
		state->update(iT,T,iQ,0); // Quality only matters for pseudo-pure fluids

		properties->Tsat = T;
		properties->psat = state->p()*1000;
		properties->dl = state->rhoL();
		properties->dv = state->rhoV();
		properties->hl = state->hL()*1000;
		properties->hv = state->hV()*1000;
		properties->dTp = state->dTdp_along_sat()/1000; //[1/kPa --> 1/Pa]
		
		properties->ddldp = state->drhodp_along_sat_liquid()/1000; //[1/kPa -- > 1/Pa]
		properties->ddvdp = state->drhodp_along_sat_vapor()/1000; //[1/kPa -- > 1/Pa]
		properties->dhldp = state->dhdp_along_sat_liquid(); // [kJ/kg/kPa --> J/kg/Pa]
		properties->dhvdp = state->dhdp_along_sat_vapor(); // [kJ/kg/kPa --> J/kg/Pa]
	}
	catch(std::exception &e)
	{
		errorMessage((char*)e.what());
	}
}

// Note: the phase input is currently not supported
void CoolPropSolver::setState_ph(double &p, double &h, int &phase, ExternalThermodynamicState *const properties){

	//Convert to CoolProp input units
	h /= 1000; // [J/kg --> kJ/kg]
	p /= 1000; // [Pa --> kPa]

	if (debug_level > 5)
		std::cout << format("setState_ph(p=%0.16e,h=%0.16e)\n",p,h);

	if (enable_TTSE)
		state->enable_TTSE_LUT();
	else
		state->disable_TTSE_LUT();
	try{
		// Update the internal variables in the state instance
		state->update(iP,p,iH,h);
		
		if (!ValidNumber(state->rho()) || !ValidNumber(state->T()))
		{
			throw ValueError(format("p-h [%g, %g] failed for update",p,h));
		}
		// Set the values in the output structure
		properties->p = p*1000;
		properties->h = h*1000;

		properties->d = state->rho();
		properties->T = state->T();
		properties->s = state->s()*1000;
		
		if (state->TwoPhase){
			properties->phase = 2; 
		}
		else{
			properties->phase = 1; 
		}
		properties->cp = state->cp()*1000;
		properties->cv = state->cv()*1000;
		properties->a = state->speed_sound();
		if (state->TwoPhase && state->Q() >= 0 && state->Q() <= twophase_derivsmoothing_xend)
		{
			// Use the smoothed derivatives between a quality of 0 and twophase_derivsmoothing_xend
            properties->ddhp = state->drhodh_constp_smoothed(twophase_derivsmoothing_xend)/1000; // [1/kPa -- > 1/Pa]
            properties->ddph = state->drhodp_consth_smoothed(twophase_derivsmoothing_xend)/1000; // [1/(kJ/kg) -- > 1/(J/kg)]
		}
		else
		{
            properties->ddhp = state->drhodh_constp()/1000; // [1/kPa -- > 1/Pa]
			properties->ddph = state->drhodp_consth()/1000; // [1/(kJ/kg) -- > 1/(J/kg)]
		}
		properties->kappa = state->isothermal_compressibility()/1000; // [1/kPa -- > 1/Pa]
		properties->beta = state->isobaric_expansion_coefficient();

		if (calc_transport)
        {
            properties->eta = state->viscosity();
            properties->lambda = state->conductivity()*1000; //[kW/m/K --> W/m/K]
            properties->Pr = properties->cp*properties->eta/properties->lambda;
        }
	}
	catch(std::exception &e)
	{
		errorMessage((char*)e.what());
	}
}

void CoolPropSolver::setState_pT(double &p, double &T, ExternalThermodynamicState *const properties){
	//Convert to CoolProp inputs
	p /= 1000; // [Pa --> kPa]

	if (debug_level > 5)
		std::cout << format("setState_pT(p=%0.16e,T=%0.16e)\n",p,T);

	if (enable_TTSE)
		state->enable_TTSE_LUT();
	else
		state->disable_TTSE_LUT();
	try{
		// Update the internal variables in the state instance
		state->update(iP,p,iT,T);

		// Set the values in the output structure
		properties->p = p*1000;
		properties->T = T;
		properties->d = state->rho();
		properties->h = state->h()*1000;
		properties->s = state->s()*1000;
		properties->phase = 1; // with pT input, always one-phase conditions!
		properties->cp = state->cp()*1000;
		properties->cv = state->cv()*1000;
		properties->a = state->speed_sound();
		if (state->TwoPhase && state->Q() >= 0 && state->Q() <= twophase_derivsmoothing_xend)
		{
			// Use the smoothed derivatives between a quality of 0 and twophase_derivsmoothing_xend
            properties->ddhp = state->drhodh_constp_smoothed(twophase_derivsmoothing_xend)/1000; // [1/kPa -- > 1/Pa]
            properties->ddph = state->drhodp_consth_smoothed(twophase_derivsmoothing_xend)/1000; // [1/(kJ/kg) -- > 1/(J/kg)]
		}
		else
		{
            properties->ddhp = state->drhodh_constp()/1000; // [1/kPa -- > 1/Pa]
			properties->ddph = state->drhodp_consth()/1000; // [1/(kJ/kg) -- > 1/(J/kg)]
		}
		properties->kappa = properties->d/properties->p/state->dpdrho_constT()/1000; // [1/kPa -- > 1/Pa]  
		properties->beta = -1/properties->d*state->drhodT_constp();
		if (calc_transport)
        {
            properties->eta = state->viscosity();
            properties->lambda = state->conductivity()*1000; //[kW/m/K --> W/m/K]
            properties->Pr = properties->cp*properties->eta/properties->lambda;
        }
	}
	catch(std::exception &e)
	{
		errorMessage((char*)e.what());
	}
}

// Note: the phase input is currently not supported
void CoolPropSolver::setState_dT(double &d, double &T, int &phase, ExternalThermodynamicState *const properties)
{
	// Variables are already in the correct units

	if (debug_level > 5)
		std::cout << format("setState_dT(d=%0.16e,T=%0.16e)\n",d,T);

	if (enable_TTSE)
		state->enable_TTSE_LUT();
	else
		state->disable_TTSE_LUT();
	try{

		// Update the internal variables in the state instance
		state->update(iD,d,iT,T);

		// Set the values in the output structure
		properties->d = d;
		properties->T = T;
		properties->h = state->h()*1000;
		properties->p = state->p()*1000;
		properties->s = state->s()*1000;
		if (state->TwoPhase){
			properties->phase = 2; 
		}
		else{
			properties->phase = 1; 
		}
		properties->cp = state->cp()*1000;
		properties->cv = state->cv()*1000;
		properties->a = state->speed_sound();
		if (state->TwoPhase && state->Q() >= 0 && state->Q() <= twophase_derivsmoothing_xend)
		{
			// Use the smoothed derivatives between a quality of 0 and twophase_derivsmoothing_xend
            properties->ddhp = state->drhodh_constp_smoothed(twophase_derivsmoothing_xend)/1000; // [1/kPa -- > 1/Pa]
            properties->ddph = state->drhodp_consth_smoothed(twophase_derivsmoothing_xend)/1000; // [1/(kJ/kg) -- > 1/(J/kg)]
		}
		else
		{
            properties->ddhp = state->drhodh_constp()/1000; // [1/kPa -- > 1/Pa]
			properties->ddph = state->drhodp_consth()/1000; // [1/(kJ/kg) -- > 1/(J/kg)]
		}
		properties->kappa = properties->d/properties->p/state->dpdrho_constT()/1000; // [1/kPa -- > 1/Pa]  
		properties->beta = -1/properties->d*state->drhodT_constp();
		if (calc_transport)
        {
            properties->eta = state->viscosity();
            properties->lambda = state->conductivity()*1000; //[kW/m/K --> W/m/K]
            properties->Pr = properties->cp*properties->eta/properties->lambda;
        }
	}
	catch(std::exception &e)
	{
		errorMessage((char*)e.what());
	}	
}

// Note: the phase input is currently not supported
void CoolPropSolver::setState_ps(double &p, double &s, int &phase, ExternalThermodynamicState *const properties){
	//Convert to CoolProp inputs
	s /= 1000; // [J/kg/K --> kJ/kg/K]
	p /= 1000; // [Pa --> kPa] 

	if (debug_level > 5)
		std::cout << format("setState_ps(p=%0.16e,s=%0.16e)\n",p,s);

	if (enable_TTSE)
		state->enable_TTSE_LUT();
	else
		state->disable_TTSE_LUT();
	try{
		// Update the internal variables in the state instance
		state->update(iP,p,iS,s);

		// Set the values in the output structure
		properties->p = p*1000; // [kPa --> Pa]
		properties->s = s*1000; // [kJ/kg/K --> J/kg/K]
		properties->d = state->rho();
		properties->T = state->T();
		properties->h = state->h()*1000; //[kJ/kg/K --> J/kg/K]
		properties->phase = 1; // with pT input, always one-phase conditions!
		properties->cp = state->cp()*1000;
		properties->cv = state->cv()*1000;
		properties->a = state->speed_sound();
        if (state->TwoPhase && state->Q() >= 0 && state->Q() <= twophase_derivsmoothing_xend)
		{
			// Use the smoothed derivatives between a quality of 0 and twophase_derivsmoothing_xend
            properties->ddhp = state->drhodh_constp_smoothed(twophase_derivsmoothing_xend)/1000; // [1/kPa -- > 1/Pa]
            properties->ddph = state->drhodp_consth_smoothed(twophase_derivsmoothing_xend)/1000; // [1/(kJ/kg) -- > 1/(J/kg)]
		}
		else
		{
            properties->ddhp = state->drhodh_constp()/1000; // [1/kPa -- > 1/Pa]
			properties->ddph = state->drhodp_consth()/1000; // [1/(kJ/kg) -- > 1/(J/kg)]
		}
		properties->kappa = properties->d/properties->p/state->dpdrho_constT()/1000; // [1/kPa -- > 1/Pa]  
		properties->beta = -1/properties->d*state->drhodT_constp();
        if (calc_transport)
        {
            properties->eta = state->viscosity();
            properties->lambda = state->conductivity()*1000; //[kW/m/K --> W/m/K]
            properties->Pr = properties->cp*properties->eta/properties->lambda;
        }
	}
	catch(std::exception &e)
	{
		errorMessage((char*)e.what());
	}
}

