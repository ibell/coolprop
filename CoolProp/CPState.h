#include <iostream>
#include "FluidClass.h"
#include "CoolProp.h"

#ifndef CPSTATE_H
#define CPSTATE_H

class CoolPropStateClass
{
protected:
	std::string _Fluid;
	Fluid * pFluid;
	bool flag_SinglePhase;
	bool SaturatedL,SaturatedV;

	// Saturation values
	double rhosatL, rhosatV, psatL, psatV, TsatL, TsatV;

	// Bulk values
	double _rho,_T,_p,_Q;

	// To be used to update internal variables if you know that your parameters are T,Q or P,Q
	void update_twophase(long iInput1, double Value1, long iInput2, double Value2);

	// To be used to update internal variables if you know that your parameters are T,D
	void update_Trho(long iInput1, double Value1, long iInput2, double Value2);

	// To be used to update internal variables if you know that your parameters are T,P
	void update_Tp(long iInput1, double Value1, long iInput2, double Value2);

	// To be used to update internal variables if you know that your parameters are P,H
	void update_ph(long iInput1, double Value1, long iInput2, double Value2);

	// Check whether the quality corresponds to saturated liquid or vapor
	void check_saturated_quality(double Q);
public:

	// Phase flags
	bool TwoPhase,SinglePhase;
	
	// Constructor with fluid name
	CoolPropStateClass(std::string FluidName);

	// Constructor with fluid pointer
	CoolPropStateClass(Fluid *pFluid);

	// Property updater
	// Uses the indices in CoolProp for the input parameters
	void update(long iInput1, double Value1, long iInput2, double Value2);

	// Property accessors for saturation parameters directly
	// These all must be calculated every time if the state is saturated or two-phase
	double rhoL(void){return rhosatL;};
	double rhoV(void){return rhosatV;};
	double pL(void){return psatL;};
	double pV(void){return psatV;};
	double TL(void){return TsatL;};
	double TV(void){return TsatV;};
	// Derived parameters for the saturation states
	double hL(void){return pFluid->enthalpy_Trho(TsatL,rhosatL);};
	double hV(void){return pFluid->enthalpy_Trho(TsatV,rhosatV);};
	double sL(void){return pFluid->entropy_Trho(TsatL,rhosatL);};
	double sV(void){return pFluid->entropy_Trho(TsatV,rhosatV);};
	double cpL(void){return pFluid->specific_heat_p_Trho(TsatL,rhosatL);};
	double cpV(void){return pFluid->specific_heat_p_Trho(TsatV,rhosatV);};
	// Derivatives along the saturation curve
	void dvdp_dhdp_sat(double T, double *dvdpL, double *dvdpV, double *dhdpL, double *dhdpV);

	// Bulk properties accessors - temperature and density are directly calculated every time
	// All other parameters are calculated on an as-needed basis
	// If single-phase, just plug into the EOS, otherwise need to do two-phase analysis
	double T(void){return _T;};
	double rho(void){return _rho;};
	double p(void){return _p;};
	double h(void);
	double s(void);
	double cp(void);
	double cv(void);
	double speed_sound(void);

	double drhodT_constp(void);
	double drhodh_constp(void);
	double drhodp_consth(void);
	double dpdrho_constT(void);

};

#endif