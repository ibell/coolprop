#include <iostream>
#include "FluidClass.h"
#include "CoolProp.h"

#ifndef CPSTATE_H
#define CPSTATE_H

bool match_pair(long iI1, long iI2, long I1, long I2);
void sort_pair(long *iInput1, double *Value1, long *iInput2, double *Value2, long I1, long I2);

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
	double _rho,_T,_p,_Q,_h,_s;

	// To be used to update internal variables if you know that your parameters are T,Q or P,Q
	void update_twophase(long iInput1, double Value1, long iInput2, double Value2);

	// To be used to update internal variables if you know that your parameters are T,D
	void update_Trho(long iInput1, double Value1, long iInput2, double Value2);

	// To be used to update internal variables if you know that your parameters are T,P
	void update_Tp(long iInput1, double Value1, long iInput2, double Value2);

	/// To be used to update internal variables if you know that your parameters are P,H
	/// If T0 and rho0 are included, start at this value for the solution
	void update_ph(long iInput1, double Value1, long iInput2, double Value2, double T0 = -1, double rho0 = -1);

	// To be used to update internal variables if you know that your parameters are P,S
	void update_ps(long iInput1, double Value1, long iInput2, double Value2);

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
	double hL(void);
	double hV(void);
	double sL(void);
	double sV(void);
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