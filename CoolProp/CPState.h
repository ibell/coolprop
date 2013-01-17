#include <iostream>
#include "FluidClass.h"
#include "CoolProp.h"
#include "TTSE.h"

#ifndef CPSTATE_H
#define CPSTATE_H

bool match_pair(long iI1, long iI2, long I1, long I2);
void sort_pair(long *iInput1, double *Value1, long *iInput2, double *Value2, long I1, long I2);

class CoolPropStateClass
{
protected:
	// Helmholtz derivative cache flags
	bool cached_phi0, cached_dphi0_dTau, cached_dphi0_dDelta,  cached_d2phi0_dTau2, cached_d2phi0_dDelta_dTau, cached_d2phi0_dDelta2, cached_d3phi0_dTau3, cached_d3phi0_dDelta_dTau2, cached_d3phi0_dDelta2_dTau, cached_d3phi0_dDelta3; 
	double cachedval_phi0, cachedval_dphi0_dTau, cachedval_dphi0_dDelta,  cachedval_d2phi0_dTau2, cachedval_d2phi0_dDelta_dTau, cachedval_d2phi0_dDelta2, cachedval_d3phi0_dTau3,  cachedval_d3phi0_dDelta_dTau2, cachedval_d3phi0_dDelta2_dTau, cachedval_d3phi0_dDelta3;

	bool cached_phir, cached_dphir_dTau, cached_dphir_dDelta,  cached_d2phir_dTau2, cached_d2phir_dDelta_dTau, cached_d2phir_dDelta2, cached_d3phir_dTau3,  cached_d3phir_dDelta_dTau2, cached_d3phir_dDelta2_dTau, cached_d3phir_dDelta3;
	double cachedval_phir, cachedval_dphir_dTau, cachedval_dphir_dDelta,  cachedval_d2phir_dTau2, cachedval_d2phir_dDelta_dTau, cachedval_d2phir_dDelta2, cachedval_d3phir_dTau3,  cachedval_d3phir_dDelta_dTau2, cachedval_d3phir_dDelta2_dTau, cachedval_d3phir_dDelta3;

	std::string _Fluid;
	
	bool SaturatedL,SaturatedV;
	bool enabled_TTSE_LUT, built_TTSE_LUT;

	// Saturation values
	double rhosatL, rhosatV, psatL, psatV, TsatL, TsatV;

	// Pointers to the Liquid and Vapor classes
	CoolPropStateClass *SatL, *SatV;

	void add_saturation_states(void);
	void remove_saturation_states(void);

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

	// Update using the TTSE lookup tables
	void update_TTSE_LUT(long iInput1, double Value1, long iInput2, double Value2);

	TTSETwoPhaseTableClass TTSESatL;
	TTSETwoPhaseTableClass TTSESatV;
	TTSESinglePhaseTableClass TTSESinglePhase;
public:
	Fluid * pFluid;

	bool flag_SinglePhase, flag_TwoPhase;

	// Bulk values
	double _rho,_T,_p,_Q,_h,_s, tau, delta;

	// Phase flags
	bool TwoPhase, SinglePhase;
	
	// Default Constructor
	CoolPropStateClass(){};

	// Constructor with fluid name
	CoolPropStateClass(std::string FluidName);

	// Constructor with fluid pointer
	CoolPropStateClass(Fluid *pFluid);

	// Destructor to clear SatL and SatV
	~CoolPropStateClass();

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
	void dvdp_dhdp_sat(double T, double *dvdpL, double *dvdpV, double *dhdpL, double *dhdpV, double *);

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

	// ----------------------------------------	
	// TTSE LUT things
	// ----------------------------------------

	/// Enable the TTSE
	void enable_TTSE_LUT(void);
	/// Check if TTSE is enabled
	bool isenabled_TTSE_LUT(void);
	/// Disable the TTSE
	void disable_TTSE_LUT(void);
	/// Over-ride the default size of both of the saturation LUT
	void set_TTSESat_LUT_size(int);
	/// Over-ride the default size of the single-phase LUT
	void set_TTSESinglePhase_LUT_size(int Np, int Nh);
	/// Over-ride the default range of the single-phase LUT
	void set_TTSESinglePhase_LUT_range(double hmin, double hmax, double pmin, double pmax);

	/// Build the TTSE LUT
	bool build_TTSE_LUT();
	/// Interpolate within the TTSE LUT
	double interpolate_in_TTSE_LUT(long iParam, long iInput1, double Input1, long iInput2, double Input2);

	// ----------------------------------------	
	// Derivatives of properties
	// ----------------------------------------

	double dvdp_constT(void);
	double dvdT_constp(void);

	double drhodT_constp(void);
	double drhodh_constp(void);
	double drhodp_consth(void);
	double drhodp_constT(void);
	double d2rhodp2_constT(void);
	double d2rhodTdp(void);
	double d2rhodT2_constp(void);
	
	double dpdrho_constT(void);
	double dpdrho_consth(void);
	double dpdT_constrho(void);
	double dpdT_consth(void);
	double d2pdrho2_constT(void);
	double d2pdrhodT(void);
	double d2pdT2_constrho(void);

	double dhdrho_constT(void);
	double dhdrho_constp(void);
	double dhdT_constrho(void);
	double dhdT_constp(void);
	double dhdp_constT(void);
	double d2hdrho2_constT(void);
	double d2hdrhodT(void);
	double d2hdT2_constrho(void);
	double d2hdT2_constp(void);
	double d2hdp2_constT(void);
	double d2hdTdp(void);

	double dsdrho_constT(void);
	double dsdT_constrho(void);
	double dsdrho_constp(void);
	double dsdT_constp(void);
	double dsdp_constT(void);
	double d2sdrho2_constT(void);
	double d2sdrhodT(void);
	double d2sdT2_constrho(void);
	double d2sdT2_constp(void);
	double d2sdp2_constT(void);
	double d2sdTdp(void);

	// ----------------------------------------	
	// Derivatives along the saturation curve
	// ----------------------------------------
	
	/// Derivative of temperature w.r.t. pressure along saturation curve
	double dTdp_along_sat(void);
	/// Second derivative of temperature w.r.t. pressure along saturation curve
	double d2Tdp2_along_sat(void);
	/// Partial derivative w.r.t. pressure of dTdp along saturation curve
	double ddp_dTdp_along_sat(void);
	/// Partial derivative w.r.t. temperature of dTdp along saturation curve
	double ddT_dTdp_along_sat(void);

	double dhdp_along_sat_vapor(void);
	double dhdp_along_sat_liquid(void);
	double d2hdp2_along_sat_vapor(void);
	double d2hdp2_along_sat_liquid(void);

	double dsdp_along_sat_vapor(void);
	double dsdp_along_sat_liquid(void);
	double d2sdp2_along_sat_vapor(void);
	double d2sdp2_along_sat_liquid(void);

	double drhodp_along_sat_vapor(void);
	double drhodp_along_sat_liquid(void);
	double d2rhodp2_along_sat_vapor(void);
	double d2rhodp2_along_sat_liquid(void);

	double drhodT_along_sat_vapor(void);
	double drhodT_along_sat_liquid(void);

	/// Clear out all the cached values
	void clear_cache(void);

	// ----------------------------------------	
	// Helmholtz Energy Derivatives
	// ----------------------------------------

	double phi0(double tau, double delta);
	double dphi0_dDelta(double tau, double delta);
	double dphi0_dTau(double tau, double delta);
	double d2phi0_dDelta2(double tau, double delta);
	double d2phi0_dDelta_dTau(double tau, double delta);
	double d2phi0_dTau2(double tau, double delta);
	double d3phi0_dDelta3(double tau, double delta);
	double d3phi0_dDelta2_dTau(double tau, double delta);
	double d3phi0_dDelta_dTau2(double tau, double delta);
	double d3phi0_dTau3(double tau, double delta);

	double phir(double tau, double delta);
	double dphir_dDelta(double tau, double delta);
	double dphir_dTau(double tau, double delta);
	double d2phir_dDelta2(double tau, double delta);
	double d2phir_dDelta_dTau(double tau, double delta);
	double d2phir_dTau2(double tau, double delta);
	double d3phir_dDelta3(double tau, double delta);
	double d3phir_dDelta2_dTau(double tau, double delta);
	double d3phir_dDelta_dTau2(double tau, double delta);
	double d3phir_dTau3(double tau, double delta);
};

#endif