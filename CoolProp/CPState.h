
#ifndef CPSTATE_H
#define CPSTATE_H

#include <iostream>
#include "FluidClass.h"
#include "CoolProp.h"
#include "TTSE.h"

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
	
	bool SaturatedL,SaturatedV,_noSatLSatV;

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

	// To be used to update internal variables if you know that your parameters are P,D
	void update_prho(long iInput1, double Value1, long iInput2, double Value2);

	// To be used to update internal variables if you know that your parameters are H,S
	void update_hs(long iInput1, double Value1, long iInput2, double Value2);

	// Update using the TTSE lookup tables
	void update_TTSE_LUT(long iInput1, double Value1, long iInput2, double Value2);

	// Check whether the quality corresponds to saturated liquid or vapor
	void check_saturated_quality(double Q);

	/// Check whether within the TTSE range
	bool within_TTSE_range(long iInput1, double Value1, long iInput2, double Value2);

public:
	Fluid * pFluid;

	/// Temporarily set a flag to tell that the next call to update should be special 
	/// cased as though it is single-phase even if it isn't
	bool flag_SinglePhase;

	/// Temporarily set a flag to tell that the next call to update should be special 
	/// cased as though it is two-phase even if it isn't
	bool flag_TwoPhase;

	// Bulk values
	double _rho,_T,_p,_Q,_h,_s,_logp, _logrho, tau, delta;

	// Phase flags
	bool TwoPhase, SinglePhase, s_cached, h_cached;
	
	// Default Constructor
	CoolPropStateClass(){SatL = NULL; SatV = NULL;};

	// Constructor with fluid name
	CoolPropStateClass(std::string FluidName);

	// Constructor with fluid pointer
	CoolPropStateClass(Fluid *pFluid);

	// Destructor to clear SatL and SatV
	~CoolPropStateClass();

	/// Stop it from adding the SatL and SatV class pointers
	void no_SatLSatV(void){_noSatLSatV = true;};

	// Property updater
	// Uses the indices in CoolProp for the input parameters
	void update(long iInput1, double Value1, long iInput2, double Value2);

	// Returns an output based on the key provided
	// where iInput is one of iT,iP,iH,iS,....
	double keyed_output(long iInput);

	// Property accessors for saturation parameters directly
	// These all are calculated every time if the state is saturated or two-phase
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
	double cpL(void);
	double cpV(void);
	double viscL(void);
	double viscV(void);
	double condL(void);
	double condV(void);

	// The phase as an integer flag
	long phase(void);

	// Bulk properties accessors - temperature and density are directly calculated every time
	// All other parameters are calculated on an as-needed basis
	// If single-phase, just plug into the EOS, otherwise need to do two-phase analysis
	double T(void){return _T;};
	double rho(void){return _rho;};
	double p(void){return _p;};
	double Q(void){return _Q;};
	double h(void);
	double s(void);
	double cp(void);
	double cv(void);
	double speed_sound(void);
	double isothermal_compressibility(void);
	double isobaric_expansion_coefficient(void);
	double drhodh_constp(void);
	double drhodp_consth(void);
	/// A smoothed version of the derivative using a spline curve in the region of x=0 to x=xend
	double drhodh_constp_smoothed(double xend);
	/// A smoothed version of the derivative using a spline curve in the region of x=0 to x=xend
	double drhodp_consth_smoothed(double xend);
	/// Density corresponding to the smoothed derivatives in the region of x=0 to x=xend
	void rho_smoothed(double xend, double *rho_spline, double *dsplinedh, double *dsplinedp);

	double viscosity(void);
	double conductivity(void);

	double surface_tension(void);

	// ----------------------------------------	
	// TTSE LUT things
	// ----------------------------------------

	/// Enable the TTSE
	void enable_TTSE_LUT(void);
	/// Check if TTSE is enabled
	bool isenabled_TTSE_LUT(void);
	/// Disable the TTSE
	void disable_TTSE_LUT(void);
	/// Enable the writing of TTSE tables to file
	void enable_TTSE_LUT_writing(void);
	/// Check if the writing of TTSE tables to file is enabled
	bool isenabled_TTSE_LUT_writing(void);
	/// Disable the writing of TTSE tables to file
	void disable_TTSE_LUT_writing(void);
	/// Over-ride the default size of both of the saturation LUT
	void set_TTSESat_LUT_size(int N);
	/// Over-ride the default size of the single-phase LUT
	void set_TTSESinglePhase_LUT_size(int Np, int Nh);
	/// Over-ride the default range of the single-phase LUT
	void set_TTSESinglePhase_LUT_range(double hmin, double hmax, double pmin, double pmax);
	/// Get the current range of the single-phase LUT
	void get_TTSESinglePhase_LUT_range(double *hmin, double *hmax, double *pmin, double *pmax);

	/// Evaluate the B term from TTSE method
	double B_TTSE(double _p, double _h);
	/// Evaluate the D term from TTSE method
	double D_TTSE(double _p, double _h);
	/// Get the ratio directly which is just a bit faster
	double B_over_D_TTSE(double _p, double _h);

	/// Interpolate within the TTSE LUT
	double interpolate_in_TTSE_LUT(long iParam, long iInput1, double Input1, long iInput2, double Input2);

	// ----------------------------------------	
	// Derivatives of properties
	// ----------------------------------------

	double dvdp_constT(void);
	double dvdT_constp(void);

	double drhodT_constp(void);
	double drhodp_constT(void);
	double d2rhodp2_constT(void);
	double d2rhodTdp(void);
	double d2rhodT2_constp(void);
	double d2rhodhdQ(void);
	double d2rhodpdQ(void);
	double d2rhodhdp(void);
	double d2rhodh2_constp(void);
	
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

	// ----------------------------------------	
	// Helmholtz Energy Derivatives
	// ----------------------------------------
	
	/// Clear out all the cached values
	void clear_cache(void);

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
