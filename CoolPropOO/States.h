/*
 * StateImpl.h
 *
 *  Created on: 23 Dec 2013
 *      Author: jowr
 */

#ifndef STATES_H_
#define STATES_H_

#include "AbstractState.h"

namespace CoolProp {

class SinglePhaseState: public CoolProp::AbstractState {
protected:

	AbstractFluid * fluid;

public:
	// ----------------------------------------
	// Transport properties
	// ----------------------------------------
	double viscosity(void);
	double conductivity(void);
	double surface_tension(void);

protected:
	// ----------------------------------------
	// Helmholtz energy and derivatives
	// ----------------------------------------
	double phi0(void);
	double dphi0_dDelta(void);
	double dphi0_dTau(void);
	double d2phi0_dDelta2(void);
	double d2phi0_dDelta_dTau(void);
	double d2phi0_dTau2(void);
	double d3phi0_dDelta3(void);
	double d3phi0_dDelta2_dTau(void);
	double d3phi0_dDelta_dTau2(void);
	double d3phi0_dTau3(void);

	double phir(void);
	double dphir_dDelta(void);
	double dphir_dTau(void);
	double d2phir_dDelta2(void);
	double d2phir_dDelta_dTau(void);
	double d2phir_dTau2(void);
	double d3phir_dDelta3(void);
	double d3phir_dDelta2_dTau(void);
	double d3phir_dDelta_dTau2(void);
	double d3phir_dTau3(void);

	double dphir_dDelta_lim(void);
	double d2phir_dDelta2_lim(void);
	double d2phir_dDelta_dTau_lim(void);
	double d3phir_dDelta2_dTau_lim(void);

public:
	SinglePhaseState();
	virtual ~SinglePhaseState();
};

class TwoPhaseState: public CoolProp::SinglePhaseState {
protected:

	SinglePhaseState *satL, *satV;

public:
	// ----------------------------------------
	// Smoothing functions for density // TODO: Not implemented yet
	// ----------------------------------------
	/// A smoothed version of the derivative using a spline curve in the region of x=0 to x=xend
	virtual double drhodh_constp_smoothed(double xend);
	/// A smoothed version of the derivative using a spline curve in the region of x=0 to x=xend
	virtual double drhodp_consth_smoothed(double xend);
	/// Density corresponding to the smoothed derivatives in the region of x=0 to x=xend
	virtual void rho_smoothed(double xend, double *rho_spline, double *dsplinedh, double *dsplinedp);

	// ----------------------------------------
	// Derivatives along the saturation curve
	// ----------------------------------------
	/// Derivative of temperature w.r.t. pressure along saturation curve
	virtual double dTdp_along_sat(void);
	/// Second derivative of temperature w.r.t. pressure along saturation curve
	virtual double d2Tdp2_along_sat(void);
	/// Partial derivative w.r.t. pressure of dTdp along saturation curve
	virtual double ddp_dTdp_along_sat(void);
	/// Partial derivative w.r.t. temperature of dTdp along saturation curve
	virtual double ddT_dTdp_along_sat(void);

	virtual double dhdp_along_sat_vapor(void);
	virtual double dhdp_along_sat_liquid(void);
	virtual double d2hdp2_along_sat_vapor(void);
	virtual double d2hdp2_along_sat_liquid(void);

	virtual double dsdp_along_sat_vapor(void);
	virtual double dsdp_along_sat_liquid(void);
	virtual double d2sdp2_along_sat_vapor(void);
	virtual double d2sdp2_along_sat_liquid(void);

	virtual double drhodp_along_sat_vapor(void);
	virtual double drhodp_along_sat_liquid(void);
	virtual double d2rhodp2_along_sat_vapor(void);
	virtual double d2rhodp2_along_sat_liquid(void);

	/*virtual double dsdT_along_sat_vapor(void);
	virtual double dsdT_along_sat_liquid(void);

	virtual double dhdT_along_sat_vapor(void);
	virtual double dhdT_along_sat_liquid(void);*/

	virtual double drhodT_along_sat_vapor(void);
	virtual double drhodT_along_sat_liquid(void);

public:
	TwoPhaseState();
	virtual ~TwoPhaseState();
};

class CoolPropState: public CoolProp::TwoPhaseState {
public:
	CoolPropState();
	virtual ~CoolPropState();
};

class RefpropState: public CoolProp::TwoPhaseState {
public:
	RefpropState();
	virtual ~RefpropState();
};

class IncompressibleState: public CoolProp::SinglePhaseState {
protected:
	Incompressible * fluid;
public:
	IncompressibleState();
	virtual ~IncompressibleState();
};



} /* namespace CoolProp */
#endif /* STATES_H_ */
