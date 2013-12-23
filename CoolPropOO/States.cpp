/*
 * StateImpl.cpp
 *
 *  Created on: 23 Dec 2013
 *      Author: jowr
 */

#include "States.h"

namespace CoolProp {

SinglePhaseState::SinglePhaseState() {
	// TODO Auto-generated constructor stub

}

SinglePhaseState::~SinglePhaseState() {
	// TODO Auto-generated destructor stub
}

// ----------------------------------------
// Transport properties // TODO: Implement them
// ----------------------------------------
double SinglePhaseState::viscosity(void){return -_HUGE;}
double SinglePhaseState::conductivity(void){return -_HUGE;}
double SinglePhaseState::surface_tension(void){return -_HUGE;}

// ----------------------------------------
// Helmholtz energy and derivatives
// ----------------------------------------
double SinglePhaseState::phi0(void){
	if (!_phi0) _phi0 = fluid->phi0(tau,delta);
	return _phi0;
}
double SinglePhaseState::dphi0_dDelta(void);
double SinglePhaseState::dphi0_dTau(void);
double SinglePhaseState::d2phi0_dDelta2(void);
double SinglePhaseState::d2phi0_dDelta_dTau(void);
double SinglePhaseState::d2phi0_dTau2(void);
double SinglePhaseState::d3phi0_dDelta3(void);
double SinglePhaseState::d3phi0_dDelta2_dTau(void);
double SinglePhaseState::d3phi0_dDelta_dTau2(void);
double SinglePhaseState::d3phi0_dTau3(void);

double SinglePhaseState::phir(void);
double SinglePhaseState::dphir_dDelta(void);
double SinglePhaseState::dphir_dTau(void);
double SinglePhaseState::d2phir_dDelta2(void);
double SinglePhaseState::d2phir_dDelta_dTau(void);
double SinglePhaseState::d2phir_dTau2(void);
double SinglePhaseState::d3phir_dDelta3(void);
double SinglePhaseState::d3phir_dDelta2_dTau(void);
double SinglePhaseState::d3phir_dDelta_dTau2(void);
double SinglePhaseState::d3phir_dTau3(void);

double SinglePhaseState::dphir_dDelta_lim(void);
double SinglePhaseState::d2phir_dDelta2_lim(void);
double SinglePhaseState::d2phir_dDelta_dTau_lim(void);
double SinglePhaseState::d3phir_dDelta2_dTau_lim(void);


// ----------------------------------------
// Overwrite some functions to make them aware
// of the two-phase region // TODO: Not fully implemented yet
// ----------------------------------------


// ----------------------------------------
// Smoothing functions for density // TODO: Not implemented yet
// ----------------------------------------
/// A smoothed version of the derivative using a spline curve in the region of x=0 to x=xend
virtual double TwoPhaseState::drhodh_constp_smoothed(double xend){return -_HUGE;}
/// A smoothed version of the derivative using a spline curve in the region of x=0 to x=xend
virtual double TwoPhaseState::drhodp_consth_smoothed(double xend){return -_HUGE;}
/// Density corresponding to the smoothed derivatives in the region of x=0 to x=xend
virtual void TwoPhaseState::rho_smoothed(double xend, double *rho_spline, double *dsplinedh, double *dsplinedp){}

// ----------------------------------------
// Derivatives along the saturation curve
// ----------------------------------------
/// Derivative of temperature w.r.t. pressure along saturation curve
virtual double TwoPhaseState::dTdp_along_sat(void){
	return _T*(1.0/satV->rho()-1.0/satL->rho())/(satV->h()-satL->h());
}
/// Second derivative of temperature w.r.t. pressure along saturation curve
virtual double TwoPhaseState::d2Tdp2_along_sat(void){
	return ddp_dTdp_along_sat()+ddT_dTdp_along_sat()*dTdp_along_sat();
}
/// Partial derivative w.r.t. pressure of dTdp along saturation curve
virtual double TwoPhaseState::ddp_dTdp_along_sat(void){
	return 1.0/(satV->h()-satL->h())*(_T*(satV->dvdp_constT()-satL->dvdp_constT())-dTdp_along_sat()*(satV->dhdp_constT()-satL->dhdp_constT()));
}
/// Partial derivative w.r.t. temperature of dTdp along saturation curve
virtual double TwoPhaseState::ddT_dTdp_along_sat(void){
	return 1.0/(satV->h()-satL->h())*(_T*(satV->dvdT_constp()-satL->dvdT_constp())-dTdp_along_sat()*(satV->dhdT_constp()-satL->dhdT_constp())+(1.0/satV->rho()-1.0/satL->rho()));
}

virtual double TwoPhaseState::dhdp_along_sat_vapor(void){
	return satV->dhdp_constT()+satV->dhdT_constp()*dTdp_along_sat();
}
virtual double TwoPhaseState::dhdp_along_sat_liquid(void){
	return satL->dhdp_constT()+satL->dhdT_constp()*dTdp_along_sat();
}
virtual double TwoPhaseState::d2hdp2_along_sat_vapor(void){
	double ddp_dhdpsigmaV = satV->d2hdp2_constT()+satV->dhdT_constp()*ddp_dTdp_along_sat()+satV->d2hdTdp()      *dTdp_along_sat();
	double ddT_dhdpsigmaV = satV->d2hdTdp()      +satV->dhdT_constp()*ddT_dTdp_along_sat()+satV->d2hdT2_constp()*dTdp_along_sat();
	return ddp_dhdpsigmaV+ddT_dhdpsigmaV*dTdp_along_sat();
}
virtual double TwoPhaseState::d2hdp2_along_sat_liquid(void){
	double ddp_dhdpsigmaL = satL->d2hdp2_constT()+satL->dhdT_constp()*ddp_dTdp_along_sat()+satL->d2hdTdp()      *dTdp_along_sat();
	double ddT_dhdpsigmaL = satL->d2hdTdp()+      satL->dhdT_constp()*ddT_dTdp_along_sat()+satL->d2hdT2_constp()*dTdp_along_sat();
	return ddp_dhdpsigmaL+ddT_dhdpsigmaL*dTdp_along_sat();
}

virtual double TwoPhaseState::dsdp_along_sat_vapor(void){
	return satV->dsdp_constT()+satV->dsdT_constp()*dTdp_along_sat();
}
virtual double TwoPhaseState::dsdp_along_sat_liquid(void){
	return satL->dsdp_constT()+satL->dsdT_constp()*dTdp_along_sat();
}
virtual double TwoPhaseState::d2sdp2_along_sat_vapor(void){
	double ddp_dsdpsigmaV = satV->d2sdp2_constT()+satV->dsdT_constp()*ddp_dTdp_along_sat()+satV->d2sdTdp()      *dTdp_along_sat();
	double ddT_dsdpsigmaV = satV->d2sdTdp()      +satV->dsdT_constp()*ddT_dTdp_along_sat()+satV->d2sdT2_constp()*dTdp_along_sat();
	return ddp_dsdpsigmaV+ddT_dsdpsigmaV*dTdp_along_sat();
}
virtual double TwoPhaseState::d2sdp2_along_sat_liquid(void){
	double ddp_dsdpsigmaL = satL->d2sdp2_constT()+satL->dsdT_constp()*ddp_dTdp_along_sat()+satL->d2sdTdp()      *dTdp_along_sat();
	double ddT_dsdpsigmaL = satL->d2sdTdp()      +satL->dsdT_constp()*ddT_dTdp_along_sat()+satL->d2sdT2_constp()*dTdp_along_sat();
	return ddp_dsdpsigmaL+ddT_dsdpsigmaL*dTdp_along_sat();
}

virtual double TwoPhaseState::drhodp_along_sat_vapor(void){
	return satV->drhodp_constT()+satV->drhodT_constp()*dTdp_along_sat();
}
virtual double TwoPhaseState::drhodp_along_sat_liquid(void){
	return satL->drhodp_constT()+satL->drhodT_constp()*dTdp_along_sat();
}
virtual double TwoPhaseState::d2rhodp2_along_sat_vapor(void){
	double ddp_drhodpsigmaV = satV->d2rhodp2_constT()+satV->drhodT_constp()*ddp_dTdp_along_sat()+satV->d2rhodTdp()      *dTdp_along_sat();
	double ddT_drhodpsigmaV = satV->d2rhodTdp()      +satV->drhodT_constp()*ddT_dTdp_along_sat()+satV->d2rhodT2_constp()*dTdp_along_sat();
	return ddp_drhodpsigmaV+ddT_drhodpsigmaV*dTdp_along_sat();
}
virtual double TwoPhaseState::d2rhodp2_along_sat_liquid(void){
	double ddp_drhodpsigmaL = satL->d2rhodp2_constT()+satL->drhodT_constp()*ddp_dTdp_along_sat()+satL->d2rhodTdp()      *dTdp_along_sat();
	double ddT_drhodpsigmaL = satL->d2rhodTdp()      +satL->drhodT_constp()*ddT_dTdp_along_sat()+satL->d2rhodT2_constp()*dTdp_along_sat();
	return ddp_drhodpsigmaL+ddT_drhodpsigmaL*dTdp_along_sat();
}

/*virtual double TwoPhaseState::dsdT_along_sat_vapor(void);
virtual double TwoPhaseState::dsdT_along_sat_liquid(void);

virtual double TwoPhaseState::dhdT_along_sat_vapor(void);
virtual double TwoPhaseState::dhdT_along_sat_liquid(void);*/

virtual double TwoPhaseState::drhodT_along_sat_vapor(void){
	return satV->drhodT_constp()+satV->drhodp_constT()/dTdp_along_sat();
}
virtual double TwoPhaseState::drhodT_along_sat_liquid(void){
	return satL->drhodT_constp()+satL->drhodp_constT()/dTdp_along_sat();
}

TwoPhaseState::TwoPhaseState() {
	// TODO Auto-generated constructor stub

}

TwoPhaseState::~TwoPhaseState() {
	// TODO Auto-generated destructor stub
}

} /* namespace CoolProp */
