/*
 * AbstractState.cpp
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#include "AbstractState.h"
#include "math.h"

namespace CoolProp {

AbstractState::AbstractState() {
	// TODO Auto-generated constructor stub
}

AbstractState::~AbstractState() {
	// TODO Auto-generated destructor stub
}

bool AbstractState::clear() {
	// Reset all instances of CachedElement and overwrite
	// the internal double values with -_HUGE
	this->fluid_type = FLUID_TYPE_UNDEFINED;
	this->phase = iUnknown;
	this->forceSinglePhase = false;
	this->forceTwoPhase = false;

	this->critical.T = -_HUGE;
	this->critical.h = -_HUGE;
	this->critical.p = -_HUGE;
	this->critical.rho = -_HUGE;
	this->critical.s = -_HUGE;

	this->reducing.T = -_HUGE;
	this->reducing.h = -_HUGE;
	this->reducing.p = -_HUGE;
	this->reducing.rho = -_HUGE;
	this->reducing.s = -_HUGE;

	/// Bulk values
	this->_rho = -_HUGE;
	this->_T = -_HUGE;
	this->_p = -_HUGE;
	this->_Q = -_HUGE;
	this->tau = -_HUGE;
	this->delta = -_HUGE;
	this->_h.clear();
	this->_s.clear();
	this->_logp.clear();
	this->_logrho.clear();

	/// Smoothing values
	this->rhospline = -_HUGE;
	this->dsplinedp = -_HUGE;
	this->dsplinedh = -_HUGE;

	/// Cached low-level elements for in-place calculation of other properties
	this->_phi0.clear();
	this->_dphi0_dTau.clear();
	this->_dphi0_dDelta.clear();
	this->_d2phi0_dTau2.clear();
	this->_d2phi0_dDelta_dTau.clear();
	this->_d2phi0_dDelta2.clear();
	this->_d3phi0_dTau3.clear();
	this->_d3phi0_dDelta_dTau2.clear();
	this->_d3phi0_dDelta2_dTau.clear();
	this->_d3phi0_dDelta3.clear();
	this->_phir.clear();
	this->_dphir_dTau.clear();
	this->_dphir_dDelta.clear();
	this->_d2phir_dTau2.clear();
	this->_d2phir_dDelta_dTau.clear();
	this->_d2phir_dDelta2.clear();
	this->_d3phir_dTau3.clear();
	this->_d3phir_dDelta_dTau2.clear();
	this->_d3phir_dDelta2_dTau.clear();
	this->_d3phir_dDelta3.clear();

	this->_dphir_dDelta_lim.clear();
	this->_d2phir_dDelta2_lim.clear();
	this->_d2phir_dDelta_dTau_lim.clear();
	this->_d3phir_dDelta2_dTau_lim.clear();

	return true;
}



	virtual double AbstractState::h(void){
		if (!_h) _h = _R*_T*(1.0+tau*(dphi0_dTau()+dphir_dTau())+delta*dphir_dDelta());
		return _h;
	}
	virtual double AbstractState::s(void){
		if (!_s) _s = _R*(tau*(dphi0_dTau()+dphir_dTau())-phi0()-phir());
		return _s;
	}
	virtual double AbstractState::cp(void){
		double c1 = pow(1.0+delta*dphir_dDelta()-delta*tau*d2phir_dDelta_dTau(),2);
		double c2 = (1.0+2.0*delta*dphir_dDelta()+pow(delta,2)*d2phir_dDelta2());
		double val = _R*(-pow(tau,2)*(d2phi0_dTau2()+d2phir_dTau2())+c1/c2);
		return val;
	}
	virtual double AbstractState::cv(void){
		double val = -_R*pow(tau,2)*(d2phi0_dTau2()+d2phir_dTau2());
		return val;
	}
	virtual double AbstractState::speed_sound(void){
		double c1 = pow(tau,2)*(d2phi0_dTau2()+d2phir_dTau2());
		double c2 = (1.0+2.0*delta*dphir_dDelta()+pow(delta,2)*d2phir_dDelta2());
		return sqrt(-c2*_T*cp()/c1);
	}
	virtual double AbstractState::isothermal_compressibility(void){
		return 1.0/(_rho*dpdrho_constT());
	}
	virtual double AbstractState::isobaric_expansion_coefficient(void){
		return -1.0/(_rho*_rho)*drhodT_constp();
	}
//
//
//	// ----------------------------------------
//	// Smoothing functions for density
//	// ----------------------------------------
//	/// A smoothed version of the derivative using a spline curve in the region of x=0 to x=xend
//	virtual double AbstractState::drhodh_constp_smoothed(double xend);
//	/// A smoothed version of the derivative using a spline curve in the region of x=0 to x=xend
//	virtual double AbstractState::drhodp_consth_smoothed(double xend);
//	/// Density corresponding to the smoothed derivatives in the region of x=0 to x=xend
//	virtual void AbstractState::rho_smoothed(double xend, double *rho_spline, double *dsplinedh, double *dsplinedp);
//
//
//	// ----------------------------------------
//	// Transport properties // TODO: Fix it!
//	// ----------------------------------------
//	virtual double AbstractState::viscosity(void);
//	virtual double AbstractState::conductivity(void);
//
//	virtual double AbstractState::surface_tension(void);
//
//
	// ----------------------------------------
	// Derivatives of properties
	// ----------------------------------------
	virtual double AbstractState::dvdp_constT(void){
		return -1/(_rho*_rho)/dpdrho_constT();
	}
	virtual double AbstractState::dvdT_constp(void){
		return -1/(_rho*_rho)*drhodT_constp();
	}

	// Density
	virtual double AbstractState::drhodh_constp(void){
		return 1.0/(dhdrho_constT()-dhdT_constrho()*dpdrho_constT()/dpdT_constrho());
	}
	virtual double AbstractState::drhodp_consth(void){
		return 1.0/(dpdrho_constT()-dpdT_constrho()*dhdrho_constT()/dhdT_constrho());
	}
	virtual double AbstractState::drhodp_constT(void){
		return -dpdT_constrho()/dpdrho_constT();
	}
	virtual double AbstractState::drhodT_constp(void){
		return -dpdT_constrho()/dpdrho_constT();
	}
	virtual double AbstractState::d2rhodh2_constp(void){
		//double A = dpdT_constrho()*dhdrho_constT()-dpdrho_constT()*dhdT_constrho();
		//double dAdT_constrho = d2pdT2_constrho()*dhdrho_constT()+dpdT_constrho()*d2hdrhodT()-d2pdrhodT()*dhdT_constrho()-dpdrho_constT()*d2hdT2_constrho();
		//double dAdrho_constT = d2pdrhodT()*dhdrho_constT()+dpdT_constrho()*d2hdrho2_constT()-d2pdrho2_constT()*dhdT_constrho()-dpdrho_constT()*d2hdrhodT();
		double ddT_drhodh_p_constrho = 1.0/A()*d2pdT2_constrho()-1.0/(A()*A())*dAdT_constrho()*dpdT_constrho();
		double ddrho_drhodh_p_constT = 1.0/A()*d2pdrhodT()      -1.0/(A()*A())*dAdrho_constT()*dpdT_constrho();
		return ddT_drhodh_p_constrho/dhdT_constp()+ddrho_drhodh_p_constT/dhdrho_constp();
	}
	virtual double AbstractState::d2rhodhdp(void){
		//double A = dpdT_constrho()*dhdrho_constT()-dpdrho_constT()*dhdT_constrho();
		//double dAdT_constrho = d2pdT2_constrho()*dhdrho_constT()+dpdT_constrho()*d2hdrhodT()-d2pdrhodT()*dhdT_constrho()-dpdrho_constT()*d2hdT2_constrho();
		//double dAdrho_constT = d2pdrhodT()*dhdrho_constT()+dpdT_constrho()*d2hdrho2_constT()-d2pdrho2_constT()*dhdT_constrho()-dpdrho_constT()*d2hdrhodT();
		double ddT_drhodp_h_constrho = -1.0/A()*d2hdT2_constrho()+1.0/(A()*A())*dAdT_constrho()*dhdT_constrho();
		double ddrho_drhodp_h_constT = -1.0/A()*d2hdrhodT()      +1.0/(A()*A())*dAdrho_constT()*dhdT_constrho();
		return ddT_drhodp_h_constrho/dhdT_constp()+ddrho_drhodp_h_constT/dhdrho_constp();
	}
	//virtual double AbstractState::d2rhodhdQ(void);//TODO: only two-phase, not here
	virtual double AbstractState::d2rhodp2_constT(void){
		return -d2pdrho2_constT()/pow(dpdrho_constT(),3);
	}
	//virtual double AbstractState::d2rhodpdQ(void);//TODO: only two-phase, not here
	virtual double AbstractState::d2rhodT2_constp(void){
		double ddrho_drhodT_p_constT = (dpdT_constrho()*d2pdrho2_constT()-dpdrho_constT()*d2pdrhodT())/pow(dpdrho_constT(),2);
		double ddT_drhodT_p_constrho = (dpdT_constrho()*d2pdrhodT()-dpdrho_constT()*d2pdT2_constrho())/pow(dpdrho_constT(),2);
		return ddT_drhodT_p_constrho+ddrho_drhodT_p_constT*drhodT_constp();
	}
	virtual double AbstractState::d2rhodTdp(void){
		return (dpdT_constrho()*d2pdrho2_constT()-dpdrho_constT()*d2pdrhodT())/pow(dpdrho_constT(),3);
	}

	// Pressure
	virtual double AbstractState::dpdrho_consth(void){
		return dpdrho_constT() - dpdT_constrho()*dhdrho_constT()/dhdT_constrho();
	}
	virtual double AbstractState::dpdrho_constT(void){
		return _R*_T*(1+2*delta*dphir_dDelta()+delta*delta*d2phir_dDelta2());
	}
	virtual double AbstractState::dpdT_consth(void){
		return dpdT_constrho() - dpdrho_constT()*dhdT_constrho()/dhdrho_constT();
	}
	virtual double AbstractState::dpdT_constrho(void){
		return _R*_rho*(1+delta*dphir_dDelta()-delta*tau*d2phir_dDelta_dTau());
	}
	virtual double AbstractState::d2pdrho2_constT(void){
		return _R*_T/_rho*(2*delta*dphir_dDelta()+4*delta*delta*d2phir_dDelta2()+delta*delta*delta*d3phir_dDelta3());
	}
	virtual double AbstractState::d2pdrhodT(void){
		return _R*((1+2*delta*dphir_dDelta()+delta*delta*d2phir_dDelta2())+_T*(2*delta*d2phir_dDelta_dTau()+delta*delta*d3phir_dDelta2_dTau())*(-tau/_T));
	}
	virtual double AbstractState::d2pdT2_constrho(void){
		return _R*_rho*delta*tau*tau/_T*d3phir_dDelta_dTau2();
	}

	// Enthalpy
	virtual double AbstractState::dhdp_constrho(void){
		//	double dphir_dDelta = pFluid->dphir_dDelta(tau,delta);
		//	double d2phir_dDelta_dTau = pFluid->d2phir_dDelta_dTau(tau,delta);
		//	double d2phir_dDelta2 = pFluid->d2phir_dDelta2(tau,delta);
		//	double dpdrho = R*T*(1+2*delta*dphir_dDelta+delta*delta*d2phir_dDelta2);
		//	double dpdT = R*rho*(1+delta*dphir_dDelta-delta*tau*d2phir_dDelta_dTau);
		//	double cp = -tau*tau*R*(pFluid->d2phi0_dTau2(tau,delta)+pFluid->d2phir_dTau2(tau,delta))+T/rho/rho*(dpdT*dpdT)/dpdrho;
		double dpdrho = dpdrho_constT();
		double dpdT   = dpdT_constrho();
		double drhodT = -dpdT/dpdrho;
		return -cp()/dpdrho/drhodT-_T*drhodT*(-1/_rho/_rho)+1/_rho;
	}
	virtual double AbstractState::dhdp_constT(void){
		return dhdrho_constT()/dpdrho_constT();
	}
	virtual double AbstractState::dhdrho_constp(void){
		return dhdrho_constT() - dhdT_constrho()*dpdrho_constT()/dpdT_constrho();
	}
	virtual double AbstractState::dhdrho_constT(void){
		return _T*_R/_rho*(tau*delta*d2phir_dDelta_dTau()+delta*dphir_dDelta()+delta*delta*d2phir_dDelta2());
	}
	virtual double AbstractState::dhdT_constp(void){
		return dhdT_constrho() - dhdrho_constT()*dpdT_constrho()/dpdrho_constT();
	}
	virtual double AbstractState::dhdT_constrho(void){
		return _R*(-tau*tau*(d2phi0_dTau2()+d2phir_dTau2())+1+delta*dphir_dDelta()-delta*tau*d2phir_dDelta_dTau());
	}
	virtual double AbstractState::d2hdp2_constT(void){
		return (d2hdrho2_constT()-dhdp_constT()*d2pdrho2_constT())/pow(dpdrho_constT(),2);
	}
	virtual double AbstractState::d2hdrho2_constT(void){
		return _T*_R/_rho*(tau*delta*d3phir_dDelta2_dTau()+tau*d2phir_dDelta_dTau()+delta*d2phir_dDelta2()+dphir_dDelta()+delta*delta*d3phir_dDelta3()+2*delta*d2phir_dDelta2())/reducing.rho - dhdrho_constT()/_rho;
	}
	virtual double AbstractState::d2hdrhodT(void){
		return _R*(-tau*tau*d3phir_dDelta_dTau2()+delta*d2phir_dDelta2()+dphir_dDelta()-delta*tau*d3phir_dDelta2_dTau()-tau*d2phir_dDelta_dTau())/reducing.rho;
	}
	virtual double AbstractState::d2hdT2_constp(void){
		double ddT_dhdT = d2hdT2_constrho()-1/pow(dpdrho_constT(),2)*(dpdrho_constT()*(dhdrho_constT()*d2pdT2_constrho()+d2hdrhodT()*dpdT_constrho())-dhdrho_constT()*dpdT_constrho()*d2pdrhodT());
		double drho_dhdT = d2hdrhodT()-1/pow(dpdrho_constT(),2)*(dpdrho_constT()*(dhdrho_constT()*d2pdrhodT()+d2hdrho2_constT()*dpdT_constrho())-dhdrho_constT()*dpdT_constrho()*d2pdrho2_constT());
		return ddT_dhdT-drho_dhdT*dpdT_constrho()/dpdrho_constT();
	}
	virtual double AbstractState::d2hdT2_constrho(void){
		return _R*(-tau*tau*(d3phi0_dTau3()+d3phir_dTau3())-2*tau*(d2phi0_dTau2()+d2phir_dTau2())-delta*tau*d3phir_dDelta_dTau2())*(-tau/_T);
	}
	virtual double AbstractState::d2hdTdp(void){
		return 1.0/dpdrho_constT()*(d2hdrhodT()-dhdp_constT()*(drhodT_constp()*d2pdrho2_constT()+d2pdrhodT())+d2hdrho2_constT()*drhodT_constp());
	}

	// Entropy
	virtual double AbstractState::dsdp_constT(void){
		return dsdrho_constT()/dpdrho_constT();
	}
	virtual double AbstractState::dsdrho_constp(void){
		return dsdrho_constT() - dsdT_constrho()*dpdrho_constT()/dpdT_constrho();
	}
	virtual double AbstractState::dsdrho_constT(void){
		return -_R/_rho*(1+delta*dphir_dDelta()-delta*tau*d2phir_dDelta_dTau());
	}
	virtual double AbstractState::dsdT_constp(void){
		return dsdT_constrho() - dsdrho_constT()*dpdT_constrho()/dpdrho_constT();
	}
	virtual double AbstractState::dsdT_constrho(void){
		return -_R*tau*tau/_T*(d2phi0_dTau2()+d2phir_dTau2());
	}
	virtual double AbstractState::d2sdp2_constT(void){
		return (d2sdrho2_constT()-dsdp_constT()*d2pdrho2_constT())/pow(dpdrho_constT(),2);
	}
	virtual double AbstractState::d2sdrho2_constT(void){
		return -_R/_rho*(delta*d2phir_dDelta2()+dphir_dDelta()-tau*delta*d3phir_dDelta2_dTau()-tau*d2phir_dDelta_dTau())/reducing.rho+_R/_rho/_rho*(1+delta*dphir_dDelta()-delta*tau*d2phir_dDelta_dTau());
	}
	virtual double AbstractState::d2sdrhodT(void){
		// d2phi0_dDelta_dTau2() is zero by definition
		return -_R*tau*tau/_T*d3phir_dDelta_dTau2()/reducing.rho;
	}
	virtual double AbstractState::d2sdT2_constp(void){
		double ddT_dsdT  = d2sdT2_constrho()-1/pow(dpdrho_constT(),2)*(dpdrho_constT()*(dsdrho_constT()*d2pdT2_constrho()+d2sdrhodT()*dpdT_constrho())-dsdrho_constT()*dpdT_constrho()*d2pdrhodT());
		double drho_dsdT = d2sdrhodT()-1/pow(dpdrho_constT(),2)*(dpdrho_constT()*(dsdrho_constT()*d2pdrhodT()+d2sdrho2_constT()*dpdT_constrho())-dsdrho_constT()*dpdT_constrho()*d2pdrho2_constT());
		return ddT_dsdT-drho_dsdT*dpdT_constrho()/dpdrho_constT();
	}
	virtual double AbstractState::d2sdT2_constrho(void){
		return -_R/_T*(tau*tau*(d3phi0_dTau3()+d3phir_dTau3())+2*tau*(d2phi0_dTau2()+d2phir_dTau2()))*(-tau/_T)+_R*tau*tau/_T/_T*(d2phi0_dTau2()+d2phir_dTau2());
	}
	virtual double AbstractState::d2sdTdp(void){
		return 1.0/dpdrho_constT()*(d2sdrhodT()-dsdp_constT()*(drhodT_constp()*d2pdrho2_constT()+d2pdrhodT())+d2sdrho2_constT()*drhodT_constp());
	}

	// Fundamental derivative of gas dynamics
	virtual double AbstractState::fundamental_derivative_of_gas_dynamics(void){
		return d2pdv2_consts()/pow(speed_sound(),2)/2.0/pow(_rho,3);
	}
	virtual double AbstractState::d2pdv2_consts(void){
		double cv = this->cv();
		// Convert each derivative in terms of volume rather than density
		// Note drhodv = -rho^2
		double d2pdv2_constT = _rho*_rho*_rho*_rho*d2pdrho2_constT()+2.0*_rho*_rho*_rho*dpdrho_constT();
		double dpdT_constv = dpdT_constrho();
		double d2pdvdT = -_rho*_rho*d2pdrhodT();
		double d2pdT2_constv = d2pdT2_constrho();
		double dcv_dT_constv = _R*tau/_T*(2.0*tau*(d2phi0_dTau2()+d2phir_dTau2())+tau*tau*(d3phi0_dTau3()+d3phir_dTau3()));
		double LAMBDA1 = d2pdv2_constT;
		double LAMBDA2 = -3.0*_T/cv*dpdT_constv*d2pdvdT;
		double LAMBDA3 = +pow(_T/cv*dpdT_constv,2)*(3.0*d2pdT2_constv+1/_T*dpdT_constv*(1.0-_T/cv*dcv_dT_constv));
		return LAMBDA1 + LAMBDA2 + LAMBDA3;
	}

	// Other functions and derivatives
	virtual double AbstractState::A(void){
		return dpdT_constrho()*dhdrho_constT()-dpdrho_constT()*dhdT_constrho();
	}
	virtual double AbstractState::B(void){
		// given by B*rhoc=lim(delta --> 0) [dphir_ddelta(tau)]
		return 1.0/reducing.rho*dphir_dDelta_lim();
	}
	virtual double AbstractState::C(void){
		// given by C*rhoc^2=lim(delta --> 0) [d2phir_dDelta2(tau)]
		return 1.0/(reducing.rho*reducing.rho)*d2phir_dDelta2_lim();
	}
	virtual double AbstractState::Z(void){
		return 1+delta*dphir_dDelta();
	}
	virtual double AbstractState::dAdT_constrho(void){
		return d2pdT2_constrho()*dhdrho_constT()+dpdT_constrho()*d2hdrhodT()-d2pdrhodT()*dhdT_constrho()-dpdrho_constT()*d2hdT2_constrho();
	}
	virtual double AbstractState::dAdrho_constT(void){
		return d2pdrhodT()*dhdrho_constT()+dpdT_constrho()*d2hdrho2_constT()-d2pdrho2_constT()*dhdT_constrho()-dpdrho_constT()*d2hdrhodT();
	}
	// TODO: Add constXX qualifier
	virtual double AbstractState::dBdT(void){
		return 1.0/reducing.rho*d2phir_dDelta_dTau_lim()*-reducing.T/_T/_T;
	}
	virtual double AbstractState::dCdT(void){
		return 1.0/(reducing.rho*reducing.rho)*d3phir_dDelta2_dTau_lim()*-reducing.T/_T/_T;
	}
	virtual double AbstractState::dZdDelta(void){
		return delta*d2phir_dDelta2()+dphir_dDelta();
	}
	virtual double AbstractState::dZdTau(void){
		return delta*d2phir_dDelta_dTau();
	}

//	// ----------------------------------------
//	// Helmholtz energy and derivatives
//	// ----------------------------------------
//	virtual double AbstractState::phi0(void);
//	virtual double AbstractState::dphi0_dDelta(void);
//	virtual double AbstractState::dphi0_dTau(void);
//	virtual double AbstractState::d2phi0_dDelta2(void);
//	virtual double AbstractState::d2phi0_dDelta_dTau(void);
//	virtual double AbstractState::d2phi0_dTau2(void);
//	virtual double AbstractState::d3phi0_dDelta3(void);
//	virtual double AbstractState::d3phi0_dDelta2_dTau(void);
//	virtual double AbstractState::d3phi0_dDelta_dTau2(void);
//	virtual double AbstractState::d3phi0_dTau3(void);
//
//	virtual double AbstractState::phir(void);
//	virtual double AbstractState::dphir_dDelta(void);
//	virtual double AbstractState::dphir_dTau(void);
//	virtual double AbstractState::d2phir_dDelta2(void);
//	virtual double AbstractState::d2phir_dDelta_dTau(void);
//	virtual double AbstractState::d2phir_dTau2(void);
//	virtual double AbstractState::d3phir_dDelta3(void);
//	virtual double AbstractState::d3phir_dDelta2_dTau(void);
//	virtual double AbstractState::d3phir_dDelta_dTau2(void);
//	virtual double AbstractState::d3phir_dTau3(void);
//
//	// TODO: Add call back to calculator;
//	virtual double AbstractState::dphir_dDelta_lim(void);
//	virtual double AbstractState::d2phir_dDelta2_lim(void);
//	virtual double AbstractState::d2phir_dDelta_dTau_lim(void);
//	virtual double AbstractState::d3phir_dDelta2_dTau_lim(void);

} /* namespace CoolProp */
