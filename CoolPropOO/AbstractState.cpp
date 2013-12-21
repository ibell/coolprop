/*
 * AbstractState.cpp
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#include "AbstractState.h"

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

	return true;
}

} /* namespace CoolProp */
