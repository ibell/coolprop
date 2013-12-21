/*
 * DataStructures.h
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#ifndef DATASTRUCTURES_H_
#define DATASTRUCTURES_H_

namespace CoolProp {

struct SimpleState
{
	double rho, T, p, h, s;
};

// --------------------------------------------------
// Define some constants that will be used throughout
// --------------------------------------------------

//TODO: Remove this dummy entry.
#define _HUGE 1e13

// These are constants for the input and output parameters
// The structure is taken directly from the AbstractState class.
const enum params {
	// Bulk properties
	iT, irho, ip, iQ, ih, is, icp, icv, ispeed_sound, iisothermal_compressibility, iisobaric_expansion_coefficient,

	// Smoothing functions for density
	idrhodh_constp_smoothed, idrhodp_consth_smoothed, irho_smoothed,

	// Transport properties
	iviscosity, iconductivity, isurface_tension,

	// Derivatives of properties
	idvdp_constT, idvdT_constp,
	// Density
	idrhodh_constp, idrhodp_consth, idrhodp_constT, idrhodT_constp, id2rhodh2_constp,
	id2rhodhdp, id2rhodhdQ, id2rhodp2_constT, id2rhodpdQ, id2rhodT2_constp, id2rhodTdp,
	// Pressure
	idpdrho_consth, idpdrho_constT, idpdT_consth, idpdT_constrho, id2pdrho2_constT,
	id2pdrhodT, id2pdT2_constrho,
	// Enthalpy
	idhdp_constrho, idhdp_constT, idhdrho_constp, idhdrho_constT, idhdT_constp,
	idhdT_constrho, id2hdp2_constT, id2hdrho2_constT, id2hdrhodT, id2hdT2_constp,
	id2hdT2_constrho, id2hdTdp,
	// Entropy
	idsdp_constT, idsdrho_constp, idsdrho_constT, idsdT_constp, idsdT_constrho,
	id2sdp2_constT,	id2sdrho2_constT, id2sdrhodT, id2sdT2_constp, id2sdT2_constrho,
	id2sdTdp,
	// Fundamental derivative of gas dynamics
	ifundamental_derivative_of_gas_dynamics, id2pdv2_consts,

	// Other functions and derivatives
	iB, iC, iZ, idBdT, idCdT, idZdDelta, idZdTau
};

// These are constants for the phases of the fluid
const enum phases {iLiquid, iSupercritical, iGas, iTwoPhase, iUnknown};

// These are constants for the units
const enum unit_constants{UNIT_KPA, UNIT_PA, UNIT_BAR, UNIT_KG_M3, UNIT_KG_L};

// These are constants for the unit systems (currently only SI and KSI are supported)
const enum unit_systems{UNIT_SYSTEM_SI, UNIT_SYSTEM_KSI, UNIT_SYSTEM_KSI_MOLAR, UNIT_SYSTEM_SI_MOLAR};

// These are unit types for the fluid
const enum fluid_types{FLUID_TYPE_PURE, FLUID_TYPE_PSEUDOPURE, FLUID_TYPE_REFPROP, FLUID_TYPE_INCOMPRESSIBLE_LIQUID, FLUID_TYPE_INCOMPRESSIBLE_SOLUTION, FLUID_TYPE_UNDEFINED};

} /* namespace CoolProp */
#endif /* DATASTRUCTURES_H_ */
