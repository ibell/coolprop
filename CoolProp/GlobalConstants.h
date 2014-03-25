
#ifndef GLOBALCONSTANTS_H
#define GLOBALCONSTANTS_H

/* --------------------------------------------------
// Define some constants that will be used throughout
// --------------------------------------------------
*/

/* These are constants for the input and output parameters */
enum params {
	/* Properties and others */
	iB, iT, iP, iD, iC, iC0, iO, iU, iH, iS, iA, iG, iQ, iV, iL, iTfreeze, iPsat, iI, iDpdT, iDrhodT_p, iCritSplineT, iPrandtl,
	/* Constants */
	iMM, iTmax, iTmin, iAccentric, iDipole, iODP, iGWP20, iGWP100, iGWP500,
	/* Reduced quantities; triple point and critical point */
	iRhoreduce, iTreduce, iPtriple, iTtriple, iHcrit, iPcrit, iRhocrit, iScrit, iTcrit,
	/* Phase identifiers */
	iPhase, iPHASE_LIQUID, iPHASE_GAS, iPHASE_SUPERCRITICAL, iPHASE_TWOPHASE,
    /* Ancillary equations */
    iRhosatLanc, iRhosatVanc, iPsatLanc, iPsatVanc,
	/* Derivatives */
	iDERdh_dp__rho, iDERdh_dp__v, iDERZ, iDERdZ_dDelta, iDERdZ_dTau, iDERB, iDERdB_dT, iDERC, iDERdC_dT, iDERphir,
	iDERdphir_dTau, iDERdphir_dDelta, iDERd2phir_dTau2, iDERd2phir_dDelta2, iDERd2phir_dDelta_dTau,	iDERd3phir_dDelta3,
	iDERd3phir_dDelta2_dTau, iDERd3phir_dDelta_dTau2, iDERd3phir_dTau3, iDERphi0, iDERdphi0_dTau, iDERd2phi0_dTau2,
	iDERdphi0_dDelta, iDERd2phi0_dDelta2, iDERd2phi0_dDelta_dTau, iDERd3phi0_dTau3, iDERdp_dT__rho,
	iDERdp_drho__T, iDERdh_dT__rho, iDERdh_drho__T, iDERdrho_dT__p, iDERdrho_dh__p, iDERdh_dp__T,
	iDERdrho_dp__h, iDERrho_smoothed, iDERdrho_smoothed_dh, iDERdrho_smoothed_dp, iDERdrhodh_constp_smoothed,
	iDERdrhodp_consth_smoothed, iDERIsothermalCompressibility,
};

/* These are constants for the phases of the fluid */
enum phases {iLiquid, iSupercritical, iGas, iTwoPhase};

/* These are constants for the units */
enum unit_constants{UNIT_KPA, UNIT_PA, UNIT_BAR, UNIT_KG_M3, UNIT_KG_L};

/* These are constants for the unit systems (currently only SI and KSI are supported)  */
enum unit_systems{UNIT_SYSTEM_SI, UNIT_SYSTEM_KSI, UNIT_SYSTEM_KSI_MOLAR, UNIT_SYSTEM_SI_MOLAR};

/* These are unit types for the fluid */
enum fluid_types{FLUID_TYPE_PURE, FLUID_TYPE_PSEUDOPURE, FLUID_TYPE_REFPROP, FLUID_TYPE_INCOMPRESSIBLE_LIQUID, FLUID_TYPE_INCOMPRESSIBLE_SOLUTION};

#endif
