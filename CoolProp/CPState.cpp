
#include "CPExceptions.h"
#include "CPState.h"

// Constructor with fluid name
CoolPropStateClass::CoolPropStateClass(std::string Fluid){
	// Try to get the index of the fluid
	long iFluid = get_Fluid_index(Fluid);
	// If iFluid is greater than -1, it is a CoolProp Fluid, otherwise not
	if (iFluid > -1)
	{
		// Get a pointer to the fluid object
		pFluid = get_fluid(iFluid);
	}
	else
	{
		throw ValueError("Bad Fluid name - not a CoolProp fluid");
	}
}

// Constructor with pointer to fluid
CoolPropStateClass::CoolPropStateClass(Fluid * pFluid){
	this->pFluid = pFluid;
}

bool match_pair(long iI1, long iI2, long I1, long I2)
{
	return ((iI1 == I1 || iI1 == I2) && (iI2 == I1 || iI2 == I2) && iI1!=iI2);
}
void sort_pair(long *iInput1, double *Value1, long *iInput2, double *Value2, long I1, long I2)
{
	if (!(*iInput1 == I1) || !(*iInput2 == I2)){
		std::swap(*iInput1,*iInput2);
		std::swap(*Value1,*Value2);
	}
}
void CoolPropStateClass::check_saturated_quality(double Q){
	double mach_eps = 10*DBL_EPSILON;

	if (fabs(Q-1) < mach_eps){
		SaturatedL = true; SaturatedV = false;
	}
	else if (fabs(Q) < mach_eps){
		SaturatedL = false; SaturatedV = true;
	}
	else{
		SaturatedL = false; SaturatedV = false;
	}
}
// Main updater function
void CoolPropStateClass::update(long iInput1, double Value1, long iInput2, double Value2){
	/* Options for inputs (in either order) are:
	|  T,P
	|  T,D
	|  H,P
	|  S,P
	|  P,Q
	|  T,Q
	|
	*/

	// Reset all the internal variables to _HUGE
	_T = _HUGE;
	_p = _HUGE;
	_h = _HUGE;
	_s = _HUGE;
	_rho = _HUGE;
	_Q = _HUGE;

	// If flag_SinglePhase is true, it will always assume that it is not in the two-phase region
	// Can be over-written by changing the flag to true
	flag_SinglePhase = false;

	// Don't know if it is single phase or not, so assume it isn't
	SinglePhase = false;
	
	// If the inputs are P,Q or T,Q , it is guaranteed to require a call to the saturation routine
	if (match_pair(iInput1,iInput2,iP,iQ) || match_pair(iInput1,iInput2,iT,iQ)){
		update_twophase(iInput1,Value1,iInput2,Value2);
	}
	else if (match_pair(iInput1,iInput2,iT,iD)){
		update_Trho(iInput1,Value1,iInput2,Value2);
	}
	else if (match_pair(iInput1,iInput2,iT,iP)){
		update_Tp(iInput1,Value1,iInput2,Value2);
	}
	else if (match_pair(iInput1,iInput2,iP,iH)){
		update_ph(iInput1,Value1,iInput2,Value2);
	}
	else if (match_pair(iInput1,iInput2,iP,iS)){
		update_ps(iInput1,Value1,iInput2,Value2);
	}
	else
	{
		throw ValueError(format("Sorry your inputs didn't work"));
	}
}
void CoolPropStateClass::update_twophase(long iInput1, double Value1, long iInput2, double Value2)
{
	// This function handles setting internal variables when the state is known to be saturated
	// Either T,Q or P,Q are given
	double Q;
	
	SinglePhase = false;
	TwoPhase = true;

	if (iInput1 == iQ){
		Q = Value1;
		sort_pair(&iInput1,&Value1,&iInput2,&Value2,iInput2,iQ);
	}
	else{
		Q = Value2;
		sort_pair(&iInput1,&Value1,&iInput2,&Value2,iInput1,iQ);
	}

	// Check whether saturated liquid or vapor
	check_saturated_quality(Q);

	if (match_pair(iInput1,iInput2,iP,iQ)){
		// Sort so they are in the order P, Q
		sort_pair(&iInput1,&Value1,&iInput2,&Value2,iP,iQ);
		// Carry out the saturation call to get the temperature and density for each phases
		if (pFluid->pure()){
			if (SaturationLUTStatus()){
				TsatL = pFluid->ApplySaturationLUT(pFluid->SatLUT.iT,pFluid->SatLUT.iP,Value1);
				rhosatL = pFluid->ApplySaturationLUT(pFluid->SatLUT.iDL,pFluid->SatLUT.iP,Value1);
				rhosatV = pFluid->ApplySaturationLUT(pFluid->SatLUT.iDV,pFluid->SatLUT.iP,Value1);
				TsatV = TsatL;
			}
			else{
				pFluid->TsatP_Pure(Value1, &TsatL, &rhosatL, &rhosatV);
				TsatV = TsatL;
			}
		}
		else{
			TsatL = pFluid->Tsat_anc(Value1,0);
			TsatV = pFluid->Tsat_anc(Value1,1);
			psatL = Value1;
			psatV = Value1;
			// Saturation densities
			rhosatL = pFluid->density_Tp(TsatL, psatL, pFluid->rhosatL(TsatL));
			rhosatV = pFluid->density_Tp(TsatV, psatV, pFluid->rhosatV(TsatV));
		}
	}
	else{
		// Sort so they are in the order T, Q
		sort_pair(&iInput1,&Value1,&iInput2,&Value2,iT,iQ);
		// Carry out the saturation call to get the temperature and density for each phases
		if (pFluid->pure()){
			pFluid->saturation(Value1,SaturationLUTStatus(),&psatL,&psatV,&rhosatL,&rhosatV);
			TsatL = Value1;
			TsatV = Value1;
		}
		else{
			TsatL = Value1;
			TsatV = Value1;
			// Saturation pressures
			psatL = pFluid->psatL_anc(TsatL);
			psatV = pFluid->psatV_anc(TsatV);
			// Saturation densities
			rhosatL = pFluid->density_Tp(TsatL, psatL, pFluid->rhosatL(TsatL));
			rhosatV = pFluid->density_Tp(TsatV, psatV, pFluid->rhosatV(TsatV));
		}
	}
	// Set internal variables
	_T = Q*TsatV+(1-Q)*TsatL;
	_rho = 1/(Q/rhosatV+(1-Q)/rhosatL);
	_p = Q*psatV+(1-Q)*psatL;
	_Q = Q;
}

// Updater if T,rho are inputs
void CoolPropStateClass::update_Trho(long iInput1, double Value1, long iInput2, double Value2)
{
	// Get them in the right order
	sort_pair(&iInput1,&Value1,&iInput2,&Value2,iT,iD);

	// Set internal variables
	_T = Value1;
	_rho = Value2;

	// If either SinglePhase or flag_SinglePhase is set to true, it will not make the call to the saturation routine
	// SinglePhase is set by the class routines, and flag_SinglePhase is a flag that can be set externally
	if (!SinglePhase || !flag_SinglePhase || !pFluid->phase_Trho(_T,_rho,&psatL,&psatV,&rhosatL,&rhosatV).compare("Two-Phase"))
	{
		// If it made it to the saturation routine and it is two-phase the saturation variables have been set
		TwoPhase = true;
		SinglePhase = false;

		// Get the quality and pressure
		_Q = (1/_rho-1/rhosatL)/(1/rhosatV-1/rhosatL);
		_p = _Q*psatV+(1-_Q)*psatL;
		
		check_saturated_quality(_Q);
	}
	else{
		TwoPhase = false;
		SinglePhase = true;
		SaturatedL = false;
		SaturatedV = false;
		_p = pFluid->pressure_Trho(_T,_rho);
	}
}

// Updater if T,p are inputs
void CoolPropStateClass::update_Tp(long iInput1, double Value1, long iInput2, double Value2)
{
	// Get them in the right order
	sort_pair(&iInput1,&Value1,&iInput2,&Value2,iT,iP);

	// Set internal variables
	_T = Value1;
	_p = Value2;

	// If either SinglePhase or flag_SinglePhase is set to true, it will not make the call to the saturation routine
	// SinglePhase is set by the class routines, and flag_SinglePhase is a flag that can be set externally
	if (!SinglePhase && !flag_SinglePhase && !pFluid->phase_Tp(_T,_p,&psatL,&psatV,&rhosatL,&rhosatV).compare("Two-Phase"))
	{
		// If it made it to the saturation routine and it is two-phase the saturation variables have been set
		TwoPhase = true;
		SinglePhase = false;

		// Get the quality and pressure
		_Q = (1/_rho-1/rhosatL)/(1/rhosatV-1/rhosatL);
		_p = _Q*psatV+(1-_Q)*psatL;
		
		check_saturated_quality(_Q);
		if (pFluid->pure()){
			TsatL = _T;
			TsatV = _T;
		}
		else{
			TsatL = _T;
			TsatV = _T;
		}
	}
	else{
		TwoPhase = false;
		SinglePhase = true;
		SaturatedL = false;
		SaturatedV = false;
		_rho = pFluid->density_Tp(_T,_p);
	}
}

// Updater if p,h are inputs
void CoolPropStateClass::update_ph(long iInput1, double Value1, long iInput2, double Value2, double T0, double rho0)
{
	// Get them in the right order
	sort_pair(&iInput1,&Value1,&iInput2,&Value2,iP,iH);

	// Solve for temperature and density with or without the guess values provided
	pFluid->temperature_ph(Value1, Value2,&_T,&_rho,&rhosatL,&rhosatV,&TsatL,&TsatV, T0, rho0);

	// Set internal variables
	_p = Value1;
	_h = Value2;

	// Set the phase flags
	if ( _T < pFluid->reduce.T && _rho < rhosatL && _rho > rhosatV)
	{
		TwoPhase = true;
		SinglePhase = false;
		_Q = (1/_rho-1/rhosatL)/(1/rhosatV-1/rhosatL);
		check_saturated_quality(_Q);

		psatL = _p;
		psatV = _p;
	}
	else
	{
		TwoPhase = false;
		SinglePhase = true;
		SaturatedL = false;
		SaturatedV = true;
	}
}

// Updater if p,s are inputs
void CoolPropStateClass::update_ps(long iInput1, double Value1, long iInput2, double Value2)
{
	// Get them in the right order
	sort_pair(&iInput1,&Value1,&iInput2,&Value2,iP,iS);

	// Solve for temperature and density
	pFluid->temperature_ps(Value1, Value2,&_T,&_rho,&rhosatL,&rhosatV,&TsatL,&TsatV);

	// Set internal variables
	_p = Value1;
	_s = Value2;

	// Set the phase flags
	if ( _T < pFluid->reduce.T && _rho < rhosatL && _rho > rhosatV)
	{
		TwoPhase = true;
		SinglePhase = false;
		_Q = (1/_rho-1/rhosatL)/(1/rhosatV-1/rhosatL);
		check_saturated_quality(_Q);

		psatL = _p;
		psatV = _p;
	}
	else
	{
		TwoPhase = false;
		SinglePhase = true;
		SaturatedL = false;
		SaturatedV = true;
	}
}

double CoolPropStateClass::hL(void){
	if (SaturationLUTStatus())
	{
		return pFluid->ApplySaturationLUT(pFluid->SatLUT.iHL,pFluid->SatLUT.iT,TsatL);
	}
	else
	{
		return pFluid->enthalpy_Trho(TsatL,rhosatL);
	}
}
double CoolPropStateClass::hV(void){
	if (SaturationLUTStatus())
	{
		return pFluid->ApplySaturationLUT(pFluid->SatLUT.iHV,pFluid->SatLUT.iT,TsatV);
	}
	else
	{
		return pFluid->enthalpy_Trho(TsatV,rhosatV);
	}
}
double CoolPropStateClass::sL(void){
	if (SaturationLUTStatus())
	{
		return pFluid->ApplySaturationLUT(pFluid->SatLUT.iSL,pFluid->SatLUT.iT,TsatL);
	}
	else
	{
		return pFluid->entropy_Trho(TsatL,rhosatL);
	}
}
double CoolPropStateClass::sV(void){
	if (SaturationLUTStatus())
	{
		return pFluid->ApplySaturationLUT(pFluid->SatLUT.iSV,pFluid->SatLUT.iT,TsatV);
	}
	else
	{
		return pFluid->entropy_Trho(TsatV,rhosatV);
	}
}

double CoolPropStateClass::h(void){
	if (TwoPhase){
		return _Q*hV()+(1-_Q)*hL();
	}
	else{
		if (fabs(_h)<1e90)
		{
			// Use the pre-calculated value
			return _h;
		}
		else
		{
			// Use the EOS
			return pFluid->enthalpy_Trho(_T,_rho);
		}
	}
}
double CoolPropStateClass::s(void){
	if (TwoPhase){
		return _Q*sV()+(1-_Q)*sL();
	}
	else{
		if (fabs(_s)<1e90)
		{
			// Use the pre-calculated value
			return _s;
		}
		else
		{
			return pFluid->entropy_Trho(_T,_rho);
		}
	}
}
double CoolPropStateClass::cp(void){
	return pFluid->specific_heat_p_Trho(_T,_rho);
}
double CoolPropStateClass::cv(void){
	return pFluid->specific_heat_v_Trho(_T,_rho);
}
double CoolPropStateClass::speed_sound(void){
	return pFluid->speed_sound_Trho(_T,_rho);
}
double CoolPropStateClass::drhodT_constp(void){
	return DerivTerms("drhodT|p",_T,_rho,pFluid,SinglePhase,TwoPhase);
}
double CoolPropStateClass::dpdrho_constT(void){
	return DerivTerms("dpdrho|T",_T,_rho,pFluid,SinglePhase,TwoPhase);
}
double CoolPropStateClass::drhodh_constp(void){
	return DerivTerms("drhodh|p",_T,_rho,pFluid,SinglePhase,TwoPhase);
}
double CoolPropStateClass::drhodp_consth(void){
	return DerivTerms("drhodp|h",_T,_rho,pFluid,SinglePhase,TwoPhase);
}
void CoolPropStateClass::dvdp_dhdp_sat(double T, double *dvdpL, double *dvdpV, double *dhdpL, double *dhdpV)
{
	CoolPropStateClass *sat = new CoolPropStateClass(pFluid);
	sat->update(iT,T,iQ,0.5);
	double vV =  1/sat->rhoV();
	double vL =  1/sat->rhoL();
	double hV =  sat->hV();
	double hL =  sat->hL();
	double rhoL = 1/vL;
	double rhoV = 1/vV;
	double TL = T;
	double TV = T;
	double tauL = pFluid->reduce.T/TL;
	double tauV = pFluid->reduce.T/TV;
	double deltaL = rhoL/pFluid->reduce.rho;
	double deltaV = rhoV/pFluid->reduce.rho;

	// For the liquid
	double d2phir_dDelta_dTauL = pFluid->d2phir_dDelta_dTau(tauL,deltaL);
	double dphir_dDeltaL = pFluid->dphir_dDelta(tauL,deltaL);
	double d2phir_dDelta2L = pFluid->d2phir_dDelta2(tauL,deltaL);
	double d2phi0_dTau2L = pFluid->d2phi0_dTau2(tauL,deltaL);
	double d2phir_dTau2L = pFluid->d2phir_dTau2(tauL,deltaL);
	// For the vapor
	double d2phir_dDelta_dTauV = pFluid->d2phir_dDelta_dTau(tauV,deltaV);
	double dphir_dDeltaV = pFluid->dphir_dDelta(tauV,deltaV);
	double d2phir_dDelta2V = pFluid->d2phir_dDelta2(tauV,deltaV);
	double d2phi0_dTau2V = pFluid->d2phi0_dTau2(tauV,deltaV);
	double d2phir_dTau2V = pFluid->d2phir_dTau2(tauV,deltaV);

	// Saturation derivatives at constant temperature
	double dhdrhoL_T = TL*pFluid->R()/rhoL*(tauL*deltaL*d2phir_dDelta_dTauL+deltaL*dphir_dDeltaL+deltaL*deltaL*d2phir_dDelta2L);
	double dhdrhoV_T = TV*pFluid->R()/rhoV*(tauV*deltaV*d2phir_dDelta_dTauV+deltaV*dphir_dDeltaV+deltaV*deltaV*d2phir_dDelta2V);
	double dpdrhoL_T = TL*pFluid->R()*(1+2*deltaL*dphir_dDeltaL+deltaL*deltaL*d2phir_dDelta2L);
	double dpdrhoV_T = TV*pFluid->R()*(1+2*deltaV*dphir_dDeltaV+deltaV*deltaV*d2phir_dDelta2V);
	
	// Saturation derivatives at constant density
	double dhdTL_rho = pFluid->R()*(-tauL*tauL*(d2phi0_dTau2L+d2phir_dTau2L)+1+deltaL*dphir_dDeltaL-deltaL*tauL*d2phir_dDelta_dTauL);
	double dhdTV_rho = pFluid->R()*(-tauV*tauV*(d2phi0_dTau2V+d2phir_dTau2V)+1+deltaV*dphir_dDeltaV-deltaV*tauV*d2phir_dDelta_dTauV);
	double dpdTL_rho = rhoL*pFluid->R()*(1+deltaL*dphir_dDeltaL-deltaL*tauL*d2phir_dDelta_dTauL);
	double dpdTV_rho = rhoV*pFluid->R()*(1+deltaV*dphir_dDeltaV-deltaV*tauV*d2phir_dDelta_dTauV);

	// Now get dh/dp at constant T for saturated liquid and vapor
	double dhdpL_T = dhdrhoL_T/dpdrhoL_T;
	double dhdpV_T = dhdrhoV_T/dpdrhoV_T;

	// Derivatives of enthalpy 
	double dhdTL_p = dhdTL_rho - dhdrhoL_T*dpdTL_rho/dpdrhoL_T;
	double dhdTV_p = dhdTV_rho - dhdrhoV_T*dpdTV_rho/dpdrhoV_T;

	// Derivatives of specific volume (just to make thing easier)
	double dvdrhoL = -1/(rhoL*rhoL);
	double dvdrhoV = -1/(rhoV*rhoV);

	// Derivatives of specific volume
	double dvdpL_T = 1/dpdrhoL_T*dvdrhoL;
	double dvdTL_p = -dpdTL_rho/dpdrhoL_T*dvdrhoL;
	double dvdpV_T = 1/dpdrhoV_T*dvdrhoV;
	double dvdTV_p = -dpdTV_rho/dpdrhoV_T*dvdrhoV;

	double dTsigmadp = T*(vV-vL)/(hV-hL);

	*dhdpL = dhdpL_T+dhdTL_p*dTsigmadp;
	*dhdpV = dhdpV_T+dhdTV_p*dTsigmadp;
	*dvdpL = dvdpL_T+dvdTL_p*dTsigmadp;
	*dvdpV = dvdpV_T+dvdTV_p*dTsigmadp;
}