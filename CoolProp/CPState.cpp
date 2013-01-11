
#include "CPExceptions.h"
#include "CPState.h"
#include "Solvers.h"

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
	this->clear_cache();

	// If flag_SinglePhase is true, it will always assume that it is not in the two-phase region
	// If flag_TwoPhase is true, it it always assume that you are in the two-phase region
	// Can be over-written by changing the flag to true
	flag_SinglePhase = false;
	flag_TwoPhase = false;

	SatL = NULL;
	SatV = NULL;
}

// Constructor with pointer to fluid
CoolPropStateClass::CoolPropStateClass(Fluid * pFluid){
	this->pFluid = pFluid;
	this->clear_cache();

	// If flag_SinglePhase is true, it will always assume that it is not in the two-phase region
	// If flag_TwoPhase is true, it it always assume that you are in the two-phase region
	// Can be over-written by changing the flag to true
	flag_SinglePhase = false;
	flag_TwoPhase = false;

	SatL = NULL;
	SatV = NULL;
}

CoolPropStateClass::~CoolPropStateClass()
{
	if (SatL != NULL)
	{
		delete SatL;
		SatL = NULL;
	}
	if (SatV != NULL)
	{
		delete SatV; 
		SatV = NULL;
	}
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
void CoolPropStateClass::clear_cache(void)
{
	cached_phi0 = false;
	cached_dphi0_dDelta = false;
	cached_dphi0_dTau = false;
	cached_d2phi0_dDelta2 = false;
	cached_d2phi0_dDelta_dTau = false;
	cached_d2phi0_dTau2 = false;
	cached_d3phi0_dDelta3 = false;
	cached_d3phi0_dDelta2_dTau = false;
	cached_d3phi0_dDelta_dTau2 = false;
	cached_d3phi0_dTau3 = false;

	cached_phir = false;
	cached_dphir_dDelta = false;
	cached_dphir_dTau = false;
	cached_d2phir_dDelta2 = false;
	cached_d2phir_dDelta_dTau = false;
	cached_d2phir_dTau2 = false;
	cached_d3phir_dDelta3 = false;
	cached_d3phir_dDelta2_dTau = false;
	cached_d3phir_dDelta_dTau2 = false;
	cached_d3phir_dTau3 = false;
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

	if (SatL == NULL)
	{
		SatL = new CoolPropStateClass(pFluid);
		SatV = new CoolPropStateClass(pFluid);
	}

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
	// Clear the cached derivative flags
	this->clear_cache();

	// Reduced parameters
	delta = this->_rho/pFluid->reduce.rho;
	tau = pFluid->reduce.T/this->_T;

	if (TwoPhase && !flag_TwoPhase)
	{
		add_saturation_states();
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
	if (!SinglePhase || !flag_SinglePhase || flag_TwoPhase || !pFluid->phase_Trho(_T,_rho,&psatL,&psatV,&rhosatL,&rhosatV).compare("Two-Phase"))
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

void CoolPropStateClass::add_saturation_states(void)
{
	SatL->flag_TwoPhase = true;
	SatL->update(iT,TsatL,iD,rhosatL);
	SatL->flag_TwoPhase = false;

	SatV->flag_TwoPhase = true;
	SatV->update(iT,TsatV,iD,rhosatV);
	SatV->flag_TwoPhase = false;
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
			// Use the EOS, using the cached value if possible
			return pFluid->R()*_T*(1.0+tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+delta*dphir_dDelta(tau,delta));
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
			// Use the EOS, using the cached value if possible
			return pFluid->R()*(tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))-phi0(tau,delta)-phir(tau,delta));
		}
	}
}
double CoolPropStateClass::cp(void){
	double c1 = pow(1.0+delta*dphir_dDelta(tau,delta)-delta*tau*d2phir_dDelta_dTau(tau,delta),2);
    double c2 = (1.0+2.0*delta*dphir_dDelta(tau,delta)+pow(delta,2)*d2phir_dDelta2(tau,delta));
	return pFluid->R()*(-pow(tau,2)*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta))+c1/c2);
}
double CoolPropStateClass::cv(void){
	return -pFluid->R()*pow(tau,2)*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta));
}
double CoolPropStateClass::speed_sound(void){
	double c1 = pow(tau,2)*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta));
    double c2 = (1.0+2.0*delta*dphir_dDelta(tau,delta)+pow(delta,2)*d2phir_dDelta2(tau,delta));
    return sqrt(-c2*this->_T*this->cp()*1000/c1);
}
double CoolPropStateClass::drhodT_constp(void){
	double dpdrho_T = pFluid->R()*_T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta));
	double dpdT_rho = pFluid->R()*_rho*(1+delta*dphir_dDelta(tau,delta)-delta*tau*d2phir_dDelta_dTau(tau,delta));
	return -dpdT_rho/dpdrho_T;
}
double CoolPropStateClass::dpdrho_constT(void){
	return pFluid->R()*_T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta));
}
double CoolPropStateClass::dpdT_constrho(void){
	return pFluid->R()*_rho*(1+delta*dphir_dDelta(tau,delta)-delta*tau*d2phir_dDelta_dTau(tau,delta));
}

double CoolPropStateClass::dhdrho_constT(void){
	return _T*pFluid->R()/_rho*(tau*delta*d2phir_dDelta_dTau(tau,delta)+delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta));
}
double CoolPropStateClass::dhdT_constrho(void){
	return pFluid->R()*(-tau*tau*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta))+1+delta*dphir_dDelta(tau,delta)-delta*tau*d2phir_dDelta_dTau(tau,delta));
}

double CoolPropStateClass::dhdp_constT(void){
	return dhdrho_constT()/dpdrho_constT();
}

double CoolPropStateClass::dhdT_constp(void){
	return dhdT_constrho() - dhdrho_constT()*dpdT_constrho()/dpdrho_constT();
}

double CoolPropStateClass::drhodh_constp(void){

	if (TwoPhase)
	{
		double vV = 1/rhosatV;
		double vL = 1/rhosatL;

		double dvdh_p = (vV-vL)/(SatV->h()-SatL->h());
		return -_rho*_rho*dvdh_p;
	}
	else
	{	
		return -dpdT_constrho()/dpdrho_constT()/cp();
	}
}
double CoolPropStateClass::drhodp_consth(void){
	if (TwoPhase)
	{
		return DerivTerms("drhodp|h",_T,_rho,pFluid,SinglePhase,TwoPhase);
	}
	else
	{
		return 1/(dpdrho_constT()-dpdT_constrho()*dhdrho_constT()/dhdT_constrho());
	}
}
void CoolPropStateClass::dvdp_dhdp_sat(double T, double *dvdpL, double *dvdpV, double *dhdpL, double *dhdpV, double *d2hdp2V)
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
	double d2phir_dDelta_dTauL = sat->SatL->d2phir_dDelta_dTau(tauL,deltaL);
	double dphir_dDeltaL = sat->SatL->dphir_dDelta(tauL,deltaL);
	double d2phir_dDelta2L = sat->SatL->d2phir_dDelta2(tauL,deltaL);
	double d2phi0_dTau2L = sat->SatL->d2phi0_dTau2(tauL,deltaL);
	double d2phir_dTau2L = sat->SatL->d2phir_dTau2(tauL,deltaL);
	// For the vapor
	double d2phir_dDelta_dTauV = sat->SatV->d2phir_dDelta_dTau(tauV,deltaV);
	double dphir_dDeltaV = sat->SatV->dphir_dDelta(tauV,deltaV);
	double d2phir_dDelta2V = sat->SatV->d2phir_dDelta2(tauV,deltaV);
	double d3phir_dDelta3V = sat->SatV->d3phir_dDelta3(tauV,deltaV);
	double d3phir_dDelta2_dTauV = sat->SatV->d3phir_dDelta2_dTau(tauV,deltaV);
	double d3phir_dDelta_dTau2V = sat->SatV->d3phir_dDelta_dTau2(tauV,deltaV);
	double d3phi0_dTau3V = sat->SatV->d3phi0_dTau3(tauV,deltaV);
	double d3phir_dTau3V = sat->SatV->d3phir_dTau3(tauV,deltaV);
	double d2phi0_dTau2V = sat->SatV->d2phi0_dTau2(tauV,deltaV);
	double d2phir_dTau2V = sat->SatV->d2phir_dTau2(tauV,deltaV);

	double drhodTV_p = DerivTerms("drhodT|p",TV,rhoV,(char*)pFluid->get_name().c_str());

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
	
	// Derivatives along the saturation boundary
	*dhdpL = dhdpL_T+dhdTL_p*dTsigmadp;
	*dhdpV = dhdpV_T+dhdTV_p*dTsigmadp;
	*dvdpL = dvdpL_T+dvdTL_p*dTsigmadp;
	*dvdpV = dvdpV_T+dvdTV_p*dTsigmadp;

	//double tau = Tc/T;
	double dtaudT = -pFluid->reduce.T/T/T; //-tau/T
	double d2Tsigma_dp2_T = T*((hV-hL)*(dvdpV_T-dvdpL_T)-(vV-vL)*(dhdpV_T-dhdpL_T))/pow(hV-hL,2);
	double d2Tsigma_dpdT_p = T*((hV-hL)*(dvdTV_p-dvdTL_p)-(vV-vL)*(dhdTV_p-dhdTL_p))/pow(hV-hL,2)+(vV-vL)/(hV-hL);
	//double d2Tsigmadp2 = d2Tsigma_dp2_T+d2Tsigma_dpdT_p*dTsigmadp;

	double d2pdrho2V_T = T*pFluid->R()*(2*deltaV*d2phir_dDelta2V+2*dphir_dDeltaV+2*deltaV*d2phir_dDelta2V+deltaV*deltaV*d3phir_dDelta3V)/pFluid->reduce.rho;
	double d2hdrho2V_T = TV*pFluid->R()/rhoV*(tauV*deltaV*d3phir_dDelta2_dTauV+tauV*d2phir_dDelta_dTauV+deltaV*d2phir_dDelta2V+dphir_dDeltaV+deltaV*deltaV*d3phir_dDelta3V+2*deltaV*d2phir_dDelta2V)/pFluid->reduce.rho - dhdrhoV_T/rhoV;
	// d3phi0_dDelta_dTau2V is zero by definition
	double d2hdrhodTV = pFluid->R()*(-tauV*tauV*d3phir_dDelta_dTau2V+deltaV*deltaV*d2phir_dDelta2V+deltaV*d2phir_dDelta2V+dphir_dDeltaV-deltaV*tauV*d3phir_dDelta2_dTauV-tauV*d2phir_dDelta_dTauV)/pFluid->reduce.rho;
	double d2pdrhodTV = pFluid->R()*((1+2*deltaV*dphir_dDeltaV+deltaV*deltaV*d2phir_dDelta2V)+T*(2*deltaV*d2phir_dDelta_dTauV+deltaV*deltaV*d3phir_dDelta2_dTauV)*dtaudT);

	double d2pdT2V_rho = -pFluid->R()*rhoV*deltaV*tauV*d3phir_dDelta_dTau2V*dtaudT;
	double d2hdT2V_rho = pFluid->R()*(-tauV*tauV*(d3phi0_dTau3V+d3phir_dTau3V)-2*tauV*(d2phi0_dTau2V+d2phir_dTau2V)-deltaV*tauV*d3phir_dDelta_dTau2V)*dtaudT;

	double d2hdp2V_T = (d2hdrho2V_T-dhdpV_T*d2pdrho2V_T)/pow(dpdrhoV_T,2);
	double d2hdTdpV = 1/dpdrhoV_T*(d2hdrhodTV-dhdpV_T*(drhodTV_p*d2pdrho2V_T+d2pdrhodTV)+d2hdrho2V_T*drhodTV_p);

	double ddT_dhdT = d2hdT2V_rho-1/pow(dpdrhoV_T,2)*(dpdrhoV_T*(dhdrhoV_T*d2pdT2V_rho+d2hdrhodTV*dpdTV_rho)-dhdrhoV_T*dpdTV_rho*d2pdrhodTV);
	double drho_dhdT = d2hdrhodTV-1/pow(dpdrhoV_T,2)*(dpdrhoV_T*(dhdrhoV_T*d2pdrhodTV+d2hdrho2V_T*dpdTV_rho)-dhdrhoV_T*dpdTV_rho*d2pdrho2V_T);

	// derivative of dhdp along saturation with respect to p with constant T
	double ddp_dhdpsigmaV = d2hdp2V_T+dhdTV_p*d2Tsigma_dp2_T+d2hdTdpV*dTsigmadp;

	double d2hdT2V_p = ddT_dhdT-drho_dhdT*dpdTV_rho/dpdrhoV_T;
	double ddT_dhdpsigmaV = d2hdTdpV+dhdTV_p*d2Tsigma_dpdT_p+d2hdT2V_p*dTsigmadp;

	// ------ GOOD TO THIS POINT -------------
	
	double dd = 1e-3;
	CoolPropStateClass CPS2 = CoolPropStateClass(pFluid);
	CPS2.flag_SinglePhase=true;
	CPS2.update(iT,_T,iD,rhoV+dd);

	double d2hdTdpVn = (CPS2.dhdT_constp()-this->dhdT_constp())/(dd)/dpdrho_constT();
	double d2pdT2Vn = (pFluid->pressure_Trho(TV-dd,rhoV)-2*pFluid->pressure_Trho(TV,rhoV)+pFluid->pressure_Trho(TV+dd,rhoV))/dd/dd;
	double d2hdT2Vn = (pFluid->enthalpy_Trho(TV-dd,rhoV)-2*pFluid->enthalpy_Trho(TV,rhoV)+pFluid->enthalpy_Trho(TV+dd,rhoV))/dd/dd;
	double d2pdrhodTVn = (pFluid->pressure_Trho(TV+dd,rhoV+dd)-pFluid->pressure_Trho(TV+dd,rhoV-dd)-pFluid->pressure_Trho(TV-dd,rhoV+dd)+pFluid->pressure_Trho(TV-dd,rhoV-dd))/dd/dd/4;

	// Find T(p,rho)

	/// A stub class to do the density(T,p) calculations for near the critical point using Brent solver
	class TprhoResids : public FuncWrapper1D
	{
	private:
		double p,rho;
		Fluid *pFluid;
	public:
		TprhoResids(Fluid *pFluid, double rho, double p){this->pFluid = pFluid; this->p = p; this->rho = rho;};
		~TprhoResids(){};
		
		double call(double T)
		{
			return this->p - pFluid->pressure_Trho(T,rho);
		}
	};

	std::string errstr;
	TprhoResids TPR1 = TprhoResids(pFluid,rhoV-dd,sat->pV());
	double T1 = Secant(&TPR1,TV,1e-3,1e-10,100,&errstr);

	TprhoResids TPR2 = TprhoResids(pFluid,rhoV,sat->pV());
	double T2 = Secant(&TPR2,TV,1e-3,1e-10,100,&errstr);

	TprhoResids TPR3 = TprhoResids(pFluid,rhoV+dd,sat->pV());
	double T3 = Secant(&TPR3,TV,1e-3,1e-10,100,&errstr);

	double p1 = pFluid->pressure_Trho(T1,rhoV-dd);
	double p2 = pFluid->pressure_Trho(T2,rhoV);
	double p3 = pFluid->pressure_Trho(T3,rhoV+dd);

	double h1 = pFluid->enthalpy_Trho(T1,rhoV-dd);
	double h2 = pFluid->enthalpy_Trho(T2,rhoV);
	double h3 = pFluid->enthalpy_Trho(T3,rhoV+dd);

	double dhdTV_rhon = (h2-h1)/(T2-T1);
	double dhdTV_rhon2 = (h3-h2)/(T3-T2);
	double d2hdT2V_pn = (dhdTV_rhon2-dhdTV_rhon)/(T3-T1)*2;

	std::cout<<format("%g,%g,%g\n",T1,T2,T3);
	std::cout<<format("%g,%g,%g\n",h1,h2,h3);
	
	CoolPropStateClass CPS3 = CoolPropStateClass(pFluid);
	CPS3.flag_SinglePhase=true;
	CPS3.update(iT,_T+dd,iD,rhoV);

	double dT_dhdTn = (CPS3.dhdT_constp()-this->dhdT_constp())/(dd);

	

	
	
	*d2hdp2V = ddp_dhdpsigmaV+ddT_dhdpsigmaV*dTsigmadp;

	std::cout<<format("d2hsigmadp2V = %g\n",*d2hdp2V);
}

// All the derivatives of the ideal-gas and residual Helmholtz energy
double CoolPropStateClass::phi0(double tau, double delta){
	if (cached_phi0) 
	{
		return cachedval_phi0; 
	}
	else 
	{
		cached_phi0 = true;
		cachedval_phi0 = pFluid->phi0(tau,delta);
		return pFluid->phi0(tau,delta);
	}
};
double CoolPropStateClass::dphi0_dDelta(double tau, double delta){
	if (cached_dphi0_dDelta)
	{
		return cachedval_dphi0_dDelta; 
	}
	else 
	{
		cached_dphi0_dDelta = true;
		cachedval_dphi0_dDelta = pFluid->dphi0_dDelta(tau,delta);
		return cachedval_dphi0_dDelta;
	}
};
double CoolPropStateClass::dphi0_dTau(double tau, double delta){
	if (cached_dphi0_dTau)
	{
		return cachedval_dphi0_dTau; 
	}
	else 
	{
		cached_dphi0_dTau = true;
		cachedval_dphi0_dTau = pFluid->dphi0_dTau(tau,delta);
		return cachedval_dphi0_dTau;
	}
};
double CoolPropStateClass::d2phi0_dDelta2(double tau, double delta){
	if (cached_d2phi0_dDelta2) 
	{
		return cachedval_d2phi0_dDelta2; 
	}
	else 
	{
		cached_d2phi0_dDelta2 = true;
		cachedval_d2phi0_dDelta2 = pFluid->d2phi0_dDelta2(tau,delta);
		return cachedval_d2phi0_dDelta2;
	}
};
double CoolPropStateClass::d2phi0_dDelta_dTau(double tau, double delta){
	if (cached_d2phi0_dDelta_dTau) 
	{
		return cachedval_d2phi0_dDelta_dTau; 
	}
	else 
	{
		cached_d2phi0_dDelta_dTau = true;
		cachedval_d2phi0_dDelta_dTau = pFluid->d2phi0_dDelta_dTau(tau,delta);
		return cachedval_d2phi0_dDelta_dTau;
	}
};
double CoolPropStateClass::d2phi0_dTau2(double tau, double delta){
	if (cached_d2phi0_dTau2) 
	{
		return cachedval_d2phi0_dTau2; 
	}
	else 
	{
		cached_d2phi0_dTau2 = true;
		cachedval_d2phi0_dTau2 = pFluid->d2phi0_dTau2(tau,delta);
		return cachedval_d2phi0_dTau2;
	}
};

double CoolPropStateClass::d3phi0_dTau3(double tau, double delta){
	if (cached_d3phi0_dTau3) 
	{
		return cachedval_d3phi0_dTau3; 
	}
	else 
	{
		cached_d3phi0_dTau3 = true;
		cachedval_d3phi0_dTau3 = pFluid->d3phi0_dTau3(tau,delta);
		return cachedval_d3phi0_dTau3;
	}
};

double CoolPropStateClass::d3phi0_dDelta_dTau2(double tau, double delta){
	return 0;
};

double CoolPropStateClass::d3phi0_dDelta2_dTau(double tau, double delta){
	return 0;
};

double CoolPropStateClass::d3phi0_dDelta3(double tau, double delta){
	if (cached_d3phi0_dDelta3) 
	{
		return cachedval_d3phi0_dDelta3; 
	}
	else 
	{
		cached_d3phi0_dDelta3 = true;
		cachedval_d3phi0_dDelta3 = pFluid->d3phi0_dDelta3(tau,delta);
		return cachedval_d3phi0_dDelta3;
	}
};


double CoolPropStateClass::phir(double tau, double delta){
	if (cached_phir) 
	{
		return cachedval_phir; 
	}
	else 
	{
		cached_phir = true;
		cachedval_phir = pFluid->phir(tau,delta);
		return pFluid->phir(tau,delta);
	}
};
double CoolPropStateClass::dphir_dDelta(double tau, double delta){
	if (cached_dphir_dDelta)
	{
		return cachedval_dphir_dDelta; 
	}
	else 
	{
		cached_dphir_dDelta = true;
		cachedval_dphir_dDelta = pFluid->dphir_dDelta(tau,delta);
		return cachedval_dphir_dDelta;
	}
};
double CoolPropStateClass::dphir_dTau(double tau, double delta){
	if (cached_dphir_dTau)
	{
		return cachedval_dphir_dTau; 
	}
	else 
	{
		cached_dphir_dTau = true;
		cachedval_dphir_dTau = pFluid->dphir_dTau(tau,delta);
		return cachedval_dphir_dTau;
	}
};
double CoolPropStateClass::d2phir_dDelta2(double tau, double delta){
	if (cached_d2phir_dDelta2) 
	{
		return cachedval_d2phir_dDelta2; 
	}
	else 
	{
		cached_d2phir_dDelta2 = true;
		cachedval_d2phir_dDelta2 = pFluid->d2phir_dDelta2(tau,delta);
		return cachedval_d2phir_dDelta2;
	}
};
double CoolPropStateClass::d2phir_dDelta_dTau(double tau, double delta){
	if (cached_d2phir_dDelta_dTau) 
	{
		return cachedval_d2phir_dDelta_dTau; 
	}
	else 
	{
		cached_d2phir_dDelta_dTau = true;
		cachedval_d2phir_dDelta_dTau = pFluid->d2phir_dDelta_dTau(tau,delta);
		return cachedval_d2phir_dDelta_dTau;
	}
};
double CoolPropStateClass::d2phir_dTau2(double tau, double delta){
	if (cached_d2phir_dTau2) 
	{
		return cachedval_d2phir_dTau2; 
	}
	else 
	{
		cached_d2phir_dTau2 = true;
		cachedval_d2phir_dTau2 = pFluid->d2phir_dTau2(tau,delta);
		return cachedval_d2phir_dTau2;
	}
};

double CoolPropStateClass::d3phir_dTau3(double tau, double delta){
	if (cached_d3phir_dTau3) 
	{
		return cachedval_d3phir_dTau3; 
	}
	else 
	{
		cached_d3phir_dTau3 = true;
		cachedval_d3phir_dTau3 = pFluid->d3phir_dTau3(tau,delta);
		return cachedval_d3phir_dTau3;
	}
};

double CoolPropStateClass::d3phir_dDelta_dTau2(double tau, double delta){
	if (cached_d3phir_dDelta_dTau2) 
	{
		return cachedval_d3phir_dDelta_dTau2; 
	}
	else 
	{
		cached_d3phir_dDelta_dTau2 = true;
		cachedval_d3phir_dDelta_dTau2 = pFluid->d3phir_dDelta_dTau2(tau,delta);
		return cachedval_d3phir_dDelta_dTau2;
	}
};

double CoolPropStateClass::d3phir_dDelta2_dTau(double tau, double delta){
	if (cached_d3phir_dDelta2_dTau) 
	{
		return cachedval_d3phir_dDelta2_dTau; 
	}
	else 
	{
		cached_d3phir_dDelta2_dTau = true;
		cachedval_d3phir_dDelta2_dTau = pFluid->d3phir_dDelta2_dTau(tau,delta);
		return cachedval_d3phir_dDelta2_dTau;
	}
};

double CoolPropStateClass::d3phir_dDelta3(double tau, double delta){
	if (cached_d3phir_dDelta3) 
	{
		return cachedval_d3phir_dDelta3; 
	}
	else 
	{
		cached_d3phir_dDelta3 = true;
		cachedval_d3phir_dDelta3 = pFluid->d3phir_dDelta3(tau,delta);
		return cachedval_d3phir_dDelta3;
	}
};