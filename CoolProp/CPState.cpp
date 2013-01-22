
#include "CPExceptions.h"
#include "Solvers.h"
#include "CPState.h"

#ifndef DBL_EPSILON
	#include <limits>
	#define DBL_EPSILON std::numeric_limits<double>::epsilon()
#endif

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
	
	if (pFluid->enabled_TTSE_LUT)
	{
		// Build the tables
		pFluid->build_TTSE_LUT();
		// Update using the LUT
		update_TTSE_LUT(iInput1,Value1,iInput2,Value2);
	}
	else
	{
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

		if (TwoPhase && !flag_TwoPhase)
		{
			// Update temperature and density for SatL and SatV
			add_saturation_states();
		}
	}
	
	

	// Reduced parameters
	delta = this->_rho/pFluid->reduce.rho;
	tau = pFluid->reduce.T/this->_T;	
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
				psatV = Value1;
				psatL = Value1;
			}
			else{
				pFluid->saturation_p(Value1,SaturationLUTStatus(),&TsatL,&TsatV,&rhosatL,&rhosatV);
				TsatV = TsatL;
				psatV = Value1;
				psatL = Value1;
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
		// for the given temperature
		if (pFluid->pure()){
			pFluid->saturation_T(Value1,SaturationLUTStatus(),&psatL,&psatV,&rhosatL,&rhosatV);
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
	if ((!SinglePhase && !flag_SinglePhase) || flag_TwoPhase || !pFluid->phase_Trho(_T,_rho,&psatL,&psatV,&rhosatL,&rhosatV).compare("Two-Phase"))
	{
		if (!flag_SinglePhase || flag_TwoPhase)
		{
			pFluid->phase_Trho(_T,_rho,&psatL,&psatV,&rhosatL,&rhosatV);
			TsatV = _T;
			TsatL = _T;
			
		}
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

// Updater if you are using TTSE LUT
void CoolPropStateClass::update_TTSE_LUT(long iInput1, double Value1, long iInput2, double Value2)
{
	
	// If the inputs are P,Q or T,Q , it is guaranteed to be two-phase
	if (match_pair(iInput1,iInput2,iP,iQ))
	{
		// Set phase flags
		SinglePhase = false;
		TwoPhase = true;
		
		// Sort in the right order (P,Q)
		sort_pair(&iInput1,&Value1,&iInput2,&Value2,iP,iQ);

		double p = Value1;
		double Q = Value2;

		// Check whether saturated liquid or vapor
		check_saturated_quality(Q);

		rhosatL = pFluid->TTSESatL.evaluate(iD,p);
		rhosatV = pFluid->TTSESatV.evaluate(iD,p);
		TsatL = pFluid->TTSESatL.evaluate(iT,p);
		TsatV = pFluid->TTSESatV.evaluate(iT,p);
		psatL = p;
		psatV = p;

		// Set internal variables
		_T = Q*TsatV+(1-Q)*TsatL;
		_rho = 1/(Q/rhosatV+(1-Q)/rhosatL);
		_p = Q*psatV+(1-Q)*psatL;
		_Q = Q;
	}
	else if (match_pair(iInput1,iInput2,iP,iH))
	{
		// Sort in the right order (P,H)
		sort_pair(&iInput1,&Value1,&iInput2,&Value2,iP,iH);

		double p = Value1; 
		double h = Value2;

		// If enthalpy is outside the saturation region, it is single-phase
		if (p > pFluid->reduce.p || p < pFluid->params.ptriple ||  h < pFluid->TTSESatL.evaluate(iH,p)  || h > pFluid->TTSESatV.evaluate(iH,p))
		{
			TwoPhase = false;
			SinglePhase = true;

			_rho = pFluid->TTSESinglePhase.evaluate(iD,p,h);
			_T = pFluid->TTSESinglePhase.evaluate(iT,p,h);
			_p = p;
			_h = h;
		}
		else
		{
			TwoPhase = true;
			SinglePhase = false;

			double hsatL = pFluid->TTSESatL.evaluate(iH,p);
			double hsatV = pFluid->TTSESatV.evaluate(iH,p);
			rhosatL = pFluid->TTSESatL.evaluate(iD,p);
			rhosatV = pFluid->TTSESatV.evaluate(iD,p);
			TsatL = pFluid->TTSESatL.evaluate(iT,p);
			TsatV = pFluid->TTSESatV.evaluate(iT,p);
			psatL = p;
			psatV = p;
			
			_Q = (h-hsatL)/(hsatV-hsatL);
			_rho = 1/(_Q/rhosatV+(1-_Q)/rhosatL);
			_T = _Q*TsatV+(1-_Q)*TsatL;
			_p = _Q*psatV+(1-_Q)*psatL;
			_h = h;

			check_saturated_quality(_Q);
		}
	}
	else
	{
		printf("Sorry your inputs[%d,%d] don't work for now with TTSE\n",iInput1,iInput2);
		throw ValueError(format("Sorry your inputs don't work for now with TTSE"));
	}
}

void CoolPropStateClass::add_saturation_states(void)
{
	// While SatL and SatV are technically two-phase, we consider 
	// them to be single-phase to speed up the calcs and avoid saturation calls
	SatL->flag_TwoPhase = true;
	SatL->update(iT,TsatL,iD,rhosatL);
	SatL->flag_TwoPhase = false;
	SatL->SinglePhase = true;
	SatL->TwoPhase = false;

	SatV->flag_TwoPhase = true;
	SatV->update(iT,TsatV,iD,rhosatV);
	SatV->flag_TwoPhase = false;
	SatV->SinglePhase = true;
	SatV->TwoPhase = false;
}

/// Enable the TTSE LUT
void CoolPropStateClass::enable_TTSE_LUT(void){
	pFluid->enabled_TTSE_LUT = true;
}
/// Check if TTSE is enabled
bool CoolPropStateClass::isenabled_TTSE_LUT(void){
	return pFluid->enabled_TTSE_LUT;
};
/// Disable the TTSE
void CoolPropStateClass::disable_TTSE_LUT(void){
	pFluid->enabled_TTSE_LUT = false;
};

double CoolPropStateClass::hL(void){
	if (pFluid->enabled_TTSE_LUT)
	{
		pFluid->build_TTSE_LUT();
		return pFluid->TTSESatL.evaluate(iH,psatL);
	}
	else
	{
		return SatL->h();
	}
}
double CoolPropStateClass::hV(void){
	if (pFluid->enabled_TTSE_LUT)
	{
		pFluid->build_TTSE_LUT();
		return pFluid->TTSESatV.evaluate(iH,psatV);
	}
	else
	{
		return SatV->h();
	}
}
double CoolPropStateClass::sL(void){
	if (pFluid->enabled_TTSE_LUT)
	{
		pFluid->build_TTSE_LUT();
		return pFluid->TTSESatL.evaluate(iH,psatL);
	}
	else
	{
		return SatL->s();
	}
}
double CoolPropStateClass::sV(void){
	if (pFluid->enabled_TTSE_LUT)
	{
		pFluid->build_TTSE_LUT();
		return pFluid->TTSESatV.evaluate(iS,psatV);
	}
	else
	{
		return SatV->s();
	}
}

double CoolPropStateClass::h(void){
	if (TwoPhase){
		// This will use the TTSE LUT if enable_TTSE_LUT() has been called
		return _Q*hV()+(1-_Q)*hL();
	}
	else{
		if (fabs(_h)<1e90){
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
		// This will use the TTSE LUT if enable_TTSE_LUT() has been called
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
	if (pFluid->enabled_TTSE_LUT)
	{
		if (TwoPhase && _Q>0 && _Q < 1)
		{
			return -1;
		}
		else
		{
			pFluid->build_TTSE_LUT();
			// cp is also given by (dh/dT)|p, or 1/((dT/dh)|p) which is tabulated
			return 1/pFluid->TTSESinglePhase.evaluate_first_derivative(iT,iH,iP,_p,_h);
		}
	}
	else
	{
		double c1 = pow(1.0+delta*dphir_dDelta(tau,delta)-delta*tau*d2phir_dDelta_dTau(tau,delta),2);
		double c2 = (1.0+2.0*delta*dphir_dDelta(tau,delta)+pow(delta,2)*d2phir_dDelta2(tau,delta));
		return pFluid->R()*(-pow(tau,2)*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta))+c1/c2);
	}
}
double CoolPropStateClass::B_TTSE(double p, double h){
	// Slightly modified, doesn't use specific volume at all, also sign switched
	double drhodh_p = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iH,iP,p,h);
	double drhodp_h = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iP,iH,p,h);
	double dsdh_p = pFluid->TTSESinglePhase.evaluate_first_derivative(iS,iH,iP,p,h);
	double dsdp_h = pFluid->TTSESinglePhase.evaluate_first_derivative(iS,iP,iH,p,h);
	double rho = pFluid->TTSESinglePhase.evaluate(iD,p,h);
	return (drhodp_h*dsdh_p-drhodh_p*dsdp_h);
}
double CoolPropStateClass::B_over_D_TTSE(double p, double h)
{
	double drhodh_p = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iH,iP,p,h);
	double drhodp_h = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iP,iH,p,h);
	double   dsdh_p = pFluid->TTSESinglePhase.evaluate_first_derivative(iS,iH,iP,p,h);
	double   dsdp_h = pFluid->TTSESinglePhase.evaluate_first_derivative(iS,iP,iH,p,h);
	double   dTdh_p = pFluid->TTSESinglePhase.evaluate_first_derivative(iT,iH,iP,p,h);
	double   dTdp_h = pFluid->TTSESinglePhase.evaluate_first_derivative(iT,iP,iH,p,h);
	return (drhodh_p*dsdp_h-drhodp_h*dsdh_p)/(dTdh_p*drhodp_h-dTdp_h*drhodh_p);
}
double CoolPropStateClass::cv(void){
	if (pFluid->enabled_TTSE_LUT)
	{
		if (TwoPhase && _Q>0 && _Q < 1)
		{
			return -1;
		}
		else
		{
			pFluid->build_TTSE_LUT();
			// cv is also given by -B*T/D which is tabulated (indirectly)
			return -B_over_D_TTSE(_p,_h)*pFluid->TTSESinglePhase.evaluate(iT,_p,_h);
		}
	}
	else
	{
		return -pFluid->R()*pow(tau,2)*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta));
	}
}
double CoolPropStateClass::speed_sound(void){
	if (pFluid->enabled_TTSE_LUT)
	{
		if (TwoPhase && _Q>0 && _Q < 1)
		{
			// Not defined for two-phase
			return -1;
		}
		else
		{
			pFluid->build_TTSE_LUT();
			// speed of sound given by sqrt(v^2/B*dsdh|p), or sqrt(dpdrho|s)
			double dsdh__p = pFluid->TTSESinglePhase.evaluate_first_derivative(iS,iH,iP,_p,_h);
			double dsdp__h = pFluid->TTSESinglePhase.evaluate_first_derivative(iS,iP,iH,_p,_h);
			double drhodh__p = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iH,iP,_p,_h);
			double drhodp__h = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iP,iH,_p,_h);
			
			/// Factor of 1000 is because units within radical need to be all in base units to result in m^2/s^3
			return 1/sqrt((drhodp__h-drhodh__p*dsdp__h/dsdh__p)/1000);
		}
	}
	else
	{
		double c1 = pow(tau,2)*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta));
		double c2 = (1.0+2.0*delta*dphir_dDelta(tau,delta)+pow(delta,2)*d2phir_dDelta2(tau,delta));
		return sqrt(-c2*this->_T*this->cp()*1000/c1);
	}
}

double CoolPropStateClass::isothermal_compressibility(void){

	if (pFluid->enabled_TTSE_LUT)
	{
		if (TwoPhase && _Q>0 && _Q < 1)
		{
			// Not defined for two-phase
			return -1;
		}
		else
		{
			pFluid->build_TTSE_LUT();
			// isothermal compressibility given by kappa = -1/v*dvdp|T = 1/rho*drhodp|T
			double dTdh__p = pFluid->TTSESinglePhase.evaluate_first_derivative(iT,iH,iP,_p,_h);
			double dTdp__h = pFluid->TTSESinglePhase.evaluate_first_derivative(iT,iP,iH,_p,_h);
			double drhodh__p = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iH,iP,_p,_h);
			double drhodp__h = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iP,iH,_p,_h);

			double rho = pFluid->TTSESinglePhase.evaluate(iD,_p,_h);
			/// 1000 is needed to convert from kJ & kPa to J & Pa
			return 1/rho*(drhodp__h-drhodh__p*dTdp__h/dTdh__p)/1000;
		}
	}
	else
	{
		return 1/(_rho*dpdrho_constT());
	}
}

double CoolPropStateClass::isobaric_expansion_coefficient(void){

	if (pFluid->enabled_TTSE_LUT)
	{
		if (TwoPhase && _Q>0 && _Q < 1)
		{
			// Not defined for two-phase
			return -1;
		}
		else
		{
			pFluid->build_TTSE_LUT();
			// isobaric expansion coefficient given by kappa = 1/v*dvdT|p = -1/rho*drhodT|p
			double dTdh__p = pFluid->TTSESinglePhase.evaluate_first_derivative(iT,iH,iP,_p,_h);
			double drhodh__p = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iH,iP,_p,_h);
			double rho = pFluid->TTSESinglePhase.evaluate(iD,_p,_h);
			return -1/rho*drhodh__p/dTdh__p;
		}
	}
	else
	{
		return -1/(_rho*_rho)*drhodT_constp();
	}
}

double CoolPropStateClass::drhodh_constp(void){

	if (pFluid->enabled_TTSE_LUT)
	{
		if (TwoPhase && _Q>0 && _Q < 1)
		{
			// equals -rho^2*dvdh_p where dvdh_p = 1/T*dTdp|sat
			return -_rho*_rho/_T*pFluid->TTSESatL.evaluate_sat_derivative(iT,_p);
		}
		else
		{
			pFluid->build_TTSE_LUT();
			return pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iH,iP,_p,_h);
		}
	}
	else
	{
		if (TwoPhase)
		{
			// equals -rho^2*dvdh_p where dvdh_p = 1/T*dTdp|sat
			return -_rho*_rho/_T*dTdp_along_sat();
		}
		else
		{
			return 1/(dhdrho_constT()-dhdT_constrho()*dpdrho_constT()/dpdT_constrho());
		}
	}
}

double CoolPropStateClass::drhodp_consth(void){

	if (pFluid->enabled_TTSE_LUT)
	{
		if (TwoPhase && _Q>0 && _Q < 1)
		{
			double hL = pFluid->TTSESatL.evaluate(iH,_p);
			double hV = pFluid->TTSESatV.evaluate(iH,_p);
			double rhoL = pFluid->TTSESatL.evaluate(iD,_p);
			double rhoV = pFluid->TTSESatV.evaluate(iD,_p);
			double dhdpL = pFluid->TTSESatL.evaluate_sat_derivative(iH,_p);
			double dhdpV = pFluid->TTSESatV.evaluate_sat_derivative(iH,_p);
			double dvdpL = -pFluid->TTSESatL.evaluate_sat_derivative(iD,_p)/rhoL/rhoL;
			double dvdpV = -pFluid->TTSESatV.evaluate_sat_derivative(iD,_p)/rhoV/rhoV;
			
			double dxdp_h = (_Q*dhdpV+(1-_Q)*(dhdpV-dhdpL))/(hL-hV);
			double dvdp_h = dvdpL+dxdp_h*(1/rhoV-1/rhoL)+_Q*(dvdpV-dvdpL);

			return -_rho*_rho*dvdp_h;
		}
		else
		{
			pFluid->build_TTSE_LUT();
			return pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iP,iH,_p,_h);
		}
	}
	else
	{
		if (TwoPhase)
		{
			double dhdpL = dhdp_along_sat_liquid();
			double dhdpV = dhdp_along_sat_vapor();
			double dvdpL = -drhodp_along_sat_liquid()/rhosatL/rhosatL;
			double dvdpV = -drhodp_along_sat_vapor()/rhosatV/rhosatV;
			
			double dxdp_h = (dhdpL+_Q*(dhdpV-dhdpL))/(hL()-hV());

			double dhdpV_2 = (Props("H",'P',p()+1e-5,'Q',1,"Propane")-Props("H",'P',p(),'Q',1,"Propane"))/1e-5;
			double dhdpL_2 = (Props("H",'P',p()+1e-5,'Q',0,"Propane")-Props("H",'P',p(),'Q',0,"Propane"))/1e-5;
			double dxdp_h2 = (Props("Q",'P',p()+1e-3,'H',h(),"Propane")-Props("Q",'P',p(),'H',h(),"Propane"))/1e-3;
			double dvdp_h2 = (1/Props("D",'P',p()+1e-3,'H',h(),"Propane")-1/Props("D",'P',p(),'H',h(),"Propane"))/1e-3;
			double drhodp_h2 = (Props("D",'P',p()+1e-3,'H',h(),"Propane")-Props("D",'P',p(),'H',h(),"Propane"))/1e-3;
			double gretresgr = Props("Q",'P',p(),'H',_h,"Propane");
			
			double dvdp_h = dvdpL+dxdp_h*(1/rhosatV-1/rhosatL)+_Q*(dvdpV-dvdpL);

			return -_rho*_rho*dvdp_h;
		}
		else
		{
			return 1/(dpdrho_constT()-dpdT_constrho()*dhdrho_constT()/dhdT_constrho());
		}
	}
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
double CoolPropStateClass::d2pdrho2_constT(void){
	return _T*pFluid->R()*(2*delta*d2phir_dDelta2(tau,delta)+2*dphir_dDelta(tau,delta)+2*delta*d2phir_dDelta2(tau,delta)+delta*delta*d3phir_dDelta3(tau,delta))/pFluid->reduce.rho;
}
double CoolPropStateClass::d2pdrhodT(void){
	return pFluid->R()*((1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta))+_T*(2*delta*d2phir_dDelta_dTau(tau,delta)+delta*delta*d3phir_dDelta2_dTau(tau,delta))*(-tau/_T));
}
double CoolPropStateClass::d2pdT2_constrho(void){
	return pFluid->R()*_rho*delta*tau*tau/_T*d3phir_dDelta_dTau2(tau,delta);
}

// DERIVATIVES OF ENTROPY FROM EOS
double CoolPropStateClass::dhdrho_constT(void){
	return _T*pFluid->R()/_rho*(tau*delta*d2phir_dDelta_dTau(tau,delta)+delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta));
}
double CoolPropStateClass::dhdT_constrho(void){
	return pFluid->R()*(-tau*tau*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta))+1+delta*dphir_dDelta(tau,delta)-delta*tau*d2phir_dDelta_dTau(tau,delta));
}
double CoolPropStateClass::d2hdrho2_constT(void){
	return _T*pFluid->R()/_rho*(tau*delta*d3phir_dDelta2_dTau(tau,delta)+tau*d2phir_dDelta_dTau(tau,delta)+delta*d2phir_dDelta2(tau,delta)+dphir_dDelta(tau,delta)+delta*delta*d3phir_dDelta3(tau,delta)+2*delta*d2phir_dDelta2(tau,delta))/pFluid->reduce.rho - dhdrho_constT()/_rho;
}
double CoolPropStateClass::d2hdrhodT(void){
	// d3phi0_dDelta_dTau2V is zero by definition
	return pFluid->R()*(-tau*tau*d3phir_dDelta_dTau2(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta)+delta*d2phir_dDelta2(tau,delta)+dphir_dDelta(tau,delta)-delta*tau*d3phir_dDelta2_dTau(tau,delta)-tau*d2phir_dDelta_dTau(tau,delta))/pFluid->reduce.rho;
}
double CoolPropStateClass::d2hdT2_constrho(void){
	return pFluid->R()*(-tau*tau*(d3phi0_dTau3(tau,delta)+d3phir_dTau3(tau,delta))-2*tau*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta))-delta*tau*d3phir_dDelta_dTau2(tau,delta))*(-tau/_T);
}

// DERIVATIVES OF ENTROPY FROM EOS
double CoolPropStateClass::dsdrho_constT(void){
	return -pFluid->R()/_rho*(1+delta*dphir_dDelta(tau,delta)-delta*tau*d2phir_dDelta_dTau(tau,delta));
}
double CoolPropStateClass::dsdT_constrho(void){
	return -pFluid->R()*tau*tau/_T*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta));
}
double CoolPropStateClass::d2sdT2_constrho(void){
	return -pFluid->R()/_T*(tau*tau*(d3phi0_dTau3(tau,delta)+d3phir_dTau3(tau,delta))+2*tau*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta)))*(-tau/_T)+pFluid->R()*tau*tau/_T/_T*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta));
}
double CoolPropStateClass::d2sdrho2_constT(void){
	return -pFluid->R()/_rho*(delta*d2phir_dDelta2(tau,delta)+dphir_dDelta(tau,delta)-tau*delta*d3phir_dDelta2_dTau(tau,delta)-tau*d2phir_dDelta_dTau(tau,delta))/pFluid->reduce.rho+pFluid->R()/_rho/_rho*(1+delta*dphir_dDelta(tau,delta)-delta*tau*d2phir_dDelta_dTau(tau,delta));
}
double CoolPropStateClass::d2sdrhodT(void){
	// d2phi0_dDelta_dTau2(tau,delta) is zero by definition
	return -pFluid->R()*tau*tau/_T*d3phir_dDelta_dTau2(tau,delta)/pFluid->reduce.rho;
}


double CoolPropStateClass::dvdp_constT(void){
	return -1/(_rho*_rho)/dpdrho_constT();
}
double CoolPropStateClass::dvdT_constp(void){
	return -1/(_rho*_rho)*drhodT_constp();
}

double CoolPropStateClass::dpdT_consth(void){
	return dpdT_constrho() - dpdrho_constT()*dhdT_constrho()/dhdrho_constT();
}
double CoolPropStateClass::dpdrho_consth(void){
	return dpdrho_constT() - dpdT_constrho()*dhdrho_constT()/dhdT_constrho();
}
// Enthalpy
double CoolPropStateClass::dhdp_constT(void){
	return dhdrho_constT()/dpdrho_constT();
}
double CoolPropStateClass::dhdT_constp(void){
	return dhdT_constrho() - dhdrho_constT()*dpdT_constrho()/dpdrho_constT();
}
double CoolPropStateClass::dhdrho_constp(void){
	return dhdrho_constT() - dhdT_constrho()*dpdrho_constT()/dpdT_constrho();
}
double CoolPropStateClass::d2hdT2_constp(void)
{
	double ddT_dhdT = d2hdT2_constrho()-1/pow(dpdrho_constT(),2)*(dpdrho_constT()*(dhdrho_constT()*d2pdT2_constrho()+d2hdrhodT()*dpdT_constrho())-dhdrho_constT()*dpdT_constrho()*d2pdrhodT());
	double drho_dhdT = d2hdrhodT()-1/pow(dpdrho_constT(),2)*(dpdrho_constT()*(dhdrho_constT()*d2pdrhodT()+d2hdrho2_constT()*dpdT_constrho())-dhdrho_constT()*dpdT_constrho()*d2pdrho2_constT());
	return ddT_dhdT-drho_dhdT*dpdT_constrho()/dpdrho_constT();
}
double CoolPropStateClass::d2hdp2_constT(void)
{
	return (d2hdrho2_constT()-dhdp_constT()*d2pdrho2_constT())/pow(dpdrho_constT(),2);
}
double CoolPropStateClass::d2hdTdp(void)
{
	return 1/dpdrho_constT()*(d2hdrhodT()-dhdp_constT()*(drhodT_constp()*d2pdrho2_constT()+d2pdrhodT())+d2hdrho2_constT()*drhodT_constp());
}

// Entropy
double CoolPropStateClass::dsdp_constT(void){
	return dsdrho_constT()/dpdrho_constT();
}
double CoolPropStateClass::dsdT_constp(void){
	return dsdT_constrho() - dsdrho_constT()*dpdT_constrho()/dpdrho_constT();
}
double CoolPropStateClass::dsdrho_constp(void){
	return dsdrho_constT() - dsdT_constrho()*dpdrho_constT()/dpdT_constrho();
}
double CoolPropStateClass::d2sdT2_constp(void)
{
	double ddT_dsdT = d2sdT2_constrho()-1/pow(dpdrho_constT(),2)*(dpdrho_constT()*(dsdrho_constT()*d2pdT2_constrho()+d2sdrhodT()*dpdT_constrho())-dsdrho_constT()*dpdT_constrho()*d2pdrhodT());
	double drho_dsdT = d2sdrhodT()-1/pow(dpdrho_constT(),2)*(dpdrho_constT()*(dsdrho_constT()*d2pdrhodT()+d2sdrho2_constT()*dpdT_constrho())-dsdrho_constT()*dpdT_constrho()*d2pdrho2_constT());
	return ddT_dsdT-drho_dsdT*dpdT_constrho()/dpdrho_constT();
}
double CoolPropStateClass::d2sdp2_constT(void)
{
	return (d2sdrho2_constT()-dsdp_constT()*d2pdrho2_constT())/pow(dpdrho_constT(),2);
}
double CoolPropStateClass::d2sdTdp(void)
{
	return 1/dpdrho_constT()*(d2sdrhodT()-dsdp_constT()*(drhodT_constp()*d2pdrho2_constT()+d2pdrhodT())+d2sdrho2_constT()*drhodT_constp());
}


double CoolPropStateClass::drhodp_constT(void)
{
	return 1/dpdrho_constT();
}
double CoolPropStateClass::d2rhodp2_constT(void)
{
	return -d2pdrho2_constT()/pow(dpdrho_constT(),3);
}
double CoolPropStateClass::d2rhodTdp(void)
{
	return (dpdT_constrho()*d2pdrho2_constT()-dpdrho_constT()*d2pdrhodT())/pow(dpdrho_constT(),3);
}
double CoolPropStateClass::d2rhodT2_constp(void)
{
	double ddrho_drhodT_p_constT = (dpdT_constrho()*d2pdrho2_constT()-dpdrho_constT()*d2pdrhodT())/pow(dpdrho_constT(),2);
	double ddT_drhodT_p_constrho = (dpdT_constrho()*d2pdrhodT()-dpdrho_constT()*d2pdT2_constrho())/pow(dpdrho_constT(),2);
	return ddT_drhodT_p_constrho+ddrho_drhodT_p_constT*drhodT_constp();
}

/// SATURATION DERIVATIVES
/// SATURATION DERIVATIVES
/// SATURATION DERIVATIVES

double CoolPropStateClass::dTdp_along_sat(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	if (pFluid->enabled_TTSE_LUT)
	{
		pFluid->build_TTSE_LUT();
		return pFluid->TTSESatL.evaluate_sat_derivative(iT,psatL);
	}
	else{
		return _T*(1/SatV->rho()-1/SatL->rho())/(SatV->h()-SatL->h());
	}
}
double CoolPropStateClass::ddp_dTdp_along_sat(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	return 1/(SatV->h()-SatL->h())*(_T*(SatV->dvdp_constT()-SatL->dvdp_constT())-dTdp_along_sat()*(SatV->dhdp_constT()-SatL->dhdp_constT()));
}
double CoolPropStateClass::ddT_dTdp_along_sat(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	return 1/(SatV->h()-SatL->h())*(_T*(SatV->dvdT_constp()-SatL->dvdT_constp())-dTdp_along_sat()*(SatV->dhdT_constp()-SatL->dhdT_constp())+(1/SatV->rho()-1/SatL->rho()));
}
double CoolPropStateClass::d2Tdp2_along_sat(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	return ddp_dTdp_along_sat()+ddT_dTdp_along_sat()*dTdp_along_sat();
}

double CoolPropStateClass::dhdp_along_sat_liquid(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	if (pFluid->enabled_TTSE_LUT){
		pFluid->build_TTSE_LUT();
		return pFluid->TTSESatL.evaluate_sat_derivative(iH,psatL);
	}
	else{
		return SatL->dhdp_constT()+SatL->dhdT_constp()*dTdp_along_sat();
	}
}
double CoolPropStateClass::dhdp_along_sat_vapor(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	if (pFluid->enabled_TTSE_LUT){
		pFluid->build_TTSE_LUT();
		return pFluid->TTSESatV.evaluate_sat_derivative(iH,psatV);
	}
	else{
		return SatV->dhdp_constT()+SatV->dhdT_constp()*dTdp_along_sat();
	}
}
double CoolPropStateClass::d2hdp2_along_sat_vapor(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	double ddp_dhdpsigmaV = SatV->d2hdp2_constT()+SatV->dhdT_constp()*ddp_dTdp_along_sat()+SatV->d2hdTdp()*dTdp_along_sat();
	double ddT_dhdpsigmaV = SatV->d2hdTdp()+SatV->dhdT_constp()*ddT_dTdp_along_sat()+SatV->d2hdT2_constp()*dTdp_along_sat();
	return ddp_dhdpsigmaV+ddT_dhdpsigmaV*dTdp_along_sat();
}
double CoolPropStateClass::d2hdp2_along_sat_liquid(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	double ddp_dhdpsigmaL = SatL->d2hdp2_constT()+SatL->dhdT_constp()*ddp_dTdp_along_sat()+SatL->d2hdTdp()*dTdp_along_sat();
	double ddT_dhdpsigmaL = SatL->d2hdTdp()+SatL->dhdT_constp()*ddT_dTdp_along_sat()+SatL->d2hdT2_constp()*dTdp_along_sat();
	return ddp_dhdpsigmaL+ddT_dhdpsigmaL*dTdp_along_sat();
}

double CoolPropStateClass::dsdp_along_sat_liquid(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	return SatL->dsdp_constT()+SatL->dsdT_constp()*dTdp_along_sat();
}
double CoolPropStateClass::dsdp_along_sat_vapor(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	return SatV->dsdp_constT()+SatV->dsdT_constp()*dTdp_along_sat();
}
double CoolPropStateClass::d2sdp2_along_sat_vapor(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	double ddp_dsdpsigmaV = SatV->d2sdp2_constT()+SatV->dsdT_constp()*ddp_dTdp_along_sat()+SatV->d2sdTdp()*dTdp_along_sat();
	double ddT_dsdpsigmaV = SatV->d2sdTdp()+SatV->dsdT_constp()*ddT_dTdp_along_sat()+SatV->d2sdT2_constp()*dTdp_along_sat();
	return ddp_dsdpsigmaV+ddT_dsdpsigmaV*dTdp_along_sat();
}
double CoolPropStateClass::d2sdp2_along_sat_liquid(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	double ddp_dsdpsigmaL = SatL->d2sdp2_constT()+SatL->dsdT_constp()*ddp_dTdp_along_sat()+SatL->d2sdTdp()*dTdp_along_sat();
	double ddT_dsdpsigmaL = SatL->d2sdTdp()+SatL->dsdT_constp()*ddT_dTdp_along_sat()+SatL->d2sdT2_constp()*dTdp_along_sat();
	return ddp_dsdpsigmaL+ddT_dsdpsigmaL*dTdp_along_sat();
}

double CoolPropStateClass::drhodp_along_sat_vapor(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	if (pFluid->enabled_TTSE_LUT)
	{
		pFluid->build_TTSE_LUT();
		return pFluid->TTSESatV.evaluate_sat_derivative(iD,psatV);
	}
	else
	{
		return SatV->drhodp_constT()+SatV->drhodT_constp()*dTdp_along_sat();
	}
}
double CoolPropStateClass::drhodp_along_sat_liquid(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	if (pFluid->enabled_TTSE_LUT)
	{
		pFluid->build_TTSE_LUT();
		return pFluid->TTSESatL.evaluate_sat_derivative(iD,psatL);
	}
	else
	{
		return SatL->drhodp_constT()+SatL->drhodT_constp()*dTdp_along_sat();
	}
}
double CoolPropStateClass::d2rhodp2_along_sat_vapor(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	double ddp_drhodpsigmaV = SatV->d2rhodp2_constT()+SatV->drhodT_constp()*ddp_dTdp_along_sat()+SatV->d2rhodTdp()*dTdp_along_sat();
	double ddT_drhodpsigmaV = SatV->d2rhodTdp()+SatV->drhodT_constp()*ddT_dTdp_along_sat()+SatV->d2rhodT2_constp()*dTdp_along_sat();
	return ddp_drhodpsigmaV+ddT_drhodpsigmaV*dTdp_along_sat();
}
double CoolPropStateClass::d2rhodp2_along_sat_liquid(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	double ddp_drhodpsigmaL = SatL->d2rhodp2_constT()+SatL->drhodT_constp()*ddp_dTdp_along_sat()+SatL->d2rhodTdp()*dTdp_along_sat();
	double ddT_drhodpsigmaL = SatL->d2rhodTdp()+SatL->drhodT_constp()*ddT_dTdp_along_sat()+SatL->d2rhodT2_constp()*dTdp_along_sat();
	return ddp_drhodpsigmaL+ddT_drhodpsigmaL*dTdp_along_sat();
}

double CoolPropStateClass::drhodT_along_sat_vapor(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	return SatV->drhodT_constp()+SatV->drhodp_constT()/dTdp_along_sat();
}
double CoolPropStateClass::drhodT_along_sat_liquid(void)
{
	if (!TwoPhase){throw ValueError(format("Saturation derivative cannot be called now.  Call update() with a two-phase set of inputs"));}
	return SatL->drhodT_constp()+SatL->drhodp_constT()/dTdp_along_sat();
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
	double dhdrhoL_T = SatL->dhdrho_constT();
	double dhdrhoV_T = SatV->dhdrho_constT();
	double dpdrhoL_T = SatL->dpdrho_constT();
	double dpdrhoV_T = SatV->dpdrho_constT();
	
	// Saturation derivatives at constant density
	double dhdTL_rho = SatL->dhdT_constrho();
	double dhdTV_rho = SatV->dhdT_constrho();
	double dpdTL_rho = SatL->dpdT_constrho();
	double dpdTV_rho = SatV->dpdT_constrho();

	// Now get dh/dp at constant T for saturated liquid and vapor
	double dhdpL_T = SatL->dhdp_constT();
	double dhdpV_T = SatV->dhdp_constT();

	// Derivatives of enthalpy 
	double dhdTL_p = SatL->dhdT_constp();
	double dhdTV_p = SatV->dhdT_constp();

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
	double d2Tsigma_dp2_T = ddp_dTdp_along_sat();  //T*((hV-hL)*(dvdpV_T-dvdpL_T)-(vV-vL)*(dhdpV_T-dhdpL_T))/pow(hV-hL,2);
	double d2Tsigma_dpdT_p = ddT_dTdp_along_sat();  //T*((hV-hL)*(dvdTV_p-dvdTL_p)-(vV-vL)*(dhdTV_p-dhdTL_p))/pow(hV-hL,2)+(vV-vL)/(hV-hL);
	double d2Tsigmadp2 = d2Tsigma_dp2_T+d2Tsigma_dpdT_p*dTsigmadp;

	double d2pdrho2V_T = SatV->d2pdrho2_constT();
	double d2pdrhodTV = SatV->d2pdrhodT();
	double d2pdT2V_rho = SatV->d2pdT2_constrho();
	double d2hdrho2V_T = SatV->d2hdrho2_constT();
	double d2hdrhodTV = SatV->d2hdrhodT();
	double d2hdT2V_rho = SatV->d2hdT2_constrho();
	double d2hdp2V_T = SatV->d2hdp2_constT();
	double d2hdTdpV = SatV->d2hdTdp();
	double d2hdT2V_p = sat->SatV->d2hdT2_constp();

	// derivative of dhdp along saturation with respect to p with constant T
	double ddp_dhdpsigmaV = d2hdp2V_T+dhdTV_p*d2Tsigma_dp2_T+d2hdTdpV*dTsigmadp;
	double ddT_dhdpsigmaV = d2hdTdpV+dhdTV_p*d2Tsigma_dpdT_p+d2hdT2V_p*dTsigmadp;
	*d2hdp2V = ddp_dhdpsigmaV+ddT_dhdpsigmaV*dTsigmadp;

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