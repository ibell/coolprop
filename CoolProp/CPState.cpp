
#include "CPExceptions.h"
#include "Solvers.h"
#include "CPState.h"
#include "float.h"
#include "math.h"
#include "Spline.h" 

#ifndef __ISWINDOWS__
	#ifndef DBL_EPSILON
		#include <limits>
		#define DBL_EPSILON std::numeric_limits<double>::epsilon()
	#endif
#endif

#include <stdio.h>

// Constructor with fluid name
CoolPropStateClass::CoolPropStateClass(std::string Fluid){
	// If a refprop fluid, add the fluid to the list of fluids
	if (Fluid.find("REFPROP-")==0){ add_REFPROP_fluid(Fluid); }

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
	_noSatLSatV = false;
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
	_noSatLSatV = false;
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
	|  H,S
	|  P,D
	*/
	bool using_EOS = true;

	if (debug()>3){
		std::cout<<__FILE__<<" update: "<<iInput1<<","<<Value1<<","<<iInput2<<","<<Value2<<","<<pFluid->get_name().c_str()<<std::endl;
	}

	// Reset all the internal variables to _HUGE
	_T = _HUGE;
	_p = _HUGE;
	_logp = _HUGE;
	_h = _HUGE;
	_s = _HUGE;
	_rho = _HUGE;
	_logrho = _HUGE;
	_Q = _HUGE;

	// Reset the cached values for _h and _s
	s_cached = false;
	h_cached = false;

	// Only build the Saturation classes if this is a top-level CPState for which no_SatLSatV() has not been called
	if (!_noSatLSatV){
		if (SatL == NULL){
			SatL = new CoolPropStateClass(pFluid);
			SatL->no_SatLSatV(); // Kill the recursive building of the saturation classes
		}
		if (SatV == NULL){
			SatV = new CoolPropStateClass(pFluid);
			SatV->no_SatLSatV(); // Kill the recursive building of the saturation classes
		}
	}

	// Pseudo-pure fluids cannot use T,Q as inputs if Q is not 0 or 1
	if (!pFluid->pure() && match_pair(iInput1,iInput2,iT,iQ))
	{
		// Sort so they are in the order T, Q
		sort_pair(&iInput1,&Value1,&iInput2,&Value2,iT,iQ);
		if (!(fabs(Value2) < 10*DBL_EPSILON) && !(fabs(Value2-1) < 10*DBL_EPSILON)){
			throw ValueError(format("Pseudo-pure fluids cannot use temperature-quality as inputs if Q is not 1 or 0"));
		}
	}

	// Don't know if it is single phase or not, so assume it isn't
	SinglePhase = false;
	
	// Determine whether the EOS or the TTSE will be used
	if (!pFluid->enabled_TTSE_LUT)
	{
		using_EOS = true;
	}
	else
	{
		// Try to build the LUT; Nothing will happen if the tables are already built
		pFluid->build_TTSE_LUT();

		// If inputs are in range of LUT, use it, otherwise just use the EOS
		if (within_TTSE_range(iInput1,Value1,iInput2,Value2)){
			using_EOS = false;
		}
		else
		{
			using_EOS = true;
		}
	}

	if (using_EOS)
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
		else if (match_pair(iInput1,iInput2,iP,iD)){
			update_prho(iInput1,Value1,iInput2,Value2);
		}
		else if (match_pair(iInput1,iInput2,iH,iS)){
			update_hs(iInput1,Value1,iInput2,Value2);
		}
		else
		{
			throw ValueError(format("Sorry your inputs didn't work; valid pairs are P,Q T,Q T,D T,P P,H P,S"));
		}
		// Clear the cached derivative flags
		this->clear_cache();	
	}
	else
	{
		// Update using the LUT
		update_TTSE_LUT(iInput1, Value1, iInput2, Value2);
		// Calculate the log of the pressure since a lot of terms need it and log() is a very slow function
		if (!ValidNumber(_logp)) _logp = log(_p);
		if (!ValidNumber(_logrho)) _logrho = log(_rho);
	}
	
	if (!_noSatLSatV && TwoPhase && !flag_TwoPhase)
	{
		// Update temperature and density for SatL and SatV
		add_saturation_states();
	}

	// Reduced parameters
	delta = this->_rho/pFluid->reduce.rho;
	tau = pFluid->reduce.T/this->_T;
}

bool CoolPropStateClass::within_TTSE_range(long iInput1, double Value1, long iInput2, double Value2)
{
	// For now, only allow p,h to use values outside of the TTSE table
	if (match_pair(iInput1,iInput2,iP,iH)){
		// Sort in the order p,h
		sort_pair(&iInput1,&Value1,&iInput2,&Value2,iP,iH);

		double hmin = 0, hmax = 0, pmin = 0, pmax = 0;
		pFluid->get_TTSESinglePhase_LUT_range(&hmin,&hmax,&pmin,&pmax);
		return (Value1 > pmin && Value1 < pmax && Value2 > hmin && Value2 < hmax);
	}
	return true;
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

		// Out-of-range checks
		if (Value1 < pFluid->params.ptriple-100*DBL_EPSILON || Value1 > pFluid->crit.p+100*DBL_EPSILON){ throw ValueError(format("Your saturation pressure [%f kPa] is out of range [%f kPa, %f kPa]",Value1,pFluid->params.ptriple,pFluid->crit.p ));}
		if (Value2 > 1+10*DBL_EPSILON || Value2 < -10*DBL_EPSILON){ throw ValueError(format("Your quality [%f] is out of range (0, 1)",Value2 )); }

		// Carry out the saturation call to get the temperature and density for each phases
		if (pFluid->pure()){
			pFluid->saturation_p(Value1,pFluid->enabled_TTSE_LUT,&TsatL,&TsatV,&rhosatL,&rhosatV);
			TsatV = TsatL;
			psatV = Value1;
			psatL = Value1;
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
		
		// Out-of-range checks
		if (Value1 < pFluid->limits.Tmin-10*DBL_EPSILON || Value1 > pFluid->crit.T+10*DBL_EPSILON){ throw ValueError(format("Your saturation temperature [%f K] is out of range [%f K, %f K]",Value1,pFluid->limits.Tmin, pFluid->reduce.T ));}
		if (Value2 > 1+10*DBL_EPSILON || Value2 < -10*DBL_EPSILON){ throw ValueError(format("Your quality [%f] is out of range (0, 1)",Value2 )); }
		
		// Carry out the saturation call to get the temperature and density for each phases
		// for the given temperature
		if (pFluid->pure()){
			pFluid->saturation_T(Value1,pFluid->enabled_TTSE_LUT,&psatL,&psatV,&rhosatL,&rhosatV);
			TsatL = Value1;
			TsatV = Value1;
		}
		else{
			TsatL = Value1;
			TsatV = Value1;
			// Saturation pressures
			psatL = pFluid->psatL_anc(TsatL);
			psatV = pFluid->psatV_anc(TsatV);

			try{
			// Saturation densities
			rhosatL = pFluid->density_Tp(TsatL, psatL, pFluid->rhosatL(TsatL));
			rhosatV = pFluid->density_Tp(TsatV, psatV, pFluid->rhosatV(TsatV));
			}
			catch (std::exception){
				// Near the critical point, the behavior is not very nice, so we will just use the ancillary near the critical point
				rhosatL = pFluid->rhosatL(TsatL);
				rhosatV = pFluid->rhosatV(TsatV);
			}
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

	if (Value1 < 0 ){ throw ValueError(format("Your temperature [%f K] is less than zero",Value1));}
	if (Value2 < 0 ){ throw ValueError(format("Your density [%f kg/m^3] is less than zero",Value2));}

	if (flag_SinglePhase && flag_TwoPhase) throw ValueError(format("Only one of flag_SinglePhase and flag_TwoPhase may be set to true"));

	// If either SinglePhase or flag_SinglePhase is set to true, it will not make the call to the saturation routine
	// SinglePhase is set by the class routines, and flag_SinglePhase is a flag that can be set externally
	bool _TwoPhase;

	if (flag_SinglePhase || SinglePhase){
		_TwoPhase = false;
	}
	else if(flag_TwoPhase){
		_TwoPhase = true;
	}
	else{
		_TwoPhase = !pFluid->phase_Trho(_T,_rho,&psatL,&psatV,&rhosatL,&rhosatV).compare("Two-Phase");
	}

	if (_TwoPhase)
	{
		if (flag_TwoPhase){
			pFluid->phase_Trho(_T,_rho,&psatL,&psatV,&rhosatL,&rhosatV);
		}
		TsatV = _T;
		TsatL = _T;
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

// Updater if p,rho are inputs
void CoolPropStateClass::update_prho(long iInput1, double Value1, long iInput2, double Value2)
{
	// Get them in the right order
	sort_pair(&iInput1,&Value1,&iInput2,&Value2,iP,iD);

	// Set internal variables
	_p = Value1;
	_rho = Value2;

	if (Value1 < 0 ){ throw ValueError(format("Your pressure [%f kPa] is less than zero",Value1));}
	if (Value2 < 0 ){ throw ValueError(format("Your density [%f kg/m^3] is less than zero",Value2));}

	if (flag_SinglePhase && flag_TwoPhase) throw ValueError(format("Only one of flag_SinglePhase and flag_TwoPhase may be set to true"));

	// If either SinglePhase or flag_SinglePhase is set to true, it will not make the call to the saturation routine
	// SinglePhase is set by the class routines, and flag_SinglePhase is a flag that can be set externally
	bool _TwoPhase;

	if (flag_SinglePhase || SinglePhase){
		_TwoPhase = false;
	}
	else if(flag_TwoPhase){
		_TwoPhase = true;
	}
	else{
		_TwoPhase = (pFluid->phase_prho_indices(_p,_rho,&_T,&TsatL,&TsatV,&rhosatL,&rhosatV) == iTwoPhase);
	}

	if (_TwoPhase)
	{
		if (flag_TwoPhase){
			pFluid->phase_prho_indices(_p,_rho,&_T,&TsatL,&TsatV,&rhosatL,&rhosatV);
		}
		// If it made it to the saturation routine and it is two-phase the saturation variables have been set
		TwoPhase = true;
		SinglePhase = false;

		// Get the quality and pressure
		_Q = (1/_rho-1/rhosatL)/(1/rhosatV-1/rhosatL);
		
		check_saturated_quality(_Q);
	}
	else{
		TwoPhase = false;
		SinglePhase = true;
		SaturatedL = false;
		SaturatedV = false;
		if (!ValidNumber(_T)){
			double T0 = pFluid->temperature_prho_PengRobinson(_p,_rho);
			_T = pFluid->temperature_prho(_p, _rho, T0);
		}
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
	bool _TwoPhase;

	if (flag_SinglePhase || SinglePhase){
		_TwoPhase = false;
	}
	else if(flag_TwoPhase){
		_TwoPhase = true;
	}
	else{
		_TwoPhase = !pFluid->phase_Tp(_T,_p,&psatL,&psatV,&rhosatL,&rhosatV).compare("Two-Phase");
	}

	if (_TwoPhase)
	{
		throw ValueError(format("TwoPhase is not possible with T,P as inputs"));
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

	if (Value1 < 0 ){ throw ValueError(format("Your pressure [%g kPa] is less than zero",Value1));}

	// Solve for temperature and density with or without the guess values provided
	pFluid->temperature_ph(Value1, Value2,&_T,&_rho,&rhosatL,&rhosatV,&TsatL,&TsatV, T0, rho0);

	// Set internal variables
	_p = Value1;
	_h = Value2;
	h_cached = true;

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

	if (Value1 < 0 ){ throw ValueError(format("Your pressure [%g kPa] is less than zero",Value1));}

	// Solve for temperature and density
	pFluid->temperature_ps(Value1, Value2,&_T,&_rho,&rhosatL,&rhosatV,&TsatL,&TsatV);

	// Set internal variables
	_p = Value1;
	_s = Value2;
	s_cached = true;

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

// Updater if h,s are inputs
void CoolPropStateClass::update_hs(long iInput1, double Value1, long iInput2, double Value2)
{
	// Get them in the right order
	sort_pair(&iInput1,&Value1,&iInput2,&Value2,iH,iS);

	// Solve for temperature and density
	pFluid->temperature_hs(Value1, Value2,&_T,&_rho,&rhosatL,&rhosatV,&TsatL,&TsatV);

	// Set internal variables
	_h = Value1;
	_s = Value2;
	h_cached = true;
	s_cached = true;
	_p = pFluid->pressure_Trho(_T,_rho); // Evaluate the EOS directly to find pressure without a saturation call

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
	else if (match_pair(iInput1,iInput2,iT,iQ))
	{
		// We know it is a pure fluid since pseudo-pure can't get this far (see update function)
		// Set phase flags
		SinglePhase = false;
		TwoPhase = true;
		
		// Sort in the right order (T,Q)
		sort_pair(&iInput1,&Value1,&iInput2,&Value2,iT,iQ);

		double T = Value1;
		double Q = Value2;

		// Check whether saturated liquid or vapor
		check_saturated_quality(Q);

		double p = pFluid->TTSESatL.evaluate_T(T);
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

		// If enthalpy is outside the saturation region or flag_SinglePhase is set, it is single-phase
		if (p > pFluid->reduce.p || p < pFluid->params.ptriple || flag_SinglePhase ||  h < pFluid->TTSESatL.evaluate(iH,p)  || h > pFluid->TTSESatV.evaluate(iH,p))
		{
			TwoPhase = false;
			SinglePhase = true;
			_logp = log(p);
			_rho = pFluid->TTSESinglePhase.evaluate(iD,p,_logp,h);
			_T = pFluid->TTSESinglePhase.evaluate(iT,p,_logp,h);
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
	else if (match_pair(iInput1,iInput2,iP,iT))
	{
		SinglePhase = true;
		TwoPhase = false;
		// Sort in the right order (P,T)
		sort_pair(&iInput1,&Value1,&iInput2,&Value2,iP,iT);

		_logp = log(Value1);
		_h = pFluid->TTSESinglePhase.evaluate_one_other_input(iP,Value1,iT,Value2);
		_rho = pFluid->TTSESinglePhase.evaluate(iD,Value1,_logp, _h);
		_p = Value1;
		_T = Value2;
	}
	else if (match_pair(iInput1,iInput2,iP,iD))
	{
		// Sort in the right order (P,D)
		sort_pair(&iInput1,&Value1,&iInput2,&Value2,iP,iD);

		double p = Value1;
		double rho = Value2;

		// If density is outside the saturation region, it is single-phase
		if (p > pFluid->reduce.p || p < pFluid->params.ptriple ||  rho < pFluid->TTSESatV.evaluate(iD,p)  || rho > pFluid->TTSESatL.evaluate(iD,p))
		{
			TwoPhase = false;
			SinglePhase = true;
			_logp = log(p);
			// Saturation calls happen again inside evaluate_one_other_input - perhaps pass values
			_h = pFluid->TTSESinglePhase.evaluate_one_other_input(iP,p,iD,rho);
			_T = pFluid->TTSESinglePhase.evaluate(iT,p,_logp,_h);
			_p = Value1;
			_rho = Value2;
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
			
			_Q = (1/rho-1/rhosatL)/(1/rhosatV-1/rhosatL);
			_rho = rho;
			_T = _Q*TsatV + (1-_Q)*TsatL;
			_p = p;
			_h = _Q*hsatV + (1-_Q)*hsatL;

			check_saturated_quality(_Q);
		}
	}
	else if (match_pair(iInput1,iInput2,iP,iS))
	{
		// Sort in the right order (P,S)
		sort_pair(&iInput1,&Value1,&iInput2,&Value2,iP,iS);
		
		double p = Value1;
		double s = Value2;

		_p = p;
		_s = s;

		// If entropy is outside the saturation region, it is single-phase
		if (p > pFluid->reduce.p || p < pFluid->params.ptriple ||  s > pFluid->TTSESatV.evaluate(iS,p)  || s < pFluid->TTSESatL.evaluate(iS,p))
		{
			TwoPhase = false;
			SinglePhase = true;
			_logp = log(p);
			_h = pFluid->TTSESinglePhase.evaluate_one_other_input(iP,p,iS,s); // Get the enthalpy
			_T = pFluid->TTSESinglePhase.evaluate(iT,p,_logp,_h);
			_rho = pFluid->TTSESinglePhase.evaluate(iD,p,_logp,_h);
		}
		else
		{
			TwoPhase = true;
			SinglePhase = false;

			double hsatL = pFluid->TTSESatL.evaluate(iH,p);
			double hsatV = pFluid->TTSESatV.evaluate(iH,p);
			double ssatL = pFluid->TTSESatL.evaluate(iS,p);
			double ssatV = pFluid->TTSESatV.evaluate(iS,p);
			rhosatL = pFluid->TTSESatL.evaluate(iD,p);
			rhosatV = pFluid->TTSESatV.evaluate(iD,p);
			TsatL = pFluid->TTSESatL.evaluate(iT,p);
			TsatV = pFluid->TTSESatV.evaluate(iT,p);
			psatL = p;
			psatV = p;
			
			_Q = (s-ssatL)/(ssatV-ssatL);
			_rho = 1/(_Q/rhosatV+(1-_Q)/rhosatL);
			_T = _Q*TsatV + (1-_Q)*TsatL;
			_h = _Q*hsatV + (1-_Q)*hsatL;

			check_saturated_quality(_Q);
		}		
	}
	else if (match_pair(iInput1,iInput2,iT,iD))
	{
		// Sort in the right order (T,D)
		sort_pair(&iInput1,&Value1,&iInput2,&Value2,iT,iD);

		_T = Value1;
		_rho = Value2;
		_logrho = log(_rho);

		// If density is outside the saturation region, it is single-phase
		if (Value1 > pFluid->crit.T || Value1 < pFluid->params.Ttriple){
			SinglePhase = true;
		}
		else{
			double psatL = pFluid->TTSESatL.evaluate_T(Value1);
			double psatV = pFluid->TTSESatV.evaluate_T(Value1);
			if (Value2 < pFluid->TTSESatV.evaluate(iD,psatV)  || Value2 > pFluid->TTSESatL.evaluate(iD,psatL)){
				SinglePhase = true;
			}
			else{
				SinglePhase = false;
			}
		}
		
		if (SinglePhase)
		{
			TwoPhase = false;
			SinglePhase = true;

			_p = pFluid->TTSESinglePhase.evaluate_Trho(iP,Value1,Value2,_logrho);
		}
		else
		{
			TwoPhase = true;
			SinglePhase = false;
			_p = pFluid->TTSESatL.evaluate_T(Value1);
			rhosatL = pFluid->TTSESatL.evaluate(iD,_p);
			rhosatV = pFluid->TTSESatV.evaluate(iD,_p);
			TsatL = pFluid->TTSESatL.evaluate(iT,_p);
			TsatV = pFluid->TTSESatV.evaluate(iT,_p);
			psatL = _p;
			psatV = _p;
			
			_Q = (1/_rho-1/rhosatL)/(1/rhosatV-1/rhosatL);
		}
	}
	else
	{
		printf("Sorry your inputs[%d,%d] don't work for now with TTSE\n",(int)iInput1,(int)iInput2);
		throw ValueError(format("Sorry your inputs don't work for now with TTSE"));
	}
}

/// Return an output based on the integer key for the term
double CoolPropStateClass::keyed_output(long iOutput)
{
	switch (iOutput)
	{
		// --------------------------
		// Fluid constants
		// --------------------------
		case iMM:
			return pFluid->params.molemass;
		case iPcrit:
			return pFluid->crit.p;
		case iTcrit:
			return pFluid->crit.T;
		case iScrit:
			return pFluid->crit.s;
		case iHcrit:
			return pFluid->crit.h;
		case iTtriple:
			return pFluid->params.Ttriple;
		case iPtriple:
			return pFluid->params.ptriple;
		case iRhocrit:
			return pFluid->crit.rho;
		case iAccentric: 
			return pFluid->params.accentricfactor;
		case iTmin:
			return pFluid->limits.Tmin;
		case iCritSplineT:
			return pFluid->CriticalSpline_T.Tend;

		// --------------------------
		// Phase Constants
		// --------------------------
		case iPHASE_LIQUID:
			return iLiquid;
		case iPHASE_GAS:
			return iGas;
		case iPHASE_SUPERCRITICAL:
			return iSupercritical;
		case iPHASE_TWOPHASE:
			return iTwoPhase;

		// --------------------------
		// Environmental properties
		// --------------------------
		case iODP:
			return pFluid->environment.ODP;
		case iGWP20:
			return pFluid->environment.GWP20;
		case iGWP100:
			return pFluid->environment.GWP100;
		case iGWP500:
			return pFluid->environment.GWP500;

		// --------------------------
		// Thermodynamic properties
		// --------------------------
		case iT:
			return _T;
		case iD:
			return _rho;
		case iP:
			return _p;
		case iC:
			return cp();
		case iC0:
			return pFluid->specific_heat_p_ideal_Trho(_T);
		case iO:
			return cv();
		case iA:
			return speed_sound();
		case iG:
			return h()-_T*s();
		case iQ:
			if (TwoPhase)
				return _Q;
			else
				return -1;
		case iH:
			return h();
		case iS:
			return s();
		case iU:
			return h()-_p/_rho;

		case iPhase:
			return phase();
		
		// --------------------------
		// Transport properties
		// --------------------------
		case iV:
			return viscosity();
		case iL:
			return conductivity();
		case iI:
			return surface_tension();
		// -----------------------------------
		// A few grandfathered derivatives
		// -----------------------------------
		case iDpdT:
			return dpdT_constrho();
		case iDrhodT_p:
			return drhodT_constp();
		default:
			throw ValueError(format("Invalid Output index to CPState function keyed_output: %d ",iOutput));
	}
}

long CoolPropStateClass::phase(void)
{
	double pL,pV,rhoL,rhoV;
	return pFluid->phase_Trho_indices(_T,_rho,&pL,&pV,&rhoL,&rhoV);
}

void CoolPropStateClass::add_saturation_states(void)
{
	// While SatL and SatV are technically two-phase, we consider 
	// them to be single-phase to speed up the calcs and avoid saturation calls
	SatL->flag_SinglePhase = true;
	SatL->flag_TwoPhase = false;
	SatL->update(iT,TsatL,iD,rhosatL);
	SatL->flag_SinglePhase = false;
	SatL->SinglePhase = true;
	SatL->TwoPhase = false;

	SatV->flag_SinglePhase = true;
	SatV->update(iT,TsatV,iD,rhosatV);
	SatV->flag_SinglePhase = false;
	SatV->SinglePhase = true;
	SatV->TwoPhase = false;
}

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
		return pFluid->TTSESatL.evaluate(iS,psatL);
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

double CoolPropStateClass::cpL(void){return SatL->cp();};
double CoolPropStateClass::cpV(void){return SatV->cp();};
double CoolPropStateClass::viscL(void){return SatL->keyed_output(iV);};
double CoolPropStateClass::viscV(void){return SatV->keyed_output(iV);};
double CoolPropStateClass::condL(void){return SatL->keyed_output(iL);};
double CoolPropStateClass::condV(void){return SatV->keyed_output(iL);};

double CoolPropStateClass::h(void){
	if (TwoPhase){
		// This will use the TTSE LUT if enable_TTSE_LUT() has been called
		return _Q*hV()+(1-_Q)*hL();
	}
	else{
		if (h_cached && ValidNumber(_h)){
			// Use the pre-calculated value
			return _h;
		}
		else
		{
			if (pFluid->enabled_TTSE_LUT)
			{
				return pFluid->TTSESinglePhase.evaluate_Trho(iH,_T,_rho,_logrho);
			}
			else
			{
				// Use the EOS, using the cached value if possible
				return pFluid->R()*_T*(1.0+tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+delta*dphir_dDelta(tau,delta));
			}
		}
	}
}
double CoolPropStateClass::s(void){
	if (TwoPhase){
		// This will use the TTSE LUT if enable_TTSE_LUT() has been called
		return _Q*sV()+(1-_Q)*sL();
	}
	else{
		if (s_cached && ValidNumber(_s) && !pFluid->enabled_TTSE_LUT)
		{
			// Use the pre-calculated value
			return _s;
		}
		else
		{
			if (pFluid->enabled_TTSE_LUT && within_TTSE_range(iP, p(), iH, h()) )
			{
				return pFluid->TTSESinglePhase.evaluate_Trho(iS,_T,_rho,_logrho);
			}
			else
			{
				// Use the EOS, using the cached value if possible
				return pFluid->R()*(tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))-phi0(tau,delta)-phir(tau,delta));
			}
		}
	}
}
double CoolPropStateClass::cp(void){
	if (TwoPhase && _Q > 0 && _Q < 1)
	{
		return _HUGE;
	}
	else if (pFluid->enabled_TTSE_LUT && within_TTSE_range(iP, p(), iH, h()) )
	{
		// cp is also given by (dh/dT)|p, or 1/((dT/dh)|p) which is tabulated
		_h = h();
		return 1/pFluid->TTSESinglePhase.evaluate_first_derivative(iT,iH,iP,_p,_logp,_h);
	}
	else
	{
		double c1 = pow(1.0+delta*dphir_dDelta(tau,delta)-delta*tau*d2phir_dDelta_dTau(tau,delta),2);
		double c2 = (1.0+2.0*delta*dphir_dDelta(tau,delta)+pow(delta,2)*d2phir_dDelta2(tau,delta));
		return pFluid->R()*(-pow(tau,2)*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta))+c1/c2);
	}
}

double CoolPropStateClass::viscosity(void){
	if (pFluid->enabled_TTSE_LUT && within_TTSE_range(iP,p(),iH,h()))
	{
		return pFluid->TTSESinglePhase.evaluate_Trho(iV,_T,_rho,_logrho);
	}
	else
	{
		return pFluid->viscosity_Trho(_T,_rho);
	}
}

double CoolPropStateClass::conductivity(void){
	if (pFluid->enabled_TTSE_LUT  && within_TTSE_range(iP,p(),iH,h()))
	{
		return pFluid->TTSESinglePhase.evaluate_Trho(iL,_T,_rho,_logrho);
	}
	else{
		return pFluid->conductivity_Trho(_T,_rho);
	}
}

double CoolPropStateClass::B_TTSE(double p, double h){
	// Slightly modified, doesn't use specific volume at all, also sign switched
	double drhodh_p = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iH,iP,_p,_logp,h);
	double drhodp_h = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iP,iH,_p,_logp,h);
	double dsdh_p = pFluid->TTSESinglePhase.evaluate_first_derivative(iS,iH,iP,_p,_logp,h);
	double dsdp_h = pFluid->TTSESinglePhase.evaluate_first_derivative(iS,iP,iH,_p,_logp,h);
	return (drhodp_h*dsdh_p-drhodh_p*dsdp_h);
}
double CoolPropStateClass::B_over_D_TTSE(double p, double h)
{
	double drhodh_p = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iH,iP,_p,_logp,h);
	double drhodp_h = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iP,iH,_p,_logp,h);
	double   dsdh_p = pFluid->TTSESinglePhase.evaluate_first_derivative(iS,iH,iP,_p,_logp,h);
	double   dsdp_h = pFluid->TTSESinglePhase.evaluate_first_derivative(iS,iP,iH,_p,_logp,h);
	double   dTdh_p = pFluid->TTSESinglePhase.evaluate_first_derivative(iT,iH,iP,_p,_logp,h);
	double   dTdp_h = pFluid->TTSESinglePhase.evaluate_first_derivative(iT,iP,iH,_p,_logp,h);
	return (drhodh_p*dsdp_h-drhodp_h*dsdh_p)/(dTdh_p*drhodp_h-dTdp_h*drhodh_p);
}
double CoolPropStateClass::cv(void){
	if (pFluid->enabled_TTSE_LUT && within_TTSE_range(iP, p(), iH, h()) )
	{
		if (TwoPhase && _Q>0 && _Q < 1)
		{
			return -1;
		}
		else
		{
			// cv is also given by -B*T/D which is tabulated (indirectly)
			_h = h();
			return -B_over_D_TTSE(_p,_h)*pFluid->TTSESinglePhase.evaluate(iT,_p,_logp,_h);
		}
	}
	else
	{
		return -pFluid->R()*pow(tau,2)*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta));
	}
}
double CoolPropStateClass::speed_sound(void){
	if (pFluid->enabled_TTSE_LUT && within_TTSE_range(iP, p(), iH, h()) )
	{
		if (TwoPhase && _Q>0 && _Q < 1)
		{
			// Not defined for two-phase
			return -1;
		}
		else
		{
			_h = h();
			// speed of sound given by sqrt(v^2/B*dsdh|p), or sqrt(dpdrho|s)
			double dsdh__p = pFluid->TTSESinglePhase.evaluate_first_derivative(iS,iH,iP,_p,_logp,_h);
			double dsdp__h = pFluid->TTSESinglePhase.evaluate_first_derivative(iS,iP,iH,_p,_logp,_h);
			double drhodh__p = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iH,iP,_p,_logp,_h);
			double drhodp__h = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iP,iH,_p,_logp,_h);
			
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

	if (pFluid->enabled_TTSE_LUT && within_TTSE_range(iP, p(), iH, h()) )
	{
		if (TwoPhase && _Q>0 && _Q < 1)
		{
			// Not defined for two-phase
			return -1;
		}
		else
		{
			_h = h();
			// isothermal compressibility given by kappa = -1/v*dvdp|T = 1/rho*drhodp|T
			double dTdh__p = pFluid->TTSESinglePhase.evaluate_first_derivative(iT,iH,iP,_p,_logp,_h);
			double dTdp__h = pFluid->TTSESinglePhase.evaluate_first_derivative(iT,iP,iH,_p,_logp,_h);
			double drhodh__p = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iH,iP,_p,_logp,_h);
			double drhodp__h = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iP,iH,_p,_logp,_h);

			double rho = pFluid->TTSESinglePhase.evaluate(iD,_p,_logp,_h);
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

	if (pFluid->enabled_TTSE_LUT && within_TTSE_range(iP, p(), iH, h()) )
	{
		if (TwoPhase && _Q>0 && _Q < 1)
		{
			// Not defined for two-phase
			return -1;
		}
		else
		{
			_h = h();
			// isobaric expansion coefficient given by kappa = 1/v*dvdT|p = -1/rho*drhodT|p
			double dTdh__p = pFluid->TTSESinglePhase.evaluate_first_derivative(iT,iH,iP,_p,_logp,_h);
			double drhodh__p = pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iH,iP,_p,_logp,_h);
			double rho = pFluid->TTSESinglePhase.evaluate(iD,_p,_logp,_h);
			return -1/rho*drhodh__p/dTdh__p;
		}
	}
	else
	{
		return -1/(_rho*_rho)*drhodT_constp();
	}
}
double CoolPropStateClass::surface_tension(void){
	return pFluid->surface_tension_T(_T);
}

double CoolPropStateClass::drhodh_constp(void){

	if (pFluid->enabled_TTSE_LUT && within_TTSE_range(iP,p(),iH,h()) )
	{
		if (TwoPhase && _Q>0 && _Q < 1)
		{
			// equals -rho^2*dvdh_p where dvdh_p = 1/T*dTdp|sat
			return -_rho*_rho/_T*pFluid->TTSESatL.evaluate_sat_derivative(iT,_p);
		}
		else
		{
			_h = h();
			return pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iH,iP,_p,_logp,_h);
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

double CoolPropStateClass::drhodp_consth_smoothed(double xend){
	// Make a state class instance
	CoolPropStateClass CPS = CoolPropStateClass(pFluid);
	SplineClass SC = SplineClass();
	double hL = this->hL();
	double hV = this->hV();
	double hend = xend*hV+(1-xend)*hL;

	// We do the interpolation in terms of enthalpy because it is a bit simpler to do
	// The values at x = 0, but faking out CoolProp by using T,rho and enforcing singlephase
	CPS.flag_SinglePhase = true;
	if (isenabled_TTSE_LUT()){
		CPS.update(iP,psatL,iH,hL);
	}
	else{
		CPS.update(iT,TsatL,iD,rhosatL);	
	}
	SC.add_value_constraint(hL, CPS.drhodp_consth());
	SC.add_derivative_constraint(hL, CPS.d2rhodhdp());
	// The values at x = xend
	CPS.flag_SinglePhase = false;
	CPS.update(iT,TsatL,iQ,xend);
	SC.add_value_constraint(hend, CPS.drhodp_consth());
	double d2vdhdp = 1/CPS.T()*CPS.d2Tdp2_along_sat()-1*pow(CPS.dTdp_along_sat()/CPS.T(),2);
	SC.add_derivative_constraint(hend, 2/CPS.rho()*CPS.drhodp_consth()*CPS.drhodh_constp()-pow(CPS.rho(),2)*d2vdhdp);
	SC.build();
	return SC.evaluate(_Q*hV+(1-_Q)*hL);
}

double CoolPropStateClass::drhodh_constp_smoothed(double xend){
	// Make a state class instance
	CoolPropStateClass CPS = CoolPropStateClass(pFluid);
	SplineClass SC = SplineClass();
	double hL = this->hL();
	double hV = this->hV();
	double hend = xend*hV+(1-xend)*hL;

	// We do the interpolation in terms of enthalpy because it is a bit simpler to do
	// The values at x = 0, but faking out CoolProp by using T,rho and enforcing singlephase
	CPS.flag_SinglePhase = true;
	// Get the value at the saturated liquid part
	if (isenabled_TTSE_LUT()){
		CPS.update(iP,psatL,iH,hL);
	}
	else{
		CPS.update(iT,TsatL,iD,rhosatL);	
	}
	SC.add_value_constraint(hL, CPS.drhodh_constp());
	SC.add_derivative_constraint(hL, CPS.d2rhodh2_constp());
	// The values at x = xend
	CPS.update(iT,TsatL,iQ,xend);
	SC.add_value_constraint(hend, CPS.drhodh_constp());
	SC.add_derivative_constraint(hend, 2/CPS.rho()*CPS.drhodh_constp()*CPS.drhodh_constp());
	SC.build();
	return SC.evaluate(_Q*hV+(1-_Q)*hL);
}


void CoolPropStateClass::rho_smoothed(double xend, double *rho_spline, double *dsplinedh, double *dsplinedp){
	// Make a state class instance in two-phase at the junction point (end):
	CoolPropStateClass CPS = CoolPropStateClass(pFluid);
	CPS.update(iT,TsatL,iQ,xend);

	// A few necessary properties:
	double h_l = this->hL();
	double h_v = this->hV();
	double rho_l = this->rhoL();
	double rho_v = this->rhoV();
	double x = _Q;
	double h_end = xend * h_v + (1-xend)*h_l;
	double rho_end = (rho_l * rho_v)/(xend * rho_l + (1-xend)*rho_v);

	// Getting the total derivatives along the saturation lines:
	double drholdp = CPS.drhodp_along_sat_liquid();
	double dhldp = CPS.dhdp_along_sat_liquid();
	double drhovdp = CPS.drhodp_along_sat_vapor();
	double dhvdp = CPS.dhdp_along_sat_vapor();

	// Getting the required partial derivatives just outside the two-phase zone (faking single-phase fluid):
	double drhodh_l = CPS.SatL->drhodh_constp();
	double drhodhdp_l = CPS.SatL->d2rhodhdp();
	double drhodh_v = CPS.SatV->drhodh_constp();

	// Partial derivatives at the junction (end):
	double drhodh_end_temp = CPS.drhodh_constp();
	double drhodp_end_temp = CPS.drhodp_consth();
	double drhodhdp_end_temp = CPS.d2rhodhdp();

	// Same as above, but detailed:
	double dxdp = ((xend - 1 )* dhldp - xend* dhvdp)/(h_v - h_l);
	double drhodh_end = pow(rho_end,2)/(rho_l*rho_v) * (rho_v - rho_l)/(h_v - h_l);
	double dvdh_end = (1/rho_v - 1/rho_l)/(h_v - h_l);
	double dvdp_end = (-1/pow(rho_l,2) * drholdp + dxdp * (1/rho_v - 1/rho_l) + xend * (-1/pow(rho_v,2) * drhovdp + 1/pow(rho_l,2) * drholdp));
	double drhodp_end = -pow(rho_end,2) * dvdp_end;
	double dvdhdp_end = 1/(h_v - h_l) * (-1/pow(rho_v,2)*drhovdp + 1/pow(rho_l,2) * drholdp) - (1/rho_v - 1/rho_l) / pow((h_v - h_l),2) * (dhvdp - dhldp);
	double drhodhdp_end = -2 * rho_end * dvdh_end * drhodp_end -pow(rho_end,2)*dvdhdp_end;

	// Derivative at constant quality (xend):
	double drhoxdp =  pow(rho_end,2)*(xend / pow(rho_v,2) * drhovdp + (1-xend)/pow(rho_l,2) * drholdp);

	//Checking that the derivatives are ok:
	//double pend = CPS._p;
	//double step = 0.1;
	//CPS.update(iH,h_end+step,iP,pend);
	//double rho_hplus = CPS._rho;
	//CPS.update(iH,h_end+step,iP,pend+step);
	//double rho_hplus_pplus = CPS._rho;
	//CPS.update(iH,h_end,iP,pend+step);
	//double rho_pplus = CPS._rho;
	//double dvdh_end_pplus = (1/CPS.rhoV() - 1/CPS.rhoL())/(CPS.hV() - CPS.hL());
	//double drhodh_check_pplus = -rho_pplus*rho_pplus * dvdh_end_pplus;
	//double drhodh_check2_pplus = CPS.drhodh_constp();
	//CPS.update(iH,h_end+step,iP,pend-step);
	//double rho_hplus_pminus = CPS._rho;
	//CPS.update(iH,h_end-step,iP,pend+step);
	//double rho_hminus_pplus = CPS._rho;
	//CPS.update(iH,h_end-step,iP,pend-step);
	//double rho_hminus_pminus = CPS._rho;
	//double dvdhdp_check = (dvdh_end_pplus - dvdh_end)/step;
	//double drhodhdp_check = -rho_end*rho_end * dvdhdp_check;
	//double drhodhdp_check2 = (drhodh_check2_pplus - drhodh_end)/step;
	//double drhodp_hend = (rho_pplus - rho_end)/step;
	//double drhodp_hplus = (rho_hplus_pplus - rho_hplus)/step;
	//double drhodh_pend = (rho_hplus - rho_end)/step;
	//double drhodh_pplus = (rho_hplus_pplus - rho_pplus)/step;
	//double drhodhdp_1 = (drhodp_hplus - drhodp_hend)/step;
	//double drhodhdp_2 = (drhodh_pplus - drhodh_pend)/step;	
	//double drhodhdp_3 = (rho_hplus_pplus + rho_hminus_pminus - rho_hplus_pminus - rho_hminus_pplus)/(4*step*step);
	//double dvdhdp_3 = (1/rho_hplus_pplus + 1/rho_hminus_pminus - 1/rho_hplus_pminus - 1/rho_hminus_pplus)/(4*step*step);

	// Arguments of the spline function (rho is smoothed as a function of h):
	double delta = x * (h_v - h_l);
	double delta_end = h_end - h_l;
	double ddeltaxdp = xend * (dhvdp - dhldp);
	double ddeltadp = -dhldp;

	// Coefficients of the spline function
	double a = 1/pow(delta_end,3) * (2*rho_l - 2*rho_end + delta_end * (drhodh_l + drhodh_end));
	double b = 3/pow(delta_end,2) * (-rho_l + rho_end) - 1/delta_end * (drhodh_end + 2 * drhodh_l);
	double c = drhodh_l;
	double d = rho_l;

	// First pressure derivative of the coefficients
	double dadp = -6/pow(delta_end,4) *  ddeltaxdp * (rho_l - rho_end) + 2/pow(delta_end,3) * (drholdp - drhoxdp) + 1/pow(delta_end,2) * (drhodhdp_l + drhodhdp_end) -2/pow(delta_end,3) * (drhodh_l + drhodh_end)*ddeltaxdp;
	double dbdp = -6/pow(delta_end,3) * ddeltaxdp * (-rho_l + rho_end) + 3/pow(delta_end,2) * (-drholdp + drhoxdp) + 1/pow(delta_end,2) * ddeltaxdp * (drhodh_end + 2 * drhodh_l) - 1/delta_end * (drhodhdp_end + 2 * drhodhdp_l);
	double dcdp = drhodhdp_l;
	double dddp = drholdp;

	// Computing the final useful values:
	*rho_spline = a * pow(delta,3) + b * pow(delta,2) + c * delta + d;
	*dsplinedp = (3*a * pow(delta,2)  +2* b * delta + c)* ddeltadp + pow(delta,3) * dadp + pow(delta,2) * dbdp + delta * dcdp + dddp;
	*dsplinedh = 3 * a * pow(delta,2) + 2*b * delta + c;
	
}



double CoolPropStateClass::drhodp_consth(void){

	if (pFluid->enabled_TTSE_LUT && within_TTSE_range(iP,p(),iH,h()) )
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
			
			double dxdp_h = (dhdpL+_Q*(dhdpV-dhdpL))/(hL-hV);
			double dvdp_h = dvdpL+dxdp_h*(1/rhoV-1/rhoL)+_Q*(dvdpV-dvdpL);

			return -_rho*_rho*dvdp_h;
		}
		else
		{
			_h = h();
			return pFluid->TTSESinglePhase.evaluate_first_derivative(iD,iP,iH,_p,_logp,_h);
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
			double dvdp_h = dvdpL+dxdp_h*(1/rhosatV-1/rhosatL)+_Q*(dvdpV-dvdpL);
			return -_rho*_rho*dvdp_h;
		}
		else
		{
			return 1/(dpdrho_constT()-dpdT_constrho()*dhdrho_constT()/dhdT_constrho());
		}
	}
}

double CoolPropStateClass::d2rhodh2_constp(void){
	if (TwoPhase) { throw ValueError("TwoPhase not supported for d2rhodh2_constp");}
	double A = dpdT_constrho()*dhdrho_constT()-dpdrho_constT()*dhdT_constrho();
	double dAdT_constrho = d2pdT2_constrho()*dhdrho_constT()+dpdT_constrho()*d2hdrhodT()-d2pdrhodT()*dhdT_constrho()-dpdrho_constT()*d2hdT2_constrho();
	double dAdrho_constT = d2pdrhodT()*dhdrho_constT()+dpdT_constrho()*d2hdrho2_constT()-d2pdrho2_constT()*dhdT_constrho()-dpdrho_constT()*d2hdrhodT();
	double ddT_drhodh_p_constrho = 1/A*d2pdT2_constrho()-1/(A*A)*dAdT_constrho*dpdT_constrho();
	double ddrho_drhodh_p_constT = 1/A*d2pdrhodT()-1/(A*A)*dAdrho_constT*dpdT_constrho();
	return ddT_drhodh_p_constrho/dhdT_constp()+ddrho_drhodh_p_constT/dhdrho_constp();
}

double CoolPropStateClass::d2rhodhdp(void){
	if (TwoPhase) { throw ValueError("TwoPhase not supported for d2rhodhdp");}
	double A = dpdT_constrho()*dhdrho_constT()-dpdrho_constT()*dhdT_constrho();
	double dAdT_constrho = d2pdT2_constrho()*dhdrho_constT()+dpdT_constrho()*d2hdrhodT()-d2pdrhodT()*dhdT_constrho()-dpdrho_constT()*d2hdT2_constrho();
	double dAdrho_constT = d2pdrhodT()*dhdrho_constT()+dpdT_constrho()*d2hdrho2_constT()-d2pdrho2_constT()*dhdT_constrho()-dpdrho_constT()*d2hdrhodT();
	double ddT_drhodp_h_constrho = -1/A*d2hdT2_constrho()+1/(A*A)*dAdT_constrho*dhdT_constrho();
	double ddrho_drhodp_h_constT = -1/A*d2hdrhodT()+1/(A*A)*dAdrho_constT*dhdT_constrho();
	return ddT_drhodp_h_constrho/dhdT_constp()+ddrho_drhodp_h_constT/dhdrho_constp();
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
	return pFluid->R()*(-tau*tau*d3phir_dDelta_dTau2(tau,delta)+delta*d2phir_dDelta2(tau,delta)+dphir_dDelta(tau,delta)-delta*tau*d3phir_dDelta2_dTau(tau,delta)-tau*d2phir_dDelta_dTau(tau,delta))/pFluid->reduce.rho;
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
double CoolPropStateClass::d2rhodhdQ(void)
{
	return 2/_rho*pow(drhodh_constp(),2)*(hV() - hL());
}
double CoolPropStateClass::d2rhodpdQ(void)
{
	double d2vdhdp = 1/_T*d2Tdp2_along_sat() - pow(dTdp_along_sat()/_T,2);
	return (2/_rho*drhodp_consth()*drhodh_constp()-pow(_rho,2)*d2vdhdp)*(hV() - hL());
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

/// Enable the TTSE
void CoolPropStateClass::enable_TTSE_LUT(void){pFluid->enable_TTSE_LUT();};
/// Check if TTSE is enabled
bool CoolPropStateClass::isenabled_TTSE_LUT(void){return pFluid->isenabled_TTSE_LUT();};
/// Disable the TTSE
void CoolPropStateClass::disable_TTSE_LUT(void){pFluid->disable_TTSE_LUT();};
/// Enable the writing of TTSE tables to file
void CoolPropStateClass::enable_TTSE_LUT_writing(void){pFluid->enable_TTSE_LUT_writing();};
/// Check if the writing of TTSE tables to file is enabled
bool CoolPropStateClass::isenabled_TTSE_LUT_writing(void){return pFluid->isenabled_TTSE_LUT_writing();};
/// Disable the writing of TTSE tables to file
void CoolPropStateClass::disable_TTSE_LUT_writing(void){pFluid->disable_TTSE_LUT_writing();};
/// Over-ride the default size of both of the saturation LUT
void CoolPropStateClass::set_TTSESat_LUT_size(int N){pFluid->set_TTSESat_LUT_size(N);};
/// Over-ride the default size of the single-phase LUT
void CoolPropStateClass::set_TTSESinglePhase_LUT_size(int Np, int Nh){pFluid->set_TTSESinglePhase_LUT_size(Np,Nh);};
/// Over-ride the default range of the single-phase LUT
void CoolPropStateClass::set_TTSESinglePhase_LUT_range(double hmin, double hmax, double pmin, double pmax){pFluid->set_TTSESinglePhase_LUT_range(hmin,hmax,pmin,pmax);};
/// Get the current range of the single-phase LUT
void CoolPropStateClass::get_TTSESinglePhase_LUT_range(double *hmin, double *hmax, double *pmin, double *pmax){pFluid->get_TTSESinglePhase_LUT_range(hmin,hmax,pmin,pmax);};
