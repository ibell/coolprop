#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <exception>
#include <iostream>
#include "Helmholtz.h"
#include "FluidClass.h"
#include "CoolProp.h"
#include "CoolPropTools.h"
#include "PengRobinson.h"
#include "Solvers.h"
#include "R134a.h"
#include <stdlib.h>
#include <stdio.h>

using namespace std;

// Destructor needs to free all dynamically allocated objects
// In this case, the phir and phi0 components
Fluid::~Fluid()
{
	while (!phirlist.empty())
	{
		delete phirlist.back();  
		phirlist.pop_back();
	}
	while (!phi0list.empty())
	{
		delete phi0list.back();  
		phi0list.pop_back();
	}
	EOSReference.clear();
	TransportReference.clear();
}
void Fluid::post_load(void)
{
	// Set the reducing values from the pointer
	reduce=*preduce;
	// Set the triple-point pressure if not set in code
	// Use the dewpoint pressure for pseudo-pure fluids
	if ( fabs(params.ptriple)>1e9 ){
		double pL, pV, rhoL, rhoV;
		saturation(params.Ttriple,false,&pL,&pV,&rhoL,&rhoV);
		params.ptriple = pV;
	}
}
//--------------------------------------------
//    Residual Part
//--------------------------------------------

double Fluid::phir(double tau, double delta)
{
	double summer = 0;
	for (list<phi_BC*>::iterator it = phirlist.begin(); it != phirlist.end(); it++)
		summer += (*it)->base(tau,delta);
	return summer;
}
double Fluid::dphir_dDelta(double tau, double delta)
{
	double summer = 0;
	for (list<phi_BC*>::iterator it = phirlist.begin(); it != phirlist.end(); it++)
		summer += (*it)->dDelta(tau,delta);
	return summer;
}
double Fluid::d2phir_dDelta2(double tau, double delta)
{
	double summer = 0;
	for (list<phi_BC*>::iterator it = phirlist.begin(); it != phirlist.end(); it++)
		summer += (*it)->dDelta2(tau,delta);
	return summer;
}
double Fluid::dphir_dTau(double tau, double delta)
{
	double summer = 0;
	for (list<phi_BC*>::iterator it = phirlist.begin(); it != phirlist.end(); it++)
		summer += (*it)->dTau(tau,delta);
	return summer;
}
double Fluid::d2phir_dTau2(double tau, double delta)
{
	double summer = 0;
	for (list<phi_BC*>::iterator it = phirlist.begin(); it != phirlist.end(); it++)
		summer += (*it)->dTau2(tau,delta);
	return summer;
}
double Fluid::d2phir_dDelta_dTau(double tau, double delta)
{
	double summer = 0;
	for (list<phi_BC*>::iterator it = phirlist.begin(); it != phirlist.end(); it++)
		summer += (*it)->dDelta_dTau(tau,delta);
	return summer;
}
double Fluid::d3phir_dDelta2_dTau(double tau, double delta)
{
	double summer = 0;
	for (list<phi_BC*>::iterator it = phirlist.begin(); it != phirlist.end(); it++)
		summer += (*it)->dDelta2_dTau(tau,delta);
	return summer;
}

//--------------------------------------------
//    Ideal Gas Part
//--------------------------------------------
double Fluid::phi0(double tau, double delta)
{
	double summer = 0;
	for (list<phi_BC*>::iterator it = phi0list.begin(); it != phi0list.end(); it++){
		double contrib=(*it)->base(tau,delta);
		summer += (*it)->base(tau,delta);
	}	
	return summer;
}
double Fluid::dphi0_dDelta(double tau, double delta)
{
	double summer = 0;
	for (list<phi_BC*>::iterator it = phi0list.begin(); it != phi0list.end(); it++)
		summer += (*it)->dDelta(tau,delta);
	return summer;
}
double Fluid::d2phi0_dDelta2(double tau, double delta)
{
	double summer = 0;
	for (list<phi_BC*>::iterator it = phi0list.begin(); it != phi0list.end(); it++)
		summer += (*it)->dDelta2(tau,delta);
	return summer;
}
double Fluid::dphi0_dTau(double tau, double delta)
{
	double summer = 0;
	for (list<phi_BC*>::iterator it = phi0list.begin(); it != phi0list.end(); it++){
		double contribution  =	(*it)->dTau(tau,delta);
		summer += (*it)->dTau(tau,delta);
	}
	return summer;
}
double Fluid::d2phi0_dTau2(double tau, double delta)
{
	double summer = 0;
	for (list<phi_BC*>::iterator it = phi0list.begin(); it != phi0list.end(); it++){
		double contrib = (*it)->dTau2(tau,delta);
		summer += (*it)->dTau2(tau,delta);
	}
	return summer;
}
double Fluid::d2phi0_dDelta_dTau(double tau, double delta)
{
	double summer = 0;
	for (list<phi_BC*>::iterator it = phi0list.begin(); it != phi0list.end(); it++)
		summer += (*it)->dDelta_dTau(tau,delta);
	return summer;
}

bool Fluid::isAlias(std::string name)
{
	// Returns true if name is an alias of the fluid
	for (list<std::string>::iterator it = aliases.begin(); it != aliases.end(); it++)
		if (name.compare((*it))==0)
		{
			return true;
		}
	return false;
}
// ----------------------------------------
//             Properties
// ----------------------------------------

double Fluid::pressure_Trho(double T, double rho)
{
    double delta,tau,R;
	R=params.R_u/params.molemass;
	tau=reduce.T/T;
	delta=rho/reduce.rho;
	return R*T*rho*(1.0+delta*dphir_dDelta(tau,delta));
}

double Fluid::enthalpy_Trho(double T, double rho)
{
	double delta,tau,R;
    R=params.R_u/params.molemass;
	tau=reduce.T/T;
	delta=rho/reduce.rho;
    return R*T*(1.0+tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))+delta*dphir_dDelta(tau,delta));
}

double Fluid::internal_energy_Trho(double T, double rho)
{
    double delta,tau,R;
    R=params.R_u/params.molemass;
	tau=reduce.T/T;
	delta=rho/reduce.rho;
    return R*T*tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta));
}

double Fluid::entropy_Trho(double T, double rho)
{
    double delta,tau,R;
    R=params.R_u/params.molemass;
	tau=reduce.T/T;
	delta=rho/reduce.rho;
    return R*(tau*(dphi0_dTau(tau,delta)+dphir_dTau(tau,delta))-phi0(tau,delta)-phir(tau,delta));
}
double Fluid::specific_heat_v_Trho(double T, double rho)
{
    double delta,tau,R;
    R=params.R_u/params.molemass;
	tau=reduce.T/T;
	delta=rho/reduce.rho;
    return -R*pow(tau,2)*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta));
}
double Fluid::specific_heat_p_Trho(double T, double rho)
{
    double delta,tau,c1,c2,R,dphir_dDelta_;
    R=params.R_u/params.molemass;
	tau=reduce.T/T;
	delta=rho/reduce.rho;
	dphir_dDelta_=dphir_dDelta(tau,delta);
    c1=pow(1.0+delta*dphir_dDelta_-delta*tau*d2phir_dDelta_dTau(tau,delta),2);
    c2=(1.0+2.0*delta*dphir_dDelta_+pow(delta,2)*d2phir_dDelta2(tau,delta));
    return R*(-pow(tau,2)*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta))+c1/c2);
}
double Fluid::specific_heat_p_ideal_Trho(double T)
{
    double tau,R;
    R=params.R_u/params.molemass;
	tau=reduce.T/T;
    return R*(1-tau*tau*d2phi0_dTau2(tau, 1e-12));
}
double Fluid::speed_sound_Trho(double T, double rho)
{
    double delta,tau,c1,c2,R;
    R=params.R_u/params.molemass;
	tau=reduce.T/T;
	delta=rho/reduce.rho;

    c1=-specific_heat_v_Trho(T,rho)/R;
    c2=(1.0+2.0*delta*dphir_dDelta(tau,delta)+pow(delta,2)*d2phir_dDelta2(tau,delta));
    return sqrt(-c2*T*specific_heat_p_Trho(T,rho)*1000/c1);
}
double Fluid::gibbs_Trho(double T,double rho)
{
	double delta,tau,R;
    R=params.R_u/params.molemass;
	tau=reduce.T/T;
	delta=rho/reduce.rho;

    return R*T*(1+phi0(tau,delta)+phir(tau,delta)+delta*dphir_dDelta(tau,delta));
}
double Fluid::dpdT_Trho(double T,double rho)
{
	double delta,tau,R;
    R=params.R_u/params.molemass;
	tau=reduce.T/T;
	delta=rho/reduce.rho;

	return rho*R*(1+delta*dphir_dDelta(tau,delta)-delta*tau*d2phir_dDelta_dTau(tau,delta));
}

double Fluid::density_Tp(double T, double p)
{
	// If the guess value is not provided, calculate the guess density using Peng-Robinson
	// This overload is used to pre-calculate the guess for density using PR if possible
	if (debug()>8){
		std::cout<<__FILE__<<':'<<__LINE__<<": Fluid::density_Tp(double T, double p): "<<T<<","<<p<<std::endl;
	}
	return density_Tp(T,p,_get_rho_guess(T,p));
}

double Fluid::density_Tp(double T, double p, double rho_guess)
{
    double delta,tau,dpdrho,error=999,R,p_EOS,rho,change=999;
    R=params.R_u/params.molemass;
	tau=reduce.T/T;

	if (debug()>8){
		std::cout<<__FILE__<<':'<<__LINE__<<": Fluid::density_Tp(double T, double p, double rho_guess): "<<T<<","<<p<<","<<rho_guess<<std::endl;
	}
	
	// Start with the guess value
	rho=rho_guess;
    int iter=1;
    while (fabs(error)>1e-10 && change>1e-11)
    {
        delta=rho/reduce.rho;
		// Run and save to cut down on calculations
		double dphirdDelta = dphir_dDelta(tau,delta);
        // Use Newton's method to find the density since the derivative of pressure w.r.t. density is known from EOS
        dpdrho=R*T*(1+2*delta*dphirdDelta+delta*delta*d2phir_dDelta2(tau,delta));
		// Pressure from equation of state
		p_EOS = rho*R*T*(1+delta*dphirdDelta);
        // Update the step using Newton's method
        rho=rho-(p_EOS-p)/dpdrho;
		change = fabs((p_EOS-p)/dpdrho);
        // Residual
        error=p_EOS-p;
        iter++;
        if (iter>200)
        {
			throw SolutionError(format("Number of steps in density_TP has exceeded 100 with inputs T=%g,p=%g,rho_guess=%g for fluid %s\n",T,p,rho_guess,name.c_str()));
            return _HUGE;
        }
    }	
	//std::cout << iter << " " << rho << " " << rho_guess << std::endl;
    return rho;
}

double Fluid::viscosity_Trho( double T, double rho)
{
	// Use propane as the reference
	Fluid * ReferenceFluid = new R134aClass();
	ReferenceFluid->post_load();
	// Calculate the ECS
	double mu = viscosity_ECS_Trho(T, rho, ReferenceFluid);
	// Delete the reference fluid instance
	delete ReferenceFluid;
	return mu;
}
double Fluid::conductivity_Trho( double T, double rho)
{
	// Use propane as the reference
	Fluid * ReferenceFluid = new R134aClass();
	ReferenceFluid->post_load();
	// Calculate the ECS
	double lambda = conductivity_ECS_Trho(T, rho, ReferenceFluid);
	// Delete the reference fluid instance
	delete ReferenceFluid;
	return lambda;
}

void Fluid::set_1phase_LUT_params(int nT,int np,double Tmin, double Tmax, double pmin, double pmax)
{
	if (debug()>0){
		std::cout << __FILE__<<": Setting 1phase_LUT_params of ("<<name<<","<<nT<<","<<np<<","<<Tmin<<","<<Tmax<<","<<pmin<<","<<pmax<<")"<<std::endl;
	}
	LUT.nT=nT; 
	LUT.np=np; 
	LUT.Tmin=Tmin;
	LUT.Tmax=Tmax;
	LUT.pmin=pmin;
	LUT.pmax=pmax;
	LUT.built=false;  // Rebuild on the next call
	return;
}
void Fluid::get_1phase_LUT_params(int *nT,int *np,double *Tmin, double *Tmax, double *pmin, double *pmax)
{
	*nT=LUT.nT; 
	*np=LUT.np; 
	*Tmin=LUT.Tmin;
	*Tmax=LUT.Tmax;
	*pmin=LUT.pmin;
	*pmax=LUT.pmax;
	if (debug()>0){
		std::cout << __FILE__<<": Getting 1phase_LUT_params of ("<<","<<name<<","<<*nT<<","<<*np<<","<<*Tmin<<","<<*Tmax<<","<<*pmin<<","<<*pmax<<")"<<std::endl;
	}
	return;
}
void Fluid::saturation(double T, bool UseLUT, double *psatLout, double *psatVout, double *rhosatLout, double *rhosatVout)
{
	double rhoL, rhoV, p;
	if (isPure==true)
	{
		if (UseLUT)
		{
			// Use the saturation Lookup Table;
			*rhosatLout=ApplySaturationLUT("rhoL","T",T);
			*rhosatVout=ApplySaturationLUT("rhoV","T",T);
			*psatLout=ApplySaturationLUT("p","T",T);
			*psatVout=ApplySaturationLUT("p","T",T);
			return;
		}
		else
		{
			rhosatPure_Akasaka(T,&rhoL,&rhoV,&p);
			*rhosatLout = rhoL;
			*rhosatVout = rhoV;
			*psatLout = p;
			*psatVout = p;
			return;
		}
	}
	else
	{ 
		// Pseudo-pure fluid
		*rhosatLout = rhosatL(T);
		*rhosatVout = rhosatV(T);
		*psatLout = psatL(T);
		*psatVout = psatV(T);
		return;
	}
}

void Fluid::rhosatPure_Akasaka(double T, double *rhoLout, double *rhoVout, double *pout)
{
	double rhoL,rhoV,dphirL,dphirV,ddphirL,ddphirV,phirL,phirV,JL,JV,KL,KV,dJL,dJV,dKL,dKV;
	double DELTA, deltaL, deltaV, tau, gamma = 1.0, error, PL, PV;
	int iter=0;
	// Use the density ancillary function as the starting point for the solver
    try
	{
		rhoL=rhosatL(T);
		rhoV=rhosatV(T);
		deltaL = rhoL/reduce.rho;
		deltaV = rhoV/reduce.rho;
		tau = reduce.T/T;
	}
	catch(NotImplementedError &e)
	{
		std::cout << e.what() << "Ancillary not provided" << std::endl;
	}

	do{
		// Calculate once to save on calls to EOS
		dphirL = dphir_dDelta(tau,deltaL);
		dphirV = dphir_dDelta(tau,deltaV);
		ddphirL = d2phir_dDelta2(tau,deltaL);
		ddphirV = d2phir_dDelta2(tau,deltaV);
		phirL = phir(tau,deltaL);
		phirV = phir(tau,deltaV);

		JL = deltaL * (1 + deltaL*dphirL);
		JV = deltaV * (1 + deltaV*dphirV);
		KL = deltaL*dphirL + phirL + log(deltaL);
		KV = deltaV*dphirV + phirV + log(deltaV);
		dJL = 1 + 2*deltaL*dphirL + deltaL*deltaL*ddphirL;
		dJV = 1 + 2*deltaV*dphirV + deltaV*deltaV*ddphirV;
		dKL = 2*dphirL + deltaL*ddphirL + 1/deltaL;
		dKV = 2*dphirV + deltaV*ddphirV + 1/deltaV;
		
		DELTA = dJV*dKL-dJL*dKV;

		error = fabs(KL-KV)+fabs(JL-JV);

		//Update the step
		deltaL += gamma/DELTA*( (KV-KL)*dJV-(JV-JL)*dKV);
		deltaV += gamma/DELTA*( (KV-KL)*dJL-(JV-JL)*dKL);
		
		iter=iter+1;
	}
	while (error > 1e-10);
	*rhoLout = deltaL*reduce.rho;
	*rhoVout = deltaV*reduce.rho;
	PL = pressure_Trho(T,*rhoLout);
	PV = pressure_Trho(T,*rhoVout);
	*pout = 0.5*(PL+PV);
	return;
}
void Fluid::rhosatPure(double T, double *rhoLout, double *rhoVout, double *pout)
{
    // Only works for pure fluids (no blends)
    // At equilibrium, saturated vapor and saturated liquid are at the same pressure and the same Gibbs energy
    double rhoL,rhoV,p,error=999,x1,x2,x3,y1,y2,f,p_guess;
    int iter;

    if (T>=crit.T || T<(params.Ttriple-0.0001))
    {
		throw ValueError(format("Temperature of fluid [%g] is out of range from %g K to %g K",T,params.Ttriple,crit.T));
    }

    // Use the density ancillary function as the starting point for the secant solver
    try
	{
		rhoL=rhosatL(T);
		rhoV=rhosatV(T);
	}
	catch(NotImplementedError &e)
	{
		std::cout << e.what() << "Ancillary not provided" << std::endl;
	}
    p_guess=pressure_Trho(T,rhoV);

    iter=1;
    // Use a secant method to obtain pressure
    while ((iter<=3 || fabs(error)>1e-10) && iter<100)
    {
        if (iter==1){x1=p_guess; p=x1;}
        if (iter==2){x2=1.0001*p_guess; p=x2;}
        if (iter>2) {p=x2;}
            //Recalculate the densities based on the current pressure
            rhoL=density_Tp(T,p,rhoL);
            rhoV=density_Tp(T,p,rhoV);
            // Residual between saturated liquid and saturated vapor gibbs function
            f=gibbs_Trho(T,rhoL)-gibbs_Trho(T,rhoV);
        if (iter==1){y1=f;}
        if (iter>1)
        {
            y2=f;
            x3=x2-y2/(y2-y1)*(x2-x1);
            error=f;
            y1=y2; x1=x2; x2=x3;
        }
        iter++;
        if (iter>100)
        {
        	//ERROR
            printf("rhosatPure failed, current values of rhoL and rhoV are %g,%g\n",rhoL,rhoV);
            *rhoLout=_HUGE;
            *rhoVout=_HUGE;
            *pout=_HUGE;
            return;
        }
    }
    *rhoLout=rhoL;
    *rhoVout=rhoV;
    *pout=p;
    return;
}
bool isBetween(double x1,double x2, double x)
{
	if (x2>x1 && x>=x1 && x<=x2) return true;
	else if (x1>x2 && x<=x1 && x>=x2) return true;
	else return false;
}
void Fluid::BuildSaturationLUT()
{
    // returns the index of the LUT
    int i;
    double T_t,T_c,dT;
    
	if (SatLUT.built==true)
		return;
    // Linearly space the values from the triple point of the fluid to the just shy of the critical point
    T_t=Props(name,"Ttriple")+0.00001;
    T_c=Props(name,"Tcrit")-0.000001;
    
	SatLUT.T=std::vector<double>(SatLUT.N,0);
	SatLUT.tau=std::vector<double>(SatLUT.N,0);
	SatLUT.p=std::vector<double>(SatLUT.N,0);
	SatLUT.logp=std::vector<double>(SatLUT.N,0);
	SatLUT.rhoL=std::vector<double>(SatLUT.N,0);
	SatLUT.rhoV=std::vector<double>(SatLUT.N,0);
	SatLUT.hL=std::vector<double>(SatLUT.N,0);
	SatLUT.hV=std::vector<double>(SatLUT.N,0);

    dT=(T_c-T_t)/(SatLUT.N-1);
    for (i=0;i<SatLUT.N;i++)
    {
        // Linearly space the elements
        SatLUT.T[i]=T_t+i*dT;
		// Build the tau vector - for better interpolation behavior
		SatLUT.tau[i]=reduce.T/SatLUT.T[i];
        // Calculate the saturation properties
        rhosatPure_Akasaka(SatLUT.T[i], &(SatLUT.rhoL[i]), &(SatLUT.rhoV[i]), &(SatLUT.p[i]));
		SatLUT.logp[i]=log(SatLUT.p[i]);
        // Calculate the other saturation properties
		SatLUT.hL[i]=Props(std::string("H"),'T',SatLUT.T[i],'D',SatLUT.rhoL[i],name);
        SatLUT.hV[i]=Props(std::string("H"),'T',SatLUT.T[i],'D',SatLUT.rhoV[i],name);
    }
	SatLUT.built=true;
}   

double Fluid::ApplySaturationLUT(char *OutPropName,char *InPropName,double InPropVal)
{
    int i;
    double x0,x1,x2,y0,y1,y2,LUTVal;

    // pointers to the x and y std::vectors for later
	std::vector<double> *y,*x;
    
    // First try to build the LUT
    BuildSaturationLUT();
   
	if (!strcmp(InPropName,"T") || !strcmp(InPropName,"Tsat"))
	{
		x=&(SatLUT.tau);
		LUTVal = reduce.T/InPropVal;
	}
	else if (!strcmp(InPropName,"P") || !strcmp(InPropName,"psat"))
	{
		x=&(SatLUT.logp);
		LUTVal = log(InPropVal);
	}
	else 
	{
		throw AttributeError(format("InPropName [%s] to ApplySaturationLUT is invalid.  Valid values are T,p,Tsat,psat",InPropName));
	}

    // Then find indices that bracket the value of interest
    for(i=0;i<SatLUT.N-1;i++)
    {
        if (isBetween((*x)[i],(*x)[i+1],LUTVal)) break;
    }
    
    if (!strcmp(OutPropName,"rhoL"))
        y=&(SatLUT.rhoL);
    else if (!strcmp(OutPropName,"rhoV"))
        y=&(SatLUT.rhoV);
    else if (!strcmp(OutPropName,"p"))
        y=&(SatLUT.logp);
    else if (!strcmp(OutPropName,"hL"))
        y=&(SatLUT.hL);
    else if (!strcmp(OutPropName,"hV"))
        y=&(SatLUT.hV);
	else if (!strcmp(OutPropName,"T"))
        y=&(SatLUT.tau);

    // Need a three-point set to interpolate using a quadratic.
    if (i<SatLUT.N-2)
    {
		// Go "forwards" with the interpolation range
        x0=(*x)[i];
        x1=(*x)[i+1];
        x2=(*x)[i+2];
        y0=(*y)[i];
        y1=(*y)[i+1];
        y2=(*y)[i+2];
    }
    else
    {
        // Go "backwards" with the interpolation range
        x0=(*x)[i];
        x1=(*x)[i-1];
        x2=(*x)[i-2];
        y0=(*y)[i];
        y1=(*y)[i-1];
        y2=(*y)[i-2];
    }
    double y_LUT = QuadInterp(x0,x1,x2,y0,y1,y2,LUTVal);
	if (!strcmp(OutPropName,"p"))
		// y_LUT has the value of log(p)
		return exp(y_LUT);
	else if (!strcmp(OutPropName,"p"))
		// yLUT has the value of Tc/T
		return reduce.T/y_LUT;
	else
		return y_LUT;
}

double Fluid::Tsat(double p, double Q, double T_guess)
{
	return Tsat(p,Q,T_guess,false);
}

double Fluid::_get_rho_guess(double T, double p)
{
    double eps=1e-10,Tc,rho_guess,rho_simple;
    int counter=1;
    double pEOS;

	Tc = reduce.T;

	// Then based on phase, select the useful solution(s)
	std::string phase = Fluid::phase_Tp(T,p);

	// These are very simplistic guesses for the density, but they work ok, use them if PR fails
	if (!strcmp(phase.c_str(),"Gas") || !strcmp(phase.c_str(),"Supercritical")){
		// Guess that it is ideal gas
		rho_simple = p/(R()*T);
	}
	else if (!strcmp(phase.c_str(),"Liquid")){
		// Start at the saturation state, with a known density
		// Return subcooled liquid density
		double rhoL,rhoV,pL,pV;
		saturation(T,false,&pL,&pV,&rhoL,&rhoV);
		double delta = rhoL / reduce.rho;
		double tau = reduce.T/T;
		double dp_drho = R()*T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta));
		double drho_dp = 1/dp_drho;
		rho_simple = rhosatL(T)-drho_dp*(pL-p);
	}
	return rho_simple;

	// Try to use Peng-Robinson to get a guess value for the density
	std::vector<double> rholist= PRGuess_rho(this,T,p);
	if (rholist.size()==0){
		throw SolutionError(format("PengRobinson could not yield any solutions"));
	}
	else if (rholist.size()==1){
		rho_guess = rholist[0];
	}
	else{
		// A list of ok solutions
		std::vector<double> goodrholist;
		for (unsigned int i=0;i<rholist.size();i++){
			pEOS=Props(std::string("P"),'T',T,'D',rholist[i],name.c_str());
			if (fabs(pEOS/p-1)<0.03){
				goodrholist.push_back(rholist[i]);
			}
		}
		if (goodrholist.size()==1){
			rho_guess = goodrholist[0];
		}
		else{
			rho_guess = rho_simple;
		}
	}
	return rho_guess;
}

static void swap(double *x, double *y)
{
    double temp;
    temp = *x;
    *x = *y;
    *y = temp;
}
std::string Fluid::phase_Tp(double T, double p)
{
	double pL,pV,rhoL,rhoV;
	/*
		|         |     
		|         |    Supercritical
		|         |
	p	| Liquid (b)------------
		|        /
		|       / 
		|      /       Gas
		|     / 
		|   (a)
		|  -
		|------------------------

		           T

	   a: triple point
	   b: critical point
	   a-b: Saturation line

	*/
	if (T>crit.T && p>crit.p){
		return std::string("Supercritical");
	}
	else if (T>crit.T && p<crit.p){
		return std::string("Gas");
	}
	else if (T<crit.T && p>crit.p){
		return std::string("Liquid");
	}
	else if (p<params.ptriple){
		return std::string("Gas");
	}
	else{
		// Now start to think about the saturation stuff
		// First try to use the ancillary equations if you are far enough away
		// Ancillary equations are good to within 1% in pressure in general
		// Some industrial fluids might not be within 3%
		if (isPure && p<0.94*psat(T)){
			return std::string("Gas");
		}
		else if (isPure && p>1.06*psat(T)){
			return std::string("Liquid");
		}
		else if (!isPure && p<0.94*psatV(T)){
			return std::string("Gas");
		}
		else if (!isPure && p>1.06*psatL(T)){
			return std::string("Liquid");
		}
		else{
			// Actually have to use saturation information sadly
			// For the given temperature, find the saturation state
			// Run the saturation routines to determine the saturation densities and pressures
			saturation(T,SaturationLUTStatus(),&pL,&pV,&rhoL,&rhoV);
			if (p>pL){
				return std::string("Liquid");
			}
			else if (p<pV){
				return std::string("Gas");
			}
			else{
				return std::string("Two-Phase");
			}
		}
	}
}

double Fluid::Tsat(double p, double Q, double T_guess, bool UseLUT)
{
    double x1=0,x2=0,y1=0,y2=0;
    double Tc,Tmax,Tmin;
    
    // Do reverse interpolation in the Saturation Lookup Table
    if (UseLUT && isPure==true) { return ApplySaturationLUT("T","p",p); }

    Tc=Props(name,"Tcrit");
    Tmax=Tc-0.1;
    Tmin=Props(name,"Ttriple")+1;
    
    // Plotting Tc/T versus log(p) tends to give very close to straight line
    // Use this fact find T more easily
    
	class SatFuncClass : public FuncWrapper1D
	{
	private:
		double p,Q,Tc;
		std::string name;
	public:
		SatFuncClass(double p_, double Q_, double Tc_, std::string name_){
			p=p_;Q=Q_;Tc=Tc_,name=name_;
		};
		double call(double tau){
			return log(Props(std::string("P"),'T',Tc/tau,'Q',Q,name)/p);
		};
	} SatFunc(p,Q,reduce.T,name);

	double tau_max = Tc/Tmin;
	double tau_min = Tc/Tmax;

	// Use Brent's method to find tau = Tc/T
	std::string errstr;
    double tau = Brent(&SatFunc,tau_min,tau_max,1e-10,1e-10,50,&errstr);
	if (errstr.size()>0)
		throw SolutionError("Saturation calculation failed");
	return reduce.T/tau;

}

void Fluid::BuildLookupTable()
{
    int i,j;
	bool OldUseLUT,okval;
    double T,p,rho,Tc;
    
	// Build the Lookup matrices, but only if they ar not already built
	if (LUT.built==true)
		return;
	
	if (debug()>5){
		std::cout<<__FILE__<<": Building 1-phase lookup for "<<name<<std::endl;
		std::cout << __FILE__<<": 1phase_LUT_params of ("<<","<<name<<","<<LUT.nT
			<<","<<LUT.np<<","<<LUT.Tmin<<","<<LUT.Tmax<<","<<LUT.pmin<<","<<LUT.pmax<<")"<<std::endl;
	}
    Tc=reduce.T;

	LUT.Tvec=std::vector<double> ( LUT.nT,0 );
	LUT.pvec=std::vector<double> ( LUT.np,0 );
	LUT.kmat=std::vector< std::vector<double> >( LUT.nT, std::vector<double> ( LUT.np, 0 ) );
	LUT.kmat=std::vector< std::vector<double> >( LUT.nT, std::vector<double> ( LUT.np, 0 ) );
	LUT.cpmat=std::vector< std::vector<double> >( LUT.nT, std::vector<double> ( LUT.np, 0 ) );
	LUT.cp0mat=std::vector< std::vector<double> >( LUT.nT, std::vector<double> ( LUT.np, 0 ) );
	LUT.cvmat=std::vector< std::vector<double> >( LUT.nT, std::vector<double> ( LUT.np, 0 ) );
	LUT.pmat=std::vector< std::vector<double> >( LUT.nT, std::vector<double> ( LUT.np, 0 ) );
	LUT.rhomat=std::vector< std::vector<double> >( LUT.nT, std::vector<double> ( LUT.np, 0 ) );
	LUT.smat=std::vector< std::vector<double> >( LUT.nT, std::vector<double> ( LUT.np, 0 ) );
	LUT.umat=std::vector< std::vector<double> >( LUT.nT, std::vector<double> ( LUT.np, 0 ) );
	LUT.viscmat=std::vector< std::vector<double> >( LUT.nT, std::vector<double> ( LUT.np, 0 ) );
	LUT.hmat=std::vector< std::vector<double> >( LUT.nT, std::vector<double> ( LUT.np, 0 ) );
	LUT.rhomat=std::vector< std::vector<double> >( LUT.nT, std::vector<double> ( LUT.np, 0 ) );
	LUT.hmat=std::vector< std::vector<double> >( LUT.nT, std::vector<double> ( LUT.np, 0 ) );
	LUT.dpdTmat=std::vector< std::vector<double> >( LUT.nT, std::vector<double> ( LUT.np, 0 ) );
    
    // Properties evaluated at all points with X in the 
    // following p-h plot:
    
    /*      Supercritical
    ||X X X X X X X X X X X X X X X X X X X X X X
    ||X X X X X X X X X X X X X X X X X X X X X X
    ||X X X -------  X X X X X X X X X X X X
    ||X X  /	   \  X X X X X X X X X X X
p	||X  /		    | X X X X Superheated Gas
    ||X /	Two     / X X X X X X X X X X X
    || /  Phase	   /X X X X X X X X X X X X 
    ||/		      / X X X X X X X X X X X X 
    ||=	= = = = = = = = = = = = = = = = = = = = =
                        Enthalpy										
    */
    
	// Get the old value for the LUT flag
	OldUseLUT=SinglePhaseLUTStatus();
	
	// Turn off the LUT so that it will use EOS to determine state
    UseSinglePhaseLUT(false);

    for (i=0;i<LUT.nT;i++)
    {
        LUT.Tvec[i]=LUT.Tmin+i*(LUT.Tmax-LUT.Tmin)/(LUT.nT-1);
    }
    for (j=0;j<LUT.np;j++)
    {
        LUT.pvec[j]=LUT.pmin+j*(LUT.pmax-LUT.pmin)/(LUT.np-1);
    }
    for (i=0;i<LUT.nT;i++)
    {
        for (j=0;j<LUT.np;j++)
        {
			T = LUT.Tvec[i];
			p = LUT.pvec[j];

			/*
			For a given fluid, if the temperature is above the critical temperature, 
			add an entry to the LUT matrix.  If temperature is below critical temp, need to
			consider a few other options.

			If it is a pure fluid, the saturated vapor and saturated liquid lines
			collapse onto each other, resulting in only one pressure that yields an undefined value

			If is is a pseudo-pure fluid, the bubble point (saturated liquid) pressure will ALWAYS 
			be above that of the dew point (saturated vapor) pressure, so if the pressure is not between
			the dewpoint and bubblepoint presssures, it is single phase
			*/
			okval=false;
			if (T>Tc){
				okval=true;
			}
			else if (isPure && (p>psat(T) || p<psat(T))){
				okval=true;
			}
			else if (!isPure && (p>psatL(T) || p<psatV(T))){
				okval=true;
			}

            if (okval)
            {					
                LUT.rhomat[i][j]=density_Tp(T,p);
				rho=LUT.rhomat[i][j];
                LUT.hmat[i][j]=enthalpy_Trho(T,rho);
                LUT.smat[i][j]=entropy_Trho(T,rho);
                LUT.umat[i][j]=internal_energy_Trho(T,rho);
                LUT.cpmat[i][j]=specific_heat_p_Trho(T,rho);
				LUT.cp0mat[i][j]=specific_heat_p_ideal_Trho(T);
                LUT.cvmat[i][j]=specific_heat_v_Trho(T,rho);
                LUT.viscmat[i][j]=viscosity_Trho(T,rho);
                LUT.kmat[i][j]=conductivity_Trho(T,rho);
				LUT.dpdTmat[i][j]=dpdT_Trho(T,rho);
                LUT.pmat[i][j]=p;
            }
            else
            {
                LUT.rhomat[i][j]=_HUGE;
                LUT.hmat[i][j]=_HUGE;
                LUT.smat[i][j]=_HUGE;
                LUT.umat[i][j]=_HUGE;
                LUT.cpmat[i][j]=_HUGE;
				LUT.cp0mat[i][j]=_HUGE;
                LUT.cvmat[i][j]=_HUGE;
                LUT.viscmat[i][j]=_HUGE;
                LUT.kmat[i][j]=_HUGE;
                LUT.pmat[i][j]=_HUGE;
            }
        }
    }
    UseSinglePhaseLUT(OldUseLUT);
	LUT.built=true;
    return;
}

std::vector<std::vector<double> >* Fluid::_get_LUT_ptr(std::string Prop)
{
	std::vector<std::vector<double> > *mat;
	/* Depending on which property is desired, 
    make the matrix "mat" a pointer to the 
    desired property matrix */
	if (Prop[0] == 'C')
        mat=&LUT.cpmat;
	else if (Prop[0] == 'C' && Prop[1]=='0')
        mat=&LUT.cp0mat;
    else if (Prop[0] == 'P')
        mat=&LUT.pmat;
	else if (Prop[0] == 'D')
        mat=&LUT.rhomat;
    else if (Prop[0] == 'O')
        mat=&LUT.cvmat;
    else if (Prop[0] == 'H')
        mat=&LUT.hmat;
    else if (Prop[0] == 'S')
        mat=&LUT.smat;
    else if (Prop[0] == 'U')
        mat=&LUT.umat;
    else if (Prop[0] == 'V')
        mat=&LUT.viscmat;
    else if (Prop[0] == 'L')
        mat=&LUT.kmat;
	else if (!Prop.compare("dpdT"))
		mat=&LUT.dpdTmat;
    else
    {
    	throw ValueError(format("Invalid output type [%s] in _get_LUT_ptr\n",Prop.c_str()));
    }
	return mat;
}
double Fluid::LookupValue_TP(std::string Prop, double T, double p)
{
    int iPlow, iPhigh, iTlow, iThigh;
    double T1, T2, T3, P1, P2, P3, y1, y2, y3, a1, a2, a3;
	std::vector<std::vector<double> > *mat;

    // Build if needed
    BuildLookupTable();
    
    if (T>LUT.Tmax || T<LUT.Tmin){
		throw ValueError(format("Input temperature to LookupValue_TP(%g) is out of bounds [%g,%g]",T,LUT.Tmin,LUT.Tmax));
        return _HUGE;
    }

    if (p>LUT.pmax || p<LUT.pmin)
    {
		throw ValueError(format("Input pressure to LookupValue_TP(%g) is out of bounds [%g,%g]",p,LUT.pmin,LUT.pmax));
        return _HUGE;
    }

    iTlow=(int)floor((T-LUT.Tmin)/(LUT.Tmax-LUT.Tmin)*(LUT.nT-1));
    iThigh=iTlow+1;

    iPlow=(int)floor((p-LUT.pmin)/(LUT.pmax-LUT.pmin)*(LUT.np-1));
    iPhigh=iPlow+1;

	mat = _get_LUT_ptr(Prop);    
    
	if (iThigh==LUT.nT-1){
		iThigh--;
		iTlow--;
	}

	if (iPhigh==LUT.np-1){
		iPhigh--;
		iPlow--;
	}
    //At Low Temperature Index
    y1=(*mat)[iTlow][iPlow];
    y2=(*mat)[iTlow][iPhigh];
    y3=(*mat)[iTlow][iPhigh+1];
    P1=LUT.pvec[iPlow];
    P2=LUT.pvec[iPhigh];
    P3=LUT.pvec[iPhigh+1];
    a1=QuadInterp(P1,P2,P3,y1,y2,y3,p);

    //At High Temperature Index
    y1=(*mat)[iThigh][iPlow];
    y2=(*mat)[iThigh][iPhigh];
    y3=(*mat)[iThigh][iPhigh+1];
    a2=QuadInterp(P1,P2,P3,y1,y2,y3,p);

    //At High Temperature Index+1 (for QuadInterp() )
    y1=(*mat)[iThigh+1][iPlow];
    y2=(*mat)[iThigh+1][iPhigh];
    y3=(*mat)[iThigh+1][iPhigh+1];
    a3=QuadInterp(P1,P2,P3,y1,y2,y3,p);

    //At Final Interpolation
    T1=LUT.Tvec[iTlow];
    T2=LUT.Tvec[iThigh];
    T3=LUT.Tvec[iThigh+1];

	if (debug()>8){
		std::cout << __FILE__<<"Lookup_TP for " << Prop << " with inputs T="<< T <<" and p="<< p <<std::endl;
	}
    return QuadInterp(T1,T2,T3,a1,a2,a3,T);
}

double Fluid::LookupValue_Trho(std::string Prop, double T, double rho)
{
	// Lookup a value in the lookup table using temperature and rho as the inputs

    int irholow, irhohigh, iTlow, iThigh,L,R,M;
    double T1, T2, T3, rho1, rho2, rho3, y1, y2, y3, a1, a2, a3;
	std::vector<std::vector<double> > (*mat);
    double Tmin,Tmax,rhomin,rhomax;

    BuildLookupTable(); 
    
    Tmin=LUT.Tvec[0];
    Tmax=LUT.Tvec[LUT.nT-1];

    if (T>LUT.Tmax || T<LUT.Tmin){
		throw ValueError(format("Input temperature to LookupValue_Trho(%g K) is out of bounds [%g K,%g K]",T,LUT.Tmin,LUT.Tmax));
        return _HUGE;
    }

	iTlow=(int)floor((T-LUT.Tmin)/(LUT.Tmax-LUT.Tmin)*(LUT.nT-1));
    iThigh=iTlow+1;
	
	if (iThigh==LUT.nT-1){
		iThigh--;
		iTlow--;
	}

    rhomin=LUT.rhomat[iThigh][0];
    rhomax=LUT.rhomat[iTlow][LUT.np-1];

	if (rho>rhomax || rho<rhomin)
    {
		throw ValueError(format("Input rho to LookupValue_Trho(%g kg/m3) for given T=%g K is out of bounds [%g kg/m3,%g kg/m3]",rho,T,rhomin,rhomax));
        return _HUGE;
    }

	L=0;
	R=LUT.np-1;
	M=(L+R)/2;
	// Use interval halving to find the indices which bracket the density of interest
	while (R-L>1)
	{
		if (rho>=LUT.rhomat[iThigh][M])
		{ L=M; M=(L+R)/2; continue;}
		if (rho<LUT.rhomat[iTlow][M])
		{ R=M; M=(L+R)/2; continue;}
	}
    irholow=L;
    irhohigh=R;
	
	if (irhohigh==LUT.np-1){
		irholow--;
		irhohigh--;
	}

    mat = _get_LUT_ptr(Prop);
    
    //At Low Temperature Index
    y1=(*mat)[iTlow][irholow];
    y2=(*mat)[iTlow][irhohigh];
    y3=(*mat)[iTlow][irhohigh+1];
    rho1=LUT.rhomat[iTlow][irholow];
    rho2=LUT.rhomat[iTlow][irhohigh];
    rho3=LUT.rhomat[iTlow][irhohigh+1];
    a1=QuadInterp(rho1,rho2,rho3,y1,y2,y3,rho);

    //At High Temperature Index
    y1=(*mat)[iThigh][irholow];
    y2=(*mat)[iThigh][irhohigh];
    y3=(*mat)[iThigh][irhohigh+1];
    rho1=LUT.rhomat[iThigh][irholow];
    rho2=LUT.rhomat[iThigh][irhohigh];
    rho3=LUT.rhomat[iThigh][irhohigh+1];
    a2=QuadInterp(rho1,rho2,rho3,y1,y2,y3,rho);

    //At High Temperature Index+1 (for QuadInterp() )
    y1=(*mat)[iThigh+1][irholow];
    y2=(*mat)[iThigh+1][irhohigh];
    y3=(*mat)[iThigh+1][irhohigh+1];
    rho1=LUT.rhomat[iThigh+1][irholow];
    rho2=LUT.rhomat[iThigh+1][irhohigh];
    rho3=LUT.rhomat[iThigh+1][irhohigh+1];
    a3=QuadInterp(rho1,rho2,rho3,y1,y2,y3,rho);

    //At Final Interpolation
    T1=LUT.Tvec[iTlow];
    T2=LUT.Tvec[iThigh];
    T3=LUT.Tvec[iThigh+1];

	if (debug()>8){
		std::cout << __FILE__<<"Lookup_Trho for" << Prop << "with inputs T="<< T <<" and rho="<< rho <<std::endl;
	}
    return QuadInterp(T1,T2,T3,a1,a2,a3,T);
}

//
//int WriteLookup2File(int ILUT)
//{
//    int i,j;
//    FILE *fp_h,*fp_s,*fp_rho,*fp_u,*fp_cp,*fp_cv,*fp_visc;
//    fp_h=fopen("h.csv","w");
//    fp_s=fopen("s.csv","w");
//    fp_u=fopen("u.csv","w");
//    fp_cp=fopen("cp.csv","w");
//    fp_cv=fopen("cv.csv","w");
//    fp_rho=fopen("rho.csv","w");
//    fp_visc=fopen("visc.csv","w");
//
//    // Write the pressure header row
//    for (j=0;j<nP;j++)
//    {
//        fprintf(fp_h,",%0.12f",pvec[j][ILUT]);
//        fprintf(fp_s,",%0.12f",pvec[j][ILUT]);
//        fprintf(fp_rho,",%0.12f",pvec[j][ILUT]);
//        fprintf(fp_u,",%0.12f",pvec[j][ILUT]);
//        fprintf(fp_cp,",%0.12f",pvec[j][ILUT]);
//        fprintf(fp_cv,",%0.12f",pvec[j][ILUT]);
//        fprintf(fp_visc,",%0.12f",pvec[j][ILUT]);
//    }
//    fprintf(fp_h,"\n");
//    fprintf(fp_s,"\n");
//    fprintf(fp_rho,"\n");
//    fprintf(fp_u,"\n");
//    fprintf(fp_cp,"\n");
//    fprintf(fp_cv,"\n");
//    fprintf(fp_visc,"\n");
//    
//    for (i=1;i<nT;i++)
//    {
//        fprintf(fp_h,"%0.12f",Tvec[i][ILUT]);
//        fprintf(fp_s,"%0.12f",Tvec[i][ILUT]);
//        fprintf(fp_rho,"%0.12f",Tvec[i][ILUT]);
//        fprintf(fp_u,"%0.12f",Tvec[i][ILUT]);
//        fprintf(fp_cp,"%0.12f",Tvec[i][ILUT]);
//        fprintf(fp_cv,"%0.12f",Tvec[i][ILUT]);
//        fprintf(fp_visc,"%0.12f",Tvec[i][ILUT]);
//        for (j=0;j<nP;j++)
//        {
//            fprintf(fp_h,",%0.12f",hmat[i][j][ILUT]);
//            fprintf(fp_s,",%0.12f",smat[i][j][ILUT]);
//            fprintf(fp_rho,",%0.12f",rhomat[i][j][ILUT]);
//            fprintf(fp_u,",%0.12f",umat[i][j][ILUT]);
//            fprintf(fp_cp,",%0.12f",cpmat[i][j][ILUT]);
//            fprintf(fp_cv,",%0.12f",cvmat[i][j][ILUT]);
//            fprintf(fp_visc,",%0.12f",viscmat[i][j][ILUT]);
//        }
//        fprintf(fp_h,"\n");
//        fprintf(fp_s,"\n");
//        fprintf(fp_rho,"\n");
//        fprintf(fp_u,"\n");
//        fprintf(fp_cp,"\n");
//        fprintf(fp_cv,"\n");
//        fprintf(fp_visc,"\n");
//    }
//    return 1;
//}

double Fluid::viscosity_dilute(double T, double e_k, double sigma)
{	
	// T in [K], e_k in [K], sigma in [nm]
	// viscosity returned is in [Pa-s]
	
	/*
	Model for the Viscosity and Thermal Conductivity of Refrigerants,
	Including a New Correlation for the Viscosity of R134a,
	Marcia L. Huber, Arno Laesecke, and Richard A. Perkins
	Ind. Eng. Chem. Res. 2003, 42, 3163-3178
	*/

	double sum=0,eta_star, Tstar, OMEGA_2_2, pi=3.141592654;
	Tstar = T/e_k;
	// From Neufeld, 1972, Journal of Chemical Physics - checked coefficients
	OMEGA_2_2 = 1.16145*pow(Tstar,-0.14874)+ 0.52487*exp(-0.77320*Tstar)+2.16178*exp(-2.43787*Tstar);
	// Using the leading constant from McLinden, 2000 since the leading term from Huber 2003 gives crazy values
	eta_star = 26.69e-3*sqrt(params.molemass*T)/(pow(sigma,2)*OMEGA_2_2)/1e6;
	return eta_star;
}
double Fluid::conductivity_critical(double T, double rho)
{
	double k=1.380658e-23, //[J/K]
		R0=1.03,
		gamma=1.239,
		nu=0.63,
		Pcrit = reduce.p, //[kPa]
		Tref = 1.5*reduce.T, //[K]
		GAMMA = 0.0496,
		zeta0=1.94e-10, //[m]
		qd = 2.0e9, //[1/m]
		cp,cv,delta,num,zeta,mu,
		OMEGA_tilde,OMEGA_tilde0,pi=M_PI,tau;

	delta = rho/reduce.rho;

	tau = reduce.T/T;
	double dp_drho=R()*T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta));
	double X = Pcrit/pow(reduce.rho,2)*rho/dp_drho;
	tau = reduce.T/Tref;
	double dp_drho_ref=R()*Tref*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta));
	double Xref = Pcrit/pow(reduce.rho,2)*rho/dp_drho_ref*Tref/T;
	num=X-Xref;

	// no critical enhancement if numerator is negative
	if (num<0)
		return 0.0;
	else
		zeta=zeta0*pow(num/GAMMA,nu/gamma); //[m]

	cp=specific_heat_p_Trho(T,rho); //[kJ/kg/K]
	cv=specific_heat_v_Trho(T,rho); //[kJ/kg/K]
	mu=viscosity_Trho(T,rho)*1e6; //[uPa-s]

	OMEGA_tilde=2.0/pi*((cp-cv)/cp*atan(zeta*qd)+cv/cp*(zeta*qd)); //[-]
	OMEGA_tilde0=2.0/pi*(1.0-exp(-1.0/(1.0/(qd*zeta)+1.0/3.0*(zeta*qd)*(zeta*qd)/delta/delta))); //[-]

	double lambda=rho*cp*1e9*(R0*k*T)/(6*pi*mu*zeta)*(OMEGA_tilde-OMEGA_tilde0); //[W/m/K]
	return lambda/1e3; //[kW/m/K]
}
double Fluid::surface_tension_T(double T)
{
	/* Implements the method of Miqeu et al;, "An extended scaled equation for the temperature dependence of
    the surface tension of pure compounds inferred from an analysis of experimental data",Fluid Phase 
	Equilibria 172 (2000) 169–182

	This correlation is used as the default, in the absence of other correlation, which can be provided for fluids if desired

	sigma in [mN/m]-->[N/m]
	*/
	double N_A = 6.02214129e23,    //[-]
		k = params.R_u/N_A,         //[J/K]
		w = params.accentricfactor, //[-]
		Vc,                         //[cm^3/mol]
		Tc,                         //[K]
		t,                          //[K]
		sigma;                      //[mN/m]
	Tc=reduce.T;
	t=1-T/Tc;
	//[m3/kg]*[kg/kmol]*[0.001 kmol/mol] --> [m3/mol]
	Vc = 1/reduce.rho*params.molemass/1000;
    // N_A has units of 1/mol, Vc has units of m3/mol , (1/m3)^(2/3)->1/m^2
	// k has units of J/K
	// Tc has units of K
	// k*Tc has units of J, or N-m
	sigma = k*Tc*pow(N_A/Vc,2.0/3.0)*(4.35+4.14*w)*pow(t,1.26)*(1+0.19*sqrt(t)-0.25*t);
	return sigma;
}
class ConformalTempResids : public FuncWrapperND
{
private:
	Fluid * InterestFluid, * ReferenceFluid;
	double T,rho;
public:
	ConformalTempResids(Fluid *InterestFluid_, Fluid *ReferenceFluid_, double T_, double rho_){
		InterestFluid=InterestFluid_;
		ReferenceFluid=ReferenceFluid_;
		T=T_;
		rho=rho_;
	};
	~ConformalTempResids(){};
	Eigen::VectorXd call(Eigen::VectorXd x)
	{
		double T0=x[0]; double rho0=x[1];
		double alpha_j = DerivTerms("phir",T,rho,InterestFluid->get_namec());
		double alpha_0 = DerivTerms("phir",T0,rho0,ReferenceFluid->get_namec());
		double Z_j = DerivTerms("Z",T,rho,InterestFluid->get_namec());
		double Z_0 = DerivTerms("Z",T0,rho0,ReferenceFluid->get_namec());
		Eigen::Vector2d out;
		out(0)=alpha_j-alpha_0;
		out(1)=Z_j-Z_0;
		return out;
	}
	Eigen::MatrixXd Jacobian(Eigen::VectorXd x)
	{
		double T0=x[0]; double rho0=x[1];
		double dtau_dT = -ReferenceFluid->reduce.T/T/T;
		double ddelta_drho = 1/ReferenceFluid->reduce.rho;
		Eigen::MatrixXd out(x.rows(),x.rows());
		// Terms for the fluid of interest drop out
		double dalpha_dT0 = -DerivTerms("dphir_dTau",T0,rho0,ReferenceFluid->get_namec())*dtau_dT;
		out(0,0) = dalpha_dT0;
		double dalpha_drho0 = -DerivTerms("dphir_dDelta",T0,rho0,ReferenceFluid->get_namec())*ddelta_drho;
		out(0,1) = dalpha_drho0;
		double dZ_dT0 = -DerivTerms("dZ_dTau",T0,rho0,ReferenceFluid->get_namec())*dtau_dT;
		out(1,0) = dZ_dT0;
		double dZ_drho0 = -DerivTerms("dZ_dDelta",T0,rho0,ReferenceFluid->get_namec())*ddelta_drho;
		out(1,1) = dZ_drho0;

		return out;
	}
};

Eigen::Vector2d Fluid::ConformalTemperature(Fluid *InterestFluid, Fluid *ReferenceFluid,double T0, double rho0, std::string *errstring)
{
	int iter=0;
	double error,v0,v1,T,rho,delta,tau,dp_drho;
	Eigen::Vector2d f0,v;
	Eigen::MatrixXd J;
	ConformalTempResids CTR = ConformalTempResids(InterestFluid,ReferenceFluid,T0,rho0);
	Eigen::Vector2d x0;
	x0 << T0, rho0;
	
	// Check whether the starting guess is already pretty good
	error = sqrt(CTR.call(x0).array().square().sum());
	if (fabs(error)<1e-3){
		return x0;
	}
	// Make a copy so that if the calculations fail, we can return the original values
	Eigen::Vector2d x0_initial=x0;

	try{
		// First naively try to just use the Newton-Raphson solver without any
		// special checking of the values.
		x0=NDNewtonRaphson_Jacobian(&CTR,x0_initial,1e-10,30,errstring);
		error = sqrt(CTR.call(x0).array().square().sum());
		if (fabs(error)>1e-2 || x0(0)<0.0  || x0(1)<0.0 || !ValidNumber(x0(0)) || !ValidNumber(x0(1))){
			throw ValueError("");
		}
		return x0;
	}
	catch(std::exception &e){}
	// Ok, that didn't work, so we need to try something more interesting
	// Local Newton-Raphson solver with bounds checking on the step values
	error=999;
	iter=0;
	x0=x0_initial; // Start back at unity shape factors
	while (iter<30 && fabs(error)>1e-6)
	{
		T = x0(0);
		rho = x0(1);
		tau = reduce.T/T;
		delta = rho/reduce.rho;
		dp_drho=R()*T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta));

		//Try to take a step
		f0 = CTR.call(x0);
		J = CTR.Jacobian(x0);
		v = J.inverse()*f0;
		v0 = v(0); 
		v1 = v(1);
		if (x0(0)-v(0)>0 && x0(1)-v(1)>0)
			x0 -= v;
		else
			x0 -= 1.05*x0;

		error = sqrt(f0.array().square().sum());
		iter+=1;
	}
	return x0_initial;
};

double Fluid::viscosity_ECS_Trho(double T, double rho, Fluid * ReferenceFluid)
{
	/*
	Implements the method of
	Marcia L. Huber, Arno Laesecke, and Richard A. Perkins
	"Model for the Viscosity and Thermal Conductivity of Refrigerants,
	Including a New Correlation for the Viscosity of R134a"
	Ind. Eng. Chem. Res. 2003, 42, 3163-3178
	*/
	double e_k,sigma,e0_k, sigma0, Tc0,rhoc0,T0,rho0,rhoc,Tc,
		eta_dilute,theta,phi,f,h,eta_resid,M,M0,F_eta,eta,psi,
		rhoc0bar,rhobar,rhocbar,rho0bar;
	Eigen::Vector2d x0;

	Tc0=ReferenceFluid->reduce.T;
	rhoc0=ReferenceFluid->reduce.rho;
	M0=ReferenceFluid->params.molemass;
	rhoc0bar = rhoc0/M0;
	Tc = reduce.T;
	rhoc = reduce.rho;
	M = params.molemass;
	rhocbar=rhoc/M;
	rhobar=rho/M;
	try{
		// Get the ECS params for the fluid if it has them
		ECSParams(&e_k,&sigma);
	}
	catch(NotImplementedError &e){
		try{
			// Get the ECS parameters from the reference fluid
			ReferenceFluid->ECSParams(&e0_k,&sigma0);
		}
		catch (NotImplementedError &e){
			// Doesn't have e_k and sigma for reference fluid
			throw NotImplementedError(format("Your reference fluid for ECS [%s] does not have an implementation of ECSParams",(char *)ReferenceFluid->get_name().c_str()));
		}
		//Estimate the ECS parameters from Huber and Ely, 2003
		e_k = e0_k*Tc/Tc0;
		sigma = sigma0*pow(rhoc/rhoc0,1.0/3.0);
	}

	// The dilute portion is for the fluid of interest, not for the reference fluid
	// It is the viscosity in the limit of zero density
	eta_dilute = viscosity_dilute(T,e_k,sigma);

	if (1)//(T>reduce.T)
	{
		// Get the conformal temperature.  To start out here, assume that the shape factors are unity
		theta=1;
		phi=1;
	}
	else
	{
		/* Use the method from
		Isabel M. Marruchoa, James F. Ely, "Extended corresponding states for pure polar
		and non-polar fluids: an improved method for component shape factor prediction",
		Fluid Phase Equilibria 150–151 1998 215–223

		Reynes and Thodos,
		APPLICATION OF A REDUCED VAPOR
		PRESSURE EQUATION TO
		NONHYDROROCARBON SUBSTANCES
		beta = 5/9*gamma-40/27
		gamma = 9/5*beta + 9/5*40/27
		gamma = 9/5*beta + 8/3
		*/
		double Bstar,Bstar0,Cstar,Cstar0,DELTABstar,DELTACstar,Zc,Zc0;
		double omega = params.accentricfactor;
		double omega0 = ReferenceFluid->params.accentricfactor;
		double Tstar = T/reduce.T;
		Zc = reduce.p/(reduce.rho*R()*reduce.T);
		Zc0 = ReferenceFluid->reduce.p/(ReferenceFluid->reduce.rho*ReferenceFluid->R()*ReferenceFluid->reduce.T);
		Bstar = -6.207612-15.37641*omega-0.574946*pow(10,-omega);
		Bstar0 = -6.207612-15.37641*omega0-0.574946*pow(10,-omega0);
		Cstar = 8/3+9*Bstar/5.0/log(10.0);
		Cstar0 = 8/3+9*Bstar0/5.0/log(10.0);
		DELTABstar = Bstar-Bstar0;
		DELTACstar = Cstar-Cstar0;
		theta = (1-Cstar0+2*pow(1-Tstar,2.0/7.0)*log(Zc/Zc0)-DELTABstar+DELTACstar*log(Tstar)+Bstar/Tstar)/(1-Cstar0+Bstar0/Tstar);
		phi = pow(Zc,pow(1-Tstar,2.0/7.0))/pow(Zc0,pow(1-Tstar/theta,2.0/7.0));
	}

	try{
		psi = ECS_psi_viscosity(rho/reduce.rho);
	}
	catch(NotImplementedError &e){
		psi = 1.0;
	}
	f=Tc/Tc0*theta;
	//Must be the ratio of MOLAR densities!!
	h=rhoc0bar/rhocbar*phi;
	T0=T/f;
	rho0bar=rhobar*h;
	//Get the mass density for the reference fluid [kg/m3]
	rho0=rho0bar*M0;
	std::string errstring;
	// First check whether you should use this code in the first place.
	// Implementing the method of TRNECS from REFPROP
	double delta = rho/reduce.rho;
	double tau = reduce.T/T;
	double Z = 1+delta*dphir_dDelta(tau,delta);
	double p0 = Z*R()*T0*rho0;
	if (Z<0.3 || p0>1.1*ReferenceFluid->reduce.p || rho0>ReferenceFluid->reduce.rho){
		// Use the code to calculate the conformal state
		x0=ConformalTemperature(this,ReferenceFluid,T0,rho0,&errstring);
		T0=x0(0);
		rho0=x0(1);
	}
	rho0bar = rho0/M0;
	h=rho0bar/rhobar;
	f=T/T0;

	eta_resid = ReferenceFluid->viscosity_background(T0,rho0*psi);
	F_eta = sqrt(f)*pow(h,-2.0/3.0)*sqrt(M/M0);
	eta = eta_dilute+eta_resid*F_eta;
	return eta;
}

double Fluid::conductivity_ECS_Trho(double T, double rho, Fluid * ReferenceFluid)
{
	/*
	Implements the method of
	Marcia L. Huber, Arno Laesecke, and Richard A. Perkins
	"Model for the Viscosity and Thermal Conductivity of Refrigerants,
	Including a New Correlation for the Viscosity of R134a"
	Ind. Eng. Chem. Res. 2003, 42, 3163-3178
	*/
	double e_k,sigma,e0_k, sigma0, Tc0,rhoc0,T0,rho0,rhoc,Tc,
		eta_dilute,theta,phi,f,h,lambda_resid,M,M0,F_lambda,lambda,chi,
		f_int,lambda_int,lambda_crit,lambda_star,rhoc0bar,rhobar,rhocbar,rho0bar;
	Eigen::Vector2d x0;

	// Properties for the reference fluid
	Tc0=ReferenceFluid->reduce.T;
	rhoc0=ReferenceFluid->reduce.rho;
	M0=ReferenceFluid->params.molemass;
	rhoc0bar = rhoc0/M0;

	// Properties for the given fluid
	Tc = reduce.T;
	rhoc = reduce.rho;
	M = params.molemass;
	rhocbar=rhoc/M;
	rhobar=rho/M;
	try{
		// Get the ECS params for the fluid if it has them
		ECSParams(&e_k,&sigma);
	}
	catch(NotImplementedError &e){
		try{
			//Estimate the ECS parameters from Huber and Ely, 2003
			ReferenceFluid->ECSParams(&e0_k,&sigma0);
		}
		catch (NotImplementedError &e){
			// Doesn't have e_k and sigma for reference fluid
			throw NotImplementedError(format("Your reference fluid for ECS [%s] does not have an implementation of ECSParams",(char *)ReferenceFluid->get_name().c_str()));
		}
		e_k = e0_k*Tc/Tc0;
		sigma = sigma0*pow(rhoc/rhoc0,1.0/3.0);
	}
	// The dilute portion is for the fluid of interest, not for the reference fluid
	// It is the viscosity in the limit of zero density
	// It has units of Pa-s here
	eta_dilute = viscosity_dilute(T,e_k,sigma);
	
	try{
		chi = ECS_chi_conductivity(rho/reduce.rho);
	}
	catch(NotImplementedError){
		chi = 1.0;
	}

	try{
		f_int = ECS_f_int(T);
	}
	catch(NotImplementedError){
		f_int = 1.32e-3;
	}

	// Get the conformal temperature.  To start out here, assume that the shape factors are unity
	theta=1.0;
	phi=1.0;
	f=Tc/Tc0*theta;
	h=rhoc0bar/rhocbar*phi;
	T0=T/f;
	rho0bar=rhobar*h;
	//Get the mass density for the reference fluid [kg/m3]
	rho0=rho0bar*M0;
	
	std::string errstring;
	// First check whether you should use this code in the first place.
	// Implementing the method of TRNECS from REFPROP
	double delta = rho/reduce.rho;
	double tau = reduce.T/T;
	double Z = 1+delta*dphir_dDelta(tau,delta);
	double p0 = Z*R()*T0*rho0;
	if (Z<0.3 || p0>1.1*ReferenceFluid->reduce.p || rho0>ReferenceFluid->reduce.rho){
		// Use the code to calculate the conformal state
		x0=ConformalTemperature(this,ReferenceFluid,T0,rho0,&errstring);
		T0=x0(0);
		rho0=x0(1);
	}
	
	rho0bar = rho0/M0;
	h=rho0bar/rhobar;
	f=T/T0;
	
	// Ideal-gas specific heat in the limit of zero density
	double cpstar = specific_heat_p_ideal_Trho(T);
	lambda_star = 15e-3*R()*(eta_dilute*1e6)/4.0; //[W/m/K]
	double rrrrr2 = f_int*(eta_dilute*1e6);
	lambda_int = f_int*(eta_dilute*1e6)*(cpstar-5.0/2.0*R() ); //[W/m/K]
	F_lambda = sqrt(f)*pow(h,-2.0/3.0)*sqrt(M0/M);
	
	//Get the background viscosity from the reference fluid
	lambda_resid = ReferenceFluid->conductivity_background(T0,rho0*chi)*1e3 ;//[W/m/K]
	lambda_crit = conductivity_critical(T,rho);
	lambda = lambda_int+lambda_star+lambda_resid*F_lambda+lambda_crit;
	return lambda/1e3; //[kW/m/K]
}

// ------------------------------
// FluidsContainer Implementation
// ------------------------------

FluidsContainer::FluidsContainer()
{
	FluidsList.push_back(new AirClass());
	FluidsList.push_back(new WaterClass());	
	FluidsList.push_back(new R134aClass());
	
	// If the proprocessor key ONLY_AIR_WATER is defined, only air and water will be included
	// This is to speed up compilation of humid air package since many fewer files will be included
	#if !defined(ONLY_AIR_WATER)
	// The pure fluids
	FluidsList.push_back(new ArgonClass());
	FluidsList.push_back(new R744Class());
	FluidsList.push_back(new NitrogenClass());
	FluidsList.push_back(new R290Class());
	FluidsList.push_back(new R717Class());
	FluidsList.push_back(new R1234yfClass());
	FluidsList.push_back(new R32Class());
	FluidsList.push_back(new R22Class());
	// The industrial fluids
	FluidsList.push_back(new CarbonMonoxideClass());
	FluidsList.push_back(new CarbonylSulfideClass());
	FluidsList.push_back(new DecaneClass());
	FluidsList.push_back(new HydrogenSulfideClass());
	FluidsList.push_back(new IsopentaneClass());
	FluidsList.push_back(new NeopentaneClass());
	FluidsList.push_back(new IsohexaneClass());
	FluidsList.push_back(new KryptonClass());
	FluidsList.push_back(new NonaneClass());
	FluidsList.push_back(new TolueneClass());
	FluidsList.push_back(new XenonClass());
	FluidsList.push_back(new R116Class());
	FluidsList.push_back(new AcetoneClass());
	FluidsList.push_back(new NitrousOxideClass());
	FluidsList.push_back(new SulfurDioxideClass());
	FluidsList.push_back(new R141bClass());
	FluidsList.push_back(new R142bClass());
	FluidsList.push_back(new R218Class());
	FluidsList.push_back(new R245faClass());
	FluidsList.push_back(new R41Class());
	// The pseudo-pure fluids
	FluidsList.push_back(new R404AClass());
	FluidsList.push_back(new R410AClass());
	FluidsList.push_back(new R407CClass());
	FluidsList.push_back(new R507AClass());
	#endif

	// Build the map of fluid names mapping to pointers to the Fluid class instances
	for (std::list<Fluid*>::iterator it = FluidsList.begin(); it != FluidsList.end(); it++)
	{
		// Call the post_load routine
		(*it)->post_load();
		// Load up entry in map
		fluid_name_map[(*it)->get_name()] = *it;
	}
}

// Destructor
FluidsContainer::~FluidsContainer()
{
	while (!FluidsList.empty())
	{
		delete FluidsList.back();
		FluidsList.pop_back();
	}
}

Fluid * FluidsContainer::get_fluid(std::string name)
{
	std::map<std::string,Fluid*>::iterator it;
	// Try to find using the map if Fluid name is provided
	it = fluid_name_map.find(name);
	// If it is found the iterator will not be equal to end
	if (it != fluid_name_map.end() )
	{
		// Return a pointer to the class
		return (*it).second;
	}

	// Wasn't found, now we need to check for an alias
	for (std::list<Fluid*>::iterator it = FluidsList.begin(); it != FluidsList.end(); it++)
	{
		if ( (*it)->isAlias(name) )
		{
			return *it;
		}
	}
	throw NotImplementedError(format("Fluid [%s] not allowed",name.c_str()));
	return NULL;
}

std::string FluidsContainer::FluidList()
{
	std::string FL;
	for (std::list<Fluid*>::iterator it = FluidsList.begin(); it != FluidsList.end(); it++)
	{
		FL+=(*it)->get_name();
		FL+=",";
	}
	//Remove the tailing comma
	FL = FL.substr (0,FL.length()-1);
	return FL;
}
