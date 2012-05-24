#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include <string>
#include <vector>
#include <list>
#include <exception>
#include <iostream>
#include "Helmholtz.h"
#include "FluidClass.h"
#include "CoolProp.h"
#include "CoolPropTools.h"
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
	for (list<phi_BC*>::iterator it = phi0list.begin(); it != phi0list.end(); it++)
		summer += (*it)->base(tau,delta);
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
	for (list<phi_BC*>::iterator it = phi0list.begin(); it != phi0list.end(); it++)
		summer += (*it)->dTau(tau,delta);
	return summer;
}
double Fluid::d2phi0_dTau2(double tau, double delta)
{
	double summer = 0;
	for (list<phi_BC*>::iterator it = phi0list.begin(); it != phi0list.end(); it++)
		summer += (*it)->dTau2(tau,delta);
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
    double delta,tau,c1,c2,R;
    R=params.R_u/params.molemass;
	tau=reduce.T/T;
	delta=rho/reduce.rho;

    c1=pow(1.0+delta*dphir_dDelta(tau,delta)-delta*tau*d2phir_dDelta_dTau(tau,delta),2);
    c2=(1.0+2.0*delta*dphir_dDelta(tau,delta)+pow(delta,2)*d2phir_dDelta2(tau,delta));
	double a1=d2phi0_dTau2(tau,delta);
	double a2 = d2phir_dTau2(tau,delta);
	double a3 =  R*(-pow(tau,2)*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta))+c1/c2);
    return R*(-pow(tau,2)*(d2phi0_dTau2(tau,delta)+d2phir_dTau2(tau,delta))+c1/c2);
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

double Fluid::density_Tp(double T, double p, double rho)
{
    double delta,tau,dpdrho,error=999,R,p_EOS;
    R=params.R_u/params.molemass;
	tau=reduce.T/T;
	
    int iter=1;
    while (fabs(error)>1e-6)
    {
        delta=rho/reduce.rho;
        // Use Newton's method to find the saturation density since the derivative of pressure w.r.t. density is known from EOS
        dpdrho=R*T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta));
		// Pressure from equation of state
		p_EOS = pressure_Trho(T,rho);
        // Update the step using Newton's method
        rho=rho-(p_EOS-p)/dpdrho;
        // Residual
        error=p_EOS-p;
        iter++;
        if (iter>100)
        {
            printf("Number of steps in Density_TP has exceeded 100 with inputs %g,%g,%g\n",T,p,rho);
            return _HUGE;
        }
    }		
    return rho;
}


double Fluid::gibbs_Trho(double T,double rho)
{
    return enthalpy_Trho(T,rho)-T*entropy_Trho(T,rho);
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
			rhosatPure(T,&rhoL,&rhoV,&p);
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
        rhosatPure(SatLUT.T[i], &(SatLUT.rhoL[i]), &(SatLUT.rhoV[i]), &(SatLUT.p[i]));
		SatLUT.logp[i]=log(SatLUT.p[i]);
        // Calculate the other saturation properties
        SatLUT.hL[i]=Props('H','T',SatLUT.T[i],'D',SatLUT.rhoL[i],name);
        SatLUT.hV[i]=Props('H','T',SatLUT.T[i],'D',SatLUT.rhoV[i],name);
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
        x0=(*x)[i+1];
        x1=(*x)[i];
        x2=(*x)[i-1];
        y0=(*y)[i+1];
        y1=(*y)[i];
        y2=(*y)[i-1];
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

static void swap(double *x, double *y)
{
    double temp;
    temp = *x;
    *x = *y;
    *y = temp;
}

double Fluid::_Dekker_Tsat(double T_min, double T_max, double eps, double p, double Q)
{
    double a_k,b_k,f_ak,f_bk,f_bkn1,error,T,p_EOS, 
        b_kn1,b_kp1,s,m,a_kp1,f_akp1,f_bkp1,x,Tc,x_min,x_max;

    int iter=1;

    Tc=reduce.T;
    x_min=Tc/T_max;
    x_max=Tc/T_min;

    // Loop for the solver method
    while ((iter <= 1 || fabs(error) > eps) && iter < 100)
    {
        // Start with the maximum value
        if (iter == 1)
        {
            a_k = x_max;
            x = a_k;
        }
        // End with the minimum value
        if (iter == 2)
        {
            b_k = x_min;
            x = b_k;
        }
        if (iter > 2)
            x = b_k;

        // Evaluate residual
		T=Tc / x;
		p_EOS = Props('P','T',T,'Q',Q,name.c_str());
        error = log(p_EOS/p);

        // First time through, store the outputs
        if(iter == 1)
        {
            f_ak = error;
            b_kn1 = a_k;
            f_bkn1 = error;
        }
        if (iter > 1)
        {
            f_bk = error;
            //Secant solution
            s = b_k - (b_k - b_kn1) / (f_bk - f_bkn1) * f_bk;
            //Midpoint solution
            m = (a_k + b_k) / 2.0;

            if (s > b_k && s < m)
            {
                //Use the secant solution
                b_kp1 = s;
            }
            else
            {
                //Use the midpoint solution
                b_kp1 = m;
            }

            //See if the signs of iterate and contrapoint are the same
			T=Tc/b_kp1;
			p_EOS = Props('P','T',T,'Q',Q,name.c_str());
            f_bkp1 = log(p_EOS/p);

            if (f_ak / fabs(f_ak) != f_bkp1 / fabs(f_bkp1))
            {
                // If a and b have opposite signs, 
                //  keep the same contrapoint
                a_kp1 = a_k;
                f_akp1 = f_ak;
            }
            else
            {
                //Otherwise, keep the iterate
                a_kp1 = b_k;
                f_akp1 = f_bk;
            }

            if (fabs(f_akp1) < fabs(f_bkp1))
            {
                //a_k+1 is a better guess than b_k+1, so swap a and b values
                swap(&a_kp1, &b_kp1);
                swap(&f_akp1, &f_bkp1);
            }

            //Update variables
            //Old values
            b_kn1 = b_k;
            f_bkn1 = f_bk;
            //values at this iterate
            b_k = b_kp1;
            a_k = a_kp1;
            f_ak = f_akp1;
            f_bk = f_bkp1;
        }
        iter++;
        if (iter>90 && fabs(error)>eps)
        {
            printf("Dekker for Tsat has failed with inputs (%g,%g,%g,%g,%g,%s); value: %g\n",x_min,x_max,eps,p,Q,name,error);
        }
    }
    return Tc/b_k;
}

double Fluid::Tsat(double p, double Q, double T_guess, bool UseLUT)
{
    double x1=0,x2=0,y1=0,y2=0,tau1,tau3,tau2,logp1,logp2,logp3,T1,T2,T3;
    int i;
    double Tc,Tmax,Tmin;
    // pointers to the x and y std::vectors for later
	std::vector<double> (*y),(*x);
	double x0,y0;
    
    // Do reverse interpolation in the Saturation Lookup Table
    if (UseLUT && isPure==true) { return ApplySaturationLUT("T","p",p); }

    Tc=Props(name,"Tcrit");
    Tmax=Tc-0.000001;
    Tmin=Props(name,"Ttriple")+1;
    
    // Plotting Tc/T versus log(p) tends to give very close to straight line
    // Use this fact to figure out a reasonable guess temperature
    
    if (isPure==true)
    {
        T1=Tmin+1;
        T3=Tc-1;
        T2=(T1+T3)/2;
        tau1=Tc/T1;
        tau2=Tc/T2;
        tau3=Tc/T3;
        logp1=log(psat(T1));
        logp2=log(psat(T2));
        logp3=log(psat(T3));

        T_guess=Tc/QuadInterp(logp1,logp2,logp3,tau1,tau2,tau3,log(p));
        if (T_guess+5<Tmax)
            Tmax=T_guess+5;
        if (T_guess-5>Tmin)
            Tmin=T_guess-5;
    }
    return _Dekker_Tsat(Tmin,Tmax,0.0001,p,Q);
}