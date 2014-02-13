#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#include "float.h"
#else
#include <stdlib.h>
#include <limits>
#define DBL_EPSILON std::numeric_limits<double>::epsilon()
#endif

#include "time.h"
#include <vector>
#include <iostream>
#include "math.h"
#include "Helmholtz.h"
#include "CoolPropTools.h"

#ifndef DISABLE_CATCH
#include "Catch/catch.hpp"
#endif

void check_derivatives(phi_BC * phi, double tau, double delta, double ddelta, double dtau)
{
	double dphir_dDelta_a = (phi->base(tau,delta+ddelta)-phi->base(tau,delta-ddelta))/(2*ddelta);
	double dphir_dDelta_e = phi->dDelta(tau,delta);

	double d2phir_dDelta2_a = (phi->dDelta(tau,delta+ddelta)-phi->dDelta(tau,delta-ddelta))/(2*ddelta);
	double d2phir_dDelta2_e = phi->dDelta2(tau,delta);

	double d2phir_dDelta_dTau_a = (phi->dDelta(tau+dtau,delta)-phi->dDelta(tau-dtau,delta))/(2*dtau);
	double d2phir_dDelta_dTau_e = phi->dDelta_dTau(tau,delta);

	double d3phir_dDelta2_dTau_a = (phi->dDelta_dTau(tau,delta+ddelta)-phi->dDelta_dTau(tau,delta-ddelta))/(2*ddelta);
	double d3phir_dDelta2_dTau_e = phi->dDelta2_dTau(tau,delta);

	double d3phir_dDelta_dTau2_a = (phi->dDelta_dTau(tau+dtau,delta)-phi->dDelta_dTau(tau-dtau,delta))/(2*dtau);
	double d3phir_dDelta_dTau2_e = phi->dDelta_dTau2(tau,delta);

	double d3phir_dDelta3_a = (phi->dDelta2(tau,delta+ddelta)-phi->dDelta2(tau,delta-ddelta))/(2*ddelta);
	double d3phir_dDelta3_e = phi->dDelta3(tau,delta);

	double dphir_dTau_a = (phi->base(tau+dtau,delta)-phi->base(tau-dtau,delta))/(2*dtau);
	double dphir_dTau_e = phi->dTau(tau,delta);

	double d2phir_dTau2_a = (phi->dTau(tau+dtau,delta)-phi->dTau(tau-dtau,delta))/(2*dtau);
	double d2phir_dTau2_e = phi->dTau2(tau,delta);

	double d3phir_dTau3_a = (phi->dTau2(tau+dtau,delta)-phi->dTau2(tau-dtau,delta))/(2*dtau);
	double d3phir_dTau3_e = phi->dTau3(tau,delta);

	double d2phir_dDelta2 = 0;
}

// Constructors
phir_power::phir_power(std::vector<double> n_in,std::vector<double> d_in,std::vector<double> t_in, std::vector<double> l_in, int iStart_in,int iEnd_in)
{
	n=n_in;
	d=d_in;
	t=t_in;
	l=l_in;
	iStart=iStart_in;
	iEnd=iEnd_in;
	cache();
}
phir_power::phir_power(std::vector<double> n_in,std::vector<double> d_in,std::vector<double> t_in,int iStart_in,int iEnd_in)
{
	n=n_in;
	d=d_in;
	t=t_in;
	l.assign(d.size(),0.0);
	iStart=iStart_in;
	iEnd=iEnd_in;
	cache();
}
phir_power::phir_power(const double n_in[], const double d_in[], const double t_in[],int iStart_in,int iEnd_in, int N)
{
	n=std::vector<double>(n_in,n_in+N);
	d=std::vector<double>(d_in,d_in+N);
	t=std::vector<double>(t_in,t_in+N);
	l.assign(d.size(),0.0);
	iStart=iStart_in;
	iEnd=iEnd_in;
	cache();
}
phir_power::phir_power(double n_in[], double d_in[], double t_in[],int iStart_in,int iEnd_in, int N)
{
	n=std::vector<double>(n_in,n_in+N);
	d=std::vector<double>(d_in,d_in+N);
	t=std::vector<double>(t_in,t_in+N);
	l.assign(d.size(),0.0);
	iStart=iStart_in;
	iEnd=iEnd_in;
	cache();
}
phir_power::phir_power(double n_in[], double d_in[], double t_in[], double l_in[], int iStart_in,int iEnd_in, int N)
{
	n=std::vector<double>(n_in,n_in+N);
	d=std::vector<double>(d_in,d_in+N);
	t=std::vector<double>(t_in,t_in+N);
	l=std::vector<double>(l_in,l_in+N);
	iStart=iStart_in;
	iEnd=iEnd_in;
	cache();
}
phir_power::phir_power(const double n_in[], const double d_in[], const double t_in[], const double l_in[], int iStart_in,int iEnd_in, int N)
{
	n=std::vector<double>(n_in,n_in+N);
	d=std::vector<double>(d_in,d_in+N);
	t=std::vector<double>(t_in,t_in+N);
	l=std::vector<double>(l_in,l_in+N);
	iStart=iStart_in;
	iEnd=iEnd_in;
	cache();
}
void phir_power::cache()
{
	// Check that taking the integer representation of each of the powers of delta (delta^l_i) all yield integer
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if ( fabs((double)((int)(l[i])) - l[i]) > DBL_EPSILON)
		{
			throw ValueError(format("coefficient %0.16f does not round to integer within DBL_EPSILON",l[i]).c_str());
		}
	}
}

// Term and its derivatives
double phir_power::base(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i]>0)
			summer+=n[i]*exp(t[i]*log_tau+d[i]*log_delta-pow(delta,(int)l[i]));
		else
			summer+=n[i]*exp(t[i]*log_tau+d[i]*log_delta);
	}
	return summer;
}

double phir_power::dTau(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i]>0)
			summer+=n[i]*t[i]*exp((t[i]-1)*log_tau+d[i]*log_delta-pow(delta,(int)l[i]));
		else
			summer+=n[i]*t[i]*exp((t[i]-1)*log_tau+d[i]*log_delta);
	}
	return summer;
}
double phir_power::dTau2(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i]>0)
			summer+=n[i]*t[i]*(t[i]-1)*exp((t[i]-2)*log_tau+d[i]*log_delta-pow(delta,(int)l[i]));
		else
			summer+=n[i]*t[i]*(t[i]-1)*exp((t[i]-2)*log_tau+d[i]*log_delta);
	}
	return summer;
}
double phir_power::dTau3(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i]>0)
			summer+=n[i]*t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta-pow(delta,(int)l[i]));
		else
			summer+=n[i]*t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta);
	}
	return summer;
}
double phir_power::dDelta_dTau2(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i]>0){
			double pow_delta_li = pow(delta,(int)l[i]);
			summer+=n[i]*t[i]*(t[i]-1)*(d[i]-l[i]*pow_delta_li)*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta-pow_delta_li);
		}
		else
			summer+=n[i]*t[i]*(t[i]-1)*d[i]*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta);
	}
	return summer;
}
double phir_power::dDelta(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta), pow_delta_li, li, ni, di, ti;
	for (unsigned int i=iStart;i<=iEnd;++i)
	{	
		ni = n[i]; di = d[i]; ti = t[i]; li = l[i]; 
		if (li > 0){
			pow_delta_li = pow(delta,(int)li);
			summer += ni*(di-li*pow_delta_li)*exp(ti*log_tau+(di-1)*log_delta-pow_delta_li);
		}
		else
		{
			summer += ni*di*exp(ti*log_tau+(di-1)*log_delta);
		}
	}
	return summer;
}
double phir_power::A(double tau, double delta, int i) throw()
{
	double log_tau = log(tau), log_delta = log(delta);
	if (l[i]>0)
		return exp(t[i]*log_tau+d[i]*log_delta-pow(delta,(int)l[i]));
	else
		return exp(t[i]*log_tau+d[i]*log_delta);
}
double phir_power::dA_dDelta(double tau, double delta, int i) throw()
{
	double log_tau = log(tau), log_delta = log(delta), pow_delta_li, li, ni, di, ti;
	ni = n[i]; di = d[i]; ti = t[i]; li = l[i]; 
	if (li > 0){
		pow_delta_li = pow(delta,(int)li);
		return (di-li*pow_delta_li)*exp(ti*log_tau+(di-1)*log_delta-pow_delta_li);
	}
	else
	{
		return di*exp(ti*log_tau+(di-1)*log_delta);
	}
}
double phir_power::dDelta2(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		
		if (l[i]>0){
			double pow_delta_li = pow(delta,(int)l[i]);
			summer+=n[i]*((d[i]-l[i]*pow_delta_li)*(d[i]-1.0-l[i]*pow_delta_li) - l[i]*l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-2)*log_delta-pow_delta_li);
		}
		else
			summer+=n[i]*d[i]*(d[i]-1.0)*exp(t[i]*log_tau+(d[i]-2)*log_delta);
	}
	return summer;
}
double phir_power::dDelta3(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i]>0)
		{
			double pow_delta_li = pow(delta,(int)l[i]);
			double bracket = (d[i]*(d[i]-1)*(d[i]-2)+pow_delta_li*(-2*l[i]+6*d[i]*l[i]-3*d[i]*d[i]*l[i]-3*d[i]*l[i]*l[i]+3*l[i]*l[i]-l[i]*l[i]*l[i])+pow_delta_li*pow_delta_li*(3*d[i]*l[i]*l[i]-3*l[i]*l[i]+3*l[i]*l[i]*l[i])-l[i]*l[i]*l[i]*pow_delta_li*pow_delta_li*pow_delta_li);
			summer+=n[i]*bracket*exp(t[i]*log_tau+(d[i]-3)*log_delta-pow_delta_li);
		}
		else
			summer+=n[i]*d[i]*(d[i]-1.0)*(d[i]-2)*exp(t[i]*log_tau+(d[i]-3)*log_delta);
	}
	return summer;
}

double phir_power::dDelta2_dTau(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i]>0){
			double pow_delta_li = pow(delta,(int)l[i]);
			summer+=n[i]*t[i]*(((d[i]-l[i]*pow_delta_li))*(d[i]-1-l[i]*pow_delta_li)-l[i]*l[i]*pow_delta_li)*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta-pow_delta_li);
		}
		else
			summer+=n[i]*d[i]*t[i]*(d[i]-1)*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta);
	}
	return summer;
}
double phir_power::dDelta_dTau(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i]>0){
			double pow_delta_li = pow(delta,(int)l[i]);
			summer+=n[i]*t[i]*(d[i]-l[i]*pow_delta_li)*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta-pow_delta_li);
		}
		else
			summer+=n[i]*d[i]*t[i]*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta);
	}
	return summer;
}
std::vector<double> phir_power::dDeltaV(std::vector<double> tau, std::vector<double> delta) throw()
{
	std::vector<double> out = tau;
	for (int i = 0; i < (int)tau.size(); i++)
	{
		out[i] = this->dDelta(tau[i],delta[i]);
	}
	return out;
}
std::vector<double> phir_power::dDelta2V(std::vector<double> tau, std::vector<double> delta) throw()
{
	std::vector<double> out = tau;
	for (int i = 0; i < (int)tau.size(); i++)
	{
		out[i] = this->dDelta2(tau[i],delta[i]);
	}
	return out;
}
std::vector<double> phir_power::dTau2V(std::vector<double> tau, std::vector<double> delta) throw()
{
	std::vector<double> out = tau;
	for (int i = 0; i < (int)tau.size(); i++)
	{
		out[i] = this->dTau2(tau[i],delta[i]);
	}
	return out;
}
std::vector<double> phir_power::dDelta_dTauV(std::vector<double> tau, std::vector<double> delta) throw()
{
	std::vector<double> out = tau;
	for (int i = 0; i < (int)tau.size(); i++)
	{
		out[i] = this->dDelta_dTau(tau[i],delta[i]);
	}
	return out;
}

#ifndef DISABLE_CATCH
TEST_CASE("Power Helmholtz terms", "[helmholtz],[fast]")
{
	// From R134a
	double n[]={0.0, 0.5586817e-1, 0.4982230e0, 0.2458698e-1, 0.8570145e-3, 0.4788584e-3, -0.1800808e1, 0.2671641e0, -0.4781652e-1, 0.1423987e-1, 0.3324062e0, -0.7485907e-2, 0.1017263e-3, -0.5184567e+0, -0.8692288e-1, 0.2057144e+0, -0.5000457e-2, 0.4603262e-3, -0.3497836e-2, 0.6995038e-2, -0.1452184e-1, -0.1285458e-3};
	double d[]={0,2,1,3,6,6,1,1,2,5,2,2,4,1,4,1,2,4,1,5,3,10};
	double t[]={0.0,-1.0/2.0,0.0,0.0,0.0,3.0/2.0,3.0/2.0,2.0,2.0,1.0,3.0,5.0,1.0,5.0,5.0,6.0,10.0,10.0,10.0,18.0,22.0,50.0};
	double c[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,2.0,2.0,2.0,2.0,2.0,2.0,3.0,3.0,3.0,4.0};
	
	phir_power phir(n,d,t,c,1,21,22);
	double eps = sqrt(DBL_EPSILON);

	SECTION("dDelta")
	{
		double ANA = phir.dDelta(0.5, 0.5);
		double NUM = (phir.base(0.5, 0.5+eps) - phir.base(0.5,0.5-eps))/(2*eps);
		REQUIRE(abs(NUM-ANA) < 1e-6);
	}
	SECTION("dA_dDelta")
	{
		double ANA = phir.dA_dDelta(0.5, 0.5, 1);
		double NUM = (phir.A(0.5, 0.5+eps, 1) - phir.A(0.5,0.5-eps, 1))/(2*eps);
		REQUIRE(abs(NUM-ANA) < 1e-6);
	}
	SECTION("dTau")
	{
		double ANA = phir.dTau(0.5, 0.5);
		double NUM = (phir.base(0.5+eps, 0.5) - phir.base(0.5-eps,0.5))/(2*eps);
		REQUIRE(abs(NUM-ANA) < 1e-6);
	}
	SECTION("dDelta2")
	{
		double ANA = phir.dDelta2(0.5, 0.5);
		double NUM = (phir.dDelta(0.5, 0.5+eps) - phir.dDelta(0.5,0.5-eps))/(2*eps);
		REQUIRE(abs(NUM-ANA) < 1e-6);
	}
	SECTION("dTau2")
	{
		double ANA = phir.dTau2(0.5, 0.5);
		double NUM = (phir.dTau(0.5+eps, 0.5) - phir.dTau(0.5-eps,0.5))/(2*eps);
		REQUIRE(abs(NUM-ANA) < 1e-6);
	}
	SECTION("dDeltadTau")
	{
		double ANA = phir.dDelta_dTau(0.5, 0.5);
		double NUM = (phir.dTau(0.5, 0.5+eps) - phir.dTau(0.5,0.5-eps))/(2*eps);
		REQUIRE(abs(NUM-ANA) < 1e-6);
	}
}
#endif

// Constructors
phir_exponential::phir_exponential(std::vector<double> n, std::vector<double> d, std::vector<double> t, std::vector<double> l, std::vector<double> g, int iStart_in,int iEnd_in)
{
	this->n = n;
	this->d = d;
	this->t = t;
	this->l = l;
	this->g = g;
	iStart=iStart_in;
	iEnd=iEnd_in;
}
phir_exponential::phir_exponential(double n[], double d[], double t[], double l[], double g[], int iStart_in,int iEnd_in, int N)
{
	this->n=std::vector<double>(n, n+N);
	this->d=std::vector<double>(d, d+N);
	this->t=std::vector<double>(t, t+N);
	this->l=std::vector<double>(l, l+N);
	this->g=std::vector<double>(g, g+N);
	iStart=iStart_in;
	iEnd=iEnd_in;
}
phir_exponential::phir_exponential(const double n[], const double d[], const double t[], const double l[], const double g[], int iStart_in,int iEnd_in, int N)
{
	this->n = std::vector<double>(n, n+N);
	this->d = std::vector<double>(d, d+N);
	this->t = std::vector<double>(t, t+N);
	this->l = std::vector<double>(l, l+N);
	this->g = std::vector<double>(g, g+N);
	iStart = iStart_in;
	iEnd = iEnd_in;
}

// Term and its derivatives
double phir_exponential::base(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		summer += n[i]*exp(t[i]*log_tau+d[i]*log_delta-g[i]*pow(delta,l[i]));
	}
	return summer;
}
double phir_exponential::dTau(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		summer += n[i]*t[i]*exp((t[i]-1)*log_tau+d[i]*log_delta-g[i]*pow(delta,l[i]));
	}
	return summer;
}
double phir_exponential::dTau2(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		summer+=n[i]*t[i]*(t[i]-1)*exp((t[i]-2)*log_tau+d[i]*log_delta-g[i]*pow(delta,l[i]));
	}
	return summer;
}
double phir_exponential::dTau3(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		summer+=n[i]*t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta-g[i]*pow(delta,l[i]));
	}
	return summer;
}
double phir_exponential::dDelta_dTau2(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		double pow_delta_li = pow(delta,l[i]);
		summer+=n[i]*t[i]*(t[i]-1)*(d[i]-g[i]*l[i]*pow_delta_li)*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta-g[i]*pow_delta_li);
	}
	return summer;
}
double phir_exponential::dDelta(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta), pow_delta_li;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{	
		pow_delta_li = pow(delta,l[i]);
		summer += n[i]*(d[i]-g[i]*l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-1)*log_delta-g[i]*pow_delta_li);
	}
	return summer;
}
double phir_exponential::dDelta2(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta), pow_delta_li;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{	
		pow_delta_li = pow(delta,l[i]);
		// Typo in Span, 2000, re-derived from Sympy
		double bracket = d[i]*d[i] - 2*d[i]*pow(delta,l[i])*g[i]*l[i] - d[i] + pow(delta,2*l[i])*g[i]*g[i]*l[i]*l[i] - pow(delta,l[i])*g[i]*l[i]*l[i] + pow(delta,l[i])*g[i]*l[i];
		summer += n[i]*bracket*exp(t[i]*log_tau+(d[i]-2)*log_delta-g[i]*pow_delta_li);
	}
	return summer;
}
double phir_exponential::dDelta3(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		// >> n_i, tau, t_i, d_i, delta, g_i, l_i = symbols(' n_i tau t_i d_i delta g_i l_i')
		// >> phir = n_i*tau**t_i*delta**d_i*exp(-g_i*pow(delta,l_i))
		// >> simplify(diff(diff(diff(phir,delta),delta),delta))
		double pow_delta_li = pow(delta,l[i]);
		double pow_delta_2li = pow(delta,2*l[i]);
		double pow_delta_3li = pow(delta,3*l[i]);
		double bracket = d[i]*d[i]*d[i] - 3*d[i]*d[i]*pow_delta_li*g[i]*l[i] - 3*d[i]*d[i] + 3*d[i]*pow_delta_2li*g[i]*g[i]*l[i]*l[i] - 3*d[i]*pow_delta_li*g[i]*l[i]*l[i] + 6*d[i]*pow_delta_li*g[i]*l[i] + 2*d[i] - pow_delta_3li*g[i]*g[i]*g[i]*l[i]*l[i]*l[i] + 3*pow_delta_2li*g[i]*g[i]*l[i]*l[i]*l[i] - 3*pow_delta_2li*g[i]*g[i]*l[i]*l[i] - pow_delta_li*g[i]*l[i]*l[i]*l[i] + 3*pow_delta_li*g[i]*l[i]*l[i] - 2*pow_delta_li*g[i]*l[i];
		summer += n[i]*bracket*exp(t[i]*log_tau+(d[i]-3)*log_delta-g[i]*pow_delta_li);
	}
	return summer;
}
double phir_exponential::dDelta2_dTau(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		double pow_delta_li = pow(delta,l[i]);
		// Typo in Span, 2000, re-derived from Sympy
		double bracket = d[i]*d[i] - 2*d[i]*pow(delta,l[i])*g[i]*l[i] - d[i] + pow(delta,2*l[i])*g[i]*g[i]*l[i]*l[i] - pow(delta,l[i])*g[i]*l[i]*l[i] + pow(delta,l[i])*g[i]*l[i];
		summer += n[i]*t[i]*bracket*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta-g[i]*pow_delta_li);
	}
	return summer;
}
double phir_exponential::dDelta_dTau(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		double pow_delta_li = pow(delta,l[i]);
		summer+=n[i]*t[i]*(d[i]-g[i]*l[i]*pow_delta_li)*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta-g[i]*pow_delta_li);
	}
	return summer;
}

// Constructors
phir_exponential2::phir_exponential2(std::vector<double> n, std::vector<double> d, std::vector<double> l, std::vector<double> a, std::vector<double> g, int iStart_in,int iEnd_in)
{
	this->n = n;
	this->d = d;
	this->l = l;
	this->a = a;
	this->g = g;
	iStart=iStart_in;
	iEnd=iEnd_in;
}
phir_exponential2::phir_exponential2(double n[], double d[], double l[], double a[], double g[], int iStart_in,int iEnd_in, int N)
{
	this->n=std::vector<double>(n, n+N);
	this->d=std::vector<double>(d, d+N);
	this->a=std::vector<double>(a, a+N);
	this->l=std::vector<double>(l, l+N);
	this->g=std::vector<double>(g, g+N);
	iStart=iStart_in;
	iEnd=iEnd_in;
}
phir_exponential2::phir_exponential2(const double n[], const double d[], const double l[], const double a[], const double g[], int iStart_in,int iEnd_in, int N)
{
	this->n = std::vector<double>(n, n+N);
	this->d = std::vector<double>(d, d+N);
	this->a = std::vector<double>(a, a+N);
	this->l = std::vector<double>(l, l+N);
	this->g = std::vector<double>(g, g+N);
	iStart = iStart_in;
	iEnd = iEnd_in;
}

// Term and its derivatives
double phir_exponential2::base(double tau, double delta) throw()
{
	double summer=0;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		summer += n[i]*pow(delta,d[i])*exp(a[i]*tau-g[i]*pow(delta,l[i]));
	}
	return summer;
}
double phir_exponential2::dTau(double tau, double delta) throw()
{
	double summer=0;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		summer += n[i]*pow(delta,d[i])*a[i]*exp(a[i]*tau-g[i]*pow(delta,l[i]));
	}
	return summer;
}
double phir_exponential2::dTau2(double tau, double delta) throw()
{
	double summer = 0;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		summer += n[i]*pow(delta,d[i])*a[i]*a[i]*exp(a[i]*tau-g[i]*pow(delta,l[i]));
	}
	return summer;
}
double phir_exponential2::dTau3(double tau, double delta) throw()
{
	double summer = 0;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		summer += n[i]*pow(delta,d[i])*a[i]*a[i]*a[i]*exp(a[i]*tau-g[i]*pow(delta,l[i]));
	}
	return summer;
}
double phir_exponential2::dDelta_dTau2(double tau, double delta) throw()
{
	double summer = 0;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		double pow_delta_li = pow(delta,l[i]);
		summer += n[i]*a[i]*a[i]*pow(delta,d[i]-1)*(d[i]-g[i]*l[i]*pow_delta_li)*exp(a[i]*tau-g[i]*pow_delta_li);
	}
	return summer;
}
double phir_exponential2::dDelta(double tau, double delta) throw()
{
	double summer=0, pow_delta_li;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{	
		pow_delta_li = pow(delta,l[i]);
		summer += n[i]*pow(delta,d[i]-1)*(d[i]-g[i]*l[i]*pow_delta_li)*exp(a[i]*tau-g[i]*pow_delta_li);
	}
	return summer;
}
double phir_exponential2::dDelta2(double tau, double delta) throw()
{
	double summer = 0;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		double pow_delta_li = pow(delta,l[i]);
		// Typo in Span, 2000, re-derived from Sympy
		double bracket = d[i]*d[i] - 2*d[i]*pow(delta,l[i])*g[i]*l[i] - d[i] + pow(delta,2*l[i])*g[i]*g[i]*l[i]*l[i] - pow(delta,l[i])*g[i]*l[i]*l[i] + pow(delta,l[i])*g[i]*l[i];
		summer += n[i]*pow(delta,d[i]-2)*bracket*exp(a[i]*tau-g[i]*pow_delta_li);
	}
	return summer;
}
double phir_exponential2::dDelta3(double tau, double delta) throw()
{
	double summer = 0;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		// >> Code for sympy: (same result for bracketed term as phir_exponential
		// >> n_i, tau, t_i, d_i, delta, g_i, l_i, alpha_i = symbols(' n_i tau t_i d_i delta g_i l_i, alpha_i');
		// >> phir = n_i*delta**d_i*exp(alpha_i*tau-g_i*pow(delta,l_i))
		// >> simplify(diff(diff(diff(phir,delta),delta),delta))
		double pow_delta_li = pow(delta,l[i]);
		double pow_delta_2li = pow(delta,2*l[i]);
		double pow_delta_3li = pow(delta,3*l[i]);
		double bracket = d[i]*d[i]*d[i] - 3*d[i]*d[i]*pow_delta_li*g[i]*l[i] - 3*d[i]*d[i] + 3*d[i]*pow_delta_2li*g[i]*g[i]*l[i]*l[i] - 3*d[i]*pow_delta_li*g[i]*l[i]*l[i] + 6*d[i]*pow_delta_li*g[i]*l[i] + 2*d[i] - pow_delta_3li*g[i]*g[i]*g[i]*l[i]*l[i]*l[i] + 3*pow_delta_2li*g[i]*g[i]*l[i]*l[i]*l[i] - 3*pow_delta_2li*g[i]*g[i]*l[i]*l[i] - pow_delta_li*g[i]*l[i]*l[i]*l[i] + 3*pow_delta_li*g[i]*l[i]*l[i] - 2*pow_delta_li*g[i]*l[i];
		summer += n[i]*bracket*exp(a[i]*tau-g[i]*pow(delta,l[i]));
	}
	return summer;
}
double phir_exponential2::dDelta2_dTau(double tau, double delta) throw()
{
	double summer = 0;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		double pow_delta_li = pow(delta,l[i]);
		// Typo in Span, 2000, re-derived from Sympy
		double bracket = d[i]*d[i] - 2*d[i]*pow(delta,l[i])*g[i]*l[i] - d[i] + pow(delta,2*l[i])*g[i]*g[i]*l[i]*l[i] - pow(delta,l[i])*g[i]*l[i]*l[i] + pow(delta,l[i])*g[i]*l[i];
		summer += n[i]*a[i]*pow(delta,d[i]-2)*bracket*exp(a[i]*tau-g[i]*pow_delta_li);
	}
	return summer;
}
double phir_exponential2::dDelta_dTau(double tau, double delta) throw()
{
	double summer=0;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		double pow_delta_li = pow(delta,l[i]);
		summer += n[i]*a[i]*pow(delta,d[i]-1)*(d[i]-g[i]*l[i]*pow_delta_li)*exp(a[i]*tau-g[i]*pow_delta_li);
	}
	return summer;
}

// Constructors
phir_Lemmon2005::phir_Lemmon2005(std::vector<double> n,std::vector<double> d,std::vector<double> t, std::vector<double> l, std::vector< double> m, int iStart_in,int iEnd_in)
{
	this->n = n;
	this->d = d;
	this->t = t;
	this->l = l;
	this->m = m;
	iStart=iStart_in;
	iEnd=iEnd_in;
}
phir_Lemmon2005::phir_Lemmon2005(double n[], double d[], double t[], double l[], double m[], int iStart_in,int iEnd_in, int N)
{
	this->n=std::vector<double>(n,n+N);
	this->d=std::vector<double>(d,d+N);
	this->t=std::vector<double>(t,t+N);
	this->l=std::vector<double>(l,l+N);
	this->m=std::vector<double>(m,m+N);
	iStart=iStart_in;
	iEnd=iEnd_in;
}
phir_Lemmon2005::phir_Lemmon2005(const double n[], const double d[], const double t[], const double l[], const double m[], int iStart_in,int iEnd_in, int N)
{
	this->n=std::vector<double>(n,n+N);
	this->d=std::vector<double>(d,d+N);
	this->t=std::vector<double>(t,t+N);
	this->l=std::vector<double>(l,l+N);
	this->m=std::vector<double>(m,m+N);
	iStart=iStart_in;
	iEnd=iEnd_in;
}

// Term and its derivatives
double phir_Lemmon2005::base(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i] != 0 && m[i] != 0)
			summer += n[i]*exp(t[i]*log_tau+d[i]*log_delta-pow(delta,l[i])-pow(tau,m[i]));
		else if (l[i] != 0 && m[i] == 0)
			summer += n[i]*exp(t[i]*log_tau+d[i]*log_delta-pow(delta,l[i]));
		else
			summer += n[i]*exp(t[i]*log_tau+d[i]*log_delta);
	}
	return summer;
}
double phir_Lemmon2005::dTau(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta), pow_tau_mi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i] != 0 && m[i] != 0){
			pow_tau_mi = pow(tau,m[i]);
			summer += n[i]*(t[i]-m[i]*pow_tau_mi)*exp((t[i]-1)*log_tau+d[i]*log_delta-pow(delta,l[i])-pow_tau_mi);
		}
		else if (l[i] != 0 && m[i] == 0)
			summer += n[i]*t[i]*exp((t[i]-1)*log_tau+d[i]*log_delta-pow(delta,l[i]));
		else
			summer += n[i]*t[i]*exp((t[i]-1)*log_tau+d[i]*log_delta);
	}
	return summer;
}
double phir_Lemmon2005::dTau2(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta), pow_tau_mi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i] != 0 && m[i] != 0){
			pow_tau_mi = pow(tau,m[i]);
			double bracket = (t[i]-m[i]*pow_tau_mi)*(t[i]-1-m[i]*pow_tau_mi)-m[i]*m[i]*pow_tau_mi;
			summer+=n[i]*bracket*exp((t[i]-2)*log_tau+d[i]*log_delta-pow(delta,l[i])-pow_tau_mi);
		}
		else if (l[i] != 0 && m[i] == 0)
			summer+=n[i]*t[i]*(t[i]-1)*exp((t[i]-2)*log_tau+d[i]*log_delta-pow(delta,l[i]));
		else
			summer+=n[i]*t[i]*(t[i]-1)*exp((t[i]-2)*log_tau+d[i]*log_delta);
	}
	return summer;
}

/*!

\f[
\frac{{{\partial ^2}{\alpha ^r}}}{{\partial {\tau ^2}}} = {N_k}{\delta ^{{d_k}}}{\tau ^{{t_k} - 2}}\exp \left( { - {\delta ^{{l_k}}}} \right)\exp \left( { - {\tau ^{{m_k}}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
\f[
\frac{{{\partial ^2}{\alpha ^r}}}{{\partial {\tau ^2}}} = {N_k}{\delta ^{{d_k}}}\exp \left( { - {\delta ^{{l_k}}}} \right){\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
Group all the terms that don't depend on \$ \tau \$
\f[
\frac{{{\partial ^2}{\alpha ^r}}}{{\partial {\tau ^2}}} = A{\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
\f[
\frac{1}{A}\frac{{{\partial ^3}{\alpha ^r}}}{{\partial {\tau ^3}}} = {\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\frac{\partial }{{\partial \tau }}\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right] + \frac{\partial }{{\partial \tau }}\left[ {{\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)} \right]\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
\f[
\frac{\partial }{{\partial \tau }}\left[ {{\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)} \right] = ({t_k} - 2){\tau ^{{t_k} - 3}}\exp \left( { - {\tau ^{{m_k}}}} \right) + {\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)( - {m_k}{\tau ^{{m_k} - 1}}) = \exp \left( { - {\tau ^{{m_k}}}} \right)\left( {({t_k} - 2){\tau ^{{t_k} - 3}} - {\tau ^{{t_k} - 2}}{m_k}{\tau ^{{m_k} - 1}}} \right)\\
\f]
\f[
\frac{\partial }{{\partial \tau }}\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right] = \left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( { - m_k^2{\tau ^{{m_k} - 1}}} \right) + \left( { - m_k^2{\tau ^{{m_k} - 1}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^3{\tau ^{{m_k} - 1}} =  - m_k^2{\tau ^{{m_k} - 1}}\left[ {{t_k} - {m_k}{\tau ^{{m_k}}} + {t_k} - 1 - {m_k}{\tau ^{{m_k}}} + {m_k}} \right] =  - m_k^2{\tau ^{{m_k} - 1}}\left[ {2{t_k} - 2{m_k}{\tau ^{{m_k}}} - 1 + {m_k}} \right]\\
\f]
\f[
\frac{1}{A}\frac{{{\partial ^3}{\alpha ^r}}}{{\partial {\tau ^3}}} = {\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\left( { - m_k^2{\tau ^{{m_k} - 1}}\left[ {2{t_k} - 2{m_k}{\tau ^{{m_k}}} - 1 + {m_k}} \right]} \right) + \exp \left( { - {\tau ^{{m_k}}}} \right)\left( {({t_k} - 2){\tau ^{{t_k} - 3}} - {\tau ^{{t_k} - 2}}{m_k}{\tau ^{{m_k} - 1}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
\f[
\frac{1}{A}\frac{{{\partial ^3}{\alpha ^r}}}{{\partial {\tau ^3}}} = \exp \left( { - {\tau ^{{m_k}}}} \right)\left[ { - {\tau ^{{t_k} - 2}}m_k^2{\tau ^{{m_k} - 1}}\left[ {2{t_k} - 2{m_k}{\tau ^{{m_k}}} - 1 + {m_k}} \right] + \left( {({t_k} - 2){\tau ^{{t_k} - 3}} - {\tau ^{{t_k} - 2}}{m_k}{\tau ^{{m_k} - 1}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]} \right]
\f]
*/
double phir_Lemmon2005::dTau3(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta),pow_delta_li, pow_tau_mi, bracket;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i] != 0 && m[i] != 0){
			pow_delta_li = pow(delta,l[i]);
			pow_tau_mi = pow(tau,m[i]);
			bracket = -pow(tau,t[i]+m[i]-3)*m[i]*m[i]*(2*t[i]-2*m[i]*pow_tau_mi-1-m[i])+((t[i]-2)*pow(tau,t[i]-3)-pow(tau,t[i]-2)*m[i]*pow(tau,m[i]-1))*((t[i]-m[i]*pow_tau_mi)*(t[i]-1-m[i]*pow_tau_mi)-m[i]*m[i]*pow_tau_mi);
			summer += n[i]*t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta-pow_delta_li-pow_tau_mi);
		}
		else if (l[i] != 0 && m[i] == 0){
			summer += n[i]*t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta-pow(delta,l[i]));
		}
		else
			summer += n[i]*t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta);
	}
	return summer;
}
double phir_Lemmon2005::dDelta_dTau2(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta), pow_delta_li, pow_tau_mi, bracket;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i] != 0 && m[i] != 0){
			pow_delta_li = pow(delta,l[i]);
			pow_tau_mi = pow(tau,m[i]);
			// delta derivative of second tau derivative
			bracket = ((t[i]-m[i]*pow_tau_mi)*(t[i]-1-m[i]*pow_tau_mi)-m[i]*m[i]*pow_tau_mi)*(d[i]-l[i]*pow_delta_li);
			summer+=n[i]*bracket*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta-pow_delta_li-pow_tau_mi);
		}
		else if (l[i] != 0 && m[i] == 0){
			pow_delta_li = pow(delta,l[i]);
			summer+=n[i]*t[i]*(t[i]-1)*(d[i]-l[i]*pow_delta_li)*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta-pow_delta_li);
		}
		else
			summer+=n[i]*t[i]*(t[i]-1)*d[i]*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta);
	}
	return summer;
}
double phir_Lemmon2005::dDelta(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta), pow_delta_li, pow_tau_mi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{	
		if (l[i] != 0 && m[i] != 0){
			pow_delta_li = pow(delta,l[i]);
			pow_tau_mi = pow(tau,m[i]);
			summer += n[i]*(d[i]-l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-1)*log_delta-pow_delta_li-pow_tau_mi);
		}
		else if (l[i]>0 && m[i] == 0){
			pow_delta_li = pow(delta,l[i]);
			summer += n[i]*(d[i]-l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-1)*log_delta-pow_delta_li);
		}
		else
			summer += n[i]*d[i]*exp(t[i]*log_tau+(d[i]-1)*log_delta);
	}
	return summer;
}
double phir_Lemmon2005::dDelta2(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta), pow_delta_li, pow_tau_mi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i] != 0 && m[i] != 0){	
			pow_delta_li = pow(delta,l[i]);
			pow_tau_mi = pow(tau,m[i]);
			summer+=n[i]*((d[i]-l[i]*pow_delta_li)*(d[i]-1.0-l[i]*pow_delta_li) - l[i]*l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-2)*log_delta-pow_delta_li-pow_tau_mi);
		}
		
		else if (l[i] != 0 && m[i] == 0){
			pow_delta_li = pow(delta,l[i]);
			summer+=n[i]*((d[i]-l[i]*pow_delta_li)*(d[i]-1.0-l[i]*pow_delta_li) - l[i]*l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-2)*log_delta-pow_delta_li);
		}
		else
			summer+=n[i]*d[i]*(d[i]-1.0)*exp(t[i]*log_tau+(d[i]-2)*log_delta);
	}
	return summer;
}
double phir_Lemmon2005::dDelta3(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta), pow_delta_li, pow_tau_mi, bracket;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i] != 0 && m[i] != 0){
			pow_delta_li = pow(delta,l[i]);
			pow_tau_mi = pow(tau,m[i]);
			bracket = (d[i]*(d[i]-1)*(d[i]-2)+pow_delta_li*(-2*l[i]+6*d[i]*l[i]-3*d[i]*d[i]*l[i]-3*d[i]*l[i]*l[i]+3*l[i]*l[i]-l[i]*l[i]*l[i])+pow_delta_li*pow_delta_li*(3*d[i]*l[i]*l[i]-3*l[i]*l[i]+3*l[i]*l[i]*l[i])-l[i]*l[i]*l[i]*pow_delta_li*pow_delta_li*pow_delta_li);
			summer+=n[i]*bracket*exp(t[i]*log_tau+(d[i]-3)*log_delta-pow_delta_li-pow_tau_mi);
		}
		else if (l[i] != 0 && m[i] == 0)
		{
			pow_delta_li = pow(delta,l[i]);
			bracket = (d[i]*(d[i]-1)*(d[i]-2)+pow_delta_li*(-2*l[i]+6*d[i]*l[i]-3*d[i]*d[i]*l[i]-3*d[i]*l[i]*l[i]+3*l[i]*l[i]-l[i]*l[i]*l[i])+pow_delta_li*pow_delta_li*(3*d[i]*l[i]*l[i]-3*l[i]*l[i]+3*l[i]*l[i]*l[i])-l[i]*l[i]*l[i]*pow_delta_li*pow_delta_li*pow_delta_li);
			summer+=n[i]*bracket*exp(t[i]*log_tau+(d[i]-3)*log_delta-pow_delta_li);
		}
		else
			summer+=n[i]*d[i]*(d[i]-1.0)*(d[i]-2)*exp(t[i]*log_tau+(d[i]-3)*log_delta);
	}
	return summer;
}

double phir_Lemmon2005::dDelta2_dTau(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta), bracket, pow_tau_mi, pow_delta_li;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i] != 0 && m[i] != 0){
			pow_delta_li = pow(delta,l[i]);
			pow_tau_mi = pow(tau,m[i]);
			bracket = (t[i]-m[i]*pow_tau_mi)*(((d[i]-l[i]*pow_delta_li))*(d[i]-1-l[i]*pow_delta_li)-l[i]*l[i]*pow_delta_li);
			summer += n[i]*bracket*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta-pow_delta_li-pow_tau_mi);
		}
		else if (l[i] != 0 && m[i] == 0){
			pow_delta_li = pow(delta,l[i]);
			bracket = t[i]*(((d[i]-l[i]*pow_delta_li))*(d[i]-1-l[i]*pow_delta_li)-l[i]*l[i]*pow_delta_li);
			summer += n[i]*bracket*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta-pow_delta_li);
		}
		else
			summer += n[i]*d[i]*t[i]*(d[i]-1)*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta);
	}
	return summer;
}
double phir_Lemmon2005::dDelta_dTau(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta), pow_delta_li, pow_tau_mi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i] != 0 && m[i] != 0){
			pow_delta_li = pow(delta,l[i]);
			pow_tau_mi = pow(tau,m[i]);
			summer+=n[i]*(d[i]-l[i]*pow_delta_li)*(t[i]-m[i]*pow_tau_mi)*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta-pow_delta_li-pow_tau_mi);
		}
		else if (l[i] != 0 && m[i] == 0){
			double pow_delta_li = pow(delta,l[i]);
			summer+=n[i]*t[i]*(d[i]-l[i]*pow_delta_li)*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta-pow_delta_li);
		}
		else
			summer+=n[i]*d[i]*t[i]*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta);
	}
	return summer;
}

phir_gaussian::phir_gaussian(std::vector<double> n_in, std::vector<double> d_in, std::vector<double> t_in, 
							 std::vector<double> alpha_in, std::vector<double> epsilon_in, std::vector<double> beta_in, std::vector<double> gamma_in,
	unsigned int iStart_in, unsigned int iEnd_in)
{
	n=n_in;
	d=d_in;
	t=t_in;
	alpha=alpha_in;
	epsilon=epsilon_in;
	beta=beta_in;
	gamma=gamma_in;
	iStart=iStart_in;
	iEnd=iEnd_in;
}

phir_gaussian::phir_gaussian(double n_in[], double d_in[],double t_in[], double alpha_in[], 
							 double epsilon_in[], double beta_in[], double gamma_in[],
							 unsigned int iStart_in, unsigned int iEnd_in, unsigned int N)
{
	n=std::vector<double>(n_in,n_in+N);
	d=std::vector<double>(d_in,d_in+N);
	t=std::vector<double>(t_in,t_in+N);
	alpha=std::vector<double>(alpha_in,alpha_in+N);
	epsilon=std::vector<double>(epsilon_in,epsilon_in+N);
	beta=std::vector<double>(beta_in,beta_in+N);
	gamma=std::vector<double>(gamma_in,gamma_in+N);
	iStart=iStart_in;
	iEnd=iEnd_in;
}
phir_gaussian::phir_gaussian(const double n_in[], const double d_in[], const double t_in[], const double alpha_in[], 
							 const double epsilon_in[], const double beta_in[], const double gamma_in[],
							 unsigned int iStart_in, unsigned int iEnd_in, unsigned int N)
{
	n=std::vector<double>(n_in,n_in+N);
	d=std::vector<double>(d_in,d_in+N);
	t=std::vector<double>(t_in,t_in+N);
	alpha=std::vector<double>(alpha_in,alpha_in+N);
	epsilon=std::vector<double>(epsilon_in,epsilon_in+N);
	beta=std::vector<double>(beta_in,beta_in+N);
	gamma=std::vector<double>(gamma_in,gamma_in+N);
	iStart=iStart_in;
	iEnd=iEnd_in;
}

// Term and its derivatives
double phir_gaussian::base(double tau, double delta) throw()
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(delta,d[i])*pow(tau,t[i])*psi;
	}
	return summer;
}
double phir_gaussian::dTau(double tau, double delta) throw()
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(delta,d[i])*pow(tau,t[i])*psi*(t[i]/tau-2.0*beta[i]*(tau-gamma[i]));
	}
	return summer;
}
double phir_gaussian::dTau2(double tau, double delta) throw()
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(delta,d[i])*pow(tau,t[i])*psi*(pow(t[i]/tau-2.0*beta[i]*(tau-gamma[i]),2)-t[i]/pow(tau,2)-2.0*beta[i]);
	}
	return summer;
}
double phir_gaussian::dTau3(double tau, double delta) throw()
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		// triple derivative product rule (a*b*c)' = a'*b*c+a*b'*c+a*b*c'
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		double dpsi_dTau = exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2))*(-2*beta[i]*(tau-gamma[i]));

		double bracket = pow(t[i]/tau-2.0*beta[i]*(tau-gamma[i]),2)-t[i]/pow(tau,2)-2.0*beta[i];
		double dbracket_dTau = 2*(t[i]/tau-2.0*beta[i]*(tau-gamma[i]))*(-t[i]/tau/tau-2*beta[i])+2*t[i]/pow(tau,3);
		summer+=n[i]*pow(delta,d[i])*(t[i]*pow(tau,t[i]-1)*psi*bracket+pow(tau,t[i])*dpsi_dTau*bracket+pow(tau,t[i])*psi*dbracket_dTau);
	}
	return summer;
}
double phir_gaussian::dDelta(double tau, double delta) throw()
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(delta,d[i])*pow(tau,t[i])*psi*(d[i]/delta-2.0*alpha[i]*(delta-epsilon[i]));
	}
	return summer;
}
double phir_gaussian::dDelta2(double tau, double delta) throw()
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(tau,t[i])*psi*(-2.0*alpha[i]*pow(delta,d[i])+4.0*pow(alpha[i],2)*pow(delta,d[i])*pow(delta-epsilon[i],2)-4.0*d[i]*alpha[i]*pow(delta,d[i]-1)*(delta-epsilon[i])+d[i]*(d[i]-1.0)*pow(delta,d[i]-2));
	}
	return summer;
}
double phir_gaussian::dDelta3(double tau, double delta) throw()
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		double bracket = (pow(d[i]-2*alpha[i]*delta*(delta-epsilon[i]),3)-3*d[i]*d[i]+2*d[i]-6*d[i]*alpha[i]*delta*delta+6*alpha[i]*delta*(delta-epsilon[i])*(d[i]+2*alpha[i]*delta*delta));
		summer+=n[i]*pow(tau,t[i])*pow(delta,d[i]-3)*psi*bracket;
	}
	return summer;
}
double phir_gaussian::dDelta2_dTau(double tau, double delta) throw()
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(tau,t[i])*psi*(t[i]/tau-2.0*beta[i]*(tau-gamma[i]))*(-2.0*alpha[i]*pow(delta,d[i])+4.0*pow(alpha[i],2)*pow(delta,d[i])*pow(delta-epsilon[i],2)-4.0*d[i]*alpha[i]*pow(delta,d[i]-1)*(delta-epsilon[i])+d[i]*(d[i]-1.0)*pow(delta,d[i]-2));
	}
	return summer;
}
double phir_gaussian::dDelta_dTau(double tau, double delta) throw()
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(delta,d[i])*pow(tau,t[i])*psi*(d[i]/delta-2.0*alpha[i]*(delta-epsilon[i]))*(t[i]/tau-2.0*beta[i]*(tau-gamma[i]));
	}
	return summer;
}
double phir_gaussian::dDelta_dTau2(double tau, double delta) throw()
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(delta,d[i])*pow(tau,t[i])*psi*(d[i]/delta-2.0*alpha[i]*(delta-epsilon[i]))*(pow(t[i]-2.0*beta[i]*tau*(tau-gamma[i]),2)-t[i]-2*beta[i]*tau*tau)/tau/tau;
	}
	return summer;
}

phir_GERG2008_gaussian::phir_GERG2008_gaussian(std::vector<double> n_in, std::vector<double> d_in, std::vector<double> t_in, 
											   std::vector<double> eta_in, std::vector<double> epsilon_in, std::vector<double> beta_in, std::vector<double> gamma_in,
	unsigned int iStart_in, unsigned int iEnd_in)
{
	n=n_in;
	d=d_in;
	t=t_in;
	eta=eta_in;
	epsilon=epsilon_in;
	beta=beta_in;
	gamma=gamma_in;
	iStart=iStart_in;
	iEnd=iEnd_in;
}

phir_GERG2008_gaussian::phir_GERG2008_gaussian(double n_in[], double d_in[],double t_in[], double eta_in[], 
							 double epsilon_in[], double beta_in[], double gamma_in[],
							 unsigned int iStart_in, unsigned int iEnd_in, unsigned int N)
{
	n=std::vector<double>(n_in, n_in+N);
	d=std::vector<double>(d_in, d_in+N);
	t=std::vector<double>(t_in, t_in+N);
	eta=std::vector<double>(eta_in, eta_in+N);
	epsilon=std::vector<double>(epsilon_in, epsilon_in+N);
	beta=std::vector<double>(beta_in, beta_in+N);
	gamma=std::vector<double>(gamma_in, gamma_in+N);
	iStart=iStart_in;
	iEnd=iEnd_in;
}
double phir_GERG2008_gaussian::base(double tau, double delta)
{
	double summer = 0, psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi = exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*(delta-gamma[i]));
		summer += n[i]*pow(tau,t[i])*pow(delta,d[i])*psi;
	}
	return summer;
}
double phir_GERG2008_gaussian::dDelta(double tau, double delta)
{
	double summer = 0, psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi = exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*(delta-gamma[i]));
		summer += n[i]*pow(tau,t[i])*pow(delta,d[i])*psi*(d[i]/delta-2*eta[i]*(delta-epsilon[i])-beta[i]);
	}
	return summer;
}
double phir_GERG2008_gaussian::dDelta2(double tau, double delta)
{
	double summer = 0, psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi = exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*(delta-gamma[i]));
		summer += n[i]*pow(tau,t[i])*pow(delta,d[i])*psi*(pow(d[i]/delta-2*eta[i]*(delta-epsilon[i])-beta[i],2)-d[i]/delta/delta-2*eta[i]);
	}
	return summer;
}
double phir_GERG2008_gaussian::dTau(double tau, double delta)
{
	double summer = 0, psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi = exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*(delta-gamma[i]));
		summer += n[i]*t[i]*pow(tau,t[i]-1)*pow(delta,d[i])*psi;
	}
	return summer;
}
double phir_GERG2008_gaussian::dTau2(double tau, double delta)
{
	double summer = 0, psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi = exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*(delta-gamma[i]));
		summer += n[i]*t[i]*(t[i]-1)*pow(tau,t[i]-2)*pow(delta,d[i])*psi;
	}
	return summer;
}
double phir_GERG2008_gaussian::dDelta_dTau(double tau, double delta)
{
	double summer = 0, psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi = exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*(delta-gamma[i]));
		summer += n[i]*t[i]*pow(tau,t[i]-1)*pow(delta,d[i])*psi*(d[i]/delta-2*eta[i]*(delta-epsilon[i])-beta[i]);
	}
	return summer;
}


phir_critical::phir_critical(std::vector<double> n_in, std::vector<double> d_in, std::vector<double> t_in, 
		std::vector<double> a_in, std::vector<double> b_in, std::vector<double> beta_in,
		std::vector<double> A_in, std::vector<double> B_in, std::vector<double> C_in, 
		std::vector<double> D_in, int iStart_in, int iEnd_in)
{
	n=n_in;
	d=d_in;
	t=t_in;
	a=a_in;
	b=b_in;
	beta=beta_in;
	A=A_in;
	B=B_in;
	C=C_in;
	D=D_in;
	iStart=iStart_in;
	iEnd=iEnd_in;
}

phir_critical::phir_critical(double n[], double d[], double t[], 
							 double a[], double b[], double beta[],
							 double A[], double B[], double C[], 
							 double D[], int iStart, int iEnd, 
							 int N)
{
	this->n=std::vector<double>(n,n+N);
	this->d=std::vector<double>(d,d+N);
	this->t=std::vector<double>(t,t+N);
	this->a=std::vector<double>(a,a+N);
	this->b=std::vector<double>(b,b+N);
	this->beta=std::vector<double>(beta,beta+N);
	this->A=std::vector<double>(A,A+N);
	this->B=std::vector<double>(B,B+N);
	this->C=std::vector<double>(C,C+N);
	this->D=std::vector<double>(D,D+N);
	this->iStart=iStart;
	this->iEnd=iEnd;
}

double phir_critical::base(double tau, double delta) throw()
{
	double summer=0,theta,DELTA,PSI;
	for (int i=iStart;i<=iEnd;i++)
	{
		theta=(1.0-tau)+A[i]*pow(pow(delta-1.0,2),1/(2*beta[i]));
        DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        summer+=n[i]*pow(DELTA,b[i])*delta*PSI;
	}
	return summer;
}

double phir_critical::dDelta(double tau, double delta) throw()
{
	double summer=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta;
	for (int i=iStart;i<=iEnd;i++)
	{
		theta=(1.0-tau)+A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(pow(delta-1.0,2),a[i]-1.0));
		
		// At critical point, DELTA is 0, and 1/0^n is undefined
		if (fabs(DELTA) < 10*DBL_EPSILON)
		{
			dDELTAbi_dDelta = 0;
		}
		else{
			dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
		}
        summer+=n[i]*(pow(DELTA,b[i])*(PSI+delta*dPSI_dDelta)+dDELTAbi_dDelta*delta*PSI);
	}
	return summer;
}
double phir_critical::dDelta_dTau2(double tau, double delta) throw()
{
	double summer=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta;
	double dDELTAbi_dTau,dPSI_dTau, dtheta_dDelta;
	double dPSI2_dDelta2, dDELTA2_dDelta2,dDELTAbi2_dDelta2,dPSI2_dTau2,dDELTAbi2_dTau2;
	double dDELTAbi2_dDelta_dTau,dPSI2_dDelta_dTau;
	double dDELTAbi3_dDelta_dTau2,dPSI3_dDelta_dTau2;

	for (int i=iStart;i<=iEnd;i++)
	{
		theta=(1.0-tau)+A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(pow(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
        dPSI2_dDelta2=(2.0*C[i]*pow(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
        dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(pow(delta-1.0,2),a[i]-2.0)+2.0*pow(A[i]/beta[i],2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
        dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*pow(dDELTA_dDelta,2));
        
		dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);

		dPSI2_dTau2=(2.0*D[i]*pow(tau-1.0,2)-1.0)*2.0*D[i]*PSI;
        dDELTAbi2_dTau2=2.0*b[i]*pow(DELTA,b[i]-1.0)+4.0*pow(theta,2)*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0);
		dPSI2_dDelta_dTau=4.0*C[i]*D[i]*(delta-1.0)*(tau-1.0)*PSI;
        dDELTAbi2_dDelta_dTau=-A[i]*b[i]*2.0/beta[i]*pow(DELTA,b[i]-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)-2.0*theta*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0)*dDELTA_dDelta;

		dPSI3_dDelta_dTau2 = 2*D[i]*(2*D[i]*pow(tau-1,2)-1)*dPSI_dDelta;
		dtheta_dDelta = A[i]/(2*beta[i])*pow(pow(delta-1,2),1/(2*beta[i])-1)*2*(delta-1);
		dDELTAbi3_dDelta_dTau2 = 2*b[i]*(b[i]-1)*pow(DELTA,b[i]-2)*dDELTA_dDelta+4*pow(theta,2)*b[i]*(b[i]-1)*(b[i]-2)*pow(DELTA,b[i]-3)*dDELTA_dDelta+8*theta*b[i]*(b[i]-1)*pow(DELTA,b[i]-2)*dtheta_dDelta;
		
        summer += n[i]*delta*(dDELTAbi2_dTau2*dPSI_dDelta+dDELTAbi3_dDelta_dTau2*PSI+2*dDELTAbi_dTau*dPSI2_dDelta_dTau+2.0*dDELTAbi2_dDelta_dTau*dPSI_dTau+pow(DELTA,b[i])*dPSI3_dDelta_dTau2+dDELTAbi_dDelta*dPSI2_dTau2)+n[i]*(dDELTAbi2_dTau2*PSI+2.0*dDELTAbi_dTau*dPSI_dTau+pow(DELTA,b[i])*dPSI2_dTau2);
	}
	return summer;
}

double phir_critical::dDelta2(double tau, double delta) throw()
{
	double summer=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta;
	double dPSI2_dDelta2, dDELTA2_dDelta2,dDELTAbi2_dDelta2;
	for (int i=iStart;i<=iEnd;i++)
	{
		theta=(1.0-tau)+A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(pow(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
        dPSI2_dDelta2=(2.0*C[i]*pow(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
        dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(pow(delta-1.0,2),a[i]-2.0)+2.0*pow(A[i]/beta[i],2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
        dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*pow(dDELTA_dDelta,2));
        
        summer+=n[i]*(pow(DELTA,b[i])*(2.0*dPSI_dDelta+delta*dPSI2_dDelta2)+2.0*dDELTAbi_dDelta*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta2*delta*PSI);
	}
	return summer;
}
double phir_critical::dDelta3(double tau, double delta) throw()
{
	double summer=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta;
	double dPSI2_dDelta2, dDELTA2_dDelta2,dDELTAbi2_dDelta2;
	double dPSI3_dDelta3,PI,dPI_dDelta,dDELTA3_dDelta3, dDELTAbi3_dDelta3;
	double dtheta_dDelta;
	for (int i=iStart;i<=iEnd;i++)
	{
		theta=(1.0-tau)+A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(pow(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
        dPSI2_dDelta2=(2.0*C[i]*pow(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
        dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(pow(delta-1.0,2),a[i]-2.0)+2.0*pow(A[i]/beta[i],2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
        dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*pow(dDELTA_dDelta,2));

		dPSI3_dDelta3=2.0*C[i]*PSI*(-4*C[i]*C[i]*pow(delta-1.0,3)+6*C[i]*(delta-1));
		dtheta_dDelta = A[i]/(2*beta[i])*pow(pow(delta-1,2),1/(2*beta[i])-1)*2*(delta-1);

		PI = 4*B[i]*a[i]*(a[i]-1)*pow(pow(delta-1,2),a[i]-2)+2*pow(A[i]/beta[i],2)*pow(pow(delta-1,2),1/beta[i]-2)+4*A[i]*theta/beta[i]*(1/(2*beta[i])-1)*pow(pow(delta-1,2),1/(2*beta[i])-2);
		dPI_dDelta = -8*B[i]*a[i]*(a[i]-1)*(a[i]-2)*pow(pow(delta-1,2),a[i]-5.0/2.0)-8*pow(A[i]/beta[i],2)*(1/(2*beta[i])-1)*pow(pow(delta-1,2),1/beta[i]-5.0/2.0)-(8*A[i]*theta)/beta[i]*(1/(2*beta[i])-1)*(1/(2*beta[i])-2)*pow(pow(delta-1,2),1/(2*beta[i])-5.0/2.0)+4*A[i]/beta[i]*(1/(2*beta[i])-1)*pow(pow(delta-1,2),1/(2*beta[i])-2)*dtheta_dDelta;
		dDELTA3_dDelta3 = 1/(delta-1)*dDELTA2_dDelta2-1/pow(delta-1,2)*dDELTA_dDelta+pow(delta-1,2)*dPI_dDelta+2*(delta-1)*PI;
		dDELTAbi3_dDelta3 = b[i]*(pow(DELTA,b[i]-1)*dDELTA3_dDelta3+dDELTA2_dDelta2*(b[i]-1)*pow(DELTA,b[i]-2)*dDELTA_dDelta+(b[i]-1)*(pow(DELTA,b[i]-2)*2*dDELTA_dDelta*dDELTA2_dDelta2+pow(dDELTA_dDelta,2)*(b[i]-2)*pow(DELTA,b[i]-3)*dDELTA_dDelta));
        
		summer += n[i]*(pow(DELTA,b[i])*(3.0*dPSI2_dDelta2+delta*dPSI3_dDelta3)+3.0*dDELTAbi_dDelta*(2*dPSI_dDelta+delta*dPSI2_dDelta2)+3*dDELTAbi2_dDelta2*(PSI+delta*dPSI_dDelta)+dDELTAbi3_dDelta3*PSI*delta);
	}
	return summer;
}

double phir_critical::dDelta2_dTau(double tau, double delta) throw()
{
	double summer=0;
	double theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta,dPSI2_dDelta2,dDELTAbi2_dDelta2,dDELTA2_dDelta2;
    double dPSI2_dDelta_dTau, dDELTAbi2_dDelta_dTau, dPSI_dTau, dDELTAbi_dTau;
    double Line1,Line2,Line3,dDELTA2_dDelta_dTau,dDELTA3_dDelta2_dTau,dDELTAbim1_dTau,dDELTAbim2_dTau;
    double dDELTA_dTau,dDELTAbi3_dDelta2_dTau,dPSI3_dDelta2_dTau;

	for (int i=iStart;i<=iEnd;i++)
	{
		theta=(1.0-tau)+A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(pow(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
        dPSI2_dDelta2=(2.0*C[i]*pow(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
        dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(pow(delta-1.0,2),a[i]-2.0)+2.0*pow(A[i]/beta[i],2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
        dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*pow(dDELTA_dDelta,2));
        
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        
        dPSI2_dDelta_dTau=4.0*C[i]*D[i]*(delta-1.0)*(tau-1.0)*PSI;
        dDELTAbi2_dDelta_dTau=-A[i]*b[i]*2.0/beta[i]*pow(DELTA,b[i]-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)-2.0*theta*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0)*dDELTA_dDelta;
	
		//Following Terms added for this derivative
		dPSI3_dDelta2_dTau = (2.0*C[i]*pow(delta-1.0,2)-1.0)*2.0*C[i]*dPSI_dTau;
		dDELTA_dTau = -2*theta;
		dDELTA2_dDelta_dTau = 2.0*A[i]/(beta[i])*pow(pow(delta-1,2),1.0/(2.0*beta[i])-0.5);
		dDELTA3_dDelta2_dTau = 2.0*A[i]*(beta[i]-1)/(beta[i]*beta[i])*pow(pow(delta-1,2),1/(2*beta[i])-1.0);
		
		dDELTAbim1_dTau = (b[i]-1)*pow(DELTA,b[i]-2)*dDELTA_dTau;
		dDELTAbim2_dTau = (b[i]-2)*pow(DELTA,b[i]-3)*dDELTA_dTau;
		Line1 = dDELTAbim1_dTau*dDELTA2_dDelta2 + pow(DELTA,b[i]-1)*dDELTA3_dDelta2_dTau;
		Line2 = (b[i]-1)*(dDELTAbim2_dTau*pow(dDELTA_dDelta,2)+pow(DELTA,b[i]-2)*2*dDELTA_dDelta*dDELTA2_dDelta_dTau);
		dDELTAbi3_dDelta2_dTau = b[i]*(Line1+Line2);
		
		Line1 = pow(DELTA,b[i])*(2*dPSI2_dDelta_dTau+delta*dPSI3_dDelta2_dTau)+dDELTAbi_dTau*(2*dPSI_dDelta+delta*dPSI2_dDelta2);
		Line2 = 2*dDELTAbi_dDelta*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+2*dDELTAbi2_dDelta_dTau*(PSI+delta*dPSI_dDelta);
		Line3 = dDELTAbi2_dDelta2*delta*dPSI_dTau + dDELTAbi3_dDelta2_dTau*delta*PSI;
        summer += n[i]*(Line1+Line2+Line3);
    }
	return summer;
}

double phir_critical::dDelta_dTau(double tau, double delta) throw()
{
	double summer=0,theta,DELTA,PSI,dPSI_dDelta,dDELTAbi_dDelta;
	double dPSI_dTau, dDELTAbi_dTau,dDELTA_dDelta, dPSI2_dDelta_dTau;
	double dDELTAbi2_dDelta_dTau;
	for (int i=iStart;i<=iEnd;i++)
	{
		theta=(1.0-tau)+A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(pow(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        
        dPSI2_dDelta_dTau=4.0*C[i]*D[i]*(delta-1.0)*(tau-1.0)*PSI;
        dDELTAbi2_dDelta_dTau=-A[i]*b[i]*2.0/beta[i]*pow(DELTA,b[i]-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)-2.0*theta*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0)*dDELTA_dDelta;
        
        summer+=n[i]*(pow(DELTA,b[i])*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+delta*dDELTAbi_dDelta*dPSI_dTau+ dDELTAbi_dTau*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta_dTau*delta*PSI);
	}
	return summer;
}
	
double phir_critical::dTau(double tau, double delta) throw()
{
	double summer=0,theta,DELTA,PSI,dPSI_dTau, dDELTAbi_dTau;

	for (int i=iStart;i<=iEnd;i++)
	{
		theta=(1.0-tau)+A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
		if (fabs(DELTA)<10*DBL_EPSILON)
			dDELTAbi_dTau = 0;
		else
			dDELTAbi_dTau = -2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        summer+=n[i]*delta*(dDELTAbi_dTau*PSI+pow(DELTA,b[i])*dPSI_dTau);
	}
	return summer;
}

double phir_critical::dTau2(double tau, double delta) throw()
{
	double summer=0,theta,DELTA,PSI,dPSI_dTau, dDELTAbi_dTau;
	double dPSI2_dTau2, dDELTAbi2_dTau2;
	for (int i=iStart;i<=iEnd;i++)
	{
		theta=(1.0-tau)+A[i]*pow(pow(delta-1.0,2),1/(2*beta[i]));
        DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        dPSI2_dTau2=(2.0*D[i]*pow(tau-1.0,2)-1.0)*2.0*D[i]*PSI;
        dDELTAbi2_dTau2=2.0*b[i]*pow(DELTA,b[i]-1.0)+4.0*pow(theta,2)*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0);
        summer+=n[i]*delta*(dDELTAbi2_dTau2*PSI+2.0*dDELTAbi_dTau*dPSI_dTau+pow(DELTA,b[i])*dPSI2_dTau2);
	}
	return summer;
}
double phir_critical::dTau3(double tau, double delta) throw()
{
	double summer=0,theta,DELTA,PSI,dPSI_dTau, dDELTAbi_dTau;
	double dPSI2_dTau2, dDELTAbi2_dTau2, dPSI3_dTau3, dDELTAbi3_dTau3;
	for (int i=iStart;i<=iEnd;i++)
	{
		theta=(1.0-tau)+A[i]*pow(pow(delta-1.0,2),1/(2*beta[i]));
        DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        dPSI2_dTau2=(2.0*D[i]*pow(tau-1.0,2)-1.0)*2.0*D[i]*PSI;
        dDELTAbi2_dTau2=2.0*b[i]*pow(DELTA,b[i]-1.0)+4.0*pow(theta,2)*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0);
		
		dPSI3_dTau3=2.0*D[i]*PSI*(-4*D[i]*D[i]*pow(tau-1,3)+6*D[i]*(tau-1));
		dDELTAbi3_dTau3 = -12.0*theta*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2)-8*pow(theta,3)*b[i]*(b[i]-1)*(b[i]-2)*pow(DELTA,b[i]-3.0);
        summer += n[i]*delta*(dDELTAbi3_dTau3*PSI+(3.0*dDELTAbi2_dTau2)*dPSI_dTau+(3*dDELTAbi_dTau )*dPSI2_dTau2+pow(DELTA,b[i])*dPSI3_dTau3);
	}
	return summer;
}
#ifndef DISABLE_CATCH
	TEST_CASE("Non-analytic critical point Helmholtz derivative check", "[helmholtz],[fast]")
	{
		// From CO2
		double n[] = {0,-0.666422765408E+00,0.726086323499E+00,0.550686686128E-01};
		double d[] = {0,2,3,3};
		double t[] = {0, 1.00, 3.00, 3.00};
		double a[] = {0, 3.5, 3.5, 3.0};
		double b[] = {0, 0.875, 0.925, 0.875};
		double beta[] = {9,0.300, 0.300, 0.300};
		double A[] = {0, 0.700, 0.700, 0.700};
		double B[] = {0, 0.3, 0.3, 1.0};
		double C[] = {0, 10.0, 10.0, 12.5};
		double D[] = {0, 275.0, 275.0, 275.0};
		
		phir_critical phir = phir_critical(n,d,t,a,b,beta,A,B,C,D,1,3,4);
		double eps = sqrt(DBL_EPSILON);

		SECTION("dDelta")
		{
			double ANA = phir.dDelta(0.5, 0.5);
			double NUM = (phir.base(0.5, 0.5+eps) - phir.base(0.5,0.5-eps))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
		SECTION("dTau")
		{
			double ANA = phir.dTau(0.5, 0.5);
			double NUM = (phir.base(0.5+eps, 0.5) - phir.base(0.5-eps,0.5))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
		SECTION("dDelta2")
		{
			double ANA = phir.dDelta2(0.5, 0.5);
			double NUM = (phir.dDelta(0.5, 0.5+eps) - phir.dDelta(0.5,0.5-eps))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
		SECTION("dTau2")
		{
			double ANA = phir.dTau2(0.5, 0.5);
			double NUM = (phir.dTau(0.5+eps, 0.5) - phir.dTau(0.5-eps,0.5))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
		SECTION("dDeltadTau")
		{
			double ANA = phir.dDelta_dTau(0.5, 0.5);
			double NUM = (phir.dTau(0.5, 0.5+eps) - phir.dTau(0.5,0.5-eps))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
	}
#endif

double phir_SAFT_associating::Deltabar(double tau, double delta)
{
	return this->g(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar;
}   
double phir_SAFT_associating::dDeltabar_ddelta__consttau(double tau, double delta)
{
    return this->dg_deta(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar*this->vbarn;
}
double phir_SAFT_associating::d2Deltabar_ddelta2__consttau(double tau, double delta)
{
    return this->d2g_deta2(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar*pow(this->vbarn,(int)2);
}
double phir_SAFT_associating::dDeltabar_dtau__constdelta(double tau, double delta)
{
    return this->g(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*this->epsilonbar;
}
double phir_SAFT_associating::d2Deltabar_dtau2__constdelta(double tau, double delta)
{
    return this->g(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*pow(this->epsilonbar,(int)2);
}
double phir_SAFT_associating::d2Deltabar_ddelta_dtau(double tau, double delta)
{
    return this->dg_deta(this->eta(delta))*exp(this->epsilonbar*tau)*this->epsilonbar*this->kappabar*this->vbarn;
}
double phir_SAFT_associating::X(double delta, double Deltabar)
{
	return 2/(sqrt(1+4*Deltabar*delta)+1);
}
double phir_SAFT_associating::dX_dDeltabar__constdelta(double delta, double Deltabar)
{
    double X = this->X(delta,Deltabar);
    return -delta*X*X/(2*Deltabar*delta*X+1);
}
double phir_SAFT_associating::dX_ddelta__constDeltabar(double delta, double Deltabar)
{
    double X = this->X(delta,Deltabar);
    return -Deltabar*X*X/(2*Deltabar*delta*X+1);
}
double phir_SAFT_associating::dX_dtau(double tau, double delta)
{
    double Deltabar = this->Deltabar(tau, delta);
    return this->dX_dDeltabar__constdelta(delta, Deltabar)*this->dDeltabar_dtau__constdelta(tau, delta);
}
double phir_SAFT_associating::dX_ddelta(double tau, double delta)
{
    double Deltabar = this->Deltabar(tau, delta);
    return (this->dX_ddelta__constDeltabar(delta, Deltabar)
           + this->dX_dDeltabar__constdelta(delta, Deltabar)*this->dDeltabar_ddelta__consttau(tau, delta));
}
double phir_SAFT_associating::d2X_dtau2(double tau, double delta)
{
    double Deltabar = this->Deltabar(tau, delta);
    double X = this->X(delta, Deltabar);
    double beta = this->dDeltabar_dtau__constdelta(tau, delta);
    double d_dXdtau_dbeta = -delta*X*X/(2*Deltabar*delta*X+1);
    double d_dXdtau_dDeltabar = 2*delta*delta*X*X*X/pow(2*Deltabar*delta*X+1,2)*beta;
    double d_dXdtau_dX = -2*beta*delta*X*(Deltabar*delta*X+1)/pow(2*Deltabar*delta*X+1,2);
    double dbeta_dtau = this->d2Deltabar_dtau2__constdelta(tau, delta);
    double dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    return d_dXdtau_dX*dX_dDeltabar*beta+d_dXdtau_dDeltabar*beta+d_dXdtau_dbeta*dbeta_dtau;
}
double phir_SAFT_associating::d2X_ddeltadtau(double tau, double delta)
{
    double Deltabar = this->Deltabar(tau, delta);
    double X = this->X(delta, Deltabar);
    double alpha = this->dDeltabar_ddelta__consttau(tau, delta);
    double beta = this->dDeltabar_dtau__constdelta(tau, delta);
    double dalpha_dtau = this->d2Deltabar_ddelta_dtau(tau, delta);
    double d_dXddelta_dDeltabar = X*X*(2*delta*delta*X*alpha-1)/pow(2*Deltabar*delta*X+1,2);
    double d_dXddelta_dalpha = -delta*X*X/(2*Deltabar*delta*X+1);
    double d_dXddelta_dX = -(Deltabar+delta*alpha)*2*(Deltabar*delta*X*X+X)/pow(2*Deltabar*delta*X+1,2);
    double dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    return d_dXddelta_dX*dX_dDeltabar*beta+d_dXddelta_dDeltabar*beta+d_dXddelta_dalpha*dalpha_dtau;
}
double phir_SAFT_associating::d2X_ddelta2(double tau, double delta)
{
    double Deltabar = this->Deltabar(tau, delta);
    double X = this->X(delta, Deltabar);
    double alpha = this->dDeltabar_ddelta__consttau(tau, delta);
    double dalpha_ddelta = this->d2Deltabar_ddelta2__consttau(tau, delta);
    double d_dXddelta_dDeltabar = X*X*(2*delta*delta*X*alpha-1)/pow(2*Deltabar*delta*X+1,2);
    double d_dXddelta_dalpha = -delta*X*X/(2*Deltabar*delta*X+1);
    double d_dXddelta_dX = -(Deltabar+delta*alpha)*2*(Deltabar*delta*X*X+X)/pow(2*Deltabar*delta*X+1,2);
    double dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    double dX_ddelta_constall = X*X*(2*Deltabar*Deltabar*X-alpha)/pow(2*Deltabar*delta*X+1,2);
    return (dX_ddelta_constall
            + d_dXddelta_dX*this->dX_ddelta__constDeltabar(delta, Deltabar)
            + d_dXddelta_dX*dX_dDeltabar*alpha
            + d_dXddelta_dDeltabar*alpha
            + d_dXddelta_dalpha*dalpha_ddelta);
}   
double phir_SAFT_associating::g(double eta)
{
	return 0.5*(2-eta)/pow(1-eta,(int)3);
}    
double phir_SAFT_associating::dg_deta(double eta)
{
	return 0.5*(5-2*eta)/pow(1-eta,(int)4);
}
double phir_SAFT_associating::d2g_deta2(double eta)
{
    return 3*(3-eta)/pow(1-eta,(int)5);
}   
double phir_SAFT_associating::d3g_deta3(double eta)
{
	return 6*(7-2*eta)/pow(1-eta,(int)6);
}   
double phir_SAFT_associating::eta(double delta){
	return this->vbarn*delta;
}

double phir_SAFT_associating_1::base(double tau, double delta)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
    return this->m*((log(X)-X/2.0+0.5));
}
double phir_SAFT_associating_1::dDelta(double tau, double delta)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
    return this->m*(1/X-0.5)*this->dX_ddelta(tau, delta);
}
double phir_SAFT_associating_1::dTau(double tau, double delta)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
    return this->m*(1/X-0.5)*this->dX_dtau(tau, delta);
}
double phir_SAFT_associating_1::dTau2(double tau, double delta)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
	double X_tau = this->dX_dtau(tau, delta);
    double X_tautau = this->d2X_dtau2(tau, delta);
    return this->m*((1/X-0.5)*X_tautau-pow(X_tau/X, 2));
}
double phir_SAFT_associating_1::dDelta2(double tau, double delta)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
	double X_delta = this->dX_ddelta(tau, delta);
    double X_deltadelta = this->d2X_ddelta2(tau, delta);
    return this->m*((1/X-0.5)*X_deltadelta-pow(X_delta/X,2));
}
double phir_SAFT_associating_1::dDelta_dTau(double tau, double delta)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
	double X_delta = this->dX_ddelta(tau, delta);
    double X_deltadelta = this->d2X_ddelta2(tau, delta);
    double X_tau = this->dX_dtau(tau, delta);
    double X_deltatau = this->d2X_ddeltadtau(tau, delta);
    return this->m*((-X_tau/X/X)*X_delta+X_deltatau*(1/X-0.5));
}
#ifndef DISABLE_CATCH
	TEST_CASE("SAFT 1 Helmholtz derivative check", "[helmholtz],[fast]")
	{
		double m = 1.01871348;
		double vbarn = 0.444215309e-1;
		double kappabar = 0.109117041e-4;
		double epsilonbar = 12.2735737;
		phir_SAFT_associating_1 phir = phir_SAFT_associating_1(m,epsilonbar,vbarn,kappabar);
		double eps = sqrt(DBL_EPSILON);

		SECTION("dDelta")
		{
			double ANA = phir.dDelta(0.5, 0.5);
			double NUM = (phir.base(0.5, 0.5+eps) - phir.base(0.5,0.5-eps))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
		SECTION("dTau")
		{
			double ANA = phir.dTau(0.5, 0.5);
			double NUM = (phir.base(0.5+eps, 0.5) - phir.base(0.5-eps,0.5))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
		SECTION("dDelta2")
		{
			double ANA = phir.dDelta2(0.5, 0.5);
			double NUM = (phir.dDelta(0.5, 0.5+eps) - phir.dDelta(0.5,0.5-eps))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
		SECTION("dTau2")
		{
			double ANA = phir.dTau2(0.5, 0.5);
			double NUM = (phir.dTau(0.5+eps, 0.5) - phir.dTau(0.5-eps,0.5))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
		SECTION("dDeltadTau")
		{
			double ANA = phir.dDelta_dTau(0.5, 0.5);
			double NUM = (phir.dTau(0.5, 0.5+eps) - phir.dTau(0.5,0.5-eps))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
	}
#endif
double phir_SAFT_associating_2B::base(double tau, double delta)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
    return this->m*2*((log(X)-X/2.0+0.5));
}
double phir_SAFT_associating_2B::dDelta(double tau, double delta)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
    return this->m*2*(1/X-0.5)*this->dX_ddelta(tau, delta);
}
double phir_SAFT_associating_2B::dTau(double tau, double delta)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
    return this->m*2*(1/X-0.5)*this->dX_dtau(tau, delta);
}
double phir_SAFT_associating_2B::dTau2(double tau, double delta)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
	double X_tau = this->dX_dtau(tau, delta);
    double X_tautau = this->d2X_dtau2(tau, delta);
    return this->m*(2*(1/X-0.5)*X_tautau-2*pow(X_tau/X, 2));
}
double phir_SAFT_associating_2B::dDelta2(double tau, double delta)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
	double X_delta = this->dX_ddelta(tau, delta);
    double X_deltadelta = this->d2X_ddelta2(tau, delta);
    return this->m*(2*(1/X-0.5)*X_deltadelta-2*pow(X_delta/X,2));
}
double phir_SAFT_associating_2B::dDelta_dTau(double tau, double delta)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
	double X_delta = this->dX_ddelta(tau, delta);
    double X_deltadelta = this->d2X_ddelta2(tau, delta);
    double X_tau = this->dX_dtau(tau, delta);
    double X_deltatau = this->d2X_ddeltadtau(tau, delta);
    return this->m*(2*(-X_tau/X/X)*X_delta+2*X_deltatau*(1/X-0.5));
}
#ifndef DISABLE_CATCH

	TEST_CASE("SAFT 2B Helmholtz derivative check", "[helmholtz],[fast]")
	{
		double m = 0.977118832;
		double epsilon = 5.46341463;
		double vbarn = 0.204481952;
		double kappa = 0.148852832e-2;
		phir_SAFT_associating_2B phir = phir_SAFT_associating_2B(m,epsilon,vbarn,kappa);
		double eps = sqrt(DBL_EPSILON);

		SECTION("dDelta")
		{
			double ANA = phir.dDelta(0.5, 0.5);
			double NUM = (phir.base(0.5, 0.5+eps) - phir.base(0.5,0.5-eps))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
		SECTION("dTau")
		{
			double ANA = phir.dTau(0.5, 0.5);
			double NUM = (phir.base(0.5+eps, 0.5) - phir.base(0.5-eps,0.5))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
		SECTION("dDelta2")
		{
			double ANA = phir.dDelta2(0.5, 0.5);
			double NUM = (phir.dDelta(0.5, 0.5+eps) - phir.dDelta(0.5,0.5-eps))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
		SECTION("dTau2")
		{
			double ANA = phir.dTau2(0.5, 0.5);
			double NUM = (phir.dTau(0.5+eps, 0.5) - phir.dTau(0.5-eps,0.5))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
		SECTION("dDeltadTau")
		{
			double ANA = phir.dDelta_dTau(0.5, 0.5);
			double NUM = (phir.dTau(0.5, 0.5+eps) - phir.dTau(0.5,0.5-eps))/(2*eps);
			REQUIRE(abs(NUM-ANA) < 1e-6);
		}
	}
#endif

double phi0_Planck_Einstein::base(double tau, double delta)
{
	double summer=0;
	for (int i=iStart;i<=iEnd;i++)
	{
		summer+=a[i]*log(1.0-exp(-theta[i]*tau));
	}
	return summer;
}
double phi0_Planck_Einstein::dTau(double tau, double delta)
{
	double summer=0;
	for (int i=iStart;i<=iEnd;i++)
	{
		summer+=a[i]*theta[i]*(1.0/(1.0-exp(-theta[i]*tau))-1.0);
	}
	return summer;
}
double phi0_Planck_Einstein::dTau2(double tau, double delta)
{
	double summer=0;
	for (int i=iStart;i<=iEnd;i++)
	{
		summer -= a[i]*pow(theta[i],2.0)*exp(theta[i]*tau)/pow(1.0-exp(theta[i]*tau),2.0);
	}
	return summer;
}
double phi0_Planck_Einstein::dTau3(double tau, double delta)
{
	double summer=0;
	for (int i=iStart;i<=iEnd;i++)
	{
		summer += a[i]*pow(theta[i],2.0)*theta[i]*exp(theta[i]*tau)*(exp(theta[i]*tau)+1)/pow(exp(theta[i]*tau)-1,3.0);
	}
	return summer;
}

/*
Maxima code for the term:

term:a*log(c+exp(%theta*%tau));
ratsimp(diff(term,%tau));
ratsimp(diff(%,%tau));

(%o23) a*log(c+%e^(%tau*%theta))
(%o24) (%theta*%e^(%tau*%theta)*a)/(c+%e^(%tau*%theta))
(%o25) (%theta^2*%e^(%tau*%theta)*a*c)/(c^2+2*%e^(%tau*%theta)*c+%e^(2*%tau*%theta))
*/
double phi0_Planck_Einstein2::base(double tau, double delta)
{
	double summer=0;
	for (int i=iStart;i<=iEnd;i++)
	{
		//a_0*log(c+exp(-theta_0*tau))
		summer+=a[i]*log(c[i]+exp(theta[i]*tau));
	}
	return summer;
}
double phi0_Planck_Einstein2::dTau(double tau, double delta)
{
	double summer=0;
	for (int i=iStart;i<=iEnd;i++)
	{
		summer+=a[i]*theta[i]*exp(tau*theta[i])/(c[i]+exp(theta[i]*tau));
	}
	return summer;
}
double phi0_Planck_Einstein2::dTau2(double tau, double delta)
{
	double summer=0;
	for (int i=iStart;i<=iEnd;i++)
	{
		summer+=a[i]*pow(theta[i],2)*c[i]*exp(tau*theta[i])/pow(c[i]+exp(tau*theta[i]),2);
	}
	return summer;
}
double phi0_Planck_Einstein2::dTau3(double tau, double delta)
{
	double summer=0;
	for (int i=iStart;i<=iEnd;i++)
	{
		summer += a[i]*pow(theta[i],2.0)*c[i]*(-theta[i])*exp(theta[i]*tau)*(exp(theta[i]*tau)-c[i])/pow(exp(theta[i]*tau)+c[i],3.0);
	}
	return summer;
}

double phi0_Planck_Einstein3::base(double tau, double delta)
{
	double summer=0;
	for (int i=iStart;i<=iEnd;i++)
	{
		summer += a[i]*log(exp(theta[i]*tau)-1);
	}
	return summer;
}
double phi0_Planck_Einstein3::dTau(double tau, double delta)
{
	double summer=0;
	for (int i=iStart;i<=iEnd;i++)
	{
		summer += a[i]*theta[i]*exp(theta[i]*tau)/(exp(theta[i]*tau)-1);
	}
	return summer;
}
double phi0_Planck_Einstein3::dTau2(double tau, double delta)
{
	double summer=0;
	for (int i=iStart;i<=iEnd;i++)
	{
		summer += -a[i]*pow(theta[i],2.0)*exp(theta[i]*tau)/pow(exp(theta[i]*tau)-1,2.0);
	}
	return summer;
}
double phi0_Planck_Einstein3::dTau3(double tau, double delta)
{
	double summer=0;
	for (int i=iStart;i<=iEnd;i++)
	{
		summer += a[i]*pow(theta[i],3.0)*theta[i]*exp(theta[i]*tau)*(exp(theta[i]*tau)+1)/pow(exp(theta[i]*tau)-1,3.0);
	}
	return summer;
}

/*
Maxima code for the sinh term:
part a)
((Tc*chi/tau)/sinh(Tc*chi/tau))^2;
-integrate(%,tau);
ratsimp(%);
part b)
((Tc*chi/tau)/sinh(Tc*chi/tau))^2;
integrate(%/tau,tau);
ratsimp(%);
Swap cosh for sinh and do it again
*/
double phi0_cp0_AlyLee::base(double tau, double delta)
{	
	return -tau/R_u*(anti_deriv_cp0_tau2(tau)-anti_deriv_cp0_tau2(tau0))+1/R_u*(anti_deriv_cp0_tau(tau)-anti_deriv_cp0_tau(tau0));
}
double phi0_cp0_AlyLee::dTau(double tau, double delta)
{
	// combining the integral terms for dTau yields
	// -1/Rbar*int(cp0/tau^2,dtau,tau0,tau)

	return -1/R_u*(anti_deriv_cp0_tau2(tau) - anti_deriv_cp0_tau2(tau0));
}
double phi0_cp0_AlyLee::anti_deriv_cp0_tau2(double tau)
{
	/*
	Maxima code:
	a[1]+a[2]*(a[3]*tau/Tc/sinh(a[3]*tau/Tc))^2+a[4]*(a[5]*tau/Tc/cosh(a[5]*tau/Tc))^2;
	integrate(%/tau^2,tau);
	*/
	return (4*a[4]*a[5])/(Tc*(2*exp(-(2*a[5]*tau)/Tc)+2))+(4*a[2]*a[3])/(Tc*(2*exp(-(2*a[3]*tau)/Tc)-2))-a[1]/tau;
}
double phi0_cp0_AlyLee::anti_deriv_cp0_tau(double tau)
{
	double term1;
	/*
	Maxima code:
	a[1]+a[2]*(a[3]*tau/Tc/sinh(a[3]*tau/Tc))^2+a[4]*(a[5]*tau/Tc/cosh(a[5]*tau/Tc))^2;
	integrate(%/tau,tau);
	*/
	if (a[4] == 0.0 && a[5] == 0.0)
	{
		term1 = 0;
	}
	else
	{
		term1 = (4*a[4]*a[5]*a[5]*((tau*Tc*exp((2*a[5]*tau)/Tc))/(2*a[5]*exp((2*a[5]*tau)/Tc)+2*a[5])-(Tc*Tc*log(exp((2*a[5]*tau)/Tc)+1))/(4*a[5]*a[5])))/Tc/Tc;
	}

	double term2 = (4*a[2]*a[3]*a[3]*((Tc*Tc*log(exp((a[3]*tau)/Tc)+1))/(4*a[3]*a[3])+(Tc*Tc*log(exp((a[3]*tau)/Tc)-1))/(4*a[3]*a[3])-(tau*Tc*exp((2*a[3]*tau)/Tc))/(2*a[3]*exp((2*a[3]*tau)/Tc)-2*a[3])))/Tc/Tc;
	double term3 = a[1]*log(tau);
	return term1 + term2 + term3;
}
double phi0_cp0_AlyLee::cp0(double tau)
{
	return a[1]+a[2]*pow(a[3]*tau/Tc/sinh(a[3]*tau/Tc),2)+a[4]*pow(a[5]*tau/Tc/(cosh(a[5]*tau/Tc)),2);
}
double phi0_cp0_AlyLee::dTau2(double tau, double delta)
{
	// The first integral term goes away, leaving just the second partial of the term (1/Rbar)*int(cp0/tau,dtau,tau0,tau)
	// which is equal to 1/Rbar*((tau*dcp0_dtau-cp0)/tau^2)
	return -cp0(tau)/(tau*tau*R_u);
}
double phi0_cp0_AlyLee::dTau3(double tau, double delta)
{
	// -cp0/tau/tau/R_u = -a[1]/tau^2/R_u-a[2]/R_u*(a[3]/Tc/sinh(a[3]*tau/Tc))^2-a[4]/R_u*(a[5]/Tc/cosh(a[5]*tau/Tc))^2;
	return 2*a[1]/tau/tau/tau/R_u     -a[2]/R_u*(-2)*pow(a[3]/Tc,3)*cosh(a[3]*tau/Tc)/pow(sinh(a[3]*tau/Tc),3)      -a[4]/R_u*(-2)*pow(a[5]/Tc,3)*sinh(a[5]*tau/Tc)/pow(cosh(a[5]*tau/Tc),3);
}

double phi0_cp0_poly::dTau(double tau, double delta)
{
	double sum=0;
	for (int i = iStart; i<=iEnd; i++){
		double t=tv[i];
		if (fabs(t)<10*DBL_EPSILON)
		{
			sum += a[i]/tau-a[i]/tau0;
		}
		else if (fabs(t+1) < 10*DBL_EPSILON)
		{
			sum += a[i]/Tc*log(tau0/tau);
		}
		else
		{
			sum+=a[i]*pow(Tc,t)*pow(tau,-t-1)/(t+1)-a[i]*pow(Tc,t)/(pow(tau0,t+1)*(t+1));
		}
	}
	return sum;
}

double phi0_cp0_poly::dTau2(double tau, double delta)
{
	double sum=0;
	for (int i = iStart; i<=iEnd; i++){
		double t = tv[i];
		if (fabs(t)<10*DBL_EPSILON)
		{
			sum += -a[i]/(tau*tau);
		}
		else if (fabs(t+1) < 10*DBL_EPSILON)
		{
			sum += -a[i]/(tau*Tc);
		}
		else
		{
			sum += -a[i]*pow(Tc/tau,tv[i])/(tau*tau);
		}
	}
	return sum;
}

double phi0_cp0_poly::dTau3(double tau, double delta)
{
	double sum=0;
	for (int i = iStart; i<=iEnd; i++){
		double t = tv[i];
		if (fabs(t)<10*DBL_EPSILON)
		{
			sum += 2*a[i]/(tau*tau*tau);
		}
		else if (fabs(t+1) < 10*DBL_EPSILON)
		{
			sum += a[i]/(tau*tau*Tc);
		}
		else
		{
			sum += a[i]*pow(Tc/tau,tv[i])*(tv[i]+2)/(tau*tau*tau);
		}
	}
	return sum;
}