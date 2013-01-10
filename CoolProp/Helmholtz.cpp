#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include "time.h"
#include <vector>
#include <iostream>
#include "math.h"
using namespace std;
#include "Helmholtz.h"
#include "CoolPropTools.h"

// Constructors
phir_power::phir_power(std::vector<double> n_in,std::vector<double> d_in,std::vector<double> t_in, std::vector<double> l_in, int iStart_in,int iEnd_in)
{
	n=n_in;
	d=d_in;
	t=t_in;
	l=l_in;
	iStart=iStart_in;
	iEnd=iEnd_in;
}
phir_power::phir_power(std::vector<double> n_in,std::vector<double> d_in,std::vector<double> t_in,int iStart_in,int iEnd_in)
{
	n=n_in;
	d=d_in;
	t=t_in;
	l.assign(d.size(),0.0);
	iStart=iStart_in;
	iEnd=iEnd_in;
}
phir_power::phir_power(const double n_in[], const double d_in[], const double t_in[],int iStart_in,int iEnd_in, int N)
{
	n=std::vector<double>(n_in,n_in+N);
	d=std::vector<double>(d_in,d_in+N);
	t=std::vector<double>(t_in,t_in+N);
	l.assign(d.size(),0.0);
	iStart=iStart_in;
	iEnd=iEnd_in;
}
phir_power::phir_power(double n_in[], double d_in[], double t_in[], double l_in[], int iStart_in,int iEnd_in, int N)
{
	n=std::vector<double>(n_in,n_in+N);
	d=std::vector<double>(d_in,d_in+N);
	t=std::vector<double>(t_in,t_in+N);
	l=std::vector<double>(l_in,l_in+N);
	iStart=iStart_in;
	iEnd=iEnd_in;
}

// Term and its derivatives
double phir_power::base(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		if (l[i]>0)
			summer+=n[i]*exp(t[i]*log_tau+d[i]*log_delta-pow(delta,l[i]));
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
			summer+=n[i]*t[i]*exp((t[i]-1)*log_tau+d[i]*log_delta-pow(delta,l[i]));
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
			summer+=n[i]*t[i]*(t[i]-1)*exp((t[i]-2)*log_tau+d[i]*log_delta-pow(delta,l[i]));
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
			summer+=n[i]*t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta-pow(delta,l[i]));
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
		double pow_delta_li = pow(delta,l[i]);
		if (l[i]>0)
			summer+=n[i]*t[i]*(t[i]-1)*(d[i]-l[i]*pow_delta_li)*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta-pow(delta,l[i]));
		else
			summer+=n[i]*t[i]*(t[i]-1)*d[i]*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta);
	}
	return summer;
}
double phir_power::dDelta(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		double pow_delta_li = pow(delta,l[i]);
		if (l[i]>0)
			summer+=n[i]*(d[i]-l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-1)*log_delta-pow_delta_li);
		else
			summer+=n[i]*d[i]*exp(t[i]*log_tau+(d[i]-1)*log_delta);
	}
	return summer;
}
double phir_power::dDelta2(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		double pow_delta_li = pow(delta,l[i]);
		if (l[i]>0)
			summer+=n[i]*((d[i]-l[i]*pow_delta_li)*(d[i]-1.0-l[i]*pow_delta_li) - l[i]*l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-2)*log_delta-pow(delta,l[i]));
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
		double pow_delta_li = pow(delta,l[i]);
		if (l[i]>0)
		{
			double bracket = (d[i]*(d[i]-1)*(d[i]-2)+pow_delta_li*(-2*l[i]+6*d[i]*l[i]-3*d[i]*d[i]*l[i]-3*d[i]*l[i]*l[i]+3*l[i]*l[i]-l[i]*l[i]*l[i])+pow_delta_li*pow_delta_li*(3*d[i]*l[i]*l[i]-3*l[i]*l[i]+3*l[i]*l[i]*l[i])-l[i]*l[i]*l[i]*pow_delta_li*pow_delta_li*pow_delta_li);
			summer+=n[i]*bracket*exp(t[i]*log_tau+(d[i]-3)*log_delta-pow(delta,l[i]));
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
		double pow_delta_li = pow(delta,l[i]);
		if (l[i]>0)
			summer+=n[i]*t[i]*(((d[i]-l[i]*pow_delta_li))*(d[i]-1-l[i]*pow_delta_li)-l[i]*l[i]*pow_delta_li)*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta-pow(delta,l[i]));
		else
			summer+=n[i]*d[i]*t[i]*(d[i]-1)*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta);;
	}
	return summer;
}
double phir_power::dDelta_dTau(double tau, double delta) throw()
{
	double summer=0, log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		double pow_delta_li = pow(delta,l[i]);
		if (l[i]>0)
			summer+=n[i]*t[i]*(d[i]-l[i]*pow_delta_li)*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta-pow_delta_li);
		else
			summer+=n[i]*d[i]*t[i]*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta);
	}
	return summer;
}

phir_gaussian::phir_gaussian(vector<double> n_in, vector<double> d_in,vector<double> t_in, 
	vector<double> alpha_in, vector<double> epsilon_in, vector<double> beta_in, vector<double> gamma_in,
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

// Term and its derivatives
double phir_gaussian::base(double tau, double delta)
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(delta,d[i])*pow(tau,t[i])*psi;
	}
	return summer;
}
double phir_gaussian::dTau(double tau, double delta)
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(delta,d[i])*pow(tau,t[i])*psi*(t[i]/tau-2.0*beta[i]*(tau-gamma[i]));
	}
	return summer;
}
double phir_gaussian::dTau2(double tau, double delta)
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(delta,d[i])*pow(tau,t[i])*psi*(pow(t[i]/tau-2.0*beta[i]*(tau-gamma[i]),2)-t[i]/pow(tau,2)-2.0*beta[i]);
	}
	return summer;
}
double phir_gaussian::dTau3(double tau, double delta)
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
double phir_gaussian::dDelta(double tau, double delta)
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(delta,d[i])*pow(tau,t[i])*psi*(d[i]/delta-2.0*alpha[i]*(delta-epsilon[i]));
	}
	return summer;
}
double phir_gaussian::dDelta2(double tau, double delta)
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(tau,t[i])*psi*(-2.0*alpha[i]*pow(delta,d[i])+4.0*pow(alpha[i],2)*pow(delta,d[i])*pow(delta-epsilon[i],2)-4.0*d[i]*alpha[i]*pow(delta,d[i]-1)*(delta-epsilon[i])+d[i]*(d[i]-1.0)*pow(delta,d[i]-2));
	}
	return summer;
}
double phir_gaussian::dDelta3(double tau, double delta)
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
double phir_gaussian::dDelta2_dTau(double tau, double delta)
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(tau,t[i])*psi*(t[i]/tau-2.0*beta[i]*(tau-gamma[i]))*(-2.0*alpha[i]*pow(delta,d[i])+4.0*pow(alpha[i],2)*pow(delta,d[i])*pow(delta-epsilon[i],2)-4.0*d[i]*alpha[i]*pow(delta,d[i]-1)*(delta-epsilon[i])+d[i]*(d[i]-1.0)*pow(delta,d[i]-2));
	}
	return summer;
}
double phir_gaussian::dDelta_dTau(double tau, double delta)
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(delta,d[i])*pow(tau,t[i])*psi*(d[i]/delta-2.0*alpha[i]*(delta-epsilon[i]))*(t[i]/tau-2.0*beta[i]*(tau-gamma[i]));
	}
	return summer;
}
double phir_gaussian::dDelta_dTau2(double tau, double delta)
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(delta,d[i])*pow(tau,t[i])*psi*(d[i]/delta-2.0*alpha[i]*(delta-epsilon[i]))*(pow(t[i]-2.0*beta[i]*tau*(tau-gamma[i]),2)-t[i]-2*beta[i]*tau*tau)/tau/tau;
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

double phir_critical::base(double tau, double delta)
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

double phir_critical::dDelta(double tau, double delta)
{
	double summer=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta;
	for (int i=iStart;i<=iEnd;i++)
	{
		theta=(1.0-tau)+A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(pow(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        summer+=n[i]*(pow(DELTA,b[i])*(PSI+delta*dPSI_dDelta)+dDELTAbi_dDelta*delta*PSI);
	}
	return summer;
}
double phir_critical::dDelta_dTau2(double tau, double delta)
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

double phir_critical::dDelta2(double tau, double delta)
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
double phir_critical::dDelta3(double tau, double delta)
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
		PI = 4*B[i]*a[i]*(a[i]-1)*pow(delta-1,2*a[i]-4)+2*pow(A[i]/beta[i],2)*pow(delta-1,2/beta[i]-4)+4*A[i]*theta/beta[i]*(1/(2*beta[i])-1)*pow(delta-1,1/beta[i]-4);
		dPI_dDelta = 8*B[i]*a[i]*(a[i]-1)*(a[i]-2)*pow(delta-1,2*a[i]-5)+8*pow(A[i]/beta[i],2)*(1/(2*beta[i])-2)*pow(delta-1,2/beta[i]-5)+(8*A[i]*theta)/beta[i]*(1/(2*beta[i])-1)*(1/(2*beta[i])-2)*pow(delta-1,1/beta[i]-5)+4*A[i]/beta[i]*(1/(2*beta[i])-1)*pow(delta-1,1/beta[i]-4)*dtheta_dDelta;
		dDELTA3_dDelta3 = 1/(delta-1)*dDELTA2_dDelta2-1/pow(delta-1,2)*dDELTA_dDelta+pow(delta-1,2)*dPI_dDelta+2*(delta-1)*PI;
		dDELTAbi3_dDelta3 = b[i]*(pow(DELTA,b[i]-1)*dDELTA3_dDelta3+dDELTA2_dDelta2*(b[i]-1)*pow(DELTA,b[i]-2)*dDELTA_dDelta+(b[i]-1)*(pow(DELTA,b[i]-2)*2*dDELTA_dDelta*dDELTA2_dDelta2+pow(dDELTA_dDelta,2)*(b[i]-2)*pow(DELTA,b[i]-3)*dDELTA_dDelta));
        
		summer += n[i]*(pow(DELTA,b[i])*(3.0*dPSI2_dDelta2+delta*dPSI3_dDelta3)+3.0*dDELTAbi_dDelta*(2*dPSI_dDelta+delta*dPSI2_dDelta2)+3*dDELTAbi2_dDelta2*(PSI+delta*dPSI_dDelta)+dDELTAbi3_dDelta3*PSI*delta);
	}
	return summer;
}

double phir_critical::dDelta2_dTau(double tau, double delta)
{
	double summer=0;
	double dphir3_dDelta2_dTau=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta,dPSI2_dDelta2,dDELTAbi2_dDelta2,dDELTA2_dDelta2;
    double dPSI2_dDelta_dTau, dDELTAbi2_dDelta_dTau, dPSI_dTau, dDELTAbi_dTau;
    double Line1,Line2,Line3,dDELTA2_dDelta_dTau,dDELTA3_dDelta2_dTau,dDELTAbim1_dTau,dDELTAbim2_dTau;
    double dDELTA_dTau,dDELTAbi3_dDelta2_dTau;

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
		dDELTA_dTau=-2*((1-tau)+A[i]*pow(pow(delta-1,2),1/(2*beta[i])-1)+2*B[i]*a[i]*pow(pow(delta-1,2),a[i]-1));
		dDELTA2_dDelta_dTau=-(delta-1)*A[i]*(2/beta[i])*pow(pow(delta-1,2),1/(2*beta[i])-1);
		dDELTA3_dDelta2_dTau=1/(delta-1)*dDELTA2_dDelta_dTau-pow(delta-1,2)*A[i]*(4/beta[i])*(1/(2*beta[i])-1)*pow(pow(delta-1,2),1/(2*beta[i])-2);
		
		dDELTAbim1_dTau=(b[i]-1)*pow(DELTA,b[i]-2)*dDELTA_dTau;
		dDELTAbim2_dTau=(b[i]-2)*pow(DELTA,b[i]-3)*dDELTA_dTau;
		Line1=dDELTAbim1_dTau*dDELTA2_dDelta2+pow(DELTA,b[i]-1)*dDELTA3_dDelta2_dTau;
		Line2=(b[i]-1)*(dDELTAbim2_dTau*pow(dDELTA_dDelta,2)+pow(DELTA,b[i]-2)*2*dDELTA2_dDelta_dTau*dDELTA_dDelta);
		dDELTAbi3_dDelta2_dTau=b[i]*(Line1+Line2);
		
		Line1=pow(DELTA,b[i])*(2*delta*dPSI2_dDelta_dTau+delta*dDELTA3_dDelta2_dTau)+dDELTAbi_dTau*(2*dPSI_dDelta+delta*dPSI2_dDelta2);
		Line2=2*dDELTAbi2_dDelta_dTau*(PSI+delta*dPSI_dDelta)+2*dDELTAbi_dDelta*(dPSI_dTau+delta*dPSI2_dDelta_dTau);
		Line3=dDELTAbi3_dDelta2_dTau*delta*PSI+dDELTAbi2_dDelta2*delta*dPSI_dTau;
        summer+=n[i]*(Line1+Line2+Line3);
    }
	return summer;
}

double phir_critical::dDelta_dTau(double tau, double delta)
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
	
double phir_critical::dTau(double tau, double delta)
{
	double summer=0,theta,DELTA,PSI,dPSI_dTau, dDELTAbi_dTau;

	for (int i=iStart;i<=iEnd;i++)
	{
		theta=(1.0-tau)+A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        summer+=n[i]*delta*(dDELTAbi_dTau*PSI+pow(DELTA,b[i])*dPSI_dTau);
	}
	return summer;
}

double phir_critical::dTau2(double tau, double delta)
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
double phir_critical::dTau3(double tau, double delta)
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
		sum+=a[i]*pow(Tc,t)*pow(tau,-t-1)/(t+1)-a[i]*pow(Tc,t)/(pow(tau0,t+1)*(t+1));
	}
	return sum;
}

double phi0_cp0_poly::dTau2(double tau, double delta)
{
	double sum=0;
	for (int i = iStart; i<=iEnd; i++){
		sum += -a[i]*pow(Tc/tau,tv[i])/(tau*tau);
	}
	return sum;
}

double phi0_cp0_poly::dTau3(double tau, double delta)
{
	double sum=0;
	for (int i = iStart; i<=iEnd; i++){
		sum += a[i]*pow(Tc/tau,tv[i])*(tv[i]+2)/(tau*tau*tau);
	}
	return sum;
}