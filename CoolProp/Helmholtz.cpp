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
	vector<double> alpha_in, vector<double> epsilon_in, vector<double> beta_in, vector<double> gamma_in,unsigned int iStart_in, unsigned int iEnd_in)
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
double phir_gaussian::dDelta2_dTau(double tau, double delta)
{
	double summer=0,psi;
	for (unsigned int i=iStart;i<=iEnd;i++)
	{
		psi=exp(-alpha[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
		summer+=n[i]*pow(tau,t[i])*psi*pow(delta,d[i])*(t[i]/tau-2.0*beta[i]*(tau-gamma[i]))*(-2.0*alpha[i]*(delta-epsilon[i])+4*alpha[i]*alpha[i]*pow(delta,d[i])*pow(delta-epsilon[i],2)-4*d[i]*alpha[i]*pow(delta,d[i]-1)*(delta-epsilon[i])+d[i]*(d[i]-1)*pow(delta,d[i]-2));
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
		summer-=a[i]*pow(theta[i],2)*exp(-theta[i]*tau)/pow(1.0-exp(-theta[i]*tau),2);
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