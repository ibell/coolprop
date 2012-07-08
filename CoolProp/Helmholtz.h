#ifndef HELMHOLTZ_H
#define HELMHOLTZ_H

#include <iostream>
#include <vector>
#include "math.h"

/// This is the abstract base class upon which each residual Helmholtz energy class is built
class phi_BC{
public:
	phi_BC(){};
	virtual ~phi_BC(){};
	// Term and its derivatives
	/// Returns the base, non-dimensional, Helmholtz energy term (no derivatives) [-]
	/// @param tau Reciprocal reduced temperature where tau=Tc/T
	/// @param delta Reduced pressure where delta = rho / rhoc 
	virtual double base(double tau,double delta) = 0;
	/// Returns the first partial derivative of Helmholtz energy term with respect to tau [-]
	/// @param tau Reciprocal reduced temperature where tau=Tc/T
	/// @param delta Reduced pressure where delta = rho / rhoc 
	virtual double dTau(double tau,double delta) = 0;
	/// Returns the second partial derivative of Helmholtz energy term with respect to tau [-]
	/// @param tau Reciprocal reduced temperature where tau=Tc/T
	/// @param delta Reduced pressure where delta = rho / rhoc 
	virtual double dTau2(double tau, double delta) = 0;
	/// Returns the second mixed partial derivative (delta1,dtau1) of Helmholtz energy term with respect to delta and tau [-]
	/// @param tau Reciprocal reduced temperature where tau=Tc / T
	/// @param delta Reduced pressure where delta = rho / rhoc 
	virtual double dDelta_dTau(double tau, double delta) = 0;
	/// Returns the first partial derivative of Helmholtz energy term with respect to delta [-]
	/// @param tau Reciprocal reduced temperature where tau=Tc / T
	/// @param delta Reduced pressure where delta = rho / rhoc 
	virtual double dDelta(double tau, double delta) = 0;
	/// Returns the second partial derivative of Helmholtz energy term with respect to delta [-]
	/// @param tau Reciprocal reduced temperature where tau=Tc / T
	/// @param delta Reduced pressure where delta = rho / rhoc 
	virtual double dDelta2(double tau, double delta) = 0;
	/// Returns the third mixed partial derivative (delta2,dtau1) of Helmholtz energy term with respect to delta and tau [-]
	/// @param tau Reciprocal reduced temperature where tau=Tc / T
	/// @param delta Reduced pressure where delta = rho / rhoc 
	virtual double dDelta2_dTau(double tau, double delta) = 0;
};

/// This class implements residual Helmholtz Energy terms of the form  n * delta ^d * tau^t * exp(-gamma*delta^l) if l>0 or if
///	l==0, then n * delta ^d * tau^t
class phir_power : public phi_BC{
	/*
	Terms are of the form n * delta ^d * tau^t * exp(-delta^l) if l>0 or if
	l==0, then n * delta ^d * tau^t

	Constructor must be called with std::vector instances of double type
	*/
private:
	std::vector<double> n, ///< The coefficients multiplying each term
		                d, ///< The power for the delta terms
						t, ///< The powers for the tau terms
						l; ///< The powers for delta in the exp terms
	unsigned int iStart,iEnd;
public:
	// Constructors
	phir_power(std::vector<double>,std::vector<double>,std::vector<double>,int,int);
	phir_power(std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,int,int);

	///< Destructor for the phir_power class.  No implementation
	~phir_power(){};

	double base(double tau, double delta) throw();
	double dTau(double tau, double delta) throw();
	double dTau2(double tau, double delta) throw();
	double dDelta(double tau, double delta) throw();
	double dDelta2(double tau, double delta) throw();
	double dDelta2_dTau(double tau, double delta) throw();
	double dDelta_dTau(double tau, double delta) throw();
};

class phir_gaussian : public phi_BC{
	/*
	Terms are of the form n * delta ^d * tau^t * exp(-alpha*(delta-epsilon)^2-beta*(tau-gamma)^2)

	Constructor must be called with std::vector instances of double type
	*/
private:
	std::vector<double> n,d,t,alpha,epsilon,beta,gamma;
	unsigned int iStart,iEnd;
public:
	// Constructors
	phir_gaussian(std::vector<double> a_in, 
				  std::vector<double> d_in,
				  std::vector<double> t_in, 
				  std::vector<double> alpha_in, 
				  std::vector<double> epsilon_in, 
				  std::vector<double> beta_in, 
				  std::vector<double> gamma_in,
		unsigned int iStart_in, unsigned int iEnd_in);

	// Destructor
	~phir_gaussian(){};

	// Term and its derivatives
	double base(double tau, double delta);
	double dTau(double tau, double delta);
	double dTau2(double tau, double delta);
	double dDelta(double tau, double delta);
	double dDelta2(double tau, double delta);
	double dDelta2_dTau(double tau, double delta);
	double dDelta_dTau(double tau, double delta);
};

class phir_critical : public phi_BC{
	/*
	Terms are of the form n * delta ^d * tau^t * exp(-alpha*(delta-epsilon)^2-beta*(tau-gamma)^2)
	Constructor must be called with std::vector instances of double type
	*/
private:
	std::vector<double> n,d,t,a,b,A,B,C,D,beta;
	int iStart,iEnd;
public:
	// Constructors
	phir_critical(std::vector<double> n_in, std::vector<double> d_in, std::vector<double> t_in, 
		std::vector<double> a_in, std::vector<double> b_in, std::vector<double> beta_in,
		std::vector<double> A_in, std::vector<double> B_in, std::vector<double> C_in, 
		std::vector<double> D_in, int iStart_in, int iEnd_in);

	//Destructor
	~phir_critical(){};

	// Term and its derivatives
	double base(double tau, double delta);
	double dTau(double tau, double delta);
	double dTau2(double tau, double delta);
	double dDelta(double tau, double delta);
	double dDelta2(double tau, double delta);
	double dDelta2_dTau(double tau, double delta);
	double dDelta_dTau(double tau, double delta);
};

class phi0_lead : public phi_BC{
	/*
	Term is of the form log(delta)+a1+a2*tau
	*/
private:
	double c1,c2; // Use these variables internally
public:
	// Constructor
	phi0_lead(double a1, double a2){c1=a1; c2=a2;};

	//Destructor
	~phi0_lead(){};

	// Term and its derivatives
	double base(double tau, double delta){return log(delta)+c1+c2*tau;};
	double dTau(double tau, double delta){return c2;};
	double dTau2(double tau, double delta){return 0.0;};
	double dDelta(double tau, double delta){return 1.0/delta;};
	double dDelta2(double tau, double delta){return -1.0/delta/delta;};
	double dDelta2_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau(double tau, double delta){return 0.0;};

};

class phi0_logtau : public phi_BC{
	/*
	Term is of the form a1*log(tau)
	*/
private:
	double c1; // Use these variables internally
public:
	// Constructor
	phi0_logtau(double a1){c1=a1;};

	//Destructor
	~phi0_logtau(){};

	// Term and its derivatives
	double base(double tau, double delta){return c1*log(tau);};
	double dTau(double tau, double delta){return c1/tau;};
	double dTau2(double tau, double delta){return -c1/tau/tau;};
	double dDelta(double tau, double delta){return 0.0;};
	double dDelta2(double tau, double delta){return 0.0;};
	double dDelta2_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau(double tau, double delta){return 0.0;};
};

class phi0_Planck_Einstein : public phi_BC{
	/*
	Term is of the form a_0*log(1-exp(-theta_0*tau))
	*/
private:
	std::vector<double> a,theta; // Use these variables internally
	int iStart, iEnd;
public:
	// Constructor
	phi0_Planck_Einstein(std::vector<double> a_in, std::vector<double> theta_in, int iStart_in, int iEnd_in)
	{
		a=a_in; theta=theta_in; iStart = iStart_in; iEnd = iEnd_in;
	};

	//Destructor
	~phi0_Planck_Einstein(){};

	// Term and its derivatives
	double base(double tau, double delta);
	double dTau(double tau, double delta);
	double dTau2(double tau, double delta);
	double dDelta(double tau, double delta){return 0.0;};
	double dDelta2(double tau, double delta){return 0.0;};
	double dDelta2_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau(double tau, double delta){return 0.0;};
};

class phi0_Planck_Einstein2 : public phi_BC{
	/*
	Term is of the form a_0*log(c+exp(theta_0*tau))
	*/
private:
	std::vector<double> a,theta,c; // Use these variables internally
	int iStart, iEnd;
public:
	// Constructor
	phi0_Planck_Einstein2(double a_in, double theta_in, double c_in)
	{
		a=std::vector<double> (1,a_in); 
		theta=std::vector<double> (1,theta_in); 
		c=std::vector<double> (1,c_in); 
		iStart = 0; iEnd = 0;
	};

	//Destructor
	~phi0_Planck_Einstein2(){};

	// Term and its derivatives
	double base(double tau, double delta);
	double dTau(double tau, double delta);
	double dTau2(double tau, double delta);
	double dDelta(double tau, double delta){return 0.0;};
	double dDelta2(double tau, double delta){return 0.0;};
	double dDelta2_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau(double tau, double delta){return 0.0;};
};

class phi0_power : public phi_BC{
	/*
	Term is of the form a_0*tau^b_0
	*/
private:
	std::vector<double> a,b; // Use these variables internally
	int iStart, iEnd;
public:
	// Constructor
	phi0_power(std::vector<double> a_in, std::vector<double> theta_in, int iStart_in, int iEnd_in)
	{
		a=a_in; b=theta_in; iStart = iStart_in; iEnd = iEnd_in;
	};
	// Constructor with just double values
	phi0_power(double a_in, double theta_in)
	{
		a=std::vector<double>(1,a_in); b=std::vector<double>(1,theta_in); iStart = 0; iEnd = 0;
	};

	//Destructor
	~phi0_power(){};

	// Term and its derivatives
	double base(double tau, double delta)
	{
		double summer=0;
		for (int i=iStart;i<=iEnd;i++) summer+=a[i]*pow(tau,b[i]);
		return summer;
	};
	double dTau(double tau, double delta)
	{
		double summer=0;
		
		for (int i=iStart;i<=iEnd;i++) { 
			summer+=a[i]*b[i]*pow(tau,b[i]-1);
		}
		return summer;
	};
	double dTau2(double tau, double delta)
	{
		double summer=0;
		for (int i=iStart;i<=iEnd;i++) summer+=a[i]*b[i]*(b[i]-1)*pow(tau,b[i]-2);
		return summer;
	};
	double dDelta(double tau, double delta){return 0.0;};
	double dDelta2(double tau, double delta){return 0.0;};
	double dDelta2_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau(double tau, double delta){return 0.0;};
};


/// Term in the ideal-gas specific heat equation that is constant
class phi0_cp0_constant : public phi_BC{
	/*

	Maxima code for this term:
	assume(T>0)$
	assume(T0>0)$
	assume(T-T0>0)$
	a:(1/T)*integrate(c,T,T0,T)-integrate(c/T,T,T0,T)$
	subst(Tc/tau,T,a)$
	subst(Tc/tau0,T0,%)$
	b:ratsimp(logcontract(%));
	db:ratsimp(diff(b,tau));
	db2:ratsimp(diff(%,tau));

	*/
private:
	double c,Tc,T0,tau0; // Use these variables internally
public:

	/// Constructor with just a single double value
	phi0_cp0_constant(double c_, double Tc_, double T0_) { c=c_; T0=T0_; Tc=Tc_; tau0=Tc/T0;};

	/// Destructor
	~phi0_cp0_constant(){};

	// Term and its derivatives
	double base(double tau, double delta){ 
		return c-c*tau/tau0+c*log(tau/tau0);
	};
	double dTau(double tau, double delta)
	{
		return c/tau-c/tau0;
	};
	double dTau2(double tau, double delta)
	{
		return -c/(tau*tau);
	};
	double dDelta(double tau, double delta){return 0.0;};
	double dDelta2(double tau, double delta){return 0.0;};
	double dDelta2_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau(double tau, double delta){return 0.0;};
};

/// Term in the ideal-gas specific heat equation that is constant
class phi0_cp0_exponential : public phi_BC{
private:
	std::vector<double> a,b;
	double Tc,T0,tau0; // Use these variables internally
	int iStart, iEnd;
public:
	/// Constructor with just a single double value
	phi0_cp0_exponential(std::vector<double> a_, std::vector<double> b_, double Tc_, double T0_, int iStart_, int iEnd_) { 
		a=a_; b=b_; T0=T0_; Tc=Tc_; iStart=iStart_; iEnd=iEnd_; tau0=Tc/T0;
	};

	/// Destructor
	~phi0_cp0_exponential(){};

	// Term and its derivatives
	double base(double tau, double delta){ 
		double sum=0;
		for (int i = iStart; i<=iEnd; i++){
			sum+=a[i]*(log(exp(b[i]*tau)-1)-log(exp(b[i]*tau0)-1)-b[i]*tau+b[i]*(tau0*exp(b[i]*tau0)-tau)/(exp(b[i]*tau0)-1));
		}
		return sum;
	};
	double dTau(double tau, double delta)
	{
		double sum=0;
		for (int i = iStart; i<=iEnd; i++){
			sum+=a[i]*(b[i]*exp(b[i]*tau)/(exp(b[i]*tau)-1)-b[i]-b[i]/(exp(b[i]*tau0)-1));
		}
		return sum;
	};
	double dTau2(double tau, double delta)
	{
		double sum=0;
		for (int i = iStart; i<=iEnd; i++){
			sum+=-a[i]*b[i]*b[i]*exp(b[i]*tau)/pow(exp(b[i]*tau)-1,2);
		}
		return sum;
	};
	double dDelta(double tau, double delta){return 0.0;};
	double dDelta2(double tau, double delta){return 0.0;};
	double dDelta2_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau(double tau, double delta){return 0.0;};
};

/// Term in the ideal-gas specific heat equation that is polynomial term
class phi0_cp0_poly : public phi_BC{
private:
	std::vector<double> a,tv;
	double Tc,T0,tau0; // Use these variables internally
	int iStart, iEnd;
public:
	/// Constructor with just a single double value
	phi0_cp0_poly(double a_, double t_, double Tc_, double T0_) {
		a=std::vector<double>(1,a_);
		tv=std::vector<double>(1,t_);
		Tc=Tc_; T0=T0_; iStart=0; iEnd=0; tau0=Tc/T0;
	};

	/// Constructor with std::vectors
	phi0_cp0_poly(std::vector<double> a_, std::vector<double> t_, double Tc_, double T0_, int iStart_, int iEnd_) { 
		a=a_; tv=t_; Tc = Tc_; T0=T0_; iStart=iStart_; iEnd=iEnd; tau0=Tc/T0;
	};

	/// Destructor
	~phi0_cp0_poly(){};

	// Term and its derivatives
	double base(double tau, double delta){ 
		double sum=0;
		for (int i = iStart; i<=iEnd; i++){
			double t=tv[i];
			sum+=-a[i]*pow(Tc,t)*pow(tau,-t)/(t*(t+1))-a[i]*pow(T0,t+1)*tau/(Tc*(t+1))+a[i]*pow(T0,t)/t;
		}
		return sum;
	};
	double dTau(double tau, double delta)
	{
		double sum=0;
		for (int i = iStart; i<=iEnd; i++){
			double t=tv[i];
			sum+=a[i]*pow(Tc,t)*pow(tau,-t-1)/(t+1)-a[i]*pow(Tc,t)/(pow(tau0,t+1)*(t+1));
		}
		return sum;
	};
	double dTau2(double tau, double delta)
	{
		double sum=0;
		for (int i = iStart; i<=iEnd; i++){
			double t=tv[i];
			sum+=-a[i]*pow(Tc,t)/pow(tau,t+2);
		}
		return sum;
	};
	double dDelta(double tau, double delta){return 0.0;};
	double dDelta2(double tau, double delta){return 0.0;};
	double dDelta2_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau(double tau, double delta){return 0.0;};
};

#endif