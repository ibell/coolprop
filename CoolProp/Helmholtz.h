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
	/// Returns the third mixed partial derivative (delta1,dtau2) of Helmholtz energy term with respect to delta and tau [-]
	/// @param tau Reciprocal reduced temperature where tau=Tc / T
	/// @param delta Reduced pressure where delta = rho / rhoc 
	virtual double dDelta_dTau2(double tau, double delta) = 0;
	/// Returns the third partial derivative of Helmholtz energy term with respect to tau [-]
	/// @param tau Reciprocal reduced temperature where tau=Tc / T
	/// @param delta Reduced pressure where delta = rho / rhoc 
	virtual double dTau3(double tau, double delta) = 0;
	/// Returns the third partial derivative of Helmholtz energy term with respect to delta [-]
	/// @param tau Reciprocal reduced temperature where tau=Tc / T
	/// @param delta Reduced pressure where delta = rho / rhoc 
	virtual double dDelta3(double tau, double delta) = 0;
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
	// Add a dummy value at the end of the function since compiler sees std::vector<double> and double[] as being the same type
	phir_power(const double[],const double[],const double[],int,int,int);
	phir_power(double[],double[],double[],double[],int,int,int);

	///< Destructor for the phir_power class.  No implementation
	~phir_power(){};

	double base(double tau, double delta) throw();
	double dDelta(double tau, double delta) throw();
	double dTau(double tau, double delta) throw();
	
	double dDelta2(double tau, double delta) throw();
	double dDelta_dTau(double tau, double delta) throw();
	double dTau2(double tau, double delta) throw();
	
	double dDelta3(double tau, double delta);
	double dDelta2_dTau(double tau, double delta);
	double dDelta_dTau2(double tau, double delta);
	double dTau3(double tau, double delta);
};

class phir_gaussian : public phi_BC{
	/*
	Terms are of the form a * delta ^d * tau^t * exp(-alpha*(delta-epsilon)^2-beta*(tau-gamma)^2)

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
	phir_gaussian(double a_in[], 
				  double d_in[],
				  double t_in[], 
				  double alpha_in[], 
				  double epsilon_in[], 
				  double beta_in[], 
				  double gamma_in[],
		unsigned int iStart_in, unsigned int iEnd_in, unsigned int N);

	// Destructor
	~phir_gaussian(){};

	// Term and its derivatives
	double base(double tau, double delta) throw();
	double dDelta(double tau, double delta) throw();
	double dTau(double tau, double delta) throw();
	
	double dDelta2(double tau, double delta) throw();
	double dDelta_dTau(double tau, double delta) throw();
	double dTau2(double tau, double delta) throw();
	
	double dDelta3(double tau, double delta);
	double dDelta2_dTau(double tau, double delta);
	double dDelta_dTau2(double tau, double delta);
	double dTau3(double tau, double delta);
	
};

class phir_critical : public phi_BC{
	/*
	The critical term, used for Water and Carbon Dioxide
	It is truly horrible
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
	double base(double tau, double delta) throw();
	double dDelta(double tau, double delta) throw();
	double dTau(double tau, double delta) throw();
	
	double dDelta2(double tau, double delta) throw();
	double dDelta_dTau(double tau, double delta) throw();
	double dTau2(double tau, double delta) throw();
	
	double dDelta3(double tau, double delta);
	double dDelta2_dTau(double tau, double delta);
	double dDelta_dTau2(double tau, double delta);
	double dTau3(double tau, double delta);
};

class phi0_lead : public phi_BC{
	/*
	Term is of the form log(delta)+a1+a2*tau
	constructor: phi0_lead(double a1, double a2)
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
	double dDelta_dTau2(double tau, double delta){return 0.0;};
	double dTau3(double tau, double delta){return 0.0;};
	double dDelta3(double tau, double delta){return 2/delta/delta/delta;};
};

class phi0_logtau : public phi_BC{
	/*
	Term is of the form a1*log(tau)
	Constructor: phi0_logtau(double a1)
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
	double dDelta_dTau2(double tau, double delta){return 0.0;};
	double dTau3(double tau, double delta){return 2*c1/tau/tau/tau;};
	double dDelta3(double tau, double delta){return 0;};
};

class phi0_Planck_Einstein : public phi_BC{
	/*
	Term is of the form a_0*log(1-exp(-theta_0*tau))
	Constructors: 
		phi0_Planck_Einstein(std::vector<double> a_in, std::vector<double> theta_in, int iStart_in, int iEnd_in)
		phi0_Planck_Einstein(double a_in, double theta_in)
	*/
private:
	std::vector<double> a,theta; // Use these variables internally
	int iStart, iEnd;
public:
	// Constructor with std::vector instances
	phi0_Planck_Einstein(std::vector<double> a_in, std::vector<double> theta_in, int iStart_in, int iEnd_in)
	{
		a=a_in; theta=theta_in; iStart = iStart_in; iEnd = iEnd_in;
	};
	// Constructor with doubles
	phi0_Planck_Einstein(double a_in, double theta_in)
	{
		a=std::vector<double> (1,a_in); 
		theta=std::vector<double> (1,theta_in); 
		iStart = 0; iEnd = 0;
	};

	//Destructor
	~phi0_Planck_Einstein(){};

	// Term and its derivatives
	double base(double tau, double delta);
	double dTau(double tau, double delta);
	double dTau2(double tau, double delta);
	double dTau3(double tau, double delta);
	double dDelta(double tau, double delta){return 0.0;};
	double dDelta2(double tau, double delta){return 0.0;};
	double dDelta2_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau2(double tau, double delta){return 0.0;};
	double dDelta3(double tau, double delta){return 0;};
};

class phi0_Planck_Einstein2 : public phi_BC{
	/*
	Term is of the form a_0*log(c+exp(theta_0*tau))
	Constructor:
		phi0_Planck_Einstein2(double a_in, double theta_in, double c_in)
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
	double dTau3(double tau, double delta);
	double dDelta(double tau, double delta){return 0.0;};
	double dDelta2(double tau, double delta){return 0.0;};
	double dDelta2_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau2(double tau, double delta){return 0.0;};
	double dDelta3(double tau, double delta){return 0.0;};
};

class phi0_power : public phi_BC{
	/*
	Term is of the form a_0*tau^b_0
	Constructor:
		phi0_power(std::vector<double> a_in, std::vector<double> theta_in, int iStart_in, int iEnd_in)
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
		for (int i=iStart;i<=iEnd;i++){
			summer += a[i]*pow(tau,b[i]);
		}
		return summer;
	};
	double dTau(double tau, double delta)
	{
		double summer=0;
		
		for (int i=iStart;i<=iEnd;i++) { 
			summer += a[i]*b[i]*pow(tau,b[i]-1);
		}
		return summer;
	};
	double dTau2(double tau, double delta)
	{
		double summer=0;
		for (int i=iStart; i<=iEnd; i++){
			summer += a[i]*b[i]*(b[i]-1)*pow(tau,b[i]-2);
		}
		return summer;
	};
	double dTau3(double tau, double delta)
	{
		double summer=0;
		for (int i=iStart; i<=iEnd; i++){
			summer += a[i]*b[i]*(b[i]-1)*(b[i]-2)*pow(tau,b[i]-3);
		}
		return summer;
	};
	double dDelta(double tau, double delta){return 0.0;};
	double dDelta2(double tau, double delta){return 0.0;};
	double dDelta2_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau2(double tau, double delta){return 0.0;};
	double dDelta3(double tau, double delta){return 0.0;};
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
	double dTau3(double tau, double delta)
	{
		return 2*c/(tau*tau*tau);
	};
	double dDelta(double tau, double delta){return 0.0;};
	double dDelta2(double tau, double delta){return 0.0;};
	double dDelta2_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau2(double tau, double delta){return 0.0;};
	double dDelta3(double tau, double delta){return 0.0;};
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
		a=a_; tv=t_; Tc = Tc_; T0=T0_; iStart=iStart_; iEnd=iEnd_; tau0=Tc/T0;
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
	double dTau(double tau, double delta);
	double dTau2(double tau, double delta);
	double dDelta(double tau, double delta){return 0.0;};
	double dDelta2(double tau, double delta){return 0.0;};
	double dDelta2_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau2(double tau, double delta){return 0.0;};
	double dTau3(double tau, double delta);
	double dDelta3(double tau, double delta){return 0.0;};
};

/// Term in the ideal-gas specific heat equation that is based on Aly-Lee formulation
/// All the terms are of the form 
/// cp0 = A + B*((C/T)/sinh(C/T))^2 + D*((E/T)/cosh(E/T))^2
/// Note the LHS is NOT cp0/R
class phi0_cp0_AlyLee : public phi_BC{
private:
	std::vector<double> a;
	double Tc,tau0,T0,R_u; // Use these variables internally
public:

	/// Constructor with std::vectors
	phi0_cp0_AlyLee(std::vector<double> _a, double _Tc, double _T0, double _R) { 
		a=_a; Tc = _Tc; T0=_T0; tau0=Tc/T0; R_u = _R;
	};

	/// Destructor
	~phi0_cp0_AlyLee(){};
	double base(double tau, double delta);
	double dTau(double tau, double delta);
	double dTau2(double tau, double delta);
	double dDelta(double tau, double delta){return 0.0;};
	double dDelta2(double tau, double delta){return 0.0;};
	double dDelta2_dTau(double tau, double delta){return 0.0;};
	double dDelta_dTau(double tau, double delta){return 0.0;};
	double cp0(double tau);
	double anti_deriv_cp0_tau(double tau);
	double anti_deriv_cp0_tau2(double tau);
	double dDelta_dTau2(double tau, double delta){return 0.0;};
	double dTau3(double tau, double delta);
	double dDelta3(double tau, double delta){return 0.0;};
};

#endif