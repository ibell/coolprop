#ifndef HELMHOLTZ_H
#define HELMHOLTZ_H

#include <iostream>
#include <vector>
#include "math.h"

// The abstract base class upon which each residual Helmholtz energy class is built
class phi_BC{
public:
	phi_BC(){};
	virtual ~phi_BC(){};
	virtual double base(double,double) = 0;
	virtual double dTau(double,double) = 0;
	virtual double dTau2(double,double) = 0;
	virtual double dDelta_dTau(double,double) = 0;
	virtual double dDelta(double,double) = 0;
	virtual double dDelta2(double,double) = 0;
	virtual double dDelta2_dTau(double,double) = 0;
};

class phir_power : public phi_BC{
	/*
	Terms are of the form n * delta ^d * tau^t

	Constructor must be called with std::vector instances of double type
	*/
private:
	std::vector<double> n,d,t,l;
	int iStart,iEnd;
public:
	// Constructors
	phir_power(std::vector<double>,std::vector<double>,std::vector<double>,int,int);
	phir_power(std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,int,int);

	// Destructor
	~phir_power(){};

	// Term and its derivatives
	double base(double tau, double delta);
	double dTau(double tau, double delta);
	double dTau2(double tau, double delta);
	double dDelta(double tau, double delta);
	double dDelta2(double tau, double delta);
	double dDelta2_dTau(double tau, double delta);
	double dDelta_dTau(double tau, double delta);
};

class phir_gaussian : public phi_BC{
	/*
	Terms are of the form n * delta ^d * tau^t * exp(-alpha*(delta-epsilon)^2-beta*(tau-gamma)^2)

	Constructor must be called with std::vector instances of double type
	*/
private:
	std::vector<double> n,d,t,alpha,epsilon,beta,gamma;
	int iStart,iEnd;
public:
	// Constructors
	phir_gaussian(std::vector<double> a_in, std::vector<double> d_in,std::vector<double> t_in, 
		std::vector<double> alpha_in, std::vector<double> epsilon_in, std::vector<double> beta_in, std::vector<double> gamma_in,int iStart_in, int iEnd_in);

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
		for (int i=iStart;i<=iEnd;i++) summer+=a[i]*b[i]*pow(tau,b[i]-1);
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

#endif

