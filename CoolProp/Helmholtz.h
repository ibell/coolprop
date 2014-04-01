#ifndef HELMHOLTZ_H
#define HELMHOLTZ_H

#include <iostream>
#include <vector>
#include "math.h"
#include "CoolPropTools.h"
#include "CPExceptions.h"
#include "rapidjson_CoolProp.h"

#ifdef __ISWINDOWS__
	#define _USE_MATH_DEFINES
	#include "float.h"
#else
	#ifndef DBL_EPSILON
		#include <limits>
		#define DBL_EPSILON std::numeric_limits<double>::epsilon()
	#endif
#endif



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

    virtual void to_json(rapidjson::Value &el, rapidjson::Document &doc) = 0;
};

/// Check the derivatives for a Helmholtz energy term
/// Fluid *fl = get_fluid(get_Fluid_index("Methanol"));
/// check_derivatives(fl->phirlist.at(2),0.5,0.3);
void check_derivatives(phi_BC * phi, double tau, double delta, double ddelta = 1e-10, double dtau = 1e-10);

/*!

Terms are of the form 
\f[
\phi_r = n delta ^d tau^t \exp(-delta^l)
\f]
if l>0 or 
if l==0, then 
\f[
\phi_r = n delta ^d tau^t
\f]

*/
class phir_power : public phi_BC{
	
private:
	
public:
	unsigned int iStart,iEnd;
	std::vector<double> n, ///< The coefficients multiplying each term
		                d, ///< The power for the delta terms
						t, ///< The powers for the tau terms
						l //< The powers for delta in the exp terms
						;
	// Default Constructor
	phir_power(){};
	// Constructors
	phir_power(std::vector<double>,std::vector<double>,std::vector<double>,int,int);
	phir_power(std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,int,int);
	// Add a dummy value at the end of the function since compiler sees std::vector<double> and double[] as being the same type
	phir_power(const double[],const double[],const double[],int,int,int);
	phir_power(double[], double[], double[],int,int,int);
	phir_power(const double[], const double[], const double[], const double[],int,int,int);
	phir_power(double[],double[],double[],double[],int,int,int);
    //std::string to_json();
	
	/// Cache some terms for internal use
	void cache();

	///< Destructor for the phir_power class.  No implementation
	~phir_power(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

	double base(double tau, double delta) throw();
	double dDelta(double tau, double delta) throw();
	double dTau(double tau, double delta) throw();
	
	double dDelta2(double tau, double delta) throw();
	double dDelta_dTau(double tau, double delta) throw();
	double dTau2(double tau, double delta) throw();
	
	double dDelta3(double tau, double delta) throw();
	double dDelta2_dTau(double tau, double delta) throw();
	double dDelta_dTau2(double tau, double delta) throw();
	double dTau3(double tau, double delta) throw();

	/// Vectorized form so iteration happens at c++ level
	std::vector<double> dDeltaV(std::vector<double> tau, std::vector<double> delta) throw();
	std::vector<double> dDelta2V(std::vector<double> tau, std::vector<double> delta) throw();
	std::vector<double> dTau2V(std::vector<double> tau, std::vector<double> delta) throw();
	std::vector<double> dDelta_dTauV(std::vector<double> tau, std::vector<double> delta) throw();

	/// Derivatives for a single term for use in fitter
	double A(double log_tau, double log_delta, double delta, int i) throw();
	double dA_dDelta(double log_tau, double log_delta, double delta, int i) throw();
	double dA_dTau(double log_tau, double log_delta, double delta, int i) throw();
	double d2A_dTau2(double log_tau, double log_delta, double delta, int i) throw();
	double d2A_dDelta2(double log_tau, double log_delta, double delta, int i) throw();
	double d2A_dDelta_dTau(double log_tau, double log_delta, double delta, int i) throw();
};


/*!

Terms are of the form 
\f[
\phi_r = n \delta ^d \tau^t \exp(-\gamma*\delta^l)
\f]

*/
class phir_exponential : public phi_BC{
private:
	std::vector<double> n, ///< The coefficients multiplying each term
		                d, ///< The power for the delta terms
						t, ///< The powers for the tau terms
						l, ///< The powers for delta in the exp terms
						g; ///< Gamma in the exponential term
	unsigned int iStart,iEnd;
public:
	// Constructors
	phir_exponential(std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,int,int);
	phir_exponential(const double[], const double[], const double[], const double[],const double[],int,int,int);
	phir_exponential(double[],double[],double[],double[],double[],int,int,int);
	
	///< Destructor for the phir_power class.  No implementation
	~phir_exponential(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

	double base(double tau, double delta) throw();
	double dDelta(double tau, double delta) throw();
	double dTau(double tau, double delta) throw();
	
	double dDelta2(double tau, double delta) throw();
	double dDelta_dTau(double tau, double delta) throw();
	double dTau2(double tau, double delta) throw();
	
	double dDelta3(double tau, double delta) throw();
	double dDelta2_dTau(double tau, double delta) throw();
	double dDelta_dTau2(double tau, double delta) throw();
	double dTau3(double tau, double delta) throw();
};

/*!

Terms are of the form 
\f[
\phi_r = a \delta ^d \tau^t \exp(-\alpha(\delta-\epsilon)^2-\beta(\tau-\gamma)^2)
\f]

*/
class phir_gaussian : public phi_BC{
private:
	std::vector<double> n,d,t,alpha,epsilon,beta,gamma;
	unsigned int iStart,iEnd;
public:
	// Default Constructor
	phir_gaussian(){};
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
	phir_gaussian(const double a_in[],	
				  const double d_in[],
				  const double t_in[], 
				  const double alpha_in[], 
				  const double epsilon_in[], 
				  const double beta_in[], 
				  const double gamma_in[],
		unsigned int iStart_in, unsigned int iEnd_in, unsigned int N);

	// Destructor
	~phir_gaussian(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

	// Term and its derivatives
	double base(double tau, double delta) throw();
	double dDelta(double tau, double delta) throw();
	double dTau(double tau, double delta) throw();
	
	double dDelta2(double tau, double delta) throw();
	double dDelta_dTau(double tau, double delta) throw();
	double dTau2(double tau, double delta) throw();
	
	double dDelta3(double tau, double delta) throw();
	double dDelta2_dTau(double tau, double delta) throw();
	double dDelta_dTau2(double tau, double delta) throw();
	double dTau3(double tau, double delta) throw();
	
};

/*!
The Gaussian term from the GERG 2008 mixture formulation
\f[
\phi_r = a \delta ^d \tau^t \exp(-\eta(\delta-\epsilon)^2-\beta(\delta-\gamma))
\f]

*/
class phir_GERG2008_gaussian : public phi_BC{
private:
	std::vector<double> n,d,t,eta,epsilon,beta,gamma;
	unsigned int iStart,iEnd;
public:
	// Default Constructor
	phir_GERG2008_gaussian(){};
	// Constructors
	phir_GERG2008_gaussian(std::vector<double> a_in, 
				  std::vector<double> d_in,
				  std::vector<double> t_in, 
				  std::vector<double> eta_in, 
				  std::vector<double> epsilon_in, 
				  std::vector<double> beta_in, 
				  std::vector<double> gamma_in,
				  unsigned int iStart_in, 
				  unsigned int iEnd_in);
	phir_GERG2008_gaussian(double a_in[], 
					   double d_in[],
					   double t_in[], 
					   double eta_in[], 
					   double epsilon_in[], 
					   double beta_in[], 
					   double gamma_in[],
					   unsigned int iStart_in, 
					   unsigned int iEnd_in, 
					   unsigned int N);
	// Destructor
	~phir_GERG2008_gaussian(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

	// Term and its derivatives
	double base(double tau, double delta);
	double dDelta(double tau, double delta);
	double dTau(double tau, double delta);
	
	double dDelta2(double tau, double delta);
	double dDelta_dTau(double tau, double delta);
	double dTau2(double tau, double delta);
	
	double dDelta3(double tau, double delta){throw ValueError();};
	double dDelta2_dTau(double tau, double delta){throw ValueError();};
	double dDelta_dTau2(double tau, double delta){throw ValueError();};
	double dTau3(double tau, double delta){throw ValueError();};
};

/*!

The critical term, used for Water and Carbon Dioxide
It is truly horrible to implement

*/
class phir_critical : public phi_BC{
	
private:
	std::vector<double> n,d,t,a,b,A,B,C,D,beta;
	int iStart,iEnd;
public:
	// Constructors
	phir_critical(std::vector<double> n_in, std::vector<double> d_in, std::vector<double> t_in, 
		std::vector<double> a_in, std::vector<double> b_in, std::vector<double> beta_in,
		std::vector<double> A_in, std::vector<double> B_in, std::vector<double> C_in, 
		std::vector<double> D_in, int iStart_in, int iEnd_in);

	phir_critical(double n[], double d[], double t[], 
				  double a[], double b[], double beta[],
				  double A[], double B[], double C[], 
				  double D[], int iStart, int iEnd, int N);

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

	//Destructor
	~phir_critical(){};

	// Term and its derivatives
	double base(double tau, double delta) throw();
	double dDelta(double tau, double delta) throw();
	double dTau(double tau, double delta) throw();
	
	double dDelta2(double tau, double delta) throw();
	double dDelta_dTau(double tau, double delta) throw();
	double dTau2(double tau, double delta) throw();
	
	double dDelta3(double tau, double delta) throw();
	double dDelta2_dTau(double tau, double delta) throw();
	double dDelta_dTau2(double tau, double delta) throw();
	double dTau3(double tau, double delta) throw();
};

/*! 

This class implements residual Helmholtz Energy terms of the form  

\f[
\phi_r = n  \delta^d  \tau^t  \exp(-\delta^l) \exp(-\tau^m)
\f]

*/

class phir_Lemmon2005 : public phi_BC{
	/*
	Terms are of the form 
	n * delta ^d * tau^t if l == 0 and m == 0
	n * delta ^d * tau^t * exp(-delta^l) if l != 0  and m = 0
	n * delta ^d * tau^t * exp(-delta^l) * exp(-tau^m) if l != 0 and m != 0
	Constructor must be called with std::vector instances of double type
	*/
private:
	std::vector<double> n, ///< The coefficients multiplying each term
		                d, ///< The power for the delta terms
						t, ///< The powers for the tau terms
						l, ///< The powers for delta in the exp terms
						m; ///< The powers for tau in the exp terms
	unsigned int iStart,iEnd;
public:
	// Constructors
	phir_Lemmon2005(std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,int,int);
	// Add a dummy value at the end of the function since compiler sees std::vector<double> and double[] as being the same type
	phir_Lemmon2005(const double[], const double[], const double[], const double[],const double[],int,int,int);
	phir_Lemmon2005(double[],double[],double[],double[],double[],int,int,int);

	///< Destructor for the phir_Lemmon2005 class.  No implementation
	~phir_Lemmon2005(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

	double base(double tau, double delta) throw();
	double dDelta(double tau, double delta) throw();
	double dTau(double tau, double delta) throw();
	
	double dDelta2(double tau, double delta) throw();
	double dDelta_dTau(double tau, double delta) throw();
	double dTau2(double tau, double delta) throw();
	
	double dDelta3(double tau, double delta) throw();
	double dDelta2_dTau(double tau, double delta) throw();
	double dDelta_dTau2(double tau, double delta) throw();
	double dTau3(double tau, double delta) throw();
};

class phir_SAFT_associating : public phi_BC{
	
protected:
	double m,epsilonbar, vbarn, kappabar,a;
public:
	// Constructor
	phir_SAFT_associating(){};

	//Destructor
	~phir_SAFT_associating(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    std::string to_json();

	double Deltabar(double tau, double delta);
	double dDeltabar_ddelta__consttau(double tau, double delta);
	double d2Deltabar_ddelta2__consttau(double tau, double delta);
	double dDeltabar_dtau__constdelta(double tau, double delta);
	double d2Deltabar_dtau2__constdelta(double tau, double delta);
	double d2Deltabar_ddelta_dtau(double tau, double delta);
    double d3Deltabar_dtau3__constdelta(double tau, double delta);
    double d3Deltabar_ddelta_dtau2(double tau, double delta);
    double d3Deltabar_ddelta3__consttau(double tau, double delta);
    double d3Deltabar_ddelta2_dtau(double tau, double delta);

	double X(double delta, double Deltabar);
	double dX_dDeltabar__constdelta(double delta, double Deltabar);
	double dX_ddelta__constDeltabar(double delta, double Deltabar);
	double dX_dtau(double tau, double delta);
	double dX_ddelta(double tau, double delta);
	double d2X_dtau2(double tau, double delta);
	double d2X_ddeltadtau(double tau, double delta);
	double d2X_ddelta2(double tau, double delta);

    double d3X_dtau3(double tau, double delta);
    double d3X_ddelta3(double tau, double delta);
    double d3X_ddeltadtau2(double tau, double delta);
    double d3X_ddelta2dtau(double tau, double delta);

	double g(double eta);
	double dg_deta(double eta);
	double d2g_deta2(double eta);   
	double d3g_deta3(double eta);
	double eta(double delta);

    double base(double tau, double delta);
	double dDelta(double tau, double delta);
	double dTau(double tau, double delta);
	
	double dDelta2(double tau, double delta);
	double dDelta_dTau(double tau, double delta);
	double dTau2(double tau, double delta);
	
	double dDelta3(double tau, double delta);
	double dDelta2_dTau(double tau, double delta);
    double dDelta_dTau2(double tau, double delta);
	double dTau3(double tau, double delta);
};

class phir_SAFT_associating_2B : public phir_SAFT_associating{

public:
	phir_SAFT_associating_2B(double m, double epsilonbar, double vbarn, double kappabar)
	{
        this->a = 2;
		this->m = m;
		this->epsilonbar = epsilonbar;
		this->vbarn = vbarn;
		this->kappabar = kappabar;
	};
};
class phir_SAFT_associating_1 : public phir_SAFT_associating{

public:
	phir_SAFT_associating_1(double m, double epsilonbar, double vbarn, double kappabar)
	{
        this->a = 1;
		this->m = m;
		this->epsilonbar = epsilonbar;
		this->vbarn = vbarn;
		this->kappabar = kappabar;
	};
};

/*!

\f[
\phi_0 = \log(\delta)+a_1+a_2\tau
\f]
*/
class phi0_lead : public phi_BC{
	/*
	constructor: phi0_lead(double a1, double a2)
	*/
private:
	double c1,c2; // Use these variables internally
public:
	// Constructor
	phi0_lead(double a1, double a2){c1=a1; c2=a2;};

	//Destructor
	~phi0_lead(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc){
        el.AddMember("type","IdealGasHelmholtzLead",doc.GetAllocator());
        el.AddMember("a1",c1,doc.GetAllocator());
        el.AddMember("a2",c2,doc.GetAllocator());
    };

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

/*!
\f[
\phi_0 = a_1+a_2\tau
\f]

constructor: phi0_enthalpy_entropy_offset(double a1, double a2)
*/
class phi0_enthalpy_entropy_offset : public phi_BC{
private:
	double c1,c2; // Use these variables internally
public:
	// Constructor
	phi0_enthalpy_entropy_offset(double a1, double a2){c1=a1; c2=a2;};

	//Destructor
	~phi0_enthalpy_entropy_offset(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc){
        el.AddMember("type","IdealGasEnthalpyEntropyOffset",doc.GetAllocator());
        el.AddMember("a1",c1,doc.GetAllocator());
        el.AddMember("a2",c2,doc.GetAllocator());
    };

	// Term and its derivatives
	double base(double tau, double delta){return c1+c2*tau;};
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


/*!
Term is of the form 
\f[
\phi_0 = a_1 \log(\tau)
\f]
*/
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

    void to_json(rapidjson::Value &el, rapidjson::Document &doc){
        el.AddMember("type","IdealGasHelmholtzLogTau",doc.GetAllocator());
        el.AddMember("a",c1,doc.GetAllocator());
    };

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
	phi0_Planck_Einstein(double a_in[], double theta_in[], int iStart_in, int iEnd_in, int N)
	{
		a=std::vector<double>(a_in,a_in+N);
		theta=std::vector<double>(theta_in,theta_in+N);
		iStart = iStart_in; iEnd = iEnd_in;
	};
	phi0_Planck_Einstein(const double a_in[], const double theta_in[], int iStart_in, int iEnd_in, int N)
	{
		a=std::vector<double>(a_in,a_in+N);
		theta=std::vector<double>(theta_in,theta_in+N);
		iStart = iStart_in; iEnd = iEnd_in;
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
  
    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

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
	phi0_Planck_Einstein2(double a, double theta, double c)
	{
		this->a=std::vector<double> (1,a); 
		this->theta=std::vector<double> (1,theta); 
		this->c=std::vector<double> (1,c); 
		iStart = 0; iEnd = 0;
	};

	//Destructor
	~phi0_Planck_Einstein2(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

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

/*
Term is of the form a*tau^b
Constructor:
	phi0_power(std::vector<double> a, std::vector<double> b, int iStart, int iEnd)
*/
class phi0_power : public phi_BC{
	
private:
	std::vector<double> a,b; // Use these variables internally
	int iStart, iEnd;
public:
	// Constructor
	phi0_power(std::vector<double> a, std::vector<double> b, int iStart, int iEnd)
	{
		this->a=a; 
		this->b=b; 
		this->iStart = iStart; 
		this->iEnd = iEnd;
	};
	phi0_power(double a[], double b[], int iStart, int iEnd, int N)
	{
		this->a = std::vector<double>(a,a+N);
		this->b = std::vector<double>(b,b+N);
		this->iStart = iStart; 
		this->iEnd = iEnd;
	};
	phi0_power(const double a[], const double b[], int iStart, int iEnd, int N)
	{
		this->a = std::vector<double>(a,a+N);
		this->b = std::vector<double>(b,b+N);
		this->iStart = iStart; 
		this->iEnd = iEnd;
	};
	// Constructor with just double values
	phi0_power(double a, double b)
	{
		this->a=std::vector<double>(1,a); 
		this->b=std::vector<double>(1,b); 
		iStart = 0;
		iEnd = 0;
	};

	//Destructor
	~phi0_power(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc)
    {
        el.AddMember("type","IdealGasHelmholtzPower",doc.GetAllocator());
        rapidjson::Value _n(rapidjson::kArrayType),_t(rapidjson::kArrayType);
        for (int i=iStart;i<=iEnd;i++)
	    {
            _n.PushBack(a[i],doc.GetAllocator());
            _t.PushBack(b[i],doc.GetAllocator());
	    }
        el.AddMember("n",_n,doc.GetAllocator());
        el.AddMember("t",_t,doc.GetAllocator());
    };

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
	phi0_cp0_constant(double c, double Tc, double T0) { this->c=c; this->T0=T0; this->Tc=Tc; this->tau0=Tc/T0;};

	/// Destructor
	~phi0_cp0_constant(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc)
    {
        el.AddMember("type","IdealGasHelmholtzCP0Constant",doc.GetAllocator());
        el.AddMember("cp_over_R",c,doc.GetAllocator());
        el.AddMember("Tc",Tc,doc.GetAllocator());
        el.AddMember("T0",T0,doc.GetAllocator());
    };

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

/*!

for a term of the form
\f[
    \frac{c_p^0}{R}=cT^t, t \neq 0,-1
\f]
the contribution is found from 
\f[
    \frac{1}{T}\int_{T_0}^T c T^t dT-\int_{T_0}^T \frac{c T^t}{T}dT
\f]    
\f[
    \frac{c}{T}\left(\frac{T^{t+1}}{t+1}-\frac{T_0^{t+1}}{t+1}\right)-c\left(\frac{T^{t}}{t}-\frac{T_0^{t}}{t}\right)
\f]
\f[
    cT^{t}\left(\frac{1}{t+1}-\frac{1}{t}\right)-c\frac{T_0^{t+1}}{T(t+1)}+c\frac{T_0^t}{t}
\f]
or in terms of \$$\tau\$$ 
\f[
    cT_c^{t}\tau^{-t}\left(\frac{1}{t+1}-\frac{1}{t}\right)-c\frac{T_0^{t+1}\tau}{T_c(t+1)}+c\frac{T_0^t}{t}
\f]
if t = 0
\f[
\frac{1}{T}\int_{{T_0}}^T c dT - \int_{{T_0}}^T {\frac{c}{T}} dT = \frac{{c(T - {T_0})}}{T} - c\ln \left( {\frac{T}{{{T_0}}}} \right) = c\left( {1 - \frac{\tau }{{{\tau _0}}}} \right) - c\ln \left( {\frac{{{\tau _0}}}{\tau }} \right)
\f]
if t = -1
\f[
\frac{1}{T}\int_{{T_0}}^T {\frac{c}{T}} dT - \int_{{T_0}}^T {\frac{c}{{{T^2}}}} dT = \frac{c}{T}\ln \left( {\frac{T}{{{T_0}}}} \right) + c\left( {\frac{1}{T} - \frac{1}{{{T_0}}}} \right) = \frac{{c\tau }}{{{T_c}}}\ln \left( {\frac{{{\tau _0}}}{\tau }} \right) + \frac{c}{{{T_c}}}\left( {\tau  - {\tau _0}} \right)
\f]

*/
class phi0_cp0_poly : public phi_BC{
private:
	std::vector<double> a,tv;
	double Tc,T0,tau0; // Use these variables internally
	int iStart, iEnd;
public:
	/// Constructor with just a single double value
	phi0_cp0_poly(double a, double t, double Tc, double T0) {
		this->a=std::vector<double>(1,a);
		this->tv=std::vector<double>(1,t);
		this->Tc=Tc; this->T0=T0; iStart=0; iEnd=0; tau0=Tc/T0;
	};

	/// Constructor with std::vectors
	phi0_cp0_poly(std::vector<double> a, std::vector<double> t, double Tc, double T0, int iStart, int iEnd) { 
		this->a=a; this->tv=t; this->Tc = Tc; this->T0=T0; this->iStart=iStart; this->iEnd=iEnd; tau0=Tc/T0;
	};

	/// Destructor
	~phi0_cp0_poly(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

	// Term and its derivatives
	double base(double tau, double delta){ 
		double sum=0;
		for (int i = iStart; i<=iEnd; i++){
			double t=tv[i];
			if (fabs(t)<10*DBL_EPSILON)
			{
				sum += a[i]-a[i]*tau/tau0+a[i]*log(tau/tau0);
			}
			else if (fabs(t+1) < 10*DBL_EPSILON)
			{
				sum += a[i]*tau/Tc*log(tau0/tau)+a[i]/Tc*(tau-tau0);
			}
			else
			{
				sum += -a[i]*pow(Tc,t)*pow(tau,-t)/(t*(t+1))-a[i]*pow(T0,t+1)*tau/(Tc*(t+1))+a[i]*pow(T0,t)/t;
			}
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

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

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
