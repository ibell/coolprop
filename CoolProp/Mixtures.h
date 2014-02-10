#ifndef MIXTURES_H
#define MIXTURES_H

#include "Helmholtz.h"
#include "FluidClass.h"

#include <vector>
#include <string>

typedef std::vector<std::vector<double> > STLMatrix;

/*! 
An abstract base class for the reducing function to allow for
Lemmon-Jacobsen, GERG, or other reducing function to yield the
reducing parameters \f$ \bar\rho_r \f$ and \f$ T_r \f$
*/
class ReducingFunction
{
protected:
	unsigned int N;
public:
	ReducingFunction(){};
	virtual ~ReducingFunction(){};
	/// The reduced temperature
	virtual double Tr(const std::vector<double> &x) = 0;
	/// The derivative of reduced temperature with respect to component i mole fraction
	virtual double dTrdxi__constxj(const std::vector<double> &x, int i) = 0;
	/// The molar reducing density
	virtual double rhorbar(const std::vector<double> &x) = 0;
	///Derivative of the molar reducing density with respect to component i mole fraction
	virtual double drhorbardxi__constxj(const std::vector<double> &x, int i) = 0;

	/// Set the coefficients based on reducing parameters loaded from JSON
	virtual void set_coeffs_from_map(int i, int j, std::map<std::string,double >) = 0;

	virtual double d2rhorbardxi2__constxj(const std::vector<double> &x, int i) = 0;
	virtual double d2rhorbardxidxj(const std::vector<double> &x, int i, int j) = 0;
	virtual double d2Trdxi2__constxj(const std::vector<double> &x, int i) = 0;
	virtual double d2Trdxidxj(const std::vector<double> &x, int i, int j) = 0;

	/*! GERG 2004 Monograph equation 7.56:
	\f[
	\left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial T_r}{\partial n_i} \right)_{n_j}\right)\right)_{x_i} = \left(\frac{\partial^2T_r}{\partial x_j \partial x_i}\right)-\left(\frac{\partial T_r}{\partial x_j}\right)_{x_i}-\sum_{k=1}^Nx_k\left(\frac{\partial^2T_r}{\partial x_j \partial x_k}\right)
	\f]
	*/
	double d_ndTrdni_dxj__constxi(const std::vector<double> &x, int i, int j);
	/*! GERG 2004 Monograph equation 7.55:
	\f[
	\left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial \rho_r}{\partial n_i} \right)_{n_j}\right)\right)_{x_i} = \left(\frac{\partial^2\rho_r}{\partial x_j \partial x_i}\right)-\left(\frac{\partial \rho_r}{\partial x_j}\right)_{x_i}-\sum_{k=1}^Nx_k\left(\frac{\partial^2\rho_r}{\partial x_j \partial x_k}\right)
	\f]
	*/
	double d_ndrhorbardni_dxj__constxi(const std::vector<double> &x, int i, int j);

	double ndrhorbardni__constnj(const std::vector<double> &x, int i);
	double ndTrdni__constnj(const std::vector<double> &x, int i);
};

/*! 
The Reducing parameter model used by the GERG-2008 formulation to yield the
reducing parameters \f$ \bar\rho_r \f$ and \f$ T_r \f$ and derivatives thereof
*/
class GERG2008ReducingFunction : public ReducingFunction
{
protected:
	STLMatrix v_c;
	STLMatrix T_c; //!< \f$ 
	STLMatrix beta_v; //!< \f$ \beta_{v,ij} \f$ from GERG-2008
	STLMatrix gamma_v; //!< \f$ \gamma_{v,ij} \f$ from GERG-2008
	STLMatrix beta_T; //!< \f$ \beta_{T,ij} \f$ from GERG-2008
	STLMatrix gamma_T; //!< \f$ \gamma_{T,ij} \f$ from GERG-2008
	std::vector<Fluid *> pFluids; //!< List of pointers to fluids
public:
	GERG2008ReducingFunction(std::vector<Fluid *> pFluids, STLMatrix beta_v, STLMatrix gamma_v, STLMatrix beta_T, STLMatrix gamma_T)
	{
		this->pFluids = pFluids;
		this->beta_v = beta_v;
		this->gamma_v = gamma_v;
		this->beta_T = beta_T;
		this->gamma_T = gamma_T;
		this->N = pFluids.size();
	};
	GERG2008ReducingFunction(std::vector<Fluid *> pFluids)
	{
		this->pFluids = pFluids;
		this->N = pFluids.size();
		/// Resize reducing parameter matrices to be the same size as x in both directions
		beta_v.resize(N,std::vector<double>(N,0));
		gamma_v.resize(N,std::vector<double>(N,0));
		beta_T.resize(N,std::vector<double>(N,0));
		gamma_T.resize(N,std::vector<double>(N,0));
		T_c.resize(N,std::vector<double>(N,0));
		v_c.resize(N,std::vector<double>(N,0));
		for (unsigned int i = 0; i < N; i++)
		{
			for (unsigned int j = 0; j < N; j++)
			{
				T_c[i][j] = sqrt(pFluids[i]->reduce.T*pFluids[j]->reduce.T);
				v_c[i][j] = 1.0/8.0*pow(pow(pFluids[i]->reduce.rhobar, -1.0/3.0)+pow(pFluids[j]->reduce.rhobar, -1.0/3.0),(int)3);
			}
		}
	};
	/// Default destructor
	~GERG2008ReducingFunction(){};
	/// The reduced temperature
	double Tr(const std::vector<double> &x);
	/// The derivative of reduced temperature with respect to component i mole fraction
	double dTrdxi__constxj(const std::vector<double> &x, int i);
	/// The molar reducing density
	double rhorbar(const std::vector<double> &x);
	///Derivative of the molar reducing density with respect to component i mole fraction
	double drhorbardxi__constxj(const std::vector<double> &x, int i);
	double dvrbardxi__constxj(const std::vector<double> &x, int i);

	double d2vrbardxi2__constxj(const std::vector<double> &x, int i);
	double d2rhorbardxi2__constxj(const std::vector<double> &x, int i);
	double d2vrbardxidxj(const std::vector<double> &x, int i, int j);
	double d2rhorbardxidxj(const std::vector<double> &x, int i, int j);
	double d2Trdxi2__constxj(const std::vector<double> &x, int i);
	double d2Trdxidxj(const std::vector<double> &x, int i, int j);

	/// Set the coefficients based on reducing parameters loaded from JSON
	void set_coeffs_from_map(int i, int j, std::map<std::string,double >);

	double c_Y_ij(int i, int j, std::vector< std::vector< double> > * beta, std::vector< std::vector< double> > *gamma, std::vector< std::vector< double> > *Y_c);
	double c_Y_ji(int j, int i, std::vector< std::vector< double> > * beta, std::vector< std::vector< double> > *gamma, std::vector< std::vector< double> > *Y_c);
	double f_Y_ij(const std::vector<double> &x, int i, int j, std::vector< std::vector< double> > * beta);

	double dfYkidxi__constxk(const std::vector<double> &x, int k, int i,std::vector< std::vector< double> > * beta);
	double dfYikdxi__constxk(const std::vector<double> &x, int i, int k, std::vector< std::vector< double> > * beta);
	double d2fYkidxi2__constxk(const std::vector<double> &x, int k, int i, std::vector< std::vector< double> > * beta);
	double d2fYikdxi2__constxk(const std::vector<double> &x, int i, int k, std::vector< std::vector< double> > * beta);
	double d2fYijdxidxj(const std::vector<double> &x, int i, int k, std::vector< std::vector< double> > * beta);
};

/*! From Lemmon, JPCRD, 2000 for the properties of Dry Air, and also from Lemmon, JPCRD, 2004 for the properties of R404A, R410A, etc.	
\f[
\rho_r(\bar x) = \left[ \sum_{i=1}^m\frac{x_i}{\rho_{c_i}}+\sum_{i=1}^{m-1}\sum_{j=i+1}^{m}x_ix_j\zeta_{ij}\right]^{-1}
\f]
\f[
T_r(\bar x) = \sum_{i=1}^mx_iT_{c_i}+\sum_{i=1}^{m-1}\sum_{j=i+1}^mx_ix_j\xi_{ij}
\f]

These can be converted to the form of GERG by the following equations:
\f[
\beta_T = 1\ \ \ \ \beta_v = 1 
\f]
and
\f[
    \boxed{\gamma_T = \dfrac{T_{c0}+T_{c1}+\xi_{01}}{2\sqrt{T_{c0}T_{c1}}}}
\f]
and
\f[
    \boxed{\gamma_v = \dfrac{v_{c0}+v_{c1}+\zeta_{01}}{\frac{1}{4}\left(\frac{1}{\rho_{c,i}^{1/3}}+\frac{1}{\rho_{c,j}^{1/3}}\right)^{3}}}
\f]
*/
class LemmonAirHFCReducingFunction : public GERG2008ReducingFunction
{
protected:
	std::vector<Fluid *> pFluids; //!< List of pointers to fluids	
public:
	LemmonAirHFCReducingFunction(std::vector<Fluid *> pFluids) : GERG2008ReducingFunction(pFluids)
	{
		this->pFluids = pFluids;
		this->N = pFluids.size();
	};

	/// Set the coefficients based on reducing parameters loaded from JSON
	void set_coeffs_from_map(int i, int j, std::map<std::string,double >);
	
};








/*! 
The abstract base class for departure functions for the excess part of the Helmholtz energy
*/
class DepartureFunction
{
public:
	DepartureFunction(){};
	virtual ~DepartureFunction(){};
	
	/// The excess Helmholtz energy of the binary pair
	/// Pure-virtual function (must be implemented in derived class
	virtual double phir(double tau, double delta) = 0;
	virtual double dphir_dDelta(double tau, double delta) = 0;
	virtual double d2phir_dDelta2(double tau, double delta) = 0;
	virtual double d2phir_dDelta_dTau(double tau, double delta) = 0;
	virtual double dphir_dTau(double tau, double delta) = 0;
	virtual double d2phir_dTau2(double tau, double delta) = 0;
	virtual void set_coeffs_from_map(std::map<std::string,std::vector<double> >) = 0;
};

struct DepartureFunctionCacheElement
{
	double tau, delta, cached_val;
};

struct DepartureFunctionCache
{
	DepartureFunctionCacheElement phir,dphir_dDelta,dphir_dTau,d2phir_dDelta2,d2phir_dTau2,d2phir_dDelta_dTau;
};

class GERG2008DepartureFunction : public DepartureFunction
{
protected:
	bool using_gaussian;
	DepartureFunctionCache cache;
	phir_power phi1;
	phir_GERG2008_gaussian phi2;
public:
	GERG2008DepartureFunction(){};
	~GERG2008DepartureFunction(){};
	double phir(double tau, double delta);
	double dphir_dDelta(double tau, double delta);
	double d2phir_dDelta_dTau(double tau, double delta);
	double dphir_dTau(double tau, double delta);
	double d2phir_dDelta2(double tau, double delta);
	double d2phir_dTau2(double tau, double delta);
	void set_coeffs_from_map(std::map<std::string,std::vector<double> >);
};

class LemmonHFCDepartureFunction : public DepartureFunction
{
protected:
	bool using_gaussian;
	DepartureFunctionCache cache;
	phir_power phi1;
public:
	LemmonHFCDepartureFunction(){};
	~LemmonHFCDepartureFunction(){};
	double phir(double tau, double delta);
	double dphir_dDelta(double tau, double delta);
	double d2phir_dDelta_dTau(double tau, double delta);
	double dphir_dTau(double tau, double delta);
	double d2phir_dDelta2(double tau, double delta);
	double d2phir_dTau2(double tau, double delta);
	void set_coeffs_from_map(std::map<std::string,std::vector<double> >);
};

class ExcessTerm
{
public:
	unsigned int N;
	std::vector<std::vector<DepartureFunction*> > DepartureFunctionMatrix;
	std::vector<std::vector<double> > F;
	ExcessTerm(int N);
	~ExcessTerm();
	double phir(double tau, double delta, const std::vector<double> &x);
	double dphir_dDelta(double tau, double delta, const std::vector<double> &x);
	double d2phir_dDelta2(double tau, double delta, const std::vector<double> &x);
	double d2phir_dDelta_dTau(double tau, double delta, const std::vector<double> &x);
	double dphir_dTau(double tau, double delta, const std::vector<double> &x);
	double d2phir_dTau2(double tau, double delta, const std::vector<double> &x);
	double dphir_dxi(double tau, double delta, const std::vector<double> &x, unsigned int i);
	double d2phirdxidxj(double tau, double delta, const std::vector<double> &x, unsigned int i, unsigned int j);
	double d2phir_dxi_dTau(double tau, double delta, const std::vector<double> &x, unsigned int i);
	double d2phir_dxi_dDelta(double tau, double delta, const std::vector<double> &x, unsigned int i);
	void set_coeffs_from_map(int i, int j, std::map<std::string,std::vector<double> >);
};


class ResidualIdealMixture
{
protected:
	unsigned int N;
	std::vector<Fluid*> pFluids;
public:
	ResidualIdealMixture(std::vector<Fluid*> pFluids);
	double phir(double tau, double delta, const std::vector<double> &x);
	double dphir_dDelta(double tau, double delta, const std::vector<double> &x);
	double d2phir_dDelta2(double tau, double delta, const std::vector<double> &x);
	double d2phir_dDelta_dTau(double tau, double delta, const std::vector<double> &x);
	double d2phir_dTau2(double tau, double delta, const std::vector<double> &x);
	double dphir_dTau(double tau, double delta, const std::vector<double> &x);
};

class Mixture;  // Forward declaration since some classes that are members of Mixture take pointers to Mixture

/*!
A structure to hold information at one iteration of successive substitution for later checking, plotting, etc.
By default this is not used, but it can be enabled.
*/
struct SuccessiveSubstitutionStep
{
	double T, ///< Temperature [K]
		   p; ///< Pressure [kPa]
	std::vector<double> K, ///< K-factor [-]
						lnK, ///< Natural logarithm of K-factor [-]
						ln_phi_liq, ///< Natural logarithm of fugacity coeff. in liq. phase [-]
						ln_phi_vap; ///< Natural logarithm of fugacity coeff. in vapor phase [-]
};

/*!
A class to do successive substitution given guess values for vapor-liquid equilibria.  This class will then be included in the Mixture class

A class is used rather than a function so that it is easier to store iteration histories, additional output values, etc.
*/
class SuccessiveSubstitutionVLE
{
public:
	bool useNR; ///< If true (default is false), will call Newton-Raphson after either solving to tolerance, or taking Nmax steps
	bool logging; ///< If true (default is false), intermediate steps will be stored for each call in the step_logger structure
	int Nsteps; ///< How many steps were taken to yield tolerance
	int Nstep_max; ///< The maximum number of steps to take, can be changed as needed
	Mixture *Mix; ///< Pointer to the Mixture class instance
	double rhobar_liq, ///< The molar density of the liquid phase [mol/m^3]
		   rhobar_vap; ///< The molar density of the vapor phase [mol/m^3]
	std::vector<double> K, 
					    ln_phi_liq, 
						ln_phi_vap, 
						x, 
						y;
	std::vector<SuccessiveSubstitutionStep> step_logger;

	SuccessiveSubstitutionVLE(){useNR = false; logging = false; Nstep_max = 10;};

	double call(double beta, double T, double p, const std::vector<double> &z, std::vector<double> &K);
};

/*!
A class to do newton raphson solver for VLE given guess values for vapor-liquid equilibria.  This class will then be included in the Mixture class

A class is used rather than a function so that it is easier to store iteration histories, additional output values, etc.
*/
class NewtonRaphsonVLE
{
public:
	double error_rms, rhobar_liq, rhobar_vap, T, p;
	unsigned int N;
	bool logging;
	int Nsteps;
	int Nsteps_max;
	Mixture *Mix;
	STLMatrix J;
	std::vector<double> K, x, y, phi_ij_liq, phi_ij_vap, r, dXdS, neg_dFdS;
	std::vector<SuccessiveSubstitutionStep> step_logger;

	NewtonRaphsonVLE(){Nsteps_max = 20;};

	void resize(unsigned int N);
	
	// Reset the state of all the internal variables
	void pre_call()
	{
		K.clear();
		x.clear();
		y.clear();
		phi_ij_liq.clear();
		phi_ij_vap.clear();
		Nsteps = 0;
		step_logger.clear();
		error_rms = 1e99;
		rhobar_liq = _HUGE;
		rhobar_vap = _HUGE;
		T = _HUGE;
		p = _HUGE;
	};

	/*! Call the Newton-Raphson VLE Solver

	This solver must be passed reasonable guess values for the mole fractions, 
	densities, etc.  You may want to take a few steps of successive substitution
	before you start with Newton Raphson.

	@param T Temperature [K]
	@param p Pressure [Pa]
	@param rhobar_liq Current value of molar density of the liquid [mol/m^3]
	@param rhobar_vap Current value of molar density of the vapor [mol/m^3]
	@param z Bulk mole fractions [-]
	@param K Array of K-factors [-]
	*/
	double call(double beta, double T, double p, double rhobar_liq, double rhobar_vap, const std::vector<double> &z, std::vector<double> &K, int spec_index, double spec_value);

	/*! Build the arrays for the Newton-Raphson solve

	This method builds the Jacobian matrix, the sensitivity matrix, etc.

	@param beta Void fraction [-] (0: bubble, 1: dew)
	@param T Temperature [K]
	@param p Pressure [Pa]
	@param z Bulk mole fractions [-]
	@param K Array of K-factors [-]
	*/
	void build_arrays(double beta, double T, double p, const std::vector<double> &z, std::vector<double> & K, int spec_index, double spec_value);
};

struct PhaseEnvelopeLog
{
	std::vector< std::vector<double> > K, lnK, x, y;
	std::vector<double> T, p, lnT, lnp, rhobar_liq, rhobar_vap, lnrhobar_liq, lnrhobar_vap;
	std::vector<int> iS;

	void store_variables(const double T, 
		                 const double p, 
						 const double rhobar_liq, 
						 const double rhobar_vap, 
						 const std::vector<double> & K, 
						 const int iS, 
						 const std::vector<double> & x, 
						 const std::vector<double> & y, 
						 const unsigned int N)
	{
		this->p.push_back(p);
		this->T.push_back(T);
		this->lnT.push_back(log(T));
		this->lnp.push_back(log(p));
		this->rhobar_liq.push_back(rhobar_liq);
		this->rhobar_vap.push_back(rhobar_vap);
		this->lnrhobar_liq.push_back(log(rhobar_liq));
		this->lnrhobar_vap.push_back(log(rhobar_vap));
		this->iS.push_back(iS);
		for (unsigned int i = 0; i < N; i++)
		{
			this->K[i].push_back(K[i]);
			this->x[i].push_back(x[i]);
			this->y[i].push_back(y[i]);
			this->lnK[i].push_back(log(K[i]));
		}
	};
};
class PhaseEnvelope
{
public:
	PhaseEnvelopeLog bubble, dew;
	double rhobar_liq, rhobar_vap;
	std::vector<double> K;
	Mixture *Mix;
	void build(double p0, const std::vector<double> &z, double beta_envelope);
	void to_python_files(std::string fname);
};

/*! 
This is the class that actually implements the mixture properties
*/
class Mixture
{

public:

	unsigned int N;
	Mixture(std::vector<Fluid *> pFluids);
	Mixture(std::string FluidsString);
	~Mixture();
	void initialize();

	double Rbar(const std::vector<double> &x);

	std::vector<Fluid *> pFluids;
	ReducingFunction * pReducing;
	ExcessTerm * pExcess;
	ResidualIdealMixture * pResidualIdealMix;

	PhaseEnvelope Envelope;

	SuccessiveSubstitutionVLE SS;
	NewtonRaphsonVLE NRVLE;

	/*! Returns the natural logarithm of K for component i using the method from Wilson as in
	\f[
	\ln K_i = \ln\left(\frac{p_{c,i}}{p}\right)+5.373(1+\omega_i)\left(1-\frac{T_{c,i}}{T}\right)
	\f]
	@param T Temperature [K]
	@param p Pressure [Pa]
	@param i Index of component [-]
	*/
	double Wilson_lnK_factor(double T, double p, int i);

	double phir(double tau, double delta, const std::vector<double> &x);
	double d2phir_dDelta_dTau(double tau, double delta, const std::vector<double> &x);
	double d2phir_dTau2(double tau, double delta, const std::vector<double> &x);
	double dphir_dDelta(double tau, double delta, const std::vector<double> &x);
	double d2phir_dDelta2(double tau, double delta, const std::vector<double> &x);
	double dphir_dTau(double tau, double delta, const std::vector<double> &x);
	double dphir_dxi(double tau, double delta, const std::vector<double> &x, int i);
	double d2phir_dxi_dTau(double tau, double delta, const std::vector<double> &x, int i);
	double d2phir_dxi_dDelta(double tau, double delta, const std::vector<double> &x, int i);
	double d2phirdxidxj(double tau, double delta, const std::vector<double> &x, int i, int j);

	/// Returns the fugacity for the given component for the given total reduced density and reciprocal reduced temperature
	double fugacity(double tau, double delta, const std::vector<double> &x, int i);
	
	/*! Density as a function of T,p,z
	@param T Temperature [K]
	@param p Pressure [Pa]
	@param z Bulk mole fractions [-]
	@param rhobar0 Guess value for molar density [mol/m^3]
	*/
	double rhobar_Tpz(double T, double p, const std::vector<double> &z, double rhobar0);

	/*! Temperature-pressure-bulk mole fraction flash calculation
	@param T Temperature [K]
	@param p Pressure [Pa]
	@param z Bulk mole fractions [-]
	@param rhobar Molar density [mol/m^3]
	@param x Liquid mole fractions [-] (if saturated)
	@param y Vapor mole fractions [-] (if saturated)
	*/
	void TpzFlash(double T, double p, const std::vector<double> &z, double &rhobar, std::vector<double> &x, std::vector<double> &y);

	/*! Objective function from Rachford-Rice (Not to be confused with the Gibbs function)
	\f[
	g(\beta) = \sum_i(y_i-x_i) = \sum_i z_i\frac{K_i-1}{1-\beta+\beta K_i}
	\f]
	@param z Bulk mole fractions [-]
	@param lnK Logarithm of the K factor [-]
	@param beta Molar fraction in the gas phase [-]
	*/
	double g_RachfordRice(const std::vector<double> &z, const std::vector<double> &lnK, double beta);

	/*! Objective function from Rachford-Rice (Not to be confused with the Gibbs function)
	\f[
	\frac{dg}{d\beta} = \sum_i z_i\frac{(K_i-1)^2}{(1-\beta+\beta K_i)^2};
	\f]
	@param z Bulk mole fractions [-]
	@param lnK Logarithm of the K factor [-]
	@param beta Molar fraction in the gas phase [-]
	*/
	double dgdbeta_RachfordRice(const std::vector<double> &z, const std::vector<double> &lnK, double beta);

	/*! The derivative term
	\f[
	n\left(\frac{\partial \phi^r}{\partial n_i} \right)_{T,V,n_j}
	\f]
	which is equal to
	\f{eqnarray*}{
	n\left(\frac{\partial \phi^r}{\partial n_i} \right)_{T,V,n_j} &=& \delta \phi^r_{\delta}\left[ 1-\frac{1}{\rho_r}\left[\left(\frac{\partial \rho_r}{\partial x_i}\right)_{x_j} - \sum_{k=1}^N x_k\left(\frac{\partial \rho_r}{\partial x_k}\right)_{x_j}  \right]\right]\\
	&& +\tau \phi^r_{\tau}\frac{1}{T_r}\left[\left(\frac{\partial T_r}{\partial x_i}\right)_{x_j} - \sum_{k=1}^N x_k\left(\frac{\partial T_r}{\partial x_k}\right)_{x_j}  \right]\\
	&& +\phi^r_{x_i}-\sum_{k=1}^{N}x_k\phi^r_{x_k}
	\f}
	See Table B4, Kunz, JCED, 2012 for the original term and the subsequent substitutions
	*/
	double ndphir_dni__constT_V_nj(double tau, double delta, const std::vector<double> &x, int i);

	/*! The partial molar volume
	\f[
	\hat v_i = \left( \frac{\partial V}{\partial n_i}\right)_{T,p,n_j} = \frac{-\left(\dfrac{\partial p}{\partial n_i}\right)_{T,V,n_j}}{\left(\dfrac{\partial p}{\partial V}\right)_{T,\bar n}}
	\f]
	from GERG monograph eqn 7.32
	*/
	double partial_molar_volume(double tau, double delta, const std::vector<double> &x, int i);

	/*! The derivative term
	\f[
	\left(\frac{\partial p}{\partial T} \right)_{V,\bar n} = \rho R(1+\delta \alpha_{\delta}^r-\delta \tau \alpha^r_{\delta\tau})
	\f]
	GERG 2004 Monograph equation 7.61
	*/
	double dpdT__constV_n(double tau, double delta, const std::vector<double> &x, int i);

	/*! The derivative term
	\f[
	n\left(\frac{\partial p}{\partial V} \right)_{T,\bar n} = -\rho^2 RT(1+2\delta \alpha_{\delta}^r+\delta^2\alpha^r_{\delta\delta})
	\f]
	GERG 2004 Monograph equation 7.62
	*/
	double ndpdV__constT_n(double tau, double delta, const std::vector<double> &x, int i);

	/*! The derivative term
	\f[
	n\left(\frac{\partial p}{\partial n_i} \right)_{T,V,n_j} = \rho RT\left[1+\delta\alpha_{\delta}^r\left[2- \frac{1}{\rho_r}\cdot n\left( \frac{\partial \rho_r}{\partial n_i}\right)_{n_j}\right] +\delta\cdot n\left(\frac{\partial\alpha_{\delta}^r}{\partial n_i}\right)_{T,V,n_j}\right]
	\f]
	GERG 2004 Monograph equation 7.63
	*/
	double ndpdni__constT_V_nj(double tau, double delta, const std::vector<double> &x, int i);

	/*!
	Natural logarithm of the fugacity coefficient
	*/
	double ln_fugacity_coefficient(double tau, double delta, const std::vector<double> &x, int i);

	/*!
	Derivative of the natural logarithm of the fugacity coefficient with respect to T
	*/
	double dln_fugacity_coefficient_dT__constrho(double tau, double delta, const std::vector<double> &x, int i);

	/*! The derivative term
	\f[
	\left(\frac{\partial \ln \phi_i}{\partial T} \right)_{p,\bar n} = \left(\frac{\partial^2n\alpha^r}{\partial T\partial n_i} \right)_{V,n_j} + \frac{1}{T}-\frac{\hat v}{RT}\left(\frac{\partial p}{\partial T}\right)_{V,\bar n}
	\f]
	GERG 2004 Monograph Eqn. 7.29
	*/
	double dln_fugacity_coefficient_dT__constp_n(double tau, double delta, const std::vector<double> &x, int i);

	// GERG Equation 7.42
	double dnphir_dni__constT_V_nj(double tau, double delta, const std::vector<double> &x, int i);

	/*! The derivative term
	\f[
	\left(\frac{\partial \ln \phi_i}{\partial p} \right)_{T,\bar n} = \frac{\hat v_i}{RT}-\frac{1}{p}
	\f]
	GERG 2004 Monograph Eqn. 7.30
	*/
	double dln_fugacity_coefficient_dp__constT_n(double tau, double delta, const std::vector<double> &x, int i);

	/*!
	\f[
	n\left(\frac{\partial \ln \phi_i}{\partial n_j}\right)_{T,p} = n\left(\frac{\partial^2n\alpha^r}{\partial n_j \partial n_i} \right)_{T,V}+1+\frac{n}{RT}\frac{\left(\frac{\partial p}{\partial n_j}\right)_{T,V,n_i}\left(\frac{\partial p}{\partial n_i}\right)_{T,V,n_j}}{\left(\frac{\partial p}{\partial V}\right)_{T,\bar n}}
	\f]
	which is also equal to 
	\f[
	n\left(\frac{\partial \ln \phi_i}{\partial n_j}\right)_{T,p} = n\left(\frac{\partial^2n\alpha^r}{\partial n_j \partial n_i} \right)_{T,V}+1-\frac{\hat v_i}{RT}\left[n\left(\frac{\partial p}{\partial n_j}\right)_{T,V,n_i}\right]
	\f]
	GERG 2004 Monograph Equation 7.31
	*/
	double ndln_fugacity_coefficient_dnj__constT_p(double tau, double delta, const std::vector<double> &x, int i, int j);


	/// Gernert Equation 3.115
	// Catch test provided
	double dln_fugacity_coefficient_dxj__constT_p_xi(double tau, double delta, const std::vector<double> &x, int i, int j);
	
	/// Gernert Equation 3.130
	/// Catch test provided
	double dpdxj__constT_V_xi(double tau, double delta, const std::vector<double> &x, int j);

	/// Gernert Equation 3.117
	double d2nphir_dni_dxj__constT_V(double tau, double delta, const std::vector<double> &x, int i, int j);

	/// Gernert Equation 3.119
	/// Catch test provided
	double dphir_dxj__constT_V_xi(double tau, double delta, const std::vector<double> &x, int j);

	/// Gernert Equation 3.118
	/// Catch test provided
	double d_ndphirdni_dxj__constT_V_xi(double tau, double delta, const std::vector<double> &x, int i, int j);

	/// Gernert Equation 3.134
	/// Catch test provided
	double d_dphirddelta_dxj__constT_V_xi(double tau, double delta, const std::vector<double> &x, int j);

	/// Gernert Equation 3.121
	/// Catch test provided
	double ddelta_dxj__constT_V_xi(double tau, double delta, const std::vector<double> &x, int j);

	/// Gernert Equation 3.122
	/// Catch test provided
	double dtau_dxj__constT_V_xi(double tau, double delta, const std::vector<double> &x, int j);

	/*!
	\f[
	\left(\frac{\partial^2n\alpha^r}{\partial T\partial n_i} \right)_{V,n_j} = \left( \frac{\partial}{\partial T}\left(\frac{\partial n \alpha^r}{\partial n_i}\right)_{T,V,n_j} \right)_{V,\bar n}
	\f]
	\f[
	\left(\frac{\partial^2n\alpha^r}{\partial T\partial n_i} \right)_{V,n_j} = -\frac{\tau}{T}\left[\alpha_{\tau}^r +\left( \frac{\partial}{\partial \tau}\left(n\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\bar x}\right]
	\f]
	GERG 2004 Monograph, equations 7.44 and 7.51
	*/
	double d2nphir_dni_dT(double tau, double delta, const std::vector<double> &x, int i);

	/*! The derivative term
	\f{eqnarray*}{
	\frac{\partial }{\partial \tau} \left( n\left(\frac{\partial \phi^r}{\partial n_i} \right)_{T,V,n_j} \right) &=& \delta \phi^r_{\delta\tau}\left[ 1-\frac{1}{\rho_r}\left[\left(\frac{\partial \rho_r}{\partial x_i}\right)_{x_j} - \sum_{k=1}^N x_k\left(\frac{\partial \rho_r}{\partial x_k}\right)_{x_j}  \right]\right]\\
	&& +(\tau \phi^r_{\tau\tau}+\phi^r_{\tau})\frac{1}{T_r}\left[\left(\frac{\partial T_r}{\partial x_i}\right)_{x_j} - \sum_{k=1}^N x_k\left(\frac{\partial T_r}{\partial x_k}\right)_{x_j}  \right]\\
	&& +\phi^r_{x_i\tau}-\sum_{k=1}^{N}x_k\phi^r_{x_k\tau}
	\f}
	GERG 2004 Monograph Equation 7.51 and Table B4, Kunz, JCED, 2012
	*/
	double d_ndphirdni_dTau(double tau, double delta, const std::vector<double> &x, int i);

	/*! The derivative term
	\f{eqnarray*}{
	\left(\frac{\partial }{\partial \delta} \left( n\left(\frac{\partial \phi^r}{\partial n_i} \right)_{T,V,n_j} \right)\right)_{\tau,\bar x} &=& (\alpha_{\delta}^r+\delta\alpha_{\delta\delta}^r)\left[1-\frac{1}{\rho_r}\cdot n\left(\frac{\partial \rho_r}{\partial n_i}\right)_{n_j} \right] \\
	&+&\tau\alpha^r_{\delta\tau}\frac{1}{T_r}\cdot n\left(\frac{\partial T_r}{\partial n_i}\right)_{n_j}\\
	&+&\phi^r_{\delta x_i}-\sum_{k=1}^{N}x_k\phi^r_{\delta x_k}
	\f}
	GERG 2004 Monograph Equation 7.50 and Table B4, Kunz, JCED, 2012
	*/
	double d_ndphirdni_dDelta(double tau, double delta, const std::vector<double> &x, int i);

	/*!GERG 2004 Monograph equation 7.41:
	\f[
	n\left(\frac{\partial^2n\alpha^r}{\partial n_i \partial n_j} \right)_{T,V} = n\left( \frac{\partial}{\partial n_j}\left(\frac{\partial n\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)_{T,V,n_i}
	\f]
	and
	GERG 2004 Monograph equation 7.46:
	\f[
	n\left( \frac{\partial}{\partial n_j}\left(\frac{\partial n\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)_{T,V,n_i} = n\left( \frac{\partial \alpha^r}{\partial n_j}\right)_{T,V,n_i} + n\left( \frac{\partial}{\partial n_j}\left(n\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j} \right) \right)_{T,V,n_i}
	\f]
	GERG 2004 Monograph equation 7.47:
	\f{eqnarray*}{
	n\left( \frac{\partial}{\partial n_j}\left(n\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j} \right) \right)_{T,V,n_i} &=& \left( \frac{\partial}{\partial \delta}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\tau,\bar x}\cdot n\left(\frac{\partial\delta}{\partial n_j}\right)_{T,V,n_i}\\ 
	&+& \left( \frac{\partial}{\partial \tau}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\tau,\bar x}\cdot n\left(\frac{\partial\tau}{\partial n_j}\right)_{T,V,n_i}\\
	&+& \left( \frac{\partial}{\partial x_j}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\tau,x_i}-\sum_{k=1}^{N}x_k \left( \frac{\partial}{\partial x_k}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\tau,x_i}\\
	\f}
	*/
	double nd2nphirdnidnj__constT_V(double tau, double delta, const std::vector<double> &x, int i, int j);

	/*! 
	\f[
	n\left(\frac{\partial \delta}{\partial n_i} \right)_{T,V,n_j} = \delta - \frac{\delta}{\rho_r}\cdot n\left(\frac{\partial \rho_r}{\partial n_i} \right)_{n_j}
	\f]
	GERG 2004 Monograph equation 7.48
	*/
	double nddeltadni__constT_V_nj(double tau, double delta, const std::vector<double> &x, int i);
	/*! 
	\f[
	n\left(\frac{\partial \tau}{\partial n_i} \right)_{T,V,n_j} = \frac{\tau}{T_r}\cdot n\left(\frac{\partial T_r}{\partial n_i} \right)_{n_j}
	\f]
	GERG 2004 Monograph equation 7.49
	*/
	double ndtaudni__constT_V_nj(double tau, double delta, const std::vector<double> &x, int i);

	/*!GERG 2004 Monograph equation 7.52:
	\f{eqnarray*}{
	\left( \frac{\partial}{\partial x_j}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\tau,x_i} &=& \delta\alpha_{\delta x_j}^{r}\left[ 1-\frac{1}{\rho_r}\cdot n\left(\frac{\partial \rho_r}{\partial n_i}\right)_{n_j}\right] \\
	&-& \delta\alpha_{\delta}^{r}\frac{1}{\rho_r}\left[ \left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial \rho_r}{\partial n_i}\right)_{n_j}\right)\right)_{x_i}-\frac{1}{\rho_r}\left(\frac{\partial \rho_r}{\partial x_j}\right)_{x_i}\cdot n\left(\frac{\partial \rho_r}{\partial n_i}\right)_{n_j}\right] \\
	&+& \tau\alpha_{\tau x_j}^r\frac{1}{T_r}\cdot n\left(\frac{\partial T_r}{\partial n_i}\right)_{n_j}\\
	&+& \tau\alpha_{\tau}^{r}\frac{1}{T_r}\left[ \left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial T_r}{\partial n_i}\right)_{n_j}\right)\right)_{x_i}-\frac{1}{T_r}\left(\frac{\partial T_r}{\partial x_j}\right)_{x_i}\cdot n\left(\frac{\partial T_r}{\partial n_i}\right)_{n_j}\right] \\
	&+& \alpha_{x_ix_j}^r-\alpha_{x_j}^r-\sum_{m=1}^Nx_m\alpha_{x_jx_m}^r
	\f}
	*/
	double d_ndphirdni_dxj__constdelta_tau_xi(double tau, double delta, const std::vector<double> &x, int i, int j);

	double saturation_p_preconditioner(double p, const std::vector<double> &z);

	double saturation_p_Wilson(double beta, double p, const std::vector<double> &z, double T_guess, std::vector<double> &K);

	double saturation_p(double beta, double p, const std::vector<double> &z, std::vector<double> &x, std::vector<double> &y);

	/*! Calculate the mixture molar density based on the use of the Peng-Robinson equation of state
	*/
	double rhobar_pengrobinson(double T, double p, const std::vector<double> &x, int solution);
	
	void x_and_y_from_K(double beta, const std::vector<double> &K, const std::vector<double> &z, std::vector<double> &x, std::vector<double> &y);


	/*! Load the excess parameters (departure function parameters)
	@param i 0-based index of first component
	@param j 0-based index of second component
	*/
	void load_excess_values(int i, int j);

	/*! Load the reducing parameters
	@param i 0-based index of first component
	@param j 0-based index of second component
	*/
	std::map<std::string, double> load_reducing_values(int i, int j);

	void check_MethaneEthane();
	void check_WaterEthanol();

	void test();
};



#endif


