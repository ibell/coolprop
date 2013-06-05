#ifndef MIXTURES_H
#define MIXTURES_H

#include "Helmholtz.h"
#include "FluidClass.h"

typedef std::vector<std::vector<double> > STLMatrix;

/*! 
An abstract base class for the reducing function to allow for
Lemmon-Jacobsen, GERG, or other reducing function to yield the
reducing parameters \f$ \bar\rho_r \f$ and \f$ T_r \f$
*/
class ReducingFunction
{
public:
	ReducingFunction(){};
	virtual ~ReducingFunction(){};
	/// The reduced temperature
	virtual double Tr(std::vector<double> x) = 0;
	/// The derivative of reduced temperature with respect to component i mole fraction
	virtual double dTr_dxi(std::vector<double> x, int i) = 0;
	/// The molar reducing density
	virtual double rhorbar(std::vector<double> x) = 0;
	///Derivative of the molar reducing density with respect to component i mole fraction
	virtual double drhorbar_dxi(std::vector<double> x, int i) = 0;
};

/*! 
The Reducing parameter model used by the GERG-2008 formulation to yield the
reducing parameters \f$ \bar\rho_r \f$ and \f$ T_r \f$ and derivatives thereof
*/
class GERGReducingFunction : public ReducingFunction
{
protected:
	STLMatrix beta_v; //!< \f$ \beta_{v,ij} \f$ from GERG-2008
	STLMatrix gamma_v; //!< \f$ \gamma_{v,ij} \f$ from GERG-2008
	STLMatrix beta_T; //!< \f$ \beta_{T,ij} \f$ from GERG-2008
	STLMatrix gamma_T; //!< \f$ \gamma_{T,ij} \f$ from GERG-2008
	std::vector<Fluid *> pFluids; //!< List of pointers to fluids
public:
	GERGReducingFunction(std::vector<Fluid *> pFluids, STLMatrix beta_v, STLMatrix gamma_v, STLMatrix beta_T, STLMatrix gamma_T)
	{
		this->pFluids = pFluids;
		this->beta_v = beta_v;
		this->gamma_v = gamma_v;
		this->beta_T = beta_T;
		this->gamma_T = gamma_T;
	};
	/// Default destructor
	~GERGReducingFunction(){};
	/// The reduced temperature
	double Tr(std::vector<double> x);
	/// The derivative of reduced temperature with respect to component i mole fraction
	double dTr_dxi(std::vector<double> x, int i);
	/// The molar reducing density
	double rhorbar(std::vector<double> x);
	///Derivative of the molar reducing density with respect to component i mole fraction
	double drhorbar_dxi(std::vector<double> x, int i);
};








/*! 
The abstract base clas for departure functions for the excess part of the Helmholtz energy
*/
class DepartureFunction
{
protected:
	STLMatrix F; //!< The \f$ F_{ij} \f$ matrix
public:
	DepartureFunction(){};
	/// Instantiator for the ABC for the DepartureFunction
	DepartureFunction(std::vector<Fluid *> pFluids);
	virtual ~DepartureFunction(){};
	
	/// The excess Helmholtz energy of the binary pair
	/// Pure-virtual function (must be implemented in derived class
	virtual double phir(double tau, double delta, std::vector<double> x) = 0;
	virtual double dphir_dDelta(double tau, double delta, std::vector<double> x) = 0;
	virtual double dphir_dTau(double tau, double delta, std::vector<double> x) = 0;
	virtual double dphir_dxi(double tau, double delta, std::vector<double> x, int i) = 0;
};

class GERGDepartureFunction : public DepartureFunction
{
protected:
	phir_power phi1;
	phir_GERG_gaussian phi2;
public:
	GERGDepartureFunction(STLMatrix F);
	~GERGDepartureFunction(){};
	double phir(double tau, double delta, std::vector<double> x);
	double dphir_dDelta(double tau, double delta, std::vector<double> x);
	double dphir_dTau(double tau, double delta, std::vector<double> x);
	double dphir_dxi(double tau, double delta, std::vector<double> x, int i);
};

class ResidualIdealMixture
{
protected:
	std::vector<Fluid*> pFluids;
public:
	ResidualIdealMixture(std::vector<Fluid*> pFluids);
	double phir(double tau, double delta, std::vector<double> x);
	double dphir_dDelta(double tau, double delta, std::vector<double> x);
	double dphir_dTau(double tau, double delta, std::vector<double> x);
};







/*! 
This is the class that actually implements the mixture properties
*/

class Mixture
{
protected:
	std::vector<Fluid *> pFluids;
	ReducingFunction * pReducing;
	DepartureFunction * pExcess;
	ResidualIdealMixture * pResidualIdealMix;
public:
	Mixture(std::vector<Fluid *> pFluids);
	~Mixture();

	/*! Returns the natural logarithm of K for component i as in
	\f[
	\ln K = \ln\left(\frac{p_{c,i}}{p}\right)+5.373(1+\omega_i)\left(1-\frac{T_{c,i}}{T}\right)
	\f]
	using the method from Wilson
	@param T Temperature [K]
	@param p Pressure [kPa]
	@param i Index of component [-]
	*/
	double Wilson_lnK_factor(double T, double p, int i);
	double phir(double tau, double delta, std::vector<double> x);
	double dphir_dDelta(double tau, double delta, std::vector<double> x);
	double dphir_dTau(double tau, double delta, std::vector<double> x);
	double dphir_dxi(double tau, double delta, std::vector<double> x, int i);
	/// Returns the fugacity for the given component for the given total reduced density and reciprocal reduced temperature
	double fugacity(double tau, double delta, std::vector<double> x, int i);
	
	/*!
	Temperature-pressure-bulk mole fraction flash calculation

	@param T Temperature [K]
	@param p Pressure [kPa]
	@param z Bulk mole fractions [-]
	*/
	double TpzFlash(double T, double p, std::vector<double> z);

	/*!
	Objective function from Rachford-Rice

	Not to be confused with the Gibbs function

	\f[
	g(\beta) = \sum_i(y_i-x_i) = \sum_i z_i\frac{K_i-1}{1-\beta+\beta K_i}
	\f]
	
	@param z Bulk mole fractions [-]
	@param lnK Logarithm of the K factor [-]
	@param beta Molar fraction in the gas phase [-]
	*/
	double g_RachfordRice(std::vector<double> z, std::vector<double> lnK, double beta);
};


#endif
