#ifndef MIXTURES_H
#define MIXTURES_H

#include "Helmholtz.h"
#include "FluidClass.h"

enum mix_sat_types {TYPE_BUBBLEPOINT, TYPE_DEWPOINT};

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
	virtual double Tr(std::vector<double> *x) = 0;
	/// The derivative of reduced temperature with respect to component i mole fraction
	virtual double dTr_dxi(std::vector<double> *x, int i) = 0;
	/// The molar reducing density
	virtual double rhorbar(std::vector<double> *x) = 0;
	///Derivative of the molar reducing density with respect to component i mole fraction
	virtual double drhorbardxi__constxj(std::vector<double> *x, int i) = 0;

	/// Set the coefficients based on reducing parameters loaded from JSON
	virtual void set_coeffs_from_map(int i, int j, std::map<std::string,double >) = 0;

	virtual double d2rhorbardxi2__constxj(std::vector<double> *x, int i) = 0;
	virtual double d2rhorbardxidxj(std::vector<double> *x, int i, int j) = 0;
	virtual double d2Trdxi2__constxj(std::vector<double> *x, int i) = 0;
	virtual double d2Trdxidxj(std::vector<double> *x, int i, int j) = 0;

	double dndTr_dni_dxj__constxi(std::vector<double> *x, int i);
	double dndrhorbar_dni_dxj__constxi(std::vector<double> *x, int i);

	double ndrhorbar_dni__constnj(std::vector<double> *x, int i);
	double ndTr_dni__constnj(std::vector<double> *x, int i);
};

/*! 
The Reducing parameter model used by the GERG-2008 formulation to yield the
reducing parameters \f$ \bar\rho_r \f$ and \f$ T_r \f$ and derivatives thereof
*/
class GERG2008ReducingFunction : public ReducingFunction
{
protected:
	unsigned int N;
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
	double Tr(std::vector<double> *x);
	/// The derivative of reduced temperature with respect to component i mole fraction
	double dTr_dxi(std::vector<double> *x, int i);
	/// The molar reducing density
	double rhorbar(std::vector<double> *x);
	///Derivative of the molar reducing density with respect to component i mole fraction
	double drhorbardxi__constxj(std::vector<double> *x, int i);
	double dvrbardxi__constxj(std::vector<double> *x, int i);

	double d2vrbardxi2__constxj(std::vector<double> *x, int i);
	double d2rhorbardxi2__constxj(std::vector<double> *x, int i);
	double d2vrbardxidxj(std::vector<double> *x, int i, int j);
	double d2rhorbardxidxj(std::vector<double> *x, int i, int j);
	double d2Trdxi2__constxj(std::vector<double> *x, int i);
	double d2Trdxidxj(std::vector<double> *x, int i, int j);

	/// Set the coefficients based on reducing parameters loaded from JSON
	void set_coeffs_from_map(int i, int j, std::map<std::string,double >);

	double c_Y_ij(int i, int j, std::vector< std::vector< double> > * beta, std::vector< std::vector< double> > *gamma, std::vector< std::vector< double> > *Y_c);
	double f_Y_ij(std::vector<double> *x, int i, int j, std::vector< std::vector< double> > * beta);

	//double dYrdxi__constxj(std::vector<double> *x, int i){throw 1;};
	//double d2Yrdxi2__constxj(std::vector<double> *x, int i){throw 1;};
	//double d2Yrdxidxj(std::vector<double> *x, int i){throw 1;};

	double dfYkidxi__constxk(std::vector<double> *x, int k, int i,std::vector< std::vector< double> > * beta);
	double dfYikdxi__constxk(std::vector<double> *x, int i, int k, std::vector< std::vector< double> > * beta);
	double d2fYkidxi2__constxk(std::vector<double> *x, int k, int i, std::vector< std::vector< double> > * beta);
	double d2fYikdxi2__constxk(std::vector<double> *x, int i, int k, std::vector< std::vector< double> > * beta);
	double d2fYijdxidxj(std::vector<double> *x, int i, int k, std::vector< std::vector< double> > * beta);
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
	virtual double phir(double tau, double delta, std::vector<double> *x) = 0;
	virtual double dphir_dDelta(double tau, double delta, std::vector<double> *x) = 0;
	virtual double d2phir_dDelta2(double tau, double delta, std::vector<double> *x) = 0;
	virtual double d2phir_dDelta_dTau(double tau, double delta, std::vector<double> *x) = 0;
	virtual double dphir_dTau(double tau, double delta, std::vector<double> *x) = 0;
	virtual double d2phir_dTau2(double tau, double delta, std::vector<double> *x) = 0;
	virtual void set_coeffs_from_map(std::map<std::string,std::vector<double> >) = 0;
};

class GERG2008DepartureFunction : public DepartureFunction
{
protected:
	phir_power phi1;
	phir_GERG2008_gaussian phi2;
public:
	GERG2008DepartureFunction(){};
	~GERG2008DepartureFunction(){};
	double phir(double tau, double delta, std::vector<double> *x);
	double dphir_dDelta(double tau, double delta, std::vector<double> *x);
	double d2phir_dDelta_dTau(double tau, double delta, std::vector<double> *x);
	double dphir_dTau(double tau, double delta, std::vector<double> *x);
	double d2phir_dDelta2(double tau, double delta, std::vector<double> *x);
	double d2phir_dTau2(double tau, double delta, std::vector<double> *x);
	void set_coeffs_from_map(std::map<std::string,std::vector<double> >);
};

class ExcessTerm
{
public:
	int N;
	std::vector<std::vector<DepartureFunction*> > DepartureFunctionMatrix;
	std::vector<std::vector<double> > F;
	ExcessTerm(int N);
	~ExcessTerm();
	double phir(double tau, double delta, std::vector<double> *x);
	double dphir_dDelta(double tau, double delta, std::vector<double> *x);
	double d2phir_dDelta2(double tau, double delta, std::vector<double> *x);
	double d2phir_dDelta_dTau(double tau, double delta, std::vector<double> *x);
	double dphir_dTau(double tau, double delta, std::vector<double> *x);
	double d2phir_dTau2(double tau, double delta, std::vector<double> *x);
	double dphir_dxi(double tau, double delta, std::vector<double> *x, unsigned int i);
	double d2phir_dxi_dTau(double tau, double delta, std::vector<double> *x, unsigned int i);
	double d2phir_dxi_dDelta(double tau, double delta, std::vector<double> *x, unsigned int i);
	void set_coeffs_from_map(int i, int j, std::map<std::string,std::vector<double> >);
};


class ResidualIdealMixture
{
protected:
	std::vector<Fluid*> pFluids;
public:
	ResidualIdealMixture(std::vector<Fluid*> pFluids);
	double phir(double tau, double delta, std::vector<double> *x);
	double dphir_dDelta(double tau, double delta, std::vector<double> *x);
	double d2phir_dDelta2(double tau, double delta, std::vector<double> *x);
	double d2phir_dDelta_dTau(double tau, double delta, std::vector<double> *x);
	double d2phir_dTau2(double tau, double delta, std::vector<double> *x);
	double dphir_dTau(double tau, double delta, std::vector<double> *x);
};

class Mixture;  // Forward declaration since some classes that are members of Mixture take pointers to Mixture

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
A class to do successive substitution given guess values.  This class will then be included in the Mixture class

A class is used rather than a function so that it is easier to store iteration histories, other output values, etc.
*/
class SuccessiveSubstitution
{
public:
	bool logging;
	int Nsteps;
	int Nmax;
	Mixture *Mix;
	std::vector<double> K, ln_phi_liq, ln_phi_vap;
	std::vector<SuccessiveSubstitutionStep> step_logger;

	SuccessiveSubstitution(){};
	double call(int type, double T, double p, std::vector<double> *z, std::vector<double> *x, std::vector<double> *y);
};


/*! 
This is the class that actually implements the mixture properties
*/
class Mixture
{
	
public:

	Mixture(std::vector<Fluid *> pFluids);
	~Mixture();

	double Rbar(std::vector<double> *x);

	std::vector<Fluid *> pFluids;
	ReducingFunction * pReducing;
	ExcessTerm * pExcess;
	ResidualIdealMixture * pResidualIdealMix;

	SuccessiveSubstitution SS;

	/*! Returns the natural logarithm of K for component i using the method from Wilson as in
	\f[
	\ln K_i = \ln\left(\frac{p_{c,i}}{p}\right)+5.373(1+\omega_i)\left(1-\frac{T_{c,i}}{T}\right)
	\f]
	@param T Temperature [K]
	@param p Pressure [kPa]
	@param i Index of component [-]
	*/
	double Wilson_lnK_factor(double T, double p, int i);
	double phir(double tau, double delta, std::vector<double> *x);
	double d2phir_dDelta_dTau(double tau, double delta, std::vector<double> *x);
	double d2phir_dTau2(double tau, double delta, std::vector<double> *x);
	double dphir_dDelta(double tau, double delta, std::vector<double> *x);
	double d2phir_dDelta2(double tau, double delta, std::vector<double> *x);
	double dphir_dTau(double tau, double delta, std::vector<double> *x);
	double dphir_dxi(double tau, double delta, std::vector<double> *x, int i);
	double d2phir_dxi_dTau(double tau, double delta, std::vector<double> *x, int i);
	double d2phir_dxi_dDelta(double tau, double delta, std::vector<double> *x, int i);
	/// Returns the fugacity for the given component for the given total reduced density and reciprocal reduced temperature
	double fugacity(double tau, double delta, std::vector<double> *x, int i);
	
	/*! Density as a function of T,p,z
	@param T Temperature [K]
	@param p Pressure [Pa]
	@param z Bulk mole fractions [-]
	@param rhobar0 Guess value for molar density [mol/m^3]
	*/
	double rhobar_Tpz(double T, double p, std::vector<double> *z, double rhobar0);

	/*! Temperature-pressure-bulk mole fraction flash calculation
	@param T Temperature [K]
	@param p Pressure [kPa]
	@param z Bulk mole fractions [-]
	@param rhobar Molar density [mol/m^3]
	@param x Liquid mole fractions [-] (if saturated)
	@param y Vapor mole fractions [-] (if saturated)
	*/
	void TpzFlash(double T, double p, std::vector<double> *z, double *rhobar, std::vector<double> *x, std::vector<double> *y);

	/*! Objective function from Rachford-Rice (Not to be confused with the Gibbs function)
	\f[
	g(\beta) = \sum_i(y_i-x_i) = \sum_i z_i\frac{K_i-1}{1-\beta+\beta K_i}
	\f]
	@param z Bulk mole fractions [-]
	@param lnK Logarithm of the K factor [-]
	@param beta Molar fraction in the gas phase [-]
	*/
	double g_RachfordRice(std::vector<double> *z, std::vector<double> *lnK, double beta);

	/*! Objective function from Rachford-Rice (Not to be confused with the Gibbs function)
	\f[
	\frac{dg}{d\beta} = \sum_i z_i\frac{(K_i-1)^2}{(1-\beta+\beta K_i)^2};
	\f]
	@param z Bulk mole fractions [-]
	@param lnK Logarithm of the K factor [-]
	@param beta Molar fraction in the gas phase [-]
	*/
	double dgdbeta_RachfordRice(std::vector<double> *z, std::vector<double> *lnK, double beta);

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
	double ndphir_dni(double tau, double delta, std::vector<double> *x, int i);

	/*! The partial molar volume
	\f[
	\hat v_i = \left( \frac{\partial V}{\partial n_i}\right)_{T,p,n_j} = \frac{-\left(\dfrac{\partial p}{\partial n_i}\right)_{T,V,n_j}}{\left(\dfrac{\partial p}{\partial V}\right)_{T,\bar n}}
	\f]
	from GERG monograph eqn 7.32
	*/
	double partial_molar_volume(double tau, double delta, std::vector<double> *x, int i);

	/*! The derivative term
	\f[
	\left(\frac{\partial p}{\partial T} \right)_{V,\bar n} = \rho R(1+\delta \alpha_{\delta}^r-\delta \tau \alpha^r_{\delta\tau})
	\f]
	GERG 2004 Monograph equation 7.61
	*/
	double dpdT__constV_n(double tau, double delta, std::vector<double> *x, int i);

	/*! The derivative term
	\f[
	n\left(\frac{\partial p}{\partial V} \right)_{T,\bar n} = -\rho^2 RT(1+2\delta \alpha_{\delta}^r+\delta^2\alpha^r_{\delta\delta})
	\f]
	GERG 2004 Monograph equation 7.62
	*/
	double ndpdV__constT_n(double tau, double delta, std::vector<double> *x, int i);

	/*! The derivative term
	\f[
	n\left(\frac{\partial p}{\partial n_i} \right)_{T,V,n_j} = \rho RT\left[1+\delta\alpha_{\delta}^r\left[2- \frac{1}{\rho_r}\cdot n\left( \frac{\partial \rho_r}{\partial n_i}\right)_{n_j}\right] +\delta\cdot n\left(\frac{\partial\alpha_{\delta}^r}{\partial n_i}\right)_{T,V,n_j}\right]
	\f]
	GERG 2004 Monograph equation 7.63
	*/
	double ndpdni__constT_V_nj(double tau, double delta, std::vector<double> *x, int i);

	/*!
	Natural logarithm of the fugacity coefficient
	*/
	double ln_fugacity_coefficient(double tau, double delta, std::vector<double> *x, int i);

	/*!
	Derivative of the natural logarithm of the fugacity coefficient with respect to T
	*/
	double dln_fugacity_coefficient_dT__constrho(double tau, double delta, std::vector<double> *x, int i);

	/*! The derivative term
	\f[
	\left(\frac{\partial \ln \phi_i}{\partial T} \right)_{p,\bar n} = \left(\frac{\partial^2n\alpha^r}{\partial T\partial n_i} \right)_{V,n_j} + \frac{1}{T}-\frac{\hat v}{RT}\left(\frac{\partial p}{\partial T}\right)_{V,\bar n}
	\f]
	GERG 2004 Monograph Eqn. 7.29
	*/
	double dln_fugacity_coefficient_dT__constp_n(double tau, double delta, std::vector<double> *x, int i);

	/*! The derivative term
	\f[
	\left(\frac{\partial \ln \phi_i}{\partial p} \right)_{T,\bar n} = \frac{\hat v_i}{RT}-\frac{1}{p}
	\f]
	GERG 2004 Monograph Eqn. 7.30
	*/
	double dln_fugacity_coefficient_dp__constT_n(double tau, double delta, std::vector<double> *x, int i);

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
	double ndln_fugacity_coefficient_dnj__constT_p(double tau, double delta, std::vector<double> *x, int i);

	/*!
	\f[
	\left(\frac{\partial^2n\alpha^r}{\partial T\partial n_i} \right)_{V,n_j} = \left( \frac{\partial}{\partial T}\left(\frac{\partial n \alpha^r}{\partial n_i}\right)_{T,V,n_j} \right)_{V,\bar n}
	\f]
	\f[
	\left(\frac{\partial^2n\alpha^r}{\partial T\partial n_i} \right)_{V,n_j} = -\frac{\tau}{T}\left[\alpha_{\tau}^r +\left( \frac{\partial}{\partial \tau}\left(n\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\bar x}\right]
	\f]
	GERG 2004 Monograph, equations 7.44 and 7.51
	*/
	double d2nphir_dni_dT(double tau, double delta, std::vector<double> *x, int i); // Implemented

	

	/*! The derivative term
	\f{eqnarray*}{
	\frac{\partial }{\partial \tau} \left( n\left(\frac{\partial \phi^r}{\partial n_i} \right)_{T,V,n_j} \right) &=& \delta \phi^r_{\delta\tau}\left[ 1-\frac{1}{\rho_r}\left[\left(\frac{\partial \rho_r}{\partial x_i}\right)_{x_j} - \sum_{k=1}^N x_k\left(\frac{\partial \rho_r}{\partial x_k}\right)_{x_j}  \right]\right]\\
	&& +(\tau \phi^r_{\tau\tau}+\phi^r_{\tau})\frac{1}{T_r}\left[\left(\frac{\partial T_r}{\partial x_i}\right)_{x_j} - \sum_{k=1}^N x_k\left(\frac{\partial T_r}{\partial x_k}\right)_{x_j}  \right]\\
	&& +\phi^r_{x_i\tau}-\sum_{k=1}^{N}x_k\phi^r_{x_k\tau}
	\f}
	GERG 2004 Monograph Equation 7.51 and Table B4, Kunz, JCED, 2012
	*/
	double dndphir_dni_dTau(double tau, double delta, std::vector<double> *x, int i);

	/*! The derivative term
	\f{eqnarray*}{
	\left(\frac{\partial }{\partial \delta} \left( n\left(\frac{\partial \phi^r}{\partial n_i} \right)_{T,V,n_j} \right)\right)_{\tau,\bar x} &=& (\alpha_{\delta}^r+\delta\alpha_{\delta\delta}^r)\left[1-\frac{1}{\rho_r}\cdot n\left(\frac{\partial \rho_r}{\partial n_i}\right)_{n_j} \right] \\
	&+&\tau\alpha^r_{\delta\tau}\frac{1}{T_r}\cdot n\left(\frac{\partial T_r}{\partial n_i}\right)_{n_j}\\
	&+&\phi^r_{\delta x_i}-\sum_{k=1}^{N}x_k\phi^r_{\delta x_k}
	\f}
	GERG 2004 Monograph Equation 7.50 and Table B4, Kunz, JCED, 2012
	*/
	double dndphir_dni_dDelta(double tau, double delta, std::vector<double> *x, int i);

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
	\f[
	n\left( \frac{\partial}{\partial n_j}\left(n\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j} \right) \right)_{T,V,n_i} = ....
	\f]
	*/
	double nd2nphirdnidnj__constT_V(double tau, double delta, std::vector<double> *x, int i){throw 1;};

	double nddeltadnj__constT_V_ni(double tau, double delta, std::vector<double> *x, int ){throw 1;};

	double ndtaudnj__constT_V_n(double tau, double delta, std::vector<double> *x, int ){throw 1;};

	/*!GERG 2004 Monograph equation 7.52:
	*/
	double dndphirdni_dxj__constdelta_tau_xi(double tau, double delta, std::vector<double> *x, int i){throw 1;};

	double saturation_p(int type, double p, std::vector<double> *z, std::vector<double> *x, std::vector<double> *y);

	/*! Calculate the mixture molar density based on the use of the Peng-Robinson equation of state
	*/
	double rhobar_pengrobinson(double T, double p, std::vector<double> *x, int solution);

	/*! Load the excess parameters (departure function parameters)
	@param i 0-based index of first component
	@param j 0-based index of second component
	*/
	std::map<std::string,std::vector<double> > load_excess_values(int i, int j);

	/*! Load the reducing parameters
	@param i 0-based index of first component
	@param j 0-based index of second component
	*/
	std::map<std::string, double> load_reducing_values(int i, int j);
};


#endif
