
#ifndef INCOMPRESSIBLE_SOLUTION_H
#define INCOMPRESSIBLE_SOLUTION_H

#include <string>
#include <vector>
#include <math.h>
#include "CPExceptions.h"
#include "CoolPropTools.h"
#include "IncompBase.h"
#include <stdio.h>

/**
Notes for developers:

If you want to add a fluid, add its definition to the header 
and then add it to the map in the constructor in SolutionsContainer
in IncompSolution.cpp
**/

/// Base class for simplified brine/solution models
/** Employs the base functions implemented in IncompBase.h and
 *  provides properties as function of temperature, pressure
 *  and composition. */
class IncompressibleSolution : public IncompressibleClass{

/** A few redefinitions are needed. The base class does not support
 *  composition as input. We also have to have an additional check
 *  for maximum and minimum concentration as well as the freezing point.
 */
protected:
	double xmin, xmax;
	std::vector< std::vector<double> > cRho;
	std::vector< std::vector<double> > cHeat;
	std::vector< std::vector<double> > cVisc;
	std::vector< std::vector<double> > cCond;
	std::vector< std::vector<double> > cPsat;
	             std::vector<double>   cTfreeze;

public:
	// Constructor
	IncompressibleSolution(){
		xmin  = -1.;
		xmax  = -1.;
	};

	std::vector<std::vector<double> > getcRho() const {return cRho;}
	std::vector<std::vector<double> > getcHeat() const {return cHeat;}
	std::vector<std::vector<double> > getcVisc() const {return cVisc;}
	std::vector<std::vector<double> > getcCond() const {return cCond;}
	std::vector<std::vector<double> > getcPsat() const {return cPsat;}
	            std::vector<double>   getcTfreeze() const {return cTfreeze;}

	/* All functions need T, p and x as input. Might not
	 * be necessary, but gives a clearer structure.
	 */
	/// Density as a function of temperature, pressure and composition.
    virtual double rho (double T_K, double p, double x){return -_HUGE;};
    /// Heat capacities as a function of temperature, pressure and composition.
    virtual double c   (double T_K, double p, double x){return -_HUGE;};
	virtual double cp  (double T_K, double p, double x){return c(T_K,p,x);};
	virtual double cv  (double T_K, double p, double x){return c(T_K,p,x);};
    /// Entropy as a function of temperature, pressure and composition.
    virtual double s   (double T_K, double p, double x){return -_HUGE;};
    /// Internal energy as a function of temperature, pressure and composition.
    virtual double u   (double T_K, double p, double x){return -_HUGE;};
    /// Enthalpy as a function of temperature, pressure and composition.
    virtual double h   (double T_K, double p, double x){return -_HUGE;};
    /// Viscosity as a function of temperature, pressure and composition.
    virtual double visc(double T_K, double p, double x){return -_HUGE;};
    /// Thermal conductivity as a function of temperature, pressure and composition.
    virtual double cond(double T_K, double p, double x){return -_HUGE;};
    /// Saturation pressure as a function of temperature and composition.
    virtual double psat(double T_K          , double x){return -_HUGE;};
    /// Freezing temperature as a function of pressure and composition.
    virtual double Tfreeze(         double p, double x){return -_HUGE;};

    void testInputs(double T_K, double p, double x);


protected:
    /* Define internal energy and enthalpy as functions of the
	 * other properties to provide data in case there are no
	 * coefficients.
	 */

	/// Enthalpy from x, u, p and rho.
	/** Calculate enthalpy as a function of temperature and
	 *  pressure employing functions for internal energy and
	 *  density. Provides consistent formulations. */
	double h_u(double T_K, double p, double x) {
		return u(T_K,p,x)+p/rho(T_K,p,x);
	};

	/// Internal energy from x, h, p and rho.
	/** Calculate internal energy as a function of temperature
	 *  and pressure employing functions for enthalpy and
	 *  density. Provides consistent formulations. */
	double u_h(double T_K, double p, double x) {
		return h(T_K,p,x)-p/rho(T_K,p,x);
	};


	/*
	 * Some more functions to provide a single implementation
	 * of important routines.
	 * We start with the check functions that can validate input
	 * in terms of pressure p, temperature T and composition x.
	 */

	/// Check validity of temperature input.
	/** Compares the given temperature T to the result of a
	 *  freezing point calculation. This is not necessarily
	 *  defined for all fluids, default values do not cause errors. */
	bool checkT(double T_K, double p, double x);

	/// Check validity of pressure input.
	/** Compares the given pressure p to the saturation pressure at
	 *  temperature T and throws and exception if p is lower than
	 *  the saturation conditions.
	 *  The default value for psat is -1 yielding true if psat
	 *  is not redefined in the subclass.
	 *  */
	bool checkP(double T_K, double p, double x);

	/// Check validity of composition input.
	/** Compares the given composition x to a stored minimum and
	 *  maximum value. Enforces the redefinition of xmin and
	 *  xmax since the default values cause an error. */
	bool checkX(double x);

	/// Check validity of temperature, pressure and composition input.
	bool checkTPX(double T, double p, double x);
};


/** Basic functions to access the list of incompressible fluids.
 *  Used here for convenience, but does not really contribute
 *  any functionality.
 */
bool IsIncompressibleSolution(std::string name);
double IncompSolutionSI(long iOutput, double T, double p, double x, long iFluid);
double IncompSolutionSI(long iOutput, double T, double p, double x, std::string name);
double IncompSolutionSI(long iOutput, double T, double p, std::string name);  // TODO Solutions: Remove as soon as possible
std::string getSolutionName(std::string name); // TODO Solutions: Remove as soon as possible
double getSolutionConc(std::string name); // TODO Solutions: Remove as soon as possible


/** Handle all the objects in a single list of incompressible
 *  solutions and brines. */
class SolutionsContainer {
private:
  std::vector<IncompressibleSolution*> solution_list;
  std::map<std::string,IncompressibleSolution*> solution_map;

public:
  IncompressibleSolution * get_solution(long index);
  IncompressibleSolution * get_solution(std::string name);
  void set_solutions(std::vector<IncompressibleSolution*> list);

  SolutionsContainer();
  ~SolutionsContainer();
};

/// Class to access Lithium-Bromide solutions
/** Employs some basic wrapper-like functionality
 *  to bridge the gap between the solution functions
 *  used in the paper by Pátek and Klomfar:
 *  http://dx.doi.org/10.1016/j.ijrefrig.2005.10.007
 *
 *  We owe gratitude to the authors for providing
 *  both access to the paper as well as the equations
 *  in the form of C source code. */
class LiBrSolution : public IncompressibleSolution{

protected:
	static double const M_H2O; /* kg/mol, molar mass of H2O */
	static double const M_LiBr; /* kg/mol, molar mass of LiBr */
	static double const T0; /* K, constant */

	/* Critical point of H2O */
	static double const Tc_H2O; /* K, temperature  */
	static double const pc_H2O; /* MPa, pressure */
	static double const rhoc_H2O; /* mol/m^3 (322 kg/m^3), molar density */
	static double const hc_H2O; /* J/mol, molar enthalpy */
	static double const sc_H2O; /* J/(mol.K) molar entropy*/

	/*Triple point of H2O */
	static double const Tt_H2O; /* K, temperature */
	static double const cpt_H2O; /* J/(mol.K), molar isobaric heat capacity*/

	double ps_mix(double T, double x)
	/* Equation (1) */
	{
		static double m[8] = { 3.0, 4.0, 4.0, 8.0, 1.0, 1.0, 4.0, 6.0 };
		static double n[8] = { 0.0, 5.0, 6.0, 3.0, 0.0, 2.0, 6.0, 0.0 };
		static double t[8] = { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };
		static double a[8] = { -2.41303e2, 1.91750e7, -1.75521e8, 3.25430e7,
				3.92571e2, -2.12626e3, 1.85127e8, 1.91216e3 };
		double tau, suma = 0.0;
		int i;

		tau = T / Tc_H2O;
		for (i = 0; i <= 7; i++)
			suma += a[i] * pow(x, m[i]) * pow(0.4 - x, n[i]) * pow(tau, t[i]);
		return (ps_H2O(T - suma));

	} /* end function ps_mix */

	double rho_mix(double T, double x)
	/* Equation (2) */
	{
		static double m[2] = { 1.0, 1.0 };
		static double n[2] = { 0.0, 6.0 };
		static double a[2] = { 1.746, 4.709 };

		double tau, suma = 0.0;
		int i;

		tau = T / Tc_H2O;
		for (i = 0; i <= 1; i++)
			suma += a[i] * pow(x, m[i]) * pow(tau, n[i]);

		return ((1.0 - x) * rho_H2O(T) + rhoc_H2O * suma);

	} /* end function rho_mix */

	double cp_mix(double T, double x)
	/* Equation (3) */
	{
		static double m[8] = { 2.0, 3.0, 3.0, 3.0, 3.0, 2.0, 1.0, 1.0 };
		static double n[8] = { 0.0, 0.0, 1.0, 2.0, 3.0, 0.0, 3.0, 2.0 };
		static double t[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 3.0, 4.0 };
		static double a[8] = { -1.42094e1, 4.04943e1, 1.11135e2, 2.29980e2,
				1.34526e3, -1.41010e-2, 1.24977e-2, -6.83209e-4 };

		double tau, suma = 0.0;
		int i;

		tau = Tc_H2O / (T - T0);
		for (i = 0; i <= 7; i++)
			suma += a[i] * pow(x, m[i]) * pow(0.4 - x, n[i]) * pow(tau, t[i]);

		return ((1.0 - x) * cp_H2O(T) + cpt_H2O * suma);

	} /* end function cp_mix */

	double h_mix(double T, double x)
	/* Equation (4) */
	{
		static double m[30] = { 1.0, 1.0, 2.0, 3.0, 6.0, 1.0, 3.0, 5.0, 4.0,
				5.0, 5.0, 6.0, 6.0, 1.0, 2.0, 2.0, 2.0, 5.0, 6.0, 7.0, 1.0, 1.0,
				2.0, 2.0, 2.0, 3.0, 1.0, 1.0, 1.0, 1.0 };

		static double n[30] = { 0.0, 1.0, 6.0, 6.0, 2.0, 0.0, 0.0, 4.0, 0.0,
				4.0, 5.0, 5.0, 6.0, 0.0, 3.0, 5.0, 7.0, 0.0, 3.0, 1.0, 0.0, 4.0,
				2.0, 6.0, 7.0, 0.0, 0.0, 1.0, 2.0, 3.0 };

		static double t[30] = { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0,
				2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0,
				4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 5.0 };

		static double a[30] = { 2.27431, -7.99511, 3.85239e2, -1.63940e4,
				-4.22562e2, 1.13314e-1, -8.33474, -1.73833e4, 6.49763,
				3.24552e3, -1.34643e4, 3.99322e4, -2.58877e5, -1.93046e-3,
				2.80616, -4.04479e1, 1.45342e2, -2.74873, -4.49743e2,
				-1.21794e1, -5.83739e-3, 2.33910e-1, 3.41888e-1, 8.85259,
				-1.78731e1, 7.35179e-2, -1.79430e-4, 1.84261e-3, -6.24282e-3,
				6.84765e-3 };

		double tau, suma = 0.0;
		int i;

		tau = Tc_H2O / (T - T0);
		for (i = 0; i <= 29; i++)
			suma += a[i] * pow(x, m[i]) * pow(0.4 - x, n[i]) * pow(tau, t[i]);

		return ((1.0 - x) * h_H2O(T) + hc_H2O * suma);

	} /* end function h_mix */

	double s_mix(double T, double x)
	/* Equation (5) */
	{
		static double m[29] = { 1.0, 1.0, 2.0, 3.0, 6.0, 1.0, 3.0, 5.0, 1.0,
				2.0, 2.0, 4.0, 5.0, 5.0, 6.0, 6.0, 1.0, 3.0, 5.0, 7.0, 1.0, 1.0,
				1.0, 2.0, 3.0, 1.0, 1.0, 1.0, 1.0 };

		static double n[29] = { 0.0, 1.0, 6.0, 6.0, 2.0, 0.0, 0.0, 4.0, 0.0,
				0.0, 4.0, 0.0, 4.0, 5.0, 2.0, 5.0, 0.0, 4.0, 0.0, 1.0, 0.0, 2.0,
				4.0, 7.0, 1.0, 0.0, 1.0, 2.0, 3.0 };

		static double t[29] = { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0,
				2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0,
				4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 5.0 };

		static double a[29] = { 1.53091, -4.52564, 6.98302e+2, -2.16664e+4,
				-1.47533e+3, 8.47012e-2, -6.59523, -2.95331e+4, 9.56314e-3,
				-1.88679e-1, 9.31752, 5.78104, 1.38931e+4, -1.71762e+4,
				4.15108e+2, -5.55647e+4, -4.23409e-3, 3.05242e+1, -1.67620,
				1.48283e+1, 3.03055e-3, -4.01810e-2, 1.49252e-1, 2.59240,
				-1.77421e-1, -6.99650e-5, 6.05007e-4, -1.65228e-3, 1.22966e-3 };

		double tau, suma = 0.0;
		int i;

		tau = Tc_H2O / (T - T0);
		for (i = 0; i <= 28; i++)
			suma += a[i] * pow(x, m[i]) * pow(0.4 - x, n[i]) * pow(tau, t[i]);

		return ((1.0 - x) * s_H2O(T) + sc_H2O * suma);

	} /* end function s_mix */

	double ps_H2O(double T)
	/* Equation (28) */
	{
		static double a[7] = { 0.0, -7.85951783, 1.84408259, -11.7866497,
				22.6807411, -15.9618719, 1.80122502 };

		double tau, ps;

		tau = 1 - T / Tc_H2O;

		ps = pc_H2O
				* exp(
						Tc_H2O / T
								* (a[1] * tau + a[2] * pow(tau, 1.5)
										+ a[3] * pow(tau, 3.0)
										+ a[4] * pow(tau, 3.5)
										+ a[5] * pow(tau, 4.0)
										+ a[6] * pow(tau, 7.5)));

		return (ps * 1.0e6);

	} /* end function ps_H2O */

	double rho_H2O(double T)
	/* Equation (29) */
	{
		static double b[7] = { 0.0, 1.99274064, 1.09965342, -0.510839303,
				-1.75493479, -45.5170352, -6.7469445e5 };
		double theta, rho;

		theta = 1.0 - T / Tc_H2O;

		rho = rhoc_H2O
				* (1.0 + b[1] * pow(theta, 1.0 / 3.0)
						+ b[2] * pow(theta, 2.0 / 3.0)
						+ b[3] * pow(theta, 5.0 / 3.0)
						+ b[4] * pow(theta, 16.0 / 3.0)
						+ b[5] * pow(theta, 43.0 / 3.0)
						+ b[6] * pow(theta, 110.0 / 3.0));

		return (rho);
	} /* end function rho_H2O */

	double cp_H2O(double T)
	/* Equation (30) */
	{
		static double a[5] =
				{ 1.38801, -2.95318, 3.18721, -0.645473, 9.18946e5 };
		static double b[5] = { 0.0, 2.0, 3.0, 6.0, 34.0 };
		static double c[5] = { 0.0, 2.0, 3.0, 5.0, 0.0 };

		double suma = 0.0;
		int i;

		for (i = 0; i <= 4; i++)
			suma += a[i] * exp(b[i] * log(1.0 - T / Tc_H2O))
					* exp(c[i] * log(T / Tt_H2O));

		return (cpt_H2O * suma);

	} /* end function cp_H2O */

	double h_H2O(double T)
	/* Equation (31) */
	{
		static double a[4] = { -4.37196e-1, 3.03440e-1, -1.29582, -1.76410e-1 };
		static double alpha[4] = { 1.0 / 3.0, 2.0 / 3.0, 5.0 / 6.0, 21.0 / 6.0 };

		double suma = 0.0;
		int i;

		for (i = 0; i <= 3; i++)
			suma += a[i] * exp(alpha[i] * log(1.0 - T / Tc_H2O));

		return (hc_H2O * (1.0 + suma));

	} /* end function h_H2O */

	double s_H2O(double T)
	/* Equation (32)  */
	{
		static double a[4] = { -3.34112e-1, -8.47987e-1, -9.11980e-1, -1.64046 };
		static double alpha[4] = { 1.0 / 3.0, 3.0 / 3.0, 8.0 / 3.0, 24.0 / 3.0 };

		double suma = 0.0;
		int i;

		for (i = 0; i <= 3; i++)
			suma += a[i] * exp(alpha[i] * log(1.0 - T / Tc_H2O));

		return (sc_H2O * (1.0 + suma));

	} /* end function s_H2O */


	/** Finished with the code from the paper. Now we need to
	 *  convert the molar values to mass-based units. */
	double massToMole(double w)
	/* Equation (7)  */
	{
		return (w/M_LiBr)/(w/M_LiBr+(1.-w)/M_H2O);
		//return (w*M_LiBr)/(w*M_LiBr+(1.-w)*M_H2O);
	}

	double molarToSpecific(double w, double value)
	/* Equation (7,8)  */
	{
		double x = massToMole(w);
		//return w/(x*M_LiBr) * value;
		return 1. / ( x*M_LiBr + (1.-x)*M_H2O ) * value;
	}

	static const bool debug = false;

public:

	LiBrSolution(){
        name = std::string("LiBr");
        description = std::string("Lithium-Bromide solution from Patek2006");
        reference = std::string("Patek2006");

		Tmin     = 273.00;
		Tmax     = 500.00;
		TminPsat = Tmin;

		xmin     = 0.0;
		xmax     = 1.0;

	};

	double rho(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		return 1./molarToSpecific(x, 1./rho_mix(T_K,massToMole(x)));
	}
	double c(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		return molarToSpecific(x, cp_mix(T_K,massToMole(x)));
	}
	double h(double T_K, double p, double x){
		return h_u(T_K,p,x);
	}
	double s(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		return molarToSpecific(x, s_mix(T_K,massToMole(x)));
	}
	double visc(double T_K, double p, double x){
		throw ValueError("Viscosity is not defined for LiBr-solutions.");
	}
	double cond(double T_K, double p, double x){
		throw ValueError("Thermal conductivity is not defined for LiBr-solutions.");
	}
	double u(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		return molarToSpecific(x, h_mix(T_K,massToMole(x)));
	}
	double psat(double T_K, double x){
		//checkT(T_K,p,x);
		if (debug) throw ValueError(format("Your concentration is %f in kg/kg and %f in mol/mol.",x,massToMole(x)));
		return ps_mix(T_K,massToMole(x));
	};
	double Tfreeze(double p, double x){
		if (debug) throw ValueError(format("No freezing point data available for Lithium-Bromide: p=%f, x=%f",p,x));
		return Tmin;
	}
};


/// Class to use Melinder and SecCool parameters
/** Employs some basic wrapper-like functionality
 *  to bridge the gap between the solution functions
 *  used in CoolProp and the definition used in
 *  Melinder's book and Ian's original implementation and
 *  the definition used in SecCool. Please visit:
 *  http://en.ipu.dk/Indhold/refrigeration-and-energy-technology/seccool.aspx
 *  Many thanks to Morten Juel Skovrup for providing
 *  this nice piece of software as well as the parameters
 *  needed to calculate the composition based properties. */
class BaseSolution : public IncompressibleSolution{

protected:
	double Tbase, xbase;
	/// Some more general purpose functions
	double baseFunction(std::vector<double> const& coefficients, double T_K, double p, double x);
	std::vector< std::vector<double> > makeMatrix(std::vector<double> const& coefficients);

public:

	BaseSolution(){
		Tbase = -1.;
		xbase = -1.;
	};

	double getTbase() const {return Tbase;}
	double getxbase() const {return xbase;}

	double getTInput(double curTValue){
		return curTValue-Tbase;
	}

	double getxInput(double curxValue){
		return (curxValue-xbase)*100.0;
	}

	double rho(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cRho,6,4);
		return polyval(cRho, getxInput(x), getTInput(T_K));
	}
	double c(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cHeat,6,4);
		return polyval(cHeat, getxInput(x), getTInput(T_K));
	}
	double h(double T_K, double p, double x){
		return h_u(T_K,p,x);
	}
	double s(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cHeat,6,4);
		return polyfracintcentral(cHeat, getxInput(x), T_K, Tref, Tbase);
	}
	double visc(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cVisc,6,4);
		return expval(cVisc, getxInput(x), getTInput(T_K), 2)/1e3;
	}
	double cond(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cCond,6,4);
		return polyval(cCond, getxInput(x), getTInput(T_K));
	}
	double u(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cHeat,6,4);
		return polyint(cHeat, getxInput(x), getTInput(T_K), getTInput(Tref));
	}
//	double psat(double T_K, double x){
//		//checkT(T_K,p,x);
//		Incompressible::checkCoefficients(cPsat,6,4);
//		if (T_K<TminPsat || TminPsat<0){
//			return -1.;
//		} else {
//			return expval(cPsat, getxInput(x), getTInput(T_K), 2);
//		}
//	};
//	double Tfreeze(double p, double x){
//		IncompressibleClass::checkCoefficients(cTfreeze,5);
//		std::vector<double> tmpVector(cTfreeze.begin()+1,cTfreeze.end());
//		return polyval(tmpVector, x*100.0-cTfreeze[0])+273.15;
//	}
};


class MelinderSolution : public BaseSolution{

protected:
	/// Convert pre-v4.0-style coefficient array to new format
	std::vector<std::vector<double> > convertCoeffs(double *oldestCoeffs, int A, int B);

public:
	MelinderSolution(){

        name = std::string("MelinderSolution");
        description = std::string("Test Methanol Melinder");
        reference = std::string("Melinder-BOOK-2010");

		Tmin     = -50 + 273.15;
		Tmax     =  40 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.6;

		Tbase    =  3.5359 + 273.15;
	    xbase    = 30.5128 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
		{ -26.29           	, 958.1				,3887			,   0.4175				,   1.153			},//
		{  -0.000002575    	,  -0.4151			,   7.201		,   0.0007271			,  -0.03866			},
		{  -0.000006732    	,  -0.002261		,  -0.08979		,   0.0000002823		,   0.0002779		},
		{   0.000000163		,   0.0000002998	,  -0.000439	,   0.000000009718		,  -0.000001543		},
		{  -1.187			,  -1.391			, -18.5			,  -0.004421			,   0.005448		},//
		{  -0.00001609		,  -0.0151			,   0.2984		,  -0.00002952			,   0.0001008		},
		{   0.000000342		,   0.0001113		,  -0.001865	,   0.00000007336		,  -0.000002809		},
		{   0.0000000005687	,  -0.0000003264	,  -0.00001718	,   0.0000000004328		,   0.000000009811	},
		{  -0.01218			,  -0.01105			,  -0.03769		,   0.00002044			,  -0.0005552		},//
		{   0.0000003865	,   0.0001828		,  -0.01196		,   0.0000003413		,   0.000008384		},
		{   0.000000008768	,  -0.000001641		,   0.00009801	,  -0.000000003665		,  -0.00000003997	},
		{  -0.0000000002095	,   0.0000000151	,   0.000000666	,  -0.00000000002791	,  -0.0000000003466	},
		{  -0.00006823		,  -0.0001208		,  -0.003776	,   0.0000002943		,   0.000003038		},//
		{   0.00000002137	,   0.000002992		,  -0.00005611	,  -0.0000000009646		,  -0.00000007435	},
		{  -0.0000000004271	,   0.000000001455	,  -0.0000007811,   0.00000000003174	,   0.0000000007442	},
		{   0.0000001297	,   0.000004927		,  -0.0001504	,  -0.0000000008666		,   0.00000006669	},//
		{  -0.0000000005407	,  -0.0000001325	,   0.000007373	,  -0.0000000000004573	,  -0.0000000009105	},
		{   0.00000002363	,  -0.00000007727	,   0.000006433	,  -0.0000000002033		,  -0.0000000008472	}};//

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = get_col(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };

	/// Define freezing point calculations
	double Tfreeze(double p, double x){
		IncompressibleClass::checkCoefficients(cTfreeze,6);
		return polyval(cTfreeze, getxInput(x))+273.15;
	}
};

class EGSolution : public MelinderSolution{
public:
	EGSolution(){

        name = std::string("MEG");
        description = std::string("Ethylene Glycol");
        reference = std::string("Melinder-BOOK-2010");

		Tmin     =-100 + 273.15;
		Tmax     = 100 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.6;

		Tbase    =  31.728 + 273.15;
	    xbase    = 30.8462 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
			{-15.25,1034,3.737,0.472,0.4705},
			{-0.000001566,-0.4781,0.00293,0.0008903,-0.0255},
			{-0.0000002278,-0.002692,-0.000004675,-0.000001058,0.0001782},
			{0.000000002169,0.000004725,-0.00000001389,-0.000000002789,-0.0000007669},
			{-0.808,1.311,-0.01799,-0.004286,0.02471},
			{-0.000001339,-0.006876,0.0001046,-0.00001473,-0.0001171},
			{0.00000002047,0.00004805,-0.0000004147,0.0000001059,0.000001052},
			{-0.00000000002717,0.0000000169,0.0000000001847,-0.0000000001142,-0.00000001634},
			{-0.01334,0.0000749,-0.00009933,0.00001747,0.000003328},
			{0.00000006322,0.00007855,0.0000003516,0.00000006814,0.000001086},
			{0.0000000002373,-0.0000003995,0.000000005109,-0.000000003612,0.00000001051},
			{-0.000000000002183,0.000000004982,-0.00000000007138,0.000000000002365,-0.0000000006475},
			{-0.00007293,-0.0001062,0.00000261,0.00000003017,0.000001659},
			{0.000000001764,0.000001229,-0.000000001189,-0.000000002412,0.000000003157},
			{-0.00000000002442,-0.00000001153,-0.0000000001643,0.00000000004004,0.0000000004063},
			{0.000001006,-0.0000009623,0.00000001537,-0.000000001322,0.00000003089},
			{-0.00000000007662,-0.00000007211,-0.0000000004272,0.00000000002555,0.0000000001831},
			{0.00000000114,0.00000004891,-0.000000001618,0.00000000002678,-0.000000001865}
		};

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = get_col(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);
		// fix the typo for Ethylene Glycol
		for (unsigned int i=0; i<cHeat.size(); i++){
			for (unsigned int j=0; j<cHeat[i].size(); j++){
				cHeat[i][j] *= 1000;
			}
		}

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };
};

class PGSolution : public MelinderSolution{
public:
	PGSolution(){

        name = std::string("MPG");
        description = std::string("Propylene Glycol");
        reference = std::string("Melinder-BOOK-2010");

		Tmin     =-100 + 273.15;
		Tmax     = 100 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.6;

		Tbase    = 32.7083 + 273.15;
	    xbase    = 30.7031 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
			{-13.25,1018,3882,0.4513,0.6837},
			{-0.0000382,-0.5406,2.699,0.0007955,-0.03045},
			{0.0000007865,-0.002666,-0.001659,0.00000003482,0.0002525},
			{-0.000000001733,0.00001347,-0.00001032,-0.000000005966,-0.000001399},
			{-0.6631,0.7604,-13.04,-0.004795,0.03328},
			{0.000006774,-0.00945,0.0507,-0.00001678,-0.0003984},
			{-0.00000006242,0.00005541,-0.00004752,0.00000008941,0.000004332},
			{-0.0000000007819,-0.0000001343,0.000001522,0.0000000001493,-0.0000000186},
			{-0.01094,-0.002498,-0.1598,0.00002076,0.00005453},
			{0.00000005332,0.000027,0.00009534,0.0000001563,-0.000000086},
			{-0.000000004169,-0.0000004018,0.00001167,-0.000000004615,-0.00000001593},
			{0.00000000003288,0.000000003376,-0.0000000487,0.000000000009897,-0.00000000004465},
			{-0.0002283,-0.000155,0.0003539,-0.00000009083,-0.0000039},
			{-0.00000001131,0.000002829,0.00003102,-0.000000002518,0.0000001054},
			{0.0000000001918,-0.000000007175,-0.000000295,0.00000000006543,-0.000000001589},
			{-0.000003409,-0.000001131,0.00005,-0.0000000005952,-0.00000001587},
			{0.00000000008035,-0.00000002221,-0.0000007135,-0.00000000003605,0.0000000004475},
			{0.00000001465,0.00000002342,-0.0000004959,0.00000000002104,0.000000003564}
		};

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = get_col(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };
};

class EASolution : public MelinderSolution{
public:
	EASolution(){

        name = std::string("MEA");
        description = std::string("Ethyl Alcohol (Ethanol)");
        reference = std::string("Melinder-BOOK-2010");

		Tmin     =-100 + 273.15;
		Tmax     =  40 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.6;

		Tbase    =  8.1578 + 273.15;
	    xbase    = 29.2361 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
			{-19.41,961.9,4204,0.4067,1.474},
			{-0.0003668,-0.5222,2.319,0.0006775,-0.04745},
			{-0.00004005,-0.003281,-0.03042,0.0000003105,0.0004314},
			{0.000001524,0.00001569,0.000686,-0.00000002,-0.000003023},
			{-0.954,-1.433,-21.02,-0.005008,0.01565},
			{-0.00001209,-0.01989,0.4927,-0.00002377,-0.00004106},
			{0.000002877,0.000187,-0.003072,-0.00000003216,-0.000005135},
			{-0.00000004394,-0.0000009154,-0.0000569,0.00000000008362,0.00000007004},
			{-0.002648,-0.0226,-0.3714,0.00002801,-0.0008435},
			{-0.0000003173,0.0002281,-0.002335,0.0000002669,0.0000164},
			{0.000000008652,-0.00000008581,-0.0000196,-0.000000003606,-0.0000001091},
			{-0.0000000003717,0.000000004056,0.0000007461,0.00000000001552,-0.000000001967},
			{0.0003851,-0.000169,0.01743,-0.00000002009,0.000007552},
			{0.0000000134,0.000008594,-0.0002969,-0.000000006813,-0.0000001118},
			{-0.000000002091,-0.00000009607,0.000001901,0.0000000001429,0.000000001899},
			{-0.0000002858,0.00001291,-0.00006292,-0.000000001506,0.0000001529},
			{0.0000000009312,-0.000000159,0.000005353,0.0000000001167,-0.0000000009481},
			{-0.000000167,-0.00000008318,-0.00000829,-0.00000000001653,-0.00000000413}
		};

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = get_col(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };
};

class MASolution : public MelinderSolution{
public:
	MASolution(){

        name = std::string("MMA");
        description = std::string("Methyl Alcohol (Methanol)");
        reference = std::string("Melinder-BOOK-2010");

		Tmin     =-100 + 273.15;
		Tmax     =  40 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.6;

		Tbase    =  3.5359 + 273.15;
	    xbase    = 30.5128 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
			{-26.29,958.1,3887,0.4175,1.153},
			{-0.000002575,-0.4151,7.201,0.0007271,-0.03866},
			{-0.000006732,-0.002261,-0.08979,0.0000002823,0.0002779},
			{0.000000163,0.0000002998,-0.000439,0.000000009718,-0.000001543},
			{-1.187,-1.391,-18.5,-0.004421,0.005448},
			{-0.00001609,-0.0151,0.2984,-0.00002952,0.0001008},
			{0.000000342,0.0001113,-0.001865,0.00000007336,-0.000002809},
			{0.0000000005687,-0.0000003264,-0.00001718,0.0000000004328,0.000000009811},
			{-0.01218,-0.01105,-0.03769,0.00002044,-0.0005552},
			{0.0000003865,0.0001828,-0.01196,0.0000003413,0.000008384},
			{0.000000008768,-0.000001641,0.00009801,-0.000000003665,-0.00000003997},
			{-0.0000000002095,0.0000000151,0.000000666,-0.00000000002791,-0.0000000003466},
			{-0.00006823,-0.0001208,-0.003776,0.0000002943,0.000003038},
			{0.00000002137,0.000002992,-0.00005611,-0.0000000009646,-0.00000007435},
			{-0.0000000004271,0.000000001455,-0.0000007811,0.00000000003174,0.0000000007442},
			{0.0000001297,0.000004927,-0.0001504,-0.0000000008666,0.00000006669},
			{-0.0000000005407,-0.0000001325,0.000007373,-0.0000000000004573,-0.0000000009105},
			{0.00000002363,-0.00000007727,0.000006433,-0.0000000002033,-0.0000000008472}

		};

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = get_col(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };
};

class GLSolution : public MelinderSolution{
public:
	GLSolution(){

        name = std::string("MGL");
        description = std::string("Glycerol");
        reference = std::string("Melinder-BOOK-2010");

		Tmin     =-100 + 273.15;
		Tmax     =  40 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.6;

		Tbase    =  8.9110 + 273.15;
	    xbase    = 36.1905 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
			{-13,1093,3486,0.4532,1.52},
			{-0.0008638,-0.3624,3.766,0.0009782,-0.03729},
			{0.000006895,-0.002451,-0.0001222,-0.0000001196,0.0003572},
			{0.0000005229,0.00001547,0.00003219,-0.000000008778,-0.000003648},
			{-0.5742,2.74,-22.5,-0.003223,0.0445},
			{-0.000007991,-0.008373,0.1811,-0.00001951,-0.0002688},
			{-0.0000007515,0.00009596,0.0003766,0.0000001178,0.0000008876},
			{0.0000000171,-0.0000006999,-0.0000173,0.0000000001048,-0.00000002209},
			{-0.009119,0.004081,-0.03258,0.000005539,0.0003633},
			{0.000002973,0.00001808,0.0007249,-0.00000006878,-0.000004088},
			{0.000000002379,-0.0000008516,-0.00002133,-0.000000002587,0.00000006219},
			{-0.000000001237,0.00000001441,0.0000004907,0.00000000002262,0.0000000006331},
			{-0.0001641,-0.00004744,0.002922,-0.0000002073,0.000001069},
			{-0.000000005313,0.000001833,-0.00005346,-0.0000000002235,-0.0000001248},
			{0.0000000004546,-0.00000003077,0.00000023,-0.00000000001421,0.000000003019},
			{-0.000002408,-0.000001415,0.00002238,0.000000003,0.00000005729},
			{-0.000000001682,0.00000001863,-0.0000005383,0.0000000001604,-0.00000000178},
			{-0.000000007734,-0.000000006097,-0.0000005944,0.0000000002221,0.000000002116}
		};

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = get_col(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };
};

class AMSolution : public MelinderSolution{
public:
	AMSolution(){

        name = std::string("MAM");
        description = std::string("Ammonia (NH3)");
        reference = std::string("Melinder-BOOK-2010");

		Tmin     =-100 + 273.15;
		Tmax     =  30 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.3;

		Tbase    = -4.6490 + 273.15;
	    xbase    = 16.0784 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
			{-25.76,944.5,4233,0.4551,0.9255},
			{-0.0001817,-0.2743,-1.618,0.001673,-0.03439},
			{0.00001204,-0.003113,0.0161,-0.000002214,0.0003217},
			{0.0000005567,0.000003349,0.00001662,0.0000001228,-0.000004544},
			{-2.385,-2.914,1.145,-0.005216,0.01327},
			{0.00002315,-0.02307,0.02715,0.000003544,0.0001856},
			{0.0000001415,0.0001341,0.001072,0.000001057,-0.00001646},
			{-0.00000004244,-0.0000005151,-0.00005266,-0.00000003474,0.0000003004},
			{-0.07636,0.02262,-0.001965,0.00008909,-0.0005979},
			{-0.00000229,0.00006645,0.003472,0.000003526,-0.00002184},
			{-0.000000262,-0.0000087,-0.00009051,-0.0000001782,0.000001297},
			{0.0000000001786,0.00000008999,0.000002106,0.000000001858,-0.00000001141},
			{-0.00261,-0.0006685,0.002131,0.000006807,-0.0001097},
			{-0.000000376,-0.000001002,0.00004117,-0.0000003394,0.00000167},
			{0.0000000136,0.0000003309,0.0000004446,0.000000008315,-0.00000003377},
			{-0.000073,0.000001635,0.0002136,-0.0000001898,0.000003065},
			{0.00000003524,0.0000004465,-0.00001354,0.000000006304,-0.00000006166},
			{-0.000001075,0.0000006298,-0.000008551,-0.00000001361,0.0000003244}
		};

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = get_col(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };
};

class KCSolution : public MelinderSolution{
public:
	KCSolution(){

        name = std::string("MKC");
        description = std::string("Potassium Carbonate (K2CO3)");
        reference = std::string("Melinder-BOOK-2010");

		Tmin     =-100 + 273.15;
		Tmax     =  40 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.4;

		Tbase    = 11.2422 + 273.15;
	    xbase    = 22.0833 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
			{-10.3,1216,3217,0.5622,0.8063},
			{-0.0001575,-0.4114,1.492,0.001656,-0.02362},
			{0.000003598,-0.003383,-0.001539,0.000002108,0.0001851},
			{0.000000004324,0.0000299,-0.00003546,-0.00000004482,-0.000002372},
			{-0.7786,10.64,-37.33,-0.000892,0.03624},
			{0.00001112,-0.007614,-0.01775,0.000002031,-0.00001262},
			{-0.0000003479,0.00005214,0.0005416,-0.0000002616,-0.0000003022},
			{-0.000000001244,-0.0000008087,0.00000659,0.0000000007965,-0.000000007761},
			{-0.02766,0.04413,0.2023,-0.0000005233,0.0006659},
			{0.000001616,0.00007806,0.004553,-0.0000002962,-0.00001611},
			{-0.00000001681,-0.000001173,0.00003587,-0.000000009174,0.000000153},
			{0.00000000003847,0.00000005658,-0.0000003707,0.0000000001027,0.000000001061},
			{-0.0008226,-0.0001333,0.01971,-0.0000009283,0.00001077},
			{-0.00000004913,-0.000002381,0.0001367,-0.00000001814,-0.00000009796},
			{0.000000001395,0.0000001696,-0.000003667,0.0000000008767,0.00000000307},
			{-0.000002372,-0.00001246,0.0003405,-0.00000001011,-0.0000001097},
			{-0.000000002886,0.0000002708,-0.00001676,0.0000000008471,0.000000007825},
			{0.0000003251,0.0000002035,-0.00003488,0.000000001311,-0.000000008453}
		};

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = get_col(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };
};

class CASolution : public MelinderSolution{
public:
	CASolution(){

        name = std::string("MCA");
        description = std::string("Calcium Chloride (CaCl2)");
        reference = std::string("Melinder-BOOK-2010");

		Tmin     =-100.0 + 273.15;
		Tmax     =  40.0 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.3;

		Tbase    = 7.52570 + 273.15;
	    xbase    = 18.7414 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
			{-16.21,1171,3133,0.558,0.8939},
			{-0.0001344,-0.1463,2.81,0.00146,-0.02647},
			{0.000005073,-0.001302,-0.01563,0.0000003861,0.0001718},
			{-0.0000000482,-0.0001871,-0.00001233,0.00000001307,-0.0000007918},
			{-1.555,9.847,-44.8,-0.00115,0.04389},
			{0.00002146,-0.02488,0.03271,-0.00001008,0.0002102},
			{-0.0000015,-0.000553,-0.001205,-0.0000000654,-0.0000008688},
			{0.00000002219,0.00001665,0.000009346,0.0000000004728,-0.00000004353},
			{-0.05496,0.03389,0.9511,-0.00001784,0.0009121},
			{0.0000009415,-0.002302,-0.005191,0.0000003496,0.000003993},
			{0.00000006185,0.00003526,0.0002282,-0.00000000484,0.0000003785},
			{-0.000000001723,0.0000002788,-0.000000929,-0.0000000002011,-0.000000009979},
			{-0.002624,0.001062,0.01472,-0.0000004415,0.000001963},
			{-0.0000001082,0.00006291,0.0001615,-0.000000003462,-0.0000004892},
			{0.000000003036,0.000001806,0.000005073,0.0000000003922,0.00000000001526},
			{-0.0001726,0.00002785,-0.001346,0.00000004699,0.0000002997},
			{-0.000000004396,0.000006423,-0.000009754,0.0000000006179,-0.00000003631},
			{-0.000004494,-0.000001811,-0.00006674,0.000000002929,0.00000003435}
		};

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = get_col(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };
};

class MGSolution : public MelinderSolution{
public:
	MGSolution(){

        name = std::string("MMG");
        description = std::string("(MgCl2)");
        reference = std::string("Melinder-BOOK-2010");

		Tmin     =-100.0 + 273.15;
		Tmax     =  40.0 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.3;

		Tbase    =  9.3163 + 273.15;
	    xbase    = 14.1327 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
	    	{-15.12,1124,3365,0.5461,0.9573},
			{-0.0004843,-0.3072,2.229,0.001784,-0.03065},
			{0.00001113,-0.003295,-0.004627,-0.0000008171,0.0001115},
			{0.0000001858,0.00001015,0.00009186,-0.00000006594,-0.000002923},
			{-1.885,9.071,-52.22,-0.00273,0.04968},
			{-0.00005461,-0.006513,0.1162,-0.00001483,0.0001559},
			{0.000003579,0.00004664,0.001249,0.000000385,-0.00001796},
			{-0.00000003999,0.000002287,0.000002421,-0.000000005556,0.0000003051},
			{-0.05455,0.02449,0.6202,0.000008675,-0.002722},
			{0.00001887,0.00003574,0.002337,-0.000001489,-0.00001753},
			{-0.0000003171,0.000004337,0.0000724,-0.00000003882,0.000002021},
			{-0.000000006246,0.0000006044,-0.000003613,0.0000000009282,0.00000002614},
			{-0.0007257,0.003402,0.01052,0.0000008651,0.00009351},
			{0.0000007588,-0.0001409,-0.000459,0.0000001992,-0.000008353},
			{-0.00000004102,0.000001857,-0.00002477,-0.000000001196,0.0000001901},
			{-0.0003208,0.0003344,-0.001067,-0.0000004779,0.00007364},
			{-0.0000000492,-0.00000683,-0.0001048,0.00000001797,-0.0000004014},
			{-0.00001794,0.000007239,-0.0000696,-0.00000002503,0.000003717}
		};

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = get_col(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };
};

class NASolution : public MelinderSolution{
public:
	NASolution(){

        name = std::string("MNA");
        description = std::string("Sodium Chloride (NaCl)");
        reference = std::string("Melinder-BOOK-2010");

		Tmin     =-100.0 + 273.15;
		Tmax     =  40.0 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.23;

		Tbase    = 12.6179 + 273.15;
	    xbase    = 13.3897 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
	    	{-9.383,1099,3593,0.5736,0.4369},
			{-0.00002581,-0.3758,1.669,0.001595,-0.02666},
			{0.000001423,-0.002023,-0.02019,-0.0000004003,0.0002035},
			{0,0,0,0,0},
			{-0.9039,7.723,-32.48,-0.0009383,0.02346},
			{0.000002578,-0.01426,-0.03331,-0.00001248,-0.00005368},
			{-0.00000003318,0.0001535,-0.001164,0.0000003353,0.000002871},
			{0,0,0,0,0},
			{-0.02204,0.02567,0.6453,-0.00001057,0.0004276},
			{0.0000001192,0.0003994,-0.009314,0.0000004158,-0.000004526},
			{-0.000000008993,-0.000007281,0.0002236,-0.00000002032,0.0000001838},
			{0,0,0,0,0},
			{-0.0004827,0.0001108,-0.01629,-0.0000004853,0.000007386},
			{-0.00000001467,0.0000003522,0.0007927,-0.00000001587,0.0000005437},
			{0,0,0,0,0},
			{0.000002247,-0.00001686,0.0002331,-0.000000004654,0.0000004688},
			{0,0,0,0,0},
			{0,0,0,0,0}
		};

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = get_col(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };
};

class KASolution : public MelinderSolution{
public:
	KASolution(){

        name = std::string("MKA");
        description = std::string("Potassium Acetate (CH3CO2K)");
        reference = std::string("Melinder-BOOK-2010");

		Tmin     =-100.0 + 273.15;
		Tmax     =  40.0 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.45;

		Tbase    =  6.7757 + 273.15;
	    xbase    = 25.6757 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
	    	{-17.04,1138,3327,0.4958,1.042},
			{-0.0001082,-0.3565,1.806,0.00134,-0.03071},
			{0.000006892,-0.00202,-0.001766,0.00000006837,0.0002819},
			{-0.0000001397,0.000004205,0.00004357,0.000000002637,-0.00000219},
			{-1.228,5.796,-28.95,-0.002931,0.03405},
			{0.0000002302,-0.0079,0.04131,-0.00001477,-0.0001511},
			{-0.0000008711,0.00002753,0.0004429,0.00000004659,0.000001172},
			{0.00000002016,0.0000000514,0.00001125,0.0000000002886,-0.00000002379},
			{-0.03862,0.01306,0.04663,0.00001032,0.0005017},
			{0.000001565,0.00006845,0.0007775,-0.0000002396,-0.000007779},
			{0.000000007565,-0.00000113,0.00003463,-0.000000004352,0.00000009125},
			{-0.0000000003063,0.00000002433,-0.0000007261,-0.0000000000223,-0.000000001888},
			{-0.0004571,-0.001427,-0.001249,-0.0000002024,0.000005637},
			{-0.00000003734,0.0000008304,0.00005115,-0.000000005541,0.00000002534},
			{0.000000001268,0.00000001303,-0.000002987,0.00000000008984,0.000000001596},
			{0.00002969,0.000009353,0.0001659,-0.000000002371,-0.0000002922},
			{-0.000000002817,0.00000002322,-0.000005193,0.0000000005573,0.000000004601},
			{0.000000869,0.000002285,-0.0000004612,0.0000000005515,-0.00000000796}
		};

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = get_col(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };
};

class KFSolution : public MelinderSolution{
public:
	KFSolution(){

        name = std::string("MKF");
        description = std::string("Potassium Formate (CHKO2)");
        reference = std::string("Melinder-BOOK-2010");

		Tmin     =-100.0 + 273.15;
		Tmax     =  40.0 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.48;

		Tbase    = 5.89080 + 273.15;
	    xbase    = 29.1447 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
	    	{-20.19,1189,3144,0.5253,0.8088},
			{-0.0001703,-0.3515,1.698,0.001241,-0.02556},
			{-0.000007478,-0.001918,-0.001303,-0.00000003799,0.0002195},
			{0.0000003761,0.00003132,0.00005177,-0.000000002951,-0.000001667},
			{-1.106,7.044,-29.94,-0.001972,0.01758},
			{-0.000004203,-0.00656,0.0229,-0.000006322,0.00008603},
			{-0.0000005737,0.00007018,0.000003898,0.00000002654,0.000002498},
			{0.00000001474,-0.000001459,0.000007391,-0.0000000007324,-0.00000003569},
			{-0.021,0.01889,0.2944,-0.00002596,-0.00008372},
			{-0.0000003802,0.0001034,-0.002482,-0.00000004693,-0.000001601},
			{0.00000006504,-0.000003889,0.000033,-0.000000004362,0.000000004815},
			{-0.000000001877,0.000000009225,-0.0000004871,0.00000000002899,-0.000000001861},
			{-0.0002638,-0.0002262,0.001161,-0.0000005062,-0.000001184},
			{0.00000004027,0.0000003692,0.00006758,-0.000000005229,-0.0000001331},
			{-0.0000000004789,0.00000006609,-0.000001277,0.0000000001035,-0.000000005489},
			{0.0000008327,0.00002368,-0.0001429,0.00000001501,0.000001088},
			{0.0000000009507,0.00000008766,0.0000009949,0.0000000005376,-0.000000007003},
			{0.00000006345,0.000001061,-0.000003221,0.0000000005562,0.00000003098}
		};

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = get_col(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };
};

class LISolution : public MelinderSolution{
public:
	LISolution(){

        name = std::string("MLI");
        description = std::string("Lithium Chloride (LiCl)");
        reference = std::string("Melinder-BOOK-2010");

		Tmin     =-100.0 + 273.15;
		Tmax     =  40.0 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.24;

		Tbase    =  1.4895 + 273.15;
	    xbase    = 14.8000 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
	    	{-23.29,1088,3383,0.5362,1.013},
			{0.0006555,-0.1772,3.958,0.001454,-0.03062},
			{-0.0001208,-0.002619,-0.0003035,-0.0000000326,0.000294},
			{0.000002616,0.000006209,-0.000003477,-0.0000000142,-0.000002719},
			{-3.051,6.056,-50.36,-0.001855,0.0392},
			{-0.0003972,-0.008588,0.4415,-0.00001405,0.00006246},
			{0.00003674,0.0001567,-0.0002609,-0.000000005424,-0.000001752},
			{-0.0000005569,-0.000001847,0.000003651,0.0000000009821,0.00000008346},
			{-0.179,0.02556,0.6298,0.00001017,0.000332},
			{0.00001391,0.00007194,-0.004384,0.0000006821,0.000000784},
			{-0.000001997,-0.00001053,0.0001039,-0.000000008674,-0.00000031},
			{0.00000002931,0.00000009716,-0.000001076,0.0000000001095,0.00000001381},
			{-0.002917,0.0009407,0.04544,0.0000007084,0.000002206},
			{0.00000149,-0.000007253,-0.0008787,-0.00000007434,-0.0000006011},
			{-0.00000006904,0.0000003144,-0.000008457,0.0000000006988,0.000000004023},
			{0.0005715,-0.00008105,0.002527,-0.0000001273,0.000001745},
			{0.0000001186,0.000001072,-0.0000358,-0.000000003058,-0.00000007094},
			{0.00002757,-0.000003974,0.00004058,-0.000000009124,0.00000006699}
		};

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = get_col(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };
};









class SecCoolSolution : public BaseSolution{
public:
	SecCoolSolution(){

		std::vector<double> tmpVector;

        name = std::string("SecCoolSolution");
        description = std::string("Test Methanol SecCool");
        reference = std::string("Test");

		Tmin     = -50 + 273.15;
		Tmax     =  20 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.5;

		Tbase    = -4.48 + 273.15;
	    xbase    = 31.57 / 100.0;

	    tmpVector.clear();
	    tmpVector.push_back( 960.24665800);
	    tmpVector.push_back(-1.2903839100);
	    tmpVector.push_back(-0.0161042520);
	    tmpVector.push_back(-0.0001969888);
	    tmpVector.push_back( 1.131559E-05);
	    tmpVector.push_back( 9.181999E-08);
	    tmpVector.push_back(-0.4020348270);
	    tmpVector.push_back(-0.0162463989);
	    tmpVector.push_back( 0.0001623301);
	    tmpVector.push_back( 4.367343E-06);
	    tmpVector.push_back( 1.199000E-08);
	    tmpVector.push_back(-0.0025204776);
	    tmpVector.push_back( 0.0001101514);
	    tmpVector.push_back(-2.320217E-07);
	    tmpVector.push_back( 7.794999E-08);
	    tmpVector.push_back( 9.937483E-06);
	    tmpVector.push_back(-1.346886E-06);
	    tmpVector.push_back( 4.141999E-08);
        cRho.clear();
        cRho = makeMatrix(tmpVector);

        tmpVector.clear();
        tmpVector.push_back( 3822.9712300);
        tmpVector.push_back(-23.122409500);
        tmpVector.push_back( 0.0678775826);
        tmpVector.push_back( 0.0022413893);
        tmpVector.push_back(-0.0003045332);
        tmpVector.push_back(-4.758000E-06);
        tmpVector.push_back( 2.3501449500);
        tmpVector.push_back( 0.1788839410);
        tmpVector.push_back( 0.0006828000);
        tmpVector.push_back( 0.0002101166);
        tmpVector.push_back(-9.812000E-06);
        tmpVector.push_back(-0.0004724176);
        tmpVector.push_back(-0.0003317949);
        tmpVector.push_back( 0.0001002032);
        tmpVector.push_back(-5.306000E-06);
        tmpVector.push_back( 4.242194E-05);
        tmpVector.push_back( 2.347190E-05);
        tmpVector.push_back(-1.894000E-06);
		cHeat.clear();
		cHeat = makeMatrix(tmpVector);

		tmpVector.clear();
		tmpVector.push_back( 0.4082066700);
		tmpVector.push_back(-0.0039816870);
		tmpVector.push_back( 1.583368E-05);
		tmpVector.push_back(-3.552049E-07);
		tmpVector.push_back(-9.884176E-10);
		tmpVector.push_back( 4.460000E-10);
		tmpVector.push_back( 0.0006629321);
		tmpVector.push_back(-2.686475E-05);
		tmpVector.push_back( 9.039150E-07);
		tmpVector.push_back(-2.128257E-08);
		tmpVector.push_back(-5.562000E-10);
		tmpVector.push_back( 3.685975E-07);
		tmpVector.push_back( 7.188416E-08);
		tmpVector.push_back(-1.041773E-08);
		tmpVector.push_back( 2.278001E-10);
		tmpVector.push_back( 4.703395E-08);
		tmpVector.push_back( 7.612361E-11);
		tmpVector.push_back(-2.734000E-10);
		cCond.clear();
		cCond = makeMatrix(tmpVector);

		tmpVector.clear();
		tmpVector.push_back( 1.4725525500);
		tmpVector.push_back( 0.0022218998);
		tmpVector.push_back(-0.0004406139);
		tmpVector.push_back( 6.047984E-06);
		tmpVector.push_back(-1.954730E-07);
		tmpVector.push_back(-2.372000E-09);
		tmpVector.push_back(-0.0411841566);
		tmpVector.push_back( 0.0001784479);
		tmpVector.push_back(-3.564413E-06);
		tmpVector.push_back( 4.064671E-08);
		tmpVector.push_back( 1.915000E-08);
		tmpVector.push_back( 0.0002572862);
		tmpVector.push_back(-9.226343E-07);
		tmpVector.push_back(-2.178577E-08);
		tmpVector.push_back(-9.529999E-10);
		tmpVector.push_back(-1.699844E-06);
		tmpVector.push_back(-1.023552E-07);
		tmpVector.push_back( 4.482000E-09);
		cVisc.clear();
		cVisc = makeMatrix(tmpVector);

		cTfreeze.clear();
		cTfreeze.push_back( 27.755555600); // reference concentration in per cent
		cTfreeze.push_back(-22.973221700);
		cTfreeze.push_back(-1.1040507200);
		cTfreeze.push_back(-0.0120762281);
		cTfreeze.push_back(-9.343458E-05);

    };

	/// Define freezing point calculations
	double Tfreeze(double p, double x){
		IncompressibleClass::checkCoefficients(cTfreeze,5);
		std::vector<double> tmpVector(cTfreeze.begin()+1,cTfreeze.end());
		return polyval(tmpVector, x*100.0-cTfreeze[0])+273.15;
	}
};


class ZitrecAC : public SecCoolSolution{
public:
	ZitrecAC(){

		std::vector<double> tmpVector;

        name = std::string("ZiAC");
        description = std::string("ZitrecAC in water (corrosion inhibitor)");
        reference = std::string("SecCool Software");

		Tmin     =   0 + 273.15;
		Tmax     = 100 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.05;
		xmax     = 0.50;

		Tbase    = 50.00 + 273.15;
	    xbase    = 22.75 / 100.0;

	    tmpVector.clear();
	    tmpVector.push_back( 1003.4314200);
	    tmpVector.push_back( 0.6164672840);
	    tmpVector.push_back(-0.0075340011);
	    tmpVector.push_back( 8.227043E-05);
	    tmpVector.push_back( 1.416356E-05);
	    tmpVector.push_back(-3.611589E-07);
	    tmpVector.push_back(-0.4849253070);
	    tmpVector.push_back(-0.0015163769);
	    tmpVector.push_back(-3.076387E-05);
	    tmpVector.push_back(-4.631673E-07);
	    tmpVector.push_back( 3.958918E-08);
	    tmpVector.push_back(-0.0053161289);
	    tmpVector.push_back(-3.675038E-06);
	    tmpVector.push_back( 2.245807E-06);
	    tmpVector.push_back(-5.921160E-08);
	    tmpVector.push_back(-2.222469E-05);
	    tmpVector.push_back(-1.016950E-06);
	    tmpVector.push_back( 1.391098E-08);
        cRho.clear();
        cRho = makeMatrix(tmpVector);

        tmpVector.clear();
        tmpVector.push_back( 4129.8221400);
        tmpVector.push_back(-2.4602756000);
        tmpVector.push_back( 0.0024608374);
        tmpVector.push_back(-0.0018878988);
        tmpVector.push_back(-7.525318E-05);
        tmpVector.push_back( 3.812299E-06);
        tmpVector.push_back( 0.3207990110);
        tmpVector.push_back(-0.0059380735);
        tmpVector.push_back(-0.0007210516);
        tmpVector.push_back( 1.852057E-05);
        tmpVector.push_back( 1.280127E-07);
        tmpVector.push_back( 0.0099912322);
        tmpVector.push_back( 0.0001857846);
        tmpVector.push_back(-2.266658E-06);
        tmpVector.push_back(-5.937414E-09);
        tmpVector.push_back(-0.0001641590);
        tmpVector.push_back(-1.708329E-06);
        tmpVector.push_back( 8.577111E-08);
		cHeat.clear();
		cHeat = makeMatrix(tmpVector);

		tmpVector.clear();
		tmpVector.push_back( 0.5997434470);
        tmpVector.push_back(-0.0014518043);
        tmpVector.push_back( 7.240125E-06);
        tmpVector.push_back(-6.458395E-07);
        tmpVector.push_back(-1.035513E-08);
        tmpVector.push_back( 6.862018E-10);
        tmpVector.push_back( 0.0011059284);
        tmpVector.push_back(-1.686714E-06);
        tmpVector.push_back(-1.542610E-08);
        tmpVector.push_back(-3.730822E-09);
        tmpVector.push_back( 1.286049E-11);
        tmpVector.push_back(-5.567860E-06);
        tmpVector.push_back( 5.331112E-08);
        tmpVector.push_back(-1.738845E-09);
        tmpVector.push_back( 2.641378E-11);
        tmpVector.push_back(-1.195582E-08);
        tmpVector.push_back(-3.329696E-10);
        tmpVector.push_back( 2.343560E-11);
		cCond.clear();
		cCond = makeMatrix(tmpVector);

		tmpVector.clear();
		tmpVector.push_back(-0.2833999350);
        tmpVector.push_back( 0.0113457408);
        tmpVector.push_back(-0.0002173449);
        tmpVector.push_back( 9.645169E-06);
        tmpVector.push_back( 6.299118E-07);
        tmpVector.push_back(-2.190184E-08);
        tmpVector.push_back(-0.0185725092);
        tmpVector.push_back(-8.290652E-05);
        tmpVector.push_back(-2.186512E-06);
        tmpVector.push_back( 1.126354E-07);
        tmpVector.push_back(-2.387450E-09);
        tmpVector.push_back( 0.0001256061);
        tmpVector.push_back( 1.147598E-06);
        tmpVector.push_back( 1.128849E-08);
        tmpVector.push_back(-1.723001E-09);
        tmpVector.push_back( 7.142151E-07);
        tmpVector.push_back(-5.140398E-09);
        tmpVector.push_back( 4.535194E-10);
		cVisc.clear();
		cVisc = makeMatrix(tmpVector);

		cTfreeze.clear();
		cTfreeze.push_back( 22.750000000); // reference concentration in per cent
		cTfreeze.push_back(-2.2469093100);
		cTfreeze.push_back(-0.0942887708);
		cTfreeze.push_back( 0.0002636562);
		cTfreeze.push_back( 9.008030E-07);

    };

};


class IceSlurryEA : public SecCoolSolution{
public:
	IceSlurryEA(){

		std::vector<double> tmpVector;

        name = std::string("IceEA");
        description = std::string("Ethanol-water mixture with slurry ice");
        reference = std::string("SecCool Software");

		Tmin     = -35 + 273.15;
		Tmax     = -10 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.05;
		xmax     = 0.35;

		Tbase    = -22.5 + 273.15;
	    xbase    =  20.0 / 100.0;

	    tmpVector.clear();
	    tmpVector.push_back( 959.65328700);
	    tmpVector.push_back(-0.6063295290);
	    tmpVector.push_back(-1.249095E-05);
	    tmpVector.push_back(-2.175929E-05);
	    tmpVector.push_back( 1.414141E-06);
	    tmpVector.push_back( 8.888899E-08);
	    tmpVector.push_back( 0.4964468440);
	    tmpVector.push_back(-0.0048152589);
	    tmpVector.push_back(-3.703183E-05);
	    tmpVector.push_back( 7.619048E-07);
	    tmpVector.push_back( 1.454545E-07);
	    tmpVector.push_back(-0.0105000000);
	    tmpVector.push_back( 0.0001660431);
	    tmpVector.push_back( 1.020408E-07);
	    tmpVector.push_back(-1.587301E-07);
	    tmpVector.push_back( 0.0003601411);
	    tmpVector.push_back(-4.814815E-06);
	    tmpVector.push_back( 1.763668E-08);
        cRho.clear();
        cRho = makeMatrix(tmpVector);

        tmpVector.clear();
        tmpVector.push_back( 8.241062E+04);
	    tmpVector.push_back( 4133.4319000);
	    tmpVector.push_back(-0.0525342683);
	    tmpVector.push_back(-0.0332405251);
	    tmpVector.push_back(-0.0007070707);
	    tmpVector.push_back( -9.555626E-05);
	    tmpVector.push_back(-175.73225900);
	    tmpVector.push_back(-19.750427100);
	    tmpVector.push_back(-0.5552331250);
	    tmpVector.push_back(-0.0024634921);
	    tmpVector.push_back( 0.0001094372);
	    tmpVector.push_back(-20.031632700);
	    tmpVector.push_back(-1.1448185800);
	    tmpVector.push_back(-0.0049591837);
	    tmpVector.push_back( 0.0003301587);
	    tmpVector.push_back(-1.3004938300);
	    tmpVector.push_back(-0.0241058201);
	    tmpVector.push_back( 0.0012176367);
		cHeat.clear();
		cHeat = makeMatrix(tmpVector);

		tmpVector.clear();
		tmpVector.push_back( 0.5466734710);
	    tmpVector.push_back( 0.0098512311);
	    tmpVector.push_back( 6.916942E-05);
	    tmpVector.push_back( 4.423616E-07);
	    tmpVector.push_back( 1.313131E-09);
	    tmpVector.push_back( 1.333317E-10);
	    tmpVector.push_back( 0.0064403362);
	    tmpVector.push_back( 6.261462E-05);
	    tmpVector.push_back(-1.361660E-07);
	    tmpVector.push_back(-4.190476E-09);
	    tmpVector.push_back( 9.350649E-11);
	    tmpVector.push_back( 1.914626E-05);
	    tmpVector.push_back(-4.656462E-07);
	    tmpVector.push_back(-5.544218E-09);
	    tmpVector.push_back(-4.761923E-11);
	    tmpVector.push_back( 1.093474E-08);
	    tmpVector.push_back(-1.322751E-09);
	    tmpVector.push_back( 1.763668E-11);
		cCond.clear();
		cCond = makeMatrix(tmpVector);

		tmpVector.clear();
		tmpVector.push_back( 3.9714277000);
	    tmpVector.push_back( 0.0217778490);
	    tmpVector.push_back(-0.0002520703);
	    tmpVector.push_back( 3.339905E-06);
	    tmpVector.push_back(-6.179120E-08);
	    tmpVector.push_back( 2.654128E-09);
	    tmpVector.push_back(-0.0820149233);
	    tmpVector.push_back( 3.060191E-06);
	    tmpVector.push_back(-1.789494E-07);
	    tmpVector.push_back(-3.463500E-09);
	    tmpVector.push_back(-5.507847E-11);
	    tmpVector.push_back(-0.0010965961);
	    tmpVector.push_back( 5.631692E-08);
	    tmpVector.push_back( 1.027114E-08);
	    tmpVector.push_back(-8.462427E-10);
	    tmpVector.push_back(1.855647E-05);
	    tmpVector.push_back( 6.616363E-09);
	    tmpVector.push_back( 1.185315E-10);
		cVisc.clear();
		cVisc = makeMatrix(tmpVector);

		cTfreeze.clear();
    };
	/// Define freezing point calculations
	double Tfreeze(double p, double x){
		return Tmin;
	}
};


class IceSlurryPG : public SecCoolSolution{
public:
	IceSlurryPG(){

		std::vector<double> tmpVector;

        name = std::string("IcePG");
        description = std::string("Propylene glycol-water mixture with slurry ice");
        reference = std::string("SecCool Software");

		Tmin     = -45 + 273.15;
		Tmax     = -10 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.05;
		xmax     = 0.35;

		Tbase    = -27.5 + 273.15;
	    xbase    =  20.0 / 100.0;

	    tmpVector.clear();
	    tmpVector.push_back(1026.0807100);
	    tmpVector.push_back(-1.5901706300);
	    tmpVector.push_back( 0.0023632306);
	    tmpVector.push_back(-5.555564E-05);
	    tmpVector.push_back( 3.787877E-07);
	    tmpVector.push_back( 1.500003E-07);
	    tmpVector.push_back(-0.8565823200);
	    tmpVector.push_back( 0.0158439454);
	    tmpVector.push_back(-4.885161E-05);
	    tmpVector.push_back( 5.820106E-07);
	    tmpVector.push_back(-1.298701E-08);
	    tmpVector.push_back(-0.0205963719);
	    tmpVector.push_back( 0.0002962207);
	    tmpVector.push_back(-3.287982E-07);
	    tmpVector.push_back( 2.645504E-08);
	    tmpVector.push_back(-0.0002976431);
	    tmpVector.push_back( 3.477633E-06);
	    tmpVector.push_back( 1.827802E-08);
        cRho.clear();
        cRho = makeMatrix(tmpVector);

        tmpVector.clear();
        tmpVector.push_back( 9.202807E+04);
	    tmpVector.push_back( 4259.8599100);
	    tmpVector.push_back(-13.425892900);
	    tmpVector.push_back( 0.1972224350);
	    tmpVector.push_back( 0.0016666667);
	    tmpVector.push_back( 2.333263E-05);
	    tmpVector.push_back(-1347.4964300);
	    tmpVector.push_back(-40.622180100);
	    tmpVector.push_back( 1.2466366000);
	    tmpVector.push_back(-0.0089735450);
	    tmpVector.push_back(-0.0002542569);
	    tmpVector.push_back( 9.1439909300);
	    tmpVector.push_back( 0.3821466440);
	    tmpVector.push_back(-0.0142154195);
	    tmpVector.push_back(-0.0004074074);
	    tmpVector.push_back( 0.4962000960);
	    tmpVector.push_back( 0.0010591631);
	    tmpVector.push_back(-0.0011599808);
		cHeat.clear();
		cHeat = makeMatrix(tmpVector);

		tmpVector.clear();
		tmpVector.push_back( 0.5308857010);
	    tmpVector.push_back( 0.0102779335);
	    tmpVector.push_back( 6.990287E-05);
	    tmpVector.push_back( 4.527783E-07);
	    tmpVector.push_back( 3.030303E-09);
	    tmpVector.push_back( 9.999831E-11);
	    tmpVector.push_back( 0.0037803301);
	    tmpVector.push_back( 3.128644E-05);
	    tmpVector.push_back(-1.033550E-07);
	    tmpVector.push_back(-3.809524E-09);
	    tmpVector.push_back( 4.617604E-11);
	    tmpVector.push_back( 8.588549E-05);
	    tmpVector.push_back( 3.319350E-07);
	    tmpVector.push_back(-5.147392E-09);
	    tmpVector.push_back(-1.164022E-10);
	    tmpVector.push_back( 2.443001E-06);
	    tmpVector.push_back( 8.542569E-09);
	    tmpVector.push_back(-1.558442E-10);
		cCond.clear();
		cCond = makeMatrix(tmpVector);

		tmpVector.clear();
		tmpVector.push_back( 5.4382616400);
	    tmpVector.push_back( 0.0178821732);
	    tmpVector.push_back(-0.0001873851);
	    tmpVector.push_back( 2.395101E-06);
	    tmpVector.push_back(-3.751781E-08);
	    tmpVector.push_back( 5.875223E-10);
	    tmpVector.push_back(-0.1513106190);
	    tmpVector.push_back(-6.271586E-06);
	    tmpVector.push_back( 3.964457E-07);
	    tmpVector.push_back(-1.180828E-08);
	    tmpVector.push_back( 9.672800E-11);
	    tmpVector.push_back(-0.0002684929);
	    tmpVector.push_back(-1.099919E-07);
	    tmpVector.push_back( 6.150740E-09);
	    tmpVector.push_back(-1.902903E-10);
	    tmpVector.push_back(-4.725022E-06);
	    tmpVector.push_back(-1.027292E-10);
	    tmpVector.push_back( 8.515353E-11);
		cVisc.clear();
		cVisc = makeMatrix(tmpVector);

		cTfreeze.clear();
    };
	/// Define freezing point calculations
	double Tfreeze(double p, double x){
		return Tmin;
	}
};


class IceSlurryNA : public SecCoolSolution{
public:
	IceSlurryNA(){

		std::vector<double> tmpVector;

        name = std::string("IceNA");
        description = std::string("Sodium chloride-water mixture with slurry ice");
        reference = std::string("SecCool Software");

		Tmin     = -20 + 273.15;
		Tmax     =  -5 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.05;
		xmax     = 0.35;

		Tbase    = -12.5 + 273.15;
	    xbase    =  20.0 / 100.0;

	    tmpVector.clear();
	    tmpVector.push_back( 1081.6353100);
	    tmpVector.push_back(-2.4559523700);
	    tmpVector.push_back( 0.0058152057);
	    tmpVector.push_back(-7.500013E-05);
	    tmpVector.push_back(-7.575759E-07);
	    tmpVector.push_back( 1.666671E-07);
	    tmpVector.push_back(-5.6609963900);
	    tmpVector.push_back( 0.1002726190);
	    tmpVector.push_back(-0.0004797330);
	    tmpVector.push_back( 1.333333E-06);
	    tmpVector.push_back( 3.636364E-08);
	    tmpVector.push_back(-0.0852857143);
	    tmpVector.push_back( 0.0007904762);
	    tmpVector.push_back( 1.428571E-06);
	    tmpVector.push_back( 6.666668E-07);
	    tmpVector.push_back(-0.0037650794);
	    tmpVector.push_back( 3.333333E-05);
	    tmpVector.push_back( 6.984127E-07);
        cRho.clear();
        cRho = makeMatrix(tmpVector);

        tmpVector.clear();
        tmpVector.push_back( 7.434384E+04);
	    tmpVector.push_back( 3669.8467100);
	    tmpVector.push_back(-2.0844426400);
	    tmpVector.push_back( 0.0312501929);
	    tmpVector.push_back(-0.0002727273);
	    tmpVector.push_back(-8.333396E-05);
	    tmpVector.push_back(-794.24689800);
	    tmpVector.push_back(-33.895515900);
	    tmpVector.push_back( 0.3610772000);
	    tmpVector.push_back(-0.0016888889);
	    tmpVector.push_back( 0.0001406061);
	    tmpVector.push_back( 12.209523800);
	    tmpVector.push_back( 0.3702381290);
	    tmpVector.push_back(-0.0099523810);
	    tmpVector.push_back( 0.0001333331);
	    tmpVector.push_back(-0.1358730160);
	    tmpVector.push_back( 0.0145714286);
	    tmpVector.push_back(-0.0014412698);
		cHeat.clear();
		cHeat = makeMatrix(tmpVector);

		tmpVector.clear();
		tmpVector.push_back( 0.7579141770);
	    tmpVector.push_back( 0.0124563700);
	    tmpVector.push_back( 5.749080E-05);
	    tmpVector.push_back( 2.263889E-07);
	    tmpVector.push_back(-7.575758E-10);
	    tmpVector.push_back( 1.333333E-10);
	    tmpVector.push_back( 0.0009894098);
	    tmpVector.push_back(-5.386429E-05);
	    tmpVector.push_back( 2.049928E-07);
	    tmpVector.push_back( 1.333333E-09);
	    tmpVector.push_back(-3.757576E-10);
	    tmpVector.push_back( 8.761905E-06);
	    tmpVector.push_back(-9.531746E-07);
	    tmpVector.push_back( 3.809524E-09);
	    tmpVector.push_back( 2.222222E-10);
	    tmpVector.push_back(-9.777778E-07);
	    tmpVector.push_back(-5.904762E-08);
	    tmpVector.push_back(-1.269841E-10);
		cCond.clear();
		cCond = makeMatrix(tmpVector);

		tmpVector.clear();
		tmpVector.push_back( 1.9270346900);
	    tmpVector.push_back( 0.0216118832);
	    tmpVector.push_back(-0.0002820062);
	    tmpVector.push_back( 4.476720E-06);
	    tmpVector.push_back(-9.032034E-08);
	    tmpVector.push_back( 2.020892E-09);
	    tmpVector.push_back(-0.0713269985);
	    tmpVector.push_back(-2.729779E-05);
	    tmpVector.push_back( 2.115509E-06);
	    tmpVector.push_back(-7.777973E-08);
	    tmpVector.push_back( 2.272317E-09);
	    tmpVector.push_back( 0.0003930707);
	    tmpVector.push_back(-7.171429E-07);
	    tmpVector.push_back( 4.470017E-08);
	    tmpVector.push_back(-1.194034E-09);
	    tmpVector.push_back( 2.156903E-05);
	    tmpVector.push_back(-1.127279E-08);
	    tmpVector.push_back( 3.086759E-09);
		cVisc.clear();
		cVisc = makeMatrix(tmpVector);

		cTfreeze.clear();
    };
	/// Define freezing point calculations
	double Tfreeze(double p, double x){
		return Tmin;
	}
};


class PK2000 : public SecCoolSolution{
protected:
	std::vector<double> cMaVo;
public:
	PK2000(){

		std::vector<double> tmpVector;

        name = std::string("PK2000");
        description = std::string("Pekasol 2000 in water (Potassium acetate and formate)");
        reference = std::string("SecCool Software");

		Tmin     = -62 + 273.15;
		Tmax     = 100 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.36;
		xmax     = 1.00;

		Tbase    = 33.31 + 273.15;
	    xbase    = 67.60 / 100.0; // volume percent!

	    tmpVector.clear();
	    tmpVector.push_back( 1197.8306300);
	    tmpVector.push_back( 2.7580390000);
	    tmpVector.push_back(-0.0046328716);
	    tmpVector.push_back(-2.118894E-05);
	    tmpVector.push_back( 5.174717E-08);
	    tmpVector.push_back( 2.265693E-09);
	    tmpVector.push_back(-0.5038076360);
	    tmpVector.push_back(-0.0020754530);
	    tmpVector.push_back( 7.670317E-06);
	    tmpVector.push_back( 1.443587E-08);
	    tmpVector.push_back(-1.451878E-09);
	    tmpVector.push_back(-0.0015712351);
	    tmpVector.push_back( 2.546509E-05);
	    tmpVector.push_back(-8.010506E-08);
	    tmpVector.push_back(-4.948979E-10);
	    tmpVector.push_back( 2.947414E-08);
	    tmpVector.push_back(-4.747692E-09);
	    tmpVector.push_back( 8.927986E-11);
        cRho.clear();
        cRho = makeMatrix(tmpVector);

        tmpVector.clear();
        tmpVector.push_back( 3012.5363200);
	    tmpVector.push_back(-11.345089400);
	    tmpVector.push_back( 0.0512475571);
	    tmpVector.push_back(-0.0004261743);
	    tmpVector.push_back( 5.582778E-06);
	    tmpVector.push_back( 2.339332E-08);
	    tmpVector.push_back( 1.8771264300);
	    tmpVector.push_back(-0.0024687047);
	    tmpVector.push_back(-0.0001136660);
	    tmpVector.push_back(-3.281309E-06);
	    tmpVector.push_back(-2.598223E-08);
	    tmpVector.push_back(-0.0117207757);
	    tmpVector.push_back( 0.0001728598);
	    tmpVector.push_back( 2.221588E-06);
	    tmpVector.push_back(-1.099247E-07);
	    tmpVector.push_back( 6.566134E-05);
	    tmpVector.push_back( 3.261243E-06);
	    tmpVector.push_back(-5.841138E-08);
		cHeat.clear();
		cHeat = makeMatrix(tmpVector);

		tmpVector.clear();
		tmpVector.push_back( 0.5060902210);
	    tmpVector.push_back(-0.0015058953);
	    tmpVector.push_back( 4.296707E-06);
	    tmpVector.push_back(-1.171421E-09);
	    tmpVector.push_back(-5.116965E-12);
	    tmpVector.push_back( 8.612329E-14);
	    tmpVector.push_back( 0.0009077322);
	    tmpVector.push_back(-4.052813E-06);
	    tmpVector.push_back( 8.942927E-09);
	    tmpVector.push_back(-3.126326E-11);
	    tmpVector.push_back( 1.729361E-12);
	    tmpVector.push_back( 2.222018E-09);
	    tmpVector.push_back(-4.276078E-10);
	    tmpVector.push_back(-3.160239E-11);
	    tmpVector.push_back( 1.220910E-12);
	    tmpVector.push_back( 3.600639E-10);
	    tmpVector.push_back(-3.839049E-11);
	    tmpVector.push_back( 7.540478E-13);
		cCond.clear();
		cCond = makeMatrix(tmpVector);

		tmpVector.clear();
	    tmpVector.push_back( 0.5631031820);
	    tmpVector.push_back( 0.0147645079);
	    tmpVector.push_back(-1.024523E-05);
	    tmpVector.push_back(-7.780506E-07);
	    tmpVector.push_back( 1.243805E-08);
	    tmpVector.push_back( 7.494764E-12);
	    tmpVector.push_back(-0.0199399025);
	    tmpVector.push_back(-1.895556E-05);
	    tmpVector.push_back(-8.956373E-07);
	    tmpVector.push_back(-1.780877E-08);
	    tmpVector.push_back( 9.735698E-11);
	    tmpVector.push_back( 0.0001153086);
	    tmpVector.push_back( 2.823599E-07);
	    tmpVector.push_back( 7.608202E-09);
	    tmpVector.push_back(-3.413081E-10);
	    tmpVector.push_back(-1.580556E-06);
	    tmpVector.push_back(-1.693264E-08);
	    tmpVector.push_back(-1.648529E-10);
		cVisc.clear();
		cVisc = makeMatrix(tmpVector);

		cTfreeze.clear();
		cTfreeze.push_back( 65.000000000); // reference concentration in per cent
		cTfreeze.push_back(-28.549321300);
		cTfreeze.push_back(-0.7312947390);
		cTfreeze.push_back(-0.0061085973);
		cTfreeze.push_back(-8.714598E-06);

		cMaVo.clear();
		cMaVo.push_back( 63.693536900);
		cMaVo.push_back( 1.0864767400);
		cMaVo.push_back( 0.0033972173);
		cMaVo.push_back( 1.986130E-05);

    };

	// Redefine x_m to convert from mass to volume fraction
	// M1 = M-M_m
	// Vol% = A[1] + A[2]*M1 + A[3]*M1^2 + A[4]*M1^3
	double getxInput(double curxValue){
		double xVolume = polyval(cMaVo, (curxValue*100.-69.92));
		return xVolume-xbase*100.0;
	}


};

#endif
