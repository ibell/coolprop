
#ifndef INCOMPRESSIBLE_SOLUTION_H
#define INCOMPRESSIBLE_SOLUTION_H

#include <string>
#include <vector>
#include "CPExceptions.h"
#include "CoolProp.h"
#include "IncompLiquid.h"
#include <math.h>
#include "Solvers.h"

/**
Notes for developers:

If you want to add a fluid, add its definition to the header 
and then add it to the map in the constructor in SolutionsContainer
in IncompSolution.cpp
**/

/// Base class for simplified brine/solution models
/** Employs the base functions implemented in IncompLiquid.h.
 *  Extends the functions for composition as input. */
class IncompressibleSolution : public IncompressibleLiquid{

protected:
	double xmin, xmax;

public:
	// Constructor
	IncompressibleSolution(){
		xmin = -1.;
		xmax = -1.;
	};

	// Destructor, no implementation
	virtual ~IncompressibleSolution(){};

protected:
	/// Enthalpy from x, u, p and rho.
	/** Calculate enthalpy as a function of temperature and
	 *  pressure employing functions for internal energy and
	 *  density. Provides consistent formulations. */
	double h_u(double T_K, double p) {needComposition();return -_HUGE;};
	double h_u(double T_K, double p, double x) {
		return u(T_K,p,x)+p/rho(T_K,p,x);
	};

	/// Internal energy from x, h, p and rho.
	/** Calculate internal energy as a function of temperature
	 *  and pressure employing functions for enthalpy and
	 *  density. Provides consistent formulations. */
	double u_h(double T_K, double p) {needComposition();return -_HUGE;};
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
	 *  defined for all fluids, default values do not
	 *  cause errors. */
	bool checkT(double T_K, double x){
		if( Tmin < 0. ) {
			throw ValueError("Please specify the minimum temperature.");
		} else if( Tmax < 0.) {
			throw ValueError("Please specify the maximum temperature.");
		} else if ( (Tmin>T_K) || (T_K>Tmax) ) {
			throw ValueError(format("Your temperature %f is not between %f and %f.",T_K,Tmin,Tmax));
		} else {
			return true;
		}
		return false;
	}

	/// Check validity of pressure input.
	/** Compares the given pressure p to the saturation pressure at
	 *  temperature T and throws and exception if p is lower than
	 *  the saturation conditions.
	 *  The default value for psat is -1 yielding true is psat
	 *  is not redefined in the subclass.
	 *  */
	bool checkP(double T_K, double p) {
		double ps = psat(T_K);
		if (p<ps) {
			throw ValueError(format("Equations are valid for liquid phase only: %f < %f. ",p,ps));
		} else {
			return true;
		}
	}

	/// Check validity of temperature and pressure input.
	bool checkTP(double T, double p) {
		return (checkT(T) && checkP(T,p));
	}

	/// 2D polynomial function generator.
	/** Base function to produce n-th order polynomials
	 *  based on the size of the coefficient matrix.
	 *  Starts with only the first coefficient at x^0*y^0. */
	double basePolynomial(std::vector<std::vector<double>> coefficients, double x, double y){
		double result = 0.;
		for(unsigned int i=0; i<coefficients.size();i++) {
			result += pow(x,(int)i) * IncompressibleLiquid::basePolynomial(coefficients[i],y);
		}
		return result;
	}

	/// 2D polynomial function generator with check.
	/** Base function to produce i-th order polynomials
	 *  based on the size of the coefficient matrix.
	 *  Starts with only the first coefficient at x^0*y^0
	 *  and checks the vector length against parameter n and o. */
	double basePolynomial(std::vector<std::vector<double>> coefficients, double x, double y, unsigned int n, unsigned int o){
		if (coefficients.size() == n){
			double result = 0.;
			for(unsigned int i=0; i<n; i++) {
				result += pow(x,(int)i) * IncompressibleLiquid::basePolynomial(coefficients[i],y,o);
			}
			return result;
		} else {
			throw ValueError(format("The number of rows %d does not match with %d. ",coefficients.size(),n));
		}
	}

	/// Function used to enforce the composition as parameter
	/** Overwrites the normal functions that only take 2
	 *  parameters (T,p) and throws an exception. */
	bool needComposition() {
		throw ValueError(format("The fluid %s needs an additional input for the composition.",this->name));
		return false;
	}

public:
	/// First we disable all the standard functions ...
	double rho(double T_K, double p){
		needComposition();
		return -_HUGE;
	}
	double cp(double T_K, double p){
		needComposition();
		return -_HUGE;
	}
	double h(double T_K, double p){
		needComposition();
		return -_HUGE;
	}
	double s(double T_K, double p){
		needComposition();
		return -_HUGE;
	}
	double visc(double T_K, double p){
		needComposition();
		return -_HUGE;
	}
	double cond(double T_K, double p){
		needComposition();
		return -_HUGE;
	}
	double u(double T_K, double p){
		needComposition();
		return -_HUGE;
	}
	double psat(double T_K){
		needComposition();
		return -_HUGE;
	};
	double Tsat(double p){
		needComposition();
		return -_HUGE;
	};

	/// ... and then we define new ones with a composition parameter
    virtual double rho (double T_K, double p, double x);
    virtual double cp  (double T_K, double p, double x);
    virtual double h   (double T_K, double p, double x);
    virtual double s   (double T_K, double p, double x);
    virtual double visc(double T_K, double p, double x);
    virtual double cond(double T_K, double p, double x);
    virtual double u   (double T_K, double p, double x);
    virtual double psat(double T_K          , double x);
    virtual double Tsat(            double p, double x);
};

/// Class to use the SecCool parameters
/** Employs some basic wrapper-like functionality
 *  to bridge the gap between the solution functions
 *  used in CoolProp and the definition used in
 *  SecCool. Please visit:
 *  http://en.ipu.dk/Indhold/refrigeration-and-energy-technology/seccool.aspx
 *  A big thanks to Morten Juel Skovrup for providing
 *  this nice piece of software as well as the parameters
 *  needed to calculate the composition based properties. */
class SecCoolSolution : public IncompressibleSolution{

protected:
	double Tbase, xbase;

public:
	// Constructor
	SecCoolSolution(){
		Tbase = -1.;
		xbase = -1.;
	};

	// Destructor, no implementation
	virtual ~SecCoolSolution(){};

protected:
	double getTInput(double curTValue){
		return curTValue-Tbase;
	}

	double getxInput(double curxValue){
		return curxValue-xbase;
	}

};



#endif
