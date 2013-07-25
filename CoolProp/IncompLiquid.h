
#ifndef INCOMPRESSIBLE_LIQUID_H
#define INCOMPRESSIBLE_LIQUID_H

#include <string>
#include <vector>
#include "CPExceptions.h"
#include "CoolProp.h"
#include <math.h>
#include "Solvers.h"

/**
Notes for developers:

If you want to add a fluid, add its definition to the header 
and then add it to the map in the constructor in LiquidsContainer 
in IncompLiquid.cpp
**/

/// The abstract base class for the liquids
class IncompressibleLiquid{
	
protected:
	std::string name,description,reference;
	double Tmin, TminPsat, Tmax;

public:
	// Constructor
	IncompressibleLiquid(){
		Tmin = -1.;
		Tmax = -1.;
		TminPsat = -1.;
	};

	// Destructor.  No implementation
	virtual ~IncompressibleLiquid(){};

	/* All functions need T and p as input. Might not be necessary,
	 * but gives a clearer structure.
	 */

	/// Density as a function of temperature and pressure.
	virtual double rho (double T_K, double p){return 0;};
	/// Isobaric heat capacity as a function of temperature and pressure.
	virtual double cp  (double T_K, double p){return 0;};
	/// Entropy as a function of temperature and pressure.
	virtual double s   (double T_K, double p){return 0;};
	/// Internal energy as a function of temperature and pressure.
	virtual double u   (double T_K, double p){return 0;};
	/// Enthalpy as a function of temperature and pressure.
	virtual double h   (double T_K, double p){return 0;};
	/// Viscosity as a function of temperature and pressure.
	virtual double visc(double T_K, double p){return 0;};
	/// Thermal conductivity as a function of temperature and pressure.
	virtual double cond(double T_K, double p){return 0;};
	/// Saturation pressure as a function of temperature.
	virtual double psat(double T_K){return 0;};
	/// Saturation temperature as a function of pressure.
	virtual double Tsat(double p){return 0;};

	std::string get_name() {
		return name;
	}

protected:
	/* Define internal energy and enthalpy as functions of the
	 * other properties to provide data in case there are no
	 * coefficients.
	 */

	/// Enthalpy from u,p and rho.
	/** Calculate enthalpy as a function of temperature and
	 *  pressure employing functions for internal energy and
	 *  density. Can be used if there are no coefficients
	 *  available. */
	double h_u(double T_K, double p) {
		return u(T_K,p)+p/rho(T_K,p);
	};

	/// Internal energy from h,p and rho.
	/** Calculate internal energy as a function of temperature
	 *  and pressure employing functions for enthalpy and
	 *  density. Can be used if there are no coefficients
	 *  available. */
	double u_h(double T_K, double p) {
		return h(T_K,p)-p/rho(T_K,p);
	};
	
	/*
	This function implements a 1-D bounded solver using the algorithm from Brent, R. P.,
	Algorithms for Minimization Without Derivatives.
	Englewood Cliffs, NJ: Prentice-Hall, 1973. Ch. 3-4.

	a and b must bound the solution of interest and f(a) and f(b) must have opposite signs.
	If the function is continuous, there must be at least one solution in the interval [a,b].

	@param f A pointer to an instance of the FuncWrapper1D class that must implement the class() function
	@param a The minimum bound for the solution of f=0
	@param b The maximum bound for the solution of f=0
	@param macheps The machine precision
	@param t Tolerance (absolute)
	@param maxiter Maximum numer of steps allowed.  Will throw a SolutionError if the solution cannot be found
	@param errstr A pointer to the error string returned.  If length is zero, no errors found.
	*/
//double Brent(FuncWrapper1D *f, double a, double b, double macheps, double t, int maxiter, std::string *errstr)







	/*
	 * Some more functions to provide a single implementation
	 * of important routines.
	 * We start with the check functions that can validate input
	 * in terms of pressure p and temperature T.
	 */

	/// Check validity of temperature input.
	/** Compares the given temperature T to a stored minimum and
	 *  maximum temperature. Enforces the redefinition of Tmin and
	 *  Tmax since the default values cause and error. */
	bool checkT(double T){
		if( Tmin < 0. ) {
			throw ValueError("Please specify the minimum temperature.");
		} else if( Tmax < 0.) {
			throw ValueError("Please specify the maximum temperature.");
		} else if ( (Tmin>T) or (T>Tmax) ) {
			throw ValueError(format("Your temperature %f is not between %f and %f.",T,Tmin,Tmax));
		} else {
			return true;
		}
		return false;
	}

	/// Check validity of pressure input.
	/** Compares the given pressure p to the saturation pressure at
	 *  temperature T and throws and exception if p is lower than
	 *  the saturation conditions. */
	bool checkP(double T, double p) {
		double ps = psat(T);
		if (p<ps) {
			throw ValueError(format("Equations are valid for liquid phase only: %f < %f. ",p,ps));
		} else {
			return true;
		}
	}

	/// Check validity of temperature and pressure input.
	bool checkTP(double T, double p) {
		return (checkT(T) and checkP(T,p));
	}

	/// Polynomial function generator.
	/** Base function to produce n-th order polynomials
	 *  based on the length of the coefficients vector.
	 *  Starts with only the first coefficient at x^0. */
	double basePolynomial(std::vector<double> coefficients, double x){
	    double result = 0.;
	    for(unsigned int i=0; i<coefficients.size();i++) {
	    	result += coefficients[i] * pow(x,i);
	    }
	    return result;
	}

	/// Exponential function generator.
	/** Base function to produce exponential functions
	 *  based on the length of the coefficients vector.
	 *  An extra check is performed to compare the length
	 *  of the vector to the input num. */
	double baseExponential(std::vector<double> coefficients, double x, int num = 3){
		double result = 0.;
	    if (coefficients.size()==(unsigned int) num){
	    	if (num==3) {
	    		result = exp(coefficients[0]/(x+coefficients[1]) - coefficients[2]);
	    	} else {
	        	throw ValueError(format("There is no function defined for this number of coefficients (%d). ",coefficients.size()));
	        }
	    } else {
	    	throw ValueError(format("The number of coefficients %d does not match with %d. ",coefficients.size(),num));
	    }
	    return result;
	}
};

bool IsIncompressibleLiquid(std::string name);
double IncompLiquid(long iOutput, double T, double p, long iFluid);
double IncompLiquid(long iOutput, double T, double p, std::string name);


class ExampleClass : public IncompressibleLiquid{
public:
	ExampleClass(){
        name = std::string("ABC");
        description = std::string("ABC mixture");
        reference = std::string("BibTex Key");
    };

    double rho(double T_K, double p){
        return -0.731182*T_K+1076.5;
    }
    double cp(double T_K, double p){
        return 0.00287576*T_K+0.999729;
    }
    double u(double T_K, double p){
        return 0.00287576*(T_K*T_K-298*298)/2.0+0.999729*(T_K-298);
    }
    double s(double T_K, double p){
        return 0.00287576*(T_K-298)/2.0+0.999729*log(T_K/298);
    }
    double visc(double T_K, double p){
        return exp(7.03331e-05*T_K*T_K+-0.0566396*T_K+3.5503);
    }
    double cond(double T_K, double p){
        return -2.06364e-07*T_K+0.000189132;
    }
    double h(double T_K, double p){
    	return h_u(T_K,p);
    }
};

class TherminolD12 : public IncompressibleLiquid{
public:
	TherminolD12(){
        name = std::string("TD12");
        description = std::string("Therminol D12");
        reference = std::string("BibTex Key");
    };

    double rho(double T_K, double p){
    	std::vector<double> cDensity;
    	cDensity.push_back( 776.257);
    	cDensity.push_back(-0.696982);
    	cDensity.push_back(-0.000131384);
    	cDensity.push_back(-0.00000209079);
    	return basePolynomial(cDensity, (T_K-273.15));
    }
    double cp(double T_K, double p){
        return 0.00287576*T_K+0.999729;
    }
    double u(double T_K, double p){
        return 0.00287576*(T_K*T_K-298*298)/2.0+0.999729*(T_K-298);
    }
    double s(double T_K, double p){
        return 0.00287576*(T_K-298)/2.0+0.999729*log(T_K/298);
    }
    double visc(double T_K, double p){
        return exp(7.03331e-05*T_K*T_K+-0.0566396*T_K+3.5503);
    }
    double cond(double T_K, double p){
        return -2.06364e-07*T_K+0.000189132;
    }
    double h(double T_K, double p){
    	return h_u(T_K,p);
    }
};


#endif
