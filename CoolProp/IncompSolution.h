
#ifndef INCOMPRESSIBLE_SOLUTION_H
#define INCOMPRESSIBLE_SOLUTION_H

#include <string>
#include <vector>
#include "CPExceptions.h"
#include "CoolProp.h"
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
class IncompressibleSolution{

protected:
	std::string name,description,reference;
	double Tmin, TminPsat, Tmax, Tref, xmin, xmax;

	/// Function used to enforce the composition as parameter
	/** Overwrites the normal functions that only take 2
	 *  parameters (T,p) and throws an exception. */
	bool needComposition() {
		throw NotImplementedError(format("The fluid %s needs an additional input for the composition.",this->name));
		return false;
	}

public:
	// Constructor
	IncompressibleSolution(){
		Tmin = -1.;
		Tmax = -1.;
		TminPsat = -1.;
		Tref = 273.15 + 25. ;
		xmin = -1.;
		xmax = -1.;
	};

	// Destructor, no implementation
	virtual ~IncompressibleSolution(){};

	/* All functions need T, p and x as input. Might not
	 * be necessary, but gives a clearer structure.
	 */
    virtual double rho (double T_K, double p, double x);
    virtual double cp  (double T_K, double p, double x);
    virtual double h   (double T_K, double p, double x);
    virtual double s   (double T_K, double p, double x);
    virtual double visc(double T_K, double p, double x);
    virtual double cond(double T_K, double p, double x);
    virtual double u   (double T_K, double p, double x);
    virtual double psat(double T_K          , double x);
    virtual double Tsat(            double p, double x);
    virtual double Tfreeze(         double p, double x);

protected:
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
	 *  defined for all fluids, default values do not
	 *  cause errors. */
	bool checkT(double T_K, double p, double x){
		if( Tmin < 0. ) {
			throw ValueError("Please specify the minimum temperature.");
		} else if( Tmax < 0.) {
			throw ValueError("Please specify the maximum temperature.");
		} else if ( (Tmin>T_K) || (T_K>Tmax) ) {
			throw ValueError(format("Your temperature %f is not between %f and %f.",T_K,Tmin,Tmax));
		} else if (T_K < Tfreeze(p,x)) {
			throw ValueError("Your temperature %f is below the freezing point of %f.",T_K,Tfreeze(p,x));
		} else {
			return true;
		}
		return false;
	}

	/// Check validity of pressure input.
	/** Compares the given pressure p to the saturation pressure at
	 *  temperature T and throws and exception if p is lower than
	 *  the saturation conditions.
	 *  The default value for psat is -1 yielding true if psat
	 *  is not redefined in the subclass.
	 *  */
	bool checkP(double T_K, double p, double x) {
		double ps = psat(T_K,x);
		if (p<ps) {
			throw ValueError(format("Equations are valid for liquid phase only: %f < %f. ",p,ps));
		} else {
			return true;
		}
	}

	/// Check validity of composition input.
	/** Compares the given composition x to a stored minimum and
	 *  maximum value. Enforces the redefinition of xmin and
	 *  xmax since the default values cause an error. */
	bool checkX(double x){
		if( xmin < 0. ) {
			throw ValueError("Please specify the minimum concentration.");
		} else if( xmax < 0.) {
			throw ValueError("Please specify the maximum concentration.");
		} else if ( (xmin>x) || (x>xmax) ) {
			throw ValueError(format("Your composition %f is not between %f and %f.",x,xmin,xmax));
		} else {
			return true;
		}
		return false;
	}

	/// Check validity of temperature, pressure and composition input.
	bool checkTP(double T_K, double p){needComposition();return false;}
	bool checkTPX(double T, double p, double x) {
		return (checkT(T,p,x) && checkP(T,p,x) && checkX(x));
	}
};

bool IsIncompressibleSolution(std::string name);
double IncompSolution(long iOutput, double T, double p, double x, long iFluid);
double IncompSolution(long iOutput, double T, double p, double x, std::string name);

/// Class to use the SecCool parameters
/** Employs some basic wrapper-like functionality
 *  to bridge the gap between the solution functions
 *  used in CoolProp and the definition used in
 *  SecCool. Please visit:
 *  http://en.ipu.dk/Indhold/refrigeration-and-energy-technology/seccool.aspx
 *  Many thanks to Morten Juel Skovrup for providing
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
