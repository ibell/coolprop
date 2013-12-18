
#ifndef INCOMPRESSIBLE_BASE_H
#define INCOMPRESSIBLE_BASE_H

#include <string>
#include <vector>
#include <math.h>
#include "FluidClass.h"
#include "GlobalConstants.h"
#include "CPExceptions.h"
#include "CoolPropTools.h"
#include "MatrixMath.h"


/// The base class for incompressible fluids only
class IncompressibleClass{
	
protected:
	std::string name;
	std::string description;
	std::string reference;
	double Tmin;
	double TminPsat;
	double Tmax;
	double Tref;

	bool DEBUG;

public:
	std::string getDescription() const {return description;}
	std::string getName() const {return name;}
	// For backwards-compatibility.
	std::string get_name() const {return getName();}
	std::string getReference() const {return reference;}
	double getTmax() const {return Tmax;}
	double getTmin() const {return Tmin;}
	double getTminPsat() const {return TminPsat;}
	double getTref() const {return Tref;}

	void setDebug(bool debug) {this->DEBUG = debug;}

public:
	// Constructor
	IncompressibleClass(){
		name = "";
		description = "";
		reference = "";
		Tmin = -1.;
		Tmax = -1.;
		TminPsat = -1.;
		Tref = 273.15 + 25. ;
		DEBUG = false;
	};

	// Destructor.  No implementation
	virtual ~IncompressibleClass(){};


protected:
	/// Basic checks for coefficient vectors.
	/** Starts with only the first coefficient dimension
	 *  and checks the vector length against parameter n. */
	bool checkCoefficients(std::vector<double> const& coefficients, unsigned int n);
	bool checkCoefficients(std::vector< std::vector<double> > const& coefficients, unsigned int rows, unsigned int columns);

private:
	/** The core of the polynomial wrappers are the different
	 *  implementations that follow below. In case there are
	 *  new calculation schemes available, please do not delete
	 *  the implementations, but mark them as deprecated.
	 *  The old functions are good for debugging since the
	 *  structure is easier to read than the backward Horner-scheme
	 *  or the recursive Horner-scheme.
	 */

	/// Simple polynomial function generator. <- Deprecated due to poor performance, use Horner-scheme instead
	/** Base function to produce n-th order polynomials
	 *  based on the length of the coefficient vector.
	 *  Starts with only the first coefficient at T^0. */
	DEPRECATED(double simplePolynomial(std::vector<double> const& coefficients, double T));
	DEPRECATED(double simplePolynomial(std::vector<std::vector<double> > const& coefficients, double x, double T));

	/// Simple integrated polynomial function generator.
	/** Base function to produce integrals of n-th order polynomials based on
	 *  the length of the coefficient vector.
	 *  Starts with only the first coefficient at T^0 */
	///Indefinite integral in T-direction
	double simplePolynomialInt(std::vector<double> const& coefficients, double T);
	///Definite integral from T0 to T1
	double simplePolynomialInt(std::vector<double> const& coefficients, double T1, double T0);
	///Indefinite integral in T-direction only
	double simplePolynomialInt(std::vector<std::vector<double> > const& coefficients, double x, double T);
	///Definite integral from T0 to T1
	double simplePolynomialInt(std::vector<std::vector<double> > const& coefficients, double x, double T1, double T0);

	/// Simple integrated polynomial function generator divided by independent variable.
	/** Base function to produce integrals of n-th order
	 *  polynomials based on the length of the coefficient
	 *  vector. Starts with only the first coefficient at T^0 */
	///Indefinite integral of a polynomial divided by its independent variable
	double simpleFracInt(std::vector<double> const& coefficients, double T);
	///Definite integral from T0 to T1 of a polynomial divided by its independent variable
	double simpleFracInt(std::vector<double> const& coefficients, double T1, double T0);
	///Indefinite integral of a polynomial divided by its 2nd independent variable
	double simpleFracInt(std::vector<std::vector<double> > const& coefficients, double x, double T);
	///Definite integral from T0 to T1 of a polynomial divided by its 2nd independent variable
	double simpleFracInt(std::vector<std::vector<double> > const& coefficients, double x, double T1, double T0);


	/** Simple integrated centred(!) polynomial function generator divided by independent variable.
	 *  We need to rewrite some of the functions in order to
	 *  use central fit. Having a central temperature Tbase
	 *  allows for a better fit, but requires a different
	 *  formulation of the fracInt function group. Other
	 *  functions are not affected.
	 *  Starts with only the first coefficient at T^0 */
	///Helper function to calculate the D vector:
	double factorial(double nValue);
	double binom(double nValue, double nValue2);
	std::vector<double> fracIntCentralDvector(int m, double T, double Tbase);
	std::vector<double> fracIntCentralDvector(int m, double T1, double T0, double Tbase);
	///Indefinite integral of a centred polynomial divided by its independent variable
	double fracIntCentral(std::vector<double> const& coefficients, double T, double Tbase);
	///Definite integral from T0 to T1 of a centred polynomial divided by its independent variable
	double fracIntCentral(std::vector<double> const& coefficients, double T1, double T0, double Tbase);


	/// Horner function generator implementations
	/** Represent polynomials according to Horner's scheme.
	 *  This avoids unnecessary multiplication and thus
	 *  speeds up calculation.
	 */
	double baseHorner(std::vector<double> const& coefficients, double T);
	double baseHorner(std::vector< std::vector<double> > const& coefficients, double x, double T);
	///Indefinite integral in T-direction
	double baseHornerInt(std::vector<double> const& coefficients, double T);
	///Definite integral from T0 to T1
	double baseHornerInt(std::vector<double> const& coefficients, double T1, double T0);
	///Indefinite integral in T-direction only
	double baseHornerInt(std::vector<std::vector<double> > const& coefficients, double x, double T);
	///Definite integral from T0 to T1
	double baseHornerInt(std::vector<std::vector<double> > const& coefficients, double x, double T1, double T0);
	///Indefinite integral of a polynomial divided by its independent variable
	double baseHornerFra(std::vector<double> const& coefficients, double T);
	///Definite integral from T0 to T1 of a polynomial divided by its independent variable
	double baseHornerFra(std::vector<double> const& coefficients, double T1, double T0);
	///Indefinite integral of a polynomial divided by its 2nd independent variable
	double baseHornerFra(std::vector<std::vector<double> > const& coefficients, double x, double T);
	///Definite integral from T0 to T1 of a polynomial divided by its 2nd independent variable
	double baseHornerFra(std::vector<std::vector<double> > const& coefficients, double x, double T1, double T0);


	/** Integrating coefficients for polynomials is done by dividing the
	 *  original coefficients by (i+1) and elevating the order by 1
	 *  through adding a zero as first coefficient.
	 *  Some reslicing needs to be applied to integrate along the x-axis.
	 *  In the brine/solution equations, reordering of the parameters
	 *  avoids this expensive operation. However, it is included for the
	 *  sake of completeness.
	 */
	std::vector<double> integrateCoeffs(std::vector<double> const& coefficients);
	std::vector< std::vector<double> > integrateCoeffs(std::vector< std::vector<double> > const& coefficients, bool axis);

	/** Deriving coefficients for polynomials is done by multiplying the
	 *  original coefficients with i and lowering the order by 1.
	 *
	 *  It is not really deprecated, but untested and therefore a warning
	 *  is issued. Please check this method before you use it.
	 */
	DEPRECATED(std::vector<double> deriveCoeffs(std::vector<double> const& coefficients));
	DEPRECATED(std::vector< std::vector<double> > deriveCoeffs(std::vector< std::vector<double> > const& coefficients, unsigned int axis));

	/** Alternatives
	 *  Simple functions that heavily rely on other parts of this file.
	 *  We still need to check which combinations yield the best
	 *  performance.
	 */
	///Indefinite integral in T-direction
	double integrateIn2Steps(std::vector<double> const& coefficients, double T);
	///Definite integral from T0 to T1
	double integrateIn2Steps(std::vector<double> const& coefficients, double T1, double T0);
	///Indefinite integral in terms of x(axis=true) or T(axis=false).
	double integrateIn2Steps(std::vector< std::vector<double> > const& coefficients, double x, double T, bool axis);
	///Definite integral from T0 to T1
	double integrateIn2Steps(std::vector< std::vector<double> > const& coefficients, double x, double T1, double T0);
	///Indefinite integral in T-direction of a polynomial divided by its independent variable
	double fracIntIn2Steps(std::vector<double> const& coefficients, double T);
	///Definite integral from T0 to T1 of a polynomial divided by its independent variable
	double fracIntIn2Steps(std::vector<double> const& coefficients, double T1, double T0);
	///Indefinite integral in T-direction of a polynomial divided by its 2nd independent variable
	double fracIntIn2Steps(std::vector<std::vector<double> > const& coefficients, double x, double T);
	///Definite integral from T0 to T1 of a polynomial divided by its 2nd independent variable
	double fracIntIn2Steps(std::vector<std::vector<double> > const& coefficients, double x, double T1, double T0);
	///Indefinite integral of a centred polynomial divided by its 2nd independent variable
	double fracIntCentral2Steps(std::vector<std::vector<double> > const& coefficients, double x, double T, double Tbase);
	///Definite integral from T0 to T1 of a centred polynomial divided by its 2nd independent variable
	double fracIntCentral2Steps(std::vector<std::vector<double> > const& coefficients, double x, double T1, double T0, double Tbase);

public:
	/** Here we define the functions that should be used by the
	 *  respective implementations. Please do no use any other
	 *  method since this would break the purpose of this interface.
	 *  Note that the functions below are supposed to be aliases
	 *  to implementations declared elsewhere in this file.
	 */

	/// Evaluates a one-dimensional polynomial for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input
	double polyval(std::vector<double> const& coefficients, double x){
		//return simplePolynomial(coefficients,x);
		return baseHorner(coefficients,x);
	}

	/// Evaluates a two-dimensional polynomial for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y double value that represents the current input in the 2nd dimension
	double polyval(std::vector< std::vector<double> > const& coefficients, double x, double y){
		//return simplePolynomial(coefficients,x,y);
		return baseHorner(coefficients,x,y);
	}

	/// Evaluates the indefinite integral of a one-dimensional polynomial
	/// @param coefficients vector containing the ordered coefficients
	/// @param T double value that represents the current input
	double polyint(std::vector<double> const& coefficients, double T){
		//return simplePolynomialInt(coefficients,T);
		return baseHornerInt(coefficients,T);
		//return integrateIn2Steps(coefficients,T);
	}

	/// Evaluates the definite integral of a one-dimensional polynomial
	/// @param coefficients vector containing the ordered coefficients
	/// @param T1 double value that represents the current position
	/// @param T0 double value that represents the reference state
	double polyint(std::vector<double> const& coefficients, double T1, double T0){
		//return simplePolynomialInt(coefficients,T1,T0);
		return baseHornerInt(coefficients,T1,T0);
		//return integrateIn2Steps(coefficients,T1,T0);
	}

	/// Evaluates the indefinite integral of a two-dimensional polynomial along the 2nd axis (T)
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param T double value that represents the current input in the 2nd dimension
	double polyint(std::vector< std::vector<double> > const& coefficients, double x, double T){
		//return simplePolynomialInt(coefficients,x,T);
		return baseHornerInt(coefficients,x,T);
		//return integrateIn2Steps(coefficients,x,T,false);
	}

	/// Evaluates the definite integral of a two-dimensional polynomial along the 2nd axis (T)
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param T1 double value that represents the current input in the 2nd dimension
	/// @param T0 double value that represents the reference state in the 2nd dimension
	double polyint(std::vector< std::vector<double> > const& coefficients, double x, double T1, double T0){
		//return simplePolynomialInt(coefficients,x,T1,T0);
		return baseHornerInt(coefficients,x,T1,T0);
		//return integrateIn2Steps(coefficients,x,T1,T0);
	}

	/// Evaluates the indefinite integral of a one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param T double value that represents the current position
	double polyfracint(std::vector<double> const& coefficients, double T){
		//return simpleFracInt(coefficients,T);
		return baseHornerFra(coefficients,T);
		//return fracIntIn2Steps(coefficients,T);
	}

	/// Evaluates the definite integral of a one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param T1 double value that represents the current position
	/// @param T0 double value that represents the reference state
	double polyfracint(std::vector<double> const& coefficients, double T1, double T0){
		//return simpleFracInt(coefficients,T1,T0);
		//return baseHornerFracInt(coefficients,T1,T0);
		return fracIntIn2Steps(coefficients,T1,T0);
	}

	/// Evaluates the indefinite integral of a two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param T double value that represents the current input in the 2nd dimension
	double polyfracint(std::vector< std::vector<double> > const& coefficients, double x, double T){
		//return simpleFracInt(coefficients,x,T);
		return baseHornerFra(coefficients,x,T);
		//return fracIntIn2Steps(coefficients,x,T);
	}

	/// Evaluates the definite integral of a two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param T1 double value that represents the current input in the 2nd dimension
	/// @param T0 double value that represents the reference state in the 2nd dimension
	double polyfracint(std::vector< std::vector<double> > const& coefficients, double x, double T1, double T0){ //TODO compare speed
		//return simpleFracInt(coefficients,x,T1,T0);
		return baseHornerFra(coefficients,x,T1,T0);
		//return fracIntIn2Steps(coefficients,x,T1,T0);
	}

	/// Evaluates the indefinite integral of a centred one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param T double value that represents the current position
	/// @param Tbase central temperature for fitted function
	double polyfracintcentral(std::vector<double> const& coefficients, double T, double Tbase){
		return fracIntCentral(coefficients,T,Tbase);
	}

	/// Evaluates the definite integral of a centred one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param T1 double value that represents the current position
	/// @param T0 double value that represents the reference state
	/// @param Tbase central temperature for fitted function
	double polyfracintcentral(std::vector<double> const& coefficients, double T1, double T0, double Tbase){
		return fracIntCentral(coefficients,T1,T0,Tbase);
	}

	/// Evaluates the indefinite integral of a centred two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param T double value that represents the current input in the 2nd dimension
	/// @param Tbase central temperature for fitted function
	double polyfracintcentral(std::vector< std::vector<double> > const& coefficients, double x, double T, double Tbase){
		return fracIntCentral2Steps(coefficients,x,T,Tbase);
	}

	/// Evaluates the definite integral of a centred two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param T1 double value that represents the current input in the 2nd dimension
	/// @param T0 double value that represents the reference state in the 2nd dimension
	/// @param Tbase central temperature for fitted function
	double polyfracintcentral(std::vector< std::vector<double> > const& coefficients, double x, double T1, double T0, double Tbase){
		return fracIntCentral2Steps(coefficients,x,T1,T0,Tbase);
	}



	/// Evaluates an exponential function for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param T double value that represents the current input
	/// @param n int value that determines the kind of exponential function
	double expval(std::vector<double> const& coefficients, double T, int n);

	/// Evaluates an exponential function for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param T double value that represents the current input in the 2nd dimension
	/// @param n int value that determines the kind of exponential function
	double expval(std::vector< std::vector<double> > const& coefficients, double x, double T, int n);
};


/** Multiple inheritance could be useful to merge the base classes for
 *  incompressible fluids and the normal fluids. However, It seems like
 *  a lot of work to find ways to redefine all the fluid functions while
 *  avoiding to break the functionality of other objects like the fluid
 *  lists etc...
 */
class IncompressibleFluid : public Fluid, public IncompressibleClass {
public:
	IncompressibleFluid(){};
    ~IncompressibleFluid(){};
    // Some functions need top be overwritten!
};


#endif
