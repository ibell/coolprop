
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
	 *  Starts with only the first coefficient at x^0. */
	DEPRECATED(double simplePolynomial(std::vector<double> const& coefficients, double x));

	/// Simple integrated polynomial function generator. <- Deprecated due to poor performance, use Horner-scheme instead
	/** Base function to produce integrals of n-th order
	 *  polynomials based on the length of the coefficient
	 *  vector. Integrates from x0 to x1.
	 *  Starts with only the first coefficient at x^0 */
	DEPRECATED(double simplePolynomialInt(std::vector<double> const& coefficients, double x1, double x0));


	/// Horner function generator implementations
	/** Represent polynomials according to Horner's scheme.
	 *  This avoids unnecessary multiplication and thus
	 *  speeds up calculation.
	 */
	double baseHorner(std::vector<double> const& coefficients, double x);
	double baseHorner(std::vector< std::vector<double> > const& coefficients, double x, double y);
	double baseHornerIntegrated(std::vector<double> const& coefficients, double x);
	double baseHornerIntegrated(std::vector<double> const& coefficients, double x1, double x0);
	double baseHornerIntegrated(std::vector< std::vector<double> > const& coefficients, double x, double y, unsigned int axis);
	double baseHornerIntegrated(std::vector< std::vector<double> > const& coefficients, double x, double y1, double y0);


	/** Integrating coefficients for polynomials is done by dividing the
	 *  original coefficients by (i+1) and elevating the order by 1.
	 *  Some reslicing needs to be applied to integrate along the x-axis.
	 *  In the brine/solution equations, reordering of the parameters
	 *  avoids this expensive operation. However, it is included for the
	 *  sake of completeness.
	 */
	std::vector<double> integrateCoeffs(std::vector<double> const& coefficients);
	std::vector< std::vector<double> > integrateCoeffs(std::vector< std::vector<double> > const& coefficients, unsigned int axis);

	/** Deriving coefficients for polynomials is done by multiplying the
	 *  original coefficients with i and lowering the order by 1.
	 *
	 *  It is not really deprecated, but untested and therefore a warning
	 *  is issued. Please check this method before you use it.
	 */
	DEPRECATED(std::vector<double> deriveCoeffs(std::vector<double> const& coefficients));
	DEPRECATED(std::vector< std::vector<double> > deriveCoeffs(std::vector< std::vector<double> > const& coefficients, unsigned int axis));


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
	    return baseHorner(coefficients,x);
	}

	/// Evaluates a two-dimensional polynomial for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y double value that represents the current input in the 2nd dimension
	double polyval(std::vector< std::vector<double> > const& coefficients, double x, double y){
		return baseHorner(coefficients,x,y);
	}

	/// Evaluates the indefinite integral of a one-dimensional polynomial
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input
	double polyint(std::vector<double> const& coefficients, double x){
		return baseHornerIntegrated(coefficients,x);
	}

	/// Evaluates the indefinite integral of a two-dimensional polynomial along the 2nd axis (y)
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y double value that represents the current input in the 2nd dimension
	double polyint(std::vector< std::vector<double> > const& coefficients, double x, double y){
		return baseHornerIntegrated(coefficients,x,y,(unsigned int)1);
	}

	/// Evaluates the definite integral of a one-dimensional polynomial
	/// @param coefficients vector containing the ordered coefficients
	/// @param x1 double value that represents the current position
	/// @param x0 double value that represents the reference state
	double polyint(std::vector<double> const& coefficients, double x1, double x0){
		return baseHornerIntegrated(coefficients,x1,x0);
	}

	/// Evaluates the definite integral of a two-dimensional polynomial along the 2nd axis (y)
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y1 double value that represents the current input in the 2nd dimension
	/// @param y0 double value that represents the reference state in the 2nd dimension
	double polyint(std::vector< std::vector<double> > const& coefficients, double x, double y1, double y0){
		return baseHornerIntegrated(coefficients,x,y1,y0);
	}

	/// Evaluates the definite integral of a one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x1 double value that represents the current position
	/// @param x0 double value that represents the reference state
	double fracint(std::vector<double> const& coefficients, double x1, double x0);

	/// Evaluates the definite integral of a two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y1 double value that represents the current input in the 2nd dimension
	/// @param y0 double value that represents the reference state in the 2nd dimension
	double fracint(std::vector< std::vector<double> > const& coefficients, double x, double y1, double y0);

	/// Evaluates an exponential function for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input
	/// @param n int value that determines the kind of exponential function
	double expval(std::vector<double> const& coefficients, double x, int n);

	/// Evaluates an exponential function for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y double value that represents the current input in the 2nd dimension
	/// @param n int value that determines the kind of exponential function
	double expval(std::vector< std::vector<double> > const& coefficients, double x, double y, int n);
};




/** Multiple inheritance could be useful to merge the base classes for
 *  incompressible fluids and the normal fluids. However, It seems like
 *  a lot of work to find ways to redefine all the fluid functions while
 *  avoiding to break the functionality of other objects like the fluid
 *  lists etc...
 */
class IncompressibleFluid : public Fluid, public IncompressibleClass {
public:
	IncompressibleFluid();
    ~IncompressibleFluid(){};
    // Some functions need top be overwritten!
};


#endif
