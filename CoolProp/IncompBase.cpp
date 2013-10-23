
#include <string>
#include <vector>
#include <math.h>
#include "CoolPropTools.h"
#include "CPExceptions.h"
#include "IncompBase.h"
#include "IncompLiquid.h"
#include "IncompSolution.h"

#include <stdio.h>



/// Basic checks for coefficient vectors.
/** Starts with only the first coefficient dimension
 *  and checks the vector length against parameter n. */
bool IncompressibleClass::checkCoefficients(std::vector<double> coefficients, unsigned int n){
	if (coefficients.size() == n){
		return true;
	} else {
		throw ValueError(format("The number of coefficients %d does not match with %d. ",coefficients.size(),n));
	}
	return false;
}

bool IncompressibleClass::checkCoefficients(std::vector< std::vector<double> > coefficients, unsigned int rows, unsigned int columns){
	if (coefficients.size() == rows){
		bool result = true;
		for(unsigned int i=0; i<rows; i++) {
			result = result && checkCoefficients(coefficients[i],columns);
		}
		return result;
	} else {
		throw ValueError(format("The number of rows %d does not match with %d. ",coefficients.size(),rows));
	}
	return false;
}


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
double IncompressibleClass::simplePolynomial(std::vector<double> coefficients, double x){
	double result = 0.;
	for(unsigned int i=0; i<coefficients.size();i++) {
		result += coefficients[i] * pow(x,(int)i);
	}
	return result;
}

/// Simple integrated polynomial function generator. <- Deprecated due to poor performance, use Horner-scheme instead
/** Base function to produce integrals of n-th order
 *  polynomials based on the length of the coefficient
 *  vector. Integrates from x0 to x1.
 *  Starts with only the first coefficient at x^0 */
double IncompressibleClass::simplePolynomialInt(std::vector<double> coefficients, double x1, double x0){
	double result = 0.;
	for(unsigned int i=0; i<coefficients.size();i++) {
		result += 1./(i+1.) * coefficients[i] * (pow(x1,(i+1.)) - pow(x0,(i+1.)));
	}
	return result;
}


/// Horner function generator implementations
/** Represent polynomials according to Horner's scheme.
 *  This avoids unnecessary multiplication and thus
 *  speeds up calculation.
 */
double IncompressibleClass::baseHorner(std::vector<double> coefficients, double x){
	double result = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result = result * x + coefficients[i];
	}
	return result;
}
double IncompressibleClass::baseHorner(std::vector< std::vector<double> > coefficients, double x, double y){
	double result = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result = result * x + baseHorner(coefficients[i], y);
	}
	return result;
}
double IncompressibleClass::baseHornerIntegrated(std::vector<double> coefficients, double x){
	return baseHorner(integrateCoeffs(coefficients),x);
}
double IncompressibleClass::baseHornerIntegrated(std::vector<double> coefficients, double x1, double x0){
	std::vector<double> coefficientsInt = integrateCoeffs(coefficients);
	return baseHorner(coefficientsInt,x1)-baseHorner(coefficientsInt,x0);
}
double IncompressibleClass::baseHornerIntegrated(std::vector< std::vector<double> > coefficients, double x, double y, unsigned int axis){
	return baseHorner(integrateCoeffs(coefficients,axis),x,y);
}
double IncompressibleClass::baseHornerIntegrated(std::vector< std::vector<double> > coefficients, double x, double y1, double y0){
	std::vector< std::vector<double> > coefficientsInt = integrateCoeffs(coefficients,1);
	return baseHorner(coefficientsInt,x,y1)-baseHorner(coefficientsInt,x,y0);
}


/** Integrating coefficients for polynomials is done by dividing the
 *  original coefficients by (i+1) and elevating the order by 1.
 *  Some reslicing needs to be applied to integrate along the x-axis.
 *  In the brine/solution equations, reordering of the parameters
 *  avoids this expensive operation. However, it is included for the
 *  sake of completeness.
 */
std::vector<double> IncompressibleClass::integrateCoeffs(std::vector<double> coefficients){
	std::vector<double> newCoefficients;
	unsigned int sizeX = coefficients.size();
	if (sizeX<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",sizeX));
	// pushing a zero elevates the order by 1
	newCoefficients.push_back(0.0);
	for(unsigned int i=0; i<coefficients.size(); i++) {
		newCoefficients.push_back(coefficients[i]/(i+1.));
	}
	return newCoefficients;
}

std::vector< std::vector<double> > IncompressibleClass::integrateCoeffs(std::vector< std::vector<double> > coefficients, unsigned int axis){
	std::vector< std::vector<double> > newCoefficients;
	unsigned int sizeX = coefficients.size();
	if (sizeX<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",sizeX));

	if (axis==0){
		std::vector< std::vector<double> > tmpCoefficients;
		tmpCoefficients = transpose(coefficients);
		unsigned int sizeY = tmpCoefficients.size();
		for(unsigned int i=0; i<sizeY; i++) {
			newCoefficients.push_back(integrateCoeffs(tmpCoefficients[i]));
		}
		return transpose(newCoefficients);
	} else if (axis==1){
		for(unsigned int i=0; i<sizeX; i++) {
			newCoefficients.push_back(integrateCoeffs(coefficients[i]));
		}
		return newCoefficients;
	} else {
		throw ValueError(format("You can only use x-axis (0) and y-axis (1) for integration. %d is not a valid input. ",axis));
	}
	return newCoefficients;
}

/** Deriving coefficients for polynomials is done by multiplying the
 *  original coefficients with i and lowering the order by 1.
 */
std::vector<double> IncompressibleClass::deriveCoeffs(std::vector<double> coefficients){
	std::vector<double> newCoefficients;
	unsigned int sizeX = coefficients.size();
	if (sizeX<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",sizeX));
	// skipping the first element lowers the order
	for(unsigned int i=1; i<coefficients.size(); i++) {
		newCoefficients.push_back(coefficients[i]*i);
	}
	return newCoefficients;
}

//std::vector< std::vector<double> > IncompressibleClass::deriveCoeffs(std::vector< std::vector<double> > coefficients, unsigned int axis){
//	std::vector< std::vector<double> > newCoefficients;
//	unsigned int sizeX = coefficients.size();
//	if (sizeX<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",sizeX));
//
//	if (axis==0){
//		std::vector< std::vector<double> > tmpCoefficients;
//		tmpCoefficients = transpose(coefficients);
//		unsigned int sizeY = tmpCoefficients.size();
//		for(unsigned int i=0; i<sizeY; i++) {
//			newCoefficients.push_back(deriveCoeffs(tmpCoefficients[i]));
//		}
//		return transpose(newCoefficients);
//	} else if (axis==1){
//		for(unsigned int i=0; i<sizeX; i++) {
//			newCoefficients.push_back(deriveCoeffs(coefficients[i]));
//		}
//		return newCoefficients;
//	} else {
//		throw ValueError(format("You can only use x-axis (0) and y-axis (1) for derivation. %d is not a valid input. ",axis));
//	}
//	return newCoefficients;
//}


/** Here we define the functions that should be used by the
 *  respective implementations. Please do no use any other
 *  method since this would break the purpose of this interface.
 */

/// Evaluates the definite integral of a one-dimensional polynomial divided by its independent variable
/// @param coefficients vector containing the ordered coefficients
/// @param x1 double value that represents the current position
/// @param x0 double value that represents the reference state
double IncompressibleClass::fracint(std::vector<double> coefficients, double x1, double x0){
	double result = coefficients[0] * log(x1/x0);
	if (coefficients.size() > 1) {
		std::vector<double> newCoeffs(coefficients.begin() + 1, coefficients.end());
		result += polyint(newCoeffs,x1,x0);
	}
	return result;
}

/// Evaluates the definite integral of a two-dimensional polynomial divided by its 2nd independent variable
/// @param coefficients vector containing the ordered coefficients
/// @param x double value that represents the current input in the 1st dimension
/// @param y1 double value that represents the current input in the 2nd dimension
/// @param y0 double value that represents the reference state in the 2nd dimension
double IncompressibleClass::fracint(std::vector< std::vector<double> > coefficients, double x, double y1, double y0){
	double row    = 0;
	std::vector<double> newCoeffs;
	double dlog = log(y1/y0);
	for (unsigned int i=0; i<coefficients.size(); i++){
		// process the first entry and remove it from the vector
		row  = coefficients[i][0] * dlog;
		if (coefficients[i].size() > 1) {
			std::vector<double> tmpCoeffs(coefficients[i].begin() + 1, coefficients[i].end());
			row += polyint(tmpCoeffs,y1,y0);
		}
		// save the result in a new vector entry
		newCoeffs.push_back(row);
	}
	// in the end we process the new coefficient vector as polynomial
	return polyval(newCoeffs,x);
}


/// Evaluates an exponential function for the given coefficients
/// @param coefficients vector containing the ordered coefficients
/// @param x double value that represents the current input
/// @param n int value that determines the kind of exponential function
double IncompressibleClass::expval(std::vector<double> coefficients, double x, int n){
	double result = 0.;
	if (n==1) {
		checkCoefficients(coefficients,3);
		result = exp(coefficients[0]/(x+coefficients[1]) - coefficients[2]);
	} else if (n==2) {
		result = exp(polyval(coefficients, x));
	} else {
		throw ValueError(format("There is no function defined for this input (%d). ",n));
	}
	return result;
}

/// Evaluates an exponential function for the given coefficients
/// @param coefficients vector containing the ordered coefficients
/// @param x double value that represents the current input in the 1st dimension
/// @param y double value that represents the current input in the 2nd dimension
/// @param n int value that determines the kind of exponential function
double IncompressibleClass::expval(std::vector< std::vector<double> > coefficients, double x, double y, int n){
	double result = 0.;
	if (n==2) {
		result = exp(polyval(coefficients, x, y));
	} else {
		throw ValueError(format("There is no function defined for this input (%d). ",n));
	}
	return result;
}


//int main() {
//
//	SimpleIncompressible* liquid = new DowthermQClass();
//	double AT      =  150.0 + 273.15;
//	double Ap      =  3e5;
//    liquid->testInputs(AT,Ap);
//
//
//	SecCoolSolution* obj = new MethanolSolution();
//    double x      =   0.25;
//    double T      =   5.0 + 273.15;
//    double p      =   3e5;
//
//	obj->testInputs(T+00,p,x);
//	obj->testInputs(T+05,p,x);
//	obj->testInputs(T+10,p,x);
//	obj->testInputs(T+15,p,x);
//
//
//}
