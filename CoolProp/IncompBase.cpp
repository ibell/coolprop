
#include <string>
#include <vector>
#include <math.h>
#include <stdio.h>
#include "CoolPropTools.h"
#include "CPExceptions.h"
#include "IncompBase.h"
#include "IncompLiquid.h"
#include "IncompSolution.h"


/// Basic checks for coefficient vectors.
/** Starts with only the first coefficient dimension
 *  and checks the vector length against parameter n. */
bool IncompressibleClass::checkCoefficients(std::vector<double> const& coefficients, unsigned int n){
	if (coefficients.size() == n){
		return true;
	} else {
		throw ValueError(format("The number of coefficients %d does not match with %d. ",coefficients.size(),n));
	}
	return false;
}

bool IncompressibleClass::checkCoefficients(std::vector< std::vector<double> > const& coefficients, unsigned int rows, unsigned int columns){
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
 *  Starts with only the first coefficient at T^0. */
double IncompressibleClass::simplePolynomial(std::vector<double> const& coefficients, double T){
	if (this->DEBUG) {
		std::cout << "Running simplePolynomial(std::vector, " << T << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0.;
	for(unsigned int i=0; i<coefficients.size();i++) {
		result += coefficients[i] * pow(T,(int)i);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
double IncompressibleClass::simplePolynomial(std::vector<std::vector<double> > const& coefficients, double x, double T){
	if (this->DEBUG) {
		std::cout << "Running simplePolynomial(std::vector, " << x << ", " << T << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0;
	for(unsigned int i=0; i<coefficients.size();i++) {
		result += pow(x,(int)i) * simplePolynomial(coefficients[i], T);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

/// Simple integrated polynomial function generator.
/** Base function to produce integrals of n-th order
 *  polynomials based on the length of the coefficient
 *  vector.
 *  Starts with only the first coefficient at T^0 */
///Indefinite integral in T-direction
double IncompressibleClass::simplePolynomialInt(std::vector<double> const& coefficients, double T){
	if (this->DEBUG) {
		std::cout << "Running simplePolynomialInt(std::vector, " << T << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0.;
	for(unsigned int i=0; i<coefficients.size();i++) {
		result += 1./(i+1.) * coefficients[i] * pow(T,(int)(i+1.));
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Definite integral from T0 to T1
double IncompressibleClass::simplePolynomialInt(std::vector<double> const& coefficients, double T1, double T0){
	if (this->DEBUG) {
		std::cout << "Running simplePolynomialInt(std::vector, " << T1 << ", " << T0 << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0.;
	for(unsigned int i=0; i<coefficients.size();i++) {
		result += 1./(i+1.) * coefficients[i] * ( pow(T1,(int)(i+1.)) - pow(T0,(int)(i+1.)) );
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Indefinite integral in y-direction only
double IncompressibleClass::simplePolynomialInt(std::vector<std::vector<double> > const& coefficients, double x, double T){
	if (this->DEBUG) {
		std::cout << "Running simplePolynomialInt(std::vector, " << x << ", " << T << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0.;
	for(unsigned int i=0; i<coefficients.size();i++) {
		result += pow(x,(int)i) * simplePolynomialInt(coefficients[i], T);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Definite integral from T0 to T1
double IncompressibleClass::simplePolynomialInt(std::vector<std::vector<double> > const& coefficients, double x, double T1, double T0){
	if (this->DEBUG) {
		std::cout << "Running simplePolynomialInt(std::vector, " << x << ", " << T1 << ", " << T0 << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0.;
	for(unsigned int i=0; i<coefficients.size();i++) {
		result += pow(x,(int)i) * simplePolynomialInt(coefficients[i], T1, T0);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

/// Simple integrated polynomial function generator divided by independent variable.
/** Base function to produce integrals of n-th order
 *  polynomials based on the length of the coefficient
 *  vector.
 *  Starts with only the first coefficient at T^0 */
double IncompressibleClass::simpleFracInt(std::vector<double> const& coefficients, double T){
	if (this->DEBUG) {
		std::cout << "Running      simpleFracInt(std::vector, " << T << "): ";
	}
	double result = coefficients[0] * log(T);
	if (coefficients.size() > 1) {
		for (unsigned int i=1; i<coefficients.size(); i++){
			result += 1./(i) * coefficients[i] * pow(T,(int)(i));
		}
	}
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
double IncompressibleClass::simpleFracInt(std::vector<double> const& coefficients, double T1, double T0){
	if (this->DEBUG) {
		std::cout << "Running      simpleFracInt(std::vector, " << T1 << ", " << T0 << "): ";
	}
	double result = coefficients[0] * log(T1/T0);
	if (coefficients.size() > 1) {
		for (unsigned int i=1; i<coefficients.size(); i++){
			result += 1./(i) * coefficients[i] * (pow(T1,(int)(i))-pow(T0,(int)(i)));
		}
	}
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
double IncompressibleClass::simpleFracInt(std::vector< std::vector<double> > const& coefficients, double x, double T){
	if (this->DEBUG) {
		std::cout << "Running      simpleFracInt(std::vector, " << x << ", " << T << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0;
	for (unsigned int i=0; i<coefficients.size(); i++){
		result += pow(x,(int)i) * polyfracint(coefficients[i],T);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
double IncompressibleClass::simpleFracInt(std::vector< std::vector<double> > const& coefficients, double x, double T1, double T0){
	if (this->DEBUG) {
		std::cout << "Running      simpleFracInt(std::vector, " << x << ", " << T1 << ", " << T0 <<"): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0;
	for (unsigned int i=0; i<coefficients.size(); i++){
		result += pow(x,(int)i) * polyfracint(coefficients[i],T1,T0);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}


/** Simple integrated centred(!) polynomial function generator divided by independent variable.
 *  We need to rewrite some of the functions in order to
 *  use central fit. Having a central temperature Tbase
 *  allows for a better fit, but requires a different
 *  formulation of the fracInt function group. Other
 *  functions are not affected.
 *  Starts with only the first coefficient at T^0 */

///Helper functions to calculate binomial coefficients: http://rosettacode.org/wiki/Evaluate_binomial_coefficients#C.2B.2B
//double IncompressibleClass::factorial(double nValue){
//   double result = nValue;
//   double result_next;
//   double pc = nValue;
//   do {
//	   result_next = result*(pc-1);
//	   result = result_next;
//	   pc--;
//   } while(pc>2);
//   nValue = result;
//   return nValue;
//}
//double IncompressibleClass::factorial(double nValue){
//	if (nValue == 0) return (1);
//	else return (nValue * factorial(nValue - 1));
//}
double IncompressibleClass::factorial(double nValue){
    double value = 1;
    for(int i = 2; i <= nValue; i++){
        value = value * i;
    }
    return value;
}
double IncompressibleClass::binom(double nValue, double nValue2){
   double result;
   if(nValue2 == 1) return nValue;
   result = (factorial(nValue)) / (factorial(nValue2)*factorial((nValue - nValue2)));
   nValue2 = result;
   return nValue2;
}
///Helper functions to calculate the D vector:
std::vector<double> IncompressibleClass::fracIntCentralDvector(int m, double T, double Tbase){
	std::vector<double> D;
	double tmp;
	if (m<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",m));
	for (int j=0; j<m; j++){ // loop through row
		tmp = pow(-1.0,j) * log(T) * pow(Tbase,(int)j);
		for(int k=0; k<j; k++) { // internal loop for every entry
			tmp += binom(j,k) * pow(-1.0,k) / (j-k) * pow(T,j-k) * pow(Tbase,k);
		}
		D.push_back(tmp);
	}
	return D;
}
std::vector<double> IncompressibleClass::fracIntCentralDvector(int m, double T1, double T0, double Tbase){
	std::vector<double> D;
	double tmp;
	if (m<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",m));
	for (int j=0; j<m; j++){ // loop through row
		tmp = pow(-1.0,(int)j) * log(T1/T0) * pow(Tbase,(double)j);
		for(int k=0; k<j; k++) { // internal loop for every entry
			tmp += binom(j,k) * pow(-1.0,(int)k) / (j-k) * (pow(T1,(int)(j-k))-pow(T0,(int)(j-k))) * pow(Tbase,(int)k);
		}
		D.push_back(tmp);
	}
	return D;
}
///Indefinite integral of a centred polynomial divided by its independent variable
double IncompressibleClass::fracIntCentral(std::vector<double> const& coefficients, double T, double Tbase){
	if (this->DEBUG) {
		std::cout << "Running    fracIntCentral(std::vector, " << T << ", " << Tbase << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	int m = coefficients.size();
	std::vector<double> D = fracIntCentralDvector(m, T, Tbase);
	double result = 0;
	for(int j=0; j<m; j++) {
		result += coefficients[j] * D[j];
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Definite integral from T0 to T1 of a centred polynomial divided by its independent variable
double IncompressibleClass::fracIntCentral(std::vector<double> const& coefficients, double T1, double T0, double Tbase){
	if (this->DEBUG) {
		std::cout << "Running    fracIntCentral(std::vector, " << T1 << ", " << T0 << ", " << Tbase << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	int m = coefficients.size();
	std::vector<double> D = fracIntCentralDvector(m, T1, T0, Tbase);
	double result = 0;
	for(int j=0; j<m; j++) {
		result += coefficients[j] * D[j];
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

/// Horner function generator implementations
/** Represent polynomials according to Horner's scheme.
 *  This avoids unnecessary multiplication and thus
 *  speeds up calculation.
 */
double IncompressibleClass::baseHorner(std::vector<double> const& coefficients, double T){
	if (this->DEBUG) {
		std::cout << "Running       baseHorner(std::vector, " << T << "): ";
	}
	double result = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result *= T;
		result += coefficients[i];
	}
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
double IncompressibleClass::baseHorner(std::vector< std::vector<double> > const& coefficients, double x, double T){
	if (this->DEBUG) {
		std::cout << "Running       baseHorner(std::vector, " << x << ", " << T << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result *= x;
		result += baseHorner(coefficients[i], T);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Indefinite integral in T-direction
double IncompressibleClass::baseHornerInt(std::vector<double> const& coefficients, double T){
	if (this->DEBUG) {
		std::cout << "Running       baseHornerInt(std::vector, " << T << "): ";
	}
	double result = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result *= T;
		result += coefficients[i]/(i+1.);
	}
	result = result * T;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Definite integral from T0 to T1
double IncompressibleClass::baseHornerInt(std::vector<double> const& coefficients, double T1, double T0){
	if (this->DEBUG) {
		std::cout << "Running       baseHornerInt(std::vector, " << T1 << ", " << T0 << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;

	double result1 = 0;
	double result0 = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result1 *= T1;
		result1 += coefficients[i]/(i+1.);
		result0 *= T0;
		result0 += coefficients[i]/(i+1.);
	}
	result1 *= T1;
	result0 *= T0;

	result1 -= result0;

	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result1 << std::endl;
	}
	return result1;
}
///Indefinite integral in T-direction only
double IncompressibleClass::baseHornerInt(std::vector<std::vector<double> > const& coefficients, double x, double T){
	if (this->DEBUG) {
		std::cout << "Running       baseHornerInt(std::vector, " << x << ", " << T << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result *= x;
		result += baseHornerInt(coefficients[i], T);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Definite integral from T0 to T1
double IncompressibleClass::baseHornerInt(std::vector<std::vector<double> > const& coefficients, double x, double T1, double T0){
	if (this->DEBUG) {
		std::cout << "Running       baseHornerInt(std::vector, " << x << ", " << T1 << ", " << T0 << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result *= x;
		result += baseHornerInt(coefficients[i], T1, T0);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
/////Indefinite integral of a polynomial divided by its independent variable
//double IncompressibleClass::baseHornerFra(std::vector<double> const& coefficients, double T);
/////Definite integral from T0 to T1 of a polynomial divided by its independent variable
//double IncompressibleClass::baseHornerFra(std::vector<double> const& coefficients, double T1, double T0);
/////Indefinite integral of a polynomial divided by its 2nd independent variable
//double IncompressibleClass::baseHornerFra(std::vector<std::vector<double> > const& coefficients, double x, double T);
/////Definite integral from T0 to T1 of a polynomial divided by its 2nd independent variable
//double IncompressibleClass::baseHornerFra(std::vector<std::vector<double> > const& coefficients, double x, double T1, double T0);

///Indefinite integral in T-direction of a polynomial divided by its independent variable
double IncompressibleClass::baseHornerFra(std::vector<double> const& coefficients, double T){
	if (this->DEBUG) {
		std::cout << "Running      baseHornerFra(std::vector, " << T << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0;
	if (coefficients.size() > 1) {
		for(int i=coefficients.size()-1; i>=1; i--) {
			result *= T;
			result += coefficients[i]/(i);
		}
		result *= T;
	}
	result += coefficients[0] * log(T);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Definite integral from T0 to T1 of a polynomial divided by its independent variable
double IncompressibleClass::baseHornerFra(std::vector<double> const& coefficients, double T1, double T0){
	if (this->DEBUG) {
		std::cout << "Running      baseHornerFra(std::vector, " << T1 << ", " << T0 << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;

	double result = 0;
	if (coefficients.size() > 1) {
		double result0 = 0;
		for(int i=coefficients.size()-1; i>=1; i--) {
			result  *= T1;
			result  += coefficients[i]/(i);
			result0 *= T0;
			result0 += coefficients[i]/(i);
		}
		result  *= T1;
		result0 *= T0;
		result  -= result0;
	}
	result += coefficients[0] * log(T1/T0);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Indefinite integral in T-direction of a polynomial divided by its 2nd independent variable
double IncompressibleClass::baseHornerFra(std::vector<std::vector<double> > const& coefficients, double x, double T){
	if (this->DEBUG) {
		std::cout << "Running      baseHornerFra(std::vector, " << x << ", " << T << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;

	double result = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result *= x;
		result += baseHornerFra(coefficients[i], T);
	}

	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Definite integral from T0 to T1 of a polynomial divided by its 2nd independent variable
double IncompressibleClass::baseHornerFra(std::vector<std::vector<double> > const& coefficients, double x, double T1, double T0){
	if (this->DEBUG) {
		std::cout << "Running      baseHornerFra(std::vector, " << x << ", " << T1 << ", " << T0 << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;

	double result = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result *= x;
		result += baseHornerFra(coefficients[i], T1, T0);
	}

	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}



/** Integrating coefficients for polynomials is done by dividing the
 *  original coefficients by (i+1) and elevating the order by 1.
 *  Some reslicing needs to be applied to integrate along the x-axis.
 */
std::vector<double> IncompressibleClass::integrateCoeffs(std::vector<double> const& coefficients){
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

///Integrating coefficients for polynomial in terms of x(axis=true) or y(axis=false).
std::vector< std::vector<double> > IncompressibleClass::integrateCoeffs(std::vector< std::vector<double> > const& coefficients, bool axis){
	std::vector< std::vector<double> > newCoefficients;
	unsigned int sizeX = coefficients.size();
	if (sizeX<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",sizeX));

	if (axis==true){
		std::vector< std::vector<double> > tmpCoefficients;
		tmpCoefficients = transpose(coefficients);
		unsigned int sizeY = tmpCoefficients.size();
		for(unsigned int i=0; i<sizeY; i++) {
			newCoefficients.push_back(integrateCoeffs(tmpCoefficients[i]));
		}
		return transpose(newCoefficients);
	} else if (axis==false){
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
std::vector<double> IncompressibleClass::deriveCoeffs(std::vector<double> const& coefficients){
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


/** Alternatives
 *  Simple functions that heavily rely on other parts of this file.
 *  We still need to check which combinations yield the best
 *  performance.
 */
///Indefinite integral in T-direction
double IncompressibleClass::integrateIn2Steps(std::vector<double> const& coefficients, double T){
	if (this->DEBUG) {
		std::cout << "Running   integrateIn2Steps(std::vector, " << T << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result =  polyval(integrateCoeffs(coefficients),T);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Definite integral from T0 to T1
double IncompressibleClass::integrateIn2Steps(std::vector<double> const& coefficients, double T1, double T0){
	if (this->DEBUG) {
		std::cout << "Running   integrateIn2Steps(std::vector, " << T1 << ", " << T0 << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	std::vector<double> coefficientsInt(integrateCoeffs(coefficients));
	double result = polyval(coefficientsInt,T1)-polyval(coefficientsInt,T0);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Indefinite integral in terms of x(axis=true) or T(axis=false).
double IncompressibleClass::integrateIn2Steps(std::vector< std::vector<double> > const& coefficients, double x, double T, bool axis){
	if (this->DEBUG) {
		std::cout << "Running   integrateIn2Steps(std::vector, " << x << ", " << T << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = polyval(integrateCoeffs(coefficients,axis),x,T);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Definite integral from T0 to T1
double IncompressibleClass::integrateIn2Steps(std::vector< std::vector<double> > const& coefficients, double x, double T1, double T0){
	if (this->DEBUG) {
		std::cout << "Running   integrateIn2Steps(std::vector, " << x << ", " << T1 << ", " << T0 << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	std::vector< std::vector<double> > coefficientsInt(integrateCoeffs(coefficients,false));
	double result = polyval(coefficientsInt,x,T1)-polyval(coefficientsInt,x,T0);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Indefinite integral in T-direction of a polynomial divided by its independent variable
double IncompressibleClass::fracIntIn2Steps(std::vector<double> const& coefficients, double T){
	if (this->DEBUG) {
		std::cout << "Running    fracIntIn2Steps(std::vector, " << T << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = coefficients[0] * log(T);
	if (coefficients.size() > 1) {
		std::vector<double> newCoeffs(coefficients.begin() + 1, coefficients.end());
		result += polyint(newCoeffs,T);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Definite integral from T0 to T1 of a polynomial divided by its independent variable
double IncompressibleClass::fracIntIn2Steps(std::vector<double> const& coefficients, double T1, double T0){
	if (this->DEBUG) {
		std::cout << "Running    fracIntIn2Steps(std::vector, " << T1 << ", " << T0 << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = coefficients[0] * log(T1/T0);
	if (coefficients.size() > 1) {
		std::vector<double> newCoeffs(coefficients.begin() + 1, coefficients.end());
		result += polyint(newCoeffs,T1,T0);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Indefinite integral in T-direction of a polynomial divided by its 2nd independent variable
double IncompressibleClass::fracIntIn2Steps(std::vector<std::vector<double> > const& coefficients, double x, double T){
	if (this->DEBUG) {
		std::cout << "Running    fracIntIn2Steps(std::vector, " << x << ", " << T << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	std::vector<double> newCoeffs;
	for (unsigned int i=0; i<coefficients.size(); i++){
		newCoeffs.push_back(polyfracint(coefficients[i],T));
	}
	double result = polyval(newCoeffs,x);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Definite integral from T0 to T1 of a polynomial divided by its 2nd independent variable
double IncompressibleClass::fracIntIn2Steps(std::vector<std::vector<double> > const& coefficients, double x, double T1, double T0){
	if (this->DEBUG) {
		std::cout << "Running    fracIntIn2Steps(std::vector, " << x << ", " << T1 << ", " << T0 << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	std::vector<double> newCoeffs;
	for (unsigned int i=0; i<coefficients.size(); i++){
		newCoeffs.push_back(polyfracint(coefficients[i],T1,T0));
	}
	double result = polyval(newCoeffs,x);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

double IncompressibleClass::fracIntCentral2Steps(std::vector<std::vector<double> > const& coefficients, double x, double T, double Tbase){
	if (this->DEBUG) {
		std::cout << "Running    fracIntCentral2Steps(std::vector, " << x << ", " << T << ", " << Tbase << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	std::vector<double> newCoeffs;
	for (unsigned int i=0; i<coefficients.size(); i++){
		newCoeffs.push_back(fracIntCentral(coefficients[i], T, Tbase));
	}
	double result = polyval(newCoeffs,x);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
///Definite integral from T0 to T1 of a centred polynomial divided by its 2nd independent variable
double IncompressibleClass::fracIntCentral2Steps(std::vector<std::vector<double> > const& coefficients, double x, double T1, double T0, double Tbase){
	if (this->DEBUG) {
		std::cout << "Running    fracIntCentral2Steps(std::vector, " << x << ", " << T1 << ", " << T0 << ", " << Tbase << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	std::vector<double> newCoeffs;
	for (unsigned int i=0; i<coefficients.size(); i++){
		newCoeffs.push_back(fracIntCentral(coefficients[i], T1, T0, Tbase));
	}
	double result = polyval(newCoeffs,x);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}


/** Here we define the functions that should be used by the
 *  respective implementations. Please do no use any other
 *  method since this would break the purpose of this interface.
 */

/// Evaluates an exponential function for the given coefficients
/// @param coefficients vector containing the ordered coefficients
/// @param T double value that represents the current input
/// @param n int value that determines the kind of exponential function
double IncompressibleClass::expval(std::vector<double> const& coefficients, double T, int n){
	double result = 0.;
	if (n==1) {
		checkCoefficients(coefficients,3);
		result = exp(coefficients[0]/(T+coefficients[1]) - coefficients[2]);
	} else if (n==2) {
		result = exp(polyval(coefficients, T));
	} else {
		throw ValueError(format("There is no function defined for this input (%d). ",n));
	}
	return result;
}

/// Evaluates an exponential function for the given coefficients
/// @param coefficients vector containing the ordered coefficients
/// @param x double value that represents the current input in the 1st dimension
/// @param T double value that represents the current input in the 2nd dimension
/// @param n int value that determines the kind of exponential function
double IncompressibleClass::expval(std::vector< std::vector<double> > const& coefficients, double x, double T, int n){
	double result = 0.;
	if (n==2) {
		result = exp(polyval(coefficients, x, T));
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
