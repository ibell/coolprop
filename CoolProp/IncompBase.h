
#ifndef INCOMPRESSIBLE_BASE_H
#define INCOMPRESSIBLE_BASE_H

#include <string>
#include <vector>
#include <valarray>
#include "CPExceptions.h"
#include "CoolProp.h"
#include <math.h>
#include "Solvers.h"

/// The abstract base class
class IncompressibleFluid{
	
protected:
	std::string name;
	std::string description;
	std::string reference;
	double Tmin;
	double TminPsat;
	double Tmax;
	double Tref;

public:
	// Constructor
	IncompressibleFluid(){
		name = "";
		description = "";
		reference = "";
		Tmin = -1.;
		Tmax = -1.;
		TminPsat = -1.;
		Tref = 273.15 + 25. ;
	};

	// Destructor.  No implementation
	virtual ~IncompressibleFluid(){};

public:
	std::string getDescription() const {
		return description;
	}

	std::string getName() const {
		return name;
	}

	std::string getReference() const {
		return reference;
	}

	double getTmax() const {
		return Tmax;
	}

	double getTmin() const {
		return Tmin;
	}

	double getTminPsat() const {
		return TminPsat;
	}

	double getTref() const {
		return Tref;
	}

protected:
	void setDescription(std::string description) {
		this->description = description;
	}

	void setName(std::string name) {
		this->name = name;
	}

	void setReference(std::string reference) {
		this->reference = reference;
	}

	void setTmax(double Tmax) {
		this->Tmax = Tmax;
	}

	void setTmin(double Tmin) {
		this->Tmin = Tmin;
	}

	void setTminPsat(double TminPsat) {
		this->TminPsat = TminPsat;
	}

	void setTref(double Tref) {
		this->Tref = Tref;
	}


protected:
	/// Polynomial function generator.
	/** Base function to produce n-th order polynomials
	 *  based on the length of the coefficient vector.
	 *  Starts with only the first coefficient at x^0. */
	double basePolynomial(std::vector<double> coefficients, double x){
	    double result = 0.;
	    for(unsigned int i=0; i<coefficients.size();i++) {
	    	result += coefficients[i] * pow(x,(int)i);
	    }
	    return result;
	}

	/// Polynomial function generator with check.
	/** Base function to produce n-th order polynomials
	 *  based on the length of the coefficient vector.
	 *  Starts with only the first coefficient at x^0
	 *  and checks the vector length against parameter n. */
	double basePolynomial(std::vector<double> coefficients, double x, unsigned int n){
		if (coefficients.size() == n){
			return basePolynomial(coefficients, x);
		} else {
			throw ValueError(format("The number of coefficients %d does not match with %d. ",coefficients.size(),n));
		}
	}

	/// 2D polynomial function generator.
	/** Base function to produce n-th order polynomials
	 *  based on the size of the coefficient matrix.
	 *  Starts with only the first coefficient at x^0*y^0. */
	double basePolynomial(std::vector<std::vector<double>> coefficients, double x, double y){
		double result = 0.;
		for(unsigned int i=0; i<coefficients.size();i++) {
			result += pow(x,(int)i) * basePolynomial(coefficients[i],y);
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
				result += pow(x,(int)i) * basePolynomial(coefficients[i],y,o);
			}
			return result;
		} else {
			throw ValueError(format("The number of rows %d does not match with %d. ",coefficients.size(),n));
		}
	}

	/// Integrated polynomial function generator.
	/** Base function to produce integrals of n-th order
	 *  polynomials based on the length of the coefficient
	 *  vector. Integrates from x0 to x1.
	 *  Starts with only the first coefficient at x^0 */
	double basePolynomialInt(std::vector<double> coefficients, double x1, double x0){
		double result = 0.;
		for(unsigned int i=0; i<coefficients.size();i++) {
			result += 1./(i+1.) * coefficients[i] * (pow(x1,(i+1.)) - pow(x0,(i+1.)));
		}
		return result;
	}

	/// Integrated polynomial function generator with check.
	/** Calls the base function but checks the vector
	 *  length against parameter n. */
	double basePolynomialInt(std::vector<double> coefficients, double x1, double x0, unsigned int n){
		if (coefficients.size() == n){
			return basePolynomialInt(coefficients, x1, x0);
		} else {
			throw ValueError(format("The number of coefficients %d does not match with %d. ",coefficients.size(),n));
		}
	}

	/// Integrated 2D polynomial function generator.
	/** Base function to produce integrals of i-th order
	 *  polynomials based on the size of the coefficient
	 *  matrix. Integrates from x0 to x1 along only one axis.
	 *  Starts with only the first coefficient at x^0*y^0 */
	// vector to valarray: std:valarray<double> corpX(corps_tmp[i].data(), corps_tmp[i].size());
    // valarray to vector: corps_tmp[i].assign(std::begin(corpX), std::end(corpX));
	// column major: return m_storage[std::slice(column, m_height, m_stride)][row];
    // row major:    return m_storage[std::slice(row, m_stride, m_height)][column];
//	double basePolynomialIntdx(std::valarray<std::valarray<double>> coefficients, double x1, double y, double x0){
//		double result = 0.;
//		std::vector<double> tmpCoeff;
//		for(unsigned int i=0; i<coefficients[0].size();i++) {
//			tmpCoeff = std::vector(std::begin(coefficients[i]), std::end(coefficients[i]));
//
//			result += pow(y,(int)i) * basePolynomialInt(coefficients[std::slice()][i],x1,x0);
//		}
//		return result;
//	}
	double basePolynomialIntdx(std::vector<std::vector<double>> coefficients, double x1, double y, double x0){
		double result = 0.;
		std::vector<double> tmpCoeff;
		for(unsigned int i=0; i<coefficients[0].size();i++) {
			tmpCoeff.clear();
			for (unsigned int j=0; j<coefficients.size();j++) {
				tmpCoeff.push_back(coefficients[j][i]);
			}
			result += pow(y,(int)i) * basePolynomialInt(tmpCoeff,x1,x0);
		}
		return result;
	}
	/** Base function to produce integrals of i-th order
	 *  polynomials based on the size of the coefficient
	 *  matrix. Integrates from y0 to y1 along only one axis.
	 *  Starts with only the first coefficient at x^0*y^0 */
	double basePolynomialIntdy(std::vector<std::vector<double>> coefficients, double x, double y1, double y0){
		double result = 0.;
		for(unsigned int i=0; i<coefficients.size();i++) {
			result += pow(x,(int)i) * basePolynomialInt(coefficients[i],y1,y0);
		}
		return result;
	}

	/// Integrated fraction generator.
	/** Base function to produce integrals of n-th order
	 *  polynomials divided by their variable based on
	 *  the length of the coefficients vector.
	 *  Integrates from x0 to x1, starts with only the
	 *  first coefficient at x^0 */
	double baseFractionInt(std::vector<double> coefficients, double x1, double x0){
		double result = coefficients[0] * log(x1/x0);
		if (coefficients.size() > 1) {
			std::vector<double> newCoeffs(coefficients.begin() + 1, coefficients.end());
			result += basePolynomialInt(newCoeffs,x1,x0);
		}
		return result;
	}

	/// Integrated fraction generator with check.
	/** Calls the base function but checks the vector
	 *  length against parameter n before. */
	double baseFractionInt(std::vector<double> coefficients, double x1, double x0, unsigned int n){
		if (coefficients.size() == n){
			return baseFractionInt(coefficients, x1, x0);
		} else {
			throw ValueError(format("The number of coefficients %d does not match with %d. ",coefficients.size(),n));
		}
	}

	/// Exponential function generator.
	/** Base function to produce exponential functions
	 *  based on the input n.
	 *  An extra check is performed for the length
	 *  of the vector according to the chosen function. */
	double baseExponential(std::vector<double> coefficients, double x, int n){
		double result = 0.;
		if (n==1) {
			int c = 3;
			if (coefficients.size() != c) throw ValueError(format("The number of coefficients %d does not match with %d. ",coefficients.size(),c));
			result = exp(coefficients[0]/(x+coefficients[1]) - coefficients[2]);
		} else if (n==2) {
			int c = 3;
			if (coefficients.size() != c) throw ValueError(format("The number of coefficients %d does not match with %d. ",coefficients.size(),c));
			result = exp(basePolynomial(coefficients, x, c));
		} else {
			throw ValueError(format("There is no function defined for this input (%d). ",n));
		}
	    return result;
	}
};




#endif
