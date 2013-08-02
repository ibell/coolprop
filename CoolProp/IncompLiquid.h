
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
	double Tmin, TminPsat, Tmax, Tref;

public:
	// Constructor
	IncompressibleLiquid(){
		Tmin = -1.;
		Tmax = -1.;
		TminPsat = -1.;
		Tref = 273.15 + 25. ;
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
	virtual double psat(double T_K){return -1;};
	/// Saturation temperature as a function of pressure.
	virtual double Tsat(double p){return -1;};

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
	 *  density. Provides consistent formulations. */
	double h_u(double T_K, double p) {
		return u(T_K,p)+p/rho(T_K,p);
	};

	/// Internal energy from h,p and rho.
	/** Calculate internal energy as a function of temperature
	 *  and pressure employing functions for enthalpy and
	 *  density. Provides consistent formulations. */
	double u_h(double T_K, double p) {
		return h(T_K,p)-p/rho(T_K,p);
	};

	
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
	bool checkT(double T_K){
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

	/// Polynomial function generator.
	/** Base function to produce n-th order polynomials
	 *  based on the length of the coefficients vector.
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
	 *  based on the length of the coefficients vector.
	 *  Starts with only the first coefficient at x^0
	 *  and checks the vector length against parameter n. */
	double basePolynomial(std::vector<double> coefficients, double x, unsigned int n){
		if (coefficients.size() == n){
			return basePolynomial(coefficients, x);
		} else {
			throw ValueError(format("The number of coefficients %d does not match with %d. ",coefficients.size(),n));
		}
	}

	/// Integrated polynomial function generator.
	/** Base function to produce integrals of n-th order
	 *  polynomials based on the length of the coefficients
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

bool IsIncompressibleLiquid(std::string name);
double IncompLiquid(long iOutput, double T, double p, long iFluid);
double IncompLiquid(long iOutput, double T, double p, std::string name);

/// Base class for simplified models
/** Employs the base functions implemented above, only needs a reduced
 *  set of coefficients for density and heat capacity. The other
 *  quantities are calculated from combinations of the coefficients.
 *  Additionally, extra parameters are used for viscosity and
 *  thermal conductivity. */
class SimpleIncompressible : public IncompressibleLiquid{
protected:
	std::vector<double> cRho;
	std::vector<double> cCp;
	std::vector<double> cVisc;
	std::vector<double> cCond;
	std::vector<double> cPsat;

public:
    double rho(double T_K, double p){
    	if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
    	return basePolynomial(cRho, T_K);
    }
    double cp(double T_K, double p){
    	if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
    	return basePolynomial(cCp, T_K);
    }
    double h(double T_K, double p){
    	if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
    	return basePolynomialInt(cCp, T_K, Tref);
	}
	double s(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return baseFractionInt(cCp, T_K, Tref);
	}
	double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return baseExponential(cVisc, T_K, 1);
	}
	double cond(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return basePolynomial(cCond, T_K);
	}
    double u(double T_K, double p){
    	return u_h(T_K,p);
    }
    double psat(double T_K){
    	if (!checkT(T_K)) throw ValueError(format("T=%f is out of range.",T_K));
    	if (T_K<TminPsat || TminPsat<0){
    		return -1.;
    	} else {
    		return baseExponential(cPsat, T_K, 1);
    	}
    };
};

/*
 * The next classes follow the structure initially developed by Ian
 * to fit the data from Melinder-BOOK-2010.
 */
class DEBLiquidClass : public SimpleIncompressible{
public:
	DEBLiquidClass(){
        name = std::string("DEB");
        description = std::string("Diethylbenzene mixture - Dowtherm J Dow Chemical Co.");
        reference = std::string("Melinder-BOOK-2010");

        Tmin = -80.0 + 273.15;
        Tmax = 100.0 + 273.15;

        cRho.clear();
        cRho.push_back(1076.5);
        cRho.push_back(-0.731182);

        cCp.clear();
        cCp.push_back(0.999729);
        cCp.push_back(0.00287576);

        cVisc.clear();
        cVisc.push_back(3.5503);
        cVisc.push_back(-0.0566396);
        cVisc.push_back(7.03331e-05);

        cCond.clear();
        cCond.push_back(0.000189132);
        cCond.push_back(-2.06364e-07);
    };
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return baseExponential(cVisc, T_K, 2);
	}
    double h(double T_K, double p){
		return h_u(T_K,p);
	}
};

class HCMLiquidClass : public SimpleIncompressible{
public:
	HCMLiquidClass(){
        name = std::string("HCM");
        description = std::string("Hydrocarbon mixture (synthetic) - Therminol D12 (Gilotherm D12) Solutia");
        reference = std::string("Melinder-BOOK-2010");

        Tmin = -80.0 + 273.15;
        Tmax = 100.0 + 273.15;

        cRho.clear();
        cRho.push_back(971.725);
        cRho.push_back(-0.718788);

        cCp.clear();
        cCp.push_back(0.844023);
        cCp.push_back(0.00431212);

        cVisc.clear();
        cVisc.push_back(18.3237);
        cVisc.push_back(-0.14706);
        cVisc.push_back(0.000209096);

        cCond.clear();
        cCond.push_back(0.000153716);
        cCond.push_back(-1.51212e-07);
    };
    double u(double T_K, double p){
    	if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
    	return basePolynomialInt(cCp, T_K, Tref);
	}
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return baseExponential(cVisc, T_K, 2);
	}
    double h(double T_K, double p){
		return h_u(T_K,p);
	}
};

class HFELiquidClass : public SimpleIncompressible{
public:
	HFELiquidClass(){
        name = std::string("HFE");
        description = std::string("Hydrofluoroether - HFE-7100 3M Novec");
        reference = std::string("Melinder-BOOK-2010");

        Tmin = -80.0 + 273.15;
        Tmax = 100.0 + 273.15;

        cRho.clear();
        cRho.push_back(1822.37);
        cRho.push_back(-0.918485);

        cCp.clear();
        cCp.push_back(0.871834);
        cCp.push_back(0.000858788);

        cVisc.clear();
        cVisc.push_back(-4.22878);
        cVisc.push_back(-0.0114765);
        cVisc.push_back(7.39823e-06);

        cCond.clear();
        cCond.push_back(9.92958e-05);
        cCond.push_back(-8.33333e-08);
    };
    double u(double T_K, double p){
    	if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
    	return basePolynomialInt(cCp, T_K, Tref);
	}
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return baseExponential(cVisc, T_K, 2);
	}
    double h(double T_K, double p){
		return h_u(T_K,p);
	}
};

class PMS1LiquidClass : public SimpleIncompressible{
public:
	PMS1LiquidClass(){
        name = std::string("PMS1");
        description = std::string("Polydimethylsiloxan 1. - Baysilone KT3");
        reference = std::string("Melinder-BOOK-2010");

        Tmin = -80.0 + 273.15;
        Tmax = 100.0 + 273.15;

        cRho.clear();
        cRho.push_back(1172.35);
        cRho.push_back(-0.9025);

        cCp.clear();
        cCp.push_back(1.22369);
        cCp.push_back(0.00148417);

        cVisc.clear();
        cVisc.push_back(6.36183);
        cVisc.push_back(-0.0636352);
        cVisc.push_back(7.51428e-05);

        cCond.clear();
        cCond.push_back(0.000207526);
        cCond.push_back(-2.84167e-07);
    };
    double u(double T_K, double p){
    	if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
    	return basePolynomialInt(cCp, T_K, Tref);
	}
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return baseExponential(cVisc, T_K, 2);
	}
    double h(double T_K, double p){
		return h_u(T_K,p);
	}
};

class PMS2LiquidClass : public SimpleIncompressible{
public:
	PMS2LiquidClass(){
        name = std::string("PMS2");
        description = std::string("Polydimethylsiloxan 2. - Syltherm XLT Dow Corning Co.");
        reference = std::string("Melinder-BOOK-2010");

        Tmin = -80.0 + 273.15;
        Tmax = 100.0 + 273.15;

        cRho.clear();
        cRho.push_back(1155.94);
        cRho.push_back(-1.02576);

        cCp.clear();
        cCp.push_back(1.15355);
        cCp.push_back(0.00210788);

        cVisc.clear();
        cVisc.push_back(5.66926);
        cVisc.push_back(-0.065582);
        cVisc.push_back(8.09988e-05);

        cCond.clear();
        cCond.push_back(0.000172305);
        cCond.push_back(-2.11212e-07);
    };
    double u(double T_K, double p){
    	if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
    	return basePolynomialInt(cCp, T_K, Tref);
	}
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return baseExponential(cVisc, T_K, 2);
	}
    double h(double T_K, double p){
		return h_u(T_K,p);
	}
};

class SABLiquidClass : public SimpleIncompressible{
public:
	SABLiquidClass(){
        name = std::string("SAB");
        description = std::string("Synthetic alkyl benzene - Marlotherm X");
        reference = std::string("Melinder-BOOK-2010");

        Tmin = -80.0 + 273.15;
        Tmax = 100.0 + 273.15;

        cRho.clear();
        cRho.push_back(1102.34);
        cRho.push_back(-0.801667);

        cCp.clear();
        cCp.push_back(1.36094);
        cCp.push_back(0.00151667);

        cVisc.clear();
        cVisc.push_back(5.21288);
        cVisc.push_back(-0.0665792);
        cVisc.push_back(8.5066e-05);

        cCond.clear();
        cCond.push_back(0.000208374);
        cCond.push_back(-2.61667e-07);
    };
    double u(double T_K, double p){
    	if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
    	return basePolynomialInt(cCp, T_K, Tref);
	}
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return baseExponential(cVisc, T_K, 2);
	}
    double h(double T_K, double p){
		return h_u(T_K,p);
	}
};

class HCBLiquidClass : public SimpleIncompressible{
public:
	HCBLiquidClass(){
        name = std::string("HCB");
        description = std::string("Hydrocarbon blend - Dynalene MV");
        reference = std::string("Melinder-BOOK-2010");

		Tmin = -80.0 + 273.15;
		Tmax = 100.0 + 273.15;

        cRho.clear();
        cRho.push_back(1071.78);
        cRho.push_back(-0.772024);

        cCp.clear();
        cCp.push_back(0.761393);
        cCp.push_back(0.00352976);

        cVisc.clear();
        cVisc.push_back(7.16819);
        cVisc.push_back(-0.0863212);
        cVisc.push_back(0.000130604);

        cCond.clear();
        cCond.push_back(0.000203186);
        cCond.push_back(-2.3869e-07);
    };
    double u(double T_K, double p){
    	if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
    	return basePolynomialInt(cCp, T_K, Tref);
	}
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return baseExponential(cVisc, T_K, 2);
	}
    double h(double T_K, double p){
		return h_u(T_K,p);
	}
};

class TCOLiquidClass : public SimpleIncompressible{
public:
    TCOLiquidClass(){
        name = std::string("TCO");
        description = std::string("Terpene from citrus oils - d-Limonene");
        reference = std::string("Melinder-BOOK-2010");

		Tmin = -80.0 + 273.15;
		Tmax = 100.0 + 273.15;

        cRho.clear();
        cRho.push_back(1071.02);
        cRho.push_back(-0.778166);

        cCp.clear();
        cCp.push_back(0.223775);
        cCp.push_back(0.0052159);

        cVisc.clear();
        cVisc.push_back(-3.47971);
        cVisc.push_back(-0.0107031);
        cVisc.push_back(1.14086e-06);

        cCond.clear();
        cCond.push_back(0.000174156);
        cCond.push_back(-1.85052e-07);
    };
    double u(double T_K, double p){
    	if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
    	return basePolynomialInt(cCp, T_K, Tref);
	}
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return baseExponential(cVisc, T_K, 2);
	}
    double h(double T_K, double p){
		return h_u(T_K,p);
	}
};


/// New fluids added with more coefficients
/** New data for a few fluids. Most of these models employ
 *  an extended set of parameters consisting of a
 *  3rd order polynomial for density and heat capacity, a
 *  2nd order polynomial for thermal conductivity as well as
 *  an exponential function for viscosity.
 *  I rewrote the base class to match this new form since I
 *  expect all the new fluids to follow this pattern. The
 *  "dev" folder contains Python scripts to fit functions to
 *  data. Have a look there and you will see how easy it is
 *  to extend the fluid database. */

/// Therminol Fluids
/** Data sheets for most Therminol (Solutia) fluids are
 *  available from their homepage and we will implement
 *  some of them as liquid only (!) and incompressible
 *  heat transfer media. */
class TherminolD12Class : public SimpleIncompressible{
public:
	TherminolD12Class(){
        name = std::string("TD12");
        description = std::string("Therminol D12");
        reference = std::string("Therminol2007");

		Tmin     = 188.15;
		Tmax     = 503.15;
		TminPsat = 188.15;

        cRho.clear();
        cRho.push_back(+9.8912069747E+02);
        cRho.push_back(-9.2980872967E-01);
        cRho.push_back(+1.0537501330E-03);
        cRho.push_back(-1.5498115143E-06);

        cCp.clear();
        cCp.push_back(+1.0054715109E+00);
        cCp.push_back(+3.6645441085E-03);
        cCp.push_back(-3.5878725087E-07);
        cCp.push_back(+1.7370907291E-09);

        cCond.clear();
        cCond.push_back(+1.4137409270E-04);
        cCond.push_back(-5.9980348657E-08);
        cCond.push_back(-1.6084318930E-10);

        cVisc.clear();
        cVisc.push_back(+5.6521090296E+02);
        cVisc.push_back(-1.2468427867E+02);
        cVisc.push_back(+1.0064677591E+01);

        cPsat.clear();
        cPsat.push_back(-3.4575358244E+03);
        cPsat.push_back(-8.2860659840E+01);
        cPsat.push_back(-1.3662144073E+01);
    };
};

class TherminolVP1Class : public SimpleIncompressible{
public:
	TherminolVP1Class(){
        name = std::string("TVP1");
        description = std::string("Therminol VP-1");
        reference = std::string("Therminol2007");

		Tmin     = 285.15;
		Tmax     = 670.15;
		TminPsat = 285.15;

        cRho.clear();
        cRho.push_back(+1.4027135186E+03);
        cRho.push_back(-1.6132555941E+00);
        cRho.push_back(+2.1378377260E-03);
        cRho.push_back(-1.9310694457E-06);

        cCp.clear();
        cCp.push_back(+6.8006101123E-01);
        cCp.push_back(+3.2230860459E-03);
        cCp.push_back(-1.0974690747E-06);
        cCp.push_back(+8.1520085319E-10);

        cCond.clear();
        cCond.push_back(+1.4898820998E-04);
        cCond.push_back(+7.4237593304E-09);
        cCond.push_back(-1.7298441018E-10);

        cVisc.clear();
        cVisc.push_back(+1.0739262832E+03);
        cVisc.push_back(-8.3841430000E+01);
        cVisc.push_back(+1.0616847342E+01);

        cPsat.clear();
        cPsat.push_back(-4.3138714899E+03);
        cPsat.push_back(-8.7431401773E+01);
        cPsat.push_back(-1.4358490536E+01);
    };
};

class Therminol72Class : public SimpleIncompressible{
public:
	Therminol72Class(){
        name = std::string("T72");
        description = std::string("Therminol 72");
        reference = std::string("Therminol2007");

		Tmin     = 263.15;
		Tmax     = 653.15;
		TminPsat = 263.15;

        cRho.clear();
        cRho.push_back(+1.3673153689E+03);
        cRho.push_back(-1.0611213813E+00);
        cRho.push_back(+3.3724996764E-04);
        cRho.push_back(-2.3561359575E-07);

        cCp.clear();
        cCp.push_back(+6.9155653776E-01);
        cCp.push_back(+3.1812352525E-03);
        cCp.push_back(-1.0707667973E-06);
        cCp.push_back(+7.8051070556E-10);

        cCond.clear();
        cCond.push_back(+1.7514206629E-04);
        cCond.push_back(-1.2131347168E-07);
        cCond.push_back(-4.0980786601E-14);

        cVisc.clear();
        cVisc.push_back(+6.8390609525E+02);
        cVisc.push_back(-1.8024922198E+02);
        cVisc.push_back(+1.0066296789E+01);

        cPsat.clear();
        cPsat.push_back(-1.7535987063E+06);
        cPsat.push_back(+9.8874585029E+03);
        cPsat.push_back(-1.7271828471E+02);
    };
};

class DowthermJClass : public SimpleIncompressible{
public:
	DowthermJClass(){
        name = std::string("DowJ");
        description = std::string("Dowtherm J");
        reference = std::string("Dow Chemicals data sheet");

        Tmin     = 193.15;
        Tmax     = 618.15;
        TminPsat = 323.15;

        cRho.clear();
        cRho.push_back(+1.1413344279E+03);
        cRho.push_back(-1.4313342250E+00);
        cRho.push_back(+2.4904467725E-03);
        cRho.push_back(-2.9222650181E-06);

        cCp.clear();
        cCp.push_back(+1.1465102196E+00);
        cCp.push_back(+2.1016260744E-03);
        cCp.push_back(-2.2151557961E-07);
        cCp.push_back(+3.5493927846E-09);

        cCond.clear();
        cCond.push_back(+1.8990894157E-04);
        cCond.push_back(-2.0921055014E-07);
        cCond.push_back(-3.2093835108E-12);

        cVisc.clear();
        cVisc.push_back(+7.0617903908E+02);
        cVisc.push_back(-6.4144266810E+01);
        cVisc.push_back(+1.0083412137E+01);

        cPsat.clear();
        cPsat.push_back(-3.1876260646E+03);
        cPsat.push_back(-9.7927481216E+01);
        cPsat.push_back(-1.3576473644E+01);
    };
};

class DowthermQClass : public SimpleIncompressible{
public:
	DowthermQClass(){
        name = std::string("DowQ");
        description = std::string("Dowtherm Q");
        reference = std::string("Dow Chemicals data sheet");

        Tmin     = 238.15;
        Tmax     = 633.15;
        TminPsat = 393.15;

        cRho.clear();
        cRho.push_back(+1.2033276517E+03);
        cRho.push_back(-8.6829872708E-01);
        cRho.push_back(+2.4832944517E-04);
        cRho.push_back(-1.7756195502E-07);

        cCp.clear();
        cCp.push_back(+6.5777357045E-01);
        cCp.push_back(+3.6209780444E-03);
        cCp.push_back(-8.3411429568E-07);
        cCp.push_back(+2.2428632290E-10);

        cCond.clear();
        cCond.push_back(+1.5381076522E-04);
        cCond.push_back(-9.0635332826E-08);
        cCond.push_back(-6.2655520375E-11);

        cVisc.clear();
        cVisc.push_back(+8.2860901385E+02);
        cVisc.push_back(-1.2328762577E+02);
        cVisc.push_back(+1.0441389876E+01);

        cPsat.clear();
        cPsat.push_back(-2.8419253596E+03);
        cPsat.push_back(-1.7104225035E+02);
        cPsat.push_back(-1.2287974574E+01);
    };
};


class Texatherm22Class : public SimpleIncompressible{
public:
	Texatherm22Class(){
        name = std::string("TX22");
        description = std::string("Texatherm22");
        reference = std::string("Texaco data sheet");

        Tmin     = 273.15;
        Tmax     = 623.15;
        TminPsat = 313.15;

        cRho.clear();
        cRho.push_back(+1.0828544822E+03);
        cRho.push_back(-9.4455188837E-01);
        cRho.push_back(+9.2414382570E-04);
        cRho.push_back(-9.5365432963E-07);

        cCp.clear();
        cCp.push_back(+7.5687081085E-01);
        cCp.push_back(+3.9761218594E-03);
        cCp.push_back(-5.5667939231E-07);
        cCp.push_back(+2.5261546566E-10);

        cCond.clear();
        cCond.push_back(+1.5246860039E-04);
        cCond.push_back(-5.9875212755E-08);
        cCond.push_back(-1.4202281813E-11);

        cVisc.clear();
        cVisc.push_back(+8.8295946693E+02);
        cVisc.push_back(-1.7261015850E+02);
        cVisc.push_back(+9.6466061591E+00);

        cPsat.clear();
        cPsat.push_back(-8.8967639918E+03);
        cPsat.push_back(-4.3465054898E+01);
        cPsat.push_back(-1.7472312925E+01);
    };
};


class NitrateSaltClass : public SimpleIncompressible{
public:
	NitrateSaltClass(){
        name = std::string("NaK");
        description = std::string("60% NaNO3 and 40% KNO3");
        reference = std::string("Solar Power Tower Design Basis Document,  Alexis B. Zavoico, Sandia Labs, USA");

        Tmin     = 573.15;
        Tmax     = 873.15;
        TminPsat = 873.15;

        cRho.clear();
        cRho.push_back(+2.3537711479E+03);
        cRho.push_back(-1.0161832139E+00);
        cRho.push_back(+5.2992982291E-04);
        cRho.push_back(-2.4393584453E-07);

        cCp.clear();
        cCp.push_back(+1.3960182000E+00);
        cCp.push_back(+1.7200000001E-04);
        cCp.push_back(-9.4083072373E-18);
        cCp.push_back(+4.3126550326E-21);

        cCond.clear();
        cCond.push_back(+3.9112042040E-04);
        cCond.push_back(+1.8994695528E-07);
        cCond.push_back(-1.7231213093E-14);

        cVisc.clear();
        cVisc.push_back(+4.7467257729E+02);
        cVisc.push_back(-3.3943569667E+02);
        cVisc.push_back(+7.7431109317E+00);

        cPsat.clear();
    };
};


#endif
