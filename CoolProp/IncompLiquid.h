
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
	bool checkT(double T_K){
		if( Tmin < 0. ) {
			throw ValueError("Please specify the minimum temperature.");
		} else if( Tmax < 0.) {
			throw ValueError("Please specify the maximum temperature.");
		} else if ( (Tmin>T_K) or (T_K>Tmax) ) {
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

	/// Exponential function generator.
	/** Base function to produce exponential functions
	 *  based on the length of the coefficients vector.
	 *  An extra check is performed to compare the length
	 *  of the vector to the input n. */
	double baseExponential(std::vector<double> coefficients, double x, unsigned int n){
		double result = 0.;
	    if (coefficients.size() == n){
	    	if (n==3) {
	    		result = exp(coefficients[0]/(x+coefficients[1]) - coefficients[2]);
	    	} else {
	        	throw ValueError(format("There is no function defined for this number of coefficients (%d). ",coefficients.size()));
	        }
	    } else {
	    	throw ValueError(format("The number of coefficients %d does not match with %d. ",coefficients.size(),n));
	    }
	    return result;
	}
};

bool IsIncompressibleLiquid(std::string name);
double IncompLiquid(long iOutput, double T, double p, long iFluid);
double IncompLiquid(long iOutput, double T, double p, std::string name);

/// Base class for simplified models
/** Base class following the structure initially developed by Ian
 *  to fit the data from Melinder-BOOK-2010. Only needs a reduced
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
    	return basePolynomial(cRho, T_K, 2);
    }
    double cp(double T_K, double p){
    	if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
    	return basePolynomial(cCp, T_K, 2);
    }
    double u(double T_K, double p){
    	if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
    	if (cCp.size() != 2) throw ValueError(format("Wrong number of coefficients: %d instead of %d. ",cCp.size(),2));
    	checkP(T_K, p);
        return cCp[1]*(T_K*T_K-298.15*298.15)/2.0 + cCp[0]*(T_K-298.15);
	}
	double s(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		if (cCp.size() != 2) throw ValueError(format("Wrong number of coefficients: %d instead of %d. ",cCp.size(),2));
		return cCp[1]*(T_K-298.15)/2.0 + cCp[0]*log(T_K/298.15);
	}
	double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return exp(basePolynomial(cVisc, T_K, 3));
	}
	double cond(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return basePolynomial(cCond, T_K, 2);
	}
    double h(double T_K, double p){
    	return h_u(T_K,p);
    }
    double psat(double T_K){
    	if (!checkT(T_K)) throw ValueError(format("T=%f is out of range.",T_K));
    	if (T_K<TminPsat || TminPsat<0){
    		return -1.;
    	} else {
    		return baseExponential(cPsat, T_K, 3);
    	}
    };
};

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
};

//Therminol fluids: Therminol2007

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
        cRho.push_back(-7.8062728911E-01);

        cCp.clear();
        cCp.push_back(+9.0795582424E-01);
        cCp.push_back(+4.0657921324E-03);

        cCond.clear();
        cCond.push_back(+1.5921846605E-04);
        cCond.push_back(-1.7117124580E-07);

        cVisc.clear();
        cVisc.push_back(+8.9025486215E+00);
        cVisc.push_back(-7.7561437403E-02);
        cVisc.push_back(+8.6542890469E-05);

        cPsat.clear();
        cPsat.push_back(-3.4575358244E+03);
        cPsat.push_back(-8.2860659840E+01);
        cPsat.push_back(-1.3662144073E+01);
    };
};
//
//class TherminolD12 : public IncompressibleLiquid{
//public:
//	TherminolD12(){
//        name = std::string("TD12");
//        description = std::string("Therminol D12");
//        reference = std::string("BibTex Key");
//    };
//
//    double rho(double T_K, double p){
//    	std::vector<double> cDensity;
//    	cDensity.push_back( 776.257);
//    	cDensity.push_back(-0.696982);
//    	cDensity.push_back(-0.000131384);
//    	cDensity.push_back(-0.00000209079);
//    	return basePolynomial(cDensity, (T_K-273.15));
//    }
//    double cp(double T_K, double p){
//        return 0.00287576*T_K+0.999729;
//    }
//    double u(double T_K, double p){
//        return 0.00287576*(T_K*T_K-298*298)/2.0+0.999729*(T_K-298);
//    }
//    double s(double T_K, double p){
//        return 0.00287576*(T_K-298)/2.0+0.999729*log(T_K/298);
//    }
//    double visc(double T_K, double p){
//        return exp(7.03331e-05*T_K*T_K+-0.0566396*T_K+3.5503);
//    }
//    double cond(double T_K, double p){
//        return -2.06364e-07*T_K+0.000189132;
//    }
//    double h(double T_K, double p){
//    	return h_u(T_K,p);
//    }
//};


#endif
