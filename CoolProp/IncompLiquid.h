
#ifndef INCOMPRESSIBLE_LIQUID_H
#define INCOMPRESSIBLE_LIQUID_H

#include <string>
#include <vector>
#include <math.h>
#include "CPExceptions.h"
#include "CoolPropTools.h"
#include "IncompBase.h"
#include <stdio.h>

/**
Notes for developers:

If you want to add a fluid, add its definition to the header 
and then add it to the map in the constructor in LiquidsContainer 
in IncompLiquid.cpp
**/

/// Base class for simplified brine/solution models
/** Employs the base functions implemented in IncompBase.h and
 *  provides properties as function of temperature, pressure
 *  and composition. */
class IncompressibleLiquid : public IncompressibleClass{
	
public:
	/* All functions need T and p as input. Might not be necessary,
	 * but gives a clearer structure.
	 */
	/// Density as a function of temperature and pressure.
	virtual double rho (double T_K, double p){return -_HUGE;};
	/// Heat capacities as a function of temperature and pressure.
	virtual double c   (double T_K, double p){return -_HUGE;};
	virtual double cp  (double T_K, double p){return c(T_K,p);};
	virtual double cv  (double T_K, double p){return c(T_K,p);};
	/// Entropy as a function of temperature and pressure.
	virtual double s   (double T_K, double p){return -_HUGE;};
	/// Internal energy as a function of temperature and pressure.
	virtual double u   (double T_K, double p){return -_HUGE;};
	/// Enthalpy as a function of temperature and pressure.
	virtual double h   (double T_K, double p){return -_HUGE;};
	/// Viscosity as a function of temperature and pressure.
	virtual double visc(double T_K, double p){return -_HUGE;};
	/// Thermal conductivity as a function of temperature and pressure.
	virtual double cond(double T_K, double p){return -_HUGE;};
	/// Saturation pressure as a function of temperature.
	virtual double psat(double T_K){return -_HUGE;};

    void testInputs(double T_K, double p);

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
	/** Compares the given temperature T to a stored minimum
	 *  and maximum temperature. Enforces the redefinition of
	 *  Tmin and Tmax since the default values cause an error. */
	bool checkT(double T_K);

	/// Check validity of pressure input.
	/** Compares the given pressure p to the saturation pressure at
	 *  temperature T and throws and exception if p is lower than
	 *  the saturation conditions.
	 *  The default value for psat is -1 yielding true if psat
	 *  is not redefined in the subclass.
	 *  */
	bool checkP(double T_K, double p);

	/// Check validity of temperature and pressure input.
	bool checkTP(double T, double p);
};


/** Basic functions to access the list of incompressible fluids.
 *  Used here for convenience, but does not really contribute
 *  any functionality.
 */
bool IsIncompressibleLiquid(std::string name);
double IncompLiquidSI(long iOutput, double T, double p, long iFluid);
double IncompLiquidSI(long iOutput, double T, double p, std::string name);
// only two functions needed
// no name processing
// no concentration issues


/** Handle all the objects in a single list of incompressible
 *  liquids. */
class LiquidsContainer {
private:
  std::vector<IncompressibleLiquid*> liquid_list;
  std::map<std::string,IncompressibleLiquid*> liquid_map;

public:
  IncompressibleLiquid * get_liquid(long index);
  IncompressibleLiquid * get_liquid(std::string name);
  void set_liquids(std::vector<IncompressibleLiquid*> list);

  LiquidsContainer();
  ~LiquidsContainer();
};


/// Base class for simplified models
/** Employs the base functions implemented above, only needs a reduced
 *  set of coefficients for density and heat capacity. The other
 *  quantities are calculated from combinations of the coefficients.
 *  Additionally, extra parameters are used for viscosity and
 *  thermal conductivity. */
class SimpleIncompressible : public IncompressibleLiquid{
protected:
	std::vector<double> cRho;
	std::vector<double> cHeat;
	std::vector<double> cVisc;
	std::vector<double> cCond;
	std::vector<double> cPsat;

public:
    double rho(double T_K, double p){
    	checkTP(T_K, p);
    	return polyval(cRho, T_K);
    }
    double c(double T_K, double p){
    	checkTP(T_K, p);
    	return polyval(cHeat, T_K);
    }
    double h(double T_K, double p){
    	checkTP(T_K, p);
    	return h_u(T_K,p);
	}
	double s(double T_K, double p){
		checkTP(T_K, p);
		return polyfracint(cHeat, T_K, Tref);
	}
	double visc(double T_K, double p){
		checkTP(T_K, p);
		return expval(cVisc, T_K, 1);
	}
	double cond(double T_K, double p){
		checkTP(T_K, p);
		return polyval(cCond, T_K);
	}
    double u(double T_K, double p){
    	return polyint(cHeat, T_K, Tref);
    }
    double psat(double T_K){
    	checkT(T_K);
    	if (T_K<TminPsat || TminPsat<0){
    		return -1.;
    	} else {
    		return expval(cPsat, T_K, 1);
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

        cHeat.clear();
        cHeat.push_back(999.729);
        cHeat.push_back(2.87576);

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
		return expval(cVisc, T_K, 2);
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

        cHeat.clear();
        cHeat.push_back(844.023);
        cHeat.push_back(4.31212);

        cVisc.clear();
        cVisc.push_back(18.3237);
        cVisc.push_back(-0.14706);
        cVisc.push_back(0.000209096);

        cCond.clear();
        cCond.push_back(0.000153716);
        cCond.push_back(-1.51212e-07);
    };
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return expval(cVisc, T_K, 2);
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

        cHeat.clear();
        cHeat.push_back(871.834);
        cHeat.push_back(858788);

        cVisc.clear();
        cVisc.push_back(-4.22878);
        cVisc.push_back(-0.0114765);
        cVisc.push_back(7.39823e-06);

        cCond.clear();
        cCond.push_back(9.92958e-05);
        cCond.push_back(-8.33333e-08);
    };
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return expval(cVisc, T_K, 2);
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

        cHeat.clear();
        cHeat.push_back(1223.69);
        cHeat.push_back(1.48417);

        cVisc.clear();
        cVisc.push_back(6.36183);
        cVisc.push_back(-0.0636352);
        cVisc.push_back(7.51428e-05);

        cCond.clear();
        cCond.push_back(0.000207526);
        cCond.push_back(-2.84167e-07);
    };
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return expval(cVisc, T_K, 2);
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

        cHeat.clear();
        cHeat.push_back(1153.55);
        cHeat.push_back(2.10788);

        cVisc.clear();
        cVisc.push_back(5.66926);
        cVisc.push_back(-0.065582);
        cVisc.push_back(8.09988e-05);

        cCond.clear();
        cCond.push_back(0.000172305);
        cCond.push_back(-2.11212e-07);
    };
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return expval(cVisc, T_K, 2);
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

        cHeat.clear();
        cHeat.push_back(1360.94);
        cHeat.push_back(1.51667);

        cVisc.clear();
        cVisc.push_back(5.21288);
        cVisc.push_back(-0.0665792);
        cVisc.push_back(8.5066e-05);

        cCond.clear();
        cCond.push_back(0.000208374);
        cCond.push_back(-2.61667e-07);
    };
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return expval(cVisc, T_K, 2);
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

        cHeat.clear();
        cHeat.push_back(761.393);
        cHeat.push_back(3.52976);

        cVisc.clear();
        cVisc.push_back(7.16819);
        cVisc.push_back(-0.0863212);
        cVisc.push_back(0.000130604);

        cCond.clear();
        cCond.push_back(0.000203186);
        cCond.push_back(-2.3869e-07);
    };
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return expval(cVisc, T_K, 2);
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

        cHeat.clear();
        cHeat.push_back(223.775);
        cHeat.push_back(5.2159);

        cVisc.clear();
        cVisc.push_back(-3.47971);
        cVisc.push_back(-0.0107031);
        cVisc.push_back(1.14086e-06);

        cCond.clear();
        cCond.push_back(0.000174156);
        cCond.push_back(-1.85052e-07);
    };
    double visc(double T_K, double p){
		if (!checkTP(T_K, p)) throw ValueError(format("T=%f or p=%f is out of range.",T_K,p));
		return expval(cVisc, T_K, 2);
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
        cRho.push_back(+9.8351112723E+02);
        cRho.push_back(-9.2980873989E-01);
        cRho.push_back(+1.0537525196E-03);
        cRho.push_back(-1.5498113539E-06);

        cHeat.clear();
        cHeat.push_back(+9.8364476562E+02);
        cHeat.push_back(+3.8726050231E+00);
        cHeat.push_back(-9.8766309005E-04);
        cHeat.push_back(+2.3435575994E-06);

        cCond.clear();
        cCond.push_back(+1.4137409273E-01);
        cCond.push_back(-5.9980348818E-05);
        cCond.push_back(-1.6084318907E-07);

        cVisc.clear();
        cVisc.push_back(+5.6521089774E+02);
        cVisc.push_back(-1.2468427934E+02);
        cVisc.push_back(+1.0064677576E+01);

        cPsat.clear();
        cPsat.push_back(-3.4575358232E+03);
        cPsat.push_back(-8.2860659877E+01);
        cPsat.push_back(-2.0569899350E+01);
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
        cRho.push_back(+1.4027134791E+03);
        cRho.push_back(-1.6132556454E+00);
        cRho.push_back(+2.1378385004E-03);
        cRho.push_back(-1.9310694268E-06);

        cHeat.clear();
        cHeat.push_back(+2.8811134920E+02);
        cHeat.push_back(+5.8749383495E+00);
        cHeat.push_back(-6.8565802314E-03);
        cHeat.push_back(+4.8441771885E-06);

        cCond.clear();
        cCond.push_back(+1.4898820987E-01);
        cCond.push_back(+7.4237598012E-06);
        cCond.push_back(-1.7298441066E-07);

        cVisc.clear();
        cVisc.push_back(+1.0739262698E+03);
        cVisc.push_back(-8.3841432169E+01);
        cVisc.push_back(+1.0616847324E+01);

        cPsat.clear();
        cPsat.push_back(-4.3138714911E+03);
        cPsat.push_back(-8.7431401731E+01);
        cPsat.push_back(-2.1266245816E+01);
    };
};

class Therminol66Class : public SimpleIncompressible{
public:
	Therminol66Class(){
        name = std::string("T66");
        description = std::string("Therminol 66");
        reference = std::string("Therminol2007");

        Tmin     = 273.15;
        Tmax     = 653.15;
        TminPsat = 343.15;

        cRho.clear();
        cRho.push_back(+1.1644533740E+03);
        cRho.push_back(-4.3889170000E-01);
        cRho.push_back(-3.2100000000E-04);
        cRho.push_back(+3.7806951741E-20);

        cHeat.clear();
        cHeat.push_back(+6.5799090444E+02);
        cHeat.push_back(+2.8229260154E+00);
        cHeat.push_back(+8.9707850000E-04);
        cHeat.push_back(-6.3169106168E-20);

        cCond.clear();
        cCond.push_back(+1.1611631163E-01);
        cCond.push_back(+4.8945000000E-05);
        cCond.push_back(-1.5000000000E-07);

        cVisc.clear();
        cVisc.push_back(+6.6720362621E+02);
        cVisc.push_back(-2.0480017928E+02);
        cVisc.push_back(+9.5933675483E+00);

        cPsat.clear();
        cPsat.push_back(-9.0945100000E+03);
        cPsat.push_back(+6.6850000000E+01);
        cPsat.push_back(-2.4544855279E+01);
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
        cRho.push_back(+1.3571809600E+03);
        cRho.push_back(-9.8961574320E-01);
        cRho.push_back(+1.7605076030E-04);
        cRho.push_back(-1.1893027931E-07);

        cHeat.clear();
        cHeat.push_back(+7.5732470240E+02);
        cHeat.push_back(+2.7131176015E+00);
        cHeat.push_back(-6.5480236953E-06);
        cHeat.push_back(+4.2717093140E-09);

        cCond.clear();
        cCond.push_back(+1.7514206624E-01);
        cCond.push_back(-1.2131347146E-04);
        cCond.push_back(-4.0981053641E-11);

        cVisc.clear();
        cVisc.push_back(+6.8390591135E+02);
        cVisc.push_back(-1.8024924396E+02);
        cVisc.push_back(+1.0066296341E+01);

        cPsat.clear();
        cPsat.push_back(-2.9571373614E+05);
        cPsat.push_back(+3.7936374754E+03);
        cPsat.push_back(-7.9704232489E+01);
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
        cRho.push_back(+1.1413344369E+03);
        cRho.push_back(-1.4313342989E+00);
        cRho.push_back(+2.4904469643E-03);
        cRho.push_back(-2.9222651755E-06);

        cHeat.clear();
        cHeat.push_back(+9.6069502370E+02);
        cHeat.push_back(+3.6462333255E+00);
        cHeat.push_back(-4.2068387567E-03);
        cHeat.push_back(+6.7827145865E-06);

        cCond.clear();
        cCond.push_back(+1.8990894140E-01);
        cCond.push_back(-2.0921054918E-04);
        cCond.push_back(-3.2093847177E-09);

        cVisc.clear();
        cVisc.push_back(+7.0729353166E+02);
        cVisc.push_back(-6.3966539111E+01);
        cVisc.push_back(+1.0085461875E+01);

        cPsat.clear();
        cPsat.push_back(-3.1876142878E+03);
        cPsat.push_back(-9.7928074744E+01);
        cPsat.push_back(-2.0484211718E+01);
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
        cRho.push_back(+1.2033275827E+03);
        cRho.push_back(-8.6829833766E-01);
        cRho.push_back(+2.4832881863E-04);
        cRho.push_back(-1.7756119683E-07);

        cHeat.clear();
        cHeat.push_back(+6.8024152924E+02);
        cHeat.push_back(+3.4510813383E+00);
        cHeat.push_back(-4.2748801499E-04);
        cHeat.push_back(-8.5970499813E-08);

        cCond.clear();
        cCond.push_back(+1.5381076524E-01);
        cCond.push_back(-9.0635332892E-05);
        cCond.push_back(-6.2655520296E-08);

        cVisc.clear();
        cVisc.push_back(+8.2860901780E+02);
        cVisc.push_back(-1.2328762540E+02);
        cVisc.push_back(+1.0441389885E+01);

        cPsat.clear();
        cPsat.push_back(-2.8419559652E+03);
        cPsat.push_back(-1.7104073646E+02);
        cPsat.push_back(-1.9195781229E+01);
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
        cRho.push_back(+1.0828544667E+03);
        cRho.push_back(-9.4455186919E-01);
        cRho.push_back(+9.2414399492E-04);
        cRho.push_back(-9.5365423381E-07);

        cHeat.clear();
        cHeat.push_back(+7.7399470257E+02);
        cHeat.push_back(+3.8528705501E+00);
        cHeat.push_back(-2.7313597680E-04);
        cHeat.push_back(+4.3191182489E-08);

        cCond.clear();
        cCond.push_back(+1.5246860010E-01);
        cCond.push_back(-5.9875211524E-05);
        cCond.push_back(-1.4202283025E-08);

        cVisc.clear();
        cVisc.push_back(+8.8295948920E+02);
        cVisc.push_back(-1.7261015666E+02);
        cVisc.push_back(+9.6466062231E+00);

        cPsat.clear();
        cPsat.push_back(-8.8969171641E+03);
        cPsat.push_back(-4.3461866340E+01);
        cPsat.push_back(-2.4380261252E+01);
    };
};


class NitrateSaltClass : public SimpleIncompressible{
public:
	NitrateSaltClass(){
        name = std::string("NaK");
        description = std::string("Salt mixture with 60% NaNO3 and 40% KNO3");
        reference = std::string("Zavoico2001");

        Tmin     = 573.15;
        Tmax     = 873.15;
        TminPsat = 873.15;

        cRho.clear();
        cRho.push_back(+2.2637234000E+03);
        cRho.push_back(-6.3600000000E-01);
        cRho.push_back(-4.4160301079E-16);
        cRho.push_back(+2.0083875178E-19);

        cHeat.clear();
        cHeat.push_back(+1.3960182000E+03);
        cHeat.push_back(+1.7200000000E-01);
        cHeat.push_back(-3.4511849814E-17);
        cHeat.push_back(+1.7449850640E-20);

        cCond.clear();
        cCond.push_back(+3.9110150000E-01);
        cCond.push_back(+1.9000000000E-04);
        cCond.push_back(+6.2250839227E-21);

        cVisc.clear();
        cVisc.push_back(+4.7467256848E+02);
        cVisc.push_back(-3.3943569983E+02);
        cVisc.push_back(+7.7431109204E+00);

        cPsat.clear();
    };
};


class SylthermXLTClass : public SimpleIncompressible{
public:
	SylthermXLTClass(){

        name = std::string("XLT");
        description = std::string("SylthermXLT");
        reference = std::string("Dow Chemicals data sheet");

        Tmin     = 173.15;
        Tmax     = 533.15;
        TminPsat = 533.15;

        cRho.clear();
        cRho.push_back(+1.1563685145E+03);
        cRho.push_back(-1.0269048032E+00);
        cRho.push_back(-9.3506079577E-07);
        cRho.push_back(+1.0368116627E-09);

        cHeat.clear();
        cHeat.push_back(+1.1562261074E+03);
        cHeat.push_back(+2.0994549103E+00);
        cHeat.push_back(+7.7175381057E-07);
        cHeat.push_back(-3.7008444051E-20);

        cCond.clear();
        cCond.push_back(+1.6121957379E-01);
        cCond.push_back(-1.3023781944E-04);
        cCond.push_back(-1.4395238766E-07);

        cVisc.clear();
        cVisc.push_back(+1.0337654989E+03);
        cVisc.push_back(-4.3322764383E+01);
        cVisc.push_back(+1.0715062356E+01);

        cPsat.clear();
    };
};

class HC50Class : public SimpleIncompressible{
public:
        HC50Class(){

        description = std::string("Dynalene HC-50");
        name = std::string("HC50");
		reference = std::string("Dynalene data sheet");

		Tmin     = 223.15;
		Tmax     = 483.15;
		TminPsat = 293.15;

		cRho.clear();
		cRho.push_back(+1.4989450835E+03);
		cRho.push_back(-5.2796479536E-01);
		cRho.push_back(-7.1686735997E-05);
		cRho.push_back(+6.2219602450E-08);

		cHeat.clear();
		cHeat.push_back(+2.1287827711E+03);
		cHeat.push_back(+1.9224638196E+00);
		cHeat.push_back(+1.3287279132E-04);
		cHeat.push_back(-1.2116448898E-07);

		cCond.clear();
		cCond.push_back(+2.1115985069E-01);
		cCond.push_back(+1.0044355501E-03);
		cCond.push_back(-6.8418171866E-09);

		cVisc.clear();
		cVisc.push_back(+5.1474948873E+02);
		cVisc.push_back(-1.2991405965E+02);
		cVisc.push_back(+8.8804895031E+00);

		cPsat.clear();
		cPsat.push_back(-4.1833595311E+03);
		cPsat.push_back(-3.3779925774E+01);
		cPsat.push_back(-2.3219027215E+01);
    };
};

class HC40Class : public SimpleIncompressible{
public:
        HC40Class(){

        description = std::string("Dynalene HC-40");
        name = std::string("HC40");
        reference = std::string("Dynalene data sheet");

        Tmin     = 233.15;
        Tmax     = 473.15;
        TminPsat = 293.15;

        cRho.clear();
        cRho.push_back(+1.4720776473E+03);
        cRho.push_back(-5.0388465311E-01);
        cRho.push_back(-1.4487525769E-04);
        cRho.push_back(+1.2228923117E-07);

        cHeat.clear();
        cHeat.push_back(+2.2849444547E+03);
        cHeat.push_back(+2.1363723550E+00);
        cHeat.push_back(+3.3322790115E-04);
        cHeat.push_back(-2.2099478511E-07);

        cCond.clear();
        cCond.push_back(+2.1585000000E-01);
        cCond.push_back(+1.0000000000E-03);
        cCond.push_back(-1.0218829918E-20);

        cVisc.clear();
        cVisc.push_back(+6.7794306641E+02);
        cVisc.push_back(-1.0098293303E+02);
        cVisc.push_back(+9.4355407493E+00);

        cPsat.clear();
        cPsat.push_back(-5.5466979146E+03);
        cPsat.push_back(+1.0982652329E+01);
        cPsat.push_back(-2.5451884216E+01);
    };
};

class HC30Class : public SimpleIncompressible{
public:
        HC30Class(){

        description = std::string("Dynalene HC-30");
        name = std::string("HC30");
        reference = std::string("Dynalene data sheet");

        Tmin     = 243.15;
        Tmax     = 483.15;
        TminPsat = 293.15;

        cRho.clear();
        cRho.push_back(+1.4153034908E+03);
        cRho.push_back(-4.4327434107E-01);
        cRho.push_back(-1.5443642107E-04);
        cRho.push_back(+1.1429794039E-07);

        cHeat.clear();
        cHeat.push_back(+2.3846023310E+03);
        cHeat.push_back(+2.4376868197E+00);
        cHeat.push_back(-3.5495726496E-04);
        cHeat.push_back(+3.2206119163E-07);

        cCond.clear();
        cCond.push_back(+2.2585000000E-01);
        cCond.push_back(+1.0000000000E-03);
        cCond.push_back(-7.8145179951E-21);

        cVisc.clear();
        cVisc.push_back(+1.4791319182E+03);
        cVisc.push_back(+4.3364527896E+00);
        cVisc.push_back(+1.0940969046E+01);

        cPsat.clear();
        cPsat.push_back(-3.9183978747E+03);
        cPsat.push_back(-4.4556553119E+01);
        cPsat.push_back(-2.3041842217E+01);
    };
};

class HC20Class : public SimpleIncompressible{
public:
        HC20Class(){

        description = std::string("DynaleneHC-20");
        name = std::string("HC20");
        reference = std::string("Dynalene data sheet");

        Tmin     = 253.15;
        Tmax     = 483.15;
        TminPsat = 293.15;

        cRho.clear();
        cRho.push_back(+1.3918918541E+03);
        cRho.push_back(-5.3737521962E-01);
        cRho.push_back(+4.1692735606E-05);
        cRho.push_back(-2.8527564759E-08);

        cHeat.clear();
        cHeat.push_back(+2.5186678435E+03);
        cHeat.push_back(+2.3663436544E+00);
        cHeat.push_back(-9.6739411954E-06);
        cHeat.push_back(+1.8768134708E-09);

        cCond.clear();
        cCond.push_back(+2.2985000000E-01);
        cCond.push_back(+1.0000000000E-03);
        cCond.push_back(+4.2380113949E-21);

        cVisc.clear();
        cVisc.push_back(+1.4573630736E+03);
        cVisc.push_back(+8.3287350365E+00);
        cVisc.push_back(+1.0986724162E+01);

        cPsat.clear();
        cPsat.push_back(-4.2012148268E+03);
        cPsat.push_back(-3.2285491186E+01);
        cPsat.push_back(-2.3529315156E+01);
    };
};

class HC10Class : public SimpleIncompressible{
public:
        HC10Class(){

        description = std::string("Dynalene HC-10");
        name = std::string("HC10");
        reference = std::string("Dynalene data sheet");

        Tmin     = 263.15;
        Tmax     = 491.15;
        TminPsat = 293.15;

        cRho.clear();
        cRho.push_back(+1.3210573130E+03);
        cRho.push_back(-4.2690361588E-01);
        cRho.push_back(-9.9246921529E-05);
        cRho.push_back(+1.1284336656E-07);

        cHeat.clear();
        cHeat.push_back(+2.5991164581E+03);
        cHeat.push_back(+2.4340563753E+00);
        cHeat.push_back(+1.2102227066E-04);
        cHeat.push_back(-1.0475863139E-07);

        cCond.clear();
        cCond.push_back(+2.3085000000E-01);
        cCond.push_back(+1.0000000000E-03);
        cCond.push_back(+4.8188007943E-22);

        cVisc.clear();
        cVisc.push_back(+1.3099267208E+03);
        cVisc.push_back(-5.1123036317E+00);
        cVisc.push_back(+1.0880904782E+01);

        cPsat.clear();
        cPsat.push_back(-4.0554351446E+03);
        cPsat.push_back(-3.8590604800E+01);
        cPsat.push_back(-2.3416001339E+01);
    };
};

class AS10Class : public SimpleIncompressible{
public:
        AS10Class(){

        	name = std::string("AS10");
        	description = std::string("Aspen Temper -10");
        	reference = std::string("SecCool Software");

        	Tmin     = 265.0;
        	Tmax     = 300.0;
        	TminPsat = 300.0;

        	cRho.clear();
        	cRho.push_back(+1.1446000000E+03);
        	cRho.push_back(-2.0000000000E-01);
        	cRho.push_back(+1.0520504534E-16);
        	cRho.push_back(-7.3416990543E-20);

        	cHeat.clear();
        	cHeat.push_back(+2.0019857143E+03);
        	cHeat.push_back(+9.2664285714E+00);
        	cHeat.push_back(-1.3285714286E-02);
        	cHeat.push_back(-1.3575719937E-16);

        	cCond.clear();
        	cCond.push_back(+1.2439404762E-01);
        	cCond.push_back(+1.3752380952E-03);
        	cCond.push_back(+1.9047619048E-07);

        	cVisc.clear();
        	cVisc.push_back(+4.7558729987E+02);
        	cVisc.push_back(-1.6103694674E+02);
        	cVisc.push_back(+1.0129670635E+01);

        	cPsat.clear();
    };
};

class AS20Class : public SimpleIncompressible{
public:
        AS20Class(){

        	name = std::string("AS20");
        	description = std::string("Aspen Temper -20");
        	reference = std::string("SecCool Software");

        	Tmin     = 255.0;
        	Tmax     = 300.0;
        	TminPsat = 300.0;

        	cRho.clear();
        	cRho.push_back(+1.1052186014E+03);
        	cRho.push_back(+4.9151515152E-01);
        	cRho.push_back(-1.1118881119E-03);
        	cRho.push_back(-4.6620046617E-07);

        	cHeat.clear();
        	cHeat.push_back(+1.5218579021E+03);
        	cHeat.push_back(+1.0060303030E+01);
        	cHeat.push_back(-1.4137529138E-02);
        	cHeat.push_back(+2.3310023312E-06);

        	cCond.clear();
        	cCond.push_back(+1.1390909091E-01);
        	cCond.push_back(+1.3430303030E-03);
        	cCond.push_back(-6.7271158343E-20);

        	cVisc.clear();
        	cVisc.push_back(+4.8319746232E+02);
        	cVisc.push_back(-1.5863708418E+02);
        	cVisc.push_back(+9.9041817702E+00);

        	cPsat.clear();

    };
};

class AS30Class : public SimpleIncompressible{
public:
        AS30Class(){

        	name = std::string("AS30");
        	description = std::string("Aspen Temper -30");
        	reference = std::string("SecCool Software");

        	Tmin     = 245.0;
        	Tmax     = 300.0;
        	TminPsat = 300.0;

        	cRho.clear();
        	cRho.push_back(+8.7111547342E+02);
        	cRho.push_back(+3.7745669146E+00);
        	cRho.push_back(-1.3928293928E-02);
        	cRho.push_back(+1.5747215747E-05);

        	cHeat.clear();
        	cHeat.push_back(+3.0785717616E+02);
        	cHeat.push_back(+1.7014884005E+01);
        	cHeat.push_back(-2.4269064269E-02);
        	cHeat.push_back(-3.4188034188E-06);

        	cCond.clear();
        	cCond.push_back(+1.1394827672E-01);
        	cCond.push_back(+1.2806543457E-03);
        	cCond.push_back(-4.4955044955E-08);

        	cVisc.clear();
        	cVisc.push_back(+5.1209144775E+02);
        	cVisc.push_back(-1.5626579962E+02);
        	cVisc.push_back(+9.9016607949E+00);

        	cPsat.clear();

    };
};

class AS40Class : public SimpleIncompressible{
public:
        AS40Class(){

        	name = std::string("AS40");
        	description = std::string("Aspen Temper -40");
        	reference = std::string("SecCool Software");

        	Tmin     = 235.0;
        	Tmax     = 300.0;
        	TminPsat = 300.0;

        	cRho.clear();
        	cRho.push_back(+8.8355972263E+02);
        	cRho.push_back(+4.0994176412E+00);
        	cRho.push_back(-1.5309807839E-02);
        	cRho.push_back(+1.7359111477E-05);

        	cHeat.clear();
        	cHeat.push_back(-5.1241354234E+02);
        	cHeat.push_back(+2.3184299720E+01);
        	cHeat.push_back(-3.7767937944E-02);
        	cHeat.push_back(-1.2066365007E-06);

        	cCond.clear();
        	cCond.push_back(+1.3896631868E-01);
        	cCond.push_back(+1.1308846154E-03);
        	cCond.push_back(-6.0439560440E-08);

        	cVisc.clear();
        	cVisc.push_back(+4.7403455568E+02);
        	cVisc.push_back(-1.6081056499E+02);
        	cVisc.push_back(+9.4982567299E+00);

        	cPsat.clear();

    };
};

class AS55Class : public SimpleIncompressible{
public:
        AS55Class(){

        	name = std::string("AS55");
        	description = std::string("Aspen Temper -55");
        	reference = std::string("SecCool Software");

        	Tmin     = 220.0;
        	Tmax     = 300.0;
        	TminPsat = 300.0;

        	cRho.clear();
        	cRho.push_back(+8.9605196078E+02);
        	cRho.push_back(+4.5156492948E+00);
        	cRho.push_back(-1.7100103199E-02);
        	cRho.push_back(+1.9435844513E-05);

        	cHeat.clear();
        	cHeat.push_back(+3.6600598555E+02);
        	cHeat.push_back(+1.5792254902E+01);
        	cHeat.push_back(-2.4545923633E-02);
        	cHeat.push_back(-4.1279669761E-07);

        	cCond.clear();
        	cCond.push_back(+3.3985681115E-01);
        	cCond.push_back(-3.0877966976E-04);
        	cCond.push_back(+2.2822497420E-06);

        	cVisc.clear();
        	cVisc.push_back(+6.1353626424E+02);
        	cVisc.push_back(-1.5081016752E+02);
        	cVisc.push_back(+1.0109522524E+01);

        	cPsat.clear();

    };
};

class ZS10Class : public SimpleIncompressible{
public:
        ZS10Class(){

        	name = std::string("ZS10");
        	description = std::string("Zitrec S -10");
        	reference = std::string("SecCool Software");

        	Tmin     = 265.0;
        	Tmax     = 360.0;
        	TminPsat = 360.0;

        	cRho.clear();
        	cRho.push_back(+1.2314887434E+03);
        	cRho.push_back(-5.2375274970E-01);
        	cRho.push_back(+1.9728931723E-04);
        	cRho.push_back(-2.1700379757E-07);

        	cHeat.clear();
        	cHeat.push_back(+2.2183875526E+03);
        	cHeat.push_back(+7.4287485687E+00);
        	cHeat.push_back(-9.1301807180E-03);
        	cHeat.push_back(-1.7295039503E-07);

        	cCond.clear();
        	cCond.push_back(+1.2534548872E-01);
        	cCond.push_back(+1.4936591479E-03);
        	cCond.push_back(-3.0576441103E-07);

        	cVisc.clear();
        	cVisc.push_back(+5.2976289768E+02);
        	cVisc.push_back(-1.5197832234E+02);
        	cVisc.push_back(+1.0073901959E+01);

        	cPsat.clear();

    };
};

class ZS25Class : public SimpleIncompressible{
public:
        ZS25Class(){

        	name = std::string("ZS25");
        	description = std::string("Zitrec S -25");
        	reference = std::string("SecCool Software");

        	Tmin     = 250.0;
        	Tmax     = 360.0;
        	TminPsat = 360.0;

        	cRho.clear();
        	cRho.push_back(+1.3386515650E+03);
        	cRho.push_back(-5.2656222039E-01);
        	cRho.push_back(-1.5645224340E-05);
        	cRho.push_back(+1.2161751291E-08);

        	cHeat.clear();
        	cHeat.push_back(+3.0812417591E+03);
        	cHeat.push_back(+3.7786772648E-01);
        	cHeat.push_back(-6.5881944143E-05);
        	cHeat.push_back(+3.3242120199E-07);

        	cCond.clear();
        	cCond.push_back(+3.4681146245E-02);
        	cCond.push_back(+1.9330395257E-03);
        	cCond.push_back(-1.1992094862E-06);

        	cVisc.clear();
        	cVisc.push_back(+4.1761158398E+02);
        	cVisc.push_back(-1.6295200411E+02);
        	cVisc.push_back(+9.3205072397E+00);

        	cPsat.clear();

    };
};

class ZS40Class : public SimpleIncompressible{
public:
        ZS40Class(){

        	name = std::string("ZS40");
        	description = std::string("Zitrec S -40");
        	reference = std::string("SecCool Software");

        	Tmin     = 235.0;
        	Tmax     = 360.0;
        	TminPsat = 360.0;

        	cRho.clear();
        	cRho.push_back(+1.4183190376E+03);
        	cRho.push_back(-5.8243188235E-01);
        	cRho.push_back(+5.4922721789E-05);
        	cRho.push_back(-6.5779076274E-08);

        	cHeat.clear();
        	cHeat.push_back(+9.5472807318E+02);
        	cHeat.push_back(+1.6485063272E+01);
        	cHeat.push_back(-5.3205833713E-02);
        	cHeat.push_back(+6.2103990740E-05);

        	cCond.clear();
        	cCond.push_back(+1.5446717460E-01);
        	cCond.push_back(+1.1457880342E-03);
        	cCond.push_back(-1.5579975580E-07);

        	cVisc.clear();
        	cVisc.push_back(+4.6544465976E+02);
        	cVisc.push_back(-1.5503900008E+02);
        	cVisc.push_back(+9.2536810485E+00);

        	cPsat.clear();

    };
};

class ZS45Class : public SimpleIncompressible{
public:
        ZS45Class(){

        	name = std::string("ZS45");
        	description = std::string("Zitrec S -45");
        	reference = std::string("SecCool Software");

        	Tmin     = 230.0;
        	Tmax     = 360.0;
        	TminPsat = 360.0;

        	cRho.clear();
        	cRho.push_back(+1.4317629296E+03);
        	cRho.push_back(-5.3989162562E-01);
        	cRho.push_back(-4.2356111320E-05);
        	cRho.push_back(+4.7155909223E-08);

        	cHeat.clear();
        	cHeat.push_back(+2.0449839861E+03);
        	cHeat.push_back(+2.2070679223E+00);
        	cHeat.push_back(-4.6369977429E-06);
        	cHeat.push_back(+5.2395454720E-09);

        	cCond.clear();
        	cCond.push_back(+2.0207530462E-01);
        	cCond.push_back(+8.0384160667E-04);
        	cCond.push_back(+3.6382468107E-07);

        	cVisc.clear();
        	cVisc.push_back(+4.5803862071E+02);
        	cVisc.push_back(-1.5592748501E+02);
        	cVisc.push_back(+9.1329431082E+00);

        	cPsat.clear();

    };
};

class ZS55Class : public SimpleIncompressible{
public:
        ZS55Class(){

        	name = std::string("ZS55");
        	description = std::string("Zitrec S -55");
        	reference = std::string("SecCool Software");

        	Tmin     = 220.0;
        	Tmax     = 360.0;
        	TminPsat = 360.0;

        	cRho.clear();
        	cRho.push_back(+1.4861637422E+03);
        	cRho.push_back(-5.7747246466E-01);
        	cRho.push_back(-9.0036956000E-05);
        	cRho.push_back(+9.9350655525E-08);

        	cHeat.clear();
        	cHeat.push_back(+2.3791157676E+03);
        	cHeat.push_back(-7.7643011749E-01);
        	cHeat.push_back(+4.7957532518E-03);
        	cHeat.push_back(-6.3381598423E-08);

        	cCond.clear();
        	cCond.push_back(+2.4488072144E-01);
        	cCond.push_back(+5.5666864417E-04);
        	cCond.push_back(+5.3163126578E-07);

        	cVisc.clear();
        	cVisc.push_back(+5.0238871796E+02);
        	cVisc.push_back(-1.5126152304E+02);
        	cVisc.push_back(+9.1546618222E+00);

        	cPsat.clear();

    };
};

#endif
