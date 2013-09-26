
#ifndef INCOMPRESSIBLE_LIQUID_H
#define INCOMPRESSIBLE_LIQUID_H

#include <string>
#include <vector>
#include <math.h>
#include "CPExceptions.h"
#include "CoolPropTools.h"
#include "IncompBase.h"

/**
Notes for developers:

If you want to add a fluid, add its definition to the header 
and then add it to the map in the constructor in LiquidsContainer 
in IncompLiquid.cpp
**/

/// The abstract base class for the liquids
class IncompressibleLiquid : public Incompressible{
	
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
	 *  Tmax since the default values cause an error. */
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
	 *  The default value for psat is -1 yielding true if psat
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
		return fracint(cHeat, T_K, Tref);
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
        cRho.push_back(+9.8912069747E+02);
        cRho.push_back(-9.2980872967E-01);
        cRho.push_back(+1.0537501330E-03);
        cRho.push_back(-1.5498115143E-06);

        cHeat.clear();
        cHeat.push_back(+1.0054715109E+03);
        cHeat.push_back(+3.6645441085E-00);
        cHeat.push_back(-3.5878725087E-04);
        cHeat.push_back(+1.7370907291E-06);

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

        cHeat.clear();
        cHeat.push_back(+6.8006101123E+02);
        cHeat.push_back(+3.2230860459E-00);
        cHeat.push_back(-1.0974690747E-03);
        cHeat.push_back(+8.1520085319E-07);

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

        cHeat.clear();
        cHeat.push_back(+6.9155653776E+02);
        cHeat.push_back(+3.1812352525E-00);
        cHeat.push_back(-1.0707667973E-03);
        cHeat.push_back(+7.8051070556E-07);

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

class Therminol66Class : public SimpleIncompressible{
public:
	Therminol66Class(){
        name = std::string("T66");
        description = std::string("Therminol 66");
        reference = std::string("Therminol2007");

		Tmin     = 273.15;
		Tmax     = 618.15;
		TminPsat = 273.15;

        cRho.clear();
        cRho.push_back(+1164.45337);
        cRho.push_back(-0.4388917);
        cRho.push_back(-0.000321);

        cHeat.clear();
        cHeat.push_back(+657.990904);
        cHeat.push_back(+2.82292602);
        cHeat.push_back(+0.0008970785);

        cCond.clear();
        cCond.push_back(+0.116116312);
        cCond.push_back(+0.000048945);
        cCond.push_back(-1.50000000E-07);

        cVisc.clear();
        cVisc.push_back(+586.375);
        cVisc.push_back(-210.6500E0);
        cVisc.push_back(-2.2809);

        cPsat.clear();
        cPsat.push_back(-9094.51);
        cPsat.push_back(+66.8500E0);
        cPsat.push_back(+24.54486E0);
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

        cHeat.clear();
        cHeat.push_back(+1.1465102196E+03);
        cHeat.push_back(+2.1016260744E-00);
        cHeat.push_back(-2.2151557961E-04);
        cHeat.push_back(+3.5493927846E-06);

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

        cHeat.clear();
        cHeat.push_back(+6.5777357045E+02);
        cHeat.push_back(+3.6209780444E-00);
        cHeat.push_back(-8.3411429568E-04);
        cHeat.push_back(+2.2428632290E-07);

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

        cHeat.clear();
        cHeat.push_back(+7.5687081085E+02);
        cHeat.push_back(+3.9761218594E-00);
        cHeat.push_back(-5.5667939231E-04);
        cHeat.push_back(+2.5261546566E-07);

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

        cHeat.clear();
        cHeat.push_back(+1.3960182000E+03);
        cHeat.push_back(+1.7200000001E-01);
        cHeat.push_back(-9.4083072373E-15);
        cHeat.push_back(+4.3126550326E-18);

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
        cRho.push_back(+1.1563685147E+03);
        cRho.push_back(-1.0269048053E+00);
        cRho.push_back(-9.3505449092E-07);
        cRho.push_back(+1.0368054566E-09);

        cHeat.clear();
        cHeat.push_back(+1.1122510494E+03);
        cHeat.push_back(+2.5171817500E-00);
        cHeat.push_back(-1.2392094684E-03);
        cHeat.push_back(+1.1624033307E-06);

        cCond.clear();
        cCond.push_back(+1.6121957379E-04);
        cCond.push_back(-1.3023781944E-07);
        cCond.push_back(-1.4395238759E-10);

        cVisc.clear();
        cVisc.push_back(+1.0337654975E+03);
        cVisc.push_back(-4.3322764492E+01);
        cVisc.push_back(+1.0715062353E+01);

        cPsat.clear();
    };
};


#endif
