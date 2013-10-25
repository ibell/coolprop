
#ifndef INCOMPRESSIBLE_SOLUTION_H
#define INCOMPRESSIBLE_SOLUTION_H

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
and then add it to the map in the constructor in SolutionsContainer
in IncompSolution.cpp
**/

/// Base class for simplified brine/solution models
/** Employs the base functions implemented in IncompBase.h and
 *  provides properties as function of temperature, pressure
 *  and composition. */
class IncompressibleSolution : public IncompressibleClass{

/** A few redefinitions are needed. The base class does not support
 *  composition as input. We also have to have an additional check
 *  for maximum and minimum concentration as well as the freezing point.
 */
protected:
	double xmin, xmax;
	double Tbase, xbase;
	std::vector< std::vector<double> > cRho;
	std::vector< std::vector<double> > cHeat;
	std::vector< std::vector<double> > cVisc;
	std::vector< std::vector<double> > cCond;
	std::vector< std::vector<double> > cPsat;
	             std::vector<double>   cTfreeze;

public:
	// Constructor
	IncompressibleSolution(){
		xmin  = -1.;
		xmax  = -1.;
		Tbase = -1.;
		xbase = -1.;
	};

	double getTbase() const {return Tbase;}
	double getxbase() const {return xbase;}
	std::vector<std::vector<double> > getcRho() const {return cRho;}
	std::vector<std::vector<double> > getcHeat() const {return cHeat;}
	std::vector<std::vector<double> > getcVisc() const {return cVisc;}
	std::vector<std::vector<double> > getcCond() const {return cCond;}
	std::vector<std::vector<double> > getcPsat() const {return cPsat;}
	            std::vector<double>   getcTfreeze() const {return cTfreeze;}

	/* All functions need T, p and x as input. Might not
	 * be necessary, but gives a clearer structure.
	 */
	/// Density as a function of temperature, pressure and composition.
    virtual double rho (double T_K, double p, double x){return -_HUGE;};
    /// Heat capacities as a function of temperature, pressure and composition.
    virtual double c   (double T_K, double p, double x){return -_HUGE;};
	virtual double cp  (double T_K, double p, double x){return c(T_K,p,x);};
	virtual double cv  (double T_K, double p, double x){return c(T_K,p,x);};
    /// Entropy as a function of temperature, pressure and composition.
    virtual double s   (double T_K, double p, double x){return -_HUGE;};
    /// Internal energy as a function of temperature, pressure and composition.
    virtual double u   (double T_K, double p, double x){return -_HUGE;};
    /// Enthalpy as a function of temperature, pressure and composition.
    virtual double h   (double T_K, double p, double x){return -_HUGE;};
    /// Viscosity as a function of temperature, pressure and composition.
    virtual double visc(double T_K, double p, double x){return -_HUGE;};
    /// Thermal conductivity as a function of temperature, pressure and composition.
    virtual double cond(double T_K, double p, double x){return -_HUGE;};
    /// Saturation pressure as a function of temperature and composition.
    virtual double psat(double T_K          , double x){return -_HUGE;};
    /// Freezing temperature as a function of pressure and composition.
    virtual double Tfreeze(         double p, double x){return -_HUGE;};

    void testInputs(double T_K, double p, double x);


protected:
    /* Define internal energy and enthalpy as functions of the
	 * other properties to provide data in case there are no
	 * coefficients.
	 */

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
	 *  defined for all fluids, default values do not cause errors. */
	bool checkT(double T_K, double p, double x);

	/// Check validity of pressure input.
	/** Compares the given pressure p to the saturation pressure at
	 *  temperature T and throws and exception if p is lower than
	 *  the saturation conditions.
	 *  The default value for psat is -1 yielding true if psat
	 *  is not redefined in the subclass.
	 *  */
	bool checkP(double T_K, double p, double x);

	/// Check validity of composition input.
	/** Compares the given composition x to a stored minimum and
	 *  maximum value. Enforces the redefinition of xmin and
	 *  xmax since the default values cause an error. */
	bool checkX(double x);

	/// Check validity of temperature, pressure and composition input.
	bool checkTPX(double T, double p, double x);
};


/** Basic functions to access the list of incompressible fluids.
 *  Used here for convenience, but does not really contribute
 *  any functionality.
 */
bool IsIncompressibleSolution(std::string name);
double IncompSolution(long iOutput, double T, double p, double x, long iFluid);
double IncompSolution(long iOutput, double T, double p, double x, std::string name);
double IncompSolution(long iOutput, double T, double p, std::string name);  // TODO Solutions: Remove as soon as possible
std::string getSolutionName(std::string name); // TODO Solutions: Remove as soon as possible
double getSolutionConc(std::string name); // TODO Solutions: Remove as soon as possible


/** Handle all the objects in a single list of incompressible
 *  solutions and brines. */
class SolutionsContainer {
private:
  std::vector<IncompressibleSolution*> solution_list;
  std::map<std::string,IncompressibleSolution*> solution_map;

public:
  IncompressibleSolution * get_solution(long index);
  IncompressibleSolution * get_solution(std::string name);
  void set_solutions(std::vector<IncompressibleSolution*> list);

  SolutionsContainer();
  ~SolutionsContainer();
};


/// Class to use Melinder and SecCool parameters
/** Employs some basic wrapper-like functionality
 *  to bridge the gap between the solution functions
 *  used in CoolProp and the definition used in
 *  Melinder's book and Ian's original implementation and
 *  the definition used in SecCool. Please visit:
 *  http://en.ipu.dk/Indhold/refrigeration-and-energy-technology/seccool.aspx
 *  Many thanks to Morten Juel Skovrup for providing
 *  this nice piece of software as well as the parameters
 *  needed to calculate the composition based properties. */
class BaseSolution : public IncompressibleSolution{

protected:
	/// Some more general purpose functions
	double baseFunction(std::vector<double> const& coefficients, double T_K, double p, double x);
	std::vector< std::vector<double> > makeMatrix(std::vector<double> const& coefficients);

public:

	double getTInput(double curTValue){
		return curTValue-Tbase;
	}

	double getxInput(double curxValue){
		return (curxValue-xbase)*100.0;
	}

	double rho(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cRho,6,4);
		return polyval(cRho, getxInput(x), getTInput(T_K));
	}
	double c(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cHeat,6,4);
		return polyval(cHeat, getxInput(x), getTInput(T_K));
	}
	double h(double T_K, double p, double x){
		return h_u(T_K,p,x);
	}
	double s(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cHeat,6,4);
		return polyfracint(cHeat, getxInput(x), T_K, Tref);
	}
//	double s_alt(double T_K, double p, double x){
//		checkTPX(T_K, p, x);
//		IncompressibleClass::checkCoefficients(cHeat,6,4);
//		x = getxInput(x);
//		return (-1/2)*pow(x,5)*cHeat[5][2]*pow(Tref,2)-x*cHeat[1][1]*Tref-pow(x,4)*cHeat[4][1]*Tref-pow(x,5)*cHeat[5][1]*Tref-log(Tref)*pow(x,5)*cHeat[5][0]-log(Tref)*pow(x,4)*cHeat[4][0]-log(Tref)*x*cHeat[1][0]-log(Tref)*pow(x,2)*cHeat[2][0]+(1/2)*pow(x,5)*cHeat[5][2]*pow(T_K,2)-cHeat[0][1]*Tref-log(Tref)*cHeat[0][0]+(1/3)*cHeat[0][3]*pow(T_K,3)+(1/2)*cHeat[0][2]*pow(T_K,2)+cHeat[0][1]*T_K+log(T_K)*cHeat[0][0]+(1/2)*pow(x,2)*cHeat[2][2]*pow(T_K,2)+log(T_K)*x*cHeat[1][0]+pow(x,4)*cHeat[4][1]*T_K+pow(x,5)*cHeat[5][1]*T_K+(1/3)*pow(x,2)*cHeat[2][3]*pow(T_K,3)+(1/3)*pow(x,3)*cHeat[3][3]*pow(T_K,3)+log(T_K)*pow(x,3)*cHeat[3][0]+log(T_K)*pow(x,5)*cHeat[5][0]+(1/3)*pow(x,4)*cHeat[4][3]*pow(T_K,3)+(1/3)*x*cHeat[1][3]*pow(T_K,3)+(1/3)*pow(x,5)*cHeat[5][3]*pow(T_K,3)-log(Tref)*pow(x,3)*cHeat[3][0]+x*cHeat[1][1]*T_K+pow(x,2)*cHeat[2][1]*T_K+pow(x,3)*cHeat[3][1]*T_K-(1/3)*pow(x,2)*cHeat[2][3]*pow(Tref,3)-(1/3)*pow(x,3)*cHeat[3][3]*pow(Tref,3)-pow(x,3)*cHeat[3][1]*Tref-(1/3)*pow(x,4)*cHeat[4][3]*pow(Tref,3)-(1/3)*x*cHeat[1][3]*pow(Tref,3)-(1/3)*pow(x,5)*cHeat[5][3]*pow(Tref,3)-(1/2)*x*cHeat[1][2]*pow(Tref,2)-(1/2)*pow(x,2)*cHeat[2][2]*pow(Tref,2)-(1/2)*pow(x,3)*cHeat[3][2]*pow(Tref,2)-(1/2)*pow(x,4)*cHeat[4][2]*pow(Tref,2)+log(T_K)*pow(x,4)*cHeat[4][0]+(1/2)*x*cHeat[1][2]*pow(T_K,2)+log(T_K)*pow(x,2)*cHeat[2][0]+(1/2)*pow(x,3)*cHeat[3][2]*pow(T_K,2)+(1/2)*pow(x,4)*cHeat[4][2]*pow(T_K,2)-pow(x,2)*cHeat[2][1]*Tref-(1/3)*cHeat[0][3]*pow(Tref,3)-(1/2)*cHeat[0][2]*pow(Tref,2);
//	}
	double visc(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cVisc,6,4);
		return expval(cVisc, getxInput(x), getTInput(T_K), 2)/1e3;
	}
	double cond(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cCond,6,4);
		return polyval(cCond, getxInput(x), getTInput(T_K));
	}
	double u(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cHeat,6,4);
		return polyint(cHeat, getxInput(x), T_K, Tref);
	}
//	double psat(double T_K, double x){
//		//checkT(T_K,p,x);
//		Incompressible::checkCoefficients(cPsat,6,4);
//		if (T_K<TminPsat || TminPsat<0){
//			return -1.;
//		} else {
//			return expval(cPsat, getxInput(x), getTInput(T_K), 2);
//		}
//	};
//	double Tfreeze(double p, double x){
//		IncompressibleClass::checkCoefficients(cTfreeze,5);
//		std::vector<double> tmpVector(cTfreeze.begin()+1,cTfreeze.end());
//		return polyval(tmpVector, x*100.0-cTfreeze[0])+273.15;
//	}
};


class MelinderSolution : public BaseSolution{

protected:
	/// Convert pre-v4.0-style coefficient array to new format
	std::vector<std::vector<double> > convertCoeffs(double *oldestCoeffs, int A, int B);

public:
	MelinderSolution(){

        name = std::string("MelinderSolution");
        description = std::string("Test Methanol Melinder");
        reference = std::string("Test");

		Tmin     = -50 + 273.15;
		Tmax     =  40 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.6;

		Tbase    =  3.5359 + 273.15;
	    xbase    = 30.5128 / 100.0;

	    const int lengthA = 18;
	    const int lengthB =  5;

	    double oldCoeffs[lengthA][lengthB]={
		{ -26.29           	, 958.1				,3887			,   0.4175				,   1.153			},//
		{  -0.000002575    	,  -0.4151			,   7.201		,   0.0007271			,  -0.03866			},
		{  -0.000006732    	,  -0.002261		,  -0.08979		,   0.0000002823		,   0.0002779		},
		{   0.000000163		,   0.0000002998	,  -0.000439	,   0.000000009718		,  -0.000001543		},
		{  -1.187			,  -1.391			, -18.5			,  -0.004421			,   0.005448		},//
		{  -0.00001609		,  -0.0151			,   0.2984		,  -0.00002952			,   0.0001008		},
		{   0.000000342		,   0.0001113		,  -0.001865	,   0.00000007336		,  -0.000002809		},
		{   0.0000000005687	,  -0.0000003264	,  -0.00001718	,   0.0000000004328		,   0.000000009811	},
		{  -0.01218			,  -0.01105			,  -0.03769		,   0.00002044			,  -0.0005552		},//
		{   0.0000003865	,   0.0001828		,  -0.01196		,   0.0000003413		,   0.000008384		},
		{   0.000000008768	,  -0.000001641		,   0.00009801	,  -0.000000003665		,  -0.00000003997	},
		{  -0.0000000002095	,   0.0000000151	,   0.000000666	,  -0.00000000002791	,  -0.0000000003466	},
		{  -0.00006823		,  -0.0001208		,  -0.003776	,   0.0000002943		,   0.000003038		},//
		{   0.00000002137	,   0.000002992		,  -0.00005611	,  -0.0000000009646		,  -0.00000007435	},
		{  -0.0000000004271	,   0.000000001455	,  -0.0000007811,   0.00000000003174	,   0.0000000007442	},
		{   0.0000001297	,   0.000004927		,  -0.0001504	,  -0.0000000008666		,   0.00000006669	},//
		{  -0.0000000005407	,  -0.0000001325	,   0.000007373	,  -0.0000000000004573	,  -0.0000000009105	},
		{   0.00000002363	,  -0.00000007727	,   0.000006433	,  -0.0000000002033		,  -0.0000000008472	}};//

	    std::vector<std::vector<double> > tmpVector = convertCoeffs( *oldCoeffs, lengthA, lengthB);

	    cTfreeze.clear();
	    cTfreeze = column(makeMatrix(tmpVector[0]),0); // Discard temperature coefficients.

        cRho.clear();
        cRho = makeMatrix(tmpVector[1]);

		cHeat.clear();
		cHeat = makeMatrix(tmpVector[2]);

		cCond.clear();
		cCond = makeMatrix(tmpVector[3]);

		cVisc.clear();
		cVisc = makeMatrix(tmpVector[4]);
    };

	/// Define freezing point calculations
	double Tfreeze(double p, double x){
		IncompressibleClass::checkCoefficients(cTfreeze,6);
		return polyval(cTfreeze, getxInput(x))+273.15;
	}
};


class SecCoolSolution : public BaseSolution{
public:
	SecCoolSolution(){

		std::vector<double> tmpVector;

        name = std::string("SecCoolSolution");
        description = std::string("Test Methanol SecCool");
        reference = std::string("Test");

		Tmin     = -50 + 273.15;
		Tmax     =  20 + 273.15;
		TminPsat = Tmax;

		xmin     = 0.0;
		xmax     = 0.5;

		Tbase    = -4.48 + 273.15;
	    xbase    = 31.57 / 100.0;

	    tmpVector.clear();
	    tmpVector.push_back( 960.24665800);
	    tmpVector.push_back(-1.2903839100);
	    tmpVector.push_back(-0.0161042520);
	    tmpVector.push_back(-0.0001969888);
	    tmpVector.push_back( 1.131559E-05);
	    tmpVector.push_back( 9.181999E-08);
	    tmpVector.push_back(-0.4020348270);
	    tmpVector.push_back(-0.0162463989);
	    tmpVector.push_back( 0.0001623301);
	    tmpVector.push_back( 4.367343E-06);
	    tmpVector.push_back( 1.199000E-08);
	    tmpVector.push_back(-0.0025204776);
	    tmpVector.push_back( 0.0001101514);
	    tmpVector.push_back(-2.320217E-07);
	    tmpVector.push_back( 7.794999E-08);
	    tmpVector.push_back( 9.937483E-06);
	    tmpVector.push_back(-1.346886E-06);
	    tmpVector.push_back( 4.141999E-08);
        cRho.clear();
        cRho = makeMatrix(tmpVector);

        tmpVector.clear();
        tmpVector.push_back( 3822.9712300);
        tmpVector.push_back(-23.122409500);
        tmpVector.push_back( 0.0678775826);
        tmpVector.push_back( 0.0022413893);
        tmpVector.push_back(-0.0003045332);
        tmpVector.push_back(-4.758000E-06);
        tmpVector.push_back( 2.3501449500);
        tmpVector.push_back( 0.1788839410);
        tmpVector.push_back( 0.0006828000);
        tmpVector.push_back( 0.0002101166);
        tmpVector.push_back(-9.812000E-06);
        tmpVector.push_back(-0.0004724176);
        tmpVector.push_back(-0.0003317949);
        tmpVector.push_back( 0.0001002032);
        tmpVector.push_back(-5.306000E-06);
        tmpVector.push_back( 4.242194E-05);
        tmpVector.push_back( 2.347190E-05);
        tmpVector.push_back(-1.894000E-06);
		cHeat.clear();
		cHeat = makeMatrix(tmpVector);

		tmpVector.clear();
		tmpVector.push_back( 0.4082066700);
		tmpVector.push_back(-0.0039816870);
		tmpVector.push_back( 1.583368E-05);
		tmpVector.push_back(-3.552049E-07);
		tmpVector.push_back(-9.884176E-10);
		tmpVector.push_back( 4.460000E-10);
		tmpVector.push_back( 0.0006629321);
		tmpVector.push_back(-2.686475E-05);
		tmpVector.push_back( 9.039150E-07);
		tmpVector.push_back(-2.128257E-08);
		tmpVector.push_back(-5.562000E-10);
		tmpVector.push_back( 3.685975E-07);
		tmpVector.push_back( 7.188416E-08);
		tmpVector.push_back(-1.041773E-08);
		tmpVector.push_back( 2.278001E-10);
		tmpVector.push_back( 4.703395E-08);
		tmpVector.push_back( 7.612361E-11);
		tmpVector.push_back(-2.734000E-10);
		cCond.clear();
		cCond = makeMatrix(tmpVector);

		tmpVector.clear();
		tmpVector.push_back( 1.4725525500);
		tmpVector.push_back( 0.0022218998);
		tmpVector.push_back(-0.0004406139);
		tmpVector.push_back( 6.047984E-06);
		tmpVector.push_back(-1.954730E-07);
		tmpVector.push_back(-2.372000E-09);
		tmpVector.push_back(-0.0411841566);
		tmpVector.push_back( 0.0001784479);
		tmpVector.push_back(-3.564413E-06);
		tmpVector.push_back( 4.064671E-08);
		tmpVector.push_back( 1.915000E-08);
		tmpVector.push_back( 0.0002572862);
		tmpVector.push_back(-9.226343E-07);
		tmpVector.push_back(-2.178577E-08);
		tmpVector.push_back(-9.529999E-10);
		tmpVector.push_back(-1.699844E-06);
		tmpVector.push_back(-1.023552E-07);
		tmpVector.push_back( 4.482000E-09);
		cVisc.clear();
		cVisc = makeMatrix(tmpVector);

		cTfreeze.clear();
		cTfreeze.push_back( 27.755555600); // reference concentration in per cent
		cTfreeze.push_back(-22.973221700);
		cTfreeze.push_back(-1.1040507200);
		cTfreeze.push_back(-0.0120762281);
		cTfreeze.push_back(-9.343458E-05);

    };

	/// Define freezing point calculations
	double Tfreeze(double p, double x){
		IncompressibleClass::checkCoefficients(cTfreeze,5);
		std::vector<double> tmpVector(cTfreeze.begin()+1,cTfreeze.end());
		return polyval(tmpVector, x*100.0-cTfreeze[0])+273.15;
	}
};

#endif
