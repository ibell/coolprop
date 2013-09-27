
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
/** Employs the base functions implemented in IncompBase.h.
 *  Extends the functions for composition as input. */
class IncompressibleSolution : public IncompressibleClass{

protected:
	double xmin, xmax;

//	/// Function used to enforce the composition as parameter
//	/** Overwrites the normal functions that only take 2
//	 *  parameters (T,p) and throws an exception. */
//	bool needComposition() {
//		throw NotImplementedError(format("The fluid %s needs an additional input for the composition.",this->name));
//		return false;
//	}

public:
	// Constructor
	IncompressibleSolution(){
		xmin = -1.;
		xmax = -1.;
	};

	// Destructor, no implementation
	virtual ~IncompressibleSolution(){};

	/* All functions need T, p and x as input. Might not
	 * be necessary, but gives a clearer structure.
	 */
    virtual double rho (double T_K, double p, double x){return -_HUGE;};
    virtual double cp  (double T_K, double p, double x){return -_HUGE;};
    virtual double h   (double T_K, double p, double x){return -_HUGE;};
    virtual double s   (double T_K, double p, double x){return -_HUGE;};
    virtual double visc(double T_K, double p, double x){return -_HUGE;};
    virtual double cond(double T_K, double p, double x){return -_HUGE;};
    virtual double u   (double T_K, double p, double x){return -_HUGE;};
    virtual double psat(double T_K          , double x){return -_HUGE;};
    virtual double Tfreeze(         double p, double x){return -_HUGE;};

    void testInputs(double T_K, double p, double x);

protected:
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
	 *  defined for all fluids, default values do not
	 *  cause errors. */
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

bool IsIncompressibleSolution(std::string name);
double IncompSolution(long iOutput, double T, double p, double x, long iFluid);
double IncompSolution(long iOutput, double T, double p, double x, std::string name);


/// Class to use the SecCool parameters
/** Employs some basic wrapper-like functionality
 *  to bridge the gap between the solution functions
 *  used in CoolProp and the definition used in
 *  SecCool. Please visit:
 *  http://en.ipu.dk/Indhold/refrigeration-and-energy-technology/seccool.aspx
 *  Many thanks to Morten Juel Skovrup for providing
 *  this nice piece of software as well as the parameters
 *  needed to calculate the composition based properties. */
class SecCoolSolution : public IncompressibleSolution{

protected:
	double Tbase;
	double xbase;
	std::vector< std::vector<double> > cRho;
	std::vector< std::vector<double> > cHeat;
	std::vector< std::vector<double> > cVisc;
	std::vector< std::vector<double> > cCond;
	std::vector< std::vector<double> > cPsat;
	std::vector<double> cTfreeze;

public:
	double getTbase() const {return Tbase;}
	double getxbase() const {return xbase;}
	std::vector<std::vector<double> > getcCond() const {return cCond;}
	std::vector<std::vector<double> > getcHeat() const {return cHeat;}
	std::vector<std::vector<double> > getcPsat() const {return cPsat;}
	std::vector<std::vector<double> > getcRho() const {return cRho;}
	std::vector<double> getcTfreeze() const {return cTfreeze;}
	std::vector<std::vector<double> > getcVisc() const {return cVisc;}

public:
	// Constructor
	SecCoolSolution(){
		Tbase = -1.;
		xbase = -1.;
	};
	// Destructor, no implementation
	~SecCoolSolution(){};

public:
	double baseFunction(std::vector<double> coefficients, double T_K, double p, double x);

	std::vector< std::vector<double> > makeMatrix(std::vector<double> coefficients);

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
	double cp(double T_K, double p, double x){
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
		return fracint(cHeat, getxInput(x), getTInput(T_K), getTInput(Tref));
	}
	double visc(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cVisc,6,4);
		return expval(cVisc, getxInput(x), getTInput(T_K), 2)/1e5;
	}
	double cond(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cCond,6,4);
		return polyval(cCond, getxInput(x), getTInput(T_K))/1e3;
	}
	double u(double T_K, double p, double x){
		checkTPX(T_K, p, x);
		IncompressibleClass::checkCoefficients(cHeat,6,4);
		return polyint(cHeat, getxInput(x), getTInput(T_K), getTInput(Tref));
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
	double Tfreeze(double p, double x){
		IncompressibleClass::checkCoefficients(cTfreeze,5);
		std::vector<double> tmpVector(cTfreeze.begin()+1,cTfreeze.end());
		return polyval(tmpVector, x*100.0-cTfreeze[0])+273.15;
	}
};

class MethanolSolution : public SecCoolSolution{
public:
	MethanolSolution(){

		std::vector<double> tmpVector;

        name = std::string("Test");
        description = std::string("Test Methanol");
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
		cTfreeze.push_back (27.755555600); // concentration
		cTfreeze.push_back(-22.973221700);
		cTfreeze.push_back(-1.1040507200);
		cTfreeze.push_back(-0.0120762281);
		cTfreeze.push_back(-9.343458E-05);

    };
};

//
///** Handle all the objects in a single list of incompressible liquids
// *  and a list of solutions / brines. The singleton pattern assures
// *  there is only one list.
// */
//
//class SolutionsContainer {
//public:
//  static SolutionsContainer* Instance();
//  IncompressibleSolution * get_solution(long index);
//  IncompressibleSolution * get_solution(std::string name);
//  void set_solutions(std::vector<IncompressibleSolution*> list);
//
//private:
//  std::vector<IncompressibleSolution*> solution_list;
//  std::map<std::string,IncompressibleSolution*> solution_map;
//
//  static void destruct();
//
//  SolutionsContainer();
//  ~SolutionsContainer(){};
//
//  SolutionsContainer(SolutionsContainer const&);            // no copying
//  SolutionsContainer& operator=(SolutionsContainer const&); // no overwriting
//  static SolutionsContainer* MInstance;                     // instance
//};

#endif
