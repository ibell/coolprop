
#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include "Units.h"
#include "math.h"
#include "stdio.h"
#include <string.h>
#include "Brine.h"
#include "CoolPropTools.h"
#include "CoolProp.h"
#include "IncompBase.h"
#include "IncompSolution.h"


/// Base class for simplified brine/solution models
/** Employs the base functions implemented in IncompBase.h and
 *  provides properties as function of temperature, pressure
 *  and composition. */
void IncompressibleSolution::testInputs(double T_K, double p, double x){
	double result = 0.;

	printf(" %s \n"," ");
	printf("Testing  %s \n",this->get_name().c_str());
	printf("Inputs:  T = %3.3f degC \t p = %2.4f bar \t x = %1.5f \n",T_K-273.15,p/1e5,x);

	result = this->rho(T_K,p,x);
	printf("From object:    rho = %4.2f \t kg/m3    \n",result);
	result = this->cp(T_K,p,x);
	printf("From object:     cp = %1.5f \t kJ/kg-K  \n",result/1e3);
	result = this->h(T_K,p,x);
	printf("From object:      h = %3.3f \t kJ/kg    \n",result/1e3);
	result = this->s(T_K,p,x);
	printf("From object:      s = %1.5f \t kJ/kg-K  \n",result/1e3);
	result = this->visc(T_K,p,x);
	printf("From object:    eta = %1.5f \t 1e-5 Pa-s\n",result*1e5);
	result = this->cond(T_K,p,x);
	printf("From object: lambda = %1.5f \t W/m-k    \n",result*1e3);
	result = this->u(T_K,p,x);
	printf("From object:      u = %3.3f \t kJ/kg    \n",result/1e3);
	result = this->psat(T_K,x);
	printf("From object:   psat = %2.4f \t bar      \n",result/1e5);
	result = this->Tfreeze(p,x);
	printf("From object:Tfreeze = %3.3f \t degC   \n\n",result-273.15);
}


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
bool IncompressibleSolution::checkT(double T_K, double p, double x){
	if( Tmin < 0. ) {
		throw ValueError("Please specify the minimum temperature.");
	} else if( Tmax < 0.) {
		throw ValueError("Please specify the maximum temperature.");
	} else if ( (Tmin>T_K) || (T_K>Tmax) ) {
		throw ValueError(format("Your temperature %f is not between %f and %f.",T_K,Tmin,Tmax));
	} else if (T_K < Tfreeze(p,x)) {
		throw ValueError(format("Your temperature %f is below the freezing point of %f.",T_K,Tfreeze(p,x)));
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
bool IncompressibleSolution::checkP(double T_K, double p, double x) {
	double ps = psat(T_K,x);
	if (p<ps) {
		throw ValueError(format("Equations are valid for solution phase only: %f < %f. ",p,ps));
	} else {
		return true;
	}
}

/// Check validity of composition input.
/** Compares the given composition x to a stored minimum and
 *  maximum value. Enforces the redefinition of xmin and
 *  xmax since the default values cause an error. */
bool IncompressibleSolution::checkX(double x){
	if( xmin < 0. ) {
		throw ValueError("Please specify the minimum concentration.");
	} else if( xmax < 0.) {
		throw ValueError("Please specify the maximum concentration.");
	} else if ( (xmin>x) || (x>xmax) ) {
		throw ValueError(format("Your composition %f is not between %f and %f.",x,xmin,xmax));
	} else {
		return true;
	}
	return false;
}

/// Check validity of temperature, pressure and composition input.
bool IncompressibleSolution::checkTPX(double T, double p, double x) {
	return (checkT(T,p,x) && checkP(T,p,x) && checkX(x));
}



/** Handle all the objects in a single list of incompressible
 *  solutions and brines. */
SolutionsContainer::SolutionsContainer() {
	std::vector<IncompressibleSolution*> tmpVector;

	tmpVector.push_back(new MelinderSolution());
	tmpVector.push_back(new EGSolution());
	tmpVector.push_back(new PGSolution());
	tmpVector.push_back(new EASolution());
	tmpVector.push_back(new MASolution());
	tmpVector.push_back(new GLSolution());
	tmpVector.push_back(new AMSolution());
	tmpVector.push_back(new KCSolution());
	tmpVector.push_back(new CASolution());
	tmpVector.push_back(new MGSolution());
	tmpVector.push_back(new NASolution());
	tmpVector.push_back(new KASolution());
	tmpVector.push_back(new KFSolution());
	tmpVector.push_back(new LISolution());

	tmpVector.push_back(new SecCoolSolution());
	tmpVector.push_back(new ZitrecAC());
	tmpVector.push_back(new IceSlurryEA());
	tmpVector.push_back(new IceSlurryPG());
	tmpVector.push_back(new IceSlurryNA());
	tmpVector.push_back(new PK2000());

	tmpVector.push_back(new LiBrSolution());

	// Now we store the vector in the variable
	// and overwrite the map.
	set_solutions(tmpVector);
}

SolutionsContainer::~SolutionsContainer() {
	while (!solution_list.empty()) {
		delete solution_list.back();
		solution_list.pop_back();
	}
}

IncompressibleSolution * SolutionsContainer::get_solution(long index){
	return solution_list[index];
}

IncompressibleSolution * SolutionsContainer::get_solution(std::string name){
	std::map<std::string,IncompressibleSolution*>::iterator it;
	// Try to find using the map if Fluid name is provided
	it = solution_map.find(name);
	// If it is found the iterator will not be equal to end
	if (it != solution_map.end() )
	{
		// Return a pointer to the class
		return (*it).second;
	}
	else{
		return NULL;
	}
}

void SolutionsContainer::set_solutions(std::vector<IncompressibleSolution*> list){
	solution_list = list;
	// Build the map of fluid names mapping to pointers to the solution class instances
	for (std::vector<IncompressibleSolution*>::iterator it = solution_list.begin(); it != solution_list.end(); it++)
	{
		// Load up entry in map
		solution_map[(*it)->get_name()] = *it;
	}
}

SolutionsContainer Solutions = SolutionsContainer();

bool IsIncompressibleSolution(std::string name)
{
	IncompressibleSolution *pSolution = Solutions.get_solution(getSolutionName(name));
	if (pSolution == NULL){ //Not found since NULL pointer returned
		return false;
	}
	else{
		return true;
	}
}

double pIncompSolutionSI(long iOutput, double T, double p_SI, double x, IncompressibleSolution *pSolution)
{
	double out;
	switch (iOutput)
	{
		case iT:
			out = T; break;
		case iP:
			out = p_SI; break;
		case iD:
			out = pSolution->rho(T,p_SI,x); break;
		case iC:
			out = pSolution->cp(T,p_SI,x); break;
		case iS:
			out = pSolution->s(T,p_SI,x); break;
		case iU:
			out = pSolution->u(T,p_SI,x); break;
		case iH:
			out = pSolution->h(T,p_SI,x); break;
		case iV:
			out = pSolution->visc(T,p_SI,x); break;
		case iL:
			out = pSolution->cond(T,p_SI,x); break;
		case iTmax:
			out = pSolution->getTmax(); break;
		case iTmin:
			out = pSolution->getTmin(); break;
		case iTfreeze:
			out = pSolution->Tfreeze(p_SI,x); break;
		case iPsat:
			out = pSolution->psat(T,x); break;
		default:
			throw ValueError(format("Your index [%d] is invalid for the incompressible solution %s",iOutput,pSolution->getName().c_str()));
			out=0; break;
	}
	return out;
}

double IncompSolutionSI(long iOutput, double T, double p, double x, long iFluid)
{
	return pIncompSolutionSI(iOutput,T,p,x,Solutions.get_solution(iFluid));
}

double IncompSolutionSI(long iOutput, double T, double p, double x, std::string name)
{
	return pIncompSolutionSI(iOutput,T,p,x,Solutions.get_solution(name));
}

/** Just some convenience functions to allow for backward compatibility.
 *  We might end up using them in the beginning until we can
 *  handle a third parameter properly. */
double IncompSolutionSI(long iOutput, double T, double p, std::string name){ // TODO Solutions: Remove as soon as possible
	// Split into fluid, concentration pair
	std::vector<std::string> fluid_concentration = strsplit(std::string(name),'-');
	// Check it worked
	if (fluid_concentration.size() != 2){
		throw ValueError(format("Format of incompressible solution string [%s] is invalid, should be like \"EG-20%\" or \"EG-0.2\" ", name.c_str()) );
	}
	// Convert the concentration into a string
	char* pEnd;
	double x = strtod(fluid_concentration[1].c_str(), &pEnd);
	// Check if per cent or fraction syntax is used
	if (!strcmp(pEnd,"%")){	x *= 0.01;}
	return IncompSolutionSI(iOutput, T, p, x, fluid_concentration[0]);
}
std::string getSolutionName(std::string name){ // TODO Solutions: Remove as soon as possible
	// Split into fluid, concentration pair
	std::vector<std::string> fluid_concentration = strsplit(std::string(name),'-');
	return fluid_concentration[0];
}
double getSolutionConc(std::string name){ // TODO Solutions: Remove as soon as possible
	// Split into fluid, concentration pair
	std::vector<std::string> fluid_concentration = strsplit(std::string(name),'-');
	// Check it worked
	if (fluid_concentration.size() != 2){
		throw ValueError(format("Format of incompressible solution string [%s] is invalid, should be like \"EG-20%\" or \"EG-0.2\" ", name.c_str()) );
	}
	// Convert the concentration into a string
	char* pEnd;
	double x = strtod(fluid_concentration[1].c_str(), &pEnd);
	// Check if per cent or fraction syntax is used
	if (!strcmp(pEnd,"%")){	x *= 0.01;}
	return x;
}


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
double BaseSolution::baseFunction(std::vector<double> const& coefficients, double T_K, double p, double x){
	IncompressibleClass::checkCoefficients(coefficients,18);
	return (((((
			 coefficients[17])*x
			+coefficients[16])*x
			+coefficients[15])*T_K
		 +(((coefficients[14])*x
			+coefficients[13])*x
			+coefficients[12])*x
			+coefficients[11])*T_K
		+((((coefficients[10])*x
			+coefficients[9])*x
			+coefficients[8])*x
			+coefficients[7])*x
			+coefficients[6])*T_K
	   +(((((coefficients[5])*x
			+coefficients[4])*x
			+coefficients[3])*x
			+coefficients[2])*x
			+coefficients[1])*x
			+coefficients[0];
}

std::vector< std::vector<double> > BaseSolution::makeMatrix(std::vector<double> const& coefficients){
	IncompressibleClass::checkCoefficients(coefficients,18);
	std::vector< std::vector<double> > matrix;
	std::vector<double> tmpVector;

	tmpVector.clear();
	tmpVector.push_back(coefficients[0]);
	tmpVector.push_back(coefficients[6]);
	tmpVector.push_back(coefficients[11]);
	tmpVector.push_back(coefficients[15]);
	matrix.push_back(tmpVector);

	tmpVector.clear();
	tmpVector.push_back(coefficients[1]);
	tmpVector.push_back(coefficients[7]);
	tmpVector.push_back(coefficients[12]);
	tmpVector.push_back(coefficients[16]);
	matrix.push_back(tmpVector);

	tmpVector.clear();
	tmpVector.push_back(coefficients[2]);
	tmpVector.push_back(coefficients[8]);
	tmpVector.push_back(coefficients[13]);
	tmpVector.push_back(coefficients[17]);
	matrix.push_back(tmpVector);

	tmpVector.clear();
	tmpVector.push_back(coefficients[3]);
	tmpVector.push_back(coefficients[9]);
	tmpVector.push_back(coefficients[14]);
	tmpVector.push_back(0.0);
	matrix.push_back(tmpVector);

	tmpVector.clear();
	tmpVector.push_back(coefficients[4]);
	tmpVector.push_back(coefficients[10]);
	tmpVector.push_back(0.0);
	tmpVector.push_back(0.0);
	matrix.push_back(tmpVector);

	tmpVector.clear();
	tmpVector.push_back(coefficients[5]);
	tmpVector.push_back(0.0);
	tmpVector.push_back(0.0);
	tmpVector.push_back(0.0);
	matrix.push_back(tmpVector);

	tmpVector.clear();
	return matrix;
}

/// Convert pre-v4.0-style coefficient array to new format
std::vector<std::vector<double> > MelinderSolution::convertCoeffs(double* oldestCoeffs, const int A, const int B) {
	const int lengthA = 18;
	const int lengthB =  5;
    if ((A==lengthA) && (B==lengthB) ){
    	/// Rearrange array
    	double oldCoeffs[lengthA][lengthB];
		for (int j = 0; j < lengthB; j++) {
			oldCoeffs[ 0][j] = oldestCoeffs[( 0 * B) + j ];
			oldCoeffs[ 1][j] = oldestCoeffs[( 4 * B) + j ];
			oldCoeffs[ 2][j] = oldestCoeffs[( 8 * B) + j ];
			oldCoeffs[ 3][j] = oldestCoeffs[(12 * B) + j ];
			oldCoeffs[ 4][j] = oldestCoeffs[(15 * B) + j ];
			oldCoeffs[ 5][j] = oldestCoeffs[(17 * B) + j ];
			oldCoeffs[ 6][j] = oldestCoeffs[( 1 * B) + j ];
			oldCoeffs[ 7][j] = oldestCoeffs[( 5 * B) + j ];
			oldCoeffs[ 8][j] = oldestCoeffs[( 9 * B) + j ];
			oldCoeffs[ 9][j] = oldestCoeffs[(13 * B) + j ];
			oldCoeffs[10][j] = oldestCoeffs[(16 * B) + j ];
			oldCoeffs[11][j] = oldestCoeffs[( 2 * B) + j ];
			oldCoeffs[12][j] = oldestCoeffs[( 6 * B) + j ];
			oldCoeffs[13][j] = oldestCoeffs[(10 * B) + j ];
			oldCoeffs[14][j] = oldestCoeffs[(14 * B) + j ];
			oldCoeffs[15][j] = oldestCoeffs[( 3 * B) + j ];
			oldCoeffs[16][j] = oldestCoeffs[( 7 * B) + j ];
			oldCoeffs[17][j] = oldestCoeffs[(11 * B) + j ];
		}
		std::vector<std::vector<double> > tmpVector;
		tmpVector.push_back(std::vector<double>());
		tmpVector.push_back(std::vector<double>());
		tmpVector.push_back(std::vector<double>());
		tmpVector.push_back(std::vector<double>());
		tmpVector.push_back(std::vector<double>());
		for (int i = 0; i < lengthA; i++) {
			for (int j = 0; j < lengthB; j++) {
				if (i<1) tmpVector[j].clear();
				tmpVector[j].push_back(oldCoeffs[i][j]);
			}
		}
		return tmpVector;
    }
    throw ValueError(format("Your array has the dimensions [%d,%d], but only [%d,%d] arrays are supported.",A,B,lengthA,lengthB));
}

double const LiBrSolution::M_H2O  = 0.018015268; /* kg/mol, molar mass of H2O */
double const LiBrSolution::M_LiBr = 0.08685; /* kg/mol, molar mass of LiBr */
double const LiBrSolution::T0     = 221; /* K, constant */

/* Critical point of H2O */
double const LiBrSolution::Tc_H2O   = 647.096; /* K, temperature  */
double const LiBrSolution::pc_H2O   = 22.064; /* MPa, pressure */
double const LiBrSolution::rhoc_H2O = 17873; /* mol/m^3 (322 kg/m^3), molar density */
double const LiBrSolution::hc_H2O   = 37548.5; /* J/mol, molar enthalpy */
double const LiBrSolution::sc_H2O   = 79.3933; /* J/(mol.K) molar entropy*/

/*Triple point of H2O */
double const LiBrSolution::Tt_H2O   = 273.16; /* K, temperature */
double const LiBrSolution::cpt_H2O  = 76.0226; /* J/(mol.K), molar isobaric heat capacity*/