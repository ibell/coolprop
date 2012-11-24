/*
The name of this file is a slight misnomer.  It includes the analysis
aqueous solutions as well as single-phase liquids.
*/
#include <string>

#ifndef BRINE_H
#define BRINE_H

class IncompressibleLiquid{
	
public:
	std::string name,description,reference;
	
	// Constructor
	IncompressibleLiquid();

	///< Destructor.  No implementation
	~IncompressibleLiquid(){};

	double rho(double T_K);
	double cp(double T_K);
	double h(double T_K, double p);
	double s(double T_K);
	double visc(double T_K);
	double cond(double T_K);
};

double SecFluids(char Output, double T, double p,char * Ref);
int Brine(char * Mix, double T, double C, /*in --- out */double *Tfreeze, double *Tmax, double *rho, double *cp, double *k, double *visc, double *h, double *s);
#endif