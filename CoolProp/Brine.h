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

class PAO : IncompressibleLiquid
{
public:
	PAO(){
		name.assign("PAO");
		description.assign("Poly-alpha-olefin oil");
		reference.assign("Reference for PAO");
	};

	double rho(double T_K){
		return 808.1-0.6872*(T_K-273.15);
	};
	double cp(double T_K){
		return (0.003391*(T_K-273.15)+2.133);
	};
	double cond(double T_K){
		return (-0.0001088)*(T_K-273.15)+0.1466;
	};
	double visc(double T_K){
		return ((0.004166)*exp((-0.09615)*(T_K-273.15)) + (0.01148)*exp((-0.02563)*(T_K-273.15)));
	};
};

double SecFluids(char Output, double T, double p,char * Ref);
int Brine(char * Mix, double T, double C, /*in --- out */double *Tfreeze, double *Tmax, double *rho, double *cp, double *k, double *visc, double *h, double *s);
#endif