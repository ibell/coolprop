#ifndef ETHYLENE_H
#define ETHYLENE_H

class EthyleneClass : public Fluid {

public:
    EthyleneClass();
    ~EthyleneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma)
	{
		// Poling
		*e_k = 224.7;
		*sigma = 0.4163;
	}
	double surface_tension_T(double T)
	{
		// from Mulero, 2012, JPCRD
		return 0.0477*pow(1-T/reduce.T,1.17);
	}
};


#endif
