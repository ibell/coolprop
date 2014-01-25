#ifndef ETHER_H
#define ETHER_H

class DimethylEtherClass : public Fluid {

public:
    DimethylEtherClass();
    ~DimethylEtherClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double viscosity_Trho(double T, double rho);
	double conductivity_Trho(double T, double rho);
	double surface_tension_T(double T);
	void ECSParams(double *e_k, double *sigma);
};

#endif
