#ifndef DMC_H
#define DMC_H

class DimethylCarbonateClass : public Fluid {

public:
    DimethylCarbonateClass();
    ~DimethylCarbonateClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
};

#endif
