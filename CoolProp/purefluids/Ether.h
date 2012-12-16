#ifndef ETHER_H
#define ETHER_H

class DimethylEtherClass : public Fluid {

public:
    DimethylEtherClass();
    ~DimethylEtherClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
};

#endif
