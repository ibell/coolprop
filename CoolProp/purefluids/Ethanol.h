#ifndef ETHANOL_H
#define ETHANOL_H

class EthanolClass : public Fluid {

public:
    EthanolClass();
    ~EthanolClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double);
};

#endif
