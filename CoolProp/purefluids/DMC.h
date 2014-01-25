#ifndef DMC_H
#define DMC_H

class DimethylCarbonateClass : public Fluid {

public:
    DimethylCarbonateClass();
    ~DimethylCarbonateClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

#endif
