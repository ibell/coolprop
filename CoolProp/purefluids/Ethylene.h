#ifndef ETHYLENE_H
#define ETHYLENE_H

class EthyleneClass : public Fluid {

public:
    EthyleneClass();
    ~EthyleneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};


#endif
