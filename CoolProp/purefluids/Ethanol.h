#ifndef ETHANOL_H
#define ETHANOL_H

class EthanolClass : public Fluid {

public:
    EthanolClass();
    ~EthanolClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

};

#endif
