#ifndef CYCLOPENTANE_H
#define CYCLOPENTANE_H

class CyclopentaneClass : public Fluid {

public:
    CyclopentaneClass();
    ~CyclopentaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

#endif
