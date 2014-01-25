#ifndef UNDECANE_H
#define UNDECANE_H

class UndecaneClass : public Fluid {

public:
    UndecaneClass();
    ~UndecaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

#endif
