#ifndef R1233ZDE_H
#define R1233ZDE_H

class R1233zdEClass : public Fluid {

public:
    R1233zdEClass();
    ~R1233zdEClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

#endif
