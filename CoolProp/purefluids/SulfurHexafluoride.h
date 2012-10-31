#ifndef SF6_H
#define SF6_H

class SulfurHexafluorideClass : public Fluid {

public:
    SulfurHexafluorideClass();
    ~SulfurHexafluorideClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};


#endif
