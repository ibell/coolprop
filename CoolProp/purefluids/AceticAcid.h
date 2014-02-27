#ifndef ACETICACID_H
#define ACETICACID_H

class AceticAcidClass : public Fluid {

public:
    AceticAcidClass();
    ~AceticAcidClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

#endif
