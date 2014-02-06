#ifndef HFE143M_H
#define HFE143M_H

class HFE143mClass : public Fluid {

public:
    HFE143mClass();
    ~HFE143mClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

#endif

