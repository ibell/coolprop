#ifndef R143A_H
#define R143A_H

class R143AClass : public Fluid {

public:
    R143AClass();
    ~R143AClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};


#endif
