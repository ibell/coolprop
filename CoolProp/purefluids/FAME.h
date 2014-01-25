#ifndef FAME_H
#define FAME_H

class MethylPalmitateClass : public Fluid {

public:
    MethylPalmitateClass();
    ~MethylPalmitateClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};
class MethylStearateClass : public Fluid {

public:
    MethylStearateClass();
    ~MethylStearateClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};
class MethylOleateClass : public Fluid {

public:
    MethylOleateClass();
    ~MethylOleateClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};
class MethylLinoleateClass : public Fluid {

public:
    MethylLinoleateClass();
    ~MethylLinoleateClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};
class MethylLinolenateClass : public Fluid {

public:
    MethylLinolenateClass();
    ~MethylLinolenateClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

#endif
