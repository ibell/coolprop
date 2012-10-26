#ifndef ALKANES_H
#define ALKANES_H

class MethaneClass : public Fluid {

public:
    MethaneClass();
    ~MethaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class EthaneClass : public Fluid {

public:
    EthaneClass();
    ~EthaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class nButaneClass : public Fluid {

public:
    nButaneClass();
    ~nButaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class IsoButaneClass : public Fluid {

public:
    IsoButaneClass();
    ~IsoButaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};


#endif
