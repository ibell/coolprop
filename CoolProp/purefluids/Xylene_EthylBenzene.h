#ifndef XYLENE_ETHYLBENZENE_H
#define XYLENE_ETHYLBENZENE_H

class oXyleneClass : public Fluid {

public:
    oXyleneClass();
    ~oXyleneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class mXyleneClass : public Fluid {

public:
    mXyleneClass();
    ~mXyleneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class pXyleneClass : public Fluid {

public:
    pXyleneClass();
    ~pXyleneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class EthylBenzeneClass : public Fluid {

public:
    EthylBenzeneClass();
    ~EthylBenzeneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

#endif
