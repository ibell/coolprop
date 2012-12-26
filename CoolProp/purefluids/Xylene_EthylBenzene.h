#ifndef XYLENE_ETHYLBENZENE_H
#define XYLENE_ETHYLBENZENE_H

class oXyleneClass : public Fluid {

public:
    oXyleneClass();
    ~oXyleneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
};

class mXyleneClass : public Fluid {

public:
    mXyleneClass();
    ~mXyleneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
};

class pXyleneClass : public Fluid {

public:
    pXyleneClass();
    ~pXyleneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
};

class EthylBenzeneClass : public Fluid {

public:
    EthylBenzeneClass();
    ~EthylBenzeneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
};

#endif
