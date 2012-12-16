#ifndef BUTENES_H
#define BUTENES_H

class OneButeneClass : public Fluid {

public:
    OneButeneClass();
    ~OneButeneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
};

class IsoButeneClass : public Fluid {

public:
    IsoButeneClass();
    ~IsoButeneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 400; *sigma = 0.64947;};
};

class Cis2ButeneClass : public Fluid {

public:
    Cis2ButeneClass();
    ~Cis2ButeneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 399.3; *sigma = 0.5949;};
};

class Trans2ButeneClass : public Fluid {

public:
    Trans2ButeneClass();
    ~Trans2ButeneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 452.09; *sigma = 0.63617;};
};

#endif
