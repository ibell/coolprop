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
	void ECSParams(double *e_k, double *sigma){*e_k = 280.51; *sigma = 0.573;}; // From Chichester NISTIR 6650}
	//double viscosity_Trho(double, double);
	//double conductivity_Trho(double, double);
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
