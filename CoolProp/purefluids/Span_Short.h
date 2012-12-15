#ifndef SPAN_SHORT_H
#define SPAN_SHORT_H

class nPentaneClass : public Fluid {

public:
    nPentaneClass();
    ~nPentaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
};

class nHexaneClass : public Fluid {

public:
    nHexaneClass();
    ~nHexaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){*e_k = 399.3; *sigma = 0.5949;};
};

class nHeptaneClass : public Fluid {

public:
    nHeptaneClass();
    ~nHeptaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){*e_k = 400; *sigma = 0.64947;};
};

class nOctaneClass : public Fluid {

public:
    nOctaneClass();
    ~nOctaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){*e_k = 452.09; *sigma = 0.63617;};
};

class CyclohexaneClass : public Fluid {

public:
    CyclohexaneClass();
    ~CyclohexaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){*e_k = 297.1; *sigma = 0.6182;};
};


#endif
