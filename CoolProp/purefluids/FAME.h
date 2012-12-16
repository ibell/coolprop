#ifndef FAME_H
#define FAME_H

class MethylPalmitateClass : public Fluid {

public:
    MethylPalmitateClass();
    ~MethylPalmitateClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
};
class MethylStearateClass : public Fluid {

public:
    MethylStearateClass();
    ~MethylStearateClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
};
class MethylOleateClass : public Fluid {

public:
    MethylOleateClass();
    ~MethylOleateClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
};
class MethylLinoleateClass : public Fluid {

public:
    MethylLinoleateClass();
    ~MethylLinoleateClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
};
class MethylLinolenateClass : public Fluid {

public:
    MethylLinolenateClass();
    ~MethylLinolenateClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
};

#endif
