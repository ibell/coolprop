#ifndef R236FA_H
#define R236FA_H

class R236FAClass : public Fluid {

public:
    R236FAClass();
    ~R236FAClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){// From Chichester NISTIR report 6650
		*e_k = 307.24; *sigma = 0.5644;};
	//double surface_tension_T(double T);
};

#endif
