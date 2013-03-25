#ifndef R236EA_H
#define R236EA_H

class R236EAClass : public Fluid {

public:
    R236EAClass();
    ~R236EAClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){
		// From Chichester NISTIR report 6650
		*e_k = 318.33; *sigma = 0.5604;};
	//double surface_tension_T(double T);
};

#endif
