#ifndef BENZENE_H
#define BENZENE_H

class BenzeneClass : public Fluid {

public:
    BenzeneClass();
    ~BenzeneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double conductivity_Trho(double T, double rho);
	double surface_tension_T(double T){
		// From Mulero, 2012, JPCRD
		return 0.07298*pow(1-T/reduce.T,1.232)-0.0007802*pow(1-T/reduce.T,0.8635)-0.0001756*pow(1-T/reduce.T,0.3065);
	};
	void ECSParams(double *e_k, double *sigma)
	{
		// From Poling (2001)
		*e_k = 412.3; *sigma = 0.5349;
	}
};

#endif
