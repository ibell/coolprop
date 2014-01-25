#ifndef METHANOL_H
#define METHANOL_H

class MethanolClass : public Fluid {

public:
    MethanolClass();
    ~MethanolClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	/*double viscosity_Trho(double T, double rho);
	double conductivity_Trho(double T, double rho);*/
	//void ECSParams(double *e_k, double *sigma)
	//{
	//	// from Kiselev 2005
	//	*e_k = 362.6;
	//	*sigma = 0.453;
	//};
	double surface_tension_T(double T);
};

#endif
