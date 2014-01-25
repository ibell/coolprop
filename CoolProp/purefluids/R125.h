#ifndef R125_H
#define R125_H

class R125Class : public Fluid {

public:
    R125Class();
    ~R125Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double viscosity_Trho(double T, double rho);
	double conductivity_Trho(double T, double rho);
	double surface_tension_T(double T);
	void ECSParams(double *e_k, double *sigma);
};

#endif
