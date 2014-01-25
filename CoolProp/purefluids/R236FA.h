#ifndef R236FA_H
#define R236FA_H

class R236FAClass : public Fluid {

public:
    R236FAClass();
    ~R236FAClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T){
		// From Mulero, 2012, JPCRD
		return 0.05389*pow(1-T/reduce.T,1.249);
	};
	void ECSParams(double *e_k, double *sigma)
	{
		// From Huber (2003)
		*e_k = 307.24; *sigma = 0.5644;
	};
	double ECS_f_int(double T)
	{
		// From Huber (2003)
		return 1.00946e-3+1.21255e-6*T;
	};
	double ECS_psi_viscosity(double rhor)
	{
		// From Huber (2003)
		return 1.10195-2.94253e-2*rhor;
	};
	double ECS_chi_conductivity(double rhor)
	{
		// From Huber (2003)
		return 1.1627-4.32746e-2*rhor;
	};
};

#endif
