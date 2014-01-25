#ifndef R236EA_H
#define R236EA_H

class R236EAClass : public Fluid {

public:
    R236EAClass();
    ~R236EAClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T){
		// From Mulero, 2012, JPCRD
		return 0.306974*pow(1-T/reduce.T,1.12614)+-0.247277*pow(1-T/reduce.T,1.09899);
	};
	void ECSParams(double *e_k, double *sigma)
	{
		// From Huber (2003)
		*e_k = 318.33; *sigma = 0.5604;
	}
	double ECS_f_int(double T)
	{
		// From Huber (2003)
		return 1.70267e-3-4.91063e-7*T;
	}
	double ECS_psi_viscosity(double rhor)
	{
		// From Huber (2003)
		return 1.12216-2.73101e-2*rhor;
	}
	double ECS_chi_conductivity(double rhor)
	{
		// From Huber (2003)
		return 0.9617+3.37897e-2*rhor;
	}
};

#endif
