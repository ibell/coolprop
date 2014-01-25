#ifndef PROPYLENE_H
#define PROPYLENE_H

class PropyleneClass : public Fluid {

public:
    PropyleneClass();
    ~PropyleneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T);
	void ECSParams(double *e_k, double *sigma)
	{
		// From Huber (2003)
		*e_k = 298.9; *sigma = 0.4678;
	};
	double ECS_f_int(double T)
	{
		// From Huber (2003)
		return 1.09939e-3+3.72538e-7*T;
	};
	double ECS_psi_viscosity(double rhor)
	{
		// From Huber (2003)
		return 1.33962-0.256307*rhor+4.68211e-2*rhor*rhor;
	};
	double ECS_chi_conductivity(double rhor)
	{
		// From Huber (2003)
		return 1.3521-0.123177*rhor;
	};
};

#endif
