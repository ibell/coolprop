#ifndef R124_H
#define R124_H

class R124Class : public Fluid {

public:
    R124Class();
    ~R124Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T)
	{
		// Mulero, JPCRD 2012
		return 0.05175*pow(1-T/reduce.T,1.197);
	}
	void ECSParams(double *e_k, double *sigma)
	{
		// From Huber (2003)
		*e_k = 275.80;
		*sigma = 0.5501;
	}
	double ECS_psi_viscosity(double rhor)
	{
		// From Huber (2003)
		return 1.04253+1.38528e-3*rhor;
	}
	double ECS_f_int(double T)
	{
		// From Huber (2003)
		return 1.17690e-3+6.78397e-7*T;
	}
	double ECS_chi_conductivity(double rhor)
	{
		// From Huber (2003)
		return 1.0898-1.54229e-2*rhor;
	}
};

#endif
