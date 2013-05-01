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
		// Chichester
		*e_k = 275.80;
		*sigma = 0.5501;
	}
};

#endif
