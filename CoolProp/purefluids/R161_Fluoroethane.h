#ifndef R161EA_H
#define R161EA_H

class R161Class : public Fluid {

public:
    R161Class();
    ~R161Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T)
	{
		// Mulero, 2012
		return 0.05385*pow(1-T/reduce.T,1.111);
	}
};

#endif
