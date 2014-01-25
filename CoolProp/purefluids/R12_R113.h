#ifndef R12_R113_H
#define R12_R113_H

class R12Class : public Fluid {

public:
    R12Class();
    ~R12Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T)
	{
		// Mulero, JPCRD, 2012
		return -0.000124*pow(1-T/reduce.T,0.4318)+0.05662*pow(1-T/reduce.T,1.263);
	}
};

class R113Class : public Fluid {

public:
    R113Class();
    ~R113Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T)
	{
		// Mulero, JPCRD, 2012
		return 0.0556*pow(1-T/reduce.T,1.24);
	}
};

#endif
