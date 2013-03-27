#ifndef R143A_H
#define R143A_H

class R143AClass : public Fluid {

public:
    R143AClass();
    ~R143AClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.05416*pow(1-T/reduce.T,1.255);
	};
};


#endif
