#ifndef CYCLOPROPANECLASS_H
#define CYCLOPROPANECLASS_H

class CycloPropaneClass : public Fluid {

public:
    CycloPropaneClass();
    ~CycloPropaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class PropyneClass : public Fluid {

public:
    PropyneClass();
    ~PropyneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T)
	{
		// Mulero, JPCRD, 2012
		return 0.05801*pow(1-T/reduce.T,1.205);
	}
};

#endif
