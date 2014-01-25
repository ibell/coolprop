#ifndef BUTENES_H
#define BUTENES_H

class OneButeneClass : public Fluid {

public:
    OneButeneClass();
    ~OneButeneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.05644*pow(1-T/reduce.T,1.248);
	};
};

class IsoButeneClass : public Fluid {

public:
    IsoButeneClass();
    ~IsoButeneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.0545*pow(1-T/reduce.T,1.23);
	};
};

class Cis2ButeneClass : public Fluid {

public:
    Cis2ButeneClass();
    ~Cis2ButeneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class Trans2ButeneClass : public Fluid {

public:
    Trans2ButeneClass();
    ~Trans2ButeneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

#endif
