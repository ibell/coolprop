#ifndef RC318_R21_R114_R13_R14_H
#define RC318_R21_R114_R13_R14_H

class RC318Class : public Fluid {

public:
    RC318Class();
    ~RC318Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T)
	{
		// Mulero, JPCRD, 2012
		return 0.0507*pow(1-T/reduce.T,1.25);
	}
};

class R21Class : public Fluid {

public:
    R21Class();
    ~R21Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T)
	{
		// Mulero, JPCRD, 2012
		return 0.06924*pow(1-T/reduce.T,1.259);
	}
};

class R114Class : public Fluid {

public:
    R114Class();
    ~R114Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T)
	{
		// Mulero, JPCRD, 2012
		return 0.05239*pow(1-T/reduce.T,1.258);
	}
};

class R13Class : public Fluid {

public:
    R13Class();
    ~R13Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T)
	{
		// Mulero, JPCRD, 2012
		return 0.05045*pow(1-T/reduce.T,1.269);
	}
};

class R14Class : public Fluid {

public:
    R14Class();
    ~R14Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T)
	{
		// Mulero, JPCRD, 2012
		return 0.0423*pow(1-T/reduce.T,1.24);
	}
};

#endif
