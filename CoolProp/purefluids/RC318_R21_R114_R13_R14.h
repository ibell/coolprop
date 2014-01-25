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
	void ECSParams(double *e_k, double *sigma)
	{
		// From Huber (2003)
		*e_k = 299.76; *sigma = 0.5947;
	}
	double ECS_f_int(double T)
	{
		// From Huber (2003)
		return 1.35697e-3-1.11635e-7*T;
	}
	double ECS_psi_viscosity(double rhor)
	{
		// From Huber (2003)
		return 1.21141-3.37573e-2*rhor;
	}
	double ECS_chi_conductivity(double rhor)
	{
		// From Huber (2003)
		return 1.5249-0.147564*rhor;
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
	void ECSParams(double *e_k, double *sigma)
	{
		// From Huber (2003)
		*e_k = 204.00; *sigma = 0.4971;
	}
	double ECS_f_int(double T)
	{
		// From Huber (2003)
		return 1.07447e-3-6.42373e-7*T;
	}
	double ECS_psi_viscosity(double rhor)
	{
		// From Huber (2003)
		return 0.97618+1.48047e-2*rhor;
	}
	double ECS_chi_conductivity(double rhor)
	{
		// From Huber (2003)
		return 1.1394-3.65562e-2*rhor;
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
	void ECSParams(double *e_k, double *sigma)
	{
		// From Huber (2003)
		*e_k = 164.44; *sigma = 0.4543;
	}
	double ECS_psi_viscosity(double rhor)
	{
		// From Huber (2003)
		return 1.10941-0.630268e-1*rhor;
	}
	double ECS_f_int(double T)
	{
		// From Huber (2003)
		return 1.19864e-3-1.90048e-7*T;
	}
	double ECS_chi_conductivity(double rhor)
	{
		// From Huber (2003)
		return 1.0442;
	}
};

#endif
