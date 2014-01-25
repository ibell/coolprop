#ifndef ALKANES_H
#define ALKANES_H

class MethaneClass : public Fluid {

public:
    MethaneClass();
    ~MethaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double viscosity_dilute(double T);
	double viscosity_residual(double T, double rho);
	double viscosity_Trho(double T, double rho);
	double conductivity_dilute(double T);
	double conductivity_residual(double T, double rho);
	double conductivity_Trho(double T, double rho);
	void ECSParams(double *e_k, double *sigma)
	{
		*e_k = 174;
		*sigma = 0.36652;
	}
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.03825*pow(1-T/reduce.T,1.191)+-0.006024*pow(1-T/reduce.T,5.422)-0.0007065*pow(1-T/reduce.T,0.6161);
	}
};

class EthaneClass : public Fluid {

public:
    EthaneClass();
    ~EthaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double viscosity_dilute(double T);
	double viscosity_residual(double T, double rho);
	double viscosity_Trho(double T, double rho);
	double conductivity_dilute(double T);
	double conductivity_residual(double T, double rho);
	double conductivity_Trho(double T, double rho);
	void ECSParams(double *e_k, double *sigma)
	{
		*e_k = 245.0;
		*sigma = 0.43682;
	}
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.07602*pow(1-T/reduce.T,1.32)+-0.02912*pow(1-T/reduce.T,1.676);
	}
};

class nButaneClass : public Fluid {

public:
    nButaneClass();
    ~nButaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){
		// From Vogel HTHP 1999
		*e_k = 280.51; *sigma = 0.57335;
	};
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.05138*pow(1-T/reduce.T,1.209);
	}
	double viscosity_Trho(double, double);
	double conductivity_Trho(double, double);
};

class IsoButaneClass : public Fluid {

public:
    IsoButaneClass();
    ~IsoButaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma);
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return -0.01639*pow(1-T/reduce.T,2.102)+0.06121*pow(1-T/reduce.T,1.304);
	}
	double viscosity_Trho(double, double);
	double conductivity_Trho(double, double);
};


#endif
