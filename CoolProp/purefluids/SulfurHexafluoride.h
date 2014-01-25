#ifndef SF6_H
#define SF6_H

class SulfurHexafluorideClass : public Fluid {

public:
    SulfurHexafluorideClass();
    ~SulfurHexafluorideClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double viscosity_Trho(double T, double rho);
	double conductivity_Trho(double T, double rho);
	void ECSParams(double *e_k, double *sigma)
	{
		// QuinonesCisneros-JPCRD-2012
		*e_k = 215;
		*sigma = 0.5205;
	}
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.0538*pow(1-T/reduce.T,1.271)+-0.00004064*pow(1-T/reduce.T,0.2116);
	};
};


#endif
