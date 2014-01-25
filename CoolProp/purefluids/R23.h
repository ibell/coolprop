#ifndef R23_H
#define R23_H

class R23Class : public Fluid {

public:
    R23Class();
    ~R23Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double viscosity_Trho(double, double);
	double conductivity_Trho(double, double);
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return -0.32359*pow(1-T/reduce.T,1.6055)+0.37702*pow(1-T/reduce.T,1.5232);
	}
	void ECSParams(double *e_k, double *sigma){
		// Chichester
		*e_k = 243.91; *sigma = 0.4278;
	};
};

#endif
