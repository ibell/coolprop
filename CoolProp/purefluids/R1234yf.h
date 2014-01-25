#ifndef R1234yf_H
#define R1234yf_H

	class R1234yfClass : public Fluid{

	public:
		R1234yfClass();
		~R1234yfClass(){};
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		double conductivity_Trho(double T, double rho);
		double surface_tension_T(double T);
	};
#endif
