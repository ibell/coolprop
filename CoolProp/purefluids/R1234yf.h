#ifndef R1234yf_H
#define R1234yf_H

	class R1234yfClass : public Fluid{

	public:
		R1234yfClass();
		~R1234yfClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		double surface_tension_T(double T);
	};
#endif
