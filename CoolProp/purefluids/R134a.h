#ifndef R134a_H
#define R134a_H

	class R134aClass : public Fluid{

	public:
		R134aClass();
		~R134aClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
	};
#endif
