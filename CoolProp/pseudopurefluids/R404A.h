
#ifndef R404A_H
#define R404A_H

	class R404AClass : public Fluid{

	public:
		R404AClass();
		~R404AClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psatL(double);
		double psatV(double);
		double rhosatL(double);
		double rhosatV(double);
		double surface_tension_T(double T);
	};
#endif
