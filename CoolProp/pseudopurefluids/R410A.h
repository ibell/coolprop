
#ifndef R410_H
#define R410_H

	class R410AClass : public Fluid{

	public:
		R410AClass();
		~R410AClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psatL(double);
		double psatV(double);
		double rhosatL(double);
		double rhosatV(double);
	};
#endif
