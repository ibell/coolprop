
#ifndef R507A_H
#define R507A_H

	class R507AClass : public Fluid{

	public:
		R507AClass();
		~R507AClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psatL(double);
		double psatV(double);
		double rhosatL(double);
		double rhosatV(double);
	};

#endif
