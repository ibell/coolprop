#ifndef R32_H
#define R32_H

	class R32Class : public Fluid{

	public:
		R32Class();
		~R32Class(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
	};
#endif
