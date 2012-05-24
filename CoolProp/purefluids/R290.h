#ifndef R290_H
#define R290_H

	class R290Class : public Fluid{

	public:
		R290Class();
		~R290Class(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
	};

#endif
