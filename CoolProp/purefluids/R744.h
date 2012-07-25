#ifndef R744_H
#define R744_H

	class R744Class : public Fluid{

	public:
		R744Class();
		~R744Class(){};
		virtual double conductivity_Trho(double, double);
		double conductivity_critical(double, double);
		virtual double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
	};

#endif
