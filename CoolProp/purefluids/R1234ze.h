#ifndef R1234ZE_H
#define R1234ZE_H

	class R1234zeClass : public Fluid{

	public:
		R1234zeClass();
		~R1234zeClass(){};
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		double conductivity_Trho(double T, double rho);
	};

#endif
