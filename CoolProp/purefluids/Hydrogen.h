#ifndef Hydrogen_H
#define Hydrogen_H

	class HydrogenClass : public Fluid{

	public:
		HydrogenClass();
		~HydrogenClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		double surface_tension_T(double);
	};

#endif