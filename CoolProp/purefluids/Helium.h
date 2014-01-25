#include "FluidClass.h"
#ifndef HELIUM_H
#define HELIUM_H

	class HeliumClass : public Fluid{

	private:
		double X_tilde(double T,double tau,double delta);
	public:
		HeliumClass();
		~HeliumClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		double surface_tension_T(double T);
	};
#endif
