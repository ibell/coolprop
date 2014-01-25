
#ifndef ARGON_H
#define ARGON_H
	#include "FluidClass.h"

	class ArgonClass : public Fluid{

	private:
		double X_tilde(double T,double tau,double delta);
	public:
		ArgonClass();
		~ArgonClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		double surface_tension_T(double T);
	};
#endif
