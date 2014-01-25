#include "FluidClass.h"
#ifndef R407C_H
#define R407C_H

	class R407CClass : public Fluid{

	public:
		R407CClass();
		~R407CClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psatL(double);
		double psatV(double);
		double rhosatL(double);
		double rhosatV(double);
		double surface_tension_T(double T);
	};

#endif
