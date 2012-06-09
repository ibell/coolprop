#include "FluidClass.h"
#ifndef R717_H
#define R717_H

	class R717Class : public Fluid{

	public:
		R717Class();
		~R717Class(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		double surface_tension_T(double T);
	};

#endif
