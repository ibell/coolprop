
#ifndef NEON_H
#define NEON_H
	#include "FluidClass.h"

	class NeonClass : public Fluid{

	public:
		NeonClass();
		~NeonClass(){};
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		void ECSParams(double *e_k, double *sigma)
		{
			// From Poling, 2001
			*e_k = 32.8;
			*sigma = 0.282;
		}
		double surface_tension_T(double T);
	};
#endif
