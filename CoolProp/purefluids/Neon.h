
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
		double surface_tension_T(double T)
		{
			// From Mulero, 2012, JPCRD
			return 0.012254*pow(1-T/reduce.T,1.4136)+ 0.02728*pow(1-T/reduce.T,1.4517) - 0.025715*pow(1-T/reduce.T,1.6567);
		}
	};
#endif
