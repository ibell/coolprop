#ifndef R744_H
#define R744_H

	class R744Class : public Fluid{

	public:
		R744Class();
		~R744Class(){};
		double conductivity_Trho(double, double);
		double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		double surface_tension_T(double T)
		{
			// From Mulero, 2012, JPCRD
			return 0.07863*pow(1-T/reduce.T,1.254);
		};
	};

#endif
