#ifndef R22_H
#define R22_H

	class R22Class : public Fluid{

	public:
		R22Class();
		~R22Class(){};
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);

		void ECSParams(double *e_k, double *sigma);
		double ECS_chi_conductivity(double rhor);
		double ECS_f_int(double T);
		double ECS_psi_viscosity(double rhor);
		double surface_tension_T(double T)
		{
			// From Mulero, 2012, JPCRD
			return 3.0587*pow(1-T/reduce.T,1.41809)+-2.99856*pow(1-T/reduce.T,1.42291);
		};
	};
#endif
