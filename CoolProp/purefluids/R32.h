#ifndef R32_H
#define R32_H

	class R32Class : public Fluid{

	public:
		R32Class();
		~R32Class(){};
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		void ECSParams(double *e_k, double *sigma);
		double ECS_psi_viscosity(double rhor);
		double ECS_chi_conductivity(double rhor);
		double ECS_f_int(double T);
		double surface_tension_T(double T);
	};
#endif
