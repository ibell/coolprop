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
	};
#endif
