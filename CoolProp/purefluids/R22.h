#ifndef R22_H
#define R22_H

	class R22Class : public Fluid{

	public:
		R22Class();
		~R22Class(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);

		void ECSParams(double *e_k, double *sigma);
		double ECS_psi_viscosity(double rhor);
	};
#endif
