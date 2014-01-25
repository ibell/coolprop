#ifndef R290_H
#define R290_H

	class R290Class : public Fluid{

	public:
		R290Class();
		~R290Class(){};
		double conductivity_Trho(double, double);
		double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);

		double viscosity_dilute(double T);
		double viscosity_dilute2(double T, double rho);
		double viscosity_residual(double T, double rho);
		double viscosity_higher_order(double T, double rho);
		double viscosity_background(double T, double rho);
		double conductivity_background(double T, double rho);
		void ECSParams(double *e_k, double *sigma);
		double surface_tension_T(double T);
	};

#endif
