#ifndef R134a_H
#define R134a_H

	class R134aClass : public Fluid{

	public:
		R134aClass();
		~R134aClass(){};
		double conductivity_Trho(double, double);
		double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);

		double viscosity_background(double T, double rho);
		double viscosity_dilute(double T);
		double viscosity_residual(double T, double rho);
		double conductivity_background(double T, double rho);
		double surface_tension_T(double T);
		double conductivity_dilute(double T);
		double conductivity_residual(double T, double rho);
		void ECSParams(double *e_k, double *sigma);
	};
#endif
