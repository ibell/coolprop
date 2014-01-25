#ifndef OXYGEN_H
#define OXYGEN_H

	class OxygenClass : public Fluid{

	public:
		OxygenClass();
		~OxygenClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);

		double X_tilde(double T,double tau,double delta);
		double surface_tension_T(double T);
	};

	double rhosatV_Oxygen(double T);
	double rhosatL_Oxygen(double T);
	double psat_Oxygen(double T);

	double Viscosity_Trho_Oxygen(double T, double rho);
	double Conductivity_Trho_Oxygen(double T, double rho);

#endif
