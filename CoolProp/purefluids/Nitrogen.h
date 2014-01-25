#ifndef NITROGEN_H
#define NITROGEN_H

	class NitrogenClass : public Fluid{

	public:
		NitrogenClass();
		~NitrogenClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		
		double viscosity_dilute(double T);
		double viscosity_residual(double T, double rho);
		double viscosity_background(double T, double rho);
		double conductivity_dilute(double T);
		double conductivity_background(double T, double rho);
		double conductivity_critical(double T, double rho);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		void ECSParams(double *e_k, double *sigma){
			// From Poling
			*e_k = 71.4;	
			*sigma = 0.3798;
		};

		double X_tilde(double T,double tau,double delta);
		double surface_tension_T(double T);
	};

	int Load_Nitrogen(struct fluidParamsVals *Fluid);

	double rhosatV_Nitrogen(double T);
	double rhosatL_Nitrogen(double T);
	double psat_Nitrogen(double T);

	double Viscosity_Trho_Nitrogen(double T, double rho);
	double Conductivity_Trho_Nitrogen(double T, double rho);

	//Residual Helmholtz formulation and derivatives
	double phir_Nitrogen(double tau, double delta);
	double dphir_dDelta_Nitrogen(double tau, double delta);
	double dphir2_dDelta2_Nitrogen(double tau, double delta);
	double dphir2_dDelta_dTau_Nitrogen(double tau, double delta);
	double dphir_dTau_Nitrogen(double tau, double delta);
	double dphir2_dTau2_Nitrogen(double tau, double delta);
	double phi0_Nitrogen(double tau, double delta);
	double dphi0_dDelta_Nitrogen(double tau, double delta);
	double dphi02_dDelta2_Nitrogen(double tau, double delta);
	double dphi0_dTau_Nitrogen(double tau, double delta);
	double dphi02_dTau2_Nitrogen(double tau, double delta);

#endif
