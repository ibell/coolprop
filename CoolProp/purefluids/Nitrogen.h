#ifndef NITROGEN_H
#define NITROGEN_H

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
