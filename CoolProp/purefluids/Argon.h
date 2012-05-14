#ifndef ARGON_H
#define ARGON_H

	int Load_Argon(struct fluidParamsVals *Fluid);

	double rhosatV_Argon(double T);
	double rhosatL_Argon(double T);
	double psat_Argon(double T);

	double Viscosity_Trho_Argon(double T, double rho);
	double Conductivity_Trho_Argon(double T, double rho);

	//Residual Helmholtz formulation and derivatives
	double phir_Argon(double tau, double delta);
	double dphir_dDelta_Argon(double tau, double delta);
	double dphir2_dDelta2_Argon(double tau, double delta);
	double dphir2_dDelta_dTau_Argon(double tau, double delta);
	double dphir_dTau_Argon(double tau, double delta);
	double dphir2_dTau2_Argon(double tau, double delta);
	double phi0_Argon(double tau, double delta);
	double dphi0_dDelta_Argon(double tau, double delta);
	double dphi02_dDelta2_Argon(double tau, double delta);
	double dphi0_dTau_Argon(double tau, double delta);
	double dphi02_dTau2_Argon(double tau, double delta);

#endif
