#ifndef R290_H
#define R290_H

	int Load_R290(struct fluidParamsVals *Fluid);

	double rhosatV_R290(double T);
	double rhosatL_R290(double T);
	double psat_R290(double T);

	double Viscosity_Trho_R290(double T, double rho);
	double Conductivity_Trho_R290(double T, double rho);

	//Residual Helmholtz formulation and derivatives
	double phir_R290(double tau, double delta);
	double dphir_dDelta_R290(double tau, double delta);
	double dphir2_dDelta2_R290(double tau, double delta);
	double dphir2_dDelta_dTau_R290(double tau, double delta);
	double dphir_dTau_R290(double tau, double delta);
	double dphir2_dTau2_R290(double tau, double delta);
	double phi0_R290(double tau, double delta);
	double dphi0_dDelta_R290(double tau, double delta);
	double dphi02_dDelta2_R290(double tau, double delta);
	double dphi0_dTau_R290(double tau, double delta);
	double dphi02_dTau2_R290(double tau, double delta);

	// Derivatives
	double dhdrho_R290(double T, double p_rho, int Types);
	double dhdT_R290(double T, double p_rho, int Types);
	double dpdrho_R290(double T, double p_rho, int Types);
	double dpdT_R290(double T, double p_rho, int Types);

#endif
