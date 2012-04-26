#ifndef R744_H
#define R744_H

	int Load_R744(struct fluidParamsVals *Fluid);

	double rhosatV_R744(double T);
	double rhosatL_R744(double T);
	double psat_R744(double T);

	double Viscosity_Trho_R744(double T,double rho);
	double Conductivity_Trho_R744(double T,double rho);

	//Residual Helmholtz formulation and derivatives
	double phir_R744(double tau, double delta);
	double dphir_dDelta_R744(double tau, double delta);
	double dphir2_dDelta2_R744(double tau, double delta);
	double dphir2_dDelta_dTau_R744(double tau, double delta);
	double dphir_dTau_R744(double tau, double delta);
	double dphir2_dTau2_R744(double tau, double delta);
	double phi0_R744(double tau, double delta);
	double dphi0_dDelta_R744(double tau, double delta);
	double dphi02_dDelta2_R744(double tau, double delta);
	double dphi0_dTau_R744(double tau, double delta);
	double dphi02_dTau2_R744(double tau, double delta);

#endif
