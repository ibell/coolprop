#ifndef R717_H
#define R717_H

	int Load_R717(struct fluidParamsVals *Fluid);

	double psat_R717(double T);
	double rhosatL_R717(double T);
	double rhosatV_R717(double T);
	
	double Viscosity_Trho_R717(double T, double rho);
	double Conductivity_Trho_R717(double T, double rho);

	//Residual Helmholtz formulation and derivatives
	double phir_R717(double tau, double delta);
	double dphir_dDelta_R717(double tau, double delta);
	double dphir2_dDelta2_R717(double tau, double delta);
	double dphir2_dDelta_dTau_R717(double tau, double delta);
	double dphir_dTau_R717(double tau, double delta);
	double dphir2_dTau2_R717(double tau, double delta);
	double phi0_R717(double tau, double delta);
	double dphi0_dDelta_R717(double tau, double delta);
	double dphi02_dDelta2_R717(double tau, double delta);
	double dphi0_dTau_R717(double tau, double delta);
	double dphi02_dTau2_R717(double tau, double delta);


#endif
