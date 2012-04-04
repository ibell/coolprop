#ifndef R134a_H
#define R134a_H
	
	int Load_R134a(struct fluidParamsVals *Fluid);

	double psat_R134a(double T);
	double rhosatL_R134a(double T);
	double rhosatV_R134a(double T);

	double Conductivity_Trho_R134a(double T, double rho);
	double Viscosity_Trho_R134a(double T, double rho);

	//Residual Helmholtz formulation and derivatives
	double phir_R134a(double tau, double delta);
	double dphir_dDelta_R134a(double tau, double delta);
	double dphir2_dDelta2_R134a(double tau, double delta);
	double dphir2_dDelta_dTau_R134a(double tau, double delta);
	double dphir_dTau_R134a(double tau, double delta);
	double dphir2_dTau2_R134a(double tau, double delta);
	double phi0_R134a(double tau, double delta);
	double dphi0_dDelta_R134a(double tau, double delta);
	double dphi02_dDelta2_R134a(double tau, double delta);
	double dphi0_dTau_R134a(double tau, double delta);
	double dphi02_dTau2_R134a(double tau, double delta);


#endif
