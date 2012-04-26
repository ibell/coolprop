#ifndef AIR_H
#define AIR_H
	// -------------------
	// Required parameters
	// -------------------

	int Load_Air(struct fluidParamsVals *Fluid);

	double pdp_Air(double T);
    double pbp_Air(double T);
    double rhosatV_Air(double T);
	double rhosatL_Air(double T);
    
	double Viscosity_Trho_Air(double T, double rho);
	double Conductivity_Trho_Air(double T, double rho);

	//Residual Helmholtz formulation and derivatives
	double phir_Air(double tau, double delta);
	double dphir_dDelta_Air(double tau, double delta);
	double dphir2_dDelta2_Air(double tau, double delta);
	double dphir2_dDelta_dTau_Air(double tau, double delta);
    double dphir3_dDelta2_dTau_Air(double tau, double delta);
	double dphir_dTau_Air(double tau, double delta);
	double dphir2_dTau2_Air(double tau, double delta);
	double phi0_Air(double tau, double delta);
	double dphi0_dDelta_Air(double tau, double delta);
	double dphi02_dDelta2_Air(double tau, double delta);
	double dphi0_dTau_Air(double tau, double delta);
	double dphi02_dTau2_Air(double tau, double delta);

	// -------------------
	// Optional parameters
	// -------------------

	// Virial terms and their derivatives (for humid air properties
	double B_Air(double tau);
	double dBdT_Air(double tau);
	double C_Air(double tau);
	double dCdT_Air(double tau);

#endif
