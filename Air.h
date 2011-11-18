#ifndef AIR_H
#define AIR_H

	void setCoeffs(void);

	double p_Air(double T, double rho);
	double rho_Air(double T, double p_rho, int Types);
	double h_Air(double T, double p_rho, int Types);
	double u_Air(double T, double p_rho, int Types);
	double s_Air(double T, double p_rho, int Types);
	double cv_Air(double T, double p_rho, int Types);
	double cp_Air(double T, double p_rho, int Types);
	double c_Air(double T, double p_rho, int Types);
	double k_Air(double T, double p_rho, int Types);
	double visc_Air(double T, double p_rho, int Types);
	double w_Air(double T, double p_rho, int Types);
	double MM_Air(void);

	double rhosatV_Air(double T);
	double rhosatL_Air(double T);

	double pcrit_Air(void);
	double Tcrit_Air(void);
	double rhocrit_Air(void);
	double Ttriple_Air(void);

	double rhosat_Air(double T,double x);
	double pdp_Air(double T);
    double pbp_Air(double T);
	double Tsat_Air(double P);
    
    // Virial terms and their derivatives
	double B_Air(double tau);
	double dBdT_Air(double tau);
	double C_Air(double tau);
	double dCdT_Air(double tau);

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

	// Derivatives
	double dhdrho_Air(double T, double p_rho, int Types);
	double dhdT_Air(double T, double p_rho, int Types);
	double dpdrho_Air(double T, double p_rho, int Types);
	double dpdT_Air(double T, double p_rho, int Types);


	int errCode_Air(void);

#endif
