#ifndef ARGON_H
#define ARGON_H

	void setCoeffs(void);

	double p_Argon(double T, double rho);
	double rho_Argon(double T, double p_rho, int Types);
	double h_Argon(double T, double p_rho, int Types);
	double u_Argon(double T, double p_rho, int Types);
	double s_Argon(double T, double p_rho, int Types);
	double cv_Argon(double T, double p_rho, int Types);
	double cp_Argon(double T, double p_rho, int Types);
	double c_Argon(double T, double p_rho, int Types);
	double k_Argon(double T, double p_rho, int Types);
	double visc_Argon(double T, double p_rho, int Types);
	double w_Argon(double T, double p_rho, int Types);
	double MM_Argon(void);

	double rhosatV_Argon(double T);
	double rhosatL_Argon(double T);

	double pcrit_Argon(void);
	double Tcrit_Argon(void);
	double rhocrit_Argon(void);
	double Ttriple_Argon(void);

	double rhosat_Argon(double T,double x);
	double psat_Argon(double T);
	double Tsat_Argon(double P);

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

	// Derivatives
	double dhdrho_Argon(double T, double p_rho, int Types);
	double dhdT_Argon(double T, double p_rho, int Types);
	double dpdrho_Argon(double T, double p_rho, int Types);
	double dpdT_Argon(double T, double p_rho, int Types);


	int errCode_Argon(void);

#endif
