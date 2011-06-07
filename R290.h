#ifndef R290_H
#define R290_H

	void setCoeffs(void);

	double p_R290(double T, double rho);
	double rho_R290(double T, double p_rho, int Types);
	double h_R290(double T, double p_rho, int Types);
	double u_R290(double T, double p_rho, int Types);
	double s_R290(double T, double p_rho, int Types);
	double cv_R290(double T, double p_rho, int Types);
	double cp_R290(double T, double p_rho, int Types);
	double c_R290(double T, double p_rho, int Types);
	double k_R290(double T, double p_rho, int Types);
	double visc_R290(double T, double p_rho, int Types);
	double w_R290(double T, double p_rho, int Types);
	double MM_R290(void);

	double rhosatV_R290(double T);
	double rhosatL_R290(double T);
	double psat_R290(double T);

	double pcrit_R290(void);
	double Tcrit_R290(void);
	double rhocrit_R290(void);
	double Ttriple_R290(void);

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


	int errCode_R290(void);

#endif
