#ifndef R717_H
#define R717_H
	
	double p_R717(double T, double rho);
	double rho_R717(double T, double p, int Types);
	double h_R717(double T, double p_rho, int Types);
	double s_R717(double T, double p_rho, int Types);
	double u_R717(double T, double p_rho, int Types);
	double cp_R717(double T, double p_rho, int Types);
	double cv_R717(double T, double p_rho, int Types);
	double visc_R717(double T, double p_rho, int Types);
	double k_R717(double T, double p_rho, int Types);
	double w_R717(double T, double p_rho, int Types);

	double MM_R717(void);
	int errCode_R717(void);

	double Tcrit_R717(void); 
	double pcrit_R717(void); 
	double rhocrit_R717(void);
	double Ttriple_R717(void);

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

	double psat_R717(double T);
	double rhosatL_R717(double T); 
	double rhosatV_R717(double T);
#endif
