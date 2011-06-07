
void setCoeffs(void);

double p_R744(double T, double rho);
double rho_R744(double T, double p_rho, int Types);
double h_R744(double T, double p_rho, int Types);
double u_R744(double T, double p_rho, int Types);
double s_R744(double T, double p_rho, int Types);
double cv_R744(double T, double p_rho, int Types);
double cp_R744(double T, double p_rho, int Types);
double c_R744(double T, double p_rho, int Types);
double visc_R744(double T, double p_rho, int Types);
double k_R744(double T, double p_rho, int Types);
double w_R744(double T, double p_rho, int Types);


double MM_R744(void);
double rhosatV_R744(double T);
double rhosatL_R744(double T);
double pcrit_R744(void);
double Tcrit_R744(void);
double rhocrit_R744(void);
double Ttriple_R744(void);
double psat_R744(double T);
double Tsat_R744(double P);

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

// Derivatives
double dhdrho_R744(double T, double p_rho, int Types);
double dhdT_R744(double T, double p_rho, int Types);
double dpdrho_R744(double T, double p_rho, int Types);
double dpdT_R744(double T, double p_rho, int Types);

int errCode_R744(void);


