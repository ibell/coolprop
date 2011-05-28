
void setCoeffs(void);

double p_Nitrogen(double T, double rho);
double rho_Nitrogen(double T, double p_rho, int Types);
double h_Nitrogen(double T, double p_rho, int Types);
double u_Nitrogen(double T, double p_rho, int Types);
double s_Nitrogen(double T, double p_rho, int Types);
double cv_Nitrogen(double T, double p_rho, int Types);
double cp_Nitrogen(double T, double p_rho, int Types);
double c_Nitrogen(double T, double p_rho, int Types);
double k_Nitrogen(double T, double p_rho, int Types);
double visc_Nitrogen(double T, double p_rho, int Types);
double w_Nitrogen(double T, double p_rho, int Types);
double MM_Nitrogen(void);

double rhosatV_Nitrogen(double T);
double rhosatL_Nitrogen(double T);

double pcrit_Nitrogen(void);
double Tcrit_Nitrogen(void);
double rhocrit_Nitrogen(void);
double Ttriple_Nitrogen(void);

double psat_Nitrogen(double T);
double Tsat_Nitrogen(double P);

//Residual Helmholtz formulation and derivatives
double phir_Nitrogen(double tau, double delta);
double dphir_dDelta_Nitrogen(double tau, double delta);
double dphir2_dDelta2_Nitrogen(double tau, double delta);
double dphir2_dDelta_dTau_Nitrogen(double tau, double delta);
double dphir_dTau_Nitrogen(double tau, double delta);
double dphir2_dTau2_Nitrogen(double tau, double delta);
double phi0_Nitrogen(double tau, double delta);
double dphi0_dDelta_Nitrogen(double tau, double delta);
double dphi02_dDelta2_Nitrogen(double tau, double delta);
double dphi0_dTau_Nitrogen(double tau, double delta);
double dphi02_dTau2_Nitrogen(double tau, double delta);

// Derivatives
double dhdrho_Nitrogen(double T, double p_rho, int Types);
double dhdT_Nitrogen(double T, double p_rho, int Types);
double dpdrho_Nitrogen(double T, double p_rho, int Types);
double dpdT_Nitrogen(double T, double p_rho, int Types);


int errCode_Nitrogen(void);


