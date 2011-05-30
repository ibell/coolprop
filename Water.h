
void setCoeffs(void);

double p_Water(double T, double rho);
double rho_Water(double T, double p_rho, int Types);
double h_Water(double T, double p_rho, int Types);
double u_Water(double T, double p_rho, int Types);
double s_Water(double T, double p_rho, int Types);
double cv_Water(double T, double p_rho, int Types);
double cp_Water(double T, double p_rho, int Types);
double c_Water(double T, double p_rho, int Types);
double visc_Water(double T, double p_rho, int Types);
double k_Water(double T, double p_rho, int Types);
double w_Water(double T, double p_rho, int Types);


double MM_Water(void);
double rhosatV_Water(double T);
double rhosatL_Water(double T);
double pcrit_Water(void);
double Tcrit_Water(void);
double rhocrit_Water(void);
double Ttriple_Water(void);
double Psat_Water(double T);
double Tsat_Water(double P);

//Residual Helmholtz formulation and derivatives
double phir_Water(double tau, double delta);
double dphir_dDelta_Water(double tau, double delta);
double dphir2_dDelta2_Water(double tau, double delta);
double dphir2_dDelta_dTau_Water(double tau, double delta);
double dphir_dTau_Water(double tau, double delta);
double dphir2_dTau2_Water(double tau, double delta);
double phi0_Water(double tau, double delta);
double dphi0_dDelta_Water(double tau, double delta);
double dphi02_dDelta2_Water(double tau, double delta);
double dphi0_dTau_Water(double tau, double delta);
double dphi02_dTau2_Water(double tau, double delta);

// Derivatives
double dhdrho_Water(double T, double p_rho, int Types);
double dhdT_Water(double T, double p_rho, int Types);
double dpdrho_Water(double T, double p_rho, int Types);
double dpdT_Water(double T, double p_rho, int Types);

int errCode_Water(void);


