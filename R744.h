
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

double hsat_R744(double T,double x);
double rhosat_R744(double T,double x);
double ssat_R744(double T,double x);
double Psat_R744(double T);
double Tsat_R744(double P);

// Derivatives
double dhdrho_R744(double T, double p_rho, int Types);
double dhdT_R744(double T, double p_rho, int Types);
double dpdrho_R744(double T, double p_rho, int Types);
double dpdT_R744(double T, double p_rho, int Types);

int errCode_R744(void);


