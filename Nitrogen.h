
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
double Ttriple_Nitrogen(void);

double hsat_Nitrogen(double T,double x);
double rhosat_Nitrogen(double T,double x);
double ssat_Nitrogen(double T,double x);
double psat_Nitrogen(double T);
double Tsat_Nitrogen(double P);

// Derivatives
double dhdrho_Nitrogen(double T, double p_rho, int Types);
double dhdT_Nitrogen(double T, double p_rho, int Types);
double dpdrho_Nitrogen(double T, double p_rho, int Types);
double dpdT_Nitrogen(double T, double p_rho, int Types);


int errCode_Nitrogen(void);


