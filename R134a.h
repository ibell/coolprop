#ifndef R134a_H
#define R134a_H
	
	double p_R134a(double T, double rho);
	double rho_R134a(double T, double p, int Types);
	double h_R134a(double T, double p_rho, int Types);
	double s_R134a(double T, double p_rho, int Types);
	double u_R134a(double T, double p_rho, int Types);
	double cp_R134a(double T, double p_rho, int Types);
	double cv_R134a(double T, double p_rho, int Types);
	double visc_R134a(double T, double p_rho, int Types);
	double k_R134a(double T, double p_rho, int Types);
	double w_R134a(double T, double p_rho, int Types);

	double MM_R134a(void);
	int errCode_R134a(void);

	double Tcrit_R134a(void); 
	double pcrit_R134a(void); 
	double Ttriple_R134a(void);

	// Derivatives
	double dhdrho_R134a(double T, double p_rho, int Types);
	double dhdT_R134a(double T, double p_rho, int Types);
	double dpdrho_R134a(double T, double p_rho, int Types);
	double dpdT_R134a(double T, double p_rho, int Types);

	double psat_R134a(double T);
	double rhosatL_R134a(double T); 
	double rhosatV_R134a(double T);
#endif
