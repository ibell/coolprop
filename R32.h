#ifndef R32_H
#define R32_H
	
	double p_R32(double T, double rho);
	double rho_R32(double T, double p, int Types);
	double h_R32(double T, double p_rho, int Types);
	double s_R32(double T, double p_rho, int Types);
	double u_R32(double T, double p_rho, int Types);
	double cp_R32(double T, double p_rho, int Types);
	double cv_R32(double T, double p_rho, int Types);
	double visc_R32(double T, double p_rho, int Types);
	double k_R32(double T, double p_rho, int Types);
	double w_R32(double T, double p_rho, int Types);

	double MM_R32(void);
	int errCode_R32(void);

	double Tcrit_R32(void); 
	double pcrit_R32(void); 
	double Ttriple_R32(void);

	double psat_R32(double T);
	double rhosatL_R32(double T); 
	double rhosatV_R32(double T);
#endif