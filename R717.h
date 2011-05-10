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
	double Ttriple_R717(void);

	double psat_R717(double T);
	double rhosatL_R717(double T); 
	double rhosatV_R717(double T);
#endif
