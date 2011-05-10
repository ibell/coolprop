
#ifndef R410_H
#define R410_H

	double rho_R410A(double T, double p, int Types);
	double p_R410A(double T, double rho);
	double h_R410A(double T, double p_rho, int Types);
	double s_R410A(double T, double p_rho, int Types);
	double u_R410A(double T, double p_rho, int Types);
	double cp_R410A(double T, double p_rho, int Types);
	double cv_R410A(double T, double p_rho, int Types);
	double visc_R410A(double T, double p_rho, int Types);
	double k_R410A(double T, double p_rho, int Types);
	double w_R410A(double T, double p_rho, int Types);

	double MM_R410A(void);

	double Tcrit_R410A(void);
	double pcrit_R410A(void);
	double Ttriple_R410A(void);

	double p_dp_R410A(double T);
	double p_bp_R410A(double T);
	double rhosatL_R410A(double T);
	double rhosatV_R410A(double T);

	// Derivatives
	double dhdrho_R410A(double T, double p_rho, int Types);
	double dhdT_R410A(double T, double p_rho, int Types);
	double dpdrho_R410A(double T, double p_rho, int Types);
	double dpdT_R410A(double T, double p_rho, int Types);

	int errCode_R410A(void);

#endif
