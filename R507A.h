
#ifndef R507A_H
#define R507A_H

	double rho_R507A(double T, double p, int Types);
	double p_R507A(double T, double rho);
	double h_R507A(double T, double p_rho, int Types);
	double s_R507A(double T, double p_rho, int Types);
	double u_R507A(double T, double p_rho, int Types);
	double cp_R507A(double T, double p_rho, int Types);
	double cv_R507A(double T, double p_rho, int Types);
	double visc_R507A(double T, double p_rho, int Types);
	double k_R507A(double T, double p_rho, int Types);
	double w_R507A(double T, double p_rho, int Types);

	double MM_R507A(void);

	double Tcrit_R507A(void);
	double pcrit_R507A(void);

	double p_dp_R507A(double T);
	double p_bp_R507A(double T);
	double rhosatL_R507A(double T);
	double rhosatV_R507A(double T);

	int errCode_R507A(void);

#endif
