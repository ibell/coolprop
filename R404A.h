
#ifndef R404A_H
#define R404A_H

	double rho_R404A(double T, double p, int Types);
	double p_R404A(double T, double rho);
	double h_R404A(double T, double p_rho, int Types);
	double s_R404A(double T, double p_rho, int Types);
	double u_R404A(double T, double p_rho, int Types);
	double cp_R404A(double T, double p_rho, int Types);
	double cv_R404A(double T, double p_rho, int Types);
	double visc_R404A(double T, double p_rho, int Types);
	double k_R404A(double T, double p_rho, int Types);
	double w_R404A(double T, double p_rho, int Types);

	double MM_R404A(void);

	double Tcrit_R404A(void);
	double pcrit_R404A(void);

	double p_dp_R404A(double T);
	double p_bp_R404A(double T);
	double rhosatL_R404A(double T);
	double rhosatV_R404A(double T);

	int errCode_R404A(void);

#endif
