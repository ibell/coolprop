
#ifndef R407C_H
#define R407C_H

	double rho_R407C(double T, double p, int Types);
	double p_R407C(double T, double rho);
	double h_R407C(double T, double p_rho, int Types);
	double s_R407C(double T, double p_rho, int Types);
	double u_R407C(double T, double p_rho, int Types);
	double cp_R407C(double T, double p_rho, int Types);
	double cv_R407C(double T, double p_rho, int Types);
	double visc_R407C(double T, double p_rho, int Types);
	double k_R407C(double T, double p_rho, int Types);
	double w_R407C(double T, double p_rho, int Types);

	double MM_R407C(void);

	double Tcrit_R407C(void);
	double pcrit_R407C(void);
	double Ttriple_R407C(void);

	double p_dp_R407C(double T);
	double p_bp_R407C(double T);
	double rhosatL_R407C(double T);
	double rhosatV_R407C(double T);

	int errCode_R407C(void);

#endif
