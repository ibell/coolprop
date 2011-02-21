#ifndef ARGON_H
#define ARGON_H

	void setCoeffs(void);

	double p_Argon(double T, double rho);
	double rho_Argon(double T, double p_rho, int Types);
	double h_Argon(double T, double p_rho, int Types);
	double u_Argon(double T, double p_rho, int Types);
	double s_Argon(double T, double p_rho, int Types);
	double cv_Argon(double T, double p_rho, int Types);
	double cp_Argon(double T, double p_rho, int Types);
	double c_Argon(double T, double p_rho, int Types);
	double k_Argon(double T, double p_rho, int Types);
	double visc_Argon(double T, double p_rho, int Types);
	double w_Argon(double T, double p_rho, int Types);
	double MM_Argon(void);

	double rhosatV_Argon(double T);
	double rhosatL_Argon(double T);

	double pcrit_Argon(void);
	double Tcrit_Argon(void);

	double hsat_Argon(double T,double x);
	double rhosat_Argon(double T,double x);
	double ssat_Argon(double T,double x);
	double psat_Argon(double T);
	double Tsat_Argon(double P);

	// Derivatives
	double dhdrho_Argon(double T, double p_rho, int Types);
	double dhdT_Argon(double T, double p_rho, int Types);
	double dpdrho_Argon(double T, double p_rho, int Types);
	double dpdT_Argon(double T, double p_rho, int Types);


	int errCode_Argon(void);

#endif
