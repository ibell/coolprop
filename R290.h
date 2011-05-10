#ifndef R290_H
#define R290_H

	void setCoeffs(void);

	double p_R290(double T, double rho);
	double rho_R290(double T, double p_rho, int Types);
	double h_R290(double T, double p_rho, int Types);
	double u_R290(double T, double p_rho, int Types);
	double s_R290(double T, double p_rho, int Types);
	double cv_R290(double T, double p_rho, int Types);
	double cp_R290(double T, double p_rho, int Types);
	double c_R290(double T, double p_rho, int Types);
	double k_R290(double T, double p_rho, int Types);
	double visc_R290(double T, double p_rho, int Types);
	double w_R290(double T, double p_rho, int Types);
	double MM_R290(void);

	double rhosatV_R290(double T);
	double rhosatL_R290(double T);

	double pcrit_R290(void);
	double Tcrit_R290(void);
	double Ttriple_R290(void);

	double hsat_R290(double T,double x);
	double rhosat_R290(double T,double x);
	double ssat_R290(double T,double x);
	double psat_R290(double T);
	double Tsat_R290(double P);

	// Derivatives
	double dhdrho_R290(double T, double p_rho, int Types);
	double dhdT_R290(double T, double p_rho, int Types);
	double dpdrho_R290(double T, double p_rho, int Types);
	double dpdT_R290(double T, double p_rho, int Types);


	int errCode_R290(void);

#endif
