#include "PropMacros.h"

#ifndef R290_H
#define R290_H

	double p_R290(double T, double rho);	
	double rho_R290(double T, double p_rho, int Type);
	double h_R290(double T, double p_rho, int Type);
	double s_R290(double T, double p_rho, int Types);
	double u_R290(double T, double p_rho, int Types);
	double cp_R290(double T, double p_rho, int Types);
	double cv_R290(double T, double p_rho, int Types);
	double visc_R290(double T, double p_rho, int Types);
	double k_R290(double T, double p_rho, int Types);
	double w_R290(double T, double p_rho, int Types);

	double MM_R290(void);
	int errCode_R290(void);

	double Tcrit_R290(void);
	double pcrit_R290(void);


	double psat_R290(double T);
	double rhosatL_R290(double T);
	double rhosatV_R290(double T);
#endif
