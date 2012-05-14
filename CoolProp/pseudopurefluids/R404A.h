
#ifndef R404A_H
#define R404A_H

	int Load_R404A(struct fluidParamsVals *Fluid);

	double p_dp_R404A(double T);
	double p_bp_R404A(double T);
	double rhosatL_R404A(double T);
	double rhosatV_R404A(double T);
	double Viscosity_Trho_R404A(double T, double rho);
	double Conductivity_Trho_R404A(double T, double rho);

	double a0_R404A(double tau, double delta);
	double da0_dtau_R404A(double tau, double delta);
	double d2a0_dtau2_R404A(double tau, double delta);
	double da0_ddelta_R404A(double tau, double delta);
	double d2a0_ddelta2_R404A(double tau, double delta);

	double ar_R404A(double tau, double delta);
	double dar_dtau_R404A(double tau,double delta);
	double d2ar_dtau2_R404A(double tau, double delta);
	double dar_ddelta_R404A(double tau,double delta);
	double d2ar_ddelta2_R404A(double tau,double delta);
	double d2ar_ddelta_dtau_R404A(double tau,double delta);

#endif
