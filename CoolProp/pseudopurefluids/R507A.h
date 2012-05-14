
#ifndef R507A_H
#define R507A_H

	int Load_R507A(struct fluidParamsVals *Fluid);

	double p_dp_R507A(double T);
	double p_bp_R507A(double T);
	double rhosatL_R507A(double T);
	double rhosatV_R507A(double T);

	double Viscosity_Trho_R507A(double T, double rho);
	double Conductivity_Trho_R507A(double T, double rho);

	double a0_R507A(double tau, double delta);
	double da0_dtau_R507A(double tau, double delta);
	double d2a0_dtau2_R507A(double tau, double delta);
	double da0_ddelta_R507A(double tau, double delta);
	double d2a0_ddelta2_R507A(double tau, double delta);

	double ar_R507A(double tau, double delta);
	double dar_dtau_R507A(double tau,double delta);
	double d2ar_dtau2_R507A(double tau, double delta);
	double dar_ddelta_R507A(double tau,double delta);
	double d2ar_ddelta2_R507A(double tau,double delta);
	double d2ar_ddelta_dtau_R507A(double tau,double delta);

#endif
