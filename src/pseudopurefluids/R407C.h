
#ifndef R407C_H
#define R407C_H

	int Load_R407C(struct fluidParamsVals *Fluid);

	double p_dp_R407C(double T);
	double p_bp_R407C(double T);
	double rhosatL_R407C(double T);
	double rhosatV_R407C(double T);

	double Viscosity_Trho_R407C(double T, double rho);
	double Conductivity_Trho_R407C(double T, double rho);

	double a0_R407C(double tau, double delta);
	double da0_dtau_R407C(double tau, double delta);
	double d2a0_dtau2_R407C(double tau, double delta);
	double da0_ddelta_R407C(double tau, double delta);
	double d2a0_ddelta2_R407C(double tau, double delta);

	double ar_R407C(double tau, double delta);
	double dar_dtau_R407C(double tau,double delta);
	double d2ar_dtau2_R407C(double tau, double delta);
	double dar_ddelta_R407C(double tau,double delta);
	double d2ar_ddelta2_R407C(double tau,double delta);
	double d2ar_ddelta_dtau_R407C(double tau,double delta);

#endif
