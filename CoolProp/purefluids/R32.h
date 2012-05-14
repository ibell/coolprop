#ifndef R32_H
#define R32_H

	double psat_R32(double T);
	double rhosatL_R32(double T); 
	double rhosatV_R32(double T);

	double Viscosity_Trho_R32(double T, double rho);
	double Conductivity_Trho_R32(double T, double rho);

	double phi0_R32(double tau,double delta);
	double dphi0_ddelta_R32(double tau,double delta);
	double d2phi0_ddelta2_R32(double tau,double delta);
	double dphi0_dtau_R32(double tau,double delta);
	double d2phi0_dtau2_R32(double tau,double delta);
	double d2phi0_ddelta_dtau_R32(double tau, double delta);

	double phir_R32(double tau,double delta);
	double dphir_ddelta_R32(double tau, double delta);
	double d2phir_ddelta2_R32(double tau, double delta);
	double dphir_dtau_R32(double tau, double delta);
	double d2phir_dtau2_R32(double tau, double delta);
	double d2phir_ddelta_dtau_R32(double tau, double delta);
#endif
