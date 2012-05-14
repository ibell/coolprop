#ifndef _WATER_H
#define WATER_H

	int Load_Water(struct fluidParamsVals *Fluid);

	double rhosatV_Water(double T);
	double rhosatL_Water(double T);
	double psat_Water(double T);

	double Viscosity_Trho_Water(double T,double rho);
	double Conductivity_Trho_Water(double T,double rho);

	//Residual Helmholtz formulation and derivatives
	double phir_Water(double tau, double delta);
	double dphir_dDelta_Water(double tau, double delta);
	double dphir2_dDelta2_Water(double tau, double delta);
	double dphir2_dDelta_dTau_Water(double tau, double delta);
	double dphir_dTau_Water(double tau, double delta);
	double dphir2_dTau2_Water(double tau, double delta);
	double phi0_Water(double tau, double delta);
	double dphi0_dDelta_Water(double tau, double delta);
	double dphi02_dDelta2_Water(double tau, double delta);
	double dphi0_dTau_Water(double tau, double delta);
	double dphi02_dTau2_Water(double tau, double delta);
	double dphir3_dDelta2_dTau_Water(double tau, double delta);

	// Virial terms and their derivatives
	double B_Water(double tau);
	double dBdT_Water(double tau);
	double C_Water(double tau);
	double dCdT_Water(double tau);

    double IsothermCompress_Water(double T, double p);

#endif
