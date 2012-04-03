#ifndef R134a_H
#define R134a_H
	
	double p_R134a(double T, double rho);
	double rho_R134a(double T, double p, int Types);
	double h_R134a(double T, double p_rho, int Types);
	double s_R134a(double T, double p_rho, int Types);
	double u_R134a(double T, double p_rho, int Types);
	double cp_R134a(double T, double p_rho, int Types);
	double cv_R134a(double T, double p_rho, int Types);
	double visc_R134a(double T, double p_rho, int Types);
	double k_R134a(double T, double p_rho, int Types);
	double w_R134a(double T, double p_rho, int Types);

	double MM_R134a(void);
	int errCode_R134a(void);

	double Tcrit_R134a(void); 
	double pcrit_R134a(void); 
	double rhocrit_R134a(void);
	double Ttriple_R134a(void);
    
    int Load_R134a(struct fluidParamsVals *Fluid);

	//Residual Helmholtz formulation and derivatives
	double phir_R134a(double tau, double delta);
	double dphir_dDelta_R134a(double tau, double delta);
	double dphir2_dDelta2_R134a(double tau, double delta);
	double dphir2_dDelta_dTau_R134a(double tau, double delta);
	double dphir_dTau_R134a(double tau, double delta);
	double dphir2_dTau2_R134a(double tau, double delta);
	double phi0_R134a(double tau, double delta);
	double dphi0_dDelta_R134a(double tau, double delta);
	double dphi02_dDelta2_R134a(double tau, double delta);
	double dphi0_dTau_R134a(double tau, double delta);
	double dphi02_dTau2_R134a(double tau, double delta);

	double psat_R134a(double T);
	double rhosatL_R134a(double T); 
	double rhosatV_R134a(double T);
#endif
