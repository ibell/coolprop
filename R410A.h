
#ifndef R410_H
#define R410_H

    double rho_R410A(double T, double p, int Types);
    double p_R410A(double T, double rho);
    double h_R410A(double T, double p_rho, int Types);
    double s_R410A(double T, double p_rho, int Types);
    double u_R410A(double T, double p_rho, int Types);
    double cp_R410A(double T, double p_rho, int Types);
    double cv_R410A(double T, double p_rho, int Types);
    double visc_R410A(double T, double p_rho, int Types);
    double k_R410A(double T, double p_rho, int Types);
    double w_R410A(double T, double p_rho, int Types);

    double MM_R410A(void);

    double Tcrit_R410A(void);
    double pcrit_R410A(void);
    double rhocrit_R410A(void);
    double Ttriple_R410A(void);

    double p_dp_R410A(double T);
    double p_bp_R410A(double T);
    double rhosatL_R410A(double T);
    double rhosatV_R410A(double T);

    // Derivatives
    double dhdrho_R410A(double T, double p_rho, int Types);
    double dhdT_R410A(double T, double p_rho, int Types);
    double dpdrho_R410A(double T, double p_rho, int Types);
    double dpdT_R410A(double T, double p_rho, int Types);

    int errCode_R410A(void);
    
    int Load_R410A(struct fluidParamsVals *Fluid);
    
    double ar_R410A(double tau, double delta);
    double dar_dtau_R410A(double tau,double delta);
    double d2ar_dtau2_R410A(double tau, double delta);
    double dar_ddelta_R410A(double tau,double delta);
    double d2ar_ddelta2_R410A(double tau,double delta);
    double d2ar_ddelta_dtau_R410A(double tau,double delta);
    
    double a0_R410A(double tau, double delta);
    double da0_dtau_R410A(double tau,double delta);
    double d2a0_dtau2_R410A(double tau, double delta);
    double da0_ddelta_R410A(double tau,double delta);
    double d2a0_ddelta2_R410A(double tau,double delta);
    double d2a0_ddelta_dtau_R410A(double tau,double delta);

#endif
