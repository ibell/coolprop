
#ifndef R410_H
#define R410_H

    double Viscosity_Trho_R410A(double T, double rho);
    double Conductivity_Trho_R410A(double T, double rho);
    double p_dp_R410A(double T);
    double p_bp_R410A(double T);
    double rhosatL_R410A(double T);
    double rhosatV_R410A(double T);
    
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
