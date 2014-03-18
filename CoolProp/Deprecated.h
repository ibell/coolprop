#ifndef DEPRECATED_H
#define DEPRECATED_H

EXPORT_CODE double CONVENTION rhosatL_anc(const char* Fluid, double T);
EXPORT_CODE double CONVENTION rhosatV_anc(const char* Fluid, double T);
EXPORT_CODE double CONVENTION psatL_anc(const char* Fluid, double T);
EXPORT_CODE double CONVENTION psatV_anc(const char* Fluid, double T);
// Expose some functions that are useful for ECS debugging
EXPORT_CODE double CONVENTION viscosity_dilute(const char* FluidName, double T);
EXPORT_CODE double CONVENTION viscosity_residual(const char* FluidName, double T, double rho);
EXPORT_CODE double CONVENTION conductivity_critical(const char* FluidName, double T, double rho);
EXPORT_CODE double CONVENTION conductivity_background(const char* FluidName, double T, double rho);
EXPORT_CODE double CONVENTION conformal_Trho(const char* FluidName, const char* ReferenceFluidName, double T, double rho, double *Tconform, double *rhoconform);

#endif