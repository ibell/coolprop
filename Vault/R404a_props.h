#ifndef TYPE_TP
#define TYPE_TP 1
#endif

#ifndef TYPE_Trho
#define TYPE_Trho 2
#endif

#ifndef R404a_H
#define R404a_H

void setCoeffs_R404a();
double get_Delta_R404a(double T, double p);


double u_R404a(double T, double param2, int typeInputs);
double h_R404a(double T, double param2, int typeInputs);
double s_R404a(double T, double param2, int typeInputs);
double cp_R404a(double T, double param2, int typeInputs);
double cv_R404a(double T, double param2, int typeInputs);
double visc_R404a(double T, double param2, int typeInputs);
double rho_R404a(double T, double p);
double Psat_R404a(double T, double x);
double MM_R404a();


#endif
