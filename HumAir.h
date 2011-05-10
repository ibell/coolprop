#ifndef GIVEN_TDP
#define GIVEN_TDP 1
#endif

#ifndef GIVEN_HUMRAT
#define GIVEN_HUMRAT 2
#endif

#ifndef GIVEN_TWB
#define GIVEN_TWB 3
#endif

#ifndef GIVEN_RH
#define GIVEN_RH 4
#endif

#ifndef GIVEN_ENTHALPY
#define GIVEN_ENTHALPY 5
#endif

#ifndef HUMAIR_H
#define HUMAIR_H

int HumAir(double tSI, double pSI, int HumInput, double xSI, /* in --- out */ double *Tdp_out, double *W_out, double *h_out, double *RH_out, double *v_out);
double cair_sat(double T);
double hair_sat(double T);
double T_hss(double h, double p, double T_guess);
double T_homega(double h, double omega, double P, double T_guess, double deltaT);
void Help();
int returnHumAirCode(char * Code);
double HumAir_Single(double T, double pSI, char *HumInputStr, double xSI, char *OutputStr);

#endif