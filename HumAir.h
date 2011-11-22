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

#ifndef GIVEN_T
#define GIVEN_T 6
#endif

#ifndef GIVEN_P
#define GIVEN_P 7
#endif


#ifndef HUMAIR_H
#define HUMAIR_H

// -----------------------
// Standard I/O functions
// -----------------------
double HAProps(char *OutputName, char *Input1Name, double Input1, char *Input2Name, double Input2, char *Input3Name, double Input3);
// HumAir is the backwards compatible call for consistency
void HumAir(double tSI, double pSI, int HumInput, double xSI, /* in --- out */ double *Tdp_out, double *W_out, double *h_out, double *RH_out, double *v_out);
// HumAir_Single is more extensible
double HumAir_Single(double T, double p, char *HumInputStr, double xSI, char *OutputStr);

// --------------
// Help functions
// --------------
void HAHelp();
int returnHumAirCode(char * Code);

// ----------------------
// Other simple functions
// ----------------------
double cair_sat(double T);

#endif