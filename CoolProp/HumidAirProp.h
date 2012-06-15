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

#ifndef GIVEN_VISC
#define GIVEN_VISC 8
#endif

#ifndef GIVEN_COND
#define GIVEN_COND 9
#endif

#ifndef HUMAIR_H
#define HUMAIR_H

// -----------------------
// Standard I/O function
// -----------------------
double HAProps(char *OutputName, char *Input1Name, double Input1, char *Input2Name, double Input2, char *Input3Name, double Input3);
// -----------------------
// Extra I/O function
// -----------------------
double HAProps_Aux(char* Name,double T, double p, double W, char *units);

//Turn on the use of virial correlations for air and water
void UseVirialCorrelations(int flag);
void UseIsothermCompressCorrelation(int flag);
void UseIdealGasEnthalpyCorrelations(int flag);

// --------------
// Help functions
// --------------
void HAHelp(void);
int returnHumAirCode(char * Code);

// ----------------------
// Other simple functions
// ----------------------
double cair_sat(double T);

#endif
