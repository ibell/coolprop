#ifndef BRINE_H
#define BRINE_H
double SecFluids(char Output, double T, double p,char * Ref);
int Brine(char * Mix, double T, double C, /*in --- out */double *Tfreeze, double *Tmax, double *rho, double *cp, double *k, double *visc, double *h, double *s);
#endif