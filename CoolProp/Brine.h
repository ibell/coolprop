/*
The name of this file is a slight misnomer.  It includes the analysis
aqueous solutions as well as single-phase liquids.
*/

#ifndef BRINE_H
#define BRINE_H

#include <string>
double SecFluidsSI(std::string Output, double T, double p, std::string Ref);
double SecFluids(std::string Output, double T, double p, std::string Ref);
int Brine(const char * Mix, double T, double C, /*in --- out */double *Tfreeze, double *Tmax, double *rho, double *cp, double *k, double *visc, double *h, double *s);


#endif
