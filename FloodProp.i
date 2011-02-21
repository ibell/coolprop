%module FloodProp
%{
#include "FloodProp.h"
%}

%constant int I_N2=0;
//, I_He=1, I_Ne=2, I_Ar=3, I_Kr=4, I_Xe=5, I_CO2=6;
//%constant const int I_Methanol=0, I_Ethanol=1, I_Propanol=2, I_Butanol=3, I_Water=4, I_NH3=5, I_Zerol=6;

%include "FloodProp.c"

