
#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include "math.h"
#include "stdio.h"
#include <string.h>
#include "CoolProp.h"
#include <vector>
#include <iostream>
#include <list>
#include "Helmholtz.h"
#include "FluidClass.h"
#include "Neon.h"
#include "Nitrogen.h"

using namespace std;

NeonClass::NeonClass()
{
	double n[] = {0, 3.532653449, -4.513954384, -0.152402796, 2.188568609, -7.44299997, 7.755627402, -3.122553128, 1.014206899, -0.052892141, 0.156684924, -0.222852705, -0.014101509, 0.070362297, -0.058820484, 0.015711727, 0.001292203, 0.000790204, -0.00037944, 0.046527993, 0.045240018, -0.238342199, 0.00629359, -0.001272314, -1.75235E-07, 0.007188419, -0.054030069, 0.075782222, -0.038085883, 0.006034022};
	double d[] = {0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 6, 6, 6, 1, 2, 2, 2, 2, 2, 4, 8, 8, 8, 8};
	double t[] = {0, 0.5, 0.75, 3.5, 0.5, 0.75, 1, 1.5, 2.5, 0.25, 0.5, 2.5, 1, 3, 4, 5, 1, 5, 6, 4, 1, 5, 8, 12, 32, 10, 6, 7, 8, 9};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 4, 6, 6, 2, 2, 2, 2, 2};

	phirlist.push_back(new phir_power(n,d,t,c,1,29,30));

	// Critical parameters
	crit.rho = 23.882*20.179;
	crit.p = 2680.0;
	crit.T = 44.4918;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 20.179;
	params.Ttriple = 24.56;
	params.ptriple = 43.432339578188873;
	params.accentricfactor = -0.038449299273685900;
	params.R_u = 8.31434;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 2000.0;
	limits.pmax = 1000000.0;
	limits.rhomax = 50.65*params.molemass;

	//Constants for ideal gas expression
	phi0list.push_back(new phi0_lead(0,0));
	phi0list.push_back(new phi0_logtau(1.5));

	name.assign("Neon");
	aliases.push_back("neon");
	REFPROPname.assign("Neon");

	BibTeXKeys.EOS = "Katti-ACE-1986";
	BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double NeonClass::psat(double T)
{
    // Maximum absolute error is 0.090446 % between 24.556001 K and 44.491790 K
    const double t[]={0, 1, 2, 3, 6};
    const double N[]={0, 0.030056278177045211, -5.9517425005604085, 1.4556263648417147, -0.72769727713194987};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1;i<=3;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p*exp(reduce.T/T*summer);
}

double NeonClass::rhosatL(double T)
{
    // Maximum absolute error is 0.184674 % between 24.556001 K and 44.491790 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.8333333333333333};
    const double N[] = {0, -7.2369686484173252, 211.94763997778651, -2354.4170363141825, 13773.297638494711, -47566.671061179397, 101604.99338894009, -134063.9401911686, 103055.345207873, -37489.598569658992, 2840.9862970955251};
    double summer=0,theta;

    theta=1-T/reduce.T;
	for (int i=1; i<=10; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double NeonClass::rhosatV(double T)
{
    // Maximum absolute error is 0.152314 % between 24.556001 K and 44.491790 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333};
    const double N[] = {0, -0.42110191619014642, 9.2190373809956672, -72.08994032423108, 243.69954118084271, -458.84607430631803, 462.45011147675018, -229.64665885151433, 41.165593531034638};
    double summer=0,theta;

    theta=1-T/reduce.T;
	for (int i=1; i<=8; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}
double NeonClass::viscosity_Trho(double T, double rho)
{
	// Use nitrogen as the reference
	Fluid * ReferenceFluid = new NitrogenClass();
	ReferenceFluid->post_load();
	// Calculate the ECS
	double mu = viscosity_ECS_Trho(T, rho, ReferenceFluid);
	// Delete the reference fluid instance
	delete ReferenceFluid;
	return mu;
}
double NeonClass::conductivity_Trho(double T, double rho)
{
	// Use nitrogen as the reference
	Fluid * ReferenceFluid = new NitrogenClass();
	ReferenceFluid->post_load();
	// Calculate the ECS
	double cond = conductivity_ECS_Trho(T, rho, ReferenceFluid);
	// Delete the reference fluid instance
	delete ReferenceFluid;
	return cond;
}