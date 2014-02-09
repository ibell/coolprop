
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

NeonClass::NeonClass()
{
	double n[] = {0, 3.532653449, -4.513954384, -0.1524027959, 2.188568609, -7.44299997, 7.755627402, -3.122553128, 1.014206899, -0.05289214086, 0.1566849239, -0.2228527050, -0.01410150942, 0.07036229719, -0.05882048367, 0.01571172741, 0.001292202769, 0.0007902035603, -0.0003794403616, 0.04652799333, 0.04524001818, -0.2383421991, 0.00629359013, -0.001272313644, -1.75235256E-07, 0.007188419232, -0.05403006914, 0.07578222187, -0.03808588254, 0.006034022431};
	double d[] = {0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 6, 6, 6, 1, 2, 2, 2, 2, 2, 4, 8, 8, 8, 8};
	double t[] = {0, 0.5, 0.75, 3.5, 0.5, 0.75, 1, 1.5, 2.5, 0.25, 0.5, 2.5, 1, 3, 4, 5, 1, 5, 6, 4, 1, 5, 8, 12, 32, 10, 6, 7, 8, 9};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 4, 6, 6, 2, 2, 2, 2, 2};

	phirlist.push_back(new phir_power(n,d,t,c,1,29,30));

	// Critical parameters
	crit.rho = 23.882*20.179;
	crit.p = PressureUnit(2680.0, UNIT_KPA);
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
	aliases.push_back(std::string("NEON"));
	REFPROPname.assign("Neon");

	ECSReferenceFluid = "Nitrogen";

	BibTeXKeys.EOS = "Katti-ACE-1986";
	BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double NeonClass::psat(double T)
{
    // Max error is  0.0338668363766 % between 24.56 and 44.491799 K
    const double t[]={0, 0.373, 0.39149999999999996, 0.39199999999999996, 0.39699999999999996, 0.8333333333333334, 4.333333333333333};
    const double N[]={0, -607.32379162883592, 125090.01843013773, -134770.7106365099, 10289.950413284598, -6.6808869319498001, -0.71410047975147495};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
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
