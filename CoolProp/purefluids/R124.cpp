#include "CoolProp.h"
#include "FluidClass.h"
#include "R124.h"

R124Class::R124Class()
{
	double n[] = {0.0, -0.1262962e-1, 0.2168373e1, -0.3330033e1, 0.1610361e0, -0.9666145e-4, 0.1191310e-1, -0.2880217e-2, 0.1681346e-2, 0.1594968e-4, 0.1289674e0, 0.1182213e-4, -0.4713997e0, -0.2412873e0, 0.6868066e0, -0.8621095e-1, 0.4728645e-5, 0.1487933e-1, -0.3001338e-1, 0.1849606e-2, 0.4126073e-3};
	double t[] = {0, 2, 0.5, 1, 0.5, 2.5, -1, 1, 0, -0.5, 1.5, 1, 2.5, -0.25, 1, 5, 2, 15, 20, 15, 45};
	double d[] = {0, 1, 1, 1, 2, 2, 3, 5, 6, 8, 2, 12, 1, 1, 1, 1, 15, 3, 3, 4, 9};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4};

	//Critical parameters
	crit.rho = 560; //[kg/m^3]
	crit.p = PressureUnit(3624.295, UNIT_KPA); //[kPa]
	crit.T = 395.425; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 136.4762;
	params.Ttriple = 75;
	params.accentricfactor = 0.28809508422142915;
	params.R_u = 8.314471;
	params.ptriple = 2.67380659624e-05;

	// Limits of EOS
	limits.Tmin = 120;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,20,21));

	double a0[] = {0, -11.669406, 9.8760443, 2.175638, -7.389735, 0.8736831, -0.1115133};
	double t0[] = {0, 0, 0, 0, -1, -2, -3};
	phi0list.push_back(new phi0_lead(a0[1],a0[2]));
	phi0list.push_back(new phi0_logtau(a0[3]));
	phi0list.push_back(new phi0_power(a0,t0,4,6,7));

	name.assign("R124");
	REFPROPname.assign("R124");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "deVries-ICR-1995";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R124Class::psat(double T)
{
 // Max error is  0.0279351679648 % between 120.0 and 395.424999 K

    const double t[]={0, 1.0, 1.1666666666666667, 2.1666666666666665, 4.0, 0.10300000000000001, 20.5};
    const double N[]={0, -8.2641394088060451, 1.9700771584419117, -0.68404895752825745, -3.8373863290050365, 0.0010940746131683551, -5.7602251358930756};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double R124Class::rhosatL(double T)
{
    // Maximum absolute error is 0.047253 % between 120.000001 K and 395.424990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.8333333333333333};
    const double N[] = {0, 0.54460398539271582, -4.4987463762100068, 251.49391057033526, -1397.1450617052847, 3683.0652090962581, -5428.9802833050198, 4444.2844628889843, -1678.1869995481597, 132.34762851690652};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=9; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}
double R124Class::rhosatV(double T)
{
    // Maximum absolute error is 0.238754 % between 120.000001 K and 395.424990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.8333333333333333, 2.1666666666666665};
    const double N[] = {0, -3.3190002541376638, 98.066603725937682, -1099.3702782361622, 6474.8477240028406, -23329.943477608645, 54044.724766470717, -80717.214878933504, 73706.624741330088, -33782.751910421401, 5453.4809807406809, -855.85806906737628};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=11; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}
