#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Propylene.h"

PropyleneClass::PropyleneClass()
{
	double n[] = {0.0, 4.341002E-02, 1.136592E+00, -8.528611E-01, 5.216669E-01, -1.382953E+00, 1.214347E-01, -5.984662E-01, -1.391883E+00, -1.008434E+00, 1.961249E-01, -3.606930E-01, -2.407175E-03, 7.432121E-01, 1.475162E-01, -2.503391E-02, -2.734409E-01, 6.378889E-03, 1.502940E-02, -3.162971E-02, -4.107194E-02, -1.190241E+00};
	double t[] = {0, 1.0, 0.205, 0.56, 0.676, 1.0, 0.5, 1.0, 1.94, 2.0, 1.0, 2.66, 0.83, 1.6, 2.5, 3.0, 2.5, 2.72, 4.0, 4.0, 1.0, 4.0};
	double d[] = {0, 4., 1., 1., 2., 2., 3., 1., 1., 3., 2., 2., 8., 1., 1., 2., 3., 3., 2., 1., 2., 3.};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.07, 0.66, 1.2, 1.12, 1.47, 1.93, 3.3, 15.4, 6.0};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.77, 0.83, 0.607, 0.4, 0.66, 0.07, 3.1, 387.0, 41.0};
	double gamma[] = {0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.21, 1.08, 0.83, 0.56, 1.22, 1.81, 1.54, 1.12, 1.4};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.78, 0.82, 1.94, 0.69, 1.96, 1.3, 0.38, 0.91, 0.7};

	//Critical parameters
	crit.rho = 5.457*42.07974; //[kg/m^3]
	crit.p = PressureUnit(4555.0, UNIT_KPA); //[kPa]
	crit.T = 364.211; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 42.07974;
	params.Ttriple = 87.953;
	params.accentricfactor = 0.146;
	params.R_u = 8.314472;
	params.ptriple = 7.5e-7;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,12,22));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 13,21,22));

	const double n5 = 4.9916462, n6 = -0.1709449, n0 = 4;
	phi0list.push_back(new phi0_lead(n5,n6));
	phi0list.push_back(new phi0_logtau(n0-1));

	const double u0[] = {0, 324/crit.T, 973/crit.T, 1932/crit.T, 4317/crit.T};
	const double v0[] = {0, 1.544, 4.013, 8.923, 6.020};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	EOSReference.assign("Lemmon, E.W., Overhoff, U., McLinden, M.O., Wagner, W. to be submitted to J. Phys. Chem. Ref. Data, 2010. Coefficients from REFPROP with permission");
	TransportReference.assign("Using ECS in fully predictive mode.  Lennard-Jones parameters from Chichester NISTIR 6650\n\nSurface Tension: Mulero \"Recommended Correlations for the Surface Tension of Common Fluids\", 2012");

	name.assign("Propylene");
	aliases.push_back(std::string("propylene"));
	aliases.push_back(std::string("PROPYLENE"));
	REFPROPname.assign("PROPYLEN");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Lemmon-PROPYLENE-2013";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double PropyleneClass::psat(double T)
{
    // Maximum absolute error is 0.396942 % between 87.953001 K and 364.210990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-6.6958563845941557, 1.4311880756878577, -0.8392648926700903, -1.404370123028005, -1.29441159795851, -0.26643717806886891 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double PropyleneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.411874 % between 87.953001 K and 364.210990 K
    const double ti[]={0,0.30707507732809863, 1.6023345456599933, 1.7791635211003165, 2.0465240028604148, 5.1411742782570071};
    const double Ni[]={0,1.3786786344106812, -2.9774720934333989, 4.5456795310163045, -1.6615812178019818, 0.069846888007405475};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double PropyleneClass::rhosatV(double T)
{
    // Maximum absolute error is 1.016248 % between 87.953001 K and 364.210990 K
    const double ti[]={0,0.16707532108933132, 0.77455491593780279, 0.37832972634970236, 3.6472145456693612, 4.8837368709837277};
    const double Ni[]={0,-0.18885062045276654, -2.9307286022811478, -1.7316143840826113, -2.3756819694192388, -1.6778871008288518};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
double PropyleneClass::surface_tension_T(double T)
{
	// Mulero 2012
	return 0.05268*pow(1-T/crit.T,1.186);
}
