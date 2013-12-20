#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "DMC.h"

static char EOSstr [] = "Yong Zhou, Jiangtao Wu and Eric W. Lemmon, \"Thermodynamic Properties of Dimethyl Carbonate\", J. Phys. Chem. Ref. Data, Vol. 40, No. 4, 2011";

DimethylCarbonateClass::DimethylCarbonateClass()
{
	double n[] = {0.0, 0.00052683187, 1.353396, -2.649283, -0.2785412, 0.1742554, 0.031606252, 0.399866, 1.178144, -0.0235281, -1.015, -0.7880436, -0.12696, 1.2198, -0.4883, -0.0033293, -0.0035387, -0.51172, -0.16882};
	double t[] = {0, 1, 0.227, 1.05, 1.06, 0.5, 0.78, 1.3, 1.347, 0.706, 2, 2.5, 4.262, 1, 2.124, 0.4, 3.5, 0.5, 2.7};
	double d[] = {0, 5, 1, 1, 2, 3, 4, 1, 2, 7, 1, 2, 3, 1, 1, 2, 2, 3, 3};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 0, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.9667, 1.5154, 1.0591, 1.6642, 12.4856, 0.9662};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.240, 0.821, 15.45, 2.210, 437.0, 0.743};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.2827, 0.4317, 1.1217, 1.1871, 1.1243, 0.4203};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6734, 0.9239, 0.8636, 1.0507, 0.8482, 0.7522};

	//Critical parameters
	crit.rho = 4*90.0779; //[kg/m^3]
	crit.p = PressureUnit(4908.8, UNIT_KPA); //[kPa]
	crit.T = 557; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 90.0779;
	params.Ttriple = 277.06;
	params.accentricfactor = 0.346;
	params.R_u = 8.314472;
	params.ptriple = 2.2265237265924869;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,12,19));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 13,18,19));

	const double n5 = 4.9916462, n6 = -0.1709449, n0 = 9.28421;
	phi0list.push_back(new phi0_lead(n5,n6));
	phi0list.push_back(new phi0_logtau(n0-1));

	const double u0[] = {0, 21/crit.T, 1340/crit.T, 1672/crit.T, 7395/crit.T};
	const double v0[] = {0, 1.48525, 0.822585, 16.2453, 1.15925};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("DimethylCarbonate");
	aliases.push_back(std::string("DMC"));
	aliases.push_back(std::string("dimethylcarbonate"));
	aliases.push_back(std::string("DIMETHYLCARBONATE"));
	REFPROPname.assign("DMC");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Zhou-JPCRD-2011";
}

double DimethylCarbonateClass::psat(double T)
{
    // Maximum absolute error is 0.015238 % between 277.060001 K and 556.999990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-8.3202660729706679, 3.4264854309711628, -3.5672563248732438, -0.4374159102683291, -3.895322642705378, 1.5398890451904521 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double DimethylCarbonateClass::rhosatL(double T)
{
    // Maximum absolute error is 0.652451 % between 277.060001 K and 556.999990 K
    const double ti[]={0,0.33284736596795433, 1.1366445658032931, 2.7906891264436053, 1.7799430494935189, 2.0875761898754157};
    const double Ni[]={0,1.7636989589272027, -2.9872459633577173, 2.7945804747539671, 12.18823328946845, -12.277305677322932};
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
double DimethylCarbonateClass::rhosatV(double T)
{
	// Max error is  0.152889150153 % between 277.06 and 556.999999 K
    const double ti[]={0, 0.384, 0.38649999999999995, 0.39049999999999996, 1.1666666666666667, 1.6666666666666667, 5.0};
    const double Ni[]={0, -8046.3851723522694, 13270.96908464885, -5230.1789590545313, 3.5069500792786417, -5.1611358518618955, -5.7960154390551759};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
