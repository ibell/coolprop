#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Cyclopentane.h"

CyclopentaneClass::CyclopentaneClass()
{
	double n[] = {0.0, 0.0536938, 1.60394, -2.41244, -0.474009, 0.203482, -0.965616, -0.344543, 0.353975, -0.231373, -0.0379099, 0.867586, -0.381827, -0.108741, -0.0976984};
	double t[] = {0, 1.0, 0.29, 0.8, 1.14, 0.5, 2.0, 1.5, 1.0, 3.36, 0.95, 1.0, 2.5, 2.5, 1.5};
	double d[] = {0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3, 3};
	double c[] = {0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.82, 1.19, 0.79, 1.52};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.15, 1.61, 0.66, 2.72};
	double gamma[] = {0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.08, 0.36, 0.09, 1.48};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.68, 0.97, 0.84, 0.66};

	//Critical parameters
	crit.rho = 3.82*70.1329; //[kg/m^3]
	crit.p = PressureUnit(4571.2, UNIT_KPA); //[kPa]
	crit.T = 511.72; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 70.1329;
	params.Ttriple = 179.7;
	params.accentricfactor = 0.20111257338059896;
	params.R_u = 8.314472;
	params.ptriple = 0.0088542726061155586;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,10,15));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 11,14,15));

	const double n5 = 3.2489131288, n6 = 2.6444166315, n0 = 1.96;
	phi0list.push_back(new phi0_lead(n5,n6));
	phi0list.push_back(new phi0_logtau(n0-1));

	const double u0[] = {0, 120/crit.T, 1300/crit.T, 2700/crit.T, 5300/crit.T};
	const double v0[] = {0, 3.34, 18.6, 13.9, 4.86};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	EOSReference.assign("Holger Gedanitz, Maria J. Davila, Eric W. Lemmon \" Speed of sound measurements and a fundamental equation of state for cyclopentane\" Preprint provided by Eric Lemmon");
	TransportReference.assign("Using ECS in fully predictive mode.");

	name.assign("Cyclopentane");
	aliases.push_back(std::string("CycloPentane"));
	aliases.push_back(std::string("cyclopentane"));
	aliases.push_back(std::string("CYCLOPENTANE"));
	REFPROPname.assign("CYCLOPEN");

	BibTeXKeys.EOS = "Gedanitz-PREPRINT-2013";

}

double CyclopentaneClass::psat(double T)
{
    // Maximum absolute error is 0.044586 % between 179.722001 K and 511.689990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-6.5559892604046448, 0.1683169483823736, 0.98984878102513452, -3.7731391793932954, 0.062343331455625178, -0.66025999749599096 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double CyclopentaneClass::rhosatL(double T)
{
   // Max error is  0.0537787674075 % between 179.7 and 511.7199 K
    const double ti[]={0, 0.127, 0.137, 0.14100000000000001, 0.14700000000000002, 0.3505, 0.3685, 0.38699999999999996, 0.39849999999999997, 1.8333333333333333, 3.8333333333333335};
    const double Ni[]={0, -2949907.6373172961, 26157071.053032953, -34291932.874714434, 11119726.461818002, -2256287.1652462273, 5780550.4968119375, -6150540.6041890942, 2591325.8956157742, -3.4860599305100841, 0.71396246842304301};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i <= 10;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*(summer+1);
}
double CyclopentaneClass::rhosatV(double T)
{
	// Max error is  0.130802597771 % between 179.7 and 511.7199 K
    const double ti[]={0, 0.053000000000000005, 0.356, 0.3605, 0.364, 0.378, 0.39149999999999996, 0.39549999999999996, 0.8333333333333334, 6.333333333333333, 7.5};
    const double Ni[]={0, 12.066536979526935, -561852563.38938332, 2137733277.8963461, -1918813779.8001397, 552902395.41400063, -465184138.2543214, 255214994.8445133, -205.13360552504494, -12.344541888181883, 10.571728743217056};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1; i<=10; i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
