#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "R1233zd(E).h"

R1233zdEClass::R1233zdEClass()
{
	double n[] = {0.0, 0.03920261, 1.639052, -1.997147, -0.6603372, 0.1498682, -1.408791, -0.7920426, 0.8549678, -0.5301920, -0.01408562, 1.335117, -0.5441797, -0.05862723, -0.04123614, -0.6619106};
	double t[] = {0, 1.0, 0.24, 0.83, 1.17, 0.6, 2.2, 2.88, 1.1, 2.0, 1.07, 1.27, 1.94, 2.0, 1.5, 1.0};
	double d[] = {0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3, 2, 3};
	double c[] = {0, 0., 0., 0., 0., 0., 2., 2., 1., 2., 1., 0, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.215, 1.5, 1.1, 2.52, 4.55};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.27, 0.82, 0.94, 20, 32};
	double gamma[] = {0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.32, 0.82, 0.66, 0.66, 1.39};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.77, 0.976, 1.08, 0.62, 0.61};

	//Critical parameters
	crit.rho = 3.67*130.4961896 ; //[kg/m^3]
	crit.p = PressureUnit(3570.9, UNIT_KPA); //[kPa]
	crit.T = 438.75; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 130.4961896;
	params.Ttriple = 195.15;
	params.accentricfactor = 0.305;
	params.R_u = 8.314472;
	params.ptriple = 0.25;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,10,16));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 11,15,16));

	const double n5 = 0, n6 = 0, n0 = 0;
	phi0list.push_back(new phi0_lead(n5,n6));
	phi0list.push_back(new phi0_logtau(n0-1));

	const double u0[] = {0, 400/crit.T, 1900/crit.T};
	const double v0[] = {0, 8.962, 11.94};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));
	phi0list.push_back(new phi0_cp0_poly(4,0,crit.T,298.15));
	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,2));

	name.assign("R1233zd(E)");
	aliases.push_back(std::string("R1233zdE"));
	aliases.push_back(std::string("R1233ZDE"));
	aliases.push_back(std::string("R1233ZD(E)"));
	REFPROPname.assign("R1233zd");
}

double R1233zdEClass::psat(double T)
{
    const double ti[]={0,1.0, 1.5, 2.0, 4.3, 14.0};
    const double Ni[]={0,-7.6021, 2.3265, -1.9771, -4.8451, -4.8762};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double R1233zdEClass::rhosatL(double T)
{
    const double ti[]={0,0.355, 0.9, 3.5, 8.0};
    const double Ni[]={0,2.13083, 0.583568, 0.247871, 0.472173};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
	double rho = reduce.rho*(summer+1);
    return rho;
}
double R1233zdEClass::rhosatV(double T)
{
    const double ti[]={0, 0.397, 1.2, 3.1, 6.6, 15.0};
    const double Ni[]={0, -3.0152, -6.5621, -19.427, -62.650, -181.64  };
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
	double rho = reduce.rho*exp(summer);
    return rho;
}
