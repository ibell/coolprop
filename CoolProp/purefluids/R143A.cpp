
#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "R143A.h"

R143AClass::R143AClass()
{
	double n [] = {0, 7.7736443, -8.7018500, -0.27779799, 0.14609220, 0.0089581616, -0.20552116, 0.10653258, 0.023270816, -0.013247542, -0.042793870, 0.36221685, -0.25671899, -0.092326113, 0.083774837, 0.017128445, -0.017256110, 0.0049080492};
	double d [] = {0, 1, 1, 1, 2, 5, 1, 3, 5, 7, 1, 2, 2, 3, 4, 2, 3, 5};
	double t [] = {0, 0.67, 0.833, 1.7, 1.82, 0.35, 3.9, 0.95, 0.0, 1.19, 7.2, 5.9, 7.65, 7.5, 7.45, 15.5, 22.0, 19.0};
	double l [] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3};

	//Critical parameters
	crit.rho = 5.12845*84.041; //[kg/m^3]
	crit.p = PressureUnit(3761, UNIT_KPA); //[kPa]
	crit.T = 345.857; //[K]
	crit.v = 1/crit.rho;

	// Other fluid parameters
	params.molemass = 84.041;
	params.Ttriple = 161.34;
	params.ptriple = 1.077;
	params.accentricfactor = 0.261489646224;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,l,1,17,18));
	
	// lead term: log(delta)+c+m*tau
	phi0list.push_back(new phi0_lead(5.903087, 7.307253));
	phi0list.push_back(new phi0_logtau(-1));
	phi0list.push_back(new phi0_cp0_poly(1.0578,0.33,crit.T,273.15));
	phi0list.push_back(new phi0_Planck_Einstein(4.4402, 1791/crit.T));
	phi0list.push_back(new phi0_Planck_Einstein(3.7515, 823/crit.T));

	EOSReference.assign("Eric W. Lemmon and Richard T. Jacobsen, \"An International Standard Formulation for the Thermodynamic Properties of 1,1,1-Trifluoroethane (HFC-143a) for Temperatures From 161 to 450 K and Pressures to 50 MPa,\" J. Phys. Chem. Ref. Data, Vol. 29, No. 4, 2000");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("R143a");
	aliases.push_back("R143A");
	REFPROPname.assign("R143A");

	BibTeXKeys.EOS = "LemmonJacobsen-JPCRD-2000";
	BibTeXKeys.ECS_LENNARD_JONES = "McLinden-IJR-2000";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R143AClass::psat(double T)
{
    const double ti[]={0,1.0,1.5,2.0,3.5,5.5};
    const double Ni[]={0,-7.3526, 1.9162,-1.2203,-2.0532,-1.7354};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(crit.T/T*summer);
}
double R143AClass::rhosatL(double T)
{
    const double ti[]={0,1.0/3.0,2.0/3.0,8.0/3.0};
    const double Ni[]={0,1.795,0.8709,0.3116};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=3;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return crit.rho*(summer+1);
}
double R143AClass::rhosatV(double T)
{
    const double ti[]={0,0.39,4.0/3.0,11.0/3.0,23.0/3.0};
    const double Ni[]={0,-3.072,-8.061,-23.74,-61.37};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
