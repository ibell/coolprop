#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "HFE143m.h"

HFE143mClass::HFE143mClass()
{
	double n[] = {0.0, 0.77715884e+01, -0.87042570e+01, -0.28095049e+00, 0.14540153e+00, 0.92291277e-02, -0.21416510e+00, 0.99475155e-01, 0.23247135e-01, -0.12873573e-01, -0.57366549e-01, 0.36504650e+00, -0.25433763e+00, -0.90896436e-01, 0.83503619e-01, 0.15477603e-01, -0.16641941e-01, 0.52410163e-02};
	double t[] = {0, 0.682, 0.851, 1.84, 1.87, 0.353, 3.92, 1.14, 0.104, 1.19, 6.58, 6.73, 7.99, 7.31, 7.45, 16.5, 24.8, 10.5};
	double d[] = {0, 1, 1, 1, 2, 5, 1, 3, 5, 7, 1, 2, 2, 3, 4, 2, 3, 5};
	double l[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3};

	//Critical parameters
	crit.rho = 4.648140744*100.04; //[kg/m^3]
	crit.p = PressureUnit(3635, UNIT_KPA); //[kPa]
	crit.T = 377.921; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 100.04;
	params.Ttriple = 240;
	params.accentricfactor = 0.28887136567003235;
	params.R_u = 8.314472;
	params.ptriple = 65.359393007;

	// Limits of EOS
	limits.Tmin = 240;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,l,1,17,18));

	double T0 = 273.15, 
		   p0 = 256.64, 
		   R_u = 8.314472,
		   R_= R_u/params.molemass,
		   rho0=p0/(R_*T0),
		   m,
		   c,
		   H0 = 40617.6, /// kJ/kmol
		   S0 = 174.758, /// kJ/kmol/K
		   tau0=crit.T/T0, 
		   delta0=rho0/crit.rho;
	
	// log(delta)+c+m*tau
	
	/// c is the constant term
	c=-S0/R_u-1+log(tau0/delta0);/*<< from the leading term*/

	/// m multiplies the tau term in the leading term (slope)
	m=H0/(R_u*crit.T); /*<< from the leading term */

	double N[] = {20.37, 0.2918, -1.950e-4, 4.650e-8};
	phi0list.push_back(new phi0_lead(c,m));
	phi0list.push_back(new phi0_logtau(-1.0));
	phi0list.push_back(new phi0_cp0_constant(N[0]/R_u,crit.T,T0));
	phi0list.push_back(new phi0_cp0_poly(N[1]/R_u,1,crit.T,T0));
	phi0list.push_back(new phi0_cp0_poly(N[2]/R_u,2,crit.T,T0));
	phi0list.push_back(new phi0_cp0_poly(N[3]/R_u,3,crit.T,T0));

	EOSReference.assign("Ryo Akasaka and Yohei Kayukawa \" A fundamental equation of state for trifluoromethyl methyl ether (HFE-143m) and its application to refrigeration cycle analysis \" International Journal of Refrigeration 35 (2012) 1003-1013 ");
	TransportReference.assign("Using ECS in fully predictive mode.");

	name.assign("HFE143m");
	aliases.push_back("HFE-143m");
	aliases.push_back("HFE143M");
	aliases.push_back("HFE-143M");
	aliases.push_back("RE143A");
	aliases.push_back("RE143a");
	REFPROPname.assign("RE143A");

	BibTeXKeys.EOS = "Akasaka-IJR-2012";
}
double HFE143mClass::psat(double T)
{
    const double ti[]={0, 1.0, 1.5, 2.5, 5};
    const double Ni[]={0, -7.44314,  1.69164, -2.27778, -4.09400};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double HFE143mClass::rhosatL(double T)
{
    const double ti[]={0, 0.33, 0.5, 1.5, 2.5};
    const double Ni[]={0, 1.20552, 1.33568, 0.0981486, 0.248917};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*(1+summer);
}
double HFE143mClass::rhosatV(double T)
{
    const double ti[]={0, 0.38, 1.24, 3.2, 6.9};
    const double Ni[]={0, -3.02576, -6.97239, -20.2601, -53.4441};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
