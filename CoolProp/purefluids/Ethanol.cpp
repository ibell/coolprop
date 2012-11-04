#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Ethanol.h"

EthanolClass::EthanolClass()
{
    double _n[]= {0, 1.14008942201e1, -3.95227128302e1, 4.13063408370e1, -1.88892923721e1, 4.72310314140, -7.78322827052e-3, 1.71707850032e-1, -1.53758307602, 1.42405508571, 1.32732097050e-1, -1.14231649761e-1, 3.27686088736e-6, 4.95699527725e-4, -7.01090149558e-5, -2.25019381648e-6, -2.55406026981e-1, -6.32036870646e-2, -3.14882729522e-2, 2.56187828185e-2, -3.08694499382e-2, 7.22046283076e-3, 2.99286406225e-3, 9.72795913095e-4};
	double _d[] = {0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 6, 7, 8, 8, 1, 3, 3, 6, 7, 8, 2, 7};
	double _t[] = {0, -0.5, 0, 0.5, 1.5, 2, 5, -0.5, 1, 2, 0, 2.5, 6, 2, 2, 4, 5, 3, 7, 5.5, 4, 1, 22, 23};
	double _l[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 4, 4};
    
    // Critical parameters
    crit.rho = 5.991*46.06844;
    crit.p = 6148.0;
    crit.T = 513.9;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 46.06844;
    params.Ttriple = 159.1;
	params.ptriple = 0.00000088;
    params.accentricfactor = 0.644;
    params.R_u = 8.314472;

    // Limits of EOS
	limits.Tmin = params.Ttriple;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

	// Residual part
    std::vector<double> n_v(_n, _n+sizeof(_n)/sizeof(double));
    std::vector<double> d_v(_d, _d+sizeof(_d)/sizeof(double));
    std::vector<double> t_v(_t, _t+sizeof(_t)/sizeof(double));
    std::vector<double> l_v(_l, _l+sizeof(_l)/sizeof(double));
    phirlist.push_back(new phir_power(n_v,d_v,t_v,l_v,1,23));

	// Ideal-gas part
	const double u0[]={0.0, 6.41129, 1.95989, 7.60084, 3.89583, 4.23238};
    const double v0[]={0.0, 0, 694, 1549, 2911, 4659};
	std::vector<double> u0_v(u0, u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0, v0+sizeof(v0)/sizeof(double));
	for (unsigned int i=1;i<v0_v.size();i++) { v0_v[i]/=crit.T; }

	double T0 = 273.15, 
		   p0 = 1.0, 
		   rho0 = p0/(params.R_u/params.molemass*T0),
		   m, c,
		   R_u = params.R_u,
		   R_ = params.R_u/params.molemass,
		   H0 = -28.200895800653548,
		   S0 = -0.13719644488861227,
		   tau0 = crit.T/T0, 
		   delta0 = rho0/crit.rho;
	
	// log(delta)+c+m*tau
	
	/// c is the constant term
	c=-S0/R_-1+log(tau0/delta0); /*<< from the leading term*/

	/// m multiplies the tau term in the leading term (slope)
	m=H0/(R_*crit.T); /*<< from the leading term */

	phi0list.push_back(new phi0_lead(c, m));
	phi0list.push_back(new phi0_logtau(-1));
	phi0list.push_back(new phi0_cp0_constant(u0_v[1],crit.T,T0));
	phi0list.push_back(new phi0_Planck_Einstein(u0_v,v0_v,2,5));

    EOSReference.assign("Dillon, H.E., and S. G. Penoncello, \"A Fundamental Equation for Calculation of the Thermodynamic Properties of Ethanol,\", International Journal of Thermophysics, Vol. 25, No. 2, March 2004.");
    TransportReference.assign("Using ECS");

    name.assign("Ethanol");
    //aliases.push_back(std::string("C2H6O")); 
    REFPROPname.assign("ETHANOL");
}
double EthanolClass::rhosatL(double T)
{
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 0.484219 %
	RHS = +0.205104*pow(theta,0.000000)-15.416759*pow(theta,0.714104)+15.243853*pow(theta,0.646579)+0.133914*pow(theta,8.236855)-0.198413*pow(theta,7.134261);
	rho = exp(RHS*reduce.T/T)*reduce.rho;
	return rho;
}
double EthanolClass::rhosatV(double T)
{
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 0.662657 %
	RHS = -0.486696*pow(theta,0.000000)-5.679398*pow(theta,0.656847)-3.890921*pow(theta,2.162034)-1.665655*pow(theta,3.649010)-0.477865*pow(theta,5.342953)+0.180460*pow(theta,6.372113)+0.563563*pow(theta,7.253718);
	rho = exp(RHS*reduce.T/T)*reduce.rho;
	return rho;
}
double EthanolClass::psat(double T)
{
	double theta = 1-T/reduce.T;
	double RHS,p;

	// Max error is 2.196317 %
	RHS = -11.227298*pow(theta,1.081329)-30.349599*pow(theta,3.427789)-8.178391*pow(theta,3.426262)-5.420801*pow(theta,3.380834)-3.447215*pow(theta,6.181072);
	p = exp(RHS)*reduce.p;
	return p;
}