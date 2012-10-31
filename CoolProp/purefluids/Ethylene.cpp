
#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Ethylene.h"

EthyleneClass::EthyleneClass()
{
	double _n [] = {0,0.18617429100670e1, -0.30913708460844e1, -0.17384817095516, 0.80370985692840e-1, 0.23682707317354, 0.21922786610247e-1, 0.11827885813193, -0.21736384396776e-1, 0.44007990661139e-1, 0.12554058863881, -0.13167945577241, -0.52116984575897e-2, 0.15236081265419e-3, -0.24505335342756e-4, 0.28970524924022, -0.18075836674288, 0.15057272878461, -0.14093151754458, 0.22755109070253e-1, 0.14026070529061e-1, 0.61697454296214e-2, -0.41286083451333e-3, 0.12885388714785e-1, -0.69128692157093e-1, 0.10936225568483, -0.81818875271794e-2, -0.56418472117170e-1, 0.16517867750633e-2, 0.95904006517001e-2, -0.26236572984886e-2, -0.50242414011355e2, 0.74846420119299e4, -0.68734299232625e4, -0.93577982814338e3, 0.94133024786113e3};
	double _d [] = {0, 1, 1, 1, 2, 2, 4, 1, 1, 3, 4, 5, 7, 10, 11, 1, 1, 2, 2, 4, 4, 6, 7, 4, 5, 6, 6, 7, 8, 9, 10, 2, 2, 2, 3, 3};
	double _t [] = {0, 0.5, 1, 2.5, 0, 2, 0.5, 1, 4, 1.25, 2.75, 2.25, 1, 0.75, 0.5, 2.5, 3.5, 4, 6, 1.5, 5, 4.5, 15, 20, 23, 22, 29, 19, 15, 13, 10, 1, 0, 1, 2, 3};
	double _l [] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0};
	double _eta [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 25, 25, 25, 25};
	double _beta [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 325, 300, 300, 300, 300};
	double _gamma [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.16, 1.19, 1.19, 1.19, 1.19};
	double _epsilon [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1};

	std::vector<double> n_v(_n,_n+sizeof(_n)/sizeof(double));
	std::vector<double> d_v(_d,_d+sizeof(_d)/sizeof(double));
	std::vector<double> t_v(_t,_t+sizeof(_t)/sizeof(double));
	std::vector<double> l_v(_l,_l+sizeof(_l)/sizeof(double));
	std::vector<double> eta_v(_eta,_eta+sizeof(_eta)/sizeof(double));
	std::vector<double> epsilon_v(_epsilon,_epsilon+sizeof(_epsilon)/sizeof(double));
	std::vector<double> beta_v(_beta,_beta+sizeof(_beta)/sizeof(double));
	std::vector<double> gamma_v(_gamma,_gamma+sizeof(_gamma)/sizeof(double));

	//Critical parameters
	crit.rho = 214.24; //[kg/m^3]
	crit.p = 5041.8; //[kPa]
	crit.T = 282.35; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 28.05376;
	params.Ttriple = 103.989;
	params.ptriple = 0.12265;
	params.accentricfactor = 0.0866;
	params.R_u = 8.31451;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n_v,d_v,t_v,l_v,1,30));
	phirlist.push_back(new phir_gaussian( n_v,d_v,t_v,eta_v,epsilon_v,beta_v,gamma_v,31,35));

	double _theta [] ={0, 0, 0, 0, 4.43266896, 5.74840149,	7.8027825, 15.5851154};
	std::vector<double> theta_v (_theta,_theta+sizeof(_theta)/sizeof(double));
	double _n0 [] ={0, 8.68815523, -4.47960564, 3.0, 2.49395851, 3.0027152, 2.5126584, 3.99064217};
	std::vector<double> n0_v (_n0,_n0+sizeof(_n0)/sizeof(double));
	
	// lead term: log(delta)+c+m*tau
	phi0list.push_back(new phi0_lead(n0_v[1], n0_v[2]));
	phi0list.push_back(new phi0_logtau(n0_v[3]));
	phi0list.push_back(new phi0_Planck_Einstein(n0_v,theta_v,4,7));

	EOSReference.assign("Smukala, J., R. Span, and W. Wagner, \"New Equation of State for Ethylene Covering the Fluid Region for Temperatures From the Melting Line to 450 K at Pressures up to 300 MPa,\" J. Phys. Chem. Ref. Data, Vol. 29, No. 5, 2000");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("Ethylene");
	REFPROPname.assign("ETHYLENE");
}
double EthyleneClass::rhosatL(double T)
{
	double theta = 1-T/reduce.T;
	double RHS,rho;
	RHS = +1.8673079*pow(theta,0.353) -0.61533892*pow(theta,1.0/2.0) -0.058973772*pow(theta,8.0/6.0) +0.10744720*pow(theta,12.0/6.0);
	rho = exp(RHS)*reduce.rho;
	return rho;
}
double EthyleneClass::rhosatV(double T)
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	RHS = -1.9034556*pow(theta,0.349)-0.75837929*pow(theta,4.0/6.0)-3.7717969*pow(theta,1.0)-8.7478586*pow(theta,14.0/6.0)-23.885296*pow(theta,29.0/6.0)-54.197979*pow(theta,56.0/6.0);
	rho = exp(RHS*reduce.T/T)*rhoc;
	return rho;
}
double EthyleneClass::psat(double T)
{
	double theta = 1-T/reduce.T;
	double p = reduce.p*exp(reduce.T/T*(-6.3905741*pow(theta,1.0)+1.4060338*pow(theta,1.5)-1.6589923*pow(theta,2.5)+1.0278028*pow(theta,3.0)-2.5071716*pow(theta,4.5)));
	return p;
}