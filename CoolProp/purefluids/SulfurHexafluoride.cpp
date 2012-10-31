
#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "SulfurHexafluoride.h"

SulfurHexafluorideClass::SulfurHexafluorideClass()
{
	double _n [] = {0, 0.54958259132835, -0.87905033269396, -0.84656969731452, 0.27692381593529, -0.49864958372345e1, 0.48879127058055e1, 0.36917081634281e-1, 0.37030130305087e-3, 0.39389132911585e-1, 0.42477413690006e-3, -0.24150013863890e-1, 0.59447650642255e-1, -0.38302880142267, 0.32606800951983, -0.29955940562031e-1, -0.86579186671173e-1, 0.41600684707562e1, -0.41398128855814e1, -0.55842159922714, 0.56531382776891, 0.82612463415545e-2, -0.10200995338080e-1, -0.21662523861406e-1, 0.34650943893908e-1, -0.28694281385812e-1, 0.84007238998053e-2, -0.26969359922498, 0.90415215646344e1, -0.37233103557977e1, -0.27524670823704e4, 0.57711861697319e4, -0.30234003119748e4, 0.22252778435360e7, -0.23056065559032e7, 0.63918852944475e7, -0.60792091415592e7};
	double _l [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double _d [] = {0, 1, 1, 1, 2, 2, 2, 3, 3, 4, 6, 1, 2, 2, 2, 3, 6, 2, 2, 4, 4, 2, 2, 1, 3, 4, 1, 1, 4, 3, 4, 4, 4, 1, 1, 3, 3};
	double _t [] = {0, 0.125, 1.25, 1.875, 0.125, 1.5, 1.625, 1.5, 5.625, 0.625, 0.25, 6, 0.25, 4.75, 5.375, 5.875, 2, 5.875, 6, 5.625, 5.75, 0, 0.5, 4, 1, 3, 2, 4, 3, 4, 1, 2, 3, 3, 4, 3, 4};
	double _eta [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 10, 11, 25, 30, 30, 30, 30, 30, 30, 30, 30};
	double _beta [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 150, 150, 150, 150, 225, 300, 350, 350, 350, 350, 400, 400, 400, 400};
	double _gamma [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.13, 1.13, 1.13, 1.16, 1.19, 1.19, 1.16, 1.16, 1.16, 1.16, 1.22, 1.22, 1.22, 1.22};
	double _epsilon [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.85, 0.85, 0.85, 0.85, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

	std::vector<double> n_v(_n,_n+sizeof(_n)/sizeof(double));
	std::vector<double> d_v(_d,_d+sizeof(_d)/sizeof(double));
	std::vector<double> t_v(_t,_t+sizeof(_t)/sizeof(double));
	std::vector<double> l_v(_l,_l+sizeof(_l)/sizeof(double));
	std::vector<double> eta_v(_eta,_eta+sizeof(_eta)/sizeof(double));
	std::vector<double> epsilon_v(_epsilon,_epsilon+sizeof(_epsilon)/sizeof(double));
	std::vector<double> beta_v(_beta,_beta+sizeof(_beta)/sizeof(double));
	std::vector<double> gamma_v(_gamma,_gamma+sizeof(_gamma)/sizeof(double));

	//Critical parameters
	crit.rho = 742.3; //[kg/m^3]
	crit.p = 3754.983; //[kPa]
	crit.T = 318.7232; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 146.0554192;
	params.Ttriple = 223.555;
	params.ptriple = 231.429;
	params.accentricfactor = 0.21;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n_v,d_v,t_v,l_v,1,22));
	phirlist.push_back(new phir_gaussian( n_v,d_v,t_v,eta_v,epsilon_v,beta_v,gamma_v,23,36));

	double _theta [] ={0, 0, 0, 0, 1.617282065, 2.747115139, 4.232907175};
	std::vector<double> theta_v (_theta,_theta+sizeof(_theta)/sizeof(double));
	double _n0 [] ={0, 11.638611086, -6.392241811, 3.000000000, 3.661182320, 7.878851030, 3.459816790};
	std::vector<double> n0_v (_n0,_n0+sizeof(_n0)/sizeof(double));
	
	// lead term: log(delta)+c+m*tau
	phi0list.push_back(new phi0_lead(n0_v[1], n0_v[2]));
	phi0list.push_back(new phi0_logtau(n0_v[3]));
	phi0list.push_back(new phi0_Planck_Einstein(n0_v,theta_v,4,6));

	EOSReference.assign("Guder C., and W. Wagner, \"A Reference Equation of State for the Thermodynamic Properties of Sulfur Hexafluoride SF6 for Temperatures from the Melting Line to 625 K and Pressures up to 150 MPa,\" J. Phys. Chem. Ref. Data, Vol. 38, No. 1, 2009");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("SulfurHexafluoride");
	aliases.push_back("SF6");
	REFPROPname.assign("SF6");
}
double SulfurHexafluorideClass::rhosatL(double T)
{
	double theta = 1-T/reduce.T;
	double RHS,rho;
	RHS = +2.31174688*pow(theta,0.355) -1.12912486*pow(theta,1.0/2.0) -1.439347*pow(theta,8.0/6.0) +0.282489982*pow(theta,10.0/6.0);
	rho = exp(RHS*reduce.T/T)*reduce.rho;
	return rho;
}
double SulfurHexafluorideClass::rhosatV(double T)
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	RHS = 23.6806342*pow(theta,0.348)+0.513062232*pow(theta,1.0/6.0)-24.4706238*pow(theta,2.0/6.0)-4.6715244*pow(theta,4.0/6.0)-1.7536843*pow(theta,16.0/6.0)-6.65585369*pow(theta,34.0/6.0);
	rho = exp(RHS*reduce.T/T)*rhoc;
	return rho;
}
double SulfurHexafluorideClass::psat(double T)
{
	double theta = 1-T/reduce.T;
	double p = reduce.p*exp(reduce.T/T*(-7.09634642*pow(theta,1.0)+1.676662*pow(theta,1.5)-2.3921599*pow(theta,2.5)+5.86078302*pow(theta,4.0)-9.02978735*pow(theta,4.5)));
	return p;
}