
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
	crit.p = PressureUnit(3754.983, UNIT_KPA); //[kPa]
	crit.T = 318.7232; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 146.0554192;
	params.Ttriple = 223.555;
	params.ptriple = 231.42447394906830;
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
	aliases.push_back(std::string("SULFURHEXAFLUORIDE"));
	aliases.push_back("SF6");
	REFPROPname.assign("SF6");
	
	BibTeXKeys.EOS = "Guder-JPCRD-2009";
	BibTeXKeys.VISCOSITY = "QuinonesCisneros-JPCRD-2012";
	BibTeXKeys.CONDUCTIVITY = "Assael-JPCRD-2012";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "QuinonesCisneros-JPCRD-2012";
}
double SulfurHexafluorideClass::psat(double T)
{
    // Maximum absolute error is 0.020570 % between 223.555001 K and 318.723190 K
    const double t[]={0, 1, 2, 3, 6};
    const double N[]={0, -0.0087978014446881223, -6.9786635371802852, 1.2429125042516707, -2.5810326181235714};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double SulfurHexafluorideClass::rhosatL(double T)
{
    // Maximum absolute error is 0.144211 % between 223.555001 K and 318.723190 K
    const double t[] = {0, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667};
    const double N[] = {0, 14.005984855994196, -44.97667319923324, 74.392487901972885, -60.452498335429304, 19.684763024871778};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
	for (i=1; i<=5; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double SulfurHexafluorideClass::rhosatV(double T)
{
    // Maximum absolute error is 0.038534 % between 223.555001 K and 318.723190 K
    const double t[] = {0, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.8333333333333333};
    const double N[] = {0, -143.71556224434678, 955.84084676133841, -2974.3911165244631, 5116.1229629656846, -4874.2285072102668, 2143.8511434952561, -231.58521696825781};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;	
	for (i=1; i<=7; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}
double SulfurHexafluorideClass::viscosity_Trho(double T, double rho)
{
	double Tr = T/crit.T;

	// Dilute
	double d[] = {0.118561, -0.378103, 0.416428, -0.165295, 0.0245381};
	double eta_0 = d[0] + d[1]* pow(Tr,0.25) + d[2]*sqrt(Tr) + d[3]*pow(Tr,0.75) + d[4]*Tr; // mPa-s
	
	//// Initial density
	double c[] = {5.38783e-5, 1.63805e-6, -2.08160e-5};
	double psi1 = exp(1/Tr)-1;
	double psi2 = exp(1/Tr/Tr)-1;
	double ki = (c[0] + c[1]*psi1 + c[2]*psi2)/Tr;
	//double B = params.R_u*T/eta_0*ki;

	// Residual part
	double a0 = -6.87811e-4, b0 =  1.72737e-4, A0 =  9.99563e-8, B0 = -8.98256e-8;
	double a1 =  8.22661e-4, b1 = -2.02448e-4, A1 = -9.64167e-9, B1 = -8.49428e-8;
	double a2 = -3.54867e-4, b2 =  1.95952e-4, A2 = -7.54196e-9, B2 = 0;
	double C0 = -8.53432e-6, D0 = 0,           E0 = 0;
	double C1 =  1.14404e-5, D1 = 0,           E1 = -5.69402e-11;
	double C2 = -5.65762e-6, D2 = 2.27980e-11, E2 = 2.92190e-11;
	double ka = (a0 + a1*psi1 + a2*psi2)/Tr;
	double kr = (b0 + b1*psi1 + b2*psi2)/Tr;
	double kaa = (A0 + A1*psi1 + A2*psi2)/Tr/Tr/Tr;
	double krr = (B0 + B1*psi1 + B2*psi2)/Tr/Tr/Tr;
	double kii = (C0 + C1*psi1 + C2*psi2)/Tr/Tr/Tr;
	double krrr = (D0 + D1*psi1 + D2*psi2)/Tr;
	double kaaa = (E0 + E1*psi1 + E2*psi2)/Tr;

	double p = this->pressure_Trho(T,rho)/1e5; // Pa -> bar
	double pr = T*this->dpdT_Trho(T,rho)/1e5; // Pa-> bar
	double pa = p - pr;
	double pid = rho * R() * T / 1e5; // kPa -> bar
	double deltapr = pr - pid;

	double eta_f = ka*pa + kr*deltapr + ki*pid + kaa*pa*pa + krr*deltapr*deltapr + kii*pid*pid + krrr*pr*pr*pr + kaaa*pa*pa*pa;

	return (eta_0 + eta_f)/1000;
}
double SulfurHexafluorideClass::conductivity_Trho(double T, double rho)
{
	double sumresid = 0;
	double B1[] = {0, -2.83746e-2, 2.07472e-2, -5.57180e-3, 5.32890e-3, -1.61688e-3};
	double B2[] = {0, 3.52768e-2, -4.33053e-2, 5.12084e-2, -2.90262e-2, 5.98438e-3};
	// Assael JPCRD 2012
	double lambda_0 = (1461860 - 18539.4*T+77.7891*T*T+0.0241059*T*T*T)/(29661.7+505.67*T+T*T)/1000; //[W/m/K]

	for (int i = 1; i <= 5; i++)
	{
		sumresid += (B1[i]+B2[i]*T/crit.T)*pow(rho/crit.rho,i);
	}
	double lambda_r = sumresid; //[W/m/K]

	double lambda_c = this->conductivity_critical(T,rho,1/(3.5e-10)); //[W/m/K]

	return lambda_0 + lambda_r + lambda_c; //[W/m/K]
}
