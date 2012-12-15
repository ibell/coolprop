
#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Alkanes.h"

MethaneClass::MethaneClass()
{
	double _n[] = {0, 4.367901028E-02, 6.709236199E-01, -1.765577859E+00, 8.582330241E-01, -1.206513052E+00, 5.120467220E-01, -4.000010791E-04, -1.247842423E-02, 3.100269701E-02, 1.754748522E-03, -3.171921605E-06, -2.240346840E-06, 2.947056156E-07, 1.830487909E-01, 1.511883679E-01, -4.289363877E-01, 6.894002446E-02, -1.408313996E-02, -3.063054830E-02, -2.969906708E-02, -1.932040831E-02, -1.105739959E-01, 9.952548995E-02, 8.548437825E-03, -6.150555662E-02, -4.291792423E-02, -1.813207290E-02, 3.445904760E-02, -2.385919450E-03, -1.159094939E-02, 6.641693602E-02, -2.371549590E-02, -3.961624905E-02, -1.387292044E-02, 3.389489599E-02, -2.927378753E-03, 9.324799946E-05, -6.287171518E+00, 1.271069467E+01, -6.423953466E+00};
	double _d[] = {0, 1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 8, 9, 10, 1, 1, 1, 2, 4, 5, 6, 1, 2, 3, 4, 4, 3, 5, 5, 8, 2, 3, 4, 4, 4, 5, 6, 2, 0, 0, 0};
	double _t[] = {0, -0.5, 0.5, 1, 0.5, 1, 1.5, 4.5, 0, 1, 3, 1, 3, 3, 0, 1, 2, 0, 0, 2, 2, 5, 5, 5, 2, 4, 12, 8, 10, 10, 10, 14, 12, 18, 22, 18, 14, 2, 0, 1, 2};
	double _l[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2};
	double _eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 20, 40, 40, 40};
	double _epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1};
	double _beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 200, 250, 250, 250};
	double _gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.07, 1.11, 1.11, 1.11};

	std::vector<double> n_v(_n,_n+sizeof(_n)/sizeof(double));
	std::vector<double> d_v(_d,_d+sizeof(_d)/sizeof(double));
	std::vector<double> t_v(_t,_t+sizeof(_t)/sizeof(double));
	std::vector<double> l_v(_l,_l+sizeof(_l)/sizeof(double));
	std::vector<double> eta_v(_eta,_eta+sizeof(_eta)/sizeof(double));
	std::vector<double> epsilon_v(_epsilon,_epsilon+sizeof(_epsilon)/sizeof(double));
	std::vector<double> beta_v(_beta,_beta+sizeof(_beta)/sizeof(double));
	std::vector<double> gamma_v(_gamma,_gamma+sizeof(_gamma)/sizeof(double));

	//Critical parameters
	crit.rho = 162.66; //[kg/m^3]
	crit.p = 4599.2; //[kPa]
	crit.T = 190.564; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 16.0428;
	params.Ttriple = 90.6941;
	params.ptriple = 11.696;
	params.accentricfactor = 0.01142;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n_v,d_v,t_v,l_v,1,36));
	phirlist.push_back(new phir_gaussian( n_v,d_v,t_v,eta_v,epsilon_v,beta_v,gamma_v,37,40));

	double _theta [] ={0,0,0,0,3.400432401,10.26951575,20.43932747,29.93744884,79.13351945};
	std::vector<double> theta_v (_theta,_theta+sizeof(_theta)/sizeof(double));
	double _n0 [] ={0,9.91243972,-6.33270087,3.0016,0.008449,4.6942,3.4865,1.6572,1.4115};
	std::vector<double> n0_v (_n0,_n0+sizeof(_n0)/sizeof(double));
	
	// lead term: log(delta)+c+m*tau
	phi0list.push_back(new phi0_lead(n0_v[1], n0_v[2]));
	phi0list.push_back(new phi0_logtau(n0_v[3]));
	phi0list.push_back(new phi0_Planck_Einstein(n0_v,theta_v,4,8));

	EOSReference.assign("Setzmann, U. and Wagner, W., \"A New Equation of State and Tables of Thermodynamic Properties for Methane Covering the Range from the Melting Line to 625 K at Pressures up to 1000 MPa,\" J. Phys. Chem. Ref. Data, 20(6):1061-1151, 1991.");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("Methane");
	aliases.push_back("CH4");
	REFPROPname.assign("METHANE");
}
double MethaneClass::rhosatL(double T)
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	RHS = +1.9906389*pow(theta,0.354)
		  -0.78756197*pow(theta,1.0/2.0)
		  +0.0369976723*pow(theta,5.0/2.0);
	rho = exp(RHS)*rhoc;
	return rho;
}
double MethaneClass::rhosatV(double T)
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	RHS = -1.8802840*pow(theta,0.354)
		  -2.8526531*pow(theta,5.0/6.0)
		  -3.0006480*pow(theta,3.0/2.0)
		  -5.2511690*pow(theta,5.0/2.0)
		  -13.191859*pow(theta,25.0/6.0)
		  -37.553961*pow(theta,47.0/6.0);
	rho = exp(RHS)*rhoc;
	return rho;
}
double MethaneClass::psat(double T)
{
	double pc = reduce.p;
	double theta = 1-T/reduce.T;
	double RHS,p;

	RHS = -6.036219*pow(theta,1.0)
		  +1.409353*pow(theta,1.5)
		  -0.4945199*pow(theta,2.0)
		  -1.443048*pow(theta,4.5);
	p = exp(reduce.T/T*RHS)*pc;
	return p;
}

EthaneClass::EthaneClass()
{

double n[] = 
{0.0,
0.83440745735241, //[1]
-1.4287360607171, //[2]
0.34430242210927, //[3]
-0.42096677920265, //[4]
0.012094500886549, //[5]
-0.57976201597341, //[6]
-0.033127037870838, //[7]
-0.1175165489413, //[8]
-0.11160957833067, //[9]
0.062181592654406, //[10]
0.098481795434443, //[11]
-0.098268582682358, //[12]
-0.00023977831007049, //[13]
0.00069885663328821, //[14]
0.000019665987803305, //[15]
-0.014586152207928, //[16]
0.046354100536781, //[17]
0.0060764622180645, //[18]
-0.0026447330147828, //[19]
-0.042931872689904, //[20]
0.0029987786517263, //[21]
0.005291933517501, //[22]
-0.0010383897798198, //[23]
-0.054260348214694, //[24]
-0.21959362918493, //[25]
0.35362456650354, //[26]
-0.12477390173714, //[27]
0.18425693591517, //[28]
-0.16192256436754, //[29]
-0.082770876149064, //[30]
0.050160758096437, //[31]
0.0093614326336655, //[32]
-0.00027839186242864, //[33]
0.000023560274071481, //[34]
0.0039238329738527, //[35]
-0.00076488325813618, //[36]
-0.004994430444073, //[37]
0.0018593386407186, //[38]
-0.00061404353331199, //[39]
-0.0023312179367924, //[40]
0.002930104790876, //[41]
-0.00026912472842883, //[42]
184.13834111814, //[43]
-10.397127984854, //[44]
};

double d[] =
{0,
1, //[1]
1, //[2]
2, //[3]
2, //[4]
4, //[5]
1, //[6]
1, //[7]
2, //[8]
2, //[9]
3, //[10]
6, //[11]
6, //[12]
7, //[13]
9, //[14]
10, //[15]
2, //[16]
4, //[17]
4, //[18]
5, //[19]
5, //[20]
6, //[21]
8, //[22]
9, //[23]
2, //[24]
3, //[25]
3, //[26]
3, //[27]
4, //[28]
4, //[29]
5, //[30]
5, //[31]
6, //[32]
11, //[33]
14, //[34]
3, //[35]
3, //[36]
4, //[37]
8, //[38]
10, //[39]
1, //[40]
1, //[41]
3, //[42]
3, //[43]
2 //[44]
};

double t[] =
{0,
0.25, //[1]
1, //[2]
0.25, //[3]
0.75, //[4]
0.75, //[5]
2, //[6]
4.25, //[7]
0.75, //[8]
2.25, //[9]
3, //[10]
1, //[11]
1.25, //[12]
2.75, //[13]
1, //[14]
2, //[15]
2.5, //[16]
5.5, //[17]
7, //[18]
0.5, //[19]
5.5, //[20]
2.5, //[21]
4, //[22]
2, //[23]
10, //[24]
16, //[25]
18, //[26]
20, //[27]
14, //[28]
18, //[29]
12, //[30]
19, //[31]
7, //[32]
15, //[33]
9, //[34]
26, //[35]
28, //[36]
28, //[37]
22, //[38]
13, //[39]
0, //[40]
3, //[41]
3, //[42]
0, //[43]
3, //[44]
};

double l[] =
{0,
0, //[1]
0, //[2]
0, //[3]
0, //[4]
0, //[5]
1, //[6]
1, //[7]
1, //[8]
1, //[9]
1, //[10]
1, //[11]
1, //[12]
1, //[13]
1, //[14]
1, //[15]
2, //[16]
2, //[17]
2, //[18]
2, //[19]
2, //[20]
2, //[21]
2, //[22]
2, //[23]
3, //[24]
3, //[25]
3, //[26]
3, //[27]
3, //[28]
3, //[29]
3, //[30]
3, //[31]
3, //[32]
3, //[33]
3, //[34]
4, //[35]
4, //[36]
4, //[37]
4, //[38]
4, //[39]
0, //[40]
0, //[41]
0, //[42]
0, //[43]
0, //[44]
};

double eta [] =
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //[indices 0-39]
15, //[40]
15, //[41]
15, //[42]
20, //[43]
20, //[44]
};

double epsilon [] =
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //[indices 0-39]
1, //[40]
1, //[41]
1, //[42]
1, //[43]
1, //[44]
};

double beta [] =
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //[indices 0-39]
150, //[40]
150, //[41]
150, //[42]
275, //[43]
400, //[44]
};

double gamma [] =
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //[indices 0-39]
1.05, //[40]
1.05, //[41]
1.05, //[42]
1.22, //[43]
1.16, //[44]
};


	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> d_v(d,d+sizeof(d)/sizeof(double));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(l,l+sizeof(l)/sizeof(double));
	std::vector<double> eta_v(eta,eta+sizeof(eta)/sizeof(double));
	std::vector<double> epsilon_v(epsilon,epsilon+sizeof(epsilon)/sizeof(double));
	std::vector<double> beta_v(beta,beta+sizeof(beta)/sizeof(double));
	std::vector<double> gamma_v(gamma,gamma+sizeof(gamma)/sizeof(double));

	//Critical parameters
	crit.rho = 206.18; //[kg/m^3]
	crit.p = 4872.2; //[kPa]
	crit.T = 305.322; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 30.06904;
	params.Ttriple = 90.368;
	params.ptriple = 0.00114240920349;
	params.accentricfactor = 0.099;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n_v,d_v,t_v,l_v,1,39));
	phirlist.push_back(new phir_gaussian( n_v,d_v,t_v,eta_v,epsilon_v,beta_v,gamma_v,40,44));

	double _theta [] ={0,0,0,0,1.409105233,4.009917071,6.596709834,13.97981027};
	std::vector<double> theta_v (_theta,_theta+sizeof(_theta)/sizeof(double));
	double _n0 [] ={0,9.212802589,-4.68224855,3.003039265,1.117433359,3.467773215,6.94194464,5.970850948};
	std::vector<double> n0_v (_n0,_n0+sizeof(_n0)/sizeof(double));
	
	// lead term: log(delta)+c+m*tau
	phi0list.push_back(new phi0_lead(n0_v[1], n0_v[2]));
	phi0list.push_back(new phi0_logtau(n0_v[3]));
	phi0list.push_back(new phi0_Planck_Einstein(n0_v,theta_v,4,7));

	EOSReference.assign("Buecker, D. and Wagner, W. \"A Reference Equation of State for the Thermodynamic Properties of Ethane for Temperatures from the Melting Line to 675 K and Pressures up to 900 MPa,\" J. Phys. Chem. Ref. Data, 35(1):205-266, 2006.");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("Ethane");
	REFPROPname.assign("ETHANE");
}
double EthaneClass::rhosatL(double T)
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	RHS = +1.56138026*pow(theta,0.329)
		  -0.381552776*pow(theta,4.0/6.0)
		  +0.0785372040*pow(theta,8.0/6.0)
		  +0.0370315089*pow(theta,19.0/6.0);
	rho = exp(RHS)*rhoc;
	return rho;
}
double EthaneClass::rhosatV(double T)
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	RHS = -1.89879145*pow(theta,0.346)
		  -3.65459262*pow(theta,5.0/6.0)
		  +0.850562745*pow(theta,1.0)
		  +0.363965487*pow(theta,2)
		  -1.50005943*pow(theta,3)
		  -2.26690389*pow(theta,5);
	rho = exp(RHS*reduce.T/T)*rhoc;
	return rho;
}
double EthaneClass::psat(double T)
{
	double pc = reduce.p;
	double theta = 1-T/reduce.T;
	double RHS,p;

	RHS = -6.48647577*pow(theta,1.0)
		  +1.47010078*pow(theta,1.5)
		  -1.66261122*pow(theta,2.5)
		  +3.57898378*pow(theta,3.5)
		  -4.79105705*pow(theta,4.0);
	p = exp(reduce.T/T*RHS)*pc;
	return p;
}



nButaneClass::nButaneClass()
{
	double _n [] = {0, 2.5536998241635E+00, -4.4585951806696E+00, 8.2425886369063E-01, 1.1215007011442E-01, -3.5910933680333E-02, 1.6790508518103E-02, 3.2734072508724E-02, 9.5571232982005E-01, -1.0003385753419E+00, 8.5581548803855E-02, -2.5147918369616E-02, -1.5202958578918E-03, 4.7060682326420E-03, -9.7845414174006E-02, -4.8317904158760E-02, 1.7841271865468E-01, 1.8173836739334E-02, -1.1399068074953E-01, 1.9329896666669E-02, 1.1575877401010E-03, 1.5253808698116E-04, -4.3688558458471E-02, -8.2403190629989E-03, -2.8390056949441E-02, 1.4904666224681E-03};
	double _l[]= {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 0, 0};
	double _d[]= {0, 1, 1, 1, 2, 3, 4, 4, 1, 1, 2, 7, 8, 8, 1, 2, 3, 3, 4, 5, 5, 10, 2, 6, 1, 2};
	double _t[]= {0, 0.5, 1, 1.5, 0, 0.5, 0.5, 0.75, 2, 2.5, 2.5, 1.5, 1, 1.5, 4, 7, 3, 7, 3, 1, 6, 0, 6, 13, 2, 0};
	double _eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10};
	double _epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.85, 1};
	double _beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 150, 200};
	double _gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.16, 1.13};
	std::vector<double> n_v(_n,_n+sizeof(_n)/sizeof(double));
	std::vector<double> d_v(_d,_d+sizeof(_d)/sizeof(double));
	std::vector<double> t_v(_t,_t+sizeof(_t)/sizeof(double));
	std::vector<double> l_v(_l,_l+sizeof(_l)/sizeof(double));
	std::vector<double> eta_v(_eta,_eta+sizeof(_eta)/sizeof(double));
	std::vector<double> epsilon_v(_epsilon,_epsilon+sizeof(_epsilon)/sizeof(double));
	std::vector<double> beta_v(_beta,_beta+sizeof(_beta)/sizeof(double));
	std::vector<double> gamma_v(_gamma,_gamma+sizeof(_gamma)/sizeof(double));

	//Critical parameters
	crit.rho = 228; //[kg/m^3]
	crit.p = 3796; //[kPa]
	crit.T = 425.125; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 58.12220;
	params.Ttriple = 134.895;
	params.ptriple = 0.000665785834101;
	params.accentricfactor = 0.099;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n_v,d_v,t_v,l_v,1,23));
	phirlist.push_back(new phir_gaussian( n_v,d_v,t_v,eta_v,epsilon_v,beta_v,gamma_v,24,25));

	double _theta [] ={0,0,0,0,0.774840445,3.340602552,4.970513096,9.975553778};
	std::vector<double> theta_v (_theta,_theta+sizeof(_theta)/sizeof(double));
	double _n0 [] ={0, 12.54882924,-5.46976878,3.24680487,5.54913289,11.4648996,7.59987584,9.66033239};
	std::vector<double> n0_v (_n0,_n0+sizeof(_n0)/sizeof(double));
	
	// lead term: log(delta)+c+m*tau
	phi0list.push_back(new phi0_lead(n0_v[1], n0_v[2]));
	phi0list.push_back(new phi0_logtau(n0_v[3]));
	phi0list.push_back(new phi0_Planck_Einstein(n0_v,theta_v,4,7));

	EOSReference.assign("Buecker, D. and Wagner, W. \"Reference Equations of State for the Thermodynamic Properties of Fluid Phase n-Butane and Isobutane,\" J. Phys. Chem. Ref. Data, Vol. 35, No. 2, 2006, 929-1019.");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("n-Butane");
	aliases.push_back("nButane");
	REFPROPname.assign("BUTANE");
}
double nButaneClass::rhosatL(double T)
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	RHS = +1.97874515*pow(theta,0.345)
		  +0.85679951*pow(theta,1.0)
		  -0.341871887*pow(theta,1.5)
		  +0.304337558*pow(theta,3.0);
	rho = (1.0+ RHS)*rhoc;
	return rho;
}
double nButaneClass::rhosatV(double T)
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	RHS = -2.07770057*pow(theta,0.345)
		  -3.0836249*pow(theta,5.0/6.0)
		  -0.485645266*pow(theta,19.0/6.0)
		  -3.83167519*pow(theta,25.0/6.0);
	rho = exp(RHS*reduce.T/T)*rhoc;
	return rho;
}
double nButaneClass::psat(double T)
{
	double pc = reduce.p;
	double theta = 1-T/reduce.T;
	double RHS,p;

	RHS = -7.17616903*pow(theta,1.0)
		  +2.53635336*pow(theta,1.5)
		  -2.07532869*pow(theta,2.0)
		  -2.82241113*pow(theta,4.5);
	p = exp(reduce.T/T*RHS)*pc;
	return p;
}

IsoButaneClass::IsoButaneClass()
{
	double _n [] = {0, 2.0686820727966E+00, -3.6400098615204E+00, 5.1968754427244E-01, 1.7745845870123E-01, -1.2361807851599E-01, 4.5145314010528E-02, 3.0476479965980E-02, 7.5508387706302E-01, -8.5885381015629E-01, 3.6324009830684E-02, -1.9548799450550E-02, -4.4452392904960E-03, 4.6410763666460E-03, -7.1444097992825E-02, -8.0765060030713E-02, 1.5560460945053E-01, 2.0318752160332E-03, -1.0624883571689E-01, 3.9807690546305E-02, 1.6371431292386E-02, 5.3212200682628E-04, -7.8681561156387E-03, -3.0981191888963E-03, -4.2276036810382E-02, -5.3001044558079E-03};
	double _l[]= {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 0, 0};
	double _d[]= {0, 1, 1, 1, 2, 3, 4, 4, 1, 1, 2, 7, 8, 8, 1, 2, 3, 3, 4, 5, 5, 10, 2, 6, 1, 2};
	double _t[]= {0, 0.5, 1, 1.5, 0, 0.5, 0.5, 0.75, 2, 2.5, 2.5, 1.5, 1, 1.5, 4, 7, 3, 7, 3, 1, 6, 0, 6, 13, 2, 0};
	double _eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10};
	double _epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.85, 1};
	double _beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 150, 200};
	double _gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.16, 1.13};
	std::vector<double> n_v(_n,_n+sizeof(_n)/sizeof(double));
	std::vector<double> d_v(_d,_d+sizeof(_d)/sizeof(double));
	std::vector<double> t_v(_t,_t+sizeof(_t)/sizeof(double));
	std::vector<double> l_v(_l,_l+sizeof(_l)/sizeof(double));
	std::vector<double> eta_v(_eta,_eta+sizeof(_eta)/sizeof(double));
	std::vector<double> epsilon_v(_epsilon,_epsilon+sizeof(_epsilon)/sizeof(double));
	std::vector<double> beta_v(_beta,_beta+sizeof(_beta)/sizeof(double));
	std::vector<double> gamma_v(_gamma,_gamma+sizeof(_gamma)/sizeof(double));

	//Critical parameters
	crit.rho = 225.5; //[kg/m^3]
	crit.p = 3629; //[kPa]
	crit.T = 407.81; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 58.12220;
	params.Ttriple = 113.73;
	params.ptriple = 2.28968758984e-05;
	params.accentricfactor = 0.099;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n_v,d_v,t_v,l_v,1,23));
	phirlist.push_back(new phir_gaussian( n_v,d_v,t_v,eta_v,epsilon_v,beta_v,gamma_v,24,25));

	double _theta [] ={0,0,0,0,0.951277902,2.387895885,4.346904269,10.36885864};
	std::vector<double> theta_v (_theta,_theta+sizeof(_theta)/sizeof(double));
	double _n0 [] ={0,11.60865546,-5.29450411,3.05956619,4.94641014,4.09475197,15.6632824,9.73918122};
	std::vector<double> n0_v (_n0,_n0+sizeof(_n0)/sizeof(double));
	
	// lead term: log(delta)+c+m*tau
	phi0list.push_back(new phi0_lead(n0_v[1], n0_v[2]));
	phi0list.push_back(new phi0_logtau(n0_v[3]));
	phi0list.push_back(new phi0_Planck_Einstein(n0_v,theta_v,4,7));

	EOSReference.assign("Buecker, D. and Wagner, W. \"Reference Equations of State for the Thermodynamic Properties of Fluid Phase n-Butane and Isobutane,\" J. Phys. Chem. Ref. Data, Vol. 35, No. 2, 2006, 929-1019.");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("IsoButane");
	aliases.push_back("isobutane");
	REFPROPname.assign("ISOBUTAN");
}
double IsoButaneClass::rhosatL(double T)
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	RHS = +2.04025104*pow(theta,0.355)
		  +0.850874089*pow(theta,1.0)
		  -0.479052281*pow(theta,4.0/3.0)
		  +0.348201252*pow(theta,7.0/3.0);
	rho = (1.0 + RHS)*rhoc;
	return rho;
}
double IsoButaneClass::rhosatV(double T)
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	RHS = -2.12933323*pow(theta,0.355)
		  -2.93790085*pow(theta,5.0/6.0)
		  -0.89441086*pow(theta,19.0/6.0)
		  -3.46343707*pow(theta,26.0/6.0);
	rho = exp(RHS*reduce.T/T)*rhoc;
	return rho;
}
double IsoButaneClass::psat(double T)
{
	double pc = reduce.p;
	double theta = 1-T/reduce.T;
	double RHS,p;

	RHS = -6.85093103*pow(theta,1.0)
		  +1.36543198*pow(theta,1.5)
		  -1.32542691*pow(theta,2.5)
		  -2.56190994*pow(theta,4.5);
	p = exp(reduce.T/T*RHS)*pc;
	return p;
}
