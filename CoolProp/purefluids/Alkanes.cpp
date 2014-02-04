
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
	crit.rho = 10.139128*16.0428; //[kg/m^3]
	crit.p = PressureUnit(4599.2,UNIT_KPA); //[kPa]
	crit.T = 190.564; //[K]
	crit.v = 1/crit.rho; 

	// Reducing parameters used in EOS
	reduce.p = PressureUnit(4599.2, UNIT_KPA);
	reduce.T = 190.564; //[K]
	reduce.rho = 10.139128*16.0428; //[kg/m^3]
	reduce.v = 1.0/reduce.rho;

	preduce = &reduce;

	// Other fluid parameters
	params.molemass = 16.0428;
	params.Ttriple = 90.6941;
	params.ptriple = 11.696;
	params.accentricfactor = 0.01142;
	params.R_u = 8.31451;

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
	aliases.push_back("methane");
	aliases.push_back("METHANE");
	REFPROPname.assign("METHANE");

	BibTeXKeys.EOS = "Setzmann-JPCRD-1991";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
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
	double pc = reduce.p.Pa;
	double theta = 1-T/reduce.T;
	double RHS,p;

	RHS = -6.036219*pow(theta,1.0)
		  +1.409353*pow(theta,1.5)
		  -0.4945199*pow(theta,2.0)
		  -1.443048*pow(theta,4.5);
	p = exp(reduce.T/T*RHS)*pc;
	return p;
}
double MethaneClass::viscosity_dilute(double T)
{
	double C[] = {0, -3.0328138281, 16.918880086, -37.189364917, 41.288861858, -24.615921140, 8.9488430959, -1.8739245042, 0.20966101390, -9.6570437074e-3};
	double OMEGA_2_2 = 0, e_k, sigma, Tstar;

	// Get the L-J parameters
	this->ECSParams(&e_k,&sigma);

	Tstar = T/e_k;
	for (int i = 1; i <= 9; i++){
		OMEGA_2_2 += C[i]*pow(Tstar,(i-1)/3.0-1);
	}

	return 10.50*sqrt(Tstar)*OMEGA_2_2/1e6; //[Pa-s]
}
double MethaneClass::viscosity_residual(double T, double rho)
{
	double r[] = {0,1,1,2,2,2,3,3,4,4,1,1};
	double s[] = {0,0,1,0,1,1.5,0,2,0,1,0,1};
	double g[] = {0, 0.41250137, -0.14390912, 0.10366993, 0.40287464, -0.24903524, -0.12953131, 0.06575776, 0.02566628, -0.03716526, -0.38798341, 0.03533815};
	
	double sum1 = 0, sum2 = 0, tau = 190.551/T, delta = rho/(10.139*16.043);

	for (int i = 1; i<= 9; i++)
	{
		sum1 += g[i]*pow(delta,r[i])*pow(tau,s[i]);
	}
	for (int i = 10; i<= 11; i++)
	{
		sum2 += g[i]*pow(delta,r[i])*pow(tau,s[i]);
	}
	return 12.149*sum1/(1+sum2)/1e6;
}
double MethaneClass::viscosity_Trho(double T, double rho)
{
	return this->viscosity_dilute(T) + this->viscosity_residual(T,rho);
}
double MethaneClass::conductivity_dilute(double T)
{
	double e_k, sigma;
	// Get the L-J parameters
	this->ECSParams(&e_k,&sigma);

	double tau = 190.551/T, Tstar = T/e_k;
	double fint = 1.458850-0.4377162/Tstar;
	return 0.51826*(this->viscosity_dilute(T)*1e6)*(3.75-fint*(tau*tau*this->d2phi0_dTau2(tau,0)+1.5))/1e3; //[W/m/K]
}
double MethaneClass::conductivity_residual(double T, double rho)
{
	double delta_sigma_star;
	double r[] = {0,1,2,3,4,5,5,2};
	double s[] = {0,0,0,0,1,0,1,0};

	if (T < 190.551 && rho < 10.139*16.043){
		delta_sigma_star = rhosatV(T)/(10.139*16.043);
	}
	else{
		delta_sigma_star = 11;
	}

	double j[] = {0, 2.4149207, 0.55166331, -0.52837734, 0.073809553, 0.24465507, -0.047613626, 1.5554612};

	double sum1 = 0, tau = 190.551/T, delta = rho/(10.139*16.043);

	for (int i = 1; i <= 6; i++)
	{
		sum1 += j[i]*pow(delta,r[i])*pow(tau,s[i]);
	}
	sum1 += j[7]*delta*delta/delta_sigma_star;
	return 6.29638*sum1/1e3; //[W/m/K]
}
double MethaneClass::conductivity_Trho(double T, double rho)
{
	return this->conductivity_dilute(T) + this->conductivity_residual(T,rho) + this->conductivity_critical(T,rho,1/(0.545e-9),0.0563,0.19e-9);
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
	crit.rho = 6.856886685*30.06904;// 206.18; //[kg/m^3]
	crit.p = PressureUnit(4872.2,UNIT_KPA); //[kPa]
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
	aliases.push_back("ethane");
	aliases.push_back("ETHANE");
	REFPROPname.assign("ETHANE");

	BibTeXKeys.EOS = "Buecker-JPCRD-2006";
	BibTeXKeys.VISCOSITY = "Friend-JPCRD-1991";
	BibTeXKeys.CONDUCTIVITY = "Friend-JPCRD-1991";
	BibTeXKeys.ECS_LENNARD_JONES = "Friend-JPCRD-1991";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
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
	// Maximum absolute error is 0.270124 % between 178.000001 K and 591.749999 K
    const double ti[]={0,0.374, 1.0, 1.5, 2.6666666666666665, 3.8333333333333335, 0.38849999999999996};
    const double Ni[]={0,-1.4796534339188399, -3.1822071897575968, 1.0445932652363741, 0.41065304624296084, -3.5816879898373259, -0.96561226454570503};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(reduce.T/T*summer);
}
double EthaneClass::psat(double T)
{
	double pc = reduce.p.Pa;
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


double EthaneClass::viscosity_dilute(double T)
{
	double C[] = {0, -3.0328138281, 16.918880086, -37.189364917, 41.288861858, -24.615921140, 8.9488430959, -1.8739245042, 0.20966101390, -9.6570437074e-3};
	double OMEGA_2_2 = 0, e_k, sigma, Tstar;

	// Get the L-J parameters
	this->ECSParams(&e_k,&sigma);

	Tstar = T/e_k;
	for (int i = 1; i<= 9; i++)
	{
		OMEGA_2_2 += C[i]*pow(Tstar,(i-1)/3.0-1);
	}

	return 12.0085*sqrt(Tstar)*OMEGA_2_2/1e6; //[Pa-s]
}
double EthaneClass::viscosity_residual(double T, double rho)
{
	double r[] = {0,1,1,2,2,2,3,3,4,4,1,1};
	double s[] = {0,0,1,0,1,1.5,0,2,0,1,0,1};
	double g[] = {0, 0.47177003, -0.23950311, 0.39808301, -0.27343335, 0.35192260, -0.21101308, -0.00478579, 0.07378129, -0.030435255, -0.30435286, 0.001215675};

	double sum1 = 0, sum2 = 0, tau = 305.33/T, delta = rho/(6.87*30.070);

	for (int i = 1; i<= 9; i++)
	{
		sum1 += g[i]*pow(delta,r[i])*pow(tau,s[i]);
	}
	for (int i = 10; i<= 11; i++)
	{
		sum2 += g[i]*pow(delta,r[i])*pow(tau,s[i]);
	}
	return 15.977*sum1/(1+sum2)/1e6;
}
double EthaneClass::viscosity_Trho(double T, double rho)
{
	return this->viscosity_dilute(T) + this->viscosity_residual(T,rho);
}
double EthaneClass::conductivity_dilute(double T)
{
	double e_k, sigma;
	// Get the L-J parameters
	this->ECSParams(&e_k,&sigma);

	double tau = 305.33/T, Tstar = T/e_k;
	double fint = 1.7104147-0.6936482/Tstar;
	return 0.276505*(this->viscosity_dilute(T)*1e6)*(3.75-fint*(tau*tau*this->d2phi0_dTau2(tau,0)+1.5))/1e3; //[W/m/K]
}
double EthaneClass::conductivity_residual(double T, double rho)
{
	double r[] = {0,1,2,3,4,5,1,3};
	double s[] = {0,0,0,0,0,0,1.5,1};
	double j[] = {0,0.96084322,2.7500235,-0.026609289,-0.078146729,0.21881339,2.3849563,-0.75113971};

	double sum1 = 0, tau = 305.33/T, delta = rho/(6.87*30.070);

	for (int i = 1; i<= 7; i++)
	{
		sum1 += j[i]*pow(delta,r[i])*pow(tau,s[i]);
	}
	return 4.41786*sum1/1e3; //[W/m/K]
}
double EthaneClass::conductivity_Trho(double T, double rho)
{
	return this->conductivity_dilute(T) + this->conductivity_residual(T,rho) + this->conductivity_critical(T,rho,1/(0.545e-9),0.0563,0.19e-9);
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
	crit.p = PressureUnit(3796,UNIT_KPA); //[kPa]
	crit.T = 425.125; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 58.12220;
	params.Ttriple = 134.895;
	params.ptriple = 0.000665785834101;
	params.accentricfactor = 0.200810094644;
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
	aliases.push_back("butane");
	aliases.push_back("BUTANE");
	aliases.push_back("N-BUTANE");
	REFPROPname.assign("BUTANE");

	BibTeXKeys.EOS = "Buecker-JPCRD-2006B";
	BibTeXKeys.VISCOSITY = "Vogel-HTHP-1999";
	BibTeXKeys.CONDUCTIVITY = "Perkins-JCED-2002A";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Vogel-HTHP-1999";
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
	// Max error is  0.121857408131 % between 134.895 and 425.124999 K
    const double ti[]={0, 0.355, 0.38799999999999996, 0.8333333333333334, 17.166666666666668, 3.5, 4.666666666666667};
    const double Ni[]={0, -1.0295542341081576, -1.2940836117556569, -2.7620656950014832, 1.8460149445757177, -2.4159729873675766, -1.9343557659300938};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(reduce.T/T*summer);
}
double nButaneClass::psat(double T)
{
	double pc = reduce.p.Pa;
	double theta = 1-T/reduce.T;
	double RHS,p;

	RHS = -7.17616903*pow(theta,1.0)
		  +2.53635336*pow(theta,1.5)
		  -2.07532869*pow(theta,2.0)
		  -2.82241113*pow(theta,4.5);
	p = exp(reduce.T/T*RHS)*pc;
	return p;
}
double nButaneClass::viscosity_Trho(double T, double rho)
{
	double a[] = {0.17067154, -0.48879666, 0.039038856};
	double b[] = {-19.572881, 219.73999, -1015.3226, 2471.01251, -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158};

	double e_k, sigma;
	this->ECSParams(&e_k,&sigma);
	double Tstar = T/e_k;
	double Gstar = exp(a[0]+a[1]*log(Tstar)+a[2]*log(Tstar)*log(Tstar));
	double eta_0 = 0.021357*sqrt(params.molemass*T)/(sigma*sigma*Gstar); // uPa-s

	//Rainwater-Friend initial density term
	double Bstar = b[0]*pow(Tstar,-0.25*0)+b[1]*pow(Tstar,-0.25*1)+b[2]*pow(Tstar,-0.25*2)+b[3]*pow(Tstar,-0.25*3)+b[4]*pow(Tstar,-0.25*4)+b[5]*pow(Tstar,-0.25*5)+b[6]*pow(Tstar,-0.25*6)+b[7]*pow(Tstar,-2.5)+b[8]*pow(Tstar,-5.5);
	double B = Bstar*0.6022137*sigma*sigma*sigma; // [L/mol]

	double e[6][2]; // init with zeros
	e[2][0] = -54.7737770846; e[2][1] = 58.0898623034;
	e[3][0] = 35.2658446259; e[3][1] = -39.6682203832;
	e[4][0] = -1.83729542151; e[4][1] = 0;
	e[5][0] = -0.833262985358; e[5][1] = 1.93837020663;
	double f1 = 188.075903903;
	double g1 = 2.30873963359, g2 = 0.881017652640;

	double sumresid = 0;
	double tau = T/crit.T, delta = rho/(3.920*58.1222);
	for (int i = 2; i<=5; i++)
	{
		for (int j = 0; j< 2; j++)
		{
			sumresid += e[i][j]*pow(delta,i)/pow(tau,j);
		}
	}

	double delta_0 = g1*(1+g2*sqrt(tau));
	double eta_r = sumresid + f1*(delta/(delta_0-delta)-delta/delta_0); // uPa-s
	
	double rhobar = rho/params.molemass; //mol/L
	return (eta_0*(1+B*rhobar)+eta_r)/1e6;
}
double nButaneClass::conductivity_Trho(double T, double rho)
{
	double lambda_0 = 1.62676e-3+9.75703e-4*(T/crit.T) + 2.89887e-2*pow(T/crit.T,2); // W/m/K

	double sumresid = 0;
	double B1[] = {0, -3.08823e-2, 1.59698e-1, -1.41629e-1, 5.03252e-2, -6.04344e-3};
	double B2[] = {0, 4.22711e-2, -1.43867e-1, 1.30043e-1, -4.73921e-2, 6.31824e-3};

	for (int i = 1; i<= 5; i++){		
		sumresid += (B1[i]+B2[i]*(T/reduce.T))*pow(rho/reduce.rho,i);
	}

	double lambda_r = sumresid; // [W/m/K]

	double lambda_c = this->conductivity_critical(T,rho,1.0/(6.12930e-10)); // [W/m/K]

	return lambda_0+lambda_r+lambda_c;
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
	crit.p = PressureUnit(3629,UNIT_KPA); //[kPa]
	crit.T = 407.817; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 58.12220;
	params.Ttriple = 113.73;
	params.ptriple = 2.28968758984e-05;
	params.accentricfactor = 0.183531783208;
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
	aliases.push_back("Isobutane");
	aliases.push_back("ISOBUTANE");
	aliases.push_back("R600A");
	aliases.push_back("R600a");
	REFPROPname.assign("ISOBUTAN");

	// Adjust to the IIR reference state (h=200 kJ/kg, s = 1 kJ/kg for sat. liq at 0C)
    params.HSReferenceState = "IIR";

	BibTeXKeys.EOS = "Buecker-JPCRD-2006B";
	BibTeXKeys.VISCOSITY = "Vogel-IJT-2000";
	BibTeXKeys.ECS_LENNARD_JONES = "Vogel-IJT-2000";
	BibTeXKeys.CONDUCTIVITY = "Perkins-JCED-2002B";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
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
	// Maximum absolute error is 0.270124 % between 178.000001 K and 591.749999 K
    const double ti[]={0,0.38949999999999996, 0.39899999999999997, 0.8333333333333334, 1.5, 2.6666666666666665, 3.5};
    const double Ni[]={0,-20.080357830155275, 18.457306494924232, -3.7935794935101383, 0.27898485736137468, 1.7770827488706802, -5.6964841742935048};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(reduce.T/T*summer);
}
double IsoButaneClass::psat(double T)
{
	// Max error is  0.0153125011953 % between 113.73 and 407.816999 K

    const double ti[]={0, 16.5, 1.0, 1.3333333333333333, 2.8333333333333335, 4.833333333333333, 7.333333333333333};
    const double Ni[]={0, -1.5783244561468581, -7.0007865713821147, 1.1998975400416103, -1.2029491215695771, -2.9725455394203251, 1.1271706682355129};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
void IsoButaneClass::ECSParams(double *e_k, double *sigma)
{
	// Vogel, 2000
	*e_k = 307.55;
	*sigma = 0.46445;
}
double IsoButaneClass::viscosity_Trho(double T, double rho)
{
	double a[] = {0.53583008, -0.45629630, 0.049911282};
	double b[] = {-19.572881, 219.73999, -1015.3226, 2471.01251, -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158};

	double e_k, sigma;
	this->ECSParams(&e_k,&sigma);
	double Tstar = T/e_k;
	double Gstar = exp(a[0]+a[1]*log(Tstar)+a[2]*log(Tstar)*log(Tstar));
	double eta_0 = 0.021357*sqrt(params.molemass*T)/(sigma*sigma*Gstar); // uPa-s

	//Rainwater-Friend
	double Bstar = b[0]*pow(Tstar,-0.25*0)+b[1]*pow(Tstar,-0.25*1)+b[2]*pow(Tstar,-0.25*2)+b[3]*pow(Tstar,-0.25*3)+b[4]*pow(Tstar,-0.25*4)+b[5]*pow(Tstar,-0.25*5)+b[6]*pow(Tstar,-0.25*6)+b[7]*pow(Tstar,-2.5)+b[8]*pow(Tstar,-5.5);
	double B = Bstar*0.602214129*sigma*sigma*sigma; // [L/mol]

	double e[6][3]; // init with zeros
	e[2][0] = 103.511763411; e[2][1] = -312.670896234;
	e[2][2] = 145.253750239; e[3][0] = -210.649894193;
	e[3][1] = 386.269696509; e[3][2] = -214.963015527;
	e[4][0] = 112.580360920; e[4][1] = -223.242033154;
	e[4][2] = 119.114788598; e[5][0] = -18.1909745900;
	e[5][1] = 36.0438957232; e[5][2] = -21.3960184050;
	double f1 = 1940.37606990;
	double g1 = 2.33859774637, g2 = 1.00596672174;

	double sumresid = 0;
	double tau = T/crit.T, delta = rho/(3.860*58.1222);
	for (int i = 2; i < 6; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			sumresid += e[i][j]*pow(delta,i)/pow(tau,j);
		}
	}

	double delta_0 = g1*(1+g2*sqrt(tau));
	double eta_r = sumresid + f1*(delta/(delta_0-delta)-delta/delta_0); // uPa-s
	
	double rhobar = rho/params.molemass;
	return (eta_0*(1+B*rhobar)+eta_r)/1e6;
}
double IsoButaneClass::conductivity_Trho(double T, double rho)
{
	double lambda_0 = -2.37901e-3+1.06601e-2*(T/crit.T) + 2.15811e-2*pow(T/crit.T,2); // W/m/K

	double sumresid = 0;
	double B1[] = {0, -3.94953e-2, 1.61607e-1, -1.38049e-1, 4.83126e-2, -5.78452e-3};
	double B2[] = {0, 4.51967e-2, -1.34395e-1, 1.15446e-1, -4.11718e-2, 5.43111e-3};

	for (int i = 1; i<= 5; i++){
		sumresid += (B1[i]+B2[i]*(T/reduce.T))*pow(rho/reduce.rho,i);
	}

	double lambda_r = sumresid; // [W/m/K]

	double lambda_c = this->conductivity_critical(T,rho,1.0/(5.37809e-10)); // [W/m/K]

	return lambda_0+lambda_r+lambda_c;
}
