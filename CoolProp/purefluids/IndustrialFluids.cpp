

// **** WARNING ******
// **** WARNING ******
// **** WARNING ******

// Do NOT modify this file.  It is created by a script in the industrialfluidsbuilder folder within the source

#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "IndustrialFluids.h"
#include "R134a.h"

static const double d_nonpolar[] =
{
0,
1.0, //[1]
1.0, //[2]
1.0, //[3]
2.0, //[4]
3.0, //[5]
7.0, //[6]
2.0, //[7]
5.0, //[8]
1.0, //[9]
4.0, //[10]
3.0, //[11]
4.0, //[12]
};

static const double t_nonpolar[] =
{
0,
0.25,  //[1]
1.125, //[2]
1.5,   //[3]
1.375, //[4]
0.25,  //[5]
0.875, //[6]
0.625, //[7]
1.75,  //[8]
3.625, //[9]
3.625, //[10]
14.5,  //[11]
12.0,  //[12]
};

static const double c_nonpolar[] =
{
0,
0.0, //[1]
0.0, //[2]
0.0, //[3]
0.0, //[4]
0.0, //[5]
0.0, //[6]
1.0, //[7]
1.0, //[8]
2.0, //[9]
2.0, //[10]
3.0, //[11]
3.0, //[12]
};

static const double d_polar[] =
{
0,
1.0, //[1]
1.0, //[2]
1.0, //[3]
3.0, //[4]
7.0, //[5]
1.0, //[6]
2.0, //[7]
5.0, //[8]
1.0, //[9]
1.0, //[10]
4.0, //[11]
2.0, //[12]
};

static const double t_polar[] =
{
0,
0.25,  //[1]
1.25,  //[2]
1.5,   //[3]
0.25,  //[4]
0.875, //[5]
2.375, //[6]
2.0,   //[7]
2.125, //[8]
3.5,   //[9]
6.5,   //[10]
4.75,  //[11]
12.5,  //[12]
};

static const double c_polar[] =
{
0,
0.0, //[1]
0.0, //[2]
0.0, //[3]
0.0, //[4]
0.0, //[5]
1.0, //[6]
1.0, //[7]
1.0, //[8]
2.0, //[9]
2.0, //[10]
2.0, //[11]
3.0, //[12]
};

CarbonMonoxideClass::CarbonMonoxideClass()
{
    const double n[]={0.0,0.905540000,-2.451500000,0.531490000,0.024173000,0.072156000,0.000188180,0.194050000,-0.043268000,-0.127780000,-0.027896000,-0.034154000,0.016329000};
    const double u0[]={0.0,3089.0};
    const double v0[]={0.0,1.0128};

    // Critical parameters
    crit.rho = 303.909585;
    crit.p = PressureUnit(3494.0, UNIT_KPA);
    crit.T = 132.86;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 28.0101;
    params.Ttriple = 68.16;
	params.ptriple = 15.5395203075;
    params.accentricfactor = 0.0497;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 68.16;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-3.3728318564,3.3683460039);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(3.5-1);
    phi0list.push_back(phi0_logtau_);

    phi_BC * phi0_power_ = new phi0_power(-2.2311e-07*pow(crit.T,1.5)/(1.5*(1.5+1)),-1.5);
	phi0list.push_back(phi0_power_);
    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("CarbonMonoxide");
    aliases.push_back(std::string("CO"));
    aliases.push_back(std::string("CARBONMONOXIDE"));
    REFPROPname.assign("CO");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";
}

CarbonylSulfideClass::CarbonylSulfideClass()
{
    const double n[]={0.0,0.943740000,-2.534800000,0.590580000,-0.021488000,0.082083000,0.000246890,0.212260000,-0.041251000,-0.223330000,-0.050828000,-0.028333000,0.016983000};
    const double u0[]={0.0,768.0,1363.0,3175.0,12829.0};
    const double v0[]={0.0,2.1651,0.93456,1.0623,0.34269};

    // Critical parameters
    crit.rho = 445.156491;
    crit.p = PressureUnit(6370.0, UNIT_KPA);
    crit.T = 378.77;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 60.0751;
    params.Ttriple = 134.3;
	params.ptriple = 0.0644440370601;
    params.accentricfactor = 0.0978;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 134.3;
    limits.Tmax = 650.0;
    limits.pmax = 50000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-3.6587449805,3.7349245016);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(3.5-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("CarbonylSulfide");
    aliases.push_back(std::string("COS"));
    aliases.push_back(std::string("CARBONYLSULFIDE"));
    REFPROPname.assign("COS");

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";
}

DecaneClass::DecaneClass()
{
    const double n[]={0.0,1.046100000,-2.480700000,0.743720000,-0.525790000,0.153150000,0.000328650,0.841780000,0.055424000,-0.735550000,-0.185070000,-0.020775000,0.012335000};
    const double u0[]={0.0,1193.0,2140.0,4763.0,10862.0};
    const double v0[]={0.0,25.685,28.233,12.417,10.035};

    // Critical parameters
    crit.rho = 233.3419552;
    crit.p = PressureUnit(2103.0, UNIT_KPA);
    crit.T = 617.7;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 142.28168;
    params.Ttriple = 243.5;
	params.ptriple = 0.00140434258288;
    params.accentricfactor = 0.4884;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 243.5;
    limits.Tmax = 675.0;
    limits.pmax = 800000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(13.9361966549,-10.5265128286);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(19.109-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("n-Decane");
    aliases.push_back("Decane"); 
	aliases.push_back("decane"); 
	aliases.push_back(std::string("DECANE"));
	aliases.push_back(std::string("N-DECANE"));
    REFPROPname.assign("decane");

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-FPE-2004";
	BibTeXKeys.VISCOSITY = "Huber-FPE-2004";
	BibTeXKeys.CONDUCTIVITY = "Huber-FPE-2005";
}
void DecaneClass::ECSParams(double *e_k, double *sigma)
{
	// From Huber 2004
	*e_k = 490.510;
	*sigma = 0.68600;
}
double DecaneClass::viscosity_Trho(double T, double rho)
{
	// Rainwater-Friend for initial density dependence
	double e_k, sigma;
	this->ECSParams(&e_k,&sigma);
	double Tstar = T/e_k;

	// Dilute gas
	double eta_0 = 0.021357*sqrt(params.molemass*T)/(sigma*sigma*exp(0.343267-0.460514*log(Tstar))); // uPa-s

	// Initial density dependence from Rainwater-Friend
	double b[] = {-19.572881, 219.73999, -1015.3226, 2471.01251, -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158};
	double Bstar = b[0]*pow(Tstar,-0.25*0)+b[1]*pow(Tstar,-0.25*1)+b[2]*pow(Tstar,-0.25*2)+b[3]*pow(Tstar,-0.25*3)+b[4]*pow(Tstar,-0.25*4)+b[5]*pow(Tstar,-0.25*5)+b[6]*pow(Tstar,-0.25*6)+b[7]*pow(Tstar,-2.5)+b[8]*pow(Tstar,-5.5);
	double B = Bstar*0.60221415*sigma*sigma*sigma; // L/mol

	double e[4][3]; // init with zeros
	e[2][1] = -0.402094e-1;
	e[2][2] =  0.404435e-1;	
	e[3][1] =  0;
	e[3][2] = -0.142063e-1;

	double c[] = {0, 0.453387, 2.55105, 1.71465, 0};

	double sumresid = 0;
	double tau = T/crit.T, delta = rho/crit.rho;
	for (int j = 2; j <= 3; j++)
	{
		for (int k = 1; k <= 2; k++)
		{
			sumresid += e[j][k]*pow(delta,j)/pow(tau,k);
		}
	}
	double delta_0 = c[2] + c[3]*sqrt(tau) + c[4]*tau;
	double eta_r = (sumresid + c[1]*(delta/(delta_0-delta)-delta/delta_0))*1000; // uPa-s

	double rhobar = rho/params.molemass; // [mol/L]
	return (eta_0*(1+B*rhobar) + eta_r)/1e6;
}
double DecaneClass::conductivity_Trho(double T, double rho)
{
	double lambda_0 = 1.05542680e-2 - 5.14530090e-2*(T/crit.T) + 1.18978971e-1*pow(T/crit.T,2) - 3.72442104e-2*pow(T/crit.T,3); // W/m/K

	double sumresid = 0;
	double B1[] = {0, -2.94394112e-2, 4.99245356e-2, -1.42700394e-2, 1.50827597e-3};
	double B2[] = {0, 1.50509474e-2, 0, -1.38857133e-2, 4.33326339e-3};

	for (int i = 1; i<= 4; i++){
		sumresid += (B1[i]+B2[i]*(T/reduce.T))*pow(rho/reduce.rho,i);
	}

	double lambda_r = sumresid; // [W/m/K]

	double lambda_c = this->conductivity_critical(T,rho,1.41115586e9); // [W/m/K]

	return lambda_0+lambda_r+lambda_c;
}

HydrogenSulfideClass::HydrogenSulfideClass()
{
    const double n[]={0.0,0.876410000,-2.036700000,0.216340000,-0.050199000,0.066994000,0.000190760,0.202270000,-0.004534800,-0.222300000,-0.034714000,-0.014885000,0.007415400};
    const double u0[]={0.0,1823.0,3965.0};
    const double v0[]={0.0,1.1364,1.9721};

    // Critical parameters
    crit.rho = 347.2841672;
    crit.p = PressureUnit(9000.0, UNIT_KPA);
    crit.T = 373.1;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 34.08088;
    params.Ttriple = 187.7;
	params.ptriple = 23.2604252601;
    params.accentricfactor = 0.1005;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 187.7;
    limits.Tmax = 760.0;
    limits.pmax = 170000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-4.0740770957,3.7632137341);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);

    phi_BC * phi0_power_ = new phi0_power(-1.4327e-06*pow(crit.T,1.5)/(1.5*(1.5+1)),-1.5);
	phi0list.push_back(phi0_power_);
    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("HydrogenSulfide");
    aliases.push_back(std::string("H2S")); 
    aliases.push_back(std::string("HYDROGENSULFIDE"));
    REFPROPname.assign("H2S");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "QuinonesCisneros-JCED-2012";
	BibTeXKeys.VISCOSITY = "QuinonesCisneros-JCED-2012";
}
void HydrogenSulfideClass::ECSParams(double *e_k, double *sigma)
{
	*e_k = 355.8;
	*sigma = 0.3565;
}
double HydrogenSulfideClass::viscosity_Trho(double T, double rho)
{	

	// Dilute
	double a[] = {0.53242, 0.93715, -0.69339, 1.16432, -0.84306, 0.20534};
	double Sstar = 0;
	for (int i = 0; i <= 5; i++) {Sstar += a[i]/pow(T/276.0,i); }
	double eta_0 = 0.87721*sqrt(T)/Sstar; // uPa-s

	// Rainwater-Friend for initial density dependence
	double e_k, sigma;
	this->ECSParams(&e_k,&sigma);
	double Tstar = T/e_k;
	double rhobar = rho/params.molemass; // mol/L

	double b[] = {-19.572881, 219.73999, -1015.3226, 2471.01251, -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158};
	double Bstar = b[0]*pow(Tstar,-0.25*0)+b[1]*pow(Tstar,-0.25*1)+b[2]*pow(Tstar,-0.25*2)+b[3]*pow(Tstar,-0.25*3)+b[4]*pow(Tstar,-0.25*4)+b[5]*pow(Tstar,-0.25*5)+b[6]*pow(Tstar,-0.25*6)+b[7]*pow(Tstar,-2.5)+b[8]*pow(Tstar,-5.5);
	double B = Bstar*0.60221415*sigma*sigma*sigma; // L/mol
	
	double eta_i = eta_0*B*rhobar; // uPa-s

	// Residual part
	double psi1 = exp(crit.T/T);
	double psi2 = exp(crit.T*crit.T/T/T);
	double a0 = 68.9659e-6, b0 = 153.406e-6, A0 = 0.782380e-9, B0 = -9.75792e-9;
	double a1 = -22.0494e-6, b1 = 8.45198e-6, A1 = -0.64717e-9, B1 = -3.19303e-9;
	double a2 = -42.6126e-6, b2 = -113.967e-6, A2 = 1.39066e-9, B2 = 12.4263e-9;
	double ka = (a0 + a1*psi1 + a2*psi2)*crit.T/T;
	double kr = (b0 + b1*psi1 + b2*psi2)*crit.T/T;
	double kaa = (A0 + A1*psi1 + A2*psi2)*crit.T/T;
	double krr = (B0 + B1*psi1 + B2*psi2)*crit.T/T;

	double p = this->pressure_Trho(T,rho)/100; // kPa -> bar
	double pr = T*this->dpdT_Trho(T,rho)/100; // kPa-> bar
	double pa = p - pr;
	double pid = rho * R() * T / 100; // kPa -> bar
	double deltapr = pr - pid;

	double eta_f = (ka*pa + kr*deltapr + kaa*pa*pa + krr*pr*pr)*1000; //mPa-s --> uPa-s

	return (eta_0 + eta_i + eta_f)/1e6;
}

IsopentaneClass::IsopentaneClass()
{
    const double n[]={0.0,1.096300000,-3.040200000,1.031700000,-0.154100000,0.115350000,0.000298090,0.395710000,-0.045881000,-0.358040000,-0.101070000,-0.035484000,0.018156000};
    const double u0[]={0.0,442.0,1109.0,2069.0,4193.0};
    const double v0[]={0.0,7.4056,9.5772,15.765,12.119};

    // Critical parameters
    crit.rho = 235.99865938;
    crit.p = PressureUnit(3378.0, UNIT_KPA);
    crit.T = 460.35;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 72.14878;
    params.Ttriple = 112.65;
	params.ptriple = 8.93808446917e-08;
    params.accentricfactor = 0.2274;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 112.65;
    limits.Tmax = 500.0;
    limits.pmax = 1000000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(2.5822330405,1.1609103419);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Isopentane");
    aliases.push_back(std::string("ipentane"));
    aliases.push_back(std::string("R601a"));
    aliases.push_back(std::string("ISOPENTANE"));
    REFPROPname.assign("ipentane");

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Chichester-NIST-2008";
}

NeopentaneClass::NeopentaneClass()
{
    const double n[]={0.0,1.113600000,-3.179200000,1.141100000,-0.104670000,0.117540000,0.000340580,0.295530000,-0.074765000,-0.314740000,-0.099401000,-0.039569000,0.023177000};
    const double u0[]={0.0,710.0,1725.0,3280.0,7787.0};
    const double v0[]={0.0,14.422,12.868,17.247,12.663};

    // Critical parameters
    crit.rho = 235.9265106;
    crit.p = PressureUnit(3196.0, UNIT_KPA);
    crit.T = 433.74;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 72.14878;
    params.Ttriple = 256.6;
	params.ptriple = 35.400947327081248;
    params.accentricfactor = 0.1961;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 256.6;
    limits.Tmax = 550.0;
    limits.pmax = 200000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(0.8702452614,1.6071746358);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Neopentane");
    aliases.push_back(std::string("neopentn"));
    aliases.push_back(std::string("NEOPENTANE"));
    REFPROPname.assign("neopentn");

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.ECS_LENNARD_JONES = "Chichester-NIST-2008";
}

IsohexaneClass::IsohexaneClass()
{
    const double n[]={0.0,1.102700000,-2.969900000,1.029500000,-0.212380000,0.118970000,0.000277380,0.401030000,-0.034238000,-0.435840000,-0.116930000,-0.019262000,0.008078300};
    const double u0[]={0.0,325.0,1150.0,2397.0,5893.0};
    const double v0[]={0.0,7.9127,16.871,19.257,14.075};

    // Critical parameters
    crit.rho = 233.9661024;
    crit.p = PressureUnit(3040.0, UNIT_KPA);
    crit.T = 497.7;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 86.17536;
    params.Ttriple = 119.6;
	params.ptriple = 7.67397444618e-09;
    params.accentricfactor = 0.2797;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 119.6;
    limits.Tmax = 550.0;
    limits.pmax = 1000000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(6.9259123919,-0.3128629679);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);

    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Isohexane");
    aliases.push_back(std::string("ihexane"));
    aliases.push_back(std::string("ISOHEXANE"));
    REFPROPname.assign("ihexane");

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Chichester-NIST-2008";
}

KryptonClass::KryptonClass()
{
    const double n[]={0.0,0.835610000,-2.372500000,0.545670000,0.014361000,0.066502000,0.000193100,0.168180000,-0.033133000,-0.150080000,-0.022897000,-0.021454000,0.006939700};
    const double u0[]={0.0,0};
    const double v0[]={0.0,0};

    // Critical parameters
    crit.rho = 909.2083;
    crit.p = PressureUnit(5525.0, UNIT_KPA);
    crit.T = 209.48;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 83.798;
    params.Ttriple = 115.77;
	params.ptriple = 73.5090071088;
    params.accentricfactor = -0.00089;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 115.77;
    limits.Tmax = 750.0;
    limits.pmax = 200000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-3.7506412806,3.7798018435);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(2.5-1);
    phi0list.push_back(phi0_logtau_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Krypton");
    aliases.push_back(std::string("krypton")); 
    aliases.push_back(std::string("KRYPTON"));
    REFPROPname.assign("krypton");

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";
}

NonaneClass::NonaneClass()
{
    const double n[]={0.0,1.115100000,-2.702000000,0.834160000,-0.388280000,0.137600000,0.000281850,0.620370000,0.015847000,-0.617260000,-0.150430000,-0.012982000,0.004432500};
    const double u0[]={0.0,1221.0,2244.0,5008.0,11724.0};
    const double v0[]={0.0,24.926,24.842,11.188,17.483};

    // Critical parameters
    crit.rho = 232.141731;
    crit.p = PressureUnit(2281.0, UNIT_KPA);
    crit.T = 594.55;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 128.2551;
    params.Ttriple = 219.7;
	params.ptriple = 0.000444543592359;
    params.accentricfactor = 0.4433;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 219.7;
    limits.Tmax = 600.0;
    limits.pmax = 800000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(10.7927224829,-8.2418318753);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(17.349-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("n-Nonane");
    aliases.push_back(std::string("nonane")); 
    aliases.push_back(std::string("NONANE"));
    aliases.push_back(std::string("N-NONANE"));
    REFPROPname.assign("nonane");

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-FPE-2004";
	BibTeXKeys.VISCOSITY = "Huber-FPE-2004";
	BibTeXKeys.CONDUCTIVITY = "Huber-FPE-2005";
}
void NonaneClass::ECSParams(double *e_k, double *sigma)
{
	// From Huber 2004
	*e_k = 472.127;
	*sigma = 0.66383;
}
double NonaneClass::viscosity_Trho(double T, double rho)
{
	// Rainwater-Friend for initial density dependence
	double e_k, sigma;
	this->ECSParams(&e_k,&sigma);
	double Tstar = T/e_k;

	// Dilute gas
	double eta_0 = 0.021357*sqrt(params.molemass*T)/(sigma*sigma*exp(0.340344-0.466455*log(Tstar))); // uPa-s

	// Initial density dependence from Rainwater-Friend
	double b[] = {-19.572881, 219.73999, -1015.3226, 2471.01251, -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158};
	double Bstar = b[0]*pow(Tstar,-0.25*0)+b[1]*pow(Tstar,-0.25*1)+b[2]*pow(Tstar,-0.25*2)+b[3]*pow(Tstar,-0.25*3)+b[4]*pow(Tstar,-0.25*4)+b[5]*pow(Tstar,-0.25*5)+b[6]*pow(Tstar,-0.25*6)+b[7]*pow(Tstar,-2.5)+b[8]*pow(Tstar,-5.5);
	double B = Bstar*0.60221415*sigma*sigma*sigma; // L/mol

	double e[4][3]; // init with zeros
	e[2][1] = -0.314367e-1;		
	e[2][2] =  0.326258e-1;	
	e[3][1] =  0.639384e-2;
	e[3][2] = -0.108922e-1;

	double c[] = {0, 0.192935, 2.66987, 1.32137, 0};

	double sumresid = 0;
	double tau = T/crit.T, delta = rho/crit.rho;
	for (int j = 2; j <= 3; j++)
	{
		for (int k = 1; k <= 2; k++)
		{
			sumresid += e[j][k]*pow(delta,j)/pow(tau,k);
		}
	}
	double delta_0 = c[2] + c[3]*sqrt(tau) + c[4]*tau;
	double eta_r = (sumresid + c[1]*(delta/(delta_0-delta)-delta/delta_0))*1000; // uPa-s

	double rhobar = rho/params.molemass; // [mol/L]
	return (eta_0*(1+B*rhobar) + eta_r)/1e6;
}
double NonaneClass::conductivity_Trho(double T, double rho)
{
	double lambda_0 = 8.7877e-3 - 4.1351e-2*(T/crit.T) + 1.0479e-1*pow(T/crit.T,2) - 3.2003e-2*pow(T/crit.T,3); // W/m/K

	double sumresid = 0;
	double B1[] = {0, 4.90087596e-3, -8.07305471e-3, 5.57430614e-3, 0};
	double B2[] = {0, 9.96486280e-3, 0 , 0, 0};

	for (int i = 1; i<= 4; i++){
		sumresid += (B1[i]+B2[i]*(T/reduce.T))*pow(rho/reduce.rho,i);
	}

	double lambda_r = sumresid; // [W/m/K]

	double lambda_c = this->conductivity_critical(T,rho,9.58722814e8); // [W/m/K]

	return lambda_0+lambda_r+lambda_c;
}

TolueneClass::TolueneClass()
{
    const double n[]={0.0,0.964640000,-2.785500000,0.867120000,-0.188600000,0.118040000,0.000251810,0.571960000,-0.029287000,-0.433510000,-0.125400000,-0.028207000,0.014076000};
    const double u0[]={0.0,190.0,797.0,1619.0,3072.0,7915.0};
    const double v0[]={0.0,1.6994,8.0577,17.059,8.4567,8.6423};

    // Critical parameters
    crit.rho = 291.98665298;
    crit.p = PressureUnit(4126.0, UNIT_KPA);
    crit.T = 591.75;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 92.13842;
    params.Ttriple = 178.0;
	params.ptriple = 3.94003300153e-05;
    params.accentricfactor = 0.2657;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 178.0;
    limits.Tmax = 700.0;
    limits.pmax = 500000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(3.5241174832,1.1360823464);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);

    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Toluene");
    aliases.push_back(std::string("toluene")); 
    aliases.push_back(std::string("TOLUENE"));
    REFPROPname.assign("toluene");

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.VISCOSITY = "";
	BibTeXKeys.CONDUCTIVITY = "ASSAEL-JPCRD-2012B";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double TolueneClass::conductivity_Trho(double T, double rho)
{
	double sumresid = 0;
	double B1[] = {0, -5.18530e-2, 1.33846e-1, -1.20446e-1, 5.30211e-2, -1.00604e-2, 6.33457e-4};
	double B2[] = {0, 5.17449e-2, -1.21902e-1, 1.37748e-1, -7.32792e-2, 1.72914e-2, -1.38585e-3};

	double lambda_0 = (5.8808-6.1693e-2*T+3.4151e-4*T*T-3.0420e-7*T*T*T+1.2868e-10*T*T*T*T-2.1303e-14*T*T*T*T*T)/1000;

	for (int i = 1; i <= 6; i++)
	{
		sumresid += (B1[i]+B2[i]*T/crit.T)*pow(rho/crit.rho,i);
	}
	double lambda_r = sumresid;

	double lambda_c = this->conductivity_critical(T,rho,1/(6.2e-10)); //[W/m/K]

	return lambda_0 + lambda_r + lambda_c; //[W/m/K]
}

XenonClass::XenonClass()
{
    const double n[]={0.0,0.831150000,-2.355300000,0.539040000,0.014382000,0.066309000,0.000196490,0.149960000,-0.035319000,-0.159290000,-0.027521000,-0.023305000,0.008694100};
    const double u0[]={0.0,0};
    const double v0[]={0.0,0};

    // Critical parameters
    crit.rho = 1102.8612;
    crit.p = PressureUnit(5842.0, UNIT_KPA);
    crit.T = 289.733;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 131.293;
    params.Ttriple = 161.4;
	params.ptriple = 81.747799073227597;
    params.accentricfactor = 0.00363;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 161.4;
    limits.Tmax = 750.0;
    limits.pmax = 700000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-3.8227178129,3.8416395351);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(2.5-1);
    phi0list.push_back(phi0_logtau_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Xenon");
    aliases.push_back(std::string("Xe"));
    aliases.push_back(std::string("xenon"));
    aliases.push_back(std::string("XENON"));
    REFPROPname.assign("xenon");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";

}

R116Class::R116Class()
{
    const double n[]={0.0,1.163200000,-2.812300000,0.772020000,-0.143310000,0.102270000,0.000246290,0.308930000,-0.028499000,-0.303430000,-0.068793000,-0.027218000,0.010665000};
    const double u0[]={0.0,190.0,622.0,1470.0};
    const double v0[]={0.0,2.4818,7.0622,7.9951};

    // Critical parameters
    crit.rho = 613.32452808;
    crit.p = PressureUnit(3048.0, UNIT_KPA);
    crit.T = 293.03;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 138.01182;
    params.Ttriple = 173.1;
	params.ptriple = 26.0855875793;
    params.accentricfactor = 0.2566;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 173.1;
    limits.Tmax = 425.0;
    limits.pmax = 50000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-10.7088650331,8.9148979056);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("R116");

    REFPROPname.assign("R116");

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
void R116Class::ECSParams(double *e_k, double *sigma)
{
    *e_k = 226.16;
    *sigma = 0.5249;
}
double R116Class::ECS_f_int(double T)
{
    return 0.00132+0*T;
}
double R116Class::ECS_psi_viscosity(double rhor)
{
    return 1.21996-0.0647835*rhor+0*rhor*rhor;
}
double R116Class::ECS_chi_conductivity(double rhor)
{
    return 1.1804-0.0539975*rhor;
}

AcetoneClass::AcetoneClass()
{
    const double n[]={0.0,0.90041000000,-2.12670000000,-0.08340900000,0.06568300000,0.00016527000,-0.03966300000,0.72085000000,0.00923180000,-0.17217000000,-0.14961000000,-0.07612400000,-0.01816600000};
    const double u0[]={0.0,310.0,3480.0,1576.0};
    const double v0[]={0.0,3.7072,7.0675,11.012};

    // Critical parameters
    crit.rho = 272.971958;
    crit.p = PressureUnit(4700.0, UNIT_KPA);
    crit.T = 508.1;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 58.07914;
    params.Ttriple = 178.5;
	params.ptriple = 0.00232681797023;
    params.accentricfactor = 0.3071;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 178.5;
    limits.Tmax = 550.0;
    limits.pmax = 700000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-9.4883659997,7.1422719708);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);

    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Acetone");
    aliases.push_back(std::string("acetone")); 
    aliases.push_back(std::string("ACETONE"));
    REFPROPname.assign("acetone");

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";
}


NitrousOxideClass::NitrousOxideClass()
{
    const double n[]={0.0,0.88045000000,-2.42350000000,0.38237000000,0.06891700000,0.00020367000,0.13122000000,0.46032000000,-0.00369850000,-0.23263000000,-0.00042859000,-0.04281000000,-0.02303800000};
    const double u0[]={0.0,879.0,2372.0,5447.0};
    const double v0[]={0.0,2.1769,1.6145,0.48393};

    // Critical parameters
    crit.rho = 452.011456;
    crit.p = PressureUnit(7245.0, UNIT_KPA);
    crit.T = 309.52;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 44.0128;
    params.Ttriple = 182.33;
	params.ptriple = 87.837439225777985;
    params.accentricfactor = 0.1613;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 182.33;
    limits.Tmax = 525.0;
    limits.pmax = 50000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-4.4262736272,4.3120475243);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(3.5-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("NitrousOxide");
    aliases.push_back(std::string("N2O")); 
    aliases.push_back(std::string("NITROUSOXIDE"));
    REFPROPname.assign("N2O");

	ECSReferenceFluid = "Nitrogen";

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";
}


SulfurDioxideClass::SulfurDioxideClass()
{
    const double n[]={0.0,0.93061000000,-1.95280000000,-0.17467000000,0.06152400000,0.00017711000,0.21615000000,0.51353000000,0.01041900000,-0.25286000000,-0.05472000000,-0.05985600000,-0.01652300000};
    const double u0[]={0.0,775.0,1851.0};
    const double v0[]={0.0,1.0620,1.9401};

    // Critical parameters
    crit.rho = 525.002841;
    crit.p = PressureUnit(7884.0, UNIT_KPA);
    crit.T = 430.64;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 64.0638;
    params.Ttriple = 197.7;
	params.ptriple = 1.66036590338;
    params.accentricfactor = 0.2557;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 197.7;
    limits.Tmax = 525.0;
    limits.pmax = 35000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-4.5328346436,4.4777967379);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);

    phi_BC * phi0_power_ = new phi0_power(-7.2453e-05*pow(crit.T,1.0)/(1.0*(1.0+1)),-1.0);
	phi0list.push_back(phi0_power_);
    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("SulfurDioxide");
    aliases.push_back(std::string("SO2")); 
    aliases.push_back(std::string("SULFURDIOXIDE"));
    REFPROPname.assign("SO2");

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";
}


R141bClass::R141bClass()
{
    const double n[]={0.0,1.14690000000,-3.67990000000,1.34690000000,0.08332900000,0.00025137000,0.32720000000,0.46946000000,-0.02982900000,-0.31621000000,-0.02621900000,-0.07804300000,-0.02049800000};
    const double u0[]={0.0,502.0,1571.0,4603.0};
    const double v0[]={0.0,6.8978,7.8157,3.2039};

    // Critical parameters
    crit.rho = 458.55946002;
    crit.p = PressureUnit(4212.0,UNIT_KPA);
    crit.T = 477.5;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 116.94962;
    params.Ttriple = 169.68;
	params.ptriple = 0.00649365146247;
    params.accentricfactor = 0.2195;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 169.68;
    limits.Tmax = 500.0;
    limits.pmax = 400000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-15.5074814985,9.1871858933);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("R141b");
    aliases.push_back(std::string("R141B"));

    REFPROPname.assign("R141b");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
void R141bClass::ECSParams(double *e_k, double *sigma)
{
    *e_k = 370.44;
    *sigma = 0.5493;
}
double R141bClass::ECS_f_int(double T)
{
    return 0.000521722+0.00000292456*T;
}
double R141bClass::ECS_psi_viscosity(double rhor)
{
    return 0.92135+0.041091*rhor+0*rhor*rhor;
}
double R141bClass::ECS_chi_conductivity(double rhor)
{
    return 1.0867-0.0216469*rhor;
}


R142bClass::R142bClass()
{
    const double n[]={0.0,1.00380000000,-2.76620000000,0.42921000000,0.08136300000,0.00024174000,0.48246000000,0.75542000000,-0.00743000000,-0.41460000000,-0.01655800000,-0.10644000000,-0.02170400000};
    const double u0[]={0.0,473.0,1256.0,2497.0,6840.0};
    const double v0[]={0.0,5.0385,6.8356,4.0591,2.8136};

    // Critical parameters
    crit.rho = 445.99694314;
    crit.p = PressureUnit(4055.0, UNIT_KPA);
    crit.T = 410.26;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 100.49503;
    params.Ttriple = 142.72;
	params.ptriple = 0.00363327066489;
    params.accentricfactor = 0.2321;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 142.72;
    limits.Tmax = 470.0;
    limits.pmax = 60000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-12.6016527149,8.3160183265);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("R142b");
    aliases.push_back(std::string("R142B"));

    REFPROPname.assign("R142b");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
void R142bClass::ECSParams(double *e_k, double *sigma)
{
    *e_k = 278.2;
    *sigma = 0.5362;
}
double R142bClass::ECS_f_int(double T)
{
    return 0.000940725+0.000000988196*T;
}
double R142bClass::ECS_psi_viscosity(double rhor)
{
    return 0.9716+0.019181*rhor+0*rhor*rhor;
}
double R142bClass::ECS_chi_conductivity(double rhor)
{
    return 1.0749-0.0177916*rhor;
}


R218Class::R218Class()
{
    const double n[]={0.0,1.32700000000,-3.84330000000,0.92200000000,0.11360000000,0.00036195000,1.10010000000,1.18960000000,-0.02514700000,-0.65923000000,-0.02796900000,-0.18330000000,-0.02163000000};
    const double u0[]={0.0,326.0,595.0,1489.0};
    const double v0[]={0.0,7.2198,7.2692,11.599};

    // Critical parameters
    crit.rho = 627.9845622;
    crit.p = PressureUnit(2640.0, UNIT_KPA);
    crit.T = 345.02;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 188.01933;
    params.Ttriple = 125.45;
	params.ptriple = 0.00201898352904;
    params.accentricfactor = 0.3172;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 125.45;
    limits.Tmax = 440.0;
    limits.pmax = 20000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-15.6587335175,11.4531412796);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("R218");
    REFPROPname.assign("R218");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

void R218Class::ECSParams(double *e_k, double *sigma)
{
    *e_k = 266.35;
    *sigma = 0.58;
}
double R218Class::ECS_f_int(double T)
{
    return 0.000892659+0.00000114912*T;
}
double R218Class::ECS_psi_viscosity(double rhor)
{
    return 1.10225-0.00550442*rhor+0*rhor*rhor;
}
double R218Class::ECS_chi_conductivity(double rhor)
{
    return 1.2877-0.0758811*rhor;
}


R245faClass::R245faClass()
{
    const double n[]={0.0,1.2904,
		                 -3.2154,
						  0.50693,
						  0.093148,
						  0.00027638,
						  0.71458,
						  0.87252,
						 -0.015077,
						 -0.40645,
						 -0.11701,
						 -0.13062,
						 -0.022952};
    const double u0[]={0.0,222.0,1010.0,2450.0};
    const double v0[]={0.0,5.5728,10.385,12.554};

    // Critical parameters
    crit.rho = 516.084569;
    crit.p = PressureUnit(3651.0, UNIT_KPA);
    crit.T = 427.16;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 134.04794;
    params.Ttriple = 171.05;
	params.ptriple = 0.0125165917597;
    params.accentricfactor = 0.3776;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 171.05;
    limits.Tmax = 440.0;
    limits.pmax = 200000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0; i < u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-13.4283638514,9.87236538);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);

    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
	TransportReference.assign("Using ECS\n\nSurface Tension:\nJames W Schmidt, Ernesto Carrillo-Nava, Michael R Moldover \"Partially halogenated hydrocarbons CHFCl-CF3, CF3-CH3, CF3-CHF-CHF2, CF3-CH2-CF3, CHF2-CF2-CH2F, CF3-CH2-CHF2, CF3-O-CHF2: critical temperature, refractive indices, surface tension and estimates of liquid, vapor and critical densities\" Fluid Phase Equilibria, Volume 122, Issues 1�2, 31 July 1996, Pages 187�206 http://dx.doi.org/10.1016/0378-3812(96)03044-0");

    name.assign("R245fa");
	aliases.push_back("R245FA");
	REFPROPname.assign("R245fa");

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";

}
//double R245faClass::viscosity_Trho(double T, double rho)
//{ 
//    /*
//    Fitting of shape factor curves to R134a data. This method is employed because solving
//    for the shape factors is computationally very expensive and not very nice
//    convergence behavior is experienced.  Thus we can use the ECS method,
//    but with about the execution time of a conventional viscosity correlation.
//
//    This function code was automatically generated by the fit_shape_factor.py 
//    script in the dev/ folder on Wednesday, 27. February 2013
//
//    Mean absolute errors of shape factor prediction:
//    theta = 0.292256 %
//    phi = 2.72008 %
//    */
//
//    double e_k, sigma, tau, delta, A1, A2, A3, A4, theta, Tc, Tc0, T0, rho0; 
//    double DELTA, PSI_theta, psi, f, h, F_eta, M, M0, delta_omega, rho0bar;
//    double B1, B2, B3, B4, PSI_phi, Zc, Zc0, rhoc0, rhoc, log_tau, phi, rhobar;
//
//    double c[] = {-0.2154409562088915, -0.9709822762884642, 158.5843595248966, -74.52983550684725, 73672.6309941645, 73.74421821646304, -73827.2063003716, 47.22333447114561, -0.9472042408595509, 0.4964305453767688, 2.397163632251555, 2.396439386998316};
//    double d[] = {1.11465305331564, 0.2066346831910978, -212.4889952190534, 257.5843857690007, 285.3352731323313, -257.5533921176973, -82.24255997887354, -82.22887164856553, 0.9772522047882048, 0.5746074769208611, 2.248907052144205, 2.674375937148288};
//
//    tau = reduce.T/T;
//    delta = rho/reduce.rho;
//
//	long iR134a = get_Fluid_index(std::string("R134a"));
//    Fluid* R134a = get_fluid(iR134a);
//    delta_omega = params.accentricfactor-R134a->params.accentricfactor;
//
//    Zc = reduce.p/(reduce.rho*R()*reduce.T);
//    Zc0 = R134a->reduce.p/(R134a->reduce.rho*R134a->R()*R134a->reduce.T);
//    Tc = reduce.T;
//    Tc0 = R134a->reduce.T;
//    rhoc = reduce.rho;
//    rhoc0 = R134a->reduce.rho;
//    M = params.molemass;
//    M0 = R134a->params.molemass;
//
//    rhobar = rho/M;
//
//    if (rho > 40.858)
//    {
//        DELTA = pow(delta-1,2)+pow(tau-1,2);
//        log_tau = log(tau);
//
//        A1 = c[0]-c[1]*log_tau;
//        A2 = c[2]-c[3]*log_tau;
//        A3 = c[4]-c[5]*log_tau;
//        A4 = c[6]-c[7]*pow(log_tau,2);
//        PSI_theta = c[8]*delta*exp(-c[9]*pow(DELTA,2));
//        theta = 1+(delta_omega)*(A1+A2*exp(-pow(delta,2))+A3*exp(-pow(delta,c[10]))+A4*exp(-pow(delta,c[11]))+PSI_theta);
//
//        B1 = d[0]-d[1]*log_tau;
//        B2 = d[2]-d[3]*log_tau;
//        B3 = d[4]-d[5]*log_tau;
//        B4 = d[6]-d[7]*pow(log_tau,2);
//        PSI_phi = d[8]*delta*exp(-d[9]*pow(DELTA,2));
//        phi = Zc0/Zc*(1+(delta_omega)*(B1+B2*exp(-pow(delta,2))+B3*exp(-pow(delta,d[10]))+B4*exp(-pow(delta,d[11]))+PSI_phi));
//    }
//    else
//    {
//        theta = 1.0; phi = 1.0;
//    }
//    T0 = T*Tc0/theta/Tc;
//    h = M/M0*rhoc0/rhoc*phi;
//    rho0bar = rhobar*h;
//    rho0 = M0*rho0bar;
//
//    psi = ECS_psi_viscosity(delta);
//    f = T/T0;
//    F_eta = sqrt(f)*pow(h,-2.0/3.0)*sqrt(M/M0);
//    ECSParams(&e_k,&sigma);
//    return viscosity_dilute(T,e_k,sigma) + R134a->viscosity_background(T0,rho0*psi)*F_eta;
//} 
void R245faClass::ECSParams(double *e_k, double *sigma)
{
	// From Huber (2003)
    *e_k = 329.72; *sigma = 0.5529;
}
double R245faClass::ECS_f_int(double T)
{
	// From Huber (2003)
    return 1.64999e-3 - 3.28868e-7*T;
}
double R245faClass::ECS_psi_viscosity(double rhor)
{
	// From Huber (2003)
    return 1.1529 - 4.41540e-2*rhor;
}
double R245faClass::ECS_chi_conductivity(double rhor)
{
	// From Huber (2003)
    return 1.1627-0.0473491*rhor;
}
double R245faClass::surface_tension_T(double T)
{
	// From Mulero, 2012, JPCRD
	return 0.073586*pow(1-T/reduce.T,1.0983)+0.0103*pow(1-T/reduce.T,0.60033)-0.02663*pow(1-T/reduce.T,0.72765);
}

R41Class::R41Class()
{
    const double n[]={0.0,
					0.85316,
					-2.6366,
					0.69129,
					0.054681,
					0.00012796,
					-0.37093,
					0.33920,
					-0.0017413,
					-0.095417,
					-0.078852,
					-0.030729,
					-0.011497};
    const double u0[]={0.0,1841.0,4232.0};
    const double v0[]={0.0,5.6936,2.9351};

    // Critical parameters
    crit.rho = 316.506156;
    crit.p = PressureUnit(5897.0, UNIT_KPA);
    crit.T = 317.28;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 34.03292;
    params.Ttriple = 129.82;
	params.ptriple = 0.344246534824;
    params.accentricfactor = 0.2004;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 129.82;
    limits.Tmax = 425.0;
    limits.pmax = 70000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-4.867644116,4.2527951258);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);

    phi_BC * phi0_power_ = new phi0_power(-0.00016937*pow(crit.T,1.0)/(1.0*(1.0+1)),-1.0);
	phi0list.push_back(phi0_power_);
    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("R41");

    REFPROPname.assign("R41");

	BibTeXKeys.EOS = "Lemmon-JCED-2006";
	BibTeXKeys.ECS_LENNARD_JONES = "Chichester-NIST-2008";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

// --------------------
// ANCILLARY EQUATIONS
// --------------------
// All of these ancillary equations were generated with the Python script in the dev folder
// REFPROP saturation data (since it could go pretty much all the way (within 1 uK) to the critical point) was fit
// Powell method used to guess the exponents, within the function ODRPack from scipy.odr used to find the coefficients
// Works quite well

double R141bClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 509.965+1091.12*pow(THETA,0.439413)+185.865*pow(THETA,2.64809);
}
double IsohexaneClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 238.626+578.165*pow(THETA,0.398393)+109.839*pow(THETA,2.87839);
}
double R245faClass::rhosatL(double T)
{
double theta = 1-T/reduce.T;
double RHS,rho;

// Max error is 0.941705 %
RHS = +1.934009*pow(theta,0.333333)-1.048842*pow(theta,0.666667)+0.705519*pow(theta,1.333333)-0.564935*pow(theta,3.000000)+1.000168*pow(theta,6.166667);
rho = exp(RHS)*reduce.rho;
return rho;
}

double AcetoneClass::psat(double T)
{
    // Maximum absolute error is 0.161624 % between 178.500001 K and 508.099999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.7420933491645938, 2.2551227728042176, -2.2097584451804844, -0.80807019454396622, -3.0130241547201058, 1.1288610270689932 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double AcetoneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.479021 % between 178.500001 K and 508.099999 K
    const double ti[]={0,0.26459772663939418, 1.4091452766461503, 1.4202544602432625, 11.775226488952327, 1.4009455029990348};
    const double Ni[]={0,1.2349786870374411, -3379.0389866312494, 1436.4764562223259, 0.16448509887442936, 1942.7394651540521};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double AcetoneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.166474 % between 178.500001 K and 508.099999 K
    const double ti[]={0,0.32486954786675754, 1.1136505233132008, 1.397124100888042, 1.4206102336811905, 4.1261683284501736};
    const double Ni[]={0,-2.0054679941104143, -27.132764899927583, 304.9671207804098, -282.29487966244773, -4.08636533228474};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double NitrousOxideClass::psat(double T)
{
    // Maximum absolute error is 0.008121 % between 182.330001 K and 309.519999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-6.868769329455664, 2.044618171539057, -2.6329617846296989, 2.5383085882276508, -9.5683854957182639, 10.373557988632387 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double NitrousOxideClass::rhosatL(double T)
{
    // Maximum absolute error is 1.073634 % between 182.330001 K and 309.519999 K
    const double ti[]={0,0.60557790700716319, 0.70939607919869829, 2.5963911853523869, 0.68042926228856082, 0.67155910920440876};
    const double Ni[]={0,389.77875979761046, -1767.6976331438359, 0.20467533761004733, 10361.905167320929, -8982.8239481801156};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double NitrousOxideClass::rhosatV(double T)
{
    // Maximum absolute error is 0.610881 % between 182.330001 K and 309.519999 K
    const double ti[]={0,0.39081966725676398, 1.1522466533120135, 2.3180840594503693, 2.6358258974386621, 8.1917571767482542};
    const double Ni[]={0,-2.7894545753779387, -3.4744736926456077, 9.951411057303936, -11.242955195846788, -9.9192629053042864};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double SulfurDioxideClass::psat(double T)
{
    // Maximum absolute error is 0.110395 % between 197.700001 K and 430.639999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.4278303501775333, 2.4740630627264735, -3.5299527654256115, 2.6573726552982837, -8.7531566071245539, 6.5589471788208602 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double SulfurDioxideClass::rhosatL(double T)
{
    // Maximum absolute error is 1.149548 % between 197.700001 K and 430.639999 K
    const double ti[]={0,0.70876204235579221, 0.7229678532566155, 0.83486766193741779, 0.73141949820423335, 1.6738565688812563};
    const double Ni[]={0,12156.409552558353, -37170.600921873636, -391.93020566978561, 25406.018688280838, 1.5016221145583379};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double SulfurDioxideClass::rhosatV(double T)
{
    // Maximum absolute error is 0.772197 % between 197.700001 K and 430.639999 K
    const double ti[]={0,0.40597757408572105, 1.0078158710335801, 1.5968694256951719, 4.007078507058492, -0.16145660992043773};
    const double Ni[]={0,-3.0032294578395011, -2.4869182335281512, -0.13184442384899894, -4.7263242155221183, -0.00033424139516086214};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double R141bClass::psat(double T)
{
    // Maximum absolute error is 0.175974 % between 169.680001 K and 477.499999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.2204561290947709, 2.2758349671865776, -2.7502808754918209, 0.45527436710325281, -3.6366202389018829, 0.018941343442678414 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
//double R141bClass::rhosatL(double T)
//{
//    // Maximum absolute error is 100.000000 % between 169.680001 K and 477.499999 K
//    const double ti[]={0,0.2828129519231074, -0.25200521806435544, -0.72225671323569796, 3.1715066594344008, 3.3149603333488011};
//    const double Ni[]={0,1.30652890669813, 0.00055365778690386611, -8.3776570451762927e-05, -0.8565619843893828, 0.95359468578186291};
//    double summer=0;
//    int i;
//    double theta;
//    theta=1-T/reduce.T;
//    for (i=1;i<=5;i++)
//    {
//        summer+=Ni[i]*pow(theta,ti[i]);
//    }
//    return reduce.rho*exp(summer);
//}
double R141bClass::rhosatV(double T)
{
    // Maximum absolute error is 0.284442 % between 169.680001 K and 477.499999 K
    const double ti[]={0,0.37287817354647435, 0.91276608962657058, 4.7497423008783093, 2.6894254703994478, 18.919444023057686};
    const double Ni[]={0,-2.6251044774898529, -2.4960108200116555, -3.9765002486540135, -1.4797257697369941, -13.683339624267827};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double R142bClass::psat(double T)
{
    // Maximum absolute error is 0.104651 % between 142.720001 K and 410.259999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.3079519349991271, 2.2608201769441369, -2.5771971051912246, 0.3158163336195175, -4.025672461018349, 1.4340055620931298 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double R142bClass::rhosatL(double T)
{
    // Maximum absolute error is 0.806206 % between 142.720001 K and 410.259999 K
    const double ti[]={0,0.27967246985799865, 2.2595545100614993, 2.2688822238233723, 2.2789373329095373, 15.311342436654607};
    const double Ni[]={0,1.3105419105059686, 819.04365787414179, -1594.0402800866095, 775.09185686337332, 0.41293568648240681};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double R142bClass::rhosatV(double T)
{
    // Maximum absolute error is 0.361112 % between 142.720001 K and 410.259999 K
    const double ti[]={0,0.37340891699198825, 0.83800239688015932, 2.7771025614151861, 4.7454937920345799, 9.028952964225164};
    const double Ni[]={0,-2.4392795173014732, -2.6846161393904771, -1.6167477168461746, -3.8473658756292215, 0.93171399760554896};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double R218Class::psat(double T)
{
    // Maximum absolute error is 0.098033 % between 125.450001 K and 345.019999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.8102674790768001, 2.6559310101692217, -3.3120240937942338, -0.1760076204907364, -3.3661799062965576, 0.76921236609229326 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double R218Class::rhosatL(double T)
{
    // Maximum absolute error is 1.122651 % between 125.450001 K and 345.019999 K
    const double ti[]={0,0.36912361663072235, 0.71091417114172772, 1.276402066765205, 2.1923155395132605, 2.874407938057681};
    const double Ni[]={0,2.1608111314789578, -1.4692977923406709, 0.94803741116228957, -0.70505497641399617, 0.47144925458636899};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double R218Class::rhosatV(double T)
{
    // Maximum absolute error is 0.066811 % between 125.450001 K and 345.019999 K
    const double ti[]={0,0.45990648911676563, 0.72073155438993686, 2.8543906505788574, 3.1271008410006345, 3.1486689815091231};
    const double Ni[]={0,-3.5237161790276961, -1.1802834239391882, -79.601787512253651, 1180.1487663875571, -1107.036693472731};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double R245faClass::psat(double T)
{
    // Maximum absolute error is 0.064905 % between 171.050001 K and 427.159999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.8747833396154343, 2.1483917502708074, -3.0988907177869018, -0.75337618597081146, -3.6008657396087567, 0.51065419034037485 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double R245faClass::rhosatV(double T)
{
    // Maximum absolute error is 0.236320 % between 171.050001 K and 427.159999 K
    const double ti[]={0,0.34115317949378737, 0.79355129238839972, 2.8038127551966383, 2.7353987486677047, 4.822133836311937};
    const double Ni[]={0,-2.0810252672377687, -3.4639909310829626, 1.9189382275384959, -5.1121366494430251, -4.0159777625303965};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double R41Class::psat(double T)
{
    // Maximum absolute error is 0.134323 % between 129.820001 K and 317.279999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.1501947108542199, 1.7685510326093266, -0.91167320234282201, -1.5386487669225608, -0.93016709851120405, -0.1884522931161588 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double R41Class::rhosatL(double T)
{
    // Maximum absolute error is 1.080025 % between 129.820001 K and 317.279999 K
    const double ti[]={0,0.34900292856984289, 0.62896434437355531, 1.1347417598669201, 1.9506054552664505, 2.8649776878900806};
    const double Ni[]={0,1.9328233782786501, -0.86609333320547166, 0.46954483433043431, -0.39626336460294198, 0.24164742846392567};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double R41Class::rhosatV(double T)
{
    // Maximum absolute error is 0.595603 % between 129.820001 K and 317.279999 K
    const double ti[]={0,0.42460168461232123, 1.0384822760759682, 3.4244642508329721, 5.7279055805230801, 2.533996044180983};
    const double Ni[]={0,-3.3953933375589633, -2.0551088914483819, -4.7134737225801491, -0.51245543874877442, 1.9051073656048443};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double CarbonMonoxideClass::psat(double T)
{
    // Maximum absolute error is 0.106378 % between 68.160001 K and 132.859999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-6.1148803182354188, 1.0194126886426906, 0.056580021372273664, -2.418806695158302, 2.035082448741151, -4.7228196982300803 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double CarbonMonoxideClass::rhosatL(double T)
{
    // Maximum absolute error is 0.898555 % between 68.160001 K and 132.859999 K
    const double ti[]={0,0.29742254451341815, 1.0554099042533132, 1.489758454772462, 3.9719137526372501, 3.961448705224472};
    const double Ni[]={0,1.2676977992280518, 0.09697196291781153, -0.13852431420071057, 0.058524336717331704, 0.055612102409353124};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double CarbonMonoxideClass::rhosatV(double T)
{
    // Maximum absolute error is 0.369651 % between 68.160001 K and 132.859999 K
    const double ti[]={0,0.3854343522767737, 1.0104603372467613, 2.1283755916162046, 2.9193957820543872, 6.6868135430695865};
    const double Ni[]={0,-2.3723744833990232, -2.8592470973387676, 2.7678810236798252, -3.7675026842177273, -2.5221867301156755};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double CarbonylSulfideClass::psat(double T)
{
    // Maximum absolute error is 0.083057 % between 134.300001 K and 378.769999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-6.6144340815954452, 2.0534257788168455, -2.0722088131457848, 1.1149167490331107, -3.8513264261829465, 1.1808252156160568 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double CarbonylSulfideClass::rhosatL(double T)
{
    // Maximum absolute error is 0.677071 % between 134.300001 K and 378.769999 K
    const double ti[]={0,0.39284393894021785, 0.38390015433790292, 2.2386517189152535, 2.2115659405613348, 27.119517181027163};
    const double Ni[]={0,-12.475292718905436, 13.76513349301672, 5.9184468000342676, -5.9019500519058772, 22.562403946162679};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double CarbonylSulfideClass::rhosatV(double T)
{
    // Maximum absolute error is 0.670265 % between 134.300001 K and 378.769999 K
    const double ti[]={0,0.4002371604431198, 1.0280542888849111, 4.7960989516998014, 2.1665471855290961, 2.6645774028673594};
    const double Ni[]={0,-2.8333597180190568, -2.2058837335661083, -2.8226141763428818, 1.6742917835056423, -2.2132884573899081};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double DecaneClass::psat(double T)
{
    // Maximum absolute error is 0.095897 % between 243.500001 K and 617.699999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-8.6751621424501106, 2.8804901074112048, -3.9783693016893227, -0.78550488926699347, -5.4428041796076343, 3.0822748884469573 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double DecaneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.397401 % between 243.500001 K and 617.699999 K
    const double ti[]={0,0.39698081733802165, 0.99375746839456702, 1.3824414058838479, 2.0870109314455987, 4.0369621174956682};
    const double Ni[]={0,2.306748788146411, -3.2750228326598099, 3.2390179121827032, -1.0158340181331245, 0.22333903274150235};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double DecaneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.192674 % between 243.500001 K and 617.699999 K
    const double ti[]={0,0.50769604931055112, 5.5518179080415102, 1.2630033323438727, 5.1843447790955253, 1.7840587243093533};
    const double Ni[]={0,-5.317247348628082, 26.915460654851557, 1.6082159690036035, -31.527804141299679, -4.4909283644660496};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double HydrogenSulfideClass::psat(double T)
{
    // Maximum absolute error is 0.014061 % between 187.700001 K and 373.099999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-6.5723352326794293, 1.896430781127326, -1.9355862944610984, 1.2531037764282309, -4.795728167798825, 2.787887611896148 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double HydrogenSulfideClass::rhosatL(double T)
{
    // Maximum absolute error is 0.566702 % between 187.700001 K and 373.099999 K
    const double ti[]={0,0.79106386489177882, 0.88939381375621784, 1.0266714393290504, 0.91099311495488788, 1.0191641121567934};
    const double Ni[]={0,738.4889085778741, -12085.6445805242, 16689.69047856304, 14107.568427822798, -19448.779316224445};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double HydrogenSulfideClass::rhosatV(double T)
{
    // Maximum absolute error is 0.601824 % between 187.700001 K and 373.099999 K
    const double ti[]={0,0.53435963955142052, 1.0318813376408911, 1.0297572598900648, 1.0392102012903406, 4.3435423458517368};
    const double Ni[]={0,-5.7597583855079861, -23677.45852297284, 18468.602400239459, 5209.7832429802875, -3.4496094719154797};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double IsopentaneClass::psat(double T)
{
    // Maximum absolute error is 0.859211 % between 112.650001 K and 460.349999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.1614631214691062, 1.6427675982493801, -1.0765148777064892, -2.4350359463228068, -0.00077830587755741939, -1.3607605842862582 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double IsopentaneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.412680 % between 112.650001 K and 460.349999 K
    const double ti[]={0,0.2910137977462911, 2.0204585261825168, 4.1882444640378527, 3.0704844993207829, 9.097169407695036};
    const double Ni[]={0,1.3160559423567124, -0.24235357370003638, -0.2239827299027437, 0.45658057492355336, 0.10114204215110814};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double IsopentaneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.314179 % between 112.650001 K and 460.349999 K
    const double ti[]={0,0.31270615359363657, 0.74616825916589047, 13.438145185545892, 3.7712973823560603, 13.112073494833108};
    const double Ni[]={0,-1.6095085105632891, -3.5410442212702162, -2.4025730674501111, -4.4705198913818949, 0.70459696994457377};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double NeopentaneClass::psat(double T)
{
    // Maximum absolute error is 0.004526 % between 256.600001 K and 433.739999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.027004544708106, 1.9793386200106815, -2.360943757602449, 1.1636391968268969, -6.9110737915478353, 6.6905802185479404 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double NeopentaneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.308171 % between 256.600001 K and 433.739999 K
    const double ti[]={0,0.30715908929085167, 0.02408980306294933, 2.2951165419376474, 4.2880092045078628, 3.5718062728106448};
    const double Ni[]={0,1.3278912346796319, -7.1433703381307663e-07, -0.73279809856370448, -2.8725817580875623, 3.0804108412532227};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double NeopentaneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.331471 % between 256.600001 K and 433.739999 K
    const double ti[]={0,0.34885244286586209, 0.96026335525006201, 1.7026350507007288, 2.6663611861234626, 8.6409603088688787};
    const double Ni[]={0,-2.2227708856840378, -4.1741890159176576, 3.3061543791847985, -4.7204205710823963, -14.871952577080259};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double IsohexaneClass::psat(double T)
{
    // Maximum absolute error is 0.867155 % between 119.600001 K and 497.699999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.3903287973234741, 1.5420871354898538, -0.84149185676025828, -3.6210436902930638, 0.99970651583741965, -1.9923902966792624 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double IsohexaneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.605575 % between 119.600001 K and 497.699999 K
    const double ti[]={0,0.47580852830565717, 0.80170994863058465, 4.1456753389418717, 2.0754191791227719, 22.863511521748787};
    const double Ni[]={0,-3.7956367210298922, -1.0172685447655876, -4.3120611056906561, -1.4909297795395249, -10.739883495485737};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double KryptonClass::psat(double T)
{
    // Maximum absolute error is 0.017436 % between 115.775001 K and 209.479999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-5.9739950665662702, 1.3051794187375625, -0.31592680234061876, -1.0365678166366274, 0.63351382066291007, -2.8912024898026756 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double KryptonClass::rhosatL(double T)
{
    // Maximum absolute error is 1.463881 % between 115.775001 K and 209.479999 K
    const double ti[]={0,0.36780895915580791, 0.70248203166996614, 1.172529422523406, 6.2635912317680793, 1.7846540318478548};
    const double Ni[]={0,1.9755226451868064, -1.3056130307532876, 0.8586443600865078, 0.12193064253237056, -0.30634286614817541};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double KryptonClass::rhosatV(double T)
{
    // Maximum absolute error is 0.283547 % between 115.775001 K and 209.479999 K
    const double ti[]={0,0.453854278772161, 1.6977667091236079, 2.0472835465021171, 2.0651987529414351, 7.5969412041939615};
    const double Ni[]={0,-3.1505671898070267, -18.65209521960794, 400.42457764318993, -383.96236538350985, -3.3751290883620086};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double NonaneClass::psat(double T)
{
    // Maximum absolute error is 0.087879 % between 219.700001 K and 594.549999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-8.4747418596031192, 2.8227667935992362, -3.5458599565380906, -1.4242947475868946, -3.5254946604872481, 1.1334065831048079 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double NonaneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.867150 % between 219.700001 K and 594.549999 K
    const double ti[]={0,0.45434458789525356, 0.72748787625567024, 0.925063834736367, 1.3573478933523566, 3.0112206747240773};
    const double Ni[]={0,5.1572118610805857, -11.028370080657231, 8.9470309229999039, -1.8357217413258595, 0.20683170461251787};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double NonaneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.312530 % between 219.700001 K and 594.549999 K
    const double ti[]={0,0.50655477629124324, 2.4876505608371025, 2.2041597822504291, 1.9430824560248301, 2.5830185280832492};
    const double Ni[]={0,-5.2194783202879194, 861.04350631411342, -342.56220842639533, 83.891172356485143, -610.13190227150676};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double TolueneClass::psat(double T)
{
    // Maximum absolute error is 0.187697 % between 178.000001 K and 591.749999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.4184042521750229, 1.7705795922329421, -1.0620069658107303, -2.8961055438497469, -0.35601288602065423, -1.1203832074545026 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double TolueneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.659779 % between 178.000001 K and 591.749999 K
    const double ti[]={0,0.32703234036345474, 0.86030837877972832, 1.9855638765309758, 4.3159296915209184, 5.9759915024218611};
    const double Ni[]={0,1.585216471761786, -0.42546736843083033, 0.23353047360987148, -0.1034086810731789, 0.13104469387313702};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double TolueneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.270124 % between 178.000001 K and 591.749999 K
    const double ti[]={0,0.26453133585430771, 0.59439611774229872, 3.2399544120904609, 5.3573873604955704, 13.640514113563155};
    const double Ni[]={0,-0.56008642910342632, -4.465640983522734, -3.7664413431830686, -1.9809221774094345, -0.062332971201866233};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}


double XenonClass::psat(double T)
{
    // Maximum absolute error is 0.007007 % between 161.405001 K and 289.732999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-6.0179736298621336, 1.4531608028462166, -0.75746733733446314, -0.05721415196027424, -1.4372963787937452, -0.28767071914740688 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double XenonClass::rhosatL(double T)
{
    // Maximum absolute error is 1.642288 % between 161.405001 K and 289.732999 K
    const double ti[]={0,0.31434776766743999, 0.67805601917916636, 5.8038130756111892, 2.3240256223897999, 3.0208547327754598};
    const double Ni[]={0,1.4226874845199544, -0.21979298723304491, 0.32922435816070406, 0.26160336954409047, -0.32636485566039475};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double XenonClass::rhosatV(double T)
{
    // Maximum absolute error is 0.384221 % between 161.405001 K and 289.732999 K
    const double ti[]={0,0.42456694647604343, 1.0929357777332234, 2.7453780242937951, 2.5165708451263624, 11.065026375451158};
    const double Ni[]={0,-2.8060599998539235, -2.0481529361896662, -9.9196250401113062, 9.0109861593982856, -19.388336288784551};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

double R116Class::psat(double T)
{
    // Maximum absolute error is 0.005177 % between 173.100001 K and 293.029999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.3955589482937381, 2.1874345047641488, -2.6407260335629017, 0.83237779144881274, -6.480800530357623, 6.0527091022412938 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double R116Class::rhosatL(double T)
{
    // Maximum absolute error is 0.976437 % between 173.100001 K and 293.029999 K
    const double ti[]={0,0.29011688981180672, 0.34284809214776468, 2.1525665789650317, 3.0547917391285706, 3.2820127321141621};
    const double Ni[]={0,1.5085971874674042, -0.19169905028904163, -0.1004169669731468, 0.068745519590789839, 0.10395251177877873};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double R116Class::rhosatV(double T)
{
    // Maximum absolute error is 0.309692 % between 173.100001 K and 293.029999 K
    const double ti[]={0,0.38604749962219331, 1.0088719781777016, 1.6277698849200228, 2.720738012316323, 6.8941744996057599};
    const double Ni[]={0,-2.6038986210545034, -4.0396795352002455, 2.7069117508619911, -4.2648317254667782, -5.93238493168551};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
