


#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Span_Short.h"

static double d_nonpolar_SpanShort[] =
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

static double t_nonpolar_SpanShort[] =
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

static double c_nonpolar_SpanShort[] =
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

static double d_polar_SpanShort[] =
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

static double t_polar_SpanShort[] =
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

static double c_polar_SpanShort[] =
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

nPentaneClass::nPentaneClass()
{
const double n[] = {0.0, 1.0968643, -2.9988888, 0.99516887, -0.16170709, 0.1133446, 0.00026760595, 0.40979882, -0.040876423, -0.38169482, -0.10931957, -0.032073223, 0.016877016};
std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
std::vector<double> t_v(t_nonpolar_SpanShort,t_nonpolar_SpanShort+sizeof(t_nonpolar_SpanShort)/sizeof(double));
std::vector<double> d_v(d_nonpolar_SpanShort,d_nonpolar_SpanShort+sizeof(d_nonpolar_SpanShort)/sizeof(double));
std::vector<double> l_v(c_nonpolar_SpanShort,c_nonpolar_SpanShort+sizeof(c_nonpolar_SpanShort)/sizeof(double));

//Critical parameters
crit.rho = 232; //[kg/m^3]
crit.p = PressureUnit(3370, UNIT_KPA); //[kPa]
crit.T = 469.7; //[K]
crit.v = 1/crit.rho; 

// Other fluid parameters
params.molemass = 72.14878;
params.Ttriple = 143.47;
params.accentricfactor = 0.251;
params.R_u = 8.314510;
params.ptriple = 7.63221842256e-05;

// Limits of EOS
limits.Tmin = params.Ttriple;
limits.Tmax = 500.0;
limits.pmax = 100000.0;
limits.rhomax = 1000000.0*params.molemass;

phirlist.push_back(new phir_power( n_v,d_v,t_v,l_v,1,12));

double T0 = 309.213623675,
R_ = 8.314510/params.molemass,
rho0 = 609.710850486,
m,
c,
H0 = 0.0, // kJ/kmol
S0 = 0.0, /// kJ/kmol/K
tau0=crit.T/T0,
delta0=rho0/crit.rho;

// log(delta)+c+m*tau

/// c is the constant term
c=-S0/R_-1+log(tau0/delta0);/*<< from the leading term*/

/// m multiplies the tau term in the leading term (slope)
m=H0/(R_*crit.T); /*<< from the leading term */

phi0list.push_back(new phi0_lead(c,m));
phi0list.push_back(new phi0_logtau(-1.0));

const double a0[] = {0,4*params.R_u,8.95043*params.R_u,178.67, 21.836*params.R_u, 840.538};
std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
const double a1[] = {0,0,33.4032*params.R_u, 1774.25, 0, 0};
std::vector<double> a1_v(a1,a1+sizeof(a1)/sizeof(double));

phi0list.push_back(new phi0_cp0_AlyLee(a0_v,crit.T,T0,params.R_u));
phi0list.push_back(new phi0_cp0_AlyLee(a1_v,crit.T,T0,params.R_u));

EOSReference.assign("Span, R. and W. Wagner, \"Equations of State for Technical Applications. II. Results for Nonpolar Fluids\", International Journal of Thermophysics, Vol. 24, No. 1, January 2003. \n\nCp0: Jaeschke, M. and P. Schley,\"Ideal-Gas Thermodynamic Properties for Natural-Gas Applications\", Int. J. Thermophys. v 16, n6 1995");
TransportReference.assign("Using ECS in fully predictive mode");

name.assign("n-Pentane");
aliases.push_back(std::string("nPentane"));
aliases.push_back(std::string("Pentane"));
aliases.push_back(std::string("PENTANE"));
aliases.push_back(std::string("N-PENTANE"));
aliases.push_back(std::string("R601"));
REFPROPname.assign("PENTANE");

BibTeXKeys.EOS = "Span-IJT-2003B";
BibTeXKeys.CP0 = "Jaeschke-IJT-1995";
BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";
BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";

}

double nPentaneClass::psat(double T)
{
    // Maximum absolute error is 0.214979 % between 143.470001 K and 469.699990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.3186030873160792, 1.8296494289307599, -1.5622713200315748, -1.6154955648222435, -1.5250180952272199, -0.17887604235363186 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double nPentaneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.735195 % between 143.470001 K and 469.699990 K
    const double ti[]={0,0.29637004389880028, 1.5209547193578159, 1.8091465611070348, 5.2452029430809954};
    const double Ni[]={0,1.3756520855776841, -0.56869128829842586, 0.5342216133477542, 0.035012032864392216};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double nPentaneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.595193 % between 143.470001 K and 469.699990 K
    const double ti[]={0,0.4021494191330745, 0.91494671969127384, 38.233907328347136, 3.8604517063013857};
    const double Ni[]={0,-2.8582545419999192, -2.5045842450420941, -2713.0646430695333, -4.7011651958462872};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

nHexaneClass::nHexaneClass()
{
const double n[] = {0.0, 1.0553238, -2.6120616, 0.76613883, -0.29770321, 0.11879908, 0.00027922861, 0.4634759, 0.011433197, -0.48256969, -0.093750559, -0.0067273247, -0.0051141584};
std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
std::vector<double> t_v(t_nonpolar_SpanShort,t_nonpolar_SpanShort+sizeof(t_nonpolar_SpanShort)/sizeof(double));
std::vector<double> d_v(d_nonpolar_SpanShort,d_nonpolar_SpanShort+sizeof(d_nonpolar_SpanShort)/sizeof(double));
std::vector<double> l_v(c_nonpolar_SpanShort,c_nonpolar_SpanShort+sizeof(c_nonpolar_SpanShort)/sizeof(double));

//Critical parameters
crit.rho = 233.18; //[kg/m^3]
crit.p = PressureUnit(3034, UNIT_KPA); //[kPa]
crit.T = 507.82; //[K]
crit.v = 1/crit.rho; 

// Other fluid parameters
params.molemass = 86.17536;
params.Ttriple = 177.83;
params.accentricfactor = 0.299;
params.R_u = 8.314510;
params.ptriple = 0.00127714043311;

// Limits of EOS
limits.Tmin = params.Ttriple;
limits.Tmax = 500.0;
limits.pmax = 100000.0;
limits.rhomax = 1000000.0*params.molemass;

phirlist.push_back(new phir_power( n_v,d_v,t_v,l_v,1,12));

double T0 = 341.864512243,
R_ = 8.314510/params.molemass,
rho0 = 613.01116119,
m,
c,
H0 = 0.0, // kJ/kmol
S0 = 0.0, /// kJ/kmol/K
tau0=crit.T/T0,
delta0=rho0/crit.rho;

// log(delta)+c+m*tau

/// c is the constant term
c=-S0/R_-1+log(tau0/delta0);/*<< from the leading term*/

/// m multiplies the tau term in the leading term (slope)
m=H0/(R_*crit.T); /*<< from the leading term */

phi0list.push_back(new phi0_lead(c,m));
phi0list.push_back(new phi0_logtau(-1.0));

const double a0[] = {0,4*params.R_u,11.6977*params.R_u,182.326, 26.8142*params.R_u, 859.207};
std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
const double a1[] = {0,0,38.6164*params.R_u, 1826.59, 0, 0};
std::vector<double> a1_v(a1,a1+sizeof(a1)/sizeof(double));

phi0list.push_back(new phi0_cp0_AlyLee(a0_v,crit.T,T0,params.R_u));
phi0list.push_back(new phi0_cp0_AlyLee(a1_v,crit.T,T0,params.R_u));

EOSReference.assign("Span, R. and W. Wagner, \"Equations of State for Technical Applications. II. Results for Nonpolar Fluids\", International Journal of Thermophysics, Vol. 24, No. 1, January 2003. \n\nCp0: Jaeschke, M. and P. Schley,\"Ideal-Gas Thermodynamic Properties for Natural-Gas Applications\", Int. J. Thermophys. v 16, n6 1995");
TransportReference.assign("Using ECS in fully predictive mode");

name.assign("n-Hexane");
aliases.push_back("nHexane");
aliases.push_back("Hexane");
aliases.push_back(std::string("HEXANE"));
aliases.push_back(std::string("N-HEXANE"));
REFPROPname.assign("HEXANE");

BibTeXKeys.EOS = "Span-IJT-2003B";
BibTeXKeys.CP0 = "Jaeschke-IJT-1995";
BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";
BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
BibTeXKeys.VISCOSITY = "Michailidou-JPCRD-2013";
BibTeXKeys.CONDUCTIVITY = "Assael-JPCRD-2013B";

}
double nHexaneClass::psat(double T)
{
	// Max error is  0.105432451163 % between 177.83 and 507.819999 K
    const double t[]={0, 0.051000000000000004, 1.0, 1.3333333333333333, 3.1666666666666665, 8.5, 9.666666666666666};
    const double N[]={0, 0.0052535021674954144, -7.8062705612610275, 1.6043453368140383, -3.550385816986847, -7.245501435072093, 7.4818124089530356};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double nHexaneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.266613 % between 177.830000 K and 507.820000 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.6666666666666667, 1.8333333333333333, 2.0, 2.3333333333333335, 2.6666666666666665, 3.0, 3.3333333333333335, 4.0, 4.333333333333333};
    const double N[] = {0, -0.36042485296131643, 43.664803867359119, -1804.2336876177826, 36769.690581511182, -434905.12293345085, 3273782.0693872692, -16583157.118874608, 58286580.722277559, -143415988.18061239, 243088854.89131206, -267312869.47766161, 156613144.04323533, -51245576.753353722, 26550511.264813472, -11851772.552882895, 3248032.5556567763, -331324.18203971104, 79683.643452126213};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=18; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double nHexaneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.369573 % between 177.830000 K and 507.820000 K
    const double t[] = {0, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.6666666666666667, 1.8333333333333333, 2.0, 2.3333333333333335, 2.6666666666666665, 3.0, 3.3333333333333335, 4.0, 4.333333333333333, 5.0, 5.5};
    const double N[] = {0, -13.32656487738959, 1246.0726749247744, -42990.264799784949, 757508.7800032777, -7876385.2128088353, 52425760.036635056, -234043333.88894936, 714747253.67619896, -1480735989.3985567, 1970815159.663389, -1389664970.180351, 655235798.12866044, -492889866.23758823, 326547190.99802482, -137859585.07323754, 41116916.874772295, -20862745.129438739, 2767339.869803539, -438323.38014475687};
    double summer=0,theta;
    theta=1-T/reduce.T;	
	for (int i=1; i<=19; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

double nHexaneClass::viscosity_Trho(double T, double rho)
{
	double e_k, sigma;
	this->ECSParams(&e_k,&sigma);
	double Tstar = T/e_k, rhor = rho/233.182, Tr = T/reduce.T;

	double log_Tstar = log(Tstar);
	
	// Dilute gas
	double eta_0 = 0.021357*sqrt(params.molemass*T)/(sigma*sigma*exp(0.18760-0.48430*log_Tstar+0.04477*log_Tstar*log_Tstar)); // uPa-s

	// Initial density dependence from Rainwater-Friend
	double b[] = {-19.572881, 219.73999, -1015.3226, 2471.01251, -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158};
	double Bstar = b[0]*pow(Tstar,-0.25*0)+b[1]*pow(Tstar,-0.25*1)+b[2]*pow(Tstar,-0.25*2)+b[3]*pow(Tstar,-0.25*3)+b[4]*pow(Tstar,-0.25*4)+b[5]*pow(Tstar,-0.25*5)+b[6]*pow(Tstar,-0.25*6)+b[7]*pow(Tstar,-2.5)+b[8]*pow(Tstar,-5.5);
	double B = Bstar*0.60221415*sigma*sigma*sigma; // L/mol

	double c[] = {2.53402335, -9.724061002, 0.469437316, 158.5571631, 72.42916856, 10.60751253, 8.628373915, -6.61346441, -2.212724566};
	double eta_r = pow(rhor,2.0/3.0)*sqrt(Tr)*(c[0]/Tr+c[1]/(c[2]+Tr+c[3]*rhor*rhor)+c[4]*(1+rhor)/(c[5]+c[6]*Tr+c[7]*rhor+rhor*rhor+c[8]*rhor*Tr));

	double rhobar = rho/params.molemass; // [mol/L]
	return (eta_0*(1+B*rhobar) + eta_r)/1e6;
}
double nHexaneClass::conductivity_Trho(double T, double rho)
{
	double Tr = T/reduce.T;
	double lambda_0 = (6.6742-23.7619*Tr+72.0155*Tr*Tr-18.3714*Tr*Tr*Tr)/1000; // W/m/K

	double sumresid = 0;
	double B1[] = {0, -3.01408e-2, 1.67975e-1, -1.29739e-1, 3.82833e-2, -3.70294e-3};
	double B2[] = {0, 2.18208e-2, -1.00833e-1, 7.74180e-2, -2.15945e-2, 2.12487e-3};

	for (int i = 1; i<= 5; i++){		
		sumresid += (B1[i]+B2[i]*(T/reduce.T))*pow(rho/reduce.rho,i);
	}

	double lambda_r = sumresid; // [W/m/K]

	double lambda_c = this->conductivity_critical(T,rho,1.0/(7.37e-10)); // [W/m/K]

	return lambda_0+lambda_r+lambda_c;
}


nHeptaneClass::nHeptaneClass()
{
const double n[] = {0.0, 1.0543748, -2.6500682, 0.81730048, -0.30451391, 0.12253869, 0.00027266473, 0.49865826, -0.00071432815, -0.54236896, -0.13801822, -0.0061595287, 0.0004860251};
std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
std::vector<double> t_v(t_nonpolar_SpanShort,t_nonpolar_SpanShort+sizeof(t_nonpolar_SpanShort)/sizeof(double));
std::vector<double> d_v(d_nonpolar_SpanShort,d_nonpolar_SpanShort+sizeof(d_nonpolar_SpanShort)/sizeof(double));
std::vector<double> l_v(c_nonpolar_SpanShort,c_nonpolar_SpanShort+sizeof(c_nonpolar_SpanShort)/sizeof(double));

//Critical parameters
crit.rho = 232; //[kg/m^3]
crit.p = PressureUnit(2736, UNIT_KPA); //[kPa]
crit.T = 540.13; //[K]
crit.v = 1/crit.rho; 

// Other fluid parameters
params.molemass = 100.202;
params.Ttriple = 182.55;
params.accentricfactor = 0.349;
params.R_u = 8.314510;
params.ptriple = 0.00017548675253297369;

// Limits of EOS
limits.Tmin = params.Ttriple;
limits.Tmax = 500.0;
limits.pmax = 100000.0;
limits.rhomax = 1000000.0*params.molemass;

phirlist.push_back(new phir_power( n_v,d_v,t_v,l_v,1,12));

double T0 = 371.533277446,
R_ = 8.314510/params.molemass,
rho0 = 614.215551363,
m,
c,
H0 = 0.0, // kJ/kmol
S0 = 0.0, /// kJ/kmol/K
tau0=crit.T/T0,
delta0=rho0/crit.rho;

// log(delta)+c+m*tau

/// c is the constant term
c=-S0/R_-1+log(tau0/delta0);/*<< from the leading term*/

/// m multiplies the tau term in the leading term (slope)
m=H0/(R_*crit.T); /*<< from the leading term */

phi0list.push_back(new phi0_lead(c,m));
phi0list.push_back(new phi0_logtau(-1.0));

const double a0[] = {0,4*params.R_u,13.7266*params.R_u,169.789, 30.4707*params.R_u, 836.195};
std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
const double a1[] = {0,0,43.5561*params.R_u, 1760.46, 0, 0};
std::vector<double> a1_v(a1,a1+sizeof(a1)/sizeof(double));

phi0list.push_back(new phi0_cp0_AlyLee(a0_v,crit.T,T0,params.R_u));
phi0list.push_back(new phi0_cp0_AlyLee(a1_v,crit.T,T0,params.R_u));

EOSReference.assign("Span, R. and W. Wagner, \"Equations of State for Technical Applications. II. Results for Nonpolar Fluids\", International Journal of Thermophysics, Vol. 24, No. 1, January 2003. \n\nCp0: Jaeschke, M. and P. Schley,\"Ideal-Gas Thermodynamic Properties for Natural-Gas Applications\", Int. J. Thermophys. v 16, n6 1995");
TransportReference.assign("Using ECS in fully predictive mode");

name.assign("n-Heptane");
aliases.push_back("nHeptane");
aliases.push_back("Heptane");
aliases.push_back(std::string("HEPTANE"));
aliases.push_back(std::string("N-HEPTANE"));
REFPROPname.assign("HEPTANE");

BibTeXKeys.EOS = "Span-IJT-2003B";
BibTeXKeys.CP0 = "Jaeschke-IJT-1995";
BibTeXKeys.ECS_LENNARD_JONES = "Chichester-NIST-2008";
BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
BibTeXKeys.CONDUCTIVITY = "Assael-JPCRD-2013C";
}
double nHeptaneClass::psat(double T)
{
    // Maximum absolute error is 0.083195 % between 182.550001 K and 540.129990 K
    const double t[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 17, 19, 22, 24, 26, 29};
    const double N[]={0, -3.04368483534491, 221.51326813318423, -7529.6367436332157, 143537.33433041952, -1782600.4037331331, 15351794.790984498, -95064479.531315014, 431824604.4413681, -1449020571.7368731, 3568070938.3565359, -6279101701.5130682, 7384006525.7245655, -4740322108.6224594, 1972249451.6581755, -1368254214.5357554, 801108370.43087959, -468263339.91664481, 327751566.67380428, -109115983.60517003, 10435955.348631999};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=20;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double nHeptaneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.033330 % between 182.550001 K and 540.129990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.6666666666666667, 2.0, 2.3333333333333335};
    const double N[] = {0, 8.521366797649824, -186.43238448180225, 1821.3194470035153, -9681.71039328181, 30640.757127668767, -59194.633606025527, 67519.411370724614, -37969.577087502075, 9354.1453367821832, -2755.5351214591965, 446.92980043327628};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
	for (i=1; i<=11; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double nHeptaneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.169613 % between 182.550001 K and 540.129990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.6666666666666667, 2.0, 2.3333333333333335};
    const double N[] = {0, -7.9622143721199228, 113.40353424908179, -848.80926796934193, 3848.254633388874, -11203.055218321815, 21063.716927422713, -24590.440523030087, 14797.521198915152, -4622.1556877203093, 1870.6514096531255, -432.72469128767966};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    	
	for (i=1; i<=11; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}
double nHeptaneClass::conductivity_Trho(double T, double rho)
{
	double Tr = T/reduce.T;
	double lambda_0 = ((-1.83367 + 16.2572*Tr - 39.0996*Tr*Tr + 47.8594*Tr*Tr*Tr + 15.1925*Tr*Tr*Tr*Tr - 3.39115*Tr*Tr*Tr*Tr*Tr)/(0.250611 - 0.320871*Tr + Tr*Tr))/1000; // W/m/K

	double sumresid = 0;
	double B1[] = {0, 5.17785e-2, -9.24052e-2, 5.11484e-2, -7.76896e-3, 1.21637e-4};
	double B2[] = {0, -7.72433e-3, 2.18899e-2, 1.71725e-3, -7.91642e-3, 1.83379e-3};

	for (int i = 1; i<= 5; i++){
		sumresid += (B1[i]+B2[i]*(T/reduce.T))*pow(rho/reduce.rho,i);
	}

	double lambda_r = sumresid; // [W/m/K]

	double lambda_c = this->conductivity_critical(T,rho,1.0/(8.0e-10),0.0586,2.45e-10); // [W/m/K]

	return lambda_0+lambda_r+lambda_c;
}

nOctaneClass::nOctaneClass()
{
	const double n[] = {0.0, 1.0722545, -2.4632951, 0.65386674, -0.36324974, 0.1271327, 0.00030713573, 0.52656857, 0.019362863, -0.58939427, -0.14069964, -0.0078966331, 0.0033036598};
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> t_v(t_nonpolar_SpanShort,t_nonpolar_SpanShort+sizeof(t_nonpolar_SpanShort)/sizeof(double));
	std::vector<double> d_v(d_nonpolar_SpanShort,d_nonpolar_SpanShort+sizeof(d_nonpolar_SpanShort)/sizeof(double));
	std::vector<double> l_v(c_nonpolar_SpanShort,c_nonpolar_SpanShort+sizeof(c_nonpolar_SpanShort)/sizeof(double));

	//Critical parameters
	crit.rho = 2.0564*114.2285; //[kg/m^3]
	crit.p = PressureUnit(2497, UNIT_KPA); //[kPa]
	crit.T = 569.32; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 114.2285;
	params.Ttriple = 216.37;
	params.accentricfactor = 0.395;
	params.R_u = 8.314510;
	params.ptriple = 0.0019888768508328075;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n_v,d_v,t_v,l_v,1,12));

	double T0 = 398.773831136,
	R_ = 8.314510/params.molemass,
	rho0 = 612.212813966,
	m,
	c,
	H0 = 0.0, // kJ/kmol
	S0 = 0.0, /// kJ/kmol/K
	tau0=crit.T/T0,
	delta0=rho0/crit.rho;

	// log(delta)+c+m*tau

	/// c is the constant term
	c=-S0/R_-1+log(tau0/delta0);/*<< from the leading term*/

	/// m multiplies the tau term in the leading term (slope)
	m=H0/(R_*crit.T); /*<< from the leading term */

	phi0list.push_back(new phi0_lead(c,m));
	phi0list.push_back(new phi0_logtau(-1.0));

	const double a0[] = {0,4*params.R_u,15.6865*params.R_u,158.922, 33.8029*params.R_u, 815.064};
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	const double a1[] = {0,0,48.1731*params.R_u, 1693.07, 0, 0};
	std::vector<double> a1_v(a1,a1+sizeof(a1)/sizeof(double));

	phi0list.push_back(new phi0_cp0_AlyLee(a0_v,crit.T,T0,params.R_u));
	phi0list.push_back(new phi0_cp0_AlyLee(a1_v,crit.T,T0,params.R_u));

	EOSReference.assign("Span, R. and W. Wagner, \"Equations of State for Technical Applications. II. Results for Nonpolar Fluids\", International Journal of Thermophysics, Vol. 24, No. 1, January 2003. \n\nCp0: Jaeschke, M. and P. Schley,\"Ideal-Gas Thermodynamic Properties for Natural-Gas Applications\", Int. J. Thermophys. v 16, n6 1995");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("n-Octane");
	aliases.push_back("nOctane");
	aliases.push_back("Octane");
	aliases.push_back(std::string("OCTANE"));
	aliases.push_back(std::string("N-OCTANE"));
	REFPROPname.assign("OCTANE");

	BibTeXKeys.EOS = "Span-IJT-2003B";
	BibTeXKeys.CP0 = "Jaeschke-IJT-1995";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-FPE-2004";
	BibTeXKeys.VISCOSITY = "Huber-FPE-2004";
	BibTeXKeys.CONDUCTIVITY = "Huber-FPE-2005";
}
double nOctaneClass::psat(double T)
{
    // Maximum absolute error is 0.067304 % between 216.370001 K and 569.319990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.9979900735822733, 2.010892438708793, -2.3118225183879315, -2.2598069384110318, -2.3011711505356049, 0.17255761105191239 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double nOctaneClass::rhosatL(double T)
{
    // Maximum absolute error is 5.019225 % between 216.370001 K and 569.319990 K
    const double ti[]={0,0.67051537540156969, 1.0506718257585101, 1.3129180160656528, 1.2031682152602974, 17.405127168362831, 1.4687877516582304};
    const double Ni[]={0,17.541685698531406, -224.22847046797935, -507.63069606678488, 603.46463891664382, 1.0428855382259361, 112.266623754866};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double nOctaneClass::rhosatV(double T)
{
	// Max error is  0.491603353952 % between 216.37 and 569.319999 K
    const double ti[]={0, 0.060000000000000005, 0.375, 0.3755, 0.3765, 2.5, 5.833333333333333};
    const double Ni[]={0,-0.42445663008165679, 56405.688240783682, -83966.630380619084, 27556.548441399842, -4.6668441379799557, -4.130640999355661};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(reduce.T/T*summer);
}
void nOctaneClass::ECSParams(double *e_k, double *sigma)
{
	// From Huber 2004
	*e_k = 452.090;
	*sigma = 0.63617;
}
double nOctaneClass::viscosity_Trho(double T, double rho)
{
	// Rainwater-Friend for initial density dependence
	double e_k, sigma;
	this->ECSParams(&e_k,&sigma);
	double Tstar = T/e_k;

	// Dilute gas
	double eta_0 = 0.021357*sqrt(params.molemass*T)/(sigma*sigma*exp(0.335103-0.467898*log(Tstar))); // uPa-s

	// Initial density dependence from Rainwater-Friend
	double b[] = {-19.572881, 219.73999, -1015.3226, 2471.01251, -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158};
	double Bstar = b[0]*pow(Tstar,-0.25*0)+b[1]*pow(Tstar,-0.25*1)+b[2]*pow(Tstar,-0.25*2)+b[3]*pow(Tstar,-0.25*3)+b[4]*pow(Tstar,-0.25*4)+b[5]*pow(Tstar,-0.25*5)+b[6]*pow(Tstar,-0.25*6)+b[7]*pow(Tstar,-2.5)+b[8]*pow(Tstar,-5.5);
	double B = Bstar*0.60221415*sigma*sigma*sigma; // L/mol

	double e[4][3]; // init with zeros
	e[2][1] = -0.103924;
	e[2][2] =   9.92302e-2;	
	e[3][1] =   1.13327e-2;
	e[3][2] =  -3.22455e-2;

	double c[] = {0, 0.606122, 2.0651, 3.07843, -0.879088};

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
double nOctaneClass::conductivity_Trho(double T, double rho)
{
	double lambda_0 = 7.7293e-3-3.7114e-2*(T/crit.T) + 9.7758e-2*pow(T/crit.T,2) - 2.8871e-2*pow(T/crit.T,3); // W/m/K

	double sumresid = 0;
	double B1[] = {0, 0.285553e-1, -0.171398e-1, 0.659971e-2, 0};
	double B2[] = {0, -0.926155e-2, 0, 0.153496e-2, 0};

	for (int i = 1; i<= 4; i++){
		sumresid += (B1[i]+B2[i]*(T/reduce.T))*pow(rho/reduce.rho,i);
	}

	double lambda_r = sumresid; // [W/m/K]

	double lambda_c = this->conductivity_critical(T,rho,0.145713e10); // [W/m/K]

	return lambda_0+lambda_r+lambda_c;
}

nDodecaneClass::nDodecaneClass()
{
	double n[] = {0.0, 1.38031, -2.85352, 0.288897, -0.165993, 0.0923993, 0.000282772, 0.956627, 0.0353076, -0.445008, -0.118911, -0.0366475, 0.0184223};
	double t[] = {0,0.32,1.23,1.5,1.4,0.07,0.8,2.16,1.1,4.1,5.6,14.5,12.0};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 3, 3};
	double d[] = {0, 1, 1, 1, 2, 3, 7, 2, 5, 1, 4, 3, 4};

	//Critical parameters
	crit.rho = 1.33*170.33484; //[kg/m^3]
	crit.p = PressureUnit(1817, UNIT_KPA); //[kPa]
	crit.T = 658.1; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 170.33484;
	params.Ttriple = 263.60;
	params.accentricfactor = 0.574182221240689;
	params.R_u = 8.314472;
	params.ptriple = 0.00062621099412445;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power(n, d, t, c, 1, 12, 13));

	phi0list.push_back(new phi0_lead(0,0));
	phi0list.push_back(new phi0_logtau(23.085-1.0));

	const double a0[] = {0, 37.776, 29.369, 12.461, 7.7733};
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	const double a1[] = {0, 1280/crit.T, 2399/crit.T, 5700/crit.T, 13869/crit.T};
	std::vector<double> a1_v(a1,a1+sizeof(a1)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(a0_v,a1_v,1,4));

	EOSReference.assign("Eric W. Lemmon and Marcia L. Huber, \"Thermodynamic Properties of n-Dodecane\" Energy & Fuels 2004, 18, 960-967");
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("n-Dodecane");
	aliases.push_back("nDodecane");
	aliases.push_back("Dodecane");
	aliases.push_back(std::string("DODECANE"));
	aliases.push_back(std::string("N-DODECANE"));
	REFPROPname.assign("C12");

	BibTeXKeys.EOS = "Lemmon-EF-2004";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.VISCOSITY = "Huber-EF-2004";
	BibTeXKeys.CONDUCTIVITY = "Huber-EF-2004";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-EF-2004";
}
double nDodecaneClass::psat(double T)
{
    // Maximum absolute error is 0.138846 % between 263.600001 K and 658.099990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-8.9866724724738951, 2.5493588271985472, -3.4360168519160204, -2.8106764008676923, -4.2863137232691715, 2.0061541458847785 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double nDodecaneClass::rhosatL(double T)
{
    // Maximum absolute error is 2.000452 % between 263.600001 K and 658.099990 K
    const double ti[]={0,0.18641474693916188, 1.1243571862574657, 12.893082664340085, 1.0855027883831063, 1.0077286384840609};
    const double Ni[]={0,0.79443715883454757, 113.11123644121949, 0.67935350544099971, -172.20630557113711, 59.770329506604732};
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
double nDodecaneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.594552 % between 263.600001 K and 658.099990 K
    const double ti[]={0,0.26387014250856733, 1.159849371473102, 0.98182373858336336, 1.1708598109208148, 3.5424277704241764};
    const double Ni[]={0,-1.3944420975391763, 445.78839795908664, -36.255938474032995, -414.80248500375251, -8.5883385816933817};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
void nDodecaneClass::ECSParams(double *e_k, double *sigma)
{
	// Huber 2004
	*e_k = 522.592;
	*sigma = 0.735639;
}

double nDodecaneClass::viscosity_dilute(double T)
{
	double sigma, e_k;
	// Dilute gas
	this->ECSParams(&e_k, &sigma);
	double Tstar = T/e_k;
	double a[] = {0.382987, -0.561050, 0.313962e-1};
	double S_star = exp(a[0]+a[1]*log(Tstar)+a[2]*log(Tstar)*log(Tstar));
	double eta_0 = 0.021357*sqrt(params.molemass*T)/(sigma*sigma*S_star); //[uPa-s]
	return eta_0/1e6;
}
double nDodecaneClass::viscosity_background(double T, double rho)
{
	double sigma, e_k;
	// Dilute gas
	this->ECSParams(&e_k, &sigma);
	double Tstar = T/e_k;
	double delta = rho/crit.rho;

	// Dilute viscosity
	double eta_0 = this->viscosity_dilute(T); //[Pa-s]

	// Initial-density
	double b[] = {-19.572881, 219.73999, -1015.3226, 2471.0125, -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158};
	double t[] = {0, -0.25, -0.50, -0.75, -1.00, -1.25, -1.50, -2.50, -5.50};
	double sumBstar = 0;	
	for (int i = 0; i<= 8; i++){ sumBstar += b[i]*pow(Tstar,t[i]); }
	double Bstar = sumBstar;
	double N_A = 6.02214129e23;
	double B = N_A*pow(sigma/1e9,3)*Bstar;

	// Residual
	double a21 = -0.471703e-1, a31 = 0.827816e-2, a22 = 0.298429e-1, a32 = -0.134156e-1;
	double c[] = {0, 0.503109, 2.32661, 2.23089};
	double tau = T/crit.T;
	double delta_0 = c[2]+c[3]*sqrt(tau);
	double eta_r = 1000*(a21*pow(delta,2)/tau+a31*pow(delta,3)/tau+a22*pow(delta,2)/pow(tau,2)+a32*pow(delta,3)/pow(tau,2)+c[1]*delta*(1/(delta_0-delta)-1/delta_0));

	double rhobar = rho/params.molemass*1000;
	return eta_0*B*rhobar+eta_r/1e6; //[Pa-s]
}

double nDodecaneClass::viscosity_Trho(double T, double rho)
{
	return this->viscosity_dilute(T) + this->viscosity_background(T,rho);
}

double nDodecaneClass::conductivity_background(double T, double rho)
{
	double sumresid = 0;
	double B1[] = {0, 0.693347e-1, -0.331695e-1, 0.676165e-2};
	double B2[] = {0, -0.280792e-1, 0.173922e-2, 0.309558e-2};
	for (int i = 1; i <= 3; i++)
	{
		sumresid += (B1[i]+B2[i]*T/crit.T)*pow(rho/crit.rho,i);
	}
	double lambda_r = sumresid; // W/m/K
	return lambda_r;
}
double nDodecaneClass::conductivity_Trho(double T, double rho)
{
	double lambda_0 = 0.436343e-2-0.264054e-1*T/crit.T+0.922394e-1*pow(T/crit.T,2)-0.291756e-1*pow(T/crit.T,3); //W/m/K
	double lambda_r = this->conductivity_background(T,rho); //[W/m/K]
	double lambda_c = this->conductivity_critical(T,rho,1/(1.52e-9)); //[W/m/K]

	return (lambda_0 + lambda_r + lambda_c); //[W/m/K]
}

CyclohexaneClass::CyclohexaneClass()
{
	double n[] = {0, 0.0548405, 1.60777, -2.376, -0.51375, 0.18584, -0.90075, -0.5629, 0.2903, -0.3279, -0.03178, 0.8669, -0.1963, -0.1426, 0.004198, 0.1777, -0.04434, -0.03861, 0.0740, 0.02036, 0.0027278};
	double t[] = {0, 1, 0.37, 0.79, 1.075, 0.37, 2.4, 2.5, 0.5, 3, 1.06, 1.6, 0.37, 1.33, 2.5, 0.9, 0.5, 0.73, 0.2, 1.5, 1.5};
	double d[] = {0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3, 3, 2, 2, 3, 2, 3, 2};
	double l[] = {0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	// In REFPROP order
    double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.99, 1.43, 0.97, 1.93, 0.92, 1.27, 0.87, 0.82, 1.4, 3};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.38, 4.2, 1.2, 0.9, 1.2, 2.6, 5.3, 4.4, 4.2, 25};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.65, 0.63, 1.14, 0.09, 0.56, 0.4, 1.01, 0.45, 0.85, 0.86};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.73, 0.75, 0.48, 2.32, 0.2, 1.33, 0.68, 1.11, 1.47, 0.99};

	//Critical parameters
	crit.rho = 3.224*84.15948; //[kg/m^3]
	crit.p = PressureUnit(4082.4, UNIT_KPA); //[kPa]
	crit.T = 553.6; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 84.15948;
	params.Ttriple = 279.47;
	params.accentricfactor = 0.20926;
	params.R_u = 8.3144621;
	params.ptriple = 5.2388581733699020;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n, d, t, l, 1, 10, 21));
	phirlist.push_back(new phir_gaussian(n, d, t, eta, epsilon, beta, gamma, 11, 20, 21));

	phi0list.push_back(new phi0_lead(0.991141, 1.63455));
	phi0list.push_back(new phi0_logtau(4-1));

	double a0[] = {0, 0.83775, 16.036, 24.636, 7.1715};
	double t0[] = {0, 773/crit.T, 941/crit.T, 2185/crit.T, 4495/crit.T};
	phi0list.push_back(new phi0_Planck_Einstein(a0, t0, 1, 4, 5));

	name.assign("CycloHexane");
	aliases.push_back("Cyclohexane");
	aliases.push_back(std::string("CYCLOHEXANE"));
	REFPROPname.assign("CYCLOHEX");

	BibTeXKeys.EOS = "Zhou-PREPRINT-2014";
	BibTeXKeys.ECS_LENNARD_JONES = "Chichester-NIST-2008";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double CyclohexaneClass::psat(double T)
{
	// Max error is  0.233590271792 % between 279.47 and 553.639999 K
    const double t[]={0, 0.10300000000000001, 0.38049999999999995, 0.6666666666666666, 1.3333333333333333, 2.5, 5.833333333333333};
    const double N[]={0, -0.17233802121779696, 1.2533304819707884, -4.1248342131116029, -3.4866469454507585, -0.72673516850438313, -5.2420350342080013};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double CyclohexaneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.198166 % between 279.470000 K and 553.600000 K
    const double t[] = {0, 0.5, 0.6666666666666666, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.6666666666666667};
    const double N[] = {0, 2.0875917314370747, 22.130466870877569, -385.19349831131763, 1207.8352601443141, -1646.6470917840184, 1082.9143691231232, -280.39909465932834};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=7; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double CyclohexaneClass::rhosatV(double T)
{
	// Max error is  0.146651285796 % between 279.47 and 553.639999 K
    const double t[] = {0, 0.39149999999999996, 0.367, 0.3765, 0.38249999999999995, 0.38599999999999995, 4.0};
    const double N[] = {0, -4200293.8657020703, -1003425.4737526291, 8818731.06305697, -25089440.984318022, 21474423.99292187, -4.8624843691761548};
    double summer=0,theta;
    theta=1-T/reduce.T;	    	
	for (int i=1; i<=6; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

//#ifndef DISABLE_CATCH
//class validator_element
//{
//public:
//	std::string in1, in3, in5;
//	double in2, in4, in6;
//	validator_element(std::string in1, double in2, std::string in3, double in4, std::string in5, double in6) : in1(in1), in2(in2), in3(in3), in4(in4), in5(in5), in6(in6) {};
//};
//TEST_CASE(  "Cyclohexane",  "[cyclohexane],[validation]" ) 
//{
//	Fluid *CHEX = get_fluid(get_Fluid_index("Cyclohexane"));
//	double mm = CHEX->params.molemass;
//	validator_element data[] = {
//		validator_element("T",300.0,"D",9.4*mm,"P",24.173705*1e6),
//		validator_element("T",300.0,"D",9.4*mm,"O",115.28612/mm*1000),
//		validator_element("T",300.0,"D",9.4*mm,"C",154.76968/mm*1000),
//		validator_element("T",300.0,"D",9.4*mm,"A",1383.3876),
//		validator_element("T",500.0,"D",6.5*mm,"P",3.9246630*1e6),
//		validator_element("T",500.0,"D",6.5*mm,"O",192.52079/mm*1000),
//		validator_element("T",500.0,"D",6.5*mm,"C",255.57110/mm*1000),
//		validator_element("T",500.0,"D",6.5*mm,"A",434.13058),
//		validator_element("T",500.0,"D",0.7*mm,"P",1.9981172*1e6),
//		validator_element("T",500.0,"D",0.7*mm,"O",191.96468/mm*1000),
//		validator_element("T",500.0,"D",0.7*mm,"C",235.52304/mm*1000),
//		validator_element("T",500.0,"D",0.7*mm,"A",155.34798),
//		validator_element("T",600.0,"D",3.5*mm,"P",6.8225506*1e6),  
//		validator_element("T",600.0,"D",3.5*mm,"O",232.79249/mm*1000),
//		validator_element("T",600.0,"D",3.5*mm,"C",388.55212/mm*1000),
//		validator_element("T",600.0,"D",3.5*mm,"A",150.53314),
//		validator_element("T",553.6,"D",3.3*mm,"P",4.0805433*1e6),
//		validator_element("T",553.6,"D",3.3*mm,"O",224.19580/mm*1000), 
//		validator_element("T",553.6,"D",3.3*mm,"C",199224.62/mm*1000), 
//		validator_element("T",553.6,"D",3.3*mm,"A",87.913862)
//	};
//
//	//Now actually construct the vector
//	std::vector<validator_element> elements(data, data + sizeof(data) / sizeof(validator_element));		
//	
//	for (std::vector<validator_element>::iterator it = elements.begin(); it != elements.end(); it++)
//	{
//		validator_element &el = *it;
//		double eos = PropsSI( el.in5.c_str(),  el.in1.c_str(), el.in2,  el.in3.c_str(), el.in4, "Cyclohexane");
//		double valid = el.in6;
//        CAPTURE(eos);
//        CAPTURE(valid);
//        CAPTURE(el.in1);
//        CAPTURE(el.in3);
//        CAPTURE(el.in5);
//		CHECK(fabs(valid/eos-1) < 1e-8);
//	}
//}
//
//#endif


R152AClass::R152AClass()
{
	double n[] = {0.0, 0.95702326, -2.3707196, 0.18748463, 0.063800843, 0.00016625977, 0.082208165, 0.57243518, 0.0039476701, -0.23848654, -0.080711618, -0.073103558, -0.015538724};

	double a0[] = {0,-20.78887,-0.6539092,0.03342831};
	double n0[] = {0,-1.0/4.0,-2,-4};
	//Critical parameters
	crit.rho = 368.0; //[kg/m^3]
	crit.p = PressureUnit(4520, UNIT_KPA); //[kPa]
	crit.T = 386.411; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 66.051;
	params.Ttriple = 154.56;
	params.accentricfactor = 0.27521711453209896;
	params.R_u = 8.314510;
	params.ptriple = 0.064089867936946279;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d_polar_SpanShort,t_polar_SpanShort,c_polar_SpanShort,1,12,13));

	phi0list.push_back(new phi0_lead(10.87227,6.839515));
	phi0list.push_back(new phi0_logtau(-1.0));
	phi0list.push_back(new phi0_power(a0,n0,1,3,4));

	static char EOSstr [] = "Span, R. and W. Wagner, \"Equations of State for Technical Applications. III. Results for Polar Fluids\", International Journal of Thermophysics, Vol. 24, No. 1, January 2003. \n\nCp0: R. Tillner-Roth \"A Fundamental Equation of State for 1,1-Difluoroethane (HFC-152a)\" International Journal of Thermophysics. vol. 16. No. 1 1995\n\nNote: REFPROP 9.0 uses the MWBR formulation";
	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode. ");

	name.assign("R152A");
	aliases.push_back("R152a");
	REFPROPname.assign("R152A");

	BibTeXKeys.EOS = "Span-IJT-2003C";
	BibTeXKeys.CP0 = "TillnerRoth-IJT-1995";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.VISCOSITY = "Krauss-IJT-1996";
	BibTeXKeys.CONDUCTIVITY = "Krauss-IJT-1996";
}
double R152AClass::psat(double T)
{
    // Maximum absolute error is 0.040697 % between 154.560001 K and 386.410990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3,9};
    const double Ni[]={0,-7.5183210509323537, 2.2667633168445649, -2.9745826359338392, 2.3313114019766128, -9.5143559559629267, 12.699734081556912, -9.5537989669345524 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=7;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double R152AClass::rhosatL(double T)
{
    // Maximum absolute error is 0.396356 % between 154.560001 K and 386.410990 K
    const double ti[]={0,0.36433465612787069, 1.3395245009162089, 1.1775553689273002, 1.2736011086820553, 1.39293878426836, 1.5006057434930407};
    const double Ni[]={0,2.2629063243095731, -8137.3963531047102, -739.48710570556011, 4453.3069484532643, 4994.763056241005, -572.05896753939112};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double R152AClass::rhosatV(double T)
{
	// Max error is  0.180395311787 % between 154.56 and 386.410999 K

    const double ti[]={0, 0.098, 0.365, 0.378, 0.38, 0.38749999999999996, 3.6666666666666665};
    const double Ni[]={0, -2.1665803273891777, 33129.223069368025, -600164.85585399868, 661053.48772902705, -94021.00606766867, -4.6394042521818628};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
void R152AClass::ECSParams(double *e_k, double *sigma)
{
	// Krauss 1996
	*e_k = 354.84; *sigma = 0.46115;
}
double R152AClass::viscosity_Trho(double T, double rho)
{
	// Krauss 1996
	double sigma, e_k;
	
	// Dilute gas
	this->ECSParams(&e_k, &sigma);
	double Tstar = T/e_k;
	double a[] = {0.4425728, -0.5138403,0.1547566,-0.02821844,0.001578286};
	double logTstar = log(Tstar);
	double S_star = exp(a[0]+a[1]*logTstar+a[2]*pow(logTstar,2)+a[3]*pow(logTstar,3)+a[4]*pow(logTstar,4));
	double eta_0 = 0.2169614*sqrt(T)/(sigma*sigma*S_star); //[uPa-s]

	// Initial-density
	double E[] = {0, -0.0737927,0.517924,-0.308875,0.108049,-0.408387,2.91733};
	double eta_r = 0;
	for (int i = 1; i<= 4; i++){ eta_r += E[i]*pow(rho/reduce.rho,i); }
	double Hc = 51.12;
	eta_r = Hc*(eta_r+E[5]/(rho/reduce.rho-E[6])+E[5]/E[6]); // [uPa-s]

	return (eta_0+eta_r)/1e6;
}
double R152AClass::conductivity_Trho(double T, double rho)
{
	// Krauss 1996
	double L[] = {0,9.18090,11.8577,-5.44730,1.71379};
	double lambda_0 = -14.9420+0.0973283*T; // [mW/m/K]

	// Residual part
	double summer = 0;
	for (int i = 1; i<= 4; i++){ summer += L[i]*pow(rho/reduce.rho,i); }
	double lambda_r = 1.155*summer; // [mW/m/K]

	double lambda_c = this->conductivity_critical(T,rho,1/(4.37e-10))*1e6; //[mW/m/K]

	return (lambda_0+lambda_r+lambda_c)/1e6; // [kW/m/K]
}

R123Class::R123Class()
{
	double n[] = {0.0, 1.116973, -3.074593, 0.51063873, 0.094478812, 0.00029532752, 0.66974438, 0.96438575, -0.014865424, -0.49221959, -0.022831038, -0.1407486, -0.025117301};

	double a0[] = {0,2.046009,22.231991/456.82,-11.658491/456.82/456.82,2.691665/456.82/456.82/456.82}; // Fit is in terms of Tr = T/Tc
	double n0[] = {0,0,1,2,3};
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));
	
	//Critical parameters
	crit.rho = 553.0; //[kg/m^3]
	crit.p = PressureUnit(3672, UNIT_KPA); //[kPa]
	crit.T = 456.82; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 152.931;
	params.Ttriple = 166;
	params.accentricfactor = 0.28192249703635186;
	params.R_u = 8.314510;
	params.ptriple = 0.0041534426564210506;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d_polar_SpanShort,t_polar_SpanShort,c_polar_SpanShort,1,12,13));

	phi0list.push_back(new phi0_lead(0,0));
	phi0list.push_back(new phi0_logtau(-1.0));
	phi0list.push_back(new phi0_cp0_constant(a0[1],crit.T,298));
	phi0list.push_back(new phi0_cp0_poly(a0_v,n0_v,crit.T,298,2,4));

	static char EOSstr [] = "Span, R. and W. Wagner, \"Equations of State for Technical Applications. III. Results for Polar Fluids\", International Journal of Thermophysics, Vol. 24, No. 1, January 2003. \n\nCp0: Ben A. Younglove \"An International Standard Equation of State for the Thermodynamic Properties of Refrigerant 123 (2,2-Dichloro-1,1,1-Trifluoroethane)\" J. Phys. Chem Ref. Data Vol 23, no 5, 1994\n\nNote: REFPROP 9.0 uses the less accurate MWBR formulation";
	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode. ");

	name.assign("R123");
	REFPROPname.assign("R123");

	BibTeXKeys.EOS = "Span-IJT-2003C";
	BibTeXKeys.CP0 = "Younglove-JPCRD-1994";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.CONDUCTIVITY = "Laesecke-IJR-1996";
	BibTeXKeys.VISCOSITY = "Tanaka-IJT-1996";
	BibTeXKeys.ECS_LENNARD_JONES = "Tanaka-IJT-1996";
}
double R123Class::psat(double T)
{
    // Maximum absolute error is 0.014209 % between 166.000001 K and 456.830990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3,9};
    const double Ni[]={0,-7.4849207348463391, 2.1531543385109324, -2.82304801818721, 1.4302672539730645, -9.5217923675635188, 14.152141293153656, -11.358606061743817 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=7;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double R123Class::rhosatL(double T)
{
    // Maximum absolute error is 0.654849 % between 166.000001 K and 456.830990 K
    const double ti[]={0,0.30743229169537528, 0.87399412478482197, 1.9671212656736285, 3.0665054723952703, 3.2018449236245567, 12.488775734339216};
    const double Ni[]={0,1.4776907822440624, -0.29754626202706141, 0.27387548048866328, -0.049678358639705204, -0.084295473834870918, 0.94874134293977097};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double R123Class::rhosatV(double T)
{
    // Max error is  0.0199273635887 % between 166.0 and 456.8199 K
    const double ti[]={0, 0.11800000000000001, 0.12100000000000001, 0.131, 0.066, 0.351, 0.361, 0.3695, 0.381, 2.5, 4.5};
    const double Ni[]={0, 429808.91452099587, -565461.07822199445, 139089.70648760325, -1555.8578707113704, -670815.34308265569, 1956365.6938359118, -1678623.2858405211, 391186.25477907102, -1.9994703100835092, -4.2275468110525605};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i = 1; i <= 10; i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
void R123Class::ECSParams(double *e_k, double *sigma)
{
	// From Tanaka 1996
	*e_k = 340;
	*sigma = 0.56;
}
double R123Class::viscosity_Trho(double T, double rho)
{
	// From Tanaka 1996
	double eta_0 = -2.273638 + 5.099859e-2*T - 2.402786e-5*T*T;
	
	double eta_1 = -2.226486e-2 + 5.550623e-5*T;

	double rho0 = 1.828263e3, a0 = -3.222951e5, a1 = -1.009812e-1, a2 = 6.161902e-5, a3 = -8.84048e-8;

	double eta_r = a0/(rho-rho0)+a0/rho0+a1*rho+a2*rho*rho+a3*rho*rho*rho;

	return (eta_0 + eta_1*rho + eta_r)/1e6;

}
double R123Class::conductivity_Trho(double T, double rho)
{
	double tau = 456.831/T, delta = rho/550;

	double lambda_0 = 5.695e-5*T-0.00778; // [W/m/K]

	double d[] = {0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4};
	double t[] = {0, 1.5, 2, 6, 0, 0.5, 1.5, 0, 0.5, 1.5, 0, 0.5, 1.5};
	double a[] = {0, 0.642894e-1, -0.530474e-1, 0.453522e-4, -0.139928, 0.166540, -0.162656e-1, 0.136819, -0.183291, 0.357146e-1, -0.231210e-1, 0.341945e-1, -0.757341e-2};
	double a13 = 0.486742e-2, a14 = -100, a15 = -7.08535;

	double sumresid = 0;

	for (int i = 1; i <= 12; i++)
	{
		sumresid += a[i]*pow(delta,d[i])*pow(tau,t[i]);
	}
	double lambda_r = sumresid; // W/m/K

	double lambda_c = a13*exp(a14*pow(tau-1,4)+a15*pow(delta-1,2));

	return (lambda_0 + lambda_r + lambda_c)/1000; // kW/m/K
}



R11Class::R11Class()
{
	const double n[] = {0.0, 1.0656383, -3.2495206, 0.87823894, 0.087611569, 0.00029950049, 0.42896949, 0.70828452, -0.017391823, -0.37626522, 0.011605284, -0.089550567, -0.030063991};

	//Critical parameters
	crit.rho = 565; //[kg/m^3]
	crit.p = PressureUnit(4394, UNIT_KPA); //[kPa]
	crit.T = 471.06; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 137.368;
	params.Ttriple = 162.68;
	params.accentricfactor = 0.18875064825280830;
	params.R_u = 8.314510;
	params.ptriple = 0.0066915477602794626;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d_polar_SpanShort,t_polar_SpanShort,c_polar_SpanShort,1,12,13));

	phi0list.push_back(new phi0_lead(0,0));
	phi0list.push_back(new phi0_logtau(-1.0));
	phi0list.push_back(new phi0_cp0_constant(4.0+0.0469706/(params.R_u),471.11,298));
	phi0list.push_back(new phi0_cp0_poly(0.0018532/(params.R_u),1,471.11,298));

	double a[] = {0,1,2,1,2,1,2};
	double A[] = {0,1085,847,535.2,398,349.5,241};
	for (int i = 1; i<=6; i++) {A[i] *= 1.43878/471.11; };
	std::vector<double> A_v(A,A+sizeof(A)/sizeof(double));
	std::vector<double> a_v(a,a+sizeof(a)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(a_v,A_v,1,6));

	static char EOSstr [] = "Span, R. and W. Wagner, \"Equations of State for Technical Applications. III. Results for Polar Fluids\", International Journal of Thermophysics, Vol. 24, No. 1, January 2003. \n\nCp0: Jacobsen et al. \"A Fundamental Equation for Trichlorofluoromethane (R-11) Fluid Phase Equilibria, 1992\"";
	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode. ");

	name.assign("R11");
	REFPROPname.assign("R11");

	BibTeXKeys.EOS = "Span-IJT-2003C";
	BibTeXKeys.CP0 = "Jacobsen-FPE-1992";
	BibTeXKeys.ECS_LENNARD_JONES = "McLinden-IJR-2000";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R11Class::psat(double T)
{
	// Max error is  0.101004895012 % between 162.68 and 471.059999 K
    const double ti[]={0, 0.07, 0.3755, 1.0, 2.0, 3.0, 6.166666666666667};
    const double Ni[]={0, 0.015192240512549519, -0.073596807283158483, -6.4214975426915206, 1.4226706834098106, -3.0447362413765453, -1.7705707854921364};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1; i<=6; i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double R11Class::rhosatL(double T)
{
	// Max error is  2.4331028261 % between 162.68 and 471.06 K
    const double ti[]={0, 0.07, 0.07300000000000001, 0.07400000000000001, 0.07600000000000001, 0.082, 0.085, 0.093, 0.12000000000000001};
    const double Ni[]={0, 248395134649.64792, -4147174258833.4897, 6205839573055.1846, -2666372960037.8979, 623915721315.88281, -277014946706.59399, 12442556940.240545, -30820380.285890266};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1; i<=8; i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*(summer+1);
}
double R11Class::rhosatV(double T)
{
    // Max error is  0.716817629755 % between 162.68 and 471.0599 K
	
    const double ti[]={0, 0.062, 0.066, 0.07300000000000001, 0.12100000000000001, 0.357, 0.3595, 0.3615, 0.3625, 0.3695, 5.833333333333333};
    const double Ni[]={0, -39932620.310829446, 73762127.862868786, -36143788.042476922, 2832253.8625494828, -297296943057.29657, 1497747443415.0725, -3063566982511.5029, 1892454073153.1245, -29338108979.462494, -3.5661612181452496};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=10;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
