


#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Span_Short.h"

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

static char EOSstr [] = "Span, R. and W. Wagner, \"Equations of State for Technical Applications. II. Results for Nonpolar Fluids\", International Journal of Thermophysics, Vol. 24, No. 1, January 2003. \n\nCp0: Jaeschke, M. and P. Schley,\"Ideal-Gas Thermodynamic Properties for Natural-Gas Applications\", Int. J. Thermophys. v 16, n6 1995";

nPentaneClass::nPentaneClass()
{
const double n[] = {0.0, 1.0968643, -2.9988888, 0.99516887, -0.16170709, 0.1133446, 0.00026760595, 0.40979882, -0.040876423, -0.38169482, -0.10931957, -0.032073223, 0.016877016};
std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));

//Critical parameters
crit.rho = 232; //[kg/m^3]
crit.p = 3370; //[kPa]
crit.T = 469.7; //[K]
crit.v = 1/crit.rho; 

// Other fluid parameters
params.molemass = 72.15;
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

EOSReference.assign(EOSstr);
TransportReference.assign("Using ECS in fully predictive mode");

name.assign("n-Pentane");
aliases.push_back(std::string("nPentane"));
aliases.push_back(std::string("Pentane"));
REFPROPname.assign("PENTANE");
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
    return reduce.p*exp(reduce.T/T*summer);
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
std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));

//Critical parameters
crit.rho = 233.18; //[kg/m^3]
crit.p = 3034; //[kPa]
crit.T = 507.82; //[K]
crit.v = 1/crit.rho; 

// Other fluid parameters
params.molemass = 86.177;
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

EOSReference.assign(EOSstr);
TransportReference.assign("Using ECS in fully predictive mode");

name.assign("n-Hexane");
aliases.push_back("nHexane");
aliases.push_back("Hexane");
REFPROPname.assign("HEXANE");

}
double nHexaneClass::psat(double T)
{
    // Maximum absolute error is 0.309290 % between 177.830001 K and 507.819990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.4068558733737699, 1.2751526610932369, -0.41863507831564806, -4.3671933839461419, 1.5308874611739092, -1.826117481826562 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double nHexaneClass::rhosatL(double T)
{
    // Maximum absolute error is 6.631534 % between 177.830001 K and 507.819990 K
    const double ti[]={0,0.33712689148545583, 0.81990784826075724, 1.4191327364931061, 2.6093422842418756, 3.8331703471313796};
    const double Ni[]={0,1.8072779248452011, -1.1053012310408719, 0.95401021898294835, -0.57675230346206185, 0.34981803464773442};
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
double nHexaneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.933021 % between 177.830001 K and 507.819990 K
    const double ti[]={0,0.42284371311323915, 0.60042568385779993, 2.1951394734352969, 3.2488641907097935, 3.2172188773573409};
    const double Ni[]={0,-2.2822943614933475, -2.3485033954974881, -4.4452043943305544, -174.10979741588503, 172.56055164398012};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}


nHeptaneClass::nHeptaneClass()
{
const double n[] = {0.0, 1.0543748, -2.6500682, 0.81730048, -0.30451391, 0.12253869, 0.00027266473, 0.49865826, -0.00071432815, -0.54236896, -0.13801822, -0.0061595287, 0.0004860251};
std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));

//Critical parameters
crit.rho = 232; //[kg/m^3]
crit.p = 2736; //[kPa]
crit.T = 540.13; //[K]
crit.v = 1/crit.rho; 

// Other fluid parameters
params.molemass = 100.204;
params.Ttriple = 182.55;
params.accentricfactor = 0.349;
params.R_u = 8.314510;
params.ptriple = 0.0019889200572;

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

EOSReference.assign(EOSstr);
TransportReference.assign("Using ECS in fully predictive mode");

name.assign("n-Heptane");
aliases.push_back("nHeptane");
aliases.push_back("Heptane");
REFPROPname.assign("HEPTANE");
}
double nHeptaneClass::psat(double T)
{
    // Maximum absolute error is 0.124229 % between 182.550001 K and 540.129990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.737566099753602, 1.6913411007500037, -1.4056091256177305, -3.2093008874967563, -0.64037236634384165, -0.80501812286028029 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double nHeptaneClass::rhosatL(double T)
{
    // Maximum absolute error is 6.968569 % between 182.550001 K and 540.129990 K
    const double ti[]={0,0.43144258579014039, 0.79871865359256855, 0.86763070714645008, -0.0039292068918103896};
    const double Ni[]={0,2.5212528240602059, -6.4957880126225032, 5.2404429673697379, 0.10230051516352442};
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
double nHeptaneClass::rhosatV(double T)
{
    // Maximum absolute error is 11.753383 % between 182.550001 K and 540.129990 K
    const double ti[]={0,0.054000000000000006, 0.357, 0.39149999999999996, 0.39549999999999996, 2.6666666666666665, 7.166666666666667};
    const double Ni[]={0,-0.38409500479534564, -83.909566768017612, 1131.3646812581319, -1053.3437282006555, -22.578893901804172, -83.582999448990478};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

nOctaneClass::nOctaneClass()
{
	const double n[] = {0.0, 1.0722545, -2.4632951, 0.65386674, -0.36324974, 0.1271327, 0.00030713573, 0.52656857, 0.019362863, -0.58939427, -0.14069964, -0.0078966331, 0.0033036598};
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
	std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
	std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));

	//Critical parameters
	crit.rho = 234.9; //[kg/m^3]
	crit.p = 2497; //[kPa]
	crit.T = 569.32; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 114.231;
	params.Ttriple = 216.37;
	params.accentricfactor = 0.395;
	params.R_u = 8.314510;
	params.ptriple = 0.000175490394469;

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

	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("n-Octane");
	aliases.push_back("nOctane");
	aliases.push_back("Octane");
	REFPROPname.assign("OCTANE");
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
    return reduce.p*exp(reduce.T/T*summer);
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
	double rho = reduce.rho*exp(summer);
    return reduce.rho*exp(summer);
}
double nOctaneClass::rhosatV(double T)
{
    // Maximum absolute error is 8.0572339241 % between 216.370001 K and 569.319990 K
    const double ti[]={0,0.12200000000000001, 0.3535, 0.373, 2.0, 5.0, 10.5};
    const double Ni[]={0,-0.5737417164715789, 37.14765536149423, -41.588310716565836, -14.168564457278846, -48.161300312807803, -92.480562862885606};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
	double rho = reduce.rho*exp(summer);
    return reduce.rho*exp(summer);
}

CyclohexaneClass::CyclohexaneClass()
{
	const double n[] = {0.0, 1.0232354, -2.9204964, 1.0736630, -0.19573985, 0.12228111, 0.00028943321, 0.27231767, -0.04483332, -0.38253334, -0.089835333, -0.024874965, 0.010836132};
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
	std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
	std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));

	//Critical parameters
	crit.rho = 273.02; //[kg/m^3]
	crit.p = 4075; //[kPa]
	crit.T = 553.64; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 84.1608;
	params.Ttriple = 279.47;
	params.accentricfactor = 0.20926;
	params.R_u = 8.314510;
	params.ptriple = 5.25384210955;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n_v,d_v,t_v,l_v,1,12));

	double T0 = 353.885834409,
	R_ = 8.314510/params.molemass,
	rho0 = 719.521755268,
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

	const double a0[] = {0,4*params.R_u,8.97575*params.R_u,438.27, 5.25156*params.R_u, 198.018};
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	const double a1[] = {0,0,25.1423*params.R_u, 1905.02, 16.1388*params.R_u, 893.765};
	std::vector<double> a1_v(a1,a1+sizeof(a1)/sizeof(double));

	phi0list.push_back(new phi0_cp0_AlyLee(a0_v,crit.T,T0,params.R_u));
	phi0list.push_back(new phi0_cp0_AlyLee(a1_v,crit.T,T0,params.R_u));

	EOSReference.assign(EOSstr);
    // Also see Penoncello, 1995
	TransportReference.assign("Using ECS in fully predictive mode.  ECS parameters from ");

	name.assign("CycloHexane");
	aliases.push_back("Cyclohexane");
	REFPROPname.assign("CYCLOHEX");
}
double CyclohexaneClass::psat(double T)
{
    // Maximum absolute error is 0.027225 % between 279.470001 K and 553.639990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.0355383856942639, 1.8142034452759799, -2.2917261023281044, 2.0538265028886675, -11.310993893929071, 13.09497402169962 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double CyclohexaneClass::rhosatL(double T)
{
    // Maximum absolute error is 1.990271 % between 279.470001 K and 553.639990 K
    const double ti[]={0,1.0435686062881251, 1.3253276706987631, 1.1698308973705493, 1.1318897634104805, 1.0814858893687547};
    const double Ni[]={0,24183.906537608644, 1629.0615494474951, -46379.189895301017, 91594.514971782599, -71026.798235435912};
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
double CyclohexaneClass::rhosatV(double T)
{
    // Maximum absolute error is 1.791652 % between 279.470001 K and 553.639990 K
    const double ti[]={0,0.53973386168599435, 2.3863576122670151, 2.1668615412073802, 2.3159135580489689, 2.3580298904980976};
    const double Ni[]={0,-4.6730291488786895, -13785.837810534413, 693.33747953168393, -12275.170581637027, 25362.530017757046};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
