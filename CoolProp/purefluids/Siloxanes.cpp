
#include "CoolProp.h"
#include "Siloxanes.h"
#include <vector>
#include "CPExceptions.h"

static const double d[] =
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

static const double t[] =
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

static const double c[] =
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

//MDM
OctamethyltrisiloxaneClass::OctamethyltrisiloxaneClass()
{
    const double n[]={0.0,1.19735372,-2.40380622,0.3256564,-0.19971259,0.11206277,0.00015893999,0.51234323,-0.020660361,-0.38978114,-0.1186931,-0.037203537,0.018359984};
    const double u0[]={0.0,275.1,612.9,1829.6,413.0,802.6};

	std::vector<double> d_v(d,d+sizeof(d)/sizeof(double));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(c,c+sizeof(c)/sizeof(double));
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));

    // Critical parameters
    crit.rho = 256.739908797209;
    crit.p = 1415;
    crit.T = 564.09;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 236.53146;
    params.Ttriple = 187.2;
    params.accentricfactor = 0.5297;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = params.Ttriple;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

	double T0 = 425.676967,
		   R_ = 8.314472/params.molemass,
		   rho0 = 675.2187471,
		   m,
		   c,
		   H0 = 156.46565194801141, /// kJ/kmol
		   S0 = 0.20329670191815147, /// kJ/kmol/K
		   tau0=crit.T/T0,
		   delta0=rho0/crit.rho;
	
	// log(delta)+c+m*tau
	
	/// c is the constant term
	c=-S0/R_-1+log(tau0/delta0);/*<< from the leading term*/

	/// m multiplies the tau term in the leading term (slope)
	m=H0/(R_*crit.T); /*<< from the leading term */

	phi_BC * phi0_lead_ = new phi0_lead(c,m);
	phi0list.push_back(phi0_lead_);

	phi_BC * phi0_logtau_ = new phi0_logtau(-1.0);
	phi0list.push_back(phi0_logtau_);

	phi_BC * phi0_cp0_AlyLee_ = new phi0_cp0_AlyLee(u0_v,crit.T,T0,params.R_u);
	phi0list.push_back(phi0_cp0_AlyLee_);

    EOSReference.assign("P. Colonna, N.R. Nannan, A. Guardone, \"Multiparameter equations of state for siloxanes: [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,...,3, and [O-Si-(CH3)2]6,\", Fluid Phase Equilibria 263 (2008) 115–130.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("MDM");
    aliases.push_back(std::string("Octamethyltrisiloxane")); 
    REFPROPname.assign("MDM");
}
double OctamethyltrisiloxaneClass::rhosatL(double T) 
{
	double A = reduce.p*params.molemass/params.R_u/crit.T;
	double B = reduce.p/(reduce.rho*params.R_u/params.molemass*reduce.T);
	double C = reduce.T;
	double D = 2.0/7.0;
    return A/pow(B,1+pow(1-T/C,D));
}
double OctamethyltrisiloxaneClass::rhosatV(double T) 
{
    double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 96.205848 %
	RHS = +0.740718*pow(theta,0.3)
		  -26.952900*pow(theta,0.9)
		  +37.116438*pow(theta,1.5)
		  -83.326173*pow(theta,3.2);
	rho = exp(RHS)*rhoc;
	return rho;
}
double OctamethyltrisiloxaneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -8.6693+2.2965*pow(THETA,1.5)-4.4658*pow(THETA,2.5)-8.4529*pow(THETA,5);
    return exp(crit.T/T*RHS)*crit.p;
}

//MD2M
DecamethyltetrasiloxaneClass::DecamethyltetrasiloxaneClass()
{
    const double n[]={0.0,1.33840331,-2.62939393,0.4398383,-0.53496715,0.1818844,0.00040774609,1.13444506,0.05774631,-0.5917498,-0.11020225,-0.034942635,0.007646298};
    const double u0[]={0.0,331.9,777.1,1813.8,521.4,795.1};

	std::vector<double> d_v(d,d+sizeof(d)/sizeof(double));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(c,c+sizeof(c)/sizeof(double));
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));

    // Critical parameters
    crit.rho=284.171639662027;
	crit.p=1227;
	crit.T=599.4;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
	params.molemass=310.685;
	params.Ttriple=205.2;
	params.accentricfactor=0.668;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = params.Ttriple;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

	double T0 = 467.506025761878,
		   rho0 = 666.100681528704,
		   R_ = params.R_u/params.molemass,
		   m,
		   c,
		   H0 = 133.21777124359349, /// kJ/kmol
		   S0 = 0.16466470350893656, /// kJ/kmol/K
		   tau0=crit.T/T0,
		   delta0=rho0/crit.rho;
	
	// log(delta)+c+m*tau
	
	/// c is the constant term
	c=-S0/R_-1+log(tau0/delta0);/*<< from the leading term*/

	/// m multiplies the tau term in the leading term (slope)
	m=H0/(R_*crit.T); /*<< from the leading term */

	phi_BC * phi0_lead_ = new phi0_lead(c,m);
	phi0list.push_back(phi0_lead_);

	phi_BC * phi0_logtau_ = new phi0_logtau(-1.0);
	phi0list.push_back(phi0_logtau_);

	phi_BC * phi0_cp0_AlyLee_ = new phi0_cp0_AlyLee(u0_v,crit.T,T0,params.R_u);
	phi0list.push_back(phi0_cp0_AlyLee_);

    EOSReference.assign("P. Colonna, N.R. Nannan, A. Guardone, \"Multiparameter equations of state for siloxanes: [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,...,3, and [O-Si-(CH3)2]6,\", Fluid Phase Equilibria 263 (2008) 115–130.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("MD2M");
    aliases.push_back(std::string("Decamethyltetrasiloxane")); 
    REFPROPname.assign("MD2M");
}
double DecamethyltetrasiloxaneClass::rhosatL(double T) 
{
	double A = reduce.p*params.molemass/params.R_u/crit.T;
	double B = reduce.p/(reduce.rho*params.R_u/params.molemass*reduce.T);
	double C = reduce.T;
	double D = 2.0/7.0;
    return A/pow(B,1+pow(1-T/C,D));
}
double DecamethyltetrasiloxaneClass::rhosatV(double T) 
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 81.362003 %
	RHS = +0.564413*pow(theta,0.3)
		  -26.829650*pow(theta,0.9)
		  +35.702268*pow(theta,1.5)
		  -85.660462*pow(theta,3.2);
	rho = exp(RHS)*rhoc;
	return rho;
}
double DecamethyltetrasiloxaneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -10.072+4.1849*pow(THETA,2.2965)+-6.9586*pow(THETA,-4.4658)+-4.8277*pow(THETA,-8.4529);
    return exp(crit.T/T*RHS)*crit.p;
}

//MD3M
DodecamethylpentasiloxaneClass::DodecamethylpentasiloxaneClass()
{
    const double n[]={0.0,1.20540386,-2.42914797,0.69016432,-0.69268041,0.18506046,0.00031161436,0.99862519,0.074229034,-0.80259136,-0.20865337,-0.036461791,0.019174051};
    const double u0[]={0.0,463.2,957.2,2117.1,738.3,908.5};

	std::vector<double> d_v(d,d+sizeof(d)/sizeof(double));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(c,c+sizeof(c)/sizeof(double));
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));

    // Critical parameters
    crit.rho=263.921879135305;
	crit.p=945;
	crit.T=628.36;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
	params.molemass=384.839;
	params.Ttriple=192;
	params.accentricfactor=0.722;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = params.Ttriple;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

	double T0 = 503.022623775204,
		   rho0 = 653.620754664036,
		   R_ = params.R_u/params.molemass,
		   m,
		   c,
		   H0 = 113.59777719014585, /// kJ/kmol
		   S0 = 0.13288817176242626, /// kJ/kmol/K
		   tau0=crit.T/T0,
		   delta0=rho0/crit.rho;
	
	// log(delta)+c+m*tau
	
	/// c is the constant term
	c=-S0/R_-1+log(tau0/delta0);/*<< from the leading term*/

	/// m multiplies the tau term in the leading term (slope)
	m=H0/(R_*crit.T); /*<< from the leading term */

	phi_BC * phi0_lead_ = new phi0_lead(c,m);
	phi0list.push_back(phi0_lead_);

	phi_BC * phi0_logtau_ = new phi0_logtau(-1.0);
	phi0list.push_back(phi0_logtau_);

	phi_BC * phi0_cp0_AlyLee_ = new phi0_cp0_AlyLee(u0_v,crit.T,T0,params.R_u);
	phi0list.push_back(phi0_cp0_AlyLee_);

    EOSReference.assign("P. Colonna, N.R. Nannan, A. Guardone, \"Multiparameter equations of state for siloxanes: [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,...,3, and [O-Si-(CH3)2]6,\", Fluid Phase Equilibria 263 (2008) 115–130.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("MD3M");
    aliases.push_back(std::string("Dodecamethylpentasiloxane")); 
    REFPROPname.assign("MD3M");
}
double DodecamethylpentasiloxaneClass::rhosatL(double T) 
{
	double A = reduce.p*params.molemass/params.R_u/crit.T;
	double B = reduce.p/(reduce.rho*params.R_u/params.molemass*reduce.T);
	double C = reduce.T;
	double D = 2.0/7.0;
    return A/pow(B,1+pow(1-T/C,D));
}
double DodecamethylpentasiloxaneClass::rhosatV(double T) 
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 211.128665 %
	RHS = +2.193063*pow(theta,0.3)
		  -39.614003*pow(theta,0.9)
		  +58.043624*pow(theta,1.5)
		  -112.804111*pow(theta,3.2);
	rho = exp(RHS)*rhoc;
	return rho;
}
double DodecamethylpentasiloxaneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -9.531+2.7012*pow(THETA,4.1849)+-6.9841*pow(THETA,-6.9586)+-6.5038*pow(THETA,-4.8277);
    return exp(crit.T/T*RHS)*crit.p;
}

//D6
DodecamethylcyclohexasiloxaneClass::DodecamethylcyclohexasiloxaneClass()
{
    const double n[]={0.0,1.69156186,-3.37962568,0.38609039,0.064598995,0.10589012,0.000045456825,0.74169279,-0.088102648,-0.17373336,-0.10951368,-0.062695695,0.037459986};
    const double u0[]={0.0,468.7,981.2,1792.1,686.7,786.8};

	std::vector<double> d_v(d,d+sizeof(d)/sizeof(double));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(c,c+sizeof(c)/sizeof(double));
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));

    // Critical parameters
	crit.rho=279.095729841367;
	crit.p=961;
	crit.T=645.78;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
	params.molemass=444.924;
	params.Ttriple=270.2;
	params.accentricfactor=0.736;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = params.Ttriple;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

	double T0 = 518.109977174843,
		   rho0 = 704.948896461078,
		   R_ = params.R_u/params.molemass,
		   m,
		   c,
		   H0 = 102.00241183448354, /// kJ/kmol
		   S0 = 0.11646477379595818, /// kJ/kmol/K
		   tau0=crit.T/T0,
		   delta0=rho0/crit.rho;
	
	// log(delta)+c+m*tau
	
	/// c is the constant term
	c=-S0/R_-1+log(tau0/delta0);/*<< from the leading term*/

	/// m multiplies the tau term in the leading term (slope)
	m=H0/(R_*crit.T); /*<< from the leading term */

	phi_BC * phi0_lead_ = new phi0_lead(c,m);
	phi0list.push_back(phi0_lead_);

	phi_BC * phi0_logtau_ = new phi0_logtau(-1.0);
	phi0list.push_back(phi0_logtau_);

	phi_BC * phi0_cp0_AlyLee_ = new phi0_cp0_AlyLee(u0_v,crit.T,T0,params.R_u);
	phi0list.push_back(phi0_cp0_AlyLee_);

    EOSReference.assign("P. Colonna, N.R. Nannan, A. Guardone, \"Multiparameter equations of state for siloxanes: [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,...,3, and [O-Si-(CH3)2]6,\", Fluid Phase Equilibria 263 (2008) 115–130.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("D6");
    aliases.push_back(std::string("Dodecamethylcyclohexasiloxane")); 
    REFPROPname.assign("D6");
}
double DodecamethylcyclohexasiloxaneClass::rhosatL(double T) 
{
	double A = reduce.p*params.molemass/params.R_u/crit.T;
	double B = reduce.p/(reduce.rho*params.R_u/params.molemass*reduce.T);
	double C = reduce.T;
	double D = 2.0/7.0;
    return A/pow(B,1+pow(1-T/C,D));
}
double DodecamethylcyclohexasiloxaneClass::rhosatV(double T) 
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 27.900863 %
	RHS = -0.516885*pow(theta,0.3)
		  -17.855403*pow(theta,0.9)
		  +19.163264*pow(theta,1.5)
		  -72.670303*pow(theta,3.2);
	rho = exp(RHS)*rhoc;
	return rho;
}
double DodecamethylcyclohexasiloxaneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -10.275+4.1393*pow(THETA,2.7012)+-7.7307*pow(THETA,-6.9841)+-5.9825*pow(THETA,-6.5038);
    return exp(crit.T/T*RHS)*crit.p;
}

//MM
HexamethyldisiloxaneClass::HexamethyldisiloxaneClass()
{
    const double n[]={0.0,1.01686012,-2.19713029,0.75443188,-0.68003426,0.19082162,0.0010530133,0.6284595,0.030903042,-0.83948727,-0.20262381,-0.035131597,0.025902341};
    // divided by R_u to give cp0/R_u terms like in Lemmon 2000
	const double u0[]={0.0,6.24140654992885,0.0891626070783569,-5.00332432414229E-05,8.41905535312405E-09};
	const double n0[]={0.0,0,1,2,3};

	std::vector<double> d_v(d,d+sizeof(d)/sizeof(double));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(c,c+sizeof(c)/sizeof(double));
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

    // Critical parameters
	crit.rho = 258.151840734;
	crit.p = 1939.0;
	crit.T = 518.75;
    crit.v = 1.0/crit.rho;

	// Load up a new structure with reducing parameters
	reduce.rho = 304.404388825315;
	reduce.p = 1939.39;
	reduce.T = 518.7;
    reduce.v = 1.0/reduce.rho;
	
	preduce = &reduce;

    // Other fluid parameters
	params.molemass=162.37752;
	params.Ttriple=204.93;
	params.accentricfactor=0.418;
    params.R_u = 8.314472;
	params.ptriple = 6.05765237413746;

    // Limits of EOS
    limits.Tmin = 273;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

	double T0 = 373.400735665744,
		   rho0 = 681.39621296599,
		   R_ = params.R_u/params.molemass,
		   m,
		   c,
		   H0 = -162.96488124949261, /// kJ/kmol
		   S0 = -0.52817673665136178, /// kJ/kmol/K
		   tau0=crit.T/T0,
		   delta0=rho0/crit.rho;
	
	// log(delta)+c+m*tau
	
	/// c is the constant term
	c=-S0/R_-1+log(tau0/delta0);/*<< from the leading term*/

	/// m multiplies the tau term in the leading term (slope)
	m=H0/(R_*reduce.T); /*<< from the leading term */

	phi_BC * phi0_lead_ = new phi0_lead(c,m);
	phi0list.push_back(phi0_lead_);

	phi0list.push_back(new phi0_logtau(-1.0));

	phi0list.push_back(new phi0_cp0_constant(u0[1],reduce.T,T0));
	
	phi0list.push_back(new phi0_cp0_poly(u0_v,n0_v,reduce.T,T0,2,4));

    EOSReference.assign("P. Colonna, N.R. Nannan, A. Guardone, E.W. Lemmon, \"Multiparameter equations of state for selected siloxanes\", Fluid Phase Equilibria 244 (2006) 193–211.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("MM");
    aliases.push_back(std::string("Hexamethyldisiloxane")); 
    REFPROPname.assign("MM");
}
double HexamethyldisiloxaneClass::rhosatL(double T) 
{
	double A = reduce.p*params.molemass/params.R_u/crit.T;
	double B = reduce.p/(reduce.rho*params.R_u/params.molemass*reduce.T);
	double C = reduce.T;
	double D = 2.0/7.0;
    double rho = A/pow(B,1+pow(1-T/C,D));
	return rho;
}
double HexamethyldisiloxaneClass::rhosatV(double T) 
{
	double rhoc = 304.404388825315;
	double theta = 1 - T/reduce.T;
	double RHS,rho;

	// Max error is 1.832708 % between 300K and crit temp
	RHS = -2.219346*pow(theta,0.3)
		  -5.000191*pow(theta,0.9)
		  -2.252852*pow(theta,1.5)
		  -31.605209*pow(theta,3.2);
	rho = exp(RHS)*rhoc;
	return rho;
}
double HexamethyldisiloxaneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -7.338-0.3093*pow(THETA,4.1393)+0.20125*pow(THETA,-7.7307)-13.455*pow(THETA,-5.9825);
    return exp(crit.T/T*RHS)*crit.p;
}

//MD4M
TetradecamethylhexasiloxaneClass::TetradecamethylhexasiloxaneClass()
{
    const double n[]={0.0,1.18492421,-1.87465636,-0.06571351,-0.61812689,0.19535804,0.0005067874,1.23544082,0.049462708,-0.73685283,-0.19991438,-0.055118673,0.028325885};
    // divided by R_u to give cp0/R_u terms like in Lemmon 2000
	const double u0[]={0.0,-2.41398371417933,0.268026640777671,-0.000157724988429812,3.44219091723443E-08};
	const double n0[]={0.0,0,1,2,3};

	std::vector<double> d_v(d,d+sizeof(d)/sizeof(double));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(c,c+sizeof(c)/sizeof(double));
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

    // Critical parameters
	crit.rho = 278.177742672768;
	crit.p = 877.0;
	crit.T =653.2;
    crit.v = 1.0/crit.rho;

	// Load up a new structure with reducing parameters
	reduce.rho = 285.657653221363;
	reduce.p = 877.47;
	reduce.T = 653.2;
    reduce.v = 1.0/reduce.rho;
	
	preduce = &reduce;

    // Other fluid parameters
	params.molemass=458.99328;
	params.Ttriple=214.15;
	params.accentricfactor=0.836;
    params.R_u = 8.314472;
	// Caution calculated at 300 K
	params.ptriple = 0.00109337695354083;


    // Limits of EOS
    limits.Tmin = 273;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

	double T0 = 532.723355998617,
		   rho0 = 651.644500503887,
		   R_ = params.R_u/params.molemass,
		   m,
		   c,
		   H0 = 103.76393935240931, /// kJ/kmol
		   S0 = 0.11760337078430246, /// kJ/kmol/K
		   tau0=crit.T/T0,
		   delta0=rho0/crit.rho;
	
	// log(delta)+c+m*tau
	
	/// c is the constant term
	c=-S0/R_-1+log(tau0/delta0);/*<< from the leading term*/

	/// m multiplies the tau term in the leading term (slope)
	m=H0/(R_*reduce.T); /*<< from the leading term */

	phi_BC * phi0_lead_ = new phi0_lead(c,m);
	phi0list.push_back(phi0_lead_);

	phi_BC * phi0_logtau_ = new phi0_logtau(-1.0);
	phi0list.push_back(phi0_logtau_);

	phi0list.push_back(new phi0_cp0_constant(u0_v[1],reduce.T,T0));

	phi_BC * phi0_cp0_poly_ = new phi0_cp0_poly(u0_v,n0_v,reduce.T,T0,2,4);
	phi0list.push_back(phi0_cp0_poly_);

    EOSReference.assign("P. Colonna, N.R. Nannan, A. Guardone, E.W. Lemmon, \"Multiparameter equations of state for selected siloxanes\", Fluid Phase Equilibria 244 (2006) 193–211.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("MD4M");
    aliases.push_back(std::string("Tetradecamethylhexasiloxane")); 
    REFPROPname.assign("MD4M");
}
double TetradecamethylhexasiloxaneClass::rhosatL(double T) 
{
	double A = reduce.p*params.molemass/params.R_u/crit.T;
	double B = reduce.p/(reduce.rho*params.R_u/params.molemass*reduce.T);
	double C = reduce.T;
	double D = 2.0/7.0;
    double rho = A/pow(B,1+pow(1-T/C,D));
	return rho;
}
double TetradecamethylhexasiloxaneClass::rhosatV(double T) 
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 17.785355 %
	RHS = -1.140416*pow(theta,0.3)
		  -14.921997*pow(theta,0.9)
		  +13.506731*pow(theta,1.5)
		  -70.413994*pow(theta,3.2);
	rho = exp(RHS)*rhoc;
	return rho;
}
double TetradecamethylhexasiloxaneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -10.2921+3.3035*pow(THETA,-0.3093)+-8.0592*pow(THETA,0.20125)+-2.4366*pow(THETA,-13.455);
    return exp(crit.T/T*RHS)*crit.p;
}

//D4
OctamethylcyclotetrasiloxaneClass::OctamethylcyclotetrasiloxaneClass()
{
    const double n[]={0.0,1.05392408,-2.22981918,0.77573923,-0.6937405,0.18721557,0.0004219333,0.70301835,0.047851888,-0.8025348,-0.18968872,-0.022211781,0.0060103354};
    // divided by R_u to give cp0/R_u terms like in Lemmon 2000
	const double u0[]={0.0,-2.19568963609475,0.171652511428266,-0.000119093551580906,3.60816657991031E-08};
	const double n0[]={0.0,0,1,2,3};

	std::vector<double> d_v(d,d+sizeof(d)/sizeof(double));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(c,c+sizeof(c)/sizeof(double));
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

    // Critical parameters
	crit.rho = 305.78949222528007;
	crit.p = 1332.0;
	crit.T = 586.5;
    crit.v = 1.0/crit.rho;

	// Load up a new structure with reducing parameters
	reduce.rho = 307.033590673606;
	reduce.p = 1332;
	reduce.T = 586.49127187;
    reduce.v = 1.0/reduce.rho;
	
	preduce = &reduce;

    // Other fluid parameters
	params.molemass=296.61576;
	params.Ttriple=290.25;
	params.accentricfactor=0.592;
    params.R_u = 8.314472;
	// Caution calculated at 300 K 
	params.ptriple = 0.147791735066714;

    // Limits of EOS
    limits.Tmin = 273;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

	double T0 = 448.503917912145,
		   rho0 = 764.965030584005,
		   R_ = params.R_u/params.molemass,
		   m,
		   c,
		   H0 = -230.50524687772113, /// kJ/kmol
		   S0 = -0.86432336859980707, /// kJ/kmol/K
		   tau0=crit.T/T0,
		   delta0=rho0/crit.rho;
	
	// log(delta)+c+m*tau
	
	/// c is the constant term
	c=-S0/R_-1+log(tau0/delta0);/*<< from the leading term*/

	/// m multiplies the tau term in the leading term (slope)
	m=H0/(R_*reduce.T); /*<< from the leading term */

	phi_BC * phi0_lead_ = new phi0_lead(c,m);
	phi0list.push_back(phi0_lead_);

	phi_BC * phi0_logtau_ = new phi0_logtau(-1.0);
	phi0list.push_back(phi0_logtau_);

	phi0list.push_back(new phi0_cp0_constant(u0_v[1],reduce.T,T0));

	phi_BC * phi0_cp0_poly_ = new phi0_cp0_poly(u0_v,n0_v,reduce.T,T0,2,4);
	phi0list.push_back(phi0_cp0_poly_);

    EOSReference.assign("P. Colonna, N.R. Nannan, A. Guardone, E.W. Lemmon, \"Multiparameter equations of state for selected siloxanes\", Fluid Phase Equilibria 244 (2006) 193–211.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("D4");
    aliases.push_back(std::string("Octamethylcyclotetrasiloxane")); 
    REFPROPname.assign("D4");
}
double OctamethylcyclotetrasiloxaneClass::rhosatL(double T) 
{
	double A = reduce.p*params.molemass/params.R_u/crit.T;
	double B = reduce.p/(reduce.rho*params.R_u/params.molemass*reduce.T);
	double C = reduce.T;
	double D = 2.0/7.0;
    double rho = A/pow(B,1+pow(1-T/C,D));
	return rho;
}
double OctamethylcyclotetrasiloxaneClass::rhosatV(double T) 
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 6.584655 %
	RHS = -1.507450*pow(theta,0.3)
		  -9.626093*pow(theta,0.9)
		  +4.782647*pow(theta,1.5)
		  -50.192115*pow(theta,3.2);
	rho = exp(RHS)*rhoc;
	return rho;
}
double OctamethylcyclotetrasiloxaneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -8.8952+2.5694*pow(THETA,3.3035)+-6.3275*pow(THETA,-8.0592)+-3.5483*pow(THETA,-2.4366);
    return exp(crit.T/T*RHS)*crit.p;
}

//D5
DecamethylcyclopentasiloxaneClass::DecamethylcyclopentasiloxaneClass()
{
    const double n[]={0.0,1.40844725,-2.29248044,0.42851607,-0.73506382,0.16103808,0.00029643278,0.82412481,0.15214274,-0.6849589,-0.055703624,0.013055391,-0.031853761};
    // divided by R_u to give cp0/R_u terms like in Lemmon 2000
	const double u0[]={0.0,-4.19725991019033,0.223886736283434,-0.000168790032608204,6.01361096651718E-08};
	const double n0[]={0.0,0,1,2,3};

	std::vector<double> d_v(d,d+sizeof(d)/sizeof(double));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(c,c+sizeof(c)/sizeof(double));
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

    // Critical parameters
	crit.rho = 304.90928495748;
	crit.p = 1160;
	crit.T = 619.15;
    crit.v = 1.0/crit.rho;

	// Load up a new structure with reducing parameters
	reduce.rho = 292.570762680819;
	reduce.p = 1161.46;
	reduce.T = 619.23462341;
    reduce.v = 1.0/reduce.rho;
	
	preduce = &reduce;

    // Other fluid parameters
	params.molemass=370.7697;
	params.Ttriple=226;
	params.accentricfactor=0.658;
    params.R_u = 8.314472;

	// Caution calculated at 300 K since Ttriple below min temp
	params.ptriple = 0.0275116652817263;

    // Limits of EOS
    limits.Tmin = 273;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

	double T0 = 484.050286854327,
		   rho0 = 727.981991948461,
		   R_ = params.R_u/params.molemass,
		   m,
		   c,
		   H0 = 115.63073039269375, /// kJ/kmol
		   S0 = 0.14017795179966980, /// kJ/kmol/K
		   tau0=crit.T/T0,
		   delta0=rho0/crit.rho;
	
	// log(delta)+c+m*tau
	
	/// c is the constant term
	c=-S0/R_-1+log(tau0/delta0);/*<< from the leading term*/

	/// m multiplies the tau term in the leading term (slope)
	m=H0/(R_*reduce.T); /*<< from the leading term */

	phi_BC * phi0_lead_ = new phi0_lead(c,m);
	phi0list.push_back(phi0_lead_);

	phi_BC * phi0_logtau_ = new phi0_logtau(-1.0);
	phi0list.push_back(phi0_logtau_);

	phi0list.push_back(new phi0_cp0_constant(u0_v[1],reduce.T,T0));

	phi_BC * phi0_cp0_poly_ = new phi0_cp0_poly(u0_v,n0_v,reduce.T,T0,2,4);
	phi0list.push_back(phi0_cp0_poly_);

    EOSReference.assign("P. Colonna, N.R. Nannan, A. Guardone, E.W. Lemmon, \"Multiparameter equations of state for selected siloxanes\", Fluid Phase Equilibria 244 (2006) 193–211.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("D5");
    aliases.push_back(std::string("Decamethylcyclopentasiloxane")); 
    REFPROPname.assign("D5");
}
double DecamethylcyclopentasiloxaneClass::rhosatL(double T) 
{
	double A = reduce.p*params.molemass/params.R_u/crit.T;
	double B = reduce.p/(reduce.rho*params.R_u/params.molemass*reduce.T);
	double C = reduce.T;
	double D = 2.0/7.0;
    double rho = A/pow(B,1+pow(1-T/C,D));
	return rho;
}
double DecamethylcyclopentasiloxaneClass::rhosatV(double T) 
{
	double rhoc = reduce.rho;
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 8.958681 %
	RHS = -1.163417*pow(theta,0.3)
		  -11.715446*pow(theta,0.9)
		  +7.824346*pow(theta,1.5)
		  -54.839101*pow(theta,3.2);
	rho = exp(RHS)*rhoc;
	return rho;
}
double DecamethylcyclopentasiloxaneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -9.4473+3.0697*pow(THETA,2.5694)+-6.9411*pow(THETA,-6.3275)+-2.1692*pow(THETA,-3.5483);
    return exp(crit.T/T*RHS)*crit.p;
}