
#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Siloxanes.h"

//MDM
OctamethyltrisiloxaneClass::OctamethyltrisiloxaneClass()
{
	const double d[] =
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

	const double t[] =
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

	const double l[] =
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
    const double n[]={0.0,1.19735372,-2.40380622,0.3256564,-0.19971259,0.11206277,0.00015893999,0.51234323,-0.020660361,-0.38978114,-0.1186931,-0.037203537,0.018359984};
    const double u0[]={0.0,275.1,612.9,1829.6,413.0,802.6};

	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));

    // Critical parameters
    crit.rho = 256.739908797209;
    crit.p = PressureUnit(1415, UNIT_KPA);
    crit.T = 564.09;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 236.53146;
    params.Ttriple = 187.2;
    params.accentricfactor = 0.5297;
    params.R_u = 8.314472;
	params.ptriple = 7.9911106087e-07;

    // Limits of EOS
    limits.Tmin = params.Ttriple;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phirlist.push_back(new phir_power(n,d,t,l,1,12,13));

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

    EOSReference.assign("Colonna, P., N.R. Nannan, A. Guardone, \"Multiparameter equations of state for siloxanes: [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,...,3, and [O-Si-(CH3)2]6,\", Fluid Phase Equilibria 263 (2008) 115-130.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("MDM");
    aliases.push_back(std::string("Octamethyltrisiloxane")); 
    aliases.push_back(std::string("OCTAMETHYLTRISILOXANE"));
    REFPROPname.assign("MDM");

	ECSReferenceFluid = "Nitrogen";

	BibTeXKeys.EOS = "Colonna-FPE-2008";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double OctamethyltrisiloxaneClass::rhosatL(double T) 
{
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 0.302445 %
	RHS = +0.078367*pow(theta,0.000000)+1.349838*pow(theta,0.296634)-0.005638*pow(theta,3.479858)-0.183289*pow(theta,4.132701)+0.118136*pow(theta,5.130423)+0.416073*pow(theta,6.155667);
	rho = exp(RHS)*reduce.rho;
	return rho;
}
double OctamethyltrisiloxaneClass::rhosatV(double T) 
{
	// Max error is 0.455581363816 % between 187.2 and 564.089999 K
	const double ti[]={0, 0.14700000000000002, 0.3515, 0.39699999999999996, 1.5, 3.8333333333333335, 16.166666666666668};
    const double Ni[]={0, -1.2376429620882381, 18.16254939596832, -21.7298671811246, -2.318395979212784, -9.020806554399309, -10.84321486627527};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
double OctamethyltrisiloxaneClass::psat(double T) 
{
    double theta = 1-T/reduce.T;
	double RHS,p;

	// Max error is 1.053362 %
	RHS = -7.765160*pow(theta,0.984389)-5.002622*pow(theta,3.506411)-2.478207*pow(theta,4.186938)-1.023016*pow(theta,4.808272)-0.201986*pow(theta,5.513089);
	p = exp(RHS*reduce.T/T)*reduce.p.Pa;
	return p;
}

//MD2M
DecamethyltetrasiloxaneClass::DecamethyltetrasiloxaneClass()
{
	const double d[] =
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

	const double t[] =
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

	const double l[] =
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
    const double n[]={0.0,1.33840331,-2.62939393,0.4398383,-0.53496715,0.1818844,0.00040774609,1.13444506,0.05774631,-0.5917498,-0.11020225,-0.034942635,0.007646298};
    const double u0[]={0.0,331.9,777.1,1813.8,521.4,795.1};

	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));

    // Critical parameters
    crit.rho=284.171639662027;
	crit.p= PressureUnit(1227, UNIT_KPA);
	crit.T=599.4;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
	params.molemass=310.685;
	params.Ttriple=205.2;
	params.accentricfactor=0.668;
    params.R_u = 8.314472;
	params.ptriple = 4.79509942923e-07;

    // Limits of EOS
    limits.Tmin = params.Ttriple;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phirlist.push_back(new phir_power(n,d,t,l,1,12,13));

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

    EOSReference.assign("Colonna, P., N.R. Nannan, A. Guardone, \"Multiparameter equations of state for siloxanes: [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,...,3, and [O-Si-(CH3)2]6,\", Fluid Phase Equilibria 263 (2008) 115-130.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("MD2M");
    aliases.push_back(std::string("Decamethyltetrasiloxane")); 
    aliases.push_back(std::string("DECAMETHYLTETRASILOXANE"));
    REFPROPname.assign("MD2M");

	ECSReferenceFluid = "Nitrogen";

	BibTeXKeys.EOS = "Colonna-FPE-2008";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double DecamethyltetrasiloxaneClass::rhosatL(double T) 
{
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 0.222790 %
	RHS = +1.308255*pow(theta,0.286019)+0.298564*pow(theta,2.774999)-0.249494*pow(theta,3.604257)-0.208559*pow(theta,4.311892)+0.062389*pow(theta,5.214777)+0.366547*pow(theta,6.192231);
	rho = exp(RHS)*reduce.rho;
	return rho;
}
double DecamethyltetrasiloxaneClass::rhosatV(double T) 
{
	// Max error is  0.169267097925 % between 205.2 and 599.399999 K
	const double ti[]={0, 0.39199999999999996, 0.39349999999999996, 0.39399999999999996, 1.5, 3.3333333333333335, 18.0};
    const double Ni[]={0, -52837.662064867334, 213259.90052518147, -160428.17073987963, -0.66311675635221967, -9.8035152206617866, -34.641954872452423};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
double DecamethyltetrasiloxaneClass::psat(double T) 
{
	double theta = 1-T/reduce.T;
	double RHS,p;

	// Max error is 0.653417 %
	RHS = -10.255244*pow(theta,1.003893)-14.680968*pow(theta,8.315832)-8.013482*pow(theta,2.914286)+3.386436*pow(theta,1.418176)+14.914376*pow(theta,9.632505);
	p = exp(RHS*reduce.T/T)*reduce.p.Pa;
	return p;
}

//MD3M
DodecamethylpentasiloxaneClass::DodecamethylpentasiloxaneClass()
{
	const double d[] =
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

	const double t[] =
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

	const double l[] =
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
    const double n[]={0.0,1.20540386,-2.42914797,0.69016432,-0.69268041,0.18506046,0.00031161436,0.99862519,0.074229034,-0.80259136,-0.20865337,-0.036461791,0.019174051};
    const double u0[]={0.0,463.2,957.2,2117.1,738.3,908.5};

	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));

    // Critical parameters
    crit.rho=263.921879135305;
	crit.p= PressureUnit(945, UNIT_KPA);
	crit.T=628.36;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
	params.molemass=384.839;
	params.Ttriple=192;
	params.accentricfactor=0.722;
    params.R_u = 8.314472;
	params.ptriple = 2.05774764356e-10;

    // Limits of EOS
    limits.Tmin = params.Ttriple;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phirlist.push_back(new phir_power(n,d,t,l,1,12,13));

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

    EOSReference.assign("Colonna, P., N.R. Nannan, A. Guardone, \"Multiparameter equations of state for siloxanes: [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,...,3, and [O-Si-(CH3)2]6,\", Fluid Phase Equilibria 263 (2008) 115-130.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("MD3M");
    aliases.push_back(std::string("Dodecamethylpentasiloxane")); 
    aliases.push_back(std::string("DODECAMETHYLPENTASILOXANE"));
    REFPROPname.assign("MD3M");

	ECSReferenceFluid = "Nitrogen";

	BibTeXKeys.EOS = "Colonna-FPE-2008";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double DodecamethylpentasiloxaneClass::psat(double T)
{
	// Max error is  0.0600093438283 % between 192.0 and 628.359999 K
    const double t[]={0, 1.5, 1.0, 1.8333333333333333, 3.1666666666666665, 6.666666666666667, 16.166666666666668};
    const double N[]={0, 1.3506034234657476, -9.2126060941992307, -0.28939038151557195, -9.0673717965357667, -2.9245800043788388, -10.66550052981526};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double DodecamethylpentasiloxaneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.288378 % between 192.000000 K and 628.360000 K
    const double t[] = {0, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.6666666666666667, 1.8333333333333333};
    const double N[] = {0, -2.8756454677232282, 58.702098949984546, -312.49817605073389, 1216.0740391388877, -3953.7402984826108, 9489.6576229693947, -14861.937754088141, 14106.261719203909, -7352.4250932407413, 1616.4437867458232};
    double summer=0,theta;
    theta=1-T/reduce.T;
    	
	for (int i=1; i<=10; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double DodecamethylpentasiloxaneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.165855 % between 192.000000 K and 628.360000 K
    const double t[] = {0, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.6666666666666667, 1.8333333333333333, 2.0, 2.3333333333333335, 2.6666666666666665, 3.0, 3.3333333333333335};
    const double N[] = {0, -8.1412243969392986, 465.78213307742897, -10070.395070453211, 113934.3138003942, -797658.75904416665, 3713183.847053702, -11886054.42914509, 26400426.037794076, -40022635.392553486, 38973041.838086836, -19997978.986398581, 4823720.3256814424, -1716855.6651460216, 473036.175360589, -66568.844348660248};
    double summer=0,theta;
    theta=1-T/reduce.T;	
	for (int i=1; i<=15; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

//D6
DodecamethylcyclohexasiloxaneClass::DodecamethylcyclohexasiloxaneClass()
{
	const double d[] =
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

	const double t[] =
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

	const double l[] =
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
    const double n[]={0.0,1.69156186,-3.37962568,0.38609039,0.064598995,0.10589012,0.000045456825,0.74169279,-0.088102648,-0.17373336,-0.10951368,-0.062695695,0.037459986};
    const double u0[]={0.0,468.7,981.2,1792.1,686.7,786.8};

	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));


    // Critical parameters
	crit.rho=279.095729841367;
	crit.p= PressureUnit(961, UNIT_KPA);
	crit.T=645.78;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
	params.molemass=444.924;
	params.Ttriple=270.2;
	params.accentricfactor=0.736;
    params.R_u = 8.314472;
	params.ptriple = 0.000159753038201;

    // Limits of EOS
    limits.Tmin = params.Ttriple;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phirlist.push_back(new phir_power(n,d,t,l,1,12,13));

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

    EOSReference.assign("Colonna, P., N.R. Nannan, A. Guardone, \"Multiparameter equations of state for siloxanes: [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,...,3, and [O-Si-(CH3)2]6,\", Fluid Phase Equilibria 263 (2008) 115-130.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("D6");
    aliases.push_back(std::string("Dodecamethylcyclohexasiloxane")); 
    aliases.push_back(std::string("DODECAMETHYLCYCLOHEXASILOXANE"));
    REFPROPname.assign("D6");

	ECSReferenceFluid = "Nitrogen";

	BibTeXKeys.EOS = "Colonna-FPE-2008";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double DodecamethylcyclohexasiloxaneClass::rhosatL(double T) 
{
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 0.390948 %
	RHS = +1.508669*pow(theta,0.299083)-0.367377*pow(theta,2.751619)+0.197308*pow(theta,3.522398)+0.393542*pow(theta,4.060807)+0.115383*pow(theta,4.768416)-0.125626*pow(theta,5.916593);
	rho = exp(RHS)*reduce.rho;
	return rho;
}
double DodecamethylcyclohexasiloxaneClass::rhosatV(double T) 
{
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 0.463858 %
	RHS = -2.334877*pow(theta,0.354340)-84.316983*pow(theta,8.854211)-36.472358*pow(theta,6.602777)-9.706856*pow(theta,1.066697)-29.464799*pow(theta,3.643699)-17.991230*pow(theta,3.468228)-10.449427*pow(theta,9.287216);
	rho = exp(RHS)*reduce.rho;
	return rho;
}
double DodecamethylcyclohexasiloxaneClass::psat(double T) 
{
	double theta = 1-T/reduce.T;
	double RHS,p;

	// Max error is 1.754539 %
	RHS = -10.075325*pow(theta,1.012056)-134.690578*pow(theta,9.074796)-43.425090*pow(theta,5.030680)-32.551643*pow(theta,11.673831)-27.321344*pow(theta,2.821284);
	p = exp(RHS)*reduce.p.Pa;
	return p;
}

//MM
HexamethyldisiloxaneClass::HexamethyldisiloxaneClass()
{
	const double d[] =
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

	const double t[] =
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

	const double l[] =
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
    const double n[]={0.0,1.01686012,-2.19713029,0.75443188,-0.68003426,0.19082162,0.0010530133,0.6284595,0.030903042,-0.83948727,-0.20262381,-0.035131597,0.025902341};
    // divided by R_u to give cp0/R_u terms like in Lemmon 2000
	const double u0[]={0.0,6.24140654992885,0.0891626070783569,-5.00332432414229E-05,8.41905535312405E-09};
	const double n0[]={0.0,0,1,2,3};

	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

    // Critical parameters
	crit.rho = 258.151840734;
	crit.p = PressureUnit(1939.0, UNIT_KPA);
	crit.T = 518.75;
    crit.v = 1.0/crit.rho;

	// Load up a new structure with reducing parameters
	reduce.rho = 1.87467076*162.37752;
	reduce.p = PressureUnit(1939.39, UNIT_KPA);
	reduce.T = 518.69997204;
    reduce.v = 1.0/reduce.rho;
	
	preduce = &reduce;

    // Other fluid parameters
	params.molemass = 162.37752;
	params.Ttriple = 204.93;
	params.accentricfactor = 0.418;
    params.R_u = 8.314472;
	params.ptriple = 1.3169852229396559;

    // Limits of EOS
    limits.Tmin = 273;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phirlist.push_back(new phir_power(n,d,t,l,1,12,13));

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

    EOSReference.assign("Colonna, P., N.R. Nannan, A. Guardone, E.W. Lemmon, \"Multiparameter equations of state for selected siloxanes\", Fluid Phase Equilibria 244 (2006) 193-211.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("MM");
    aliases.push_back(std::string("Hexamethyldisiloxane")); 
    aliases.push_back(std::string("HEXAMETHYLDISILOXANE"));
    REFPROPname.assign("MM");

	ECSReferenceFluid = "Nitrogen";

	BibTeXKeys.EOS = "Colonna-FPE-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double HexamethyldisiloxaneClass::rhosatL(double T) 
{
	double theta = 1-T/reduce.T;
	double RHS,rho;

	if (theta < 0) theta = 1e-6;
	// Max error is 0.177687 %
	RHS = +4.701517*pow(theta,0.429578)-3.554857*pow(theta,0.485308)+0.461167*pow(theta,4.050387)-0.438751*pow(theta,4.407893)-0.107601*pow(theta,5.383350)+0.294769*pow(theta,6.326992);
	rho = exp(RHS)*reduce.rho;
	return rho;
}
double HexamethyldisiloxaneClass::rhosatV(double T) 
{
	double theta = 1-T/reduce.T;
	double RHS,rho;

	if (theta < 0) theta = 1e-6;
	// Max error is 0.625617 %
	RHS = -3.349640*pow(theta,0.363375)-3.919482*pow(theta,1.308131)-5.046495*pow(theta,5.251634)-3.395130*pow(theta,4.732761)-1.880353*pow(theta,5.885600)-0.727824*pow(theta,6.806098);
	rho = exp(RHS*reduce.T/T)*reduce.rho;
	return rho;
}
double HexamethyldisiloxaneClass::psat(double T) 
{
    double theta = 1-T/reduce.T;
	double RHS,p;

	if (theta < 0) theta = 1e-6;
	// Max error is 0.636597 %
	RHS = -6.993179*pow(theta,0.968165)-3.963479*pow(theta,3.066906)-1.333772*pow(theta,4.415999)-0.096273*pow(theta,4.801166)+0.499673*pow(theta,5.391500);
	p = exp(RHS*reduce.T/T)*reduce.p.Pa;
	return p;
}

//MD4M
TetradecamethylhexasiloxaneClass::TetradecamethylhexasiloxaneClass()
{
	const double d[] =
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

	const double t[] =
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

	const double l[] =
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
    const double n[]={0.0,1.18492421,-1.87465636,-0.06571351,-0.61812689,0.19535804,0.0005067874,1.23544082,0.049462708,-0.73685283,-0.19991438,-0.055118673,0.028325885};
    // divided by R_u to give cp0/R_u terms like in Lemmon 2000
	const double u0[]={0.0,-2.41398371417933,0.268026640777671,-0.000157724988429812,3.44219091723443E-08};
	const double n0[]={0.0,0,1,2,3};

	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

    // Critical parameters
	crit.rho = 278.177742672768;
	crit.p = PressureUnit(877.0, UNIT_KPA);
	crit.T =653.2;
    crit.v = 1.0/crit.rho;

	// Load up a new structure with reducing parameters
	reduce.rho = 285.657653221363;
	reduce.p = PressureUnit(877.47, UNIT_KPA);
	reduce.T = 653.2;
    reduce.v = 1.0/reduce.rho;
	
	preduce = &reduce;

    // Other fluid parameters
	params.molemass=458.99328;
	params.Ttriple=214.15;
	params.accentricfactor=0.82464714726429245;
    params.R_u = 8.314472;
	// Caution calculated at Tmin
	params.ptriple = 0.00109337696739;


    // Limits of EOS
    limits.Tmin = 300;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phirlist.push_back(new phir_power(n,d,t,l,1,12,13));

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

    EOSReference.assign("Colonna, P., N.R. Nannan, A. Guardone, E.W. Lemmon, \"Multiparameter equations of state for selected siloxanes\", Fluid Phase Equilibria 244 (2006) 193-211.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("MD4M");
    aliases.push_back(std::string("Tetradecamethylhexasiloxane")); 
    aliases.push_back(std::string("TETRADECAMETHYLHEXASILOXANE"));
    REFPROPname.assign("MD4M");

	ECSReferenceFluid = "Nitrogen";

	BibTeXKeys.EOS = "Colonna-FPE-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double TetradecamethylhexasiloxaneClass::rhosatL(double T) 
{
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 0.736426 %
	RHS = +1.323853*pow(theta,0.281503)+0.147420*pow(theta,3.828555)+0.246462*pow(theta,3.896665)-0.346926*pow(theta,4.315933)-0.173105*pow(theta,5.291410)+0.187277*pow(theta,6.286236);
	rho = exp(RHS)*reduce.rho;
	return rho;
}
double TetradecamethylhexasiloxaneClass::rhosatV(double T) 
{
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 0.158764 %
	RHS = -0.411516*pow(theta,0.123844)-6.464782*pow(theta,0.633729)-11.135283*pow(theta,17.959104)-10.395902*pow(theta,3.308400)-2.873559*pow(theta,3.275238)-1.178585*pow(theta,9.838501);
	rho = exp(RHS*reduce.T/T)*reduce.rho;
	return rho;
}
double TetradecamethylhexasiloxaneClass::psat(double T) 
{
    double theta = 1-T/reduce.T;
	double RHS,p;

	// Max error is 0.104803 %
	RHS = -9.945866*pow(theta,0.993284)-113.010891*pow(theta,7.297341)-38.114082*pow(theta,24.782767)-33.140158*pow(theta,4.160418)-21.201196*pow(theta,2.574117);
	p = exp(RHS)*reduce.p.Pa;
	return p;
}

//D4
OctamethylcyclotetrasiloxaneClass::OctamethylcyclotetrasiloxaneClass()
{
	const double d[] =
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

	const double t[] =
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

	const double l[] =
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
    const double n[]={0.0,1.05392408,-2.22981918,0.77573923,-0.6937405,0.18721557,0.0004219333,0.70301835,0.047851888,-0.8025348,-0.18968872,-0.022211781,0.0060103354};
    // divided by R_u to give cp0/R_u terms like in Lemmon 2000
	const double u0[]={0.0,-2.19568963609475,0.171652511428266,-0.000119093551580906,3.60816657991031E-08};
	const double n0[]={0.0,0,1,2,3};

	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

    // Critical parameters
	crit.rho = 305.78949222528007;
	crit.p = PressureUnit(1332.0, UNIT_KPA);
	crit.T = 586.5;
    crit.v = 1.0/crit.rho;

	// Load up a new structure with reducing parameters
	reduce.rho = 307.033590673606;
	reduce.p = PressureUnit(1332, UNIT_KPA);
	reduce.T = 586.49127187;
    reduce.v = 1.0/reduce.rho;
	
	preduce = &reduce;

    // Other fluid parameters
	params.molemass=296.61576;
	params.Ttriple=290.25;
	params.accentricfactor=0.592;
    params.R_u = 8.314472;
	// Caution calculated at Tmin
	params.ptriple = 0.069609237623675213;

    // Limits of EOS
	limits.Tmin = params.Ttriple;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phirlist.push_back(new phir_power(n,d,t,l,1,12,13));

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

    EOSReference.assign("Colonna, P., N.R. Nannan, A. Guardone, E.W. Lemmon, \"Multiparameter equations of state for selected siloxanes\", Fluid Phase Equilibria 244 (2006) 193-211.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("D4");
    aliases.push_back(std::string("Octamethylcyclotetrasiloxane")); 
    aliases.push_back(std::string("OCTAMETHYLCYCLOTETRASILOXANE"));
    REFPROPname.assign("D4");

	ECSReferenceFluid = "Nitrogen";

	BibTeXKeys.EOS = "Colonna-FPE-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double OctamethylcyclotetrasiloxaneClass::rhosatL(double T) 
{
		double theta = 1-T/reduce.T;
		double RHS,rho;

		// Max error is 0.063302 %
		RHS = +0.104571*pow(theta,0.045868)+1.250508*pow(theta,0.294117)+4.015272*pow(theta,4.002063)-2.091843*pow(theta,3.489495)-2.218802*pow(theta,5.461110)-1.211494*pow(theta,6.781263);
		rho = exp(RHS)*reduce.rho;
		return rho;
}
double OctamethylcyclotetrasiloxaneClass::rhosatV(double T) 
{
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 0.539080 %
	RHS = -0.836967*pow(theta,0.246119)-4.887839*pow(theta,0.626148)-5.741624*pow(theta,2.543434)-3.669516*pow(theta,4.023274)-2.021077*pow(theta,5.824788)-0.819746*pow(theta,6.842893);
	rho = exp(RHS*reduce.T/T)*reduce.rho;
	return rho;
}
double OctamethylcyclotetrasiloxaneClass::psat(double T) 
{
	// Max error is  0.0796509289458 % between 290.25 and 586.499999 K
	const double ti[]={0, 0.352, 0.3745, 1.0, 3.5, 4.333333333333333, 5.5};
    const double Ni[]={0, 0.80415183721176908, -0.93068397535046898, -7.8477917073533066, -13.625622539026633, 13.019220440087111, -10.876350518822491};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(crit.T/T*summer);
}

//D5
DecamethylcyclopentasiloxaneClass::DecamethylcyclopentasiloxaneClass()
{
	const double d[] =
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

	const double t[] =
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

	const double l[] =
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
    const double n[]={0.0,1.40844725,-2.29248044,0.42851607,-0.73506382,0.16103808,0.00029643278,0.82412481,0.15214274,-0.6849589,-0.055703624,0.013055391,-0.031853761};
    // divided by R_u to give cp0/R_u terms like in Lemmon 2000
	const double u0[]={0.0,-4.19725991019033,0.223886736283434,-0.000168790032608204,6.01361096651718E-08};
	const double n0[]={0.0,0,1,2,3};
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));
	std::vector<double> d_v(d,d+sizeof(d)/sizeof(double));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(l,l+sizeof(l)/sizeof(double));

    // Critical parameters
	crit.rho = 304.90928495748;
	crit.p = PressureUnit(1160, UNIT_KPA);
	crit.T = 619.15;
    crit.v = 1.0/crit.rho;

	// Load up a new structure with reducing parameters
	reduce.rho = 292.570762680819;
	reduce.p = PressureUnit(1161.46, UNIT_KPA);
	reduce.T = 619.23462341;
    reduce.v = 1.0/reduce.rho;
	
	preduce = &reduce;

    // Other fluid parameters
	params.molemass=370.7697;
	params.Ttriple=226;
	params.accentricfactor=0.658;
    params.R_u = 8.314472;

	// Calculated at Tmin
	params.ptriple = 0.0275116653286;

    // Limits of EOS
    limits.Tmin = 300;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    phirlist.push_back(new phir_power(n,d,t,l,1,12,13));

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

    EOSReference.assign("Colonna, P., N.R. Nannan, A. Guardone, E.W. Lemmon, \"Multiparameter equations of state for selected siloxanes\", Fluid Phase Equilibria 244 (2006) 193-211.");
    TransportReference.assign("Using ECS in fully predictive mode");

    name.assign("D5");
    aliases.push_back(std::string("Decamethylcyclopentasiloxane")); 
    aliases.push_back(std::string("DECAMETHYLCYCLOPENTASILOXANE"));
    REFPROPname.assign("D5");

	ECSReferenceFluid = "Nitrogen";

	BibTeXKeys.EOS = "Colonna-FPE-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double DecamethylcyclopentasiloxaneClass::rhosatL(double T) 
{
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 0.260518 %
	RHS = +1.080094*pow(theta,0.198987)+0.585719*pow(theta,1.050756)-0.605915*pow(theta,3.284144)-0.279332*pow(theta,3.681234)+0.797515*pow(theta,5.186809)+1.061822*pow(theta,6.108004);
	rho = exp(RHS)*reduce.rho;
	return rho;
}
double DecamethylcyclopentasiloxaneClass::rhosatV(double T) 
{
	double theta = 1-T/reduce.T;
	double RHS,rho;

	// Max error is 0.660365 %
	RHS = -0.752767*pow(theta,0.311678)-5.032494*pow(theta,0.603043)-8.781241*pow(theta,9.277681)-7.681111*pow(theta,2.635003)-4.604776*pow(theta,6.635875)-2.743916*pow(theta,7.810733);
	rho = exp(RHS*reduce.T/T)*reduce.rho;
	return rho;
}
double DecamethylcyclopentasiloxaneClass::psat(double T) 
{
	double theta = 1-T/reduce.T;
	double RHS,p;

	// Max error is 1.040101 %
	RHS = -9.243668*pow(theta,0.996513)-98.611662*pow(theta,6.934223)-32.083059*pow(theta,2.871414)-18.793106*pow(theta,8.376677)-12.211323*pow(theta,9.708690)-6.644894*pow(theta,10.070804);
	p = exp(RHS)*reduce.p.Pa;
	return p;
}
