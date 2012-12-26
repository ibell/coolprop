#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Xylene_EthylBenzene.h"

static char EOSstr [] = "Yong Zhou, Jiangtao Wu, Eric W. Lemmon, \"Thermodynamic Properties of o-Xylene, m-Xylene, p-Xylene, and Ethylbenzene\", J. Phys. Chem. Ref. Data, Vol. 41, No. 2, 2012";

oXyleneClass::oXyleneClass()
{
	double n[] = {0.0, 0.0036765156, -0.13918171, 0.014104203, 1.5398899, -2.3600925, -0.44359159, 0.19596977, -1.0909408, -0.21890801, 1.1179223, -0.93563815, -0.018102996, 1.4172368, -0.57134695, -0.081944041, -40.682878};
	double t[] = {0, 1, 0.6, 0.91, 0.3, 0.895, 1.167, 0.435, 2.766, 3.8, 1.31, 3, 0.77, 1.41, 4.8, 1.856, 2};
	double d[] = {0, 5, 1, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3, 3};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.1723, 1.095, 1.6166, 20.4};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.442, 1.342, 3, 450};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.2655, 0.3959, 0.7789, 1.162};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.552, 0.728, 0.498, 0.894};

	//Critical parameters
	crit.rho = 2.6845*106.165; //[kg/m^3]
	crit.p = 3737.5; //[kPa]
	crit.T = 630.259; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 106.165;
	params.Ttriple = 247.985;
	params.accentricfactor = 0.312;
	params.R_u = 8.314472;
	params.ptriple = 0.022805778456334729;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,12,17));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 13,16,17));

	const double a1 = 10.137376, a2= -0.91282993, c0 = 3.748798;
	phi0list.push_back(new phi0_lead(a1,a2));
	phi0list.push_back(new phi0_logtau(c0-1));

	const double u0[] = {0, 225/crit.T, 627/crit.T, 1726/crit.T, 4941/crit.T};
	const double v0[] = {0, 4.754892, 6.915052, 25.84813, 10.93886};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("o-Xylene");
	aliases.push_back("oXylene");
	REFPROPname.assign("OXYLENE");
}
double oXyleneClass::psat(double T)
{
    // Maximum absolute error is 0.167477 % between 247.985001 K and 630.258990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.4975582403006129, 1.1122769997457862, 0.51915348595803801, -6.8673796049531219, 5.6010661099890999, -5.6960509837039917 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double oXyleneClass::rhosatL(double T)
{
    // Maximum absolute error is 1.758248 % between 247.985001 K and 630.258990 K
    const double ti[]={0,0.57845021125158513, 0.73462892915017219, 1.9704181229283253, 2.3913537297744285, 2.1442227226496251};
    const double Ni[]={0,9.389944423330256, -10.065201566952689, 36.995144332929428, 19.104153002871463, -53.998326006477313};
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
double oXyleneClass::rhosatV(double T)
{
    // Maximum absolute error is 3.076635 % between 247.985001 K and 630.258990 K
    const double ti[]={0,0.55197055685787888, 3.053372826026405, 3.0213687610765181, 2.9985688560578714, 2.8841097402929057};
    const double Ni[]={0,-5.2855633441813108, -66889.17640707799, 194634.66961501876, -134548.71734816648, 6797.2028366417326};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

mXyleneClass::mXyleneClass()
{
	double n[] = {0.0, 0.000012791017, 0.041063111, 1.505996, -2.3095875, -0.46969, 0.171031, -1.001728, -0.3945766, 0.6970578, -0.3002876, -0.024311, 0.815488, -0.330647, -0.123393, -0.54661};
	double t[] = {0, 1, 0.91, 0.231, 0.772, 1.205, 0.323, 2.7, 3.11, 0.768, 4.1, 0.818, 2, 2.9, 3.83, 0.5};
	double d[] = {0, 8, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3, 3};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0244, 1.3788, 0.9806, 6.3563};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.66, 1.9354, 1.0323, 78};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.1013, 0.6515, 0.4975, 1.26};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.713, 0.9169, 0.6897, 0.7245};

	//Critical parameters
	crit.rho = 2.665*106.165; //[kg/m^3]
	crit.p = 3534.6; //[kPa]
	crit.T = 616.89; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 106.165;
	params.Ttriple = 225.3;
	params.accentricfactor = 0.326;
	params.R_u = 8.314472;
	params.ptriple = 0.003123267841543599;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,11,16));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 12,15,16));

	const double a1 = 12.652887, a2 = 0.45975624, c0 = 2.169909;
	phi0list.push_back(new phi0_lead(a1,a2));
	phi0list.push_back(new phi0_logtau(c0-1));

	const double u0[] = {0, 160/crit.T, 190/crit.T, 1333/crit.T, 3496/crit.T};
	const double v0[] = {0, 4.44312, 2.862794, 24.83298, 16.26077};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("m-Xylene");
	aliases.push_back("mXylene");
	REFPROPname.assign("MXYLENE");
}
double mXyleneClass::psat(double T)
{
    // Maximum absolute error is 0.018452 % between 225.300001 K and 616.889990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.5970930777234598, 1.4876776684391573, -0.96616236674879985, -3.166803942530092, -0.22852828380711804, -1.1851878442164354 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double mXyleneClass::rhosatL(double T)
{
    // Maximum absolute error is 1.173679 % between 225.300001 K and 616.889990 K
    const double ti[]={0,0.35271681982885156, 0.81485029750189364, 0.016313449442244117, 2.3949502437791548, 2.4104177192823952};
    const double Ni[]={0,1.7374058871186462, -0.61674066715524523, 0.025703428745562459, 13.071354988026734, -12.862683737576914};
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
double mXyleneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.525916 % between 225.300001 K and 616.889990 K
    const double ti[]={0,0.21794065671105323, 1.4319225117247736, 0.92161626748727543, 1.0679705910886546, 4.2165026050064096};
    const double Ni[]={0,-0.87842179198146009, -11.595052043783022, -30.697840455426803, 36.548554391098151, -4.5821764606354369};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}

pXyleneClass::pXyleneClass()
{
	double n[] = {0.0, 0.0010786811, -0.103161822, 0.0421544125, 1.47865376, -2.4266, -0.46575193, 0.190290995, -1.06376565, -0.209934069, 1.25159879, -0.951328356, -0.0269980032, 1.3710318, -0.494160616, -0.0724317468, -3.69464746};
	double t[] = {0, 1, 0.83, 0.83, 0.281, 0.932, 1.1, 0.443, 2.62, 2.5, 1.2, 3, 0.778, 1.13, 4.5, 2.2, 2};
	double d[] = {0, 5, 1, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3, 3};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.179, 1.065, 1.764, 13.675};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.445, 1.483, 4.971, 413};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.267, 0.4242, 0.864, 1.1465};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.54944, 0.7234, 0.4926, 0.8459};

	//Critical parameters
	crit.rho = 2.69392*106.165; //[kg/m^3]
	crit.p = 3531.5; //[kPa]
	crit.T = 616.168; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 106.165;
	params.Ttriple = 286.4;
	params.accentricfactor = 0.324;
	params.R_u = 8.314472;
	params.ptriple = 0.580085039148721;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,12,17));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 13,16,17));

	const double a1 = 5.9815241, a2 = -0.52477835, c0 = 5.2430504;
	phi0list.push_back(new phi0_lead(a1,a2));
	phi0list.push_back(new phi0_logtau(c0-1));

	const double u0[] = {0, 414/crit.T, 1256/crit.T, 2649/crit.T, 6681/crit.T};
	const double v0[] = {0, 5.2291378, 19.549862, 16.656178, 5.9390291};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("p-Xylene");
	aliases.push_back("pXylene");
	REFPROPname.assign("PXYLENE");
}
double pXyleneClass::psat(double T)
{
    // Maximum absolute error is 0.072046 % between 286.400001 K and 616.167990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.7144829641292896, 1.5473929648972096, 0.043172647969443374, -6.6834249677972375, 6.398855445871666, -7.8572899524361315 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double pXyleneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.476482 % between 286.400001 K and 616.167990 K
    const double ti[]={0,0.36922931174880702, 1.6084186059964691, 1.685317285952024, 1.6737093375890244, 1.6358473659507031};
    const double Ni[]={0,1.9185912843031041, -9811.8226960758384, 30281.812695478082, -46654.337197534507, 26183.830801897715};
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
double pXyleneClass::rhosatV(double T)
{
    // Maximum absolute error is 2.036894 % between 286.400001 K and 616.167990 K
    const double ti[]={0,0.060186161740589092, 0.53381625316033809, 5.3043338091148762, 10.872745924810694, 3.0052865665149087};
    const double Ni[]={0,-0.067596669375935142, -4.962397082005725, -0.97523231532511534, -9.7466146643189031, -4.3322511157907018};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}






EthylBenzeneClass::EthylBenzeneClass()
{
	double n[] = {0.0, 0.0018109418, -0.076824284, 0.041823789, 1.5059649, -2.4122441, -0.47788846, 0.18814732, -1.0657412, -0.20797007, 1.1222031, -0.99300799, -0.027300984, 1.3757894, -0.44477155, -0.07769742, -2.16719};
	double t[] = {0, 1, 1, 0.92, 0.27, 0.962, 1.033, 0.513, 2.31, 3.21, 1.26, 2.29, 1, 0.6, 3.6, 2.1, 0.5};
	double d[] = {0, 5, 1, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3, 3};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.178, 1.07, 1.775, 15.45};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.437, 1.488, 4, 418.6};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.2667, 0.4237, 0.8573, 1.15};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5494, 0.7235, 0.493, 0.8566};

	//Critical parameters
	crit.rho = 2.741016*106.165; //[kg/m^3]
	crit.p = 3622.4; //[kPa]
	crit.T = 617.12; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 106.165;
	params.Ttriple = 178.2;
	params.accentricfactor = 0.304;
	params.R_u = 8.314472;
	params.ptriple = 4.0029622500330704e-006;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,12,17));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 13,16,17));

	const double a1 = 5.70409, a2 = -0.52414353, c0 = 5.2557889;
	phi0list.push_back(new phi0_lead(a1,a2));
	phi0list.push_back(new phi0_logtau(c0-1));

	const double u0[] = {0, 585/crit.T, 4420/crit.T, 1673/crit.T};
	const double v0[] = {0, 9.7329909, 11.201832, 25.440749};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,3));

	EOSReference.assign(EOSstr);
	TransportReference.assign("Using ECS in fully predictive mode");

	name.assign("EthylBenzene");
	REFPROPname.assign("EBENZENE");
}
double EthylBenzeneClass::psat(double T)
{
    // Maximum absolute error is 0.327455 % between 178.200001 K and 617.119990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.8721305383390021, 2.8514937477013724, -3.1520525089862179, -0.56815525664475219, -2.7058200340664929, 0.0060596500475032797 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double EthylBenzeneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.103278 % between 178.200001 K and 617.119990 K
    const double ti[]={0,0.46800079198111388, 1.2365861994755625, 0.76851332608867118, 0.7740262636007067, 3.774095180029045};
    const double Ni[]={0,6.0838517775182845, -2.0981534186455195, -373.30274973471006, 370.55311131275505, 0.1569655724457043};
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
double EthylBenzeneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.430690 % between 178.200001 K and 617.119990 K
    const double ti[]={0,0.4150387035308935, 0.5921334814842657, 2.4316112759502557, 3.9039286150614494, 2.4445797475912574};
    const double Ni[]={0,-2.6138578436692339, -1.8758825160832402, -250.15421395285367, -7.1651085385105544, 250.60121693432578};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
