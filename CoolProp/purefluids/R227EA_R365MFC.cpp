#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "R227EA_R365MFC.h"

R227EAClass::R227EAClass()
{
	double n[] = {0.0, 2.024341, -2.60593, 0.4957216, -0.824082, 0.06543703, -1.02461, 0.6247065, 0.2997521, -0.353917, -1.232043, -0.8824483, 0.1349661, -0.2662928, 0.1764733, 0.01536163, -0.004667185, -11.70854, 0.9114512};
	double t[] = {0, 0.34, 0.77, 0.36, 0.9, 1, 2.82, 2.1, 0.9, 1.13, 3.8, 2.75, 1.5, 2.5, 2.5, 5.4, 4, 1, 3.5};
	double d[] = {0, 1, 1, 2, 2, 4, 1, 3, 6, 6, 2, 3, 1, 2, 1, 1, 4, 2, 1};
	double c[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.83, 2.19, 2.44, 3.65, 8.88, 8.23, 2.01};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.72, 5.2, 2.31, 1.02, 5.63, 50.9, 1.56};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.414, 1.051, 1.226, 1.7, 0.904, 1.42, 0.926};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.13, 0.71, 1.2, 1.7, 0.546, 0.896, 0.747};

	//Critical parameters
	crit.rho = 3.495*170.02886; //[kg/m^3]
	crit.p = PressureUnit(2925, UNIT_KPA); //[kPa]
	crit.T = 374.9; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 170.02886;
	params.Ttriple = 146.35;
	params.accentricfactor = 0.35760476022557541;
	params.R_u = 8.314472;
	params.ptriple = 0.0073315669242655281;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,11,19));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 12,18,19));

	const double n5 = 0, n6 = 0, n0 = 4;
	phi0list.push_back(new phi0_lead(n5,n6));
	phi0list.push_back(new phi0_logtau(n0-1));

	const double u0[] = {0, 403/crit.T, 1428/crit.T};
	const double v0[] = {0, 11.43, 12.83};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,2));

	EOSReference.assign("Mark O. McLinden, Eric W. Lemmon, \" Thermodynamic Properties of R-227ea, R-365mfc, R-115, and R-13I1\" Preprint provided by Eric Lemmon");
	TransportReference.assign("Using ECS in fully predictive mode. Lennard-Jones parameters from NISTIR 6650");

	name.assign("R227EA");
	aliases.push_back(std::string("R227ea"));
	REFPROPname.assign("R227EA");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "McLinden-PREPRINT-2013";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
}
double R227EAClass::psat(double T)
{
    // Maximum absolute error is 0.014171 % between 146.350001 K and 374.899990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3,9};
    const double Ni[]={0,-7.7770183479122093, 1.9824076406295947, -2.6976862481795081, 0.098776424041554359, -7.4866660112540497, 9.792603493935939, -7.0950087986607366 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=7;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double R227EAClass::rhosatL(double T)
{
    // Maximum absolute error is 0.038984 % between 146.350001 K and 374.899990 K
    const double ti[]={0,0.33544749034418803, 0.99174918615685026, 2.1225866113354703, 2.2466438788850227, 2.4007651871441014, 2.5972300309334733, 2.7495964162340281, 5.1398348232050957};
    const double Ni[]={0,1.7639870592586082, -2.522453612659636, 1361.681410638186, -4122.6816874822398, 4805.8231801381262, -2982.2919993398209, 940.03615642523266, -0.46004659700788408};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=8;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double R227EAClass::rhosatV(double T)
{
    // Maximum absolute error is 0.252796 % between 146.350001 K and 374.899990 K
    const double ti[]={0,0.32727803601304156, 1.0663072855856532, 2.9689056156126843, 3.3281256431886144, 3.5152372361461248, 3.8121392440843458, 3.9132576223837896, 7.0803926901081651};
    const double Ni[]={0,-2.2028077488865478, -6.3911248425297638, 1542.6058303947407, -14554.384298585259, 25213.639486441465, -29654.132851429142, 17459.351007805486, -11.070145511065862};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=8;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}


R365MFCClass::R365MFCClass()
{
	double n[] = {0.0, 2.20027, -2.8624, 0.384559, -0.621227, 0.0665967, -1.19383, 0.635935, 0.461728, -0.533472, -1.07101, 0.13929, -0.385506, 0.885653, 0.226303, -0.166116};
	double t[] = {0, 0.24, 0.67, 0.5, 1.25, 1, 3.35, 2.5, 0.96, 1.07, 5.6, 6.9, 3, 3.6, 5, 1.25};
	double d[] = {0, 1, 1, 2, 2, 4, 1, 3, 6, 6, 2, 3, 1, 1, 1, 2};
	double c[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.97, 0.94, 2.15, 2.66};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.07, 1.08, 10.9, 22.6};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.48, 1.49, 1.01, 1.16};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.02, 0.62, 0.53, 0.48};

	//Critical parameters
	crit.rho = 3.2*148.07452 ; //[kg/m^3]
	crit.p = PressureUnit(3266, UNIT_KPA); //[kPa]
	crit.T = 460.0; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 148.07452 ;
	params.Ttriple = 239;
	params.accentricfactor = 0.37741704621544692;
	params.R_u = 8.314472;
	params.ptriple = 2.4784182018851801;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,11,16));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 12,15,16));

	const double n5 = 0, n6 = 0, n0 = 4;
	phi0list.push_back(new phi0_lead(n5, n6));
	phi0list.push_back(new phi0_logtau(n0-1));

	const double u0[] = {0, 569/crit.T, 2232/crit.T};
	const double v0[] = {0, 17.47, 16.29};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,2));

	EOSReference.assign("Mark O. McLinden, Eric W. Lemmon, \" Thermodynamic Properties of R-227ea, R-365mfc, R-115, and R-13I1\" Preprint provided by Eric Lemmon");
	TransportReference.assign("Using ECS in fully predictive mode. Lennard-Jones parameters from NISTIR 6650");

	name.assign("R365MFC");
	aliases.push_back(std::string("R365mfc"));
	REFPROPname.assign("R365MFC");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "McLinden-PREPRINT-2013";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R365MFCClass::psat(double T)
{
    // Fit from Eric Lemmon, April 8, 2013
    const double ti[]={0,1.0,1.5,3.4,4.3,5.0};
    const double Ni[]={0,-8.0955,2.0414,-13.3333,25.514,-19.967};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double R365MFCClass::rhosatL(double T)
{
    // Fit from Eric Lemmon, April 8, 2013
    const double ti[]={0,0.31,0.6,0.9,1.2,1.5};
    const double Ni[]={0,1.7933,-1.8792,9.0006,-11.669,5.6329};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*(1+summer);
}
double R365MFCClass::rhosatV(double T)
{
    // Fit from Eric Lemmon, April 8, 2013
    const double ti[]={0,0.281,0.91,3.0,5.0,8.0,10.0};
    const double Ni[]={0,-1.6120,-6.7679,-24.499,3.3398,-211.10,258.07};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
