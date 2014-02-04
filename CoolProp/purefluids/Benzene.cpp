#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Benzene.h"

BenzeneClass::BenzeneClass()
{
	double n[] = {0.0, 0.03513062, 2.229707, -3.100459, -0.5763224, 0.2504179, -0.7049091, -0.1393433, 0.8319673, -0.3310741, -0.02793578, 0.7087408, -0.3723906, -0.06267414, -0.86295};
	double t[] = {0, 1, 0.3, 0.744, 1.174, 0.68, 2.5, 3.67, 1.26, 2.6, 0.95, 1, 2.47, 3.35, 0.75};
	double d[] = {0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3, 3};
	double c[] = {0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.032, 1.423, 1.071, 14.35};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.867, 1.766, 1.824, 297.5};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.118, 0.6392, 0.6536, 1.164};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.7289, 0.9074, 0.7655, 0.8711};

	//Critical parameters
	crit.rho = 3.902*78.1118; //[kg/m^3]
	crit.p = PressureUnit(4894,UNIT_KPA); //[kPa]
	crit.T = 562.02; //[K]
	crit.v = 1/crit.rho;

	// Other fluid parameters
	params.molemass = 78.1118;
	params.Ttriple = 278.674;
	params.accentricfactor = 0.21083697327001505;
	params.R_u = 8.314472;
	params.ptriple = 4.7837725790574392;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,10,15));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 11,14,15));

	const double a0 = -0.6740687105, a1 = 2.5560188958, n0 = 2.94645;
	phi0list.push_back(new phi0_lead(a0,a1));
	phi0list.push_back(new phi0_logtau(n0));

	const double u0[] = {0, 4116/crit.T, 1511/crit.T, 630/crit.T};
	const double v0[] = {0, 7.36374, 18.649, 4.01834};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,3));

	EOSReference.assign("Monika Thol and Eric W. Lemmon and Roland Span, \"Equation of state for benzene for temperatures from the melting line up to 725 K with pressures up to 500 MPa\", High Temperatures-High Pressures, Vol. 41, pp. 81�97");
	TransportReference.assign("");

	name.assign("Benzene");
	aliases.push_back(std::string("benzene"));
	aliases.push_back(std::string("BENZENE"));
	REFPROPname.assign("BENZENE");

	BibTeXKeys.EOS = "Thol-HTHP-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";
	BibTeXKeys.CONDUCTIVITY = "Assael-JPCRD-2012A";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double BenzeneClass::psat(double T)
{
    // Maximum absolute error is 0.007155 % between 278.674001 K and 562.019990 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3,9};
    const double Ni[]={0,-7.1461120764571957, 1.9957829398567826, -1.9350951186777388, -0.57944805848240821, -3.5619379943600871, -0.020339032269282332, 0.78993073753218157 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double BenzeneClass::rhosatL(double T)
{
    // Maximum absolute error is 0.271091 % between 278.674001 K and 562.019990 K
    const double ti[]={0,0.37662058801320508, 0.73750586561787868, 1.496536312566036, 2.719354797503057, 5.0887565894044391, 4.5762119437433926};
    const double Ni[]={0,2.1287356475388064, -1.3313426148430443, 0.79624913133294628, -0.37980073421164662, 0.25973635751163832, 0.0031066997480834807};
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
double BenzeneClass::rhosatV(double T)
{
    // Maximum absolute error is 0.212877 % between 278.674001 K and 562.019990 K
    const double ti[]={0,0.4025245240893699, 1.2389222941368063, 1.7971910824209907, 1.9285664588450286, 2.9793432570646234, 2.6542255328369846};
    const double Ni[]={0,-2.8057874888404788, -16.379051852677751, 218.63205995702904, -258.28890459277318, -81.624634476247607, 129.94924241369819};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
double BenzeneClass::conductivity_Trho(double T, double rho)
{
	// Assael JPCRD 2012
	double sumresid = 0;
	double B1[] = {0, 2.82489e-2, -7.73415e-2, 7.14001e-2, -2.36798e-2, 3.00875e-3};
	double B2[] = {0, -1.19268e-2, 8.33389e-2, -8.98176e-2, 3.63025e-2, -4.90052e-3};

	// in mW/m/K in paper, converted to W/m/K
	double lambda_0 = (101.404-521.440*T/crit.T+868.266*pow(T/crit.T,2))/(1+9.714*T/crit.T+1.467*pow(T/crit.T,2))/1000; // [W/m/K]

	for (int i = 1; i <= 5; i++)
	{
		sumresid += (B1[i]+B2[i]*T/crit.T)*pow(rho/reduce.rho,i); // [W/m/K]
	}
	double lambda_r = sumresid;

	double lambda_c = this->conductivity_critical(T,rho,1/(6.2e-10)); //[W/m/K]

	return lambda_0+lambda_r+lambda_c;
}
