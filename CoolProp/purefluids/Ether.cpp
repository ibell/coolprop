#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Ether.h"
#include "REFPROP.h"

DimethylEtherClass::DimethylEtherClass()
{
	double n[] = {0.0, 0.029814139, 1.43517, -2.64964, -0.29515532, 0.17035607, -0.94642918, -0.099250514, 1.1264071, -0.76936548, -0.020717696, 0.24527037, 1.1863438, -0.49398368, -0.16388716, -0.027583584};
	double t[] = {0, 1, 0.4366, 1.011, 1.137, 0.45, 2.83, 1.5, 1.235, 2.675, 0.7272, 1.816, 1.783, 3.779, 3.282, 1.059};
	double d[] = {0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 1, 3, 3};
	double c[] = {0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 1, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.965336, 1.50858, 0.963855, 9.72643};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.28719, 0.806235, 0.777942, 197.681};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.27772, 0.43075, 0.429607, 1.13849};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.672698, 0.924246, 0.750815, 0.800022};

	//Critical parameters
	crit.rho = 5.940*46.06844; //[kg/m^3]
	crit.p = PressureUnit(5336.8, UNIT_KPA); //[kPa]
	crit.T = 400.378; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 46.06844;
	params.Ttriple = 131.66;
	params.accentricfactor = 0.196;
	params.R_u = 8.314472;
	params.ptriple = 0.0022;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,11,16));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 12,15,16));

	const double a1 = -1.980976, a2 = 3.171218, c0 = 4.039;
	phi0list.push_back(new phi0_lead(a1,a2));
	phi0list.push_back(new phi0_logtau(c0-1));

	const double u0[] = {0, 361/crit.T, 974/crit.T, 1916/crit.T, 4150/crit.T};
	const double v0[] = {0, 2.641, 2.123, 8.992, 6.191};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	EOSReference.assign("Jiangtao Wu, Yong Zhou, Eric W. Lemmon, \"An Equation of State for the Thermodynamic Properties of Dimethyl Ether\", J. Phys. Chem. Ref. Data, Vol. 40, No. 2, 2011");
	TransportReference.assign("Viscosity: Xianyang Meng, Jianbo Zhang, Jiangtao Wu, and Zhigang Liu, \" Experimental Measurement and Modeling of the Viscosity of Dimethyl Ether\" J. Chem. Eng.Data 2012, 57, 988-993\n\nErratum: Limits for deltaeta_r sums should be 0-3 and 4-9.  The correct order of terms based on the original indices are 0,1,7,9,2,3,4,5,6,8\n\n"
		"Using ECS in fully predictive mode for viscosity\n\n"
		"Lennard-Jones parameters from Chichester NISTIR 6650");
	name.assign("DimethylEther");
	aliases.push_back(std::string("DIMETHYLETHER"));
	REFPROPname.assign("DME");

	BibTeXKeys.EOS = "Wu-JPCRD-2011";
	BibTeXKeys.VISCOSITY = "Meng-JCED-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Chichester-NIST-2008";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double DimethylEtherClass::psat(double T)
{
    const double ti[]={0,1.0,1.5,2.5,5};
    const double Ni[]={0, -7.112782, 1.971239, -2.276083, -2.215774};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double DimethylEtherClass::rhosatL(double T)
{
    const double ti[]={0,0.54,0.74,0.95,11.43};
    const double Ni[]={0, 7.884834, -10.516328, 5.39142, 0.404890};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*(summer+1);
}
double DimethylEtherClass::rhosatV(double T)
{
    // Maximum absolute error is 0.161887 % between 87.800001 K and 419.289990 K
    const double ti[]={0,1.467/3.0,4.2/3.0,8.0/3.0,17.0/3.0,36.0/3.0};
    const double Ni[]={0, -4.136444, -4.302025, -12.032140, -39.527936, -89.476860};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}

double DimethylEtherClass::conductivity_Trho(double T, double rho)
{
	long iPropane = get_Fluid_index(std::string("Propane"));
	// Calculate the ECS
	double lambda = conductivity_ECS_Trho(T, rho, get_fluid(iPropane));
	return lambda;
}
//double DimethylEtherClass::viscosity_Trho(double T, double rho)
//{
//	long iPropane = get_Fluid_index(std::string("Propane"));
//	// Calculate the ECS
//	double mu = viscosity_ECS_Trho(T, rho, get_fluid(iPropane));
//	return mu;
//}
void DimethylEtherClass::ECSParams(double *e_k, double *sigma)
{
	*e_k = 395; *sigma =0.4307;
}
double DimethylEtherClass::viscosity_Trho(double T, double rho)
{
	//return REFPROP(std::string("V"),std::string("T"),T,std::string("D"),rho,std::string("REFPROP-DME"));

	double sigma = 0.446704; //[nm]
	double M = 46.06844; //[kg/kmol]
	double tau = T/reduce.T;
	double delta = rho/reduce.rho;
	double Tstar = T/317.937;
	double log_Tstar = log(Tstar);
	double log_theta_star = 0.294261-0.377826*log_Tstar-0.491673*log_Tstar*log_Tstar;
	double eta_0 = 0.021375*sqrt(M*T)/(sigma*sigma*exp(log_theta_star)); //[uPa-s]
	
	double n[] = {-2.70002,4.44583,0.21302,6.50681,-104.998,78.27474,41.3751,-175.055,62.81975,112.3219};
	double t[] = {-5.92,-4.36,-5.87,-0.45,-2.93,-1.64,-7.86,-4.25,-4.79,-3.11};
	double d[] = {3,3,5,1,3,4,5,2,2,2};
	double p[] = {0,0,0,0,1,1,2,1,1,2};

	double eta_r = n[0]*pow(tau,t[0])*pow(delta,d[0])+n[1]*pow(tau,t[1])*pow(delta,d[1])+n[2]*pow(tau,t[2])*pow(delta,d[2])+n[3]*pow(tau,t[3])*pow(delta,d[3]);
	for (unsigned int i = 4; i<=9; i++)
	{
		eta_r += n[i]*pow(tau,t[i])*pow(delta,d[i])*exp(-pow(delta,p[i]));
	}
	return (eta_r+eta_0)/1e6;
}

double DimethylEtherClass::surface_tension_T(double T)
{
	return 0.063157*pow(1-T/reduce.T,1.2595);
}
