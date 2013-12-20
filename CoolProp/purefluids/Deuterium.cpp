/* 
Properties of Normal hydrogen
by Ian Bell
*/

#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include "math.h"
#include "stdio.h"
#include <string.h>
#include "CoolProp.h"
#include "FluidClass.h"
#include "Deuterium.h"

DeuteriumClass::DeuteriumClass()
{
	double n[] = {0, 0.006267958, 10.53609, -10.14149, 0.356061, 0.1824472, -1.129638, -0.0549812, -0.6791329, 1.347918, -0.8657582, 1.719146, -1.917977, 0.1233365, -0.07936891, 1.686617, -4.240326, 1.857114, -0.5903705, 1.520171, 2.361373, -2.297315};
	double t[] = {0, 1, 0.462, 0.5584, 0.627, 1.201, 0.309, 1.314, 1.1166, 1.25, 1.25, 1.395, 1.627, 1, 2.5, 0.635, 0.664, 0.7082, 2.25, 1.524, 0.67, 0.709};
	double d[] = {0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 2, 1, 1, 3, 2, 1, 1, 2, 3, 3, 1, 3};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0};
	double alpha[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.868, 0.636, 0.668, 0.65, 0.745, 0.782, 0.693}; // phi from paper
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.46, 1.7864, 1.647, 0.541, 0.969, 1.892, 1.076}; // D from paper
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.613, 0.584, 0.57, 1.056, 1.01, 1.025, 1.029};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6306, 0.711, 0.6446, 0.8226, 0.992, 1.2184, 1.203};

	//Constants for ideal gas expression
	double a0[] = {0, -3.54145, 3.0326, -3.52422, -1.73421, -3.57135, 2.14858, 6.23107, -3.30425, 6.23098, -3.57137, 3.32901, 0.97782};
	double n0[] = {0, 7174.1/38.34, 8635/38.34, 902.7/38.34, 181.1/38.34, 438.5/38.34, 5034.2/38.34, 269.9/38.34, 229.9/38.34, 666.4/38.34, 452.8/38.34, 192/38.34, 1187.6/38.34};

	phirlist.push_back(new phir_power(n,d,t,c,1,14,22));
	phirlist.push_back(new phir_gaussian(n,d,t,alpha,epsilon,beta,gamma,15,21,22));

	//lead term of the form log(delta)+a1+a2*tau
	phi0list.push_back(new phi0_lead(-2.0677351753,2.4237151502));
	phi0list.push_back(new phi0_logtau(1.5));
	phi0list.push_back(new phi0_Planck_Einstein(a0,n0,1,12,13));

	// Critical parameters
	crit.rho = 17.23*4.0282;
	crit.p = PressureUnit(1679.6, UNIT_KPA);
	crit.T = 38.34;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 4.0282;
	params.Ttriple = 18.724;
	params.ptriple = 17.202;
	params.accentricfactor = -0.136290274128;
	params.R_u = 8.3144621;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 1000.0;
	limits.pmax = 2000000.0;
	limits.rhomax = 102.0*params.molemass;

	name.assign("Deuterium");
	aliases.push_back("deuterium");
	aliases.push_back("DEUTERIUM");
	aliases.push_back("D2");
	REFPROPname.assign("D2");

	BibTeXKeys.EOS = "Richardson-JPCRD-2013";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double DeuteriumClass::psat(double T)
{
	const double ti[]={0,1.0,1.5,2.83,4.06,5.4};
    const double Ni[]={0,-5.5706,1.7631,-0.5458,1.2154,-1.1556};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double DeuteriumClass::rhosatL(double T)
{
	const double ti[]={0,0.512,1.12,1.8,2.55,3.4,4.4};
    const double Ni[]={0,3.3769,-5.3693,11.943,-17.361,15.170,-6.3079};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*(1+summer);
}
double DeuteriumClass::rhosatV(double T)
{
	const double ti[]={0,0.528,2.03,3.6,5,6.5,9};
    const double Ni[]={0,-3.8111,-7.3624,2.2294,-21.443,12.796,-31.334};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double DeuteriumClass::surface_tension_T(double T)
{
	// From Mulero, JPCRD, 2012
	return 0.009376*pow(1-T/reduce.T,1.258);
}

OrthoDeuteriumClass::OrthoDeuteriumClass()
{
	double n[] = {0, 0.006267958, 10.53609, -10.14149, 0.356061, 0.1824472, -1.129638, -0.0549812, -0.6791329, 1.347918, -0.8657582, 1.719146, -1.917977, 0.1233365, -0.07936891, 1.686617, -4.240326, 1.857114, -0.5903705, 1.520171, 2.361373, -2.297315};
	double t[] = {0, 1, 0.462, 0.5584, 0.627, 1.201, 0.309, 1.314, 1.1166, 1.25, 1.25, 1.395, 1.627, 1, 2.5, 0.635, 0.664, 0.7082, 2.25, 1.524, 0.67, 0.709};
	double d[] = {0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 2, 1, 1, 3, 2, 1, 1, 2, 3, 3, 1, 3};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0};
	double alpha[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.868, 0.636, 0.668, 0.65, 0.745, 0.782, 0.693}; // phi from paper
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.46, 1.7864, 1.647, 0.541, 0.969, 1.892, 1.076}; // D from paper
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.613, 0.584, 0.57, 1.056, 1.01, 1.025, 1.029};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6306, 0.711, 0.6446, 0.8226, 0.992, 1.2184, 1.203};

	//Constants for ideal gas expression
	double a0[] = {0, 4.04482, -4.65391, -4.65342, 3.46313, -4.58637, -4.6503, -4.65124, 2.67024, 15.20455, 0.87164, -4.7608, 4.32447};
	double n0[] = {0, 1591/38.34, 481.6/38.34, 472.4/38.34, 362.2/38.34, 2038/38.34, 463.2/38.34, 491.3/38.34, 2713.4/38.34, 618.6/38.34, 8642/38.34, 961.7/38.34, 253.2/38.34};

	phirlist.push_back(new phir_power(n,d,t,c,1,14,22));
	phirlist.push_back(new phir_gaussian(n,d,t,alpha,epsilon,beta,gamma,15,21,22));

	//lead term of the form log(delta)+a1+a2*tau
	phi0list.push_back(new phi0_lead(-2.0672670563, 2.4234599781));
	phi0list.push_back(new phi0_logtau(1.5));
	phi0list.push_back(new phi0_Planck_Einstein(a0,n0,1,12,13));

	// Critical parameters
	crit.rho = 17.23*4.0282;
	crit.p = PressureUnit(1679.6, UNIT_KPA);
	crit.T = 38.34;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 4.0282;
	params.Ttriple = 18.724;
	params.ptriple = 17.202;
	params.accentricfactor = -0.136290274128;
	params.R_u = 8.3144621;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 1000.0;
	limits.pmax = 2000000.0;
	limits.rhomax = 102.0*params.molemass;

	name.assign("OrthoDeuterium");
	aliases.push_back("orthodeuterium");
	aliases.push_back("ORTHODEUTERIUM");
	REFPROPname = "N/A";

	BibTeXKeys.EOS = "Richardson-JPCRD-2013";
}

double OrthoDeuteriumClass::psat(double T)
{
	const double ti[]={0,1.0,1.5,2.83,4.06,5.4};
    const double Ni[]={0,-5.5706,1.7631,-0.5458,1.2154,-1.1556};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double OrthoDeuteriumClass::rhosatL(double T)
{
	const double ti[]={0,0.512,1.12,1.8,2.55,3.4,4.4};
    const double Ni[]={0,3.3769,-5.3693,11.943,-17.361,15.170,-6.3079};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*(1+summer);
}
double OrthoDeuteriumClass::rhosatV(double T)
{
	const double ti[]={0,0.528,2.03,3.6,5,6.5,9};
    const double Ni[]={0,-3.8111,-7.3624,2.2294,-21.443,12.796,-31.334};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}

ParaDeuteriumClass::ParaDeuteriumClass()
{
	double n[] = {0, 0.006267958, 10.53609, -10.14149, 0.356061, 0.1824472, -1.129638, -0.0549812, -0.6791329, 1.347918, -0.8657582, 1.719146, -1.917977, 0.1233365, -0.07936891, 1.686617, -4.240326, 1.857114, -0.5903705, 1.520171, 2.361373, -2.297315};
	double t[] = {0, 1, 0.462, 0.5584, 0.627, 1.201, 0.309, 1.314, 1.1166, 1.25, 1.25, 1.395, 1.627, 1, 2.5, 0.635, 0.664, 0.7082, 2.25, 1.524, 0.67, 0.709};
	double d[] = {0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 2, 1, 1, 3, 2, 1, 1, 2, 3, 3, 1, 3};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0};
	double alpha[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.868, 0.636, 0.668, 0.65, 0.745, 0.782, 0.693}; // phi from paper
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.46, 1.7864, 1.647, 0.541, 0.969, 1.892, 1.076}; // D from paper
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.613, 0.584, 0.57, 1.056, 1.01, 1.025, 1.029};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6306, 0.711, 0.6446, 0.8226, 0.992, 1.2184, 1.203};

	//Constants for ideal gas expression
	double a0[] = {0, 1.28527, 1.11376, -2.491, 6.38763, 6.17406, -3.13698, -3.14254, -2.29511, -3.37, 1.13634, 0.72512};
	double n0[] = {0, 5068/38.34, 1000.8/38.34, 261.5/38.34, 437.2/38.34, 312.3/38.34, 382.8/38.34, 356.8/38.34, 294.7/38.34, 682.4/38.34, 246/38.34, 277.1/38.34};

	phirlist.push_back(new phir_power(n,d,t,c,1,14,22));
	phirlist.push_back(new phir_gaussian(n,d,t,alpha,epsilon,beta,gamma,15,21,22));

	//lead term of the form log(delta)+a1+a2*tau
	phi0list.push_back(new phi0_lead(-2.0683998716, 2.4241000701));
	phi0list.push_back(new phi0_logtau(1.5));
	phi0list.push_back(new phi0_Planck_Einstein(a0,n0,1,11,12));

	// Critical parameters
	crit.rho = 17.23*4.0282;
	crit.p = PressureUnit(1679.6, UNIT_KPA);
	crit.T = 38.34;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 4.0282;
	params.Ttriple = 18.724;
	params.ptriple = 17.202;
	params.accentricfactor = -0.136290274128;
	params.R_u = 8.3144621;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 1000.0;
	limits.pmax = 2000000.0;
	limits.rhomax = 102.0*params.molemass;

	name.assign("ParaDeuterium");
	aliases.push_back("paradeuterium");
	aliases.push_back("PARADEUTERIUM");
	REFPROPname = "N/A";

	BibTeXKeys.EOS = "Richardson-JPCRD-2013";
}

double ParaDeuteriumClass::psat(double T)
{
	const double ti[]={0,1.0,1.5,2.83,4.06,5.4};
    const double Ni[]={0,-5.5706,1.7631,-0.5458,1.2154,-1.1556};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double ParaDeuteriumClass::rhosatL(double T)
{
	const double ti[]={0,0.512,1.12,1.8,2.55,3.4,4.4};
    const double Ni[]={0,3.3769,-5.3693,11.943,-17.361,15.170,-6.3079};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*(1+summer);
}
double ParaDeuteriumClass::rhosatV(double T)
{
	const double ti[]={0,0.528,2.03,3.6,5,6.5,9};
    const double Ni[]={0,-3.8111,-7.3624,2.2294,-21.443,12.796,-31.334};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
