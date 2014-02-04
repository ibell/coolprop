/*
Properties for R1234yf
by Ian Bell

Thermo props from
M. Richter and M.O. McLinden and E.W. Lemmon, "Thermodynamic Properties of 2,3,3,3-Tetrafluoroprop-1-ene
(R1234yf): Vapor Pressure and p-rho-T Measurements and an Equation of State"
by Reiner Tillner-Roth and Hans Dieter Baehr, J. Chem. Eng. Data, v. 56, 2011, pp 3254-3264
*/

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include <math.h>
#include "string.h"
#include "stdio.h"
#include "CoolProp.h"
#include "FluidClass.h"
#include "R1234yf.h"

R1234yfClass::R1234yfClass()
{
	static const double n[]={
		 0, //[0]
		 0.04592563, //[1]
		 1.546958, //[2]
		-2.355237, //[3]
		-0.4827835, //[4]
		 0.1758022, //[5]
		-1.210006, //[6]
		-0.6177084, //[7]
		 0.6805262, //[8]
		-0.6968555, //[9]
		-0.02695779, //[10]
		 1.389966, //[11]
		-0.4777136, //[12]
		-0.1975184, //[13]
		-1.147646, //[14]
		 0.0003428541 //[15]
	};

	static const double d[]={
		0, //[0]
		4, //[1]
		1, //[2]
		1, //[3]
		2, //[4]
		3, //[5]
		1, //[6]
		3, //[7]
		2, //[8]
		2, //[9]
		7, //[10]
		1, //[11]
		1, //[12]
		3, //[13]
		3, //[14]
		2, //[15]
	};

	static const double t[]={
		0, //[0]
		1, //[1]
		0.32, //[2]
		0.929, //[3]
		0.94, //[4]
		0.38, //[5]
		2.28, //[6]
		1.76, //[7]
		0.97, //[8]
		2.44, //[9]
		1.05, //[10]
		1.4, //[11]
		3, //[12]
		3.5, //[13]
		1, //[14]
		3.5 //[15]
	};

	static const double c[]={
	0, //[0]
	0, //[1]
	0, //[2]
	0, //[3]
	0, //[4]
	0, //[5]
	2, //[6]
	2, //[7]
	1, //[8]
	2, //[9]
	1, //[10]
	};

	static const double alpha[]={
	0, //[0]
	0, //[1]
	0, //[2]
	0, //[3]
	0, //[4]
	0, //[5]
	0, //[6]
	0, //[7]
	0, //[8]
	0, //[9]
	0, //[10]
	1.02, //[11]
	1.336, //[12]
	1.055, //[13]
	5.84, //[14]
	16.2, //[15]
	};

	static const double beta[]={
	0, //[0]
	0, //[1]
	0, //[2]
	0, //[3]
	0, //[4]
	0, //[5]
	0, //[6]
	0, //[7]
	0, //[8]
	0, //[9]
	0, //[10]
	1.42, //[11]
	2.31, //[12]
	0.89, //[13]
	80, //[14]
	108, //[15]
	};

	static const double gamma[]={
	0, //[0]
	0, //[1]
	0, //[2]
	0, //[3]
	0, //[4]
	0, //[5]
	0, //[6]
	0, //[7]
	0, //[8]
	0, //[9]
	0, //[10]
	1.13, //[11]
	0.67, //[12]
	0.46, //[13]
	1.28, //[14]
	1.2 //[15]
	};

	static const double epsilon[]={
	0, //[0]
	0, //[1]
	0, //[2]
	0, //[3]
	0, //[4]
	0, //[5]
	0, //[6]
	0, //[7]
	0, //[8]
	0, //[9]
	0, //[10]
	0.712, //[11]
	0.91, //[12]
	0.677, //[13]
	0.718, //[14]
	1.64 //[15]
	};

	static const double v0[]={
		0.0,	//[0]
		7.549,	//[1]
		1.537,	//[2]
		2.030,	//[3]
		7.455,	//[4]
	};
	static const double u0[]={
		// each of the ui terms are divided by the critical temp {fie on you Mr. EOS-maker}
		0.0,		//[0]
		718/367.85,	//[1]
		877/367.85,	//[2]
		4465/367.85,//[3]
		1755/367.85,//[4]
	};

	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));
	
	phirlist.push_back(new phir_power(n,d,t,c,1,10,11));
    phirlist.push_back(new phir_gaussian(n,d,t,alpha,epsilon,beta,gamma,11,15,16));

	phi0list.push_back(new phi0_lead(-12.837928,8.042605));
	phi0list.push_back(new phi0_logtau(4.944));
	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	// Critical parameters
	crit.rho = 475.553441976;
	crit.p = PressureUnit(3382.2, UNIT_KPA);
	crit.T = 367.85;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 114.0415928;
	params.Ttriple = 220;
	params.ptriple = 31.5093083629;
	params.accentricfactor = 0.276;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = 220;
	limits.Tmax = 410.0;
	limits.pmax = 30000.0;
	limits.rhomax = 11.64*params.molemass;
	
	EOSReference.assign("Richter, M. and M.O. McLinden and E.W. Lemmon, \"Thermodynamic Properties of 2,3,3,3-Tetrafluoroprop-1-ene"
						"(R1234yf): Vapor Pressure and p-rho-T Measurements and an Equation of State\""
						", J. Chem. Eng. Data, v. 56, 2011, pp 3254-3264");
	TransportReference.assign("Surface Tension: Katsuyuki Tanaka, Yukihiro Higashi, \"Thermodynamic properties of HFO-1234yf (2,3,3,3-tetrafluoropropene)\", International Journal of Refrigeration 33 (2010) 474-479");

	name.assign("R1234yf");
	aliases.push_back("R1234YF");

	BibTeXKeys.EOS = "Richter-JCED-2011";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.CONDUCTIVITY = "Perkins-JCED-2011";

}
double R1234yfClass::psat(double T)
{
    // Maximum absolute error is 0.012005 % between 220.000001 K and 367.849999 K
    const double ti[]={0,1.0,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.4299777557954254, 1.9583865741838868, -2.4230509106877167, 1.4167726905491065, -9.6120250676570098, 15.003017622880735 };
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return crit.p.Pa*exp(crit.T/T*summer);
}
double R1234yfClass::rhosatL(double T)
{
    // Maximum absolute error is 0.320186 % between 220.000001 K and 367.849999 K
    const double Ni[]={0,1.4756519570309901, -0.57454181759520018, 0.76682587551792647, -0.510882068727329, 0.66112133039760046};
    const double ti[]={0,0.31104118625261612, 1.1997863650360587, 1.8936812371106437, 3.0633279469212082, 6.4584770277609298};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return crit.rho*exp(summer);
}
double R1234yfClass::rhosatV(double T)
{
    // Maximum absolute error is 0.189556 % between 220.000001 K and 367.849999 K
    const double Ni[]={0,-2.4937203841988023, -4.0958221912377883, 2.598118460659435, -4.9686143098130451, -4.0779557647493059};
    const double ti[]={0,0.36409300811029677, 0.99593569215609024, 1.8193545503400304, 3.0060054082610717, 7.955467733338411};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return crit.rho*exp(crit.T/T*summer);
}
double R1234yfClass::surface_tension_T(double T)
{
	// From Mulero, 2012, JPCRD
	return 0.06274*pow(1-T/reduce.T,1.394);
}
double R1234yfClass::conductivity_Trho(double T, double rho)
{
	double sumresid = 0;
	double A[] = {-0.0102778, 0.0291098, 0.000860643};
	double B1[] = {0, -0.0368219, 0.0883226, -0.0705909, 0.0259026, -0.00322950};
	double B2[] = {0, 0.0397166, -0.0772390, 0.0664707, -0.0249071, 0.00336228};

	double lambda_0 = A[0]*pow(T/reduce.T,0)+A[1]*pow(T/reduce.T,1)+A[2]*pow(T/reduce.T,2);

	for (int i = 1; i <= 5; i++)
	{
		sumresid += (B1[i]+B2[i]*T/reduce.T)*pow(rho/reduce.rho,i);
	}
	double lambda_r = sumresid;

	double lambda_c = this->conductivity_critical(T,rho,1/(0.585e-9));

	return lambda_0+lambda_r+lambda_c;
}
