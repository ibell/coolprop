/* 
Properties for R407C.  
by Ian Bell (ihb2@cornell.edu)

Pseudo-pure fluid thermo props from 
"Pseudo-pure fluid Equations of State for the Refrigerant Blends R410A, R407C, R507C and R407C" 
by E.W. Lemmon, Int. J. Thermophys. v. 24, n4, 2003

*/

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include "math.h"
#include "stdio.h"
#include <string.h>
#include "CoolProp.h"
#include "FluidClass.h"
#include "R407C.h"

R407CClass::R407CClass()
{
	static double a[]={
		2.13194,		//[0]
		8.05008,		//[1]
		-14.3914,		//[2]
		1.4245,			//[3]
		3.9419,			//[4]
		3.1209			//[5]
	};

	static double b[]={
		0,				//[0]
		0,				//[1]
		-0.4,			//[2]
		2.40437,		//[3]
		5.25122,		//[4]
		13.3632			//[5]
	};

	static double N[]={
		 0.0,			//[0]
		 1.0588,		//[1]
		-1.12018,		//[2]
		 0.629064,		//[3]
		-0.351953,		//[4]
		 0.00455978,	//[5]
		-1.75725,		//[6]
		-1.12009,		//[7]
		 0.0277353,		//[8]
		 0.898881,		//[9]
		-1.17591,		//[10]
		 0.0818591,		//[11]
		-0.0794097,		//[12]
		-0.0000104047,	//[13]
		 0.233779,		//[14]
		-0.291790,		//[15]
		 0.0154776,		//[16]
		-0.0314579,		//[17]
		-0.00442552,	//[18]
		-0.0101254,		//[19]
		 0.00915953,	//[20]
		-0.003615		//[21]
	};

	static double t[]={
		0.0,			//[0]
		0.241,			//[1]
		0.69,			//[2]
		2.58,			//[3]
		1.15,			//[4]
		0.248,			//[5]
		2.15,			//[6]
		2.43,			//[7]
		5.3,			//[8]
		0.76,			//[9]
		1.48,			//[10]
		0.24,			//[11]
		2.86,			//[12]
		8.0,			//[13]
		3.3,			//[14]
		4.7,			//[15]
		0.45,			//[16]
		8.4,			//[17]
		16.2,			//[18]
		26.0,			//[19]
		16.0,			//[20]
		8.7				//[21]
	};

	static double d[]={
		0,				//[0]
		1,				//[1]
		1,				//[2]
		1,				//[3]
		2,				//[4]
		5,				//[5]
		1,				//[6]
		2,				//[7]
		2,				//[8]
		3,				//[9]
		3,				//[10]
		5,				//[11]
		5,				//[12]
		5,				//[13]
		1,				//[14]
		1,				//[15]
		4,				//[16]
		4,				//[17]
		2,				//[18]
		4,				//[19]
		5,				//[20]
		6				//[21]
	};

	static double l[]={
		0,				//[0]
		0,				//[1]
		0,				//[2]
		0,				//[3]
		0,				//[4]
		0,				//[5]
		1,				//[6]
		1,				//[7]
		1,				//[8]
		1,				//[9]
		1,				//[10]
		1,				//[11]
		1,				//[12]
		1,				//[13]
		2,				//[14]
		2,				//[15]
		2,				//[16]
		2,				//[17]
		3,				//[18]
		3,				//[19]
		3,				//[20]
		3				//[21]
	};

	phirlist.push_back(new phir_power(N,d,t,l,1,21,22));

	/*
	sum=log(delta)-log(tau)+a[0]+a[1]*tau+a[2]*pow(tau,b[2]);
	for(k=3;k<=5;k++)
	{
		sum+=a[k]*log(1.0-exp(-b[k]*tau));
	}
	*/
	phi0list.push_back(new phi0_lead(a[0],a[1]));
	phi0list.push_back(new phi0_logtau(-1.0));
	phi0list.push_back(new phi0_power(a[2],b[2]));
	phi0list.push_back(new phi0_Planck_Einstein(a,b,3,5,6));

	// Critical parameters (max condensing temperature)
	crit.rho = 453.43094;
	crit.p = PressureUnit(4631.7,UNIT_KPA);
	crit.T = 359.345;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 86.2036;
	params.Ttriple = 200.0;
	params.ptriple = 11.3123436377;
	params.accentricfactor = 0.363;
	params.R_u = 8.314472;
	isPure = false;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 50000.0;
	limits.rhomax = 17.04*params.molemass;
	
	EOSReference.assign("E.W. Lemmon, \"Pseudo-pure fluid Equations of State for the Refrigerant Blends R410A, R404A, R507C and R407C\"" 
						",Int. J. Thermophys. v. 24, n4, 2003");
	TransportReference.assign("Viscosity: V. Geller, \"Viscosity of Mixed Refrigerants R404A,R407C,"
							"R410A, and R507A\", 2000 Purdue Refrigeration conferences\n\n"
							"Thermal Conductivity: V.Z. Geller, B.Z. Nemzer, and U.V. Cheremnykh \"Thermal Conductivity "
							"of the Refrigerant mixtures R404A,R407C,R410A, and R507A\" "
							"Int. J. Thermophysics, v. 22, n 4 2001\n\n"
							"Surface Tension: R. Heide, \"The surface tension of HFC refrigerants and mixtures\", Int J. Refrig. Vol. 20, No. 7, pp. 496-503, 1997");

	name.assign("R407C");
	aliases.push_back("R407c");

	BibTeXKeys.EOS = "Lemmon-IJT-2003";
	BibTeXKeys.VISCOSITY = "Geller-PURDUE-2000";
	BibTeXKeys.CONDUCTIVITY = "Geller-IJT-2001";
	BibTeXKeys.SURFACE_TENSION = "Heide-IJR-1997";
}

double R407CClass::psatL(double T)
{
	//Bubble point of R407C
	double sum=0,theta;
	int k;
	static const double Nbp[]={
		0.0,			//[0]
		0.48722,		//[1]
		-6.6959,		//[2]
		-1.4165,		//[3]
		-2.5109			//[4]
	};
	    
	static const double tbp[]={
		0.0,			//[0]
		0.54,			//[1]
		0.925,			//[2]
		2.7,			//[3]
		4.7				//[4]
	};

    theta=1-T/crit.T;
    
    for (k=1;k<=4;k++)
    {
        sum+=Nbp[k]*pow(theta,tbp[k]);
    }
    return reduce.p.Pa*exp(reduce.T/T*sum);
}
    
double R407CClass::psatV(double T)
{
	//Dew point of R407C
	double sum=0,theta;
	int k;
	static const double Ndp[]={
		0.0,			//[0]
		-0.086077,		//[1]
		-6.6364,		//[2]
		-2.4648,		//[3]
		-3.4776			//[4]
	};

	static const double tdp[]={
		0.0,			//[0]
		0.4,			//[1]
		0.965,			//[2]
		3.1,			//[3]
		5.0				//[4]
	};
    theta=1-T/crit.T;

    for (k=1;k<=4;k++)
    {
        sum+=Ndp[k]*pow(theta,tdp[k]);
    }
    return reduce.p.Pa*exp(reduce.T/T*sum);
}

double R407CClass::rhosatV(double T)
{
	double summer = 0, theta = 1-T/reduce.T;
	// Max error is  0.019852053576 % between 200.0 and 359.3449 K
    const double t[] = {0, 0.352, 0.096, 0.3555, 0.357, 0.358, 1.3333333333333333, 3.5, 7.5};
    const double N[] = {0, -179729.84916615338, 0.10946690748177662, 1453310.3859794352, -2554835.880001666, 1281252.1234214078, -8.1237911813327504, -25.223828823832683, -77.955885979137349};
	for (int i=1; i<=8; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(summer);
}

double R407CClass::rhosatL(double T)
{
	double theta;
	theta=1-T/359.345;
	return exp(5.544589249+1.764403419125e+00*pow(theta,0.12337)+0.544950396285*theta-0.784102758738*theta*theta+0.741332715649*theta*theta*theta);
}

double R407CClass::viscosity_Trho(double T, double rho)
{
    // Properties taken from "Viscosity of Mixed 
	// Refrigerants R407C,R407C,R410A, and R507A" 
	// by Vladimir Geller, 
	// 2000 Purdue Refrigeration conferences

    // inputs in T [K], and p [kPa]
    // output in Pa-s

   double eta_microPa_s;

   //Set constants required
   double a_0=-1.507e0,a_1=4.894e-2,a_2=-9.305e-6,b_1=-3.038e-3,b_2=2.927e-4,
	   b_3=-9.559e-7,b_4=1.739e-9,b_5=-1.455e-12,b_6=4.756e-16;

   eta_microPa_s=a_0+a_1*T+a_2*T*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho+b_5*rho*rho*rho*rho*rho+b_6*rho*rho*rho*rho*rho*rho;
   return eta_microPa_s/1e6;
}

double R407CClass::conductivity_Trho(double T, double rho)
{
	// Properties taken from "Thermal Conductivity 
	// of the Refrigerant mixtures R407C,R407C,R410A, and R507A" 
	// by V.Z. Geller, B.Z. Nemzer, and U.V. Cheremnykh 
	// Int. J. Thermophysics, v. 22, n 4 2001

	// inputs in T [K], and p [kPa] or rho [kg/m^3]
	// output in W/m-K

	//Set constants required
	double a_0=-9.628e0,a_1=7.638e-2,b_1=2.715e-2,b_2=4.963e-5,b_3=-4.912e-8,b_4=2.884e-11;

	return (a_0+a_1*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho)/1.e3; // from mW/m-K to W/m-K
}
double R407CClass::surface_tension_T(double T)
{
	return 0.0583286*pow(1-T/reduce.T,1.237);
}
