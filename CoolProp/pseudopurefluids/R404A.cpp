/*
Properties for R404A.  
by Ian Bell

Pseudo-pure fluid thermo props from 
"Pseudo-pure fluid Equations of State for the Refrigerant Blends R410A, R404A, R507C and R407C" 
by E.W. Lemmon, Int. J. Thermophys. v. 24, n4, 2003

In order to call the exposed functions, rho_, h_, s_, cp_,...... 
there are three different ways the inputs can be passed, and this is expressed by the Types integer flag.  
These macros are defined in the PropMacros.h header file:
1) First parameter temperature, second parameter pressure ex: h_R404A(260,354.7,1)=274
	In this case, the lookup tables are built if needed and then interpolated
2) First parameter temperature, second parameter density ex: h_R404A(260,13.03,2)=274
	Density and temp plugged directly into EOS
3) First parameter temperature, second parameter pressure ex: h_R404A(260,354.7,3)=274
	Density solved for, then plugged into EOS (can be quite slow)

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
#include "R404A.h"

R404AClass::R404AClass()
{
	static const double a[]={
		7.00407,		//[0]
		7.98695,		//[1]
		-18.8664,		//[2]
		0.63078,		//[3]
		3.5979,			//[4]
		5.0335			//[5]
	};

	static const double b[]={
		0,				//[0]
		0,				//[1]
		-0.3,			//[2]
		1.19617,		//[3]
		2.32861,		//[4]
		5.00188			//[5]
	};

	static const double N[]={
		 0.0,			//[0]
		 6.10984,		//[1]
		-7.79453,		//[2]
		 0.0183377,		//[3]
		 0.262270,		//[4]
		-0.00351688,	//[5]
		 0.0116181,		//[6]
		 0.00105992,	//[7]
		 0.850922,		//[8]
		-0.520084,		//[9]
		-0.0464225,		//[10]
		 0.621190,		//[11]
		-0.195505,		//[12]
		 0.336159,		//[13]
		-0.0376062,		//[14]
		-0.00636579,	//[15]
		-0.0758262,		//[16]
		-0.0221041,		//[17]
		 0.0310441,		//[18]
		 0.0132798,		//[19]
		 0.0689437,		//[20]
		-0.0507525,		//[21]
		 0.0161382,		//[22]
	};

	static const double t[]={
		0.0,			//[0]
		0.67,			//[1]
		0.91,			//[2]
		5.96,			//[3]
		0.7,			//[4]
		6.0,			//[5]
		0.3,			//[6]
		0.7,			//[7]
		1.7,			//[8]
		3.3,			//[9]
		7.0,			//[10]
		2.05,			//[11]
		4.3,			//[12]
		2.7,			//[13]
		1.8,			//[14]
		1.25,			//[15]
		12.0,			//[16]
		6.0,			//[17]
		8.7,			//[18]
		11.6,			//[19]
		13.0,			//[20]
		17.0,			//[21]
		16.0			//[22]
	};

	static const int d[]={
		0,				//[0]
		1,				//[1]
		1,				//[2]
		1,				//[3]
		2,				//[4]
		2,				//[5]
		4,				//[6]
		6,				//[7]
		1,				//[8]
		1,				//[9]
		1,				//[10]
		2,				//[11]
		2,				//[12]
		3,				//[13]
		4,				//[14]
		7,				//[15]
		2,				//[16]
		3,				//[17]
		4,				//[18]
		4,				//[19]
		2,				//[20]
		3,				//[21]
		5				//[22]
	};

	static const int l[]={
		0,				//[0]
		0,				//[1]
		0,				//[2]
		0,				//[3]
		0,				//[4]
		0,				//[5]
		0,				//[6]
		0,				//[7]
		1,				//[8]
		1,				//[9]
		1,				//[10]
		1,				//[11]
		1,				//[12]
		1,				//[13]
		1,				//[14]
		1,				//[15]
		2,				//[16]
		2,				//[17]
		2,				//[18]
		2,				//[19]
		3,				//[20]
		3,				//[21]
		3				//[22]
	};

	std::vector<double> n_v(N,N+sizeof(N)/sizeof(double));
	std::vector<double> d_v(d,d+sizeof(d)/sizeof(int));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(l,l+sizeof(l)/sizeof(int));
	std::vector<double> a0(a,a+sizeof(a)/sizeof(double));
	std::vector<double> n0(b,b+sizeof(b)/sizeof(double));

	phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,22);
	phirlist.push_back(phir_);

	/*
	sum=log(delta)-log(tau)+a[0]+a[1]*tau+a[2]*pow(tau,b[2]);
	for(k=3;k<=5;k++)
	{
		sum+=a[k]*log(1.0-exp(-b[k]*tau));
	}
	*/
	phi_BC * phi0_lead_ = new phi0_lead(a0[0],a0[1]);
	phi_BC * phi0_logtau_ = new phi0_logtau(-1.0);
	phi_BC * phi0_power_ = new phi0_power(a0[2],n0[2]);
	phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(a0,n0,3,5);

	phi0list.push_back(phi0_lead_);
	phi0list.push_back(phi0_logtau_);
	phi0list.push_back(phi0_power_);
	phi0list.push_back(phi0_Planck_Einstein_);

	// Critical parameters (max condensing temperature)
	crit.rho = 482.162772;
	crit.p = PressureUnit(3734.8,UNIT_KPA);
	crit.T = 345.27;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 97.6038;
	params.Ttriple = 200.0;
	params.ptriple = 21.2656766151;
	params.accentricfactor = 0.293;
	params.R_u = 8.314472;
	isPure = false;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 50000.0;
	limits.rhomax = 14.21*params.molemass;
	
	EOSReference.assign("E.W. Lemmon, \"Pseudo-pure fluid Equations of State for the Refrigerant Blends R410A, R404A, R507C and R407C\"" 
						",Int. J. Thermophys. v. 24, n4, 2003");
	TransportReference.assign("Viscosity: V. Geller, \"Viscosity of Mixed Refrigerants R404A,R407C,"
							"R410A, and R507A\", 2000 Purdue Refrigeration conferences\n\n"
							"Thermal Conductivity: V.Z. Geller, B.Z. Nemzer, and U.V. Cheremnykh \"Thermal Conductivity "
							"of the Refrigerant mixtures R404A,R407C,R410A, and R507A\" "
							"Int. J. Thermophysics, v. 22, n 4 2001\n\n"
							"Surface Tension: R. Heide, \"The surface tension of HFC refrigerants and mixtures\", Int J. Refrig. Vol. 20, No. 7, pp. 496-503, 1997");

	name.assign("R404A");
	aliases.push_back("R404a");

	BibTeXKeys.EOS = "Lemmon-IJT-2003";
	BibTeXKeys.VISCOSITY = "Geller-PURDUE-2000";
	BibTeXKeys.CONDUCTIVITY = "Geller-IJT-2001";
	BibTeXKeys.SURFACE_TENSION ="Heide-IJR-1997";
}

double R404AClass::conductivity_Trho(double T, double rho)
{
	// Properties taken from "Thermal Conductivity 
	// of the Refrigerant mixtures R404A,R407C,R410A, and R507A" 
	// by V.Z. Geller, B.Z. Nemzer, and U.V. Cheremnykh 
	// Int. J. Thermophysics, v. 22, n 4 2001

	// inputs in T [K], and p [kPa] or rho [kg/m^3]
	// output in W/m-K

	//Set constants required
	double a_0=-8.624e0,a_1=7.360e-2,b_1=3.222e-2,b_2=2.569e-5,b_3=-2.693e-8,b_4=2.007e-11;

	return (a_0+a_1*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho)/1.e3; // from mW/m-K to kW/m-K
}

double R404AClass::psatL(double T)
{
	//Bubble point of R410A
	double sum=0,theta;
	int k;

	static const double Nbp[]={
		0.0,			//[0]
		0.061067,		//[1]
		-6.5646,		//[2]
		-3.6162,		//[3]
		-3.9771			//[4]
	};
	    
	static const double tbp[]={
		0.0,			//[0]
		0.54,			//[1]
		0.965,			//[2]
		3.7,			//[3]
		9.0				//[4]
	};
    theta=1-T/crit.T;
    
    for (k=1;k<=4;k++)
    {
        sum+=Nbp[k]*pow(theta,tbp[k]);
    }
    return reduce.p.Pa*exp(reduce.T/T*sum);
}
    
double R404AClass::psatV(double T)
{
	//Dew point of R410A
	double sum=0,theta;
	int k;

	static const double Ndp[]={
		0.0,			//[0]
		-0.00026863,	//[1]
		-6.5757,		//[2]
		-4.1802,		//[3]
		-7.9102			//[4]
	};

	static const double tdp[]={
		0.0,			//[0]
		0.1,			//[1]
		0.972,			//[2]
		3.8,			//[3]
		9.0				//[4]
	};

    theta=1-T/crit.T;

    for (k=1;k<=4;k++)
    {
        sum+=Ndp[k]*pow(theta,tdp[k]);
    }
    return reduce.p.Pa*exp(reduce.T/T*sum);
}

double R404AClass::rhosatV(double T)
{
	double THETA,a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,b6,b7;
	THETA=1-T/345.27;

	a1 =      -17.61; 
    a2 =      -16.61; 
    a3 =      -6.059; 
    a4 =      -16.43; 
    a5 =      -16.85; 
    a6 =      -16.26; 
    a7 =      -2.926; 
    b1 =       2.939; 
    b2 =       6.538; 
    b3 =       1.144; 
    b4 =       6.536; 
    b5 =        6.54;  
    b6 =       6.533; 
    b7 =      0.3967;

	return exp(a1*pow(THETA,b1)+a2*pow(THETA,b2)+a3*pow(THETA,b3)+a4*pow(THETA,b4)+a5*pow(THETA,b5)+a6*pow(THETA,b6)+a7*pow(THETA,b7))*482.13;
}

double R404AClass::rhosatL(double T)
{
	double THETA,a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,b6,b7;
	THETA=1-T/345.27;

	a1 =      0.3897;
	a2 =      0.2768;
	a3 =       0.546; 
	a4 =      0.3473;
	a5 =     -0.1371;
	a6 =      0.2254;
	a7 =     -0.1426;
	b1 =      0.2581;
	b2 =       2.985;
	b3 =       0.255;
	b4 =      0.2596;
	b5 =       1.526;
	b6 =      0.2661;
	b7 =      0.0742;

	return exp(a1*pow(THETA,b1)+a2*pow(THETA,b2)+a3*pow(THETA,b3)+a4*pow(THETA,b4)+a5*pow(THETA,b5)+a6*pow(THETA,b6)+a7*pow(THETA,b7))*482.163;
}

double R404AClass::viscosity_Trho(double T, double rho)
{
    // Properties taken from "Viscosity of Mixed 
	// Refrigerants R404A,R407C,R410A, and R507A" 
	// by Vladimir Geller, 
	// 2000 Purdue Refrigeration conferences

    // inputs in T [K], and rho [kg/m^3]
    // output in Pa-s

   double 
	eta_microPa_s;

   //Set constants required
   double a_0=9.766e-1,a_1=3.676e-2,a_2=2.938e-6,b_1=2.260e-3,b_2=1.786e-4,
	   b_3=-4.202e-7,b_4=8.489e-10,b_5=-8.670e-13,b_6=3.566e-16;

   eta_microPa_s=a_0+a_1*T+a_2*T*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho+b_5*rho*rho*rho*rho*rho+b_6*rho*rho*rho*rho*rho*rho;
   return eta_microPa_s/1e6;
}
double R404AClass::surface_tension_T(double T)
{
	return 0.0538534*pow(1-T/reduce.T,1.259);
}
