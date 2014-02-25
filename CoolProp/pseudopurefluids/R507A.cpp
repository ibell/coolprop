/*
Properties for R507A. 
by Ian Bell (ihb2@cornell.edu)

Pseudo-pure fluid thermo props from 
"Pseudo-pure fluid Equations of State for the Refrigerant Blends R410A, R404A, R507C and R407C" 
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
#include "R507A.h"

static const double a[]={
    9.93541,		//[0]
    7.9985,		    //[1]
    -21.6054,		//[2]
    0.95006,		//[3]
    4.1887,			//[4]
    5.5184			//[5]
};

static const double b[]={
    0,				//[0]
    0,				//[1]
    -0.25,			//[2]
    1.05886,		//[3]
    2.37081,		//[4]
    5.14305			//[5]
};

static const double N[]={
	 0.0,			//[0]
	 6.24982,		//[1]
	-8.07855,		//[2]
	 0.0264843,		//[3]
	 0.286215,		//[4]
	-0.00507076,	//[5]
	 0.0109552,		//[6]
	 0.00116124,	//[7]
	 1.38469,		//[8]
	-0.922473,		//[9]
	-0.0503562,		//[10]
	 0.822098,		//[11]
	-0.277727,		//[12]
	 0.358172,		//[13]
	-0.0126426,		//[14]
	-0.00607010,		//[15]
	-0.0815653,		//[16]
	-0.0233323,		//[17]
	 0.0352952,		//[18]
	 0.0159566,		//[19]
	 0.0755927,		//[20]
	-0.0542007,		//[21]
	 0.0170451		//[22]

};

static const double t[]={
	0.0,			//[0]
    0.692,			//[1]
    0.943,			//[2]
    5.8,			//[3]
    0.77,			//[4]
    5.84,			//[5]
    0.24,			//[6]
    0.69,			//[7]
    2.0,			//[8]
    3.0,			//[9]
    7.0,			//[10]
    2.2,			//[11]
    4.3,			//[12]
    2.7,			//[13]
    1.2,			//[14]
    1.23,			//[15]
    12.0,			//[16]
    6.0,			//[17]
    8.5,			//[18]
    11.5,			//[19]
    13.0,			//[20]
    17.0,			//[21]
	16.2			//[22]
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

static const double Nbp[]={
	0.0,			//[0]
	-7.4853,		//[1]
    2.0115,			//[2]
    -2.0141,		//[3]
    -3.7763			//[4]
};
    
static const double tbp[]={
	0.0,			//[0]
    1.0,			//[1]
    1.5,			//[2]
    2.2,			//[3]
    4.6				//[4]
};

static const double Ndp[]={
    0.0,			//[0]
	-7.5459,		//[1]
    2.3380,			//[2]
    -2.2370,		//[3]
    -4.1535			//[4]
};

static const double tdp[]={
	0.0,			//[0]
    1.0,			//[1]
    1.5,			//[2]
    2.1,			//[3]
    4.7				//[4]
};



R507AClass::R507AClass()
{
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
	crit.rho = 490.74;
	crit.p = PressureUnit(3704.9,UNIT_KPA);
	crit.T = 343.765;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 98.8592;
	params.Ttriple = 200.0;
	params.ptriple = 23.2234432758;
	params.accentricfactor = 0.286;
	params.R_u = 8.314472;
	isPure = false;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 50000.0;
	limits.rhomax = 14.13*params.molemass;
	
	EOSReference.assign("E.W. Lemmon, \"Pseudo-pure fluid Equations of State for the Refrigerant Blends R410A, R404A, R507C and R407C\"" 
						",Int. J. Thermophys. v. 24, n4, 2003");
	TransportReference.assign("Viscosity: V. Geller, \"Viscosity of Mixed Refrigerants R404A,R407C,"
							"R410A, and R507A\", 2000 Purdue Refrigeration conferences\n\n"
							"Thermal Conductivity: V.Z. Geller, B.Z. Nemzer, and U.V. Cheremnykh \"Thermal Conductivity "
							"of the Refrigerant mixtures R404A,R407C,R410A, and R507A\" "
							"Int. J. Thermophysics, v. 22, n 4 2001");

	name.assign("R507A");
	aliases.push_back("R507a");

	BibTeXKeys.EOS = "Lemmon-IJT-2003";
	BibTeXKeys.VISCOSITY = "Geller-PURDUE-2000";
	BibTeXKeys.CONDUCTIVITY = "Geller-IJT-2001";
	
}
double R507AClass::psatL(double T)
{
	//Bubble point of R410A
	double sum=0,theta;
	int k;
    theta=1-T/crit.T;
    
    for (k=1;k<=4;k++)
    {
        sum+=Nbp[k]*pow(theta,tbp[k]);
    }
    return reduce.p.Pa*exp(reduce.T/T*sum);
}
    
double R507AClass::psatV(double T)
{
	//Dew point of R410A
	double sum=0,theta;
	int k;
    theta=1-T/crit.T;

    for (k=1;k<=4;k++)
    {
        sum+=Ndp[k]*pow(theta,tdp[k]);
    }
    return reduce.p.Pa*exp(reduce.T/T*sum);
}

double R507AClass::rhosatV(double T)
{
	double theta;
	theta=1-T/343.765;
	return exp(+2.685721675+2.703462155742e+00*pow(theta,-0.031591)-12.9219949523*theta+48.3326995201*theta*theta-301.352305753*powInt(theta,3)+953.161091912*powInt(theta,4)-1667.82041758*powInt(theta,5)+1132.7963359*powInt(theta,6));
}

double R507AClass::rhosatL(double T)
{

	double theta;
	theta=1-T/343.765;
	return exp(6.059998874+1.305406634077e+00*pow(theta,0.21588)+0.585710808997*theta-2.09578513606*theta*theta+4.667514451*powInt(theta,3)-3.93369379654*powInt(theta,4));
}

double R507AClass::viscosity_Trho(double T, double rho)
{
    // Properties taken from "Viscosity of Mixed 
	// Refrigerants R404A,R407C,R410A, and R507A" 
	// by Vladimir Geller, 
	// 2000 Purdue Refrigeration conferences

    // inputs in T [K], and p [kPa]
    // output in Pa-s

   double eta_microPa_s;

   //Set constants required
   double a_0=-2.530e0,a_1=5.626e-2,a_2=-2.323e-5,b_1=5.308e-4,b_2=2.234e-4,
	   b_3=-6.742e-7,b_4=1.411e-9,b_5=-1.388e-12,b_6=5.274e-16;

   eta_microPa_s=a_0+a_1*T+a_2*T*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho+b_5*rho*rho*rho*rho*rho+b_6*rho*rho*rho*rho*rho*rho;
   return eta_microPa_s/1e6;
}

double R507AClass::conductivity_Trho(double T, double rho)
{
	// Properties taken from "Thermal Conductivity 
	// of the Refrigerant mixtures R404A,R407C,R410A, and R507A" 
	// by V.Z. Geller, B.Z. Nemzer, and U.V. Cheremnykh 
	// Int. J. Thermophysics, v. 22, n 4 2001

	// inputs in T [K], and p [kPa] or rho [kg/m^3]
	// output in W/m-K

	//Set constants required
	double a_0=-8.656e0,a_1=7.383e-2,b_1=2.799e-2,b_2=3.065e-5,b_3=-3.644e-8,b_4=2.609e-11;

	return (a_0+a_1*T+b_1*rho+b_2*rho*rho+b_3*rho*rho*rho+b_4*rho*rho*rho*rho)/1.e3; // from mW/m-K to W/m-K
}
