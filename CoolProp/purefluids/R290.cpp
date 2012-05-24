/* 
Properties of Propane (R290)
by Ian Bell
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

static const double Tc=369.89, rhoc=220.4781, Pc=4251.2, M_R290=44.09562, _Ttriple=85.525;
             //           K          kg/m^3       kPa            kg/kmol          K
static const double n[]={0,
0.042910051,
1.7313671,
-2.4516524,
0.34157466,
-0.46047898,
-0.66847295,
0.20889705,
0.19421381,
-0.22917851,
-0.60405866,
0.066680654,
0.017534618,
0.33874242,
0.22228777,
-0.23219062,
-0.09220694,
-0.47575718,
-0.017486824};

static const int d[]={0,
4, //[ 1]
1, //[ 2]
1, //[ 3]
2, //[ 4]
2, //[ 5]
1, //[ 6]
3, //[ 7]
6, //[ 8]
6, //[ 9]
2, //[10]
3, //[11]
1, //[12]
1, //[13]
1, //[14]
2, //[15]
2, //[16]
4, //[17]
1  //[18]
};

static const double t[]={0.00, //offset for natural indices
1.00,
0.33,
0.80,
0.43,
0.90,
2.46,
2.09,
0.88,
1.09,
3.25,
4.62,
0.76,
2.50,
2.75,
3.05,
2.55,
8.40,
6.75};

static const int c[]={
0,0,0,0,0,0, // indices [0-5]
1,
1,
1,
1,
2,
2,
0,0,0,0,0,0,0 // indices [12-18]
};

// alpha instead of eta is used here for consistency with the definitions in R744.c upon which R290.c is based
static const double alpha[]={
0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-11]
0.963,
1.977,
1.917,
2.307,
2.546,
3.28,
14.6};

static const double beta[]={
0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-11]
2.33,
3.47,
3.15,
3.19,
0.92,
18.8,
547.8};

static const double GAMMA[]={
0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-11]
0.684,
0.829,
1.419,
0.817,
1.5,
1.426,
1.093};

static const double epsilon[]={
0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-11]
1.283,
0.6936,
0.788,
0.473,
0.8577,
0.271,
0.948};

//Constants for ideal gas expression
static const double a0[]={0.0,
    -4.970583,
    4.29352,
    3.043,
    5.874,
    9.337,
    7.922
};

static const double b0[]={0.0,
    0,0, //[1 and 2 are not used]
    1.062478,
    3.344237,
    5.363757,
    11.762957
};

//For the viscosity
static const double tv[]={
    0.0,			//[0]
    0.0,			//[1]
    0.0,			//[2]
    0.0,			//[3]
    0.0,			//[4]
    1.0,			//[5]
    1.0,			//[6]
    2.0,			//[7]
    2.0,			//[8]
    2.0,			//[9]
    3.0,			//[10]
    4.0,			//[11]
    1.0,			//[12]
    2.0				//[13]
};

static const double dv[]={
    0.0,			//[0]
    1.0,			//[1]
    2.0,			//[2]
    3.0,			//[3]
    13.0,			//[4]
    12.0,			//[5]
    16.0,			//[6]
    0.0,			//[7]
    18.0,			//[8]
    20.0,			//[9]
    13.0,			//[10]
    4.0,			//[11]
    0.0,			//[12]
    1.0				//[13]
};

static const double nv[]={
    0.0,			//[0]
    -0.7548580e-1,	//[1]
    0.7607150,		//[2]
    -0.1665680,		//[3]
    0.1627612e-5,	//[4]
    0.1443764e-4,	//[5]
    -0.2759086e-6,	//[6]
    -0.1032756,		//[7]
    -0.2498159e-7,	//[8]
    0.4069891e-8,	//[9]
    -0.1513434e-5,	//[10]
    0.2591327e-2,	//[11]
    0.5650076,		//[12]
    0.1207253		//[13]
};

R290Class::R290Class()
{
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> d_v(d,d+sizeof(d)/sizeof(int));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(c,c+sizeof(c)/sizeof(int));
	std::vector<double> alpha_v(alpha,alpha+sizeof(alpha)/sizeof(double));
	std::vector<double> beta_v(beta,beta+sizeof(beta)/sizeof(double));
	std::vector<double> gamma_v(GAMMA,GAMMA+sizeof(GAMMA)/sizeof(double));
	std::vector<double> epsilon_v(epsilon,epsilon+sizeof(epsilon)/sizeof(double));
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> b0_v(b0,b0+sizeof(b0)/sizeof(double));

	phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,11);
	phi_BC * phirg_ = new phir_gaussian(n_v,d_v,t_v,alpha_v,epsilon_v,beta_v,gamma_v,12,18);
	phirlist.push_back(phir_);
	phirlist.push_back(phirg_);

	/* phi0=log(delta)+3*log(tau)+a0[1]+a0[2]*tau
        +a0[3]*log(1-exp(-b0[3]*tau))
        +a0[4]*log(1-exp(-b0[4]*tau))
        +a0[5]*log(1-exp(-b0[5]*tau))
        +a0[6]*log(1-exp(-b0[6]*tau));
	*/
	phi_BC * phi0_lead_ = new phi0_lead(a0[1],a0[2]);
	phi_BC * phi0_logtau_ = new phi0_logtau(3.0);
	phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(a0_v,b0_v,3,6);

	phi0list.push_back(phi0_lead_);
	phi0list.push_back(phi0_logtau_);
	phi0list.push_back(phi0_Planck_Einstein_);

	// Critical parameters
	crit.rho = 220.4781;
	crit.p = 4251.2;
	crit.T = 369.89;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 44.09562;
	params.Ttriple = 85.525;
	params.accentricfactor = 0.1521;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = 85.525;
	limits.Tmax = 625.0;
	limits.pmax = 1000000.0;
	limits.rhomax = 20.6*params.molemass;
	
	EOSReference.assign("\"A International Standard Formulation for the Thermodynamic Properties of 1,1,1,2-Tetrafluoroethane" 
					    "(HFC-134a) for Temperatures from 170 K to 455 K and Pressures up to 70 MPa\""
						"by Reiner Tillner-Roth and Hans Dieter Baehr, J. Phys. Chem. Ref. Data, v. 23, 1994, pp 657-729");
	TransportReference.assign("Viscosity: \"A Reference Multiparameter Viscosity Equation for R290"
						"with an Optimized Functional Form\""
						"by G. Scalabrin and P. Marchi, R. Span"
						"J. Phys. Chem. Ref. Data, Vol. 35, No. 2, 2006");

	name.assign("R290");
	aliases.push_back("propane");
	aliases.push_back("Propane");
	aliases.push_back("C3H8");
}
double R290Class::conductivity_Trho(double T, double rho)
{
	/*Properties taken from "Measurement and Correlation of the Thermal Conductivity of 
    Propane from 86 K to 600 K at Pressures to 70 MPa" 
    by Kenneth N. Marsh, Richard A. Perkins, and Maria L. V. Ramires
    J. Chem. Eng. Data 2002, 47, 932-940

    The empirical critical enhancement is implemented
    */
    
    // output in kW/m-K

    double delta,lambda0,lambdar,lambdac,sum=0,tau;
    double DELTAT_c,DELTArho_c;
    int i;

    //Set constants required
    double B1[]={
        0.0,			//[0]
        -3.51153e-2,	//[1]
         1.70890e-1,	//[2]
        -1.47688e-1,	//[3]
         5.19283e-2,	//[4]
        -6.18662e-3		//[5]
    };
    double B2[]={
        0.0,			//[0]
         4.69195e-2,	//[1]
        -1.48616e-1,	//[2]
         1.32457e-1,	//[3]
        -4.85636e-2,	//[4]
         6.60414e-3		//[5]
    };
    double C[]={
        0.0,			//[0]
         3.66486e-4,	//[1]
        -2.21696e-3,	//[2]
         2.64213e+0		//[3]
    };
    double A[]={
         0.0,			//[0]
        -1.24778e-3,	//[1]
         8.16371e-3,	//[2]
         1.99374e-2,	//[3]
    };
    delta=rho/rhoc;
    tau=Tc/T;
    lambda0=A[1]+A[2]/tau+A[3]/(tau*tau);
    for(i=1;i<=5;i++)
    {
        sum+=(B1[i]+B2[i]/tau)*pow(delta,(double)i);
    }
    lambdar=sum;
    DELTAT_c=(1.0/tau-1.0);
    DELTArho_c=delta-1.0;
    lambdac=C[1]/(C[2]+fabs(DELTAT_c))*exp(-(C[3]*DELTArho_c)*(C[3]*DELTArho_c));

    return (lambda0+lambdar+lambdac)/1000.0;
}
double R290Class::viscosity_Trho(double T, double rho)
{
	/* 
    Properties taken from "A Reference Multiparameter Viscosity Equation 
    for Propane with an Optimized Functional Form" 
    by G. Scalabrin and P. Marchi and R. Span
    J. Phys. Chem. Ref. Data, Vol. 35, No. 3, 2006, 1415-1442
    */

    // inputs in T [K], and p [kPa]
    // output in Pa-s

    int i;
    double Tr,rhor,Hc=17.1045,etar, sum=0;

    //Reduced Temperature
    Tr=T/Tc;
    rhor=rho/rhoc;
    for(i=1;i<=11;i++)
    {
        sum+=nv[i]*pow(Tr,tv[i])*pow(rhor,dv[i]);
    }
    for(i=12;i<=13;i++)
    {
        sum+=exp(-rhor*rhor/2.0)*nv[i]*pow(Tr,tv[i])*pow(rhor,dv[i]);
    }
    etar=sum;

    return (exp(etar)-1)*Hc/1e6;
}
double R290Class::psat(double T)
{
	const double ti[]={0,1.0,1.5,2.2,4.8,6.2};
    const double Ni[]={0,-6.7722,1.6938,-1.3341,-3.1876,0.94937};
    double summer=0,theta;
    int i;
    theta=1-T/Tc;
    for (i=1;i<=5;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return Pc*exp(Tc/T*summer);
}
double R290Class::rhosatL(double T)
{
	const double ti[]={0,0.345,0.74,2.6,7.2};
    const double Ni[]={0,1.82205,0.65802,0.21109,0.083973};
    double summer=1;
    int i;
    double theta;
    theta=1-T/Tc;
    for (i=1;i<=4;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return rhoc*summer;
}
double R290Class::rhosatV(double T)
{
	const double ti[]={0,0.3785,1.07,2.7,5.5,10,20};
    const double Ni[]={0,-2.4887,-5.1069,-12.174,-30.495,-52.192,-134.89};
    double summer=0,theta;
    int i;
    theta=1.0-T/Tc;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return rhoc*exp(summer);
}