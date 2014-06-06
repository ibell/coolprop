/* 
Properties of Propane (R290)
by Ian Bell
*/

#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include "math.h"
#include "stdio.h"
#include <string.h>
#include "CoolProp.h"
#include "FluidClass.h"
#include "R290.h"


//
////For the viscosity
//static const double tv[]={
//    0.0,			//[0]
//    0.0,			//[1]
//    0.0,			//[2]
//    0.0,			//[3]
//    0.0,			//[4]
//    1.0,			//[5]
//    1.0,			//[6]
//    2.0,			//[7]
//    2.0,			//[8]
//    2.0,			//[9]
//    3.0,			//[10]
//    4.0,			//[11]
//    1.0,			//[12]
//    2.0				//[13]
//};
//
//static const double dv[]={
//    0.0,			//[0]
//    1.0,			//[1]
//    2.0,			//[2]
//    3.0,			//[3]
//    13.0,			//[4]
//    12.0,			//[5]
//    16.0,			//[6]
//    0.0,			//[7]
//    18.0,			//[8]
//    20.0,			//[9]
//    13.0,			//[10]
//    4.0,			//[11]
//    0.0,			//[12]
//    1.0				//[13]
//};
//
//static const double nv[]={
//    0.0,			//[0]
//    -0.7548580e-1,	//[1]
//    0.7607150,		//[2]
//    -0.1665680,		//[3]
//    0.1627612e-5,	//[4]
//    0.1443764e-4,	//[5]
//    -0.2759086e-6,	//[6]
//    -0.1032756,		//[7]
//    -0.2498159e-7,	//[8]
//    0.4069891e-8,	//[9]
//    -0.1513434e-5,	//[10]
//    0.2591327e-2,	//[11]
//    0.5650076,		//[12]
//    0.1207253		//[13]
//};

R290Class::R290Class()
{

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

	static const double d[]={0,
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

	static const double c[]={
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

	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> b0_v(b0,b0+sizeof(b0)/sizeof(double));

	phirlist.push_back(new phir_power(n,d,t,c,1,11,12));
	phirlist.push_back(new phir_gaussian(n,d,t,alpha,epsilon,beta,GAMMA,12,18,19));

	/* phi0=log(delta)+3*log(tau)+a0[1]+a0[2]*tau
        +a0[3]*log(1-exp(-b0[3]*tau))
        +a0[4]*log(1-exp(-b0[4]*tau))
        +a0[5]*log(1-exp(-b0[5]*tau))
        +a0[6]*log(1-exp(-b0[6]*tau));
	*/

	phi0list.push_back(new phi0_lead(a0[1],a0[2]));
	phi0list.push_back(new phi0_logtau(3.0));
	phi0list.push_back(new phi0_Planck_Einstein(a0_v,b0_v,3,6));

	// Critical parameters
	crit.rho = 220.4781;
	crit.p = PressureUnit(4251.2, UNIT_KPA);
	crit.T = 369.89;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 44.09562;
	params.Ttriple = 85.525;
	params.ptriple = 1.71314090116e-07;
	params.accentricfactor = 0.1521;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = 85.525;
	limits.Tmax = 625.0;
	limits.pmax = 1000000.0;
	limits.rhomax = 20.6*params.molemass;
	
	EOSReference.assign("Lemmon, E.W., McLinden, M.O., Wagner, W., "
						"\"Thermodynamic Properties of Propane.  III.  A Reference Equation of State"
						"for Temperatures from the Melting Line to 650 K and Pressures up to 1000 MPa,\""
						"J. Chem. Eng. Data, 54:3141-3180, 2009");
	TransportReference.assign("Viscosity: E. Vogel, C. Kuchenmeister, and E. Bich, A. Laesecke,"
							"\"Reference Correlation of the Viscosity of Propane\""
							"J. Phys. Chem. Ref. Data, Vol. 27, No. 5, 1998\n\n"
							"Conductivity: Kenneth N. Marsh, Richard A. Perkins, and Maria L. V. Ramires "
							"\"Measurement and Correlation of the Thermal Conductivity of"
							"Propane from 86 K to 600 K at Pressures to 70 MPa\","
							"J. Chem. Eng. Data 2002, 47, 932-940 (Olchowy-Sengers)\n\n"
							"Surface Tension: A. Mulero and I. Cachadi�a and M. I. Parra"
							"\"Recommended Correlations for the Surface Tension of Common Fluids\""
							", J. Phys. Chem. Ref. Data, Vol. 41, No. 4, 2012");

	name.assign("n-Propane");
	aliases.push_back("Propane");
	aliases.push_back("propane");
	aliases.push_back("R290");
	aliases.push_back("C3H8");
	aliases.push_back(std::string("PROPANE"));
	aliases.push_back(std::string("N-PROPANE"));
	REFPROPname.assign("PROPANE");

	BibTeXKeys.EOS = "Lemmon-JCED-2009";
	BibTeXKeys.VISCOSITY = "Vogel-JPCRD-1998";
	BibTeXKeys.CONDUCTIVITY = "Marsh-JCED-2002";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Chichester-NIST-2008";
}
double R290Class::conductivity_background(double T, double rho)
{
	double sum=0, delta, tau;
	int i;
	double rhoc = reduce.rho;
	double Tc = 369.85;
	//Set constants required
    double B1[]={
        0.0,			//[0]
        -3.69500e-2,	//[1]
         1.48658e-1,	//[2]
        -1.19986e-1,	//[3]
         4.12431e-2,	//[4]
        -4.86905e-3		//[5]
    };
    double B2[]={
        0.0,			//[0]
         4.82798e-2,	//[1]
        -1.35636e-1,	//[2]
         1.17588e-1,	//[3]
        -4.36911e-2,	//[4]
         6.16079e-3		//[5]
    };
    delta=rho/rhoc;
    tau=Tc/T;
	for(i=1;i<=5;i++)
    {
        sum+=(B1[i]+B2[i]/tau)*pow(delta,(double)i);
    }
    return sum; //[W/m/K]
}
double R290Class::conductivity_Trho(double T, double rho)
{
	/*Properties taken from "Measurement and Correlation of the Thermal Conductivity of 
    Propane from 86 K to 600 K at Pressures to 70 MPa" 
    by Kenneth N. Marsh, Richard A. Perkins, and Maria L. V. Ramires
    J. Chem. Eng. Data 2002, 47, 932-940

    The empirical critical enhancement is implemented
    */
    
    // output in W/m-K

    double lambda0,lambdar,lambdac,tau = 369.85/T;
    
	// The dilute gas contribution [W/m/K]
	double A[]={0.0, -1.24778e-3, 8.16371e-3, 1.99374e-2};
    lambda0 = A[1]+A[2]/tau+A[3]/(tau*tau);
    
	// The background contribution [W/m/K]
	lambdar = this->conductivity_background(T,rho);

	// Critical term from Olchowy and Sengers using the value of qd^-1 from Marsh [W/m/K]
	lambdac = this->conductivity_critical(T,rho,1/7.16635e-10);

    return lambda0+lambdar+lambdac;
}

void R290Class::ECSParams(double *e_k, double *sigma)
{
	*e_k = 263.88;  *sigma = 0.49748;
}
double R290Class::viscosity_dilute(double T)
{	
	double eta_star, a[]={0.25104574,-0.47271238,0,0.060836515,0};
	double e_k = 263.88 /* K */, sigma = 0.49748, /* nm */ theta_star, Tstar;
	
	Tstar = T/e_k;
	theta_star = exp(a[0]*pow(log(Tstar),0)+a[1]*pow(log(Tstar),1)+a[3]*pow(log(Tstar),3));
	eta_star = 0.141824*sqrt(T)/(pow(sigma,2)*theta_star)/1e6;
	return eta_star;
}
double R290Class::viscosity_dilute2(double T, double rho)
{
	double e_k = 263.88 /* K */, sigma = 0.49748, /* nm */ Tstar, rhobar;
	double B_eta_star,B_eta,eta_1,sum=0; 
	Tstar = T/e_k;
	double b[]={-19.572881,219.73999,-1015.3226,2471.01251,-3375.1717,2491.6597,-787.26086,14.085455,-0.34664158};
	for (unsigned int i=0;i<=6;i++){
		sum += b[i]*pow(Tstar,-0.25*i);
	}
	B_eta_star = sum+b[7]*pow(Tstar,-2.5)+b[8]*pow(Tstar,-5.5); //[no units]
	B_eta = 0.602214129*pow(sigma,3)*B_eta_star; //[L/mol]
	rhobar = rho/params.molemass;
	eta_1 = viscosity_dilute(T) * B_eta;
	// B_eta*rho needs to be non-dimensional [m3/mol]*[kg/m3] so need divide by mole mass and multiply by 1000
	return eta_1*rhobar;
}
double R290Class::viscosity_higher_order(double T, double rho)
{
	double sum=0,delta_0,DELTA_H_eta,tau,delta,eta_r;
	double g1 = 2.50053938863, g2 =	0.860516059264,f1 = 1616.88405374;
	
	double e[6][3];
	e[2][0] =  35.9873030195; e[2][1] = -180.512188564; e[2][2] =  87.7124888223;
	e[3][0] = -105.773052525; e[3][1] =  205.319740877; e[3][2] = -129.210932610;
	e[4][0] =  58.9491587759; e[4][1] = -129.740033100; e[4][2] =  76.6280419971;
	e[5][0] = -9.59407868475; e[5][1] =  21.0726986598; e[5][2] = -14.3971968187;
	
	tau = T / 369.825;
	delta = rho / (44.09562*5); 

	delta_0=g1*(1+g2*sqrt(tau)); //[no units]
	sum=0;
	for (int i=2;i<=5;i++){
		for (int j=0;j<=2;j++){
			sum+=e[i][j]*pow(delta,i)/pow(tau,j);
		}
	}
	DELTA_H_eta = sum + f1*(delta/(delta_0-delta)-delta/delta_0);
	
	eta_r = DELTA_H_eta/1e6;
	return eta_r;
}
double R290Class::viscosity_residual(double T, double rho)
{
	return viscosity_background(T,rho);
}
double R290Class::viscosity_background(double T, double rho)
{
	return viscosity_dilute2(T,rho)+viscosity_higher_order(T,rho);
}
double R290Class::viscosity_Trho(double T, double rho)
{
	// inputs in T [K], and p [kPa]
    // output in Pa-s
	
	return viscosity_dilute(T)+viscosity_background(T,rho);
}
double R290Class::psat(double T)
{
	// Max error is  0.130389086256 % between 85.525 and 369.889999 K

	const double ti[]={0,0.09, 1.0, 2.1666666666666665, 2.5, 3.5, 6.666666666666667};
    const double Ni[]={0, -0.0036326603098915731, -6.4417163271195328, 4.2879114303836188, -4.3361540375495373, -1.5330938298966583, -0.8380028193117014};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double R290Class::rhosatL(double T)
{
	// Max error is  0.0403333824779 % between 85.525 and 369.889999 K
	const double ti[]={0, 0.3665, 0.3695, 0.39199999999999996, 0.39749999999999996, 1.5, 3.8333333333333335};
    const double Ni[]={0, -2544.5789267169393, 3262.5780670317858, -2279.9253725524513, 1564.3230394188856, 0.17057178175644663, 0.1801528182529071};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*(summer+1);
}
double R290Class::rhosatV(double T)
{
	// Max error is  0.119794943361 % between 85.525 and 369.889999 K
	const double ti[]={0,0.353, 0.374, 1.0, 2.0, 3.6666666666666665, 21.333333333333332};
    const double Ni[]={0,1.3470928285341803, -3.9397785053443175, -2.9376522081571999, 1.277500738496856, -4.2937206644727439, -2.2870921292077298};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
	return reduce.rho*exp(reduce.T/T*summer);
}
double R290Class::surface_tension_T(double T)
{
	return 0.05334*pow(1-T/reduce.T,1.235)-0.01748*pow(1-T/reduce.T,4.404);
}
