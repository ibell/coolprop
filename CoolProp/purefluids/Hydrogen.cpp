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
#include "Hydrogen.h"

static const double n[]={0,
-6.93643,
0.01,
2.1101,
4.52059,
0.732564,
-1.34086,
0.130985,
-0.777414,
0.351944,
-0.0211716,
0.0226312,
0.032187,
-0.0231752,
0.0557346
};

static const int d[]={0,
1, //[ 1]
4, //[ 2]
1, //[ 3]
1, //[ 4]
2, //[ 5]
2, //[ 6]
3, //[ 7]
1, //[ 8]
3, //[ 9]
2, //[10]
1, //[11]
3, //[12]
1, //[13]
1, //[14]
};

static const double t[]={0.00, //offset for natural indices
0.6844,
1.0,
0.989,
0.489,
0.803,
1.1444,
1.409,
1.754,
1.311,
4.187,
5.646,
0.791,
7.249,
2.986
};

static const int c[]={
0,0,0,0,0,0,0,0, // indices [0-7]
1,
1
};

// alpha instead of eta is used here for consistency with the definitions in R744.c upon which R290.c is based
static const double phi[]={
0,0,0,0,0,0,0,0,0,0, // indices [0-9]
1.685,
0.489,
0.103,
2.506,
1.607
};

static const double beta[]={
0,0,0,0,0,0,0,0,0,0, // indices [0-9]
0.171,
0.2245,
0.1304,
0.2785,
0.3967
};

static const double GAMMA[]={
0,0,0,0,0,0,0,0,0,0, // indices [0-9]
0.7164,
1.3444,
1.4517,
0.7204,
1.5445
};

static const double D[]={
0,0,0,0,0,0,0,0,0,0, // indices [0-9]
1.506,
0.156,
1.736,
0.67,
1.662
};

//Constants for ideal gas expression
static const double a0[]={0.0,
	-1.4579856475,
	1.888076782,
	1.616,
	-0.4117,
	-0.792,
	0.758,
	1.217
};

static const double b0[]={0.0,
    0,0, //[1 and 2 are not used]
    16.0205159149,
	22.6580178006,
	60.0090511389,
	74.9434303817,
	206.9392065168
};

HydrogenClass::HydrogenClass()
{
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> d_v(d,d+sizeof(d)/sizeof(int));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(c,c+sizeof(c)/sizeof(int));
	std::vector<double> alpha_v(phi,phi+sizeof(phi)/sizeof(double));
	std::vector<double> beta_v(beta,beta+sizeof(beta)/sizeof(double));
	std::vector<double> gamma_v(GAMMA,GAMMA+sizeof(GAMMA)/sizeof(double));
	std::vector<double> epsilon_v(D,D+sizeof(D)/sizeof(double));
	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> b0_v(b0,b0+sizeof(b0)/sizeof(double));

	phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,9);
	phi_BC * phirg_ = new phir_gaussian(n_v,d_v,t_v,alpha_v,epsilon_v,beta_v,gamma_v,10,14);
	phirlist.push_back(phir_);
	phirlist.push_back(phirg_);

	/* phi0=log(delta)+1.5*log(tau)+a0[1]+a0[2]*tau
        +a0[3]*log(1-exp(-b0[3]*tau))
        +a0[4]*log(1-exp(-b0[4]*tau))
        +a0[5]*log(1-exp(-b0[5]*tau))
        +a0[6]*log(1-exp(-b0[6]*tau))
        +a0[7]*log(1-exp(-b0[7]*tau));
	*/

	//lead term of the form log(delta)+a1+a2*tau
	phi_BC * phi0_lead_ = new phi0_lead(a0[1],a0[2]);
	phi_BC * phi0_logtau_ = new phi0_logtau(1.5);
	phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(a0_v,b0_v,3,7);

	phi0list.push_back(phi0_lead_);
	phi0list.push_back(phi0_logtau_);
	phi0list.push_back(phi0_Planck_Einstein_);

	// Critical parameters
	crit.rho = 15.508*2.01588;
	crit.p = 1296.4;
	crit.T = 33.145;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 2.01588;
	params.Ttriple = 13.957;
	params.ptriple = 7.36205784168;
	params.accentricfactor = -0.219;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = 13.957;
	limits.Tmax = 1000.0;
	limits.pmax = 2000000.0;
	limits.rhomax = 102.0*params.molemass;
	
	EOSReference.assign("\"Fundamental Equations of State for Parahydrogen, Normal Hydrogen, and Orthohydrogen\""
						"by J.W. Leachman and R.T. Jacobsen and S.G. Penoncello and E.W. Lemmon"
						", J. Phys. Chem. Ref. Data, Vol. 38, No. 3, 2009, pp 721-748");
	TransportReference.assign("Viscosity & Surface Tension: McCarty, R.D. and Weber, L.A., "
							  "\"Thermophysical properties of parahydrogen from the freezing liquid line to "
							  "5000 R for pressures to 10,000 psia,\" "
                              "Natl. Bur. Stand., Tech. Note 617, 1972.");

	name.assign("Hydrogen");
	aliases.push_back("hydrogen");
	aliases.push_back("H2");
	aliases.push_back("R702");
	REFPROPname.assign("hydrogen");
}

double HydrogenClass::psat(double T)
{
	const double ti[]={0,1.0,1.5,2.0,2.85};
    const double Ni[]={0,-4.89789,0.988558,0.349689,0.499356};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double HydrogenClass::rhosatL(double T)
{
	double theta = 1-T/reduce.T;
	double RHS,rho;
	// Max error is 1.585900 %
	RHS = +1.093049*pow(theta,0.3)
		  +0.260558*pow(theta,0.9)
		  -0.543586*pow(theta,1.5)
		  +0.333107*pow(theta,3.2);
	rho = exp(RHS)*reduce.rho;	
	return rho;

}
double HydrogenClass::rhosatV(double T)
{
	double theta = 1-T/reduce.T;
	double RHS,rho;
	// Max error is 2.044809 %
	RHS = -1.120772*pow(theta,0.3)
          -4.767213*pow(theta,0.9)
          +1.108207*pow(theta,1.5)
          -11.752668*pow(theta,3.2);
	rho = exp(RHS)*reduce.rho;
	return rho;
}

double HydrogenClass::conductivity_Trho(double T, double rho)
{
	double e_k = 59.7,
		   sigma = 0.2827;
	double eta_dilute = viscosity_dilute(T,e_k,sigma);
	return 15e-3*R()*(eta_dilute*1e6)/4.0/1000;
}
double HydrogenClass::viscosity_Trho(double T, double rho)
{
	double A,B,eta_0,eta_E,e_k,sigma;

    if (T < 100){
        // Dilute gas contribution
        eta_0 = 8.5558 * (pow(T,1.5)/(T+19.55)) * ((T+650.39)/(T+1175.9));

		//For the excess part, density is in units of g/cm3, so need to divide by 1000
		rho /= 1000;

        A = (306.4636*rho - 3350.628*rho*rho+3868.092*rho*rho*rho)/(1.0-18.47830*rho+110.915*rho*rho+25.3542*rho*rho*rho);

        B = 10.0 + 7.2*(pow(rho/0.07,6)-pow(rho/0.07,1.5))-17.63/exp(58.75*pow(rho/0.07,3));

        // Excess viscosity
        eta_E = A * exp(B/T);

		// Correlation in units of g/cm-s x 10e-6, so to get Pa-s, need to divide by 10 and divide by 1e6
		
		return (eta_0 + eta_E) / 1e7;
    }
    else{
		// Use dilute gas properties
		e_k = 59.7;
		sigma = 0.2827;
		return viscosity_dilute(T,e_k,sigma);
    }
}

double HydrogenClass::surface_tension_T(double T)
{
	return 0.005369*pow(1-T/reduce.T,1.065);
}



ParaHydrogenClass::ParaHydrogenClass()
{
	double _n [] = {0, -7.33375, 0.01, 2.60375, 4.66279, 0.68239, -1.47078, 0.135801, -1.05327, 0.328239, -0.0577833, 0.0449743, 0.0703464, -0.0401766, 0.11951};
	double _t [] = {0, 0.6855, 1, 1, 0.489, 0.774, 1.133, 1.386, 1.619, 1.162, 3.96, 5.276, 0.99, 6.791, 3.19};
	double _d [] = {0, 1, 4, 1, 1, 2, 2, 3, 1, 3, 2, 1, 3, 1, 1};
	double _c [] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
	double _a0 [] = {0, -1.4485891134, 1.884521239, 4.30256, 13.0289, -47.7365, 50.0013, -18.6261, 0.993973, 0.536078};
	double _b0 [] = {0, 0, 0, -15.14967511472, -25.0925982148, -29.4735563787, -35.4059141417, -40.724998482, -163.7925799988, -309.2173173842};
	double _phi[] = {0,0,0,0,0,0,0,0,0,0,1.7437, 0.5516, 0.0634, 2.1341, 1.777};
	double _beta[] = {0,0,0,0,0,0,0,0,0,0,0.194, 0.2019, 0.0301, 0.2383, 0.3253};
	double _GAMMA[] = {0,0,0,0,0,0,0,0,0,0,0.8048, 1.5248, 0.6648, 0.6832, 1.493};
	double _D[] = {0,0,0,0,0,0,0,0,0,0,1.5487, 0.1785, 1.28, 0.6319, 1.7104};
	std::vector<double> n_v(_n,_n+sizeof(_n)/sizeof(double));
	std::vector<double> d_v(_d,_d+sizeof(_d)/sizeof(double));
	std::vector<double> t_v(_t,_t+sizeof(_t)/sizeof(double));
	std::vector<double> l_v(_c,_c+sizeof(_c)/sizeof(double));
	std::vector<double> alpha_v(_phi,_phi+sizeof(_phi)/sizeof(double));
	std::vector<double> epsilon_v(_D,_D+sizeof(_D)/sizeof(double));
	std::vector<double> beta_v(_beta,_beta+sizeof(_beta)/sizeof(double));
	std::vector<double> gamma_v(_GAMMA,_GAMMA+sizeof(_GAMMA)/sizeof(double));
	std::vector<double> a0_v(_a0,_a0+sizeof(_a0)/sizeof(double));
	std::vector<double> b0_v(_b0,_b0+sizeof(_b0)/sizeof(double));

	phirlist.push_back(new phir_power(n_v,d_v,t_v,l_v,1,9));
	phirlist.push_back(new phir_gaussian(n_v,d_v,t_v,alpha_v,epsilon_v,beta_v,gamma_v,10,14));

	/* phi0=log(delta)+1.5*log(tau)+a0[1]+a0[2]*tau
        +a0[3]*log(1-exp(-b0[3]*tau))
        +a0[4]*log(1-exp(-b0[4]*tau))
        +a0[5]*log(1-exp(-b0[5]*tau))
        +a0[6]*log(1-exp(-b0[6]*tau))
        +a0[7]*log(1-exp(-b0[7]*tau));
	*/

	//lead term of the form log(delta)+a1+a2*tau
	phi_BC * phi0_lead_ = new phi0_lead(a0[1],a0[2]);
	phi_BC * phi0_logtau_ = new phi0_logtau(1.5);
	phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(a0_v,b0_v,3,7);

	phi0list.push_back(phi0_lead_);
	phi0list.push_back(phi0_logtau_);
	phi0list.push_back(phi0_Planck_Einstein_);

	// Critical parameters
	crit.rho = 15.538*2.01588;
	crit.p = 1285.8;
	crit.T = 32.938;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 2.01588;
	params.Ttriple = 13.8033;
	params.ptriple = 7.041;
	params.accentricfactor = -0.219;
	params.R_u = 8.314472;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 1000.0;
	limits.pmax = 2000000.0;
	limits.rhomax = 102.0*params.molemass;
	
	EOSReference.assign("\"Fundamental Equations of State for Parahydrogen, Normal Hydrogen, and Orthohydrogen\""
						"by J.W. Leachman and R.T. Jacobsen and S.G. Penoncello and E.W. Lemmon"
						", J. Phys. Chem. Ref. Data, Vol. 38, No. 3, 2009, pp 721-748");
	TransportReference.assign("Viscosity & Surface Tension: McCarty, R.D. and Weber, L.A., "
							  "\"Thermophysical properties of parahydrogen from the freezing liquid line to "
							  "5000 R for pressures to 10,000 psia,\" "
                              "Natl. Bur. Stand., Tech. Note 617, 1972.");

	name.assign("ParaHydrogen");
	REFPROPname.assign("PARAHYD");
}

double ParaHydrogenClass::psat(double T)
{
	const double ti[]={0,1.0,1.5,2.65,7.4};
    const double Ni[]={0,-4.87767,1.03359,0.82668,-0.129412};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    double rho = reduce.p*exp(reduce.T/T*summer);
	return rho;
}
double ParaHydrogenClass::rhosatL(double T)
{
	double theta = 1-T/reduce.T;
	double RHS,rho;
	// Max error is 1.585900 %
	RHS = +1.093049*pow(theta,0.3)
		  +0.260558*pow(theta,0.9)
		  -0.543586*pow(theta,1.5)
		  +0.333107*pow(theta,3.2);
	rho = exp(RHS)*reduce.rho;	
	return rho;

}
double ParaHydrogenClass::rhosatV(double T)
{
	double theta = 1-T/reduce.T;
	double RHS,rho;
	// Max error is 2.044809 %
	RHS = -1.120772*pow(theta,0.3)
          -4.767213*pow(theta,0.9)
          +1.108207*pow(theta,1.5)
          -11.752668*pow(theta,3.2);
	rho = exp(RHS)*reduce.rho;
	return rho;
}

double ParaHydrogenClass::conductivity_Trho(double T, double rho)
{
	double e_k = 59.7,
		   sigma = 0.2827;
	double eta_dilute = viscosity_dilute(T,e_k,sigma);
	return 15e-3*R()*(eta_dilute*1e6)/4.0/1000;
}
double ParaHydrogenClass::viscosity_Trho(double T, double rho)
{
	double A,B,eta_0,eta_E,e_k,sigma;

    if (T < 100){
        // Dilute gas contribution
        eta_0 = 8.5558 * (pow(T,1.5)/(T+19.55)) * ((T+650.39)/(T+1175.9));

		//For the excess part, density is in units of g/cm3, so need to divide by 1000
		rho /= 1000;

        A = (306.4636*rho - 3350.628*rho*rho+3868.092*rho*rho*rho)/(1.0-18.47830*rho+110.915*rho*rho+25.3542*rho*rho*rho);

        B = 10.0 + 7.2*(pow(rho/0.07,6)-pow(rho/0.07,1.5))-17.63/exp(58.75*pow(rho/0.07,3));

        // Excess viscosity
        eta_E = A * exp(B/T);

		// Correlation in units of g/cm-s x 10e-6, so to get Pa-s, need to divide by 10 and divide by 1e6
		
		return (eta_0 + eta_E) / 1e7;
    }
    else{
		// Use dilute gas properties
		e_k = 59.7;
		sigma = 0.2827;
		return viscosity_dilute(T,e_k,sigma);
    }
}

double ParaHydrogenClass::surface_tension_T(double T)
{
	return 0.005369*pow(1-T/reduce.T,1.065);
}

//OrthoHydrogenClass::OrthoHydrogenClass()
//{
//	double _n [] = {0, -6.83148, 0.01, 2.11505, 4.38353, 0.211292, -1.00939, 0.142086, -0.87696, 0.804927, -0.710775, 0.0639688, 0.0710858, -0.087654, 0.647088};
//	double _t [] = {0, 0.7333, 1, 1.1372, 0.5136, 0.5638, 1.6248, 1.829, 2.404, 2.105, 4.1, 7.658, 1.259, 7.589, 3.946};
//	double _d [] = {0, 1, 4, 1, 1, 2, 2, 3, 1, 3, 2, 1, 3, 1, 1};
//	double _c [] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
//	double _a0 [] = {0, -1.4675442336, 1.8845068862, 2.54151, -2.3661, 1.00365, 1.22447};
//	double _b0 [] = {0, 0, 0, -25.7676098736, -43.4677904877, -66.0445514750, -209.7531607465};
//    double _phi[] = {0,0,0,0,0,0,0,0,0,0,1.169, 0.894, 0.04, 2.072, 1.306};
//    double _beta[] = {0,0,0,0,0,0,0,0,0,0,0.4555, 0.4046, 0.0869, 0.4415, 0.5743};
//    double _GAMMA[] = {0,0,0,0,0,0,0,0,0,0,1.5444, 0.6627, 0.763, 0.6587, 1.4327};
//    double _D[] = {0,0,0,0,0,0,0,0,0,0,0.6366, 0.3876, 0.9437, 0.3976, 0.9626};
//	std::vector<double> n_v(_n,_n+sizeof(_n)/sizeof(double));
//	std::vector<double> d_v(_d,_d+sizeof(_d)/sizeof(double));
//	std::vector<double> t_v(_t,_t+sizeof(_t)/sizeof(double));
//	std::vector<double> l_v(_c,_c+sizeof(_c)/sizeof(double));
//	std::vector<double> alpha_v(_phi,_phi+sizeof(_phi)/sizeof(double));
//	std::vector<double> beta_v(_beta,_beta+sizeof(_beta)/sizeof(double));
//	std::vector<double> gamma_v(_GAMMA,_GAMMA+sizeof(_GAMMA)/sizeof(double));
//	std::vector<double> epsilon_v(_D,_D+sizeof(_D)/sizeof(double));
//	std::vector<double> a0_v(_a0,_a0+sizeof(_a0)/sizeof(double));
//	std::vector<double> b0_v(_b0,_b0+sizeof(_b0)/sizeof(double));
//
//	phirlist.push_back(new phir_power(n_v,d_v,t_v,l_v,1,9));
//	phirlist.push_back(new phir_gaussian(n_v,d_v,t_v,alpha_v,epsilon_v,beta_v,gamma_v,10,14));
//
//	/* phi0=log(delta)+1.5*log(tau)+a0[1]+a0[2]*tau
//        +a0[3]*log(1-exp(-b0[3]*tau))
//        +a0[4]*log(1-exp(-b0[4]*tau))
//        +a0[5]*log(1-exp(-b0[5]*tau))
//        +a0[6]*log(1-exp(-b0[6]*tau))
//        +a0[7]*log(1-exp(-b0[7]*tau));
//	*/
//
//	//lead term of the form log(delta)+a1+a2*tau
//	phi_BC * phi0_lead_ = new phi0_lead(a0[1],a0[2]);
//	phi_BC * phi0_logtau_ = new phi0_logtau(1.5);
//	phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(a0_v,b0_v,3,7);
//
//	phi0list.push_back(phi0_lead_);
//	phi0list.push_back(phi0_logtau_);
//	phi0list.push_back(phi0_Planck_Einstein_);
//
//	// Critical parameters
//	crit.rho = 15.445*2.01588;
//	crit.p = 1310.65;
//	crit.T = 33.22;
//	crit.v = 1.0/crit.rho;
//
//	// Other fluid parameters
//	params.molemass = 2.01594;
//	params.Ttriple = 14.008;
//	params.ptriple = 7.461;
//	params.accentricfactor = -0.219;
//	params.R_u = 8.314472;
//
//	// Limits of EOS
//	limits.Tmin = params.Ttriple;
//	limits.Tmax = 1000.0;
//	limits.pmax = 2000000.0;
//	limits.rhomax = 102.0*params.molemass;
//	
//	EOSReference.assign("\"Fundamental Equations of State for Parahydrogen, Normal Hydrogen, and Orthohydrogen\""
//						"by J.W. Leachman and R.T. Jacobsen and S.G. Penoncello and E.W. Lemmon"
//						", J. Phys. Chem. Ref. Data, Vol. 38, No. 3, 2009, pp 721-748");
//	TransportReference.assign("Viscosity & Surface Tension: McCarty, R.D. and Weber, L.A., "
//							  "\"Thermophysical properties of parahydrogen from the freezing liquid line to "
//							  "5000 R for pressures to 10,000 psia,\" "
//                              "Natl. Bur. Stand., Tech. Note 617, 1972.");
//
//	name.assign("OrthoHydrogen");
//	REFPROPname.assign("ORTHOHYD");
//}
//
//double OrthoHydrogenClass::psat(double T)
//{
//	const double ti[]={0,1.0,1.5,2.7,6.2};
//    const double Ni[]={0,-4.88684,1.05310,0.856947,-0.185355};
//    double summer=0,theta;
//    int i;
//    theta=1-T/reduce.T;
//    for (i=1;i<=4;i++)
//    {
//        summer=summer+Ni[i]*pow(theta,ti[i]);
//    }
//    return reduce.p*exp(reduce.T/T*summer);
//}
//double OrthoHydrogenClass::rhosatL(double T)
//{
//	double theta = 1-T/reduce.T;
//	double RHS,rho;
//	// Max error is 1.585900 %
//	RHS = +1.093049*pow(theta,0.3)
//		  +0.260558*pow(theta,0.9)
//		  -0.543586*pow(theta,1.5)
//		  +0.333107*pow(theta,3.2);
//	rho = exp(RHS)*reduce.rho;	
//	return rho;
//
//}
//double OrthoHydrogenClass::rhosatV(double T)
//{
//	double theta = 1-T/reduce.T;
//	double RHS,rho;
//	// Max error is 2.044809 %
//	RHS = -1.120772*pow(theta,0.3)
//          -4.767213*pow(theta,0.9)
//          +1.108207*pow(theta,1.5)
//          -11.752668*pow(theta,3.2);
//	rho = exp(RHS)*reduce.rho;
//	return rho;
//}
//
//double OrthoHydrogenClass::conductivity_Trho(double T, double rho)
//{
//	double e_k = 59.7,
//		   sigma = 0.2827;
//	double eta_dilute = viscosity_dilute(T,e_k,sigma);
//	return 15e-3*R()*(eta_dilute*1e6)/4.0/1000;
//}
//double OrthoHydrogenClass::viscosity_Trho(double T, double rho)
//{
//	double A,B,eta_0,eta_E,e_k,sigma;
//
//    if (T < 100){
//        // Dilute gas contribution
//        eta_0 = 8.5558 * (pow(T,1.5)/(T+19.55)) * ((T+650.39)/(T+1175.9));
//
//		//For the excess part, density is in units of g/cm3, so need to divide by 1000
//		rho /= 1000;
//
//        A = (306.4636*rho - 3350.628*rho*rho+3868.092*rho*rho*rho)/(1.0-18.47830*rho+110.915*rho*rho+25.3542*rho*rho*rho);
//
//        B = 10.0 + 7.2*(pow(rho/0.07,6)-pow(rho/0.07,1.5))-17.63/exp(58.75*pow(rho/0.07,3));
//
//        // Excess viscosity
//        eta_E = A * exp(B/T);
//
//		// Correlation in units of g/cm-s x 10e-6, so to get Pa-s, need to divide by 10 and divide by 1e6
//		
//		return (eta_0 + eta_E) / 1e7;
//    }
//    else{
//		// Use dilute gas properties
//		e_k = 59.7;
//		sigma = 0.2827;
//		return viscosity_dilute(T,e_k,sigma);
//    }
//}
//
//double OrthoHydrogenClass::surface_tension_T(double T)
//{
//	return 0.005369*pow(1-T/reduce.T,1.065);
//}