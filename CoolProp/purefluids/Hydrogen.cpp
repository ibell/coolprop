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

HydrogenClass::HydrogenClass()
{
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

	static const double d[]={0,
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

	static const double c[]={
	0,0,0,0,0,0,0,0, // indices [0-7]
	1,
	1
	};

	// alpha instead of eta is used here for consistency with the definitions in R744.c upon which R290.c is based
	static const double alpha[]={ // phi from paper
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

	static const double epsilon[]={ // D from paper
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

	phirlist.push_back(new phir_power(n,d,t,c,1,9,10));
	phirlist.push_back(new phir_gaussian(n,d,t,alpha,epsilon,beta,GAMMA,10,14,15));

	/* phi0=log(delta)+1.5*log(tau)+a0[1]+a0[2]*tau
        +a0[3]*log(1-exp(-b0[3]*tau))
        +a0[4]*log(1-exp(-b0[4]*tau))
        +a0[5]*log(1-exp(-b0[5]*tau))
        +a0[6]*log(1-exp(-b0[6]*tau))
        +a0[7]*log(1-exp(-b0[7]*tau));
	*/

	//lead term of the form log(delta)+a1+a2*tau
	phi0list.push_back(new phi0_lead(a0[1],a0[2]));
	phi0list.push_back(new phi0_logtau(1.5));
	phi0list.push_back(new phi0_Planck_Einstein(a0,b0,3,7,8));

	// Critical parameters
	crit.rho = 15.508*2.01588;
	crit.p = PressureUnit(1296.4, UNIT_KPA);
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
	TransportReference.assign("Conductivity: Assael, JPCRD, 2011");

	name.assign("Hydrogen");
	aliases.push_back("hydrogen");
	aliases.push_back("HYDROGEN");
	aliases.push_back("H2");
	aliases.push_back("R702");
	REFPROPname.assign("hydrogen");

	BibTeXKeys.EOS = "Leachman-JPCRD-2009";
	BibTeXKeys.VISCOSITY = "Muzny-JCED-2013";
	BibTeXKeys.CONDUCTIVITY = "Assael-JPCRD-2011";
	BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
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
    return reduce.p.Pa*exp(reduce.T/T*summer);
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
	// Max error is  0.213338271447 % between 13.957 and 33.144999 K
    const double t[] = {0, 0.3525, 0.3565, 0.367, 0.3765, 0.39249999999999996, 1.3333333333333333};
	const double N[] = {0, -426795.37901103671, 810857.2657358282, -721220.63750821573, 392092.37071915239, -54938.300252502217, 1.9072707235406241};
    double summer=0,theta;
    theta=1-T/reduce.T;    	
	for (int i=1; i<=6; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

double HydrogenClass::conductivity_Trho(double T, double rho)
{
	double sumnum = 0, sumden = 0, sumresid = 0;
	double A1[] = {-3.40976e-1, 4.58820e0, -1.45080e0, 3.26394e-1, 3.16939e-3, 1.90592e-4, -1.13900e-6};
	double A2[] = {1.38497e2, -2.21878e1, 4.57151e0, 1.00000e0};
	double B1[] = {0, 3.63081e-2, -2.07629e-2, 3.14810e-2, -1.43097e-2, 1.74980e-3};
	double B2[] = {0, 1.83370e-3, -8.86716e-3, 1.58260e-2, -1.06283e-2, 2.80673e-3};

	for (int i = 0; i<= 6; i++){		
		sumnum += A1[i]*pow(T/reduce.T,i);
	}
	for (int i = 0; i<= 3; i++){		
		sumden += A2[i]*pow(T/reduce.T,i);
	}

	double lambda_0 = sumnum/sumden; // [W/m/K]

	for (int i = 1; i<= 5; i++){		
		sumresid += (B1[i]+B2[i]*(T/reduce.T))*pow(rho/reduce.rho,i);
	}

	double lambda_r = sumresid; // [W/m/K]

	double lambda_c = this->conductivity_critical(T,rho,1.0/(4.0e-10)); // [W/m/K]

	return lambda_0+lambda_r+lambda_c;
}
double HydrogenClass::viscosity_Trho(double T, double rho)
{
	// Note, ECS L-J constants do not agree with Poling, 2001
	double sigma = 0.297, e_k = 30.41;
	
	// Dilute gas
	this->ECSParams(&e_k, &sigma);
	double Tstar = T/e_k;
	double a[] = {2.09630e-1, -4.55274e-1, 1.43602e-1, -3.35325e-2, 2.76981e-3};
	double S_star = exp(a[0]+a[1]*log(Tstar)+a[2]*log(Tstar)*log(Tstar)+a[3]*pow(log(Tstar),3)+a[4]*pow(log(Tstar),4));
	double lambda_0 = 0.021357*sqrt(params.molemass*T)/(sigma*sigma*S_star); //[uPa-s]

	// Initial-density
	double b[] = {-0.1870, 2.4871, 3.7151, -11.0972, 9.0965, -3.8292, 0.5166};
	double sumBstar = 0;
	for (int i = 0; i<= 6; i++){ sumBstar += b[i]/Tstar; }
	double Bstar = sumBstar;
	double N_A = 6.02214129e23;
	double B = N_A*pow(sigma/1e9,3)*Bstar;

	// Residual
	double c[] = {0, 6.43449673, 4.56334068e-2, 2.32797868e-1, 9.58326120e-1, 1.27941189e-1, 3.63576595e-1};
	double Tr = T/crit.T;
	double rhor = rho/90.5; 
	double lambda_r = c[1]*rhor*rhor*exp(c[2]*Tr+c[3]/Tr+(c[4]*rhor*rhor)/(c[5]+Tr)+c[6]*pow(rhor,6));

	double rhobar = rho/params.molemass*1000;
	return (lambda_0*(1+B*rhobar)+lambda_r)/1e6;

}

double HydrogenClass::surface_tension_T(double T)
{
	return -1.4165*pow(1-T/reduce.T,0.63882)+0.746383*pow(1-T/reduce.T,0.659804)+0.675625*pow(1-T/reduce.T,0.619149);
}

ParaHydrogenClass::ParaHydrogenClass()
{
	double n [] = {0, -7.33375, 0.01, 2.60375, 4.66279, 0.68239, -1.47078, 0.135801, -1.05327, 0.328239, -0.0577833, 0.0449743, 0.0703464, -0.0401766, 0.11951};
	double t [] = {0, 0.6855, 1, 1, 0.489, 0.774, 1.133, 1.386, 1.619, 1.162, 3.96, 5.276, 0.99, 6.791, 3.19};
	double d [] = {0, 1, 4, 1, 1, 2, 2, 3, 1, 3, 2, 1, 3, 1, 1};
	double c [] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
	double a0 [] = {0, -1.4485891134, 1.884521239, 4.30256, 13.0289, -47.7365, 50.0013, -18.6261, 0.993973, 0.536078};
	double b0 [] = {0, 0, 0, 15.14967511472, 25.0925982148, 29.4735563787, 35.4059141417, 40.724998482, 163.7925799988, 309.2173173842};
	double alpha[] = {0,0,0,0,0,0,0,0,0,0,1.7437, 0.5516, 0.0634, 2.1341, 1.777}; // phi from paper
	double beta[] = {0,0,0,0,0,0,0,0,0,0,0.194, 0.2019, 0.0301, 0.2383, 0.3253};
	double GAMMA[] = {0,0,0,0,0,0,0,0,0,0,0.8048, 1.5248, 0.6648, 0.6832, 1.493};
	double epsilon[] = {0,0,0,0,0,0,0,0,0,0,1.5487, 0.1785, 1.28, 0.6319, 1.7104}; // D from paper

	phirlist.push_back(new phir_power(n,d,t,c,1,9,10));
	phirlist.push_back(new phir_gaussian(n,d,t,alpha,epsilon,beta,GAMMA,10,14,15));

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
	phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(a0,b0,3,9,10);

	phi0list.push_back(phi0_lead_);
	phi0list.push_back(phi0_logtau_);
	phi0list.push_back(phi0_Planck_Einstein_);

	// Critical parameters
	crit.rho = 15.538*2.01588;
	crit.p = PressureUnit(1285.8, UNIT_KPA);
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
	aliases.push_back("Parahydrogen");
	aliases.push_back("parahydrogen");
	aliases.push_back("PARAHYDROGEN");
	REFPROPname.assign("PARAHYD");

	BibTeXKeys.EOS = "Leachman-JPCRD-2009";
	BibTeXKeys.CONDUCTIVITY = "Assael-JPCRD-2011";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double ParaHydrogenClass::psat(double T)
{
    // Max error is  0.0524435114644 % between 13.8033 and 32.937999 K

    const double t[]={0, 0.132, 0.3605, 1.0, 1.5, 3.3333333333333335, 7.0};
    const double N[]={0, -0.0031611083814973629, 0.01606103512717301, -5.0101011461385303, 1.3458439473996564, 0.82353198183584131, -0.57502774437436288};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1;i<=6;i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double ParaHydrogenClass::rhosatL(double T)
{
    // Maximum absolute error is 0.112227 % between 13.803301 K and 32.937990 K
    const double t[] = {0, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333};
    const double N[] = {0, 32.575128830119098, -105.60522773310021, 145.38945707541046, -93.098356183908166, 22.448904889155894};
    double summer=0,theta;
    theta=1-T/reduce.T;
	for (int i=1; i<=5; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double ParaHydrogenClass::rhosatV(double T)
{
	// Max error is  0.213338271447 % between 13.957 and 33.144999 K
    const double t[] = {0, 0.3525, 0.3565, 0.367, 0.3765, 0.39249999999999996, 1.3333333333333333};
	const double N[] = {0, -426795.37901103671, 810857.2657358282, -721220.63750821573, 392092.37071915239, -54938.300252502217, 1.9072707235406241};
    double summer=0,theta;
    theta=1-T/reduce.T;    	
	for (int i=1; i<=6; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

double ParaHydrogenClass::conductivity_Trho(double T, double rho)
{
	double sumnum = 0, sumden = 0, sumresid = 0;
	double A1[] = {-1.24500e0, 3.10212e2, -3.31004e2, 2.46016e2, -6.57810e1, 1.08260e1, -5.19659e-1, 1.43979e-2};
	double A2[] = {1.42304e4, -1.93922e4, 1.58379e4, -4.81812e3, 7.28639e2, -3.57365e1, 1.00000e0};
	double B1[] = {0, 2.65975e-2, -1.33826e-3, 1.30219e-2, -5.67678e-3, -9.23380e-5};
	double B2[] = {0, -1.21727e-3, 3.66663e-3, 3.88715e-3, -9.21055e-3, 4.00723e-3};

	for (int i = 0; i<= 7; i++){		
		sumnum += A1[i]*pow(T/reduce.T,i);
	}
	for (int i = 0; i<= 6; i++){		
		sumden += A2[i]*pow(T/reduce.T,i);
	}

	double lambda_0 = sumnum/sumden; // [W/m/K]

	for (int i = 1; i<= 5; i++){		
		sumresid += (B1[i]+B2[i]*(T/reduce.T))*pow(rho/reduce.rho,i);
	}

	double lambda_r = sumresid; // [W/m/K]

	double lambda_c = this->conductivity_critical(T,rho,1.0/(5.0e-10)); // [W/m/K]

	return lambda_0+lambda_r+lambda_c;

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
	// Mulero, JPCRD 2012
	return 0.005314*pow(1-T/reduce.T,1.06);
}

OrthoHydrogenClass::OrthoHydrogenClass()
{
	double n [] = {0, -6.83148, 0.01, 2.11505, 4.38353, 0.211292, -1.00939, 0.142086, -0.87696, 0.804927, -0.710775, 0.0639688, 0.0710858, -0.087654, 0.647088};
	double t [] = {0, 0.7333, 1, 1.1372, 0.5136, 0.5638, 1.6248, 1.829, 2.404, 2.105, 4.1, 7.658, 1.259, 7.589, 3.946};
	double d [] = {0, 1, 4, 1, 1, 2, 2, 3, 1, 3, 2, 1, 3, 1, 1};
	double l [] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
	
    double alpha[] = {0,0,0,0,0,0,0,0,0,0,1.169, 0.894, 0.04, 2.072, 1.306}; // Originally phi in the paper
    double beta[] = {0,0,0,0,0,0,0,0,0,0,0.4555, 0.4046, 0.0869, 0.4415, 0.5743};
    double gamma[] = {0,0,0,0,0,0,0,0,0,0,1.5444, 0.6627, 0.763, 0.6587, 1.4327};
    double epsilon[] = {0,0,0,0,0,0,0,0,0,0,0.6366, 0.3876, 0.9437, 0.3976, 0.9626}; // Originally D in the paper	

	phirlist.push_back(new phir_power(n,d,t,l,1,9,15));
	phirlist.push_back(new phir_gaussian(n,d,t,alpha,epsilon,beta,gamma,10,14,15));

	/* phi0=log(delta)+1.5*log(tau)+a0[1]+a0[2]*tau
        +a0[3]*log(1-exp(-b0[3]*tau))
        +a0[4]*log(1-exp(-b0[4]*tau))
        +a0[5]*log(1-exp(-b0[5]*tau))
        +a0[6]*log(1-exp(-b0[6]*tau))
        +a0[7]*log(1-exp(-b0[7]*tau));
	*/

	double a0 [] = {0, -1.4675442336, 1.8845068862, 2.54151, -2.3661, 1.00365, 1.22447};
	double b0 [] = {0, 0, 0, 25.7676098736, 43.4677904877, 66.0445514750, 209.7531607465};

	//lead term of the form log(delta)+a1+a2*tau
	phi0list.push_back(new phi0_lead(a0[1],a0[2]));
	phi0list.push_back(new phi0_logtau(1.5));
	phi0list.push_back(new phi0_Planck_Einstein(a0,b0,3,6,7));

	// Critical parameters
	crit.rho = 15.445*2.01588;
	crit.p = PressureUnit(1310.65, UNIT_KPA);
	crit.T = 33.22;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 2.01594;
	params.Ttriple = 14.008;
	params.ptriple = 7.5598823410394012;
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

	name.assign("OrthoHydrogen");
	aliases.push_back("Orthohydrogen");
	aliases.push_back("orthohydrogen");
	aliases.push_back("ORTHOHYDROGEN");
	REFPROPname.assign("ORTHOHYD");

	BibTeXKeys.EOS = "Leachman-JPCRD-2009";
}

double OrthoHydrogenClass::psat(double T)
{
	// Max error is  0.053161754126 % between 14.008 and 33.219999 K
    const double ti[]={0, 0.3565, 0.381, 0.6666666666666666, 0.8333333333333334, 3.3333333333333335, 3.5};
    const double Ni[]={0, 0.23308544489395641, -0.44331497371946116, 2.596446168645369, -6.3537033133836278, 6.8116055471042287, -5.9839832788614595};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double OrthoHydrogenClass::rhosatL(double T)
{
    // Maximum absolute error is 0.898525 % between 14.008001 K and 33.219990 K
    const double ti[]={0,0.44474280444194492, 0.75146340805936485, 1.4887993676087643, 2.5830154027036798, 2.6777147409416382, 3.7537296061337075};
    const double Ni[]={0,2.5493992312992853, -2.168729482271587, 1.1414905046089598, -0.45527343384150787, -0.48960090300149495, 0.4904521362247205};
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
double OrthoHydrogenClass::rhosatV(double T)
{
    // Maximum absolute error is 0.558198 % between 14.008001 K and 33.219990 K
    const double ti[]={0,0.47329684677466771, 1.1913159863668648, 1.8991992477616062, 2.494090628975616, 3.6361782580222841, 24.83499518232826};
    const double Ni[]={0,-2.8790193675821163, -0.32069213937243635, -1.6289408828391703, 3.4288996368340423, -1.5495479021680494, -303.61663521542357};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(crit.T/T*summer);
}
