/*
Properties for R717 (ammonia).  
by Ian Bell

Thermo props from
"Eine neue Fundamentalgleichung fur Ammoniak (A new Equation of State for Ammonia)"
by R. Tillner-Roth and F. Harms-Watzenberg and H.D. Baehr, Deutscher Kaelte- und Klimatechnischer Verein Tagung 1993

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
#include "R717.h"

/*
From REFPROP documentation:
The original paper has a typographical error that shows a positive coefficient
instead of negative.  The correct value should be -0.3497111e-01.
*/

R717Class::R717Class()
{

	static const double n[]={
		 0.0,			//[0]
		 0.4554431e-1, 	//[1]
		 0.7238548e+0,	//[2]
		 0.1229470e-1,	//[3]
		-0.1858814e+1,	//[4]
		 0.2141882e-10,	//[5]
		-0.1430020e-1,	//[6]
		 0.3441324e+0,	//[7]
		-0.2873571e+0,	//[8]
		 0.2352589e-4,	//[9]
		-0.3497111e-1,	//[10]
		 0.2397852e-1,	//[11]
		 0.1831117e-2,	//[12]
		-0.4085375e-1,	//[13]
		 0.2379275e+0,	//[14]
 		-0.3548972e-1,	//[15]
		-0.1823729e+0,	//[16]
		 0.2281556e-1,	//[17]
		-0.6663444e-2,	//[18]
		-0.8847486e-2,	//[19]
		 0.2272635e-2,	//[20]
		-0.5588655e-3,	//[21]
	};

	static const double d[]={
		0,			//[0]
		2, 			//[1]
		1, 			//[2]
		4, 			//[3]
		1, 			//[4]
		15, 		//[5]
		3, 			//[6]
		3, 			//[7]
		1, 			//[8]
		8, 			//[9]
		2, 			//[10]
		1, 			//[11]
		8, 			//[12]
		1, 			//[13]
		2, 			//[14]
		3, 			//[15]
		2, 			//[16]
		4, 			//[17]
		3, 			//[18]
		1, 			//[19]
		2, 			//[20]
		4 			//[21]
	};

	static const double t[]={
		0.0,		//[0]
		-1.0/2.0,	//[1]
		1.0/2.0,	//[2]
		1.0,		//[3]
		3.0/2.0,	//[4]
		3.0,		//[5]
		0.0,		//[6]
		3.0,		//[7]
		4.0,		//[8]
		4.0,		//[9]
		5.0,		//[10]
		3.0,		//[11]
		5.0,		//[12]
		6.0, 		//[13]
		8.0,		//[14]
		8.0,		//[15]
		10.0,		//[16]
		10.0,		//[17]
		5.0,		//[18]
		15.0/2.0,	//[19]
		15.0,		//[20]
		30.0		//[21]
	};

	static const double c[]={
		0.0,		//[0]
		0.0,		//[1]
		0.0,		//[2]
		0.0,		//[3]
		0.0,		//[4]
		0.0,		//[5]
		1.0,		//[6]
		1.0,		//[7]
		1.0,		//[8]
		1.0,		//[9]
		1.0,		//[10]
		2.0,		//[11]
		2.0,		//[12]
		2.0, 		//[13]
		2.0,		//[14]
		2.0,		//[15]
		2.0,		//[16]
		2.0,		//[17]
		3.0,		//[18]
		3.0,		//[19]
		3.0,		//[20]
		3.0			//[21]
	};

	static const double a0[]={
		0.0,		//[0]
		-15.815020,	//[1]
		4.255726,	//[2]
		11.474340,	//[3]
		-1.296211,	//[4]
		0.5706757	//[5]
	};
	static const double t0[]={
		0.0,		//[0]
		0.0,		//[1]
		0.0,		//[2]
		1.0/3.0,	//[3]
		-3.0/2.0,	//[4]
		-7.0/4.0	//[5]
	};

	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> t0_v(t0,t0+sizeof(t0)/sizeof(double));

	phirlist.push_back(new phir_power(n,d,t,c,1,21,22));

	// phi0=log(delta)+a0[1]+a0[2]*tau-log(tau)+a0[3]*pow(tau,1.0/3.0)+a0[4]*pow(tau,-3.0/2.0)+a0[5]*pow(tau,-7.0/4.0);
	phi0list.push_back(new phi0_lead(a0[1],a0[2]));
	phi0list.push_back(new phi0_logtau(-1));
	phi0list.push_back(new phi0_power(a0_v,t0_v,3,5));

	// Critical parameters
	crit.rho = 225;
	crit.p = PressureUnit(11333, UNIT_KPA);
	crit.T = 405.40;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 17.03026;
	params.Ttriple = 195.495;
	params.ptriple = 6.09170982378;
	params.accentricfactor = 0.25601;
	params.R_u = 8.314471;

	// Limits of EOS
	limits.Tmin = 195.495;
	limits.Tmax = 700.0;
	limits.pmax = 1000000.0;
	limits.rhomax = 52.915*params.molemass;
	
	EOSReference.assign("\"Eine neue Fundamentalgleichung fur Ammoniak (A new Equation of State for Ammonia)\""
						" R. Tillner-Roth and F. Harms-Watzenberg and H.D. Baehr, "
						"Deutscher Kaelte- und Klimatechnischer Verein Tagung 1993");
	TransportReference.assign("Viscosity: \"The Viscosity of Ammonia\", "
							"A. Fenghour and W.A. Wakeham and V. Vesovic and J.T.R. Watson and J. Millat and E. Vogel"
							"J. Phys. Chem. Ref. Data, Vol. 24, No. 5, 1995 \n\n"
							"Conductivity: \"Thermal Conductivity of Ammonia in a Large"
							"Temperature and Pressure Range Including the Critical Region\""
							"by R. Tufeu, D.Y. Ivanov, Y. Garrabos, B. Le Neindre, "
							"Bereicht der Bunsengesellschaft Phys. Chem. 88 (1984) 422-427\n\n"
							"Does not include the critical enhancement.  Comparison of EES (without enhancement) and Refprop (with enhancement) give "
							"errors in saturated liquid and saturated vapor conductivities of \n\n"
							"T < 325K, error < 0.1%\n"
							"325K < T < 355 K, error <1%\n\n"
							"Most practical conditions will be in the <325K range\n"
							"Surface tension:  Michael Kleiber and Ralph Joh, \"VDI Heat Atlas 2010 Chapter D3.1 Liquids and Gases\" ");

	name.assign("Ammonia");
	aliases.push_back("NH3");
	aliases.push_back("ammonia");
	aliases.push_back("R717");
	aliases.push_back(std::string("AMMONIA"));
	REFPROPname.assign("AMMONIA");

	BibTeXKeys.EOS = "TillnerRoth-DKV-1993";
	BibTeXKeys.VISCOSITY = "Fenghour-JPCRD-1995";
	BibTeXKeys.CONDUCTIVITY = "Tufeu-BBPC-1984";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R717Class::conductivity_Trho(double T, double rho)
{
	/* 
	From "Thermal Conductivity of Ammonia in a Large 
	Temperature and Pressure Range Including the Critical Region"
	by R. Tufeu, D.Y. Ivanov, Y. Garrabos, B. Le Neindre, 
	Bereicht der Bunsengesellschaft Phys. Chem. 88 (1984) 422-427
	*/

	/* 
	Does not include the critical enhancement.  Comparison of EES (without enhancement) and Refprop (with enhancement) give 
	errors in saturated liquid and saturated vapor conductivities of 

	T < 325K, error < 0.1%
	325K < T < 355 K, error <1%

	Nearly all practical conditions will be in the <325K range
	*/


	double lambda_0,lambda_tilde;

	double a= 0.3589e-1;
	double b=-0.1750e-3;
	double c= 0.4551e-6;
	double d= 0.1685e-9;
	double e=-0.4828e-12;
	double lambda_1= 0.16207e-3;
	double lambda_2= 0.12038e-5;
	double lambda_3=-0.23139e-8;
	double lambda_4= 0.32749e-11;

	// This variable appears to be unused in the code:
	//double a_zeta_plus=0.7;

	double LAMBDA=1.2, nu=0.63, gamma =1.24, DELTA=0.50, rhoc_visc=235,t,zeta_0_plus=1.34e-10,a_zeta=1,GAMMA_0_plus=0.423e-8;
	double pi=3.141592654,a_chi,k_B=1.3806504e-23,X_T,DELTA_lambda,dPdT,eta_B,DELTA_lambda_id,DELTA_lambda_i;

	lambda_0=a+b*T+c*T*T+d*T*T*T+e*T*T*T*T;
	lambda_tilde=lambda_1*rho+lambda_2*rho*rho+lambda_3*rho*rho*rho+lambda_4*rho*rho*rho*rho;
	
	if (0)//(T>Tc)
	{
		t=fabs((T-reduce.T)/reduce.T);
		a_chi=a_zeta/0.7;
		eta_B=(2.60*1.6*t)*1e-5;
		dPdT=(2.18-0.12/exp(17.8*t))*1e5; //[Pa-K]
		X_T=0.61*rhoc_visc+16.5*log(t);
		// Along the critical isochore (only a function of temperature) (Eq. 9)
		DELTA_lambda_i=LAMBDA*(k_B*T*T)/(6*pi*eta_B*(zeta_0_plus*pow(t,-nu)*(1+a_zeta*pow(t,DELTA))))*dPdT*dPdT*GAMMA_0_plus*pow(t,-gamma)*(1+a_chi*pow(t,DELTA));
		DELTA_lambda_id=DELTA_lambda_i*exp(-36*t*t);
		if (rho<0.6*reduce.rho)
		{
			DELTA_lambda=DELTA_lambda_id*(X_T*X_T)/(X_T*X_T+powInt(0.6*rhoc_visc-0.96*rhoc_visc,2))*powInt(rho,2)/powInt(0.6*rhoc_visc,2);
		}
		else
		{
			DELTA_lambda=DELTA_lambda_id*(X_T*X_T)/(X_T*X_T+powInt(rho-0.96*rhoc_visc,2));
		}
	}
	else
	{
		DELTA_lambda=0.0;
	}

	return lambda_0+lambda_tilde+DELTA_lambda;
}
double R717Class::viscosity_Trho(double T, double rho)
{
	/* 
	From "The Viscosity of Ammonia"
	by A. Fenghour and W.A. Wakeham and V. Vesovic and J.T.R. Watson and J. Millat and E. Vogel
	J. Phys. Chem. Ref. Data, Vol. 24, No. 5, 1995 
	*/
	double sum=0, e_k=386.0,sigma=0.2957,M=17.03026,sum2=0.0;
	int i,j;
	double T_star,G_eta_star,eta0,B_eta_star,B_eta,b_1,deltaeta_h;
	double a[]={4.99318220,-0.61122364,0.0,0.18535124,-0.11160946};
	double c[]={-0.17999496e1,0.46692621e2,-0.53460794e3,0.33604074e4,-0.13019164e5,
		0.33414230e5,-0.58711743e5,0.71426686e5,-0.59834012e5,0.33652741e5,
		-0.12027350e5,0.24348205e4,-0.20807957e3};
	// indices are backwards from paper
	double d[5][3]={{0,0.17366936e-2,0.0},
	                {0.0,-0.64250359e-2,0.0},
	                {2.19664285e-1,0.0,1.67668649e-4},
	                {0.0,0.0,-1.49710093e-4},
	                {-0.83651107e-1,0.0,0.77012274e-4}};

	// rho is units of mol/L, so convert the density from kg/m^3 to mol/L (poorly documented in paper)
	rho=rho/M;

	sum=0;
	T_star=T/e_k;
	for (i=0;i<=4;i++)
	{
		sum+=a[i]*powInt(log(T_star),i);
	}
	G_eta_star=exp(sum);

	// From REFPROP fluid file: !=0.021357*SQRT(MW)*(unknown factor of 100)  [Chapman-Enskog term]
	// Seems like there is a typo in Fenghour - or am I missing something?
	eta0=0.021357*sqrt(T*M)*100/(sigma*sigma*G_eta_star);
	
	sum=0;
	for (i=0;i<=12;i++)
	{
		sum+=c[i]*powInt(sqrt(T_star),-i);
	}
	B_eta_star=sum;
	B_eta=B_eta_star*(0.6022137*sigma*sigma*sigma);
	b_1=B_eta*eta0;

	sum=0;
	for (i=2;i<=4;i++)
	{
		sum2=0.0;
		for (j=0;j<=4;j++)
		{
			// indices of d are backwards from paper
			sum2+=d[j][i-2]/powInt(T_star,j);
		}
		sum+=sum2*powInt(rho,i);
	}
	deltaeta_h=sum;
	return (eta0+b_1*rho+deltaeta_h)/1e6;
}
double R717Class::psat(double T)
{
    // Maximum absolute error is 0.055212 % between 195.495001 K and 405.399990 K
    const double t[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 26};
    const double N[]={0, 0.04909506379170004, -8.5985591284104057, 22.98850790970738, -222.28840209935356, 1491.9401420842644, -6641.0811040366289, 19697.275454679908, -38321.157974702161, 45841.962482923009, -27382.278742352592, 6743.9692067905371, -2077.8944819332091, 1272.2022925639299, -597.49829753755705};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=14;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double R717Class::rhosatL(double T)
{
    // Maximum absolute error is 0.264196 % between 195.495001 K and 405.399990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 2.1666666666666665};
    const double N[] = {0, 3.2591904910225704, -41.687339030964814, 246.43755287753763, -837.61113396671999, 2080.7279064276067, -3935.4817672564195, 5171.1096164456239, -4027.6487405218049, 1387.0671995039477, -43.611646611054738};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;	
	for (i=1; i<=10; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double R717Class::rhosatV(double T)
{
    // Maximum absolute error is 0.123990 % between 195.495001 K and 405.399990 K
    const double t[] = {0, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.8333333333333333};
    const double N[] = {0, -62.065663192057258, 869.56245921894947, -5338.5660401310852, 18273.366873938554, -38074.430040623774, 48944.180419243385, -36862.534086957254, 13256.30649363281, -1014.8176862674514};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
	for (i=1; i<=9; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}
double R717Class::surface_tension_T(double T)
{
	return 0.1028*pow(1-T/reduce.T,1.211)-0.09453*pow(1-T/reduce.T,5.585);
}
