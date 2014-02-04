/*

Properties of Water (R718)
by Ian Bell

Themo properties from 
"The IAPWS Formulation 1995 for the Thermodynamic Properties
of Ordinary Water Substance for General and Scientific Use", 
W. Wagner and A. Pruss, J. Phys. Chem. Ref. Data, v. 31, 2000

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
#include <vector>
#include "FluidClass.h"
#include "Water.h"

WaterClass::WaterClass()
{
	static const double n[]={0,    
		 0.125335479355233e-1, //[1]
		 0.78957634722828e1, //[2]
		-0.87803203303561e1, //[3]
		 0.31802509345418, //[4]
		-0.26145533859358, //[5]
		-0.78199751687981e-2, //[6]
		 0.88089493102134e-2, //[7]
		-0.66856572307965, //[8]
		 0.20433810950965, //[9]
		-0.66212605039687e-4, //[10]
		-0.19232721156002, //[11]
		-0.25709043003438, //[12]
		 0.16074868486251, //[13]
		-0.40092828925807e-1, //[14]
		 0.39343422603254e-6, //[15]
		-0.75941377088144e-5, //[16]
		 0.56250979351888e-3, //[17]
		-0.15608652257135e-4, //[18]
		 0.11537996422951e-8, //[19]
		 0.36582165144204e-6, //[20]
		-0.13251180074668e-11, //[21]
		-0.62639586912454e-9, //[22]
		-0.10793600908932, //[23]
		 0.17611491008752e-1, //[24]
		 0.22132295167546, //[25]
		-0.40247669763528, //[26]
		 0.58083399985759, //[27]
		 0.49969146990806e-2, //[28]
		-0.31358700712549e-1, //[29]
		-0.74315929710341, //[30]
		 0.47807329915480, //[31]
		 0.20527940895948e-1, //[32]
		-0.13636435110343, //[33]
		 0.14180634400617e-1, //[34]
		 0.83326504880713e-2, //[35]
		-0.29052336009585e-1, //[36]
		 0.38615085574206e-1, //[37]
		-0.20393486513704e-1, //[38]
		-0.16554050063734e-2, //[39]
		 0.19955571979541e-2, //[40]
		 0.15870308324157e-3, //[41]
		-0.16388568342530e-4, //[42]
		 0.43613615723811e-1, //[43]
		 0.34994005463765e-1, //[44]
		-0.76788197844621e-1, //[45]
		 0.22446277332006e-1, //[46]
		-0.62689710414685e-4, //[47]
		-0.55711118565645e-9, //[48]
		-0.19905718354408, //[49]
		 0.31777497330738, //[50]
		-0.11841182425981, //[51]
		-0.31306260323435e2, //[52]
		 0.31546140237781e2, //[53]
		-0.25213154341695e4, //[54]
		-0.14874640856724, //[55]
		 0.31806110878444, //[56]
	};

	static const int d[]={0,
		1, //[1]
		1, //[2]
		1, //[3]
		2, //[4]
		2, //[5]
		3, //[6]
		4, //[7]
		1, //[8]
		1, //[9]
		1, //[10]
		2, //[11]
		2, //[12]
		3, //[13]
		4, //[14]
		4, //[15]
		5, //[16]
		7, //[17]
		9, //[18]
		10, //[19]
		11, //[20]
		13, //[21]
		15, //[22]
		1, //[23]
		2, //[24]
		2, //[25]
		2, //[26]
		3, //[27]
		4, //[28]
		4, //[29]
		4, //[30]
		5, //[31]
		6, //[32]
		6, //[33]
		7, //[34]
		9, //[35]
		9, //[36]
		9, //[37]
		9, //[38]
		9, //[39]
		10, //[40]
		10, //[41]
		12, //[42]
		3, //[43]
		4, //[44]
		4, //[45]
		5, //[46]
		14, //[47]
		3, //[48]
		6, //[49]
		6, //[50]
		6, //[51]
		3, //[52]
		3, //[53]
		3, //[54]
	};

	static const double t[]={0.00,
		-0.5, //[1]
		0.875, //[2]
		1, //[3]
		0.5, //[4]
		0.75, //[5]
		0.375, //[6]
		1, //[7]
		4, //[8]
		6, //[9]
		12, //[10]
		1, //[11]
		5, //[12]
		4, //[13]
		2, //[14]
		13, //[15]
		9, //[16]
		3, //[17]
		4, //[18]
		11, //[19]
		4, //[20]
		13, //[21]
		1, //[22]
		7, //[23]
		1, //[24]
		9, //[25]
		10, //[26]
		10, //[27]
		3, //[28]
		7, //[29]
		10, //[30]
		10, //[31]
		6, //[32]
		10, //[33]
		10, //[34]
		1, //[35]
		2, //[36]
		3, //[37]
		4, //[38]
		8, //[39]
		6, //[40]
		9, //[41]
		8, //[42]
		16, //[43]
		22, //[44]
		23, //[45]
		23, //[46]
		10, //[47]
		50, //[48]
		44, //[49]
		46, //[50]
		50, //[51]
		0, //[52]
		1, //[53]
		4, //[54]
	};

	static const int c[]={0,
		0, //[1]
		0, //[2]
		0, //[3]
		0, //[4]
		0, //[5]
		0, //[6]
		0, //[7]
		1, //[8]
		1, //[9]
		1, //[10]
		1, //[11]
		1, //[12]
		1, //[13]
		1, //[14]
		1, //[15]
		1, //[16]
		1, //[17]
		1, //[18]
		1, //[19]
		1, //[20]
		1, //[21]
		1, //[22]
		2, //[23]
		2, //[24]
		2, //[25]
		2, //[26]
		2, //[27]
		2, //[28]
		2, //[29]
		2, //[30]
		2, //[31]
		2, //[32]
		2, //[33]
		2, //[34]
		2, //[35]
		2, //[36]
		2, //[37]
		2, //[38]
		2, //[39]
		2, //[40]
		2, //[41]
		2, //[42]
		3, //[43]
		3, //[44]
		3, //[45]
		3, //[46]
		4, //[47]
		6, //[48]
		6, //[49]
		6, //[50]
		6, //[51]
		0, //[52]
		0, //[53]
		0, //[54]
	};

	static const double alpha[]={
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 51
		20, //[52]
		20, //[53]
		20, //[54]
	};

	static const double beta[]={
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 51
		150, //[52]
		150, //[53]
		250, //[54]
		0.3, //[55]
		0.3, //[56]
	};

	static const double GAMMA[]={
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 51
		1.21, //[52]
		1.21, //[53]
		1.25, //[54]
	};

	static const double epsilon[]={
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 51
		1, //[52]
		1, //[53]
		1, //[54]
	};

	static const double a[]={
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 54
		3.5, //[55]
		3.5, //[56]
	};

	static const double b[]={
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 54
		0.85, //[55]
		0.95, //[56]
	};

	static const double A[]={
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 54
		0.32, //[55]
		0.32, //[56]
	};

	static const double B[]={
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 54
		0.2, //[55]
		0.2, //[56]
	};

	static const double C[]={
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 54
		28., //[55]
		32., //[56]
	};

	static const double D[]={
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 54
		700.0, //[55]
		800.0, //[56]
	};

	// Constant for ideal gas part
	static const double n0[]={0,
		-8.3204464837497, //[1] //From Herrmann 2099 ASHRAE RP-1485
		6.6832105275932, //[2] //From Herrmann 2099 ASHRAE RP-1485
		3.00632, //[3]
		0.012436, //[4]
		0.97315, //[5]
		1.27950, //[6]
		0.96956, //[7]
		0.24873, //[8]
	};
	static const double gamma0[]={0,
		0,0,0, // 1 to 3
		1.28728967, //[4]
		3.53734222, //[5]
		7.74073708, //[6]
		9.24437796, //[7]
		27.5075105 //[8]
	};
	std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	std::vector<double> d_v(d,d+sizeof(d)/sizeof(int));
	std::vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	std::vector<double> l_v(c,c+sizeof(c)/sizeof(int));
	std::vector<double> alpha_v(alpha,alpha+sizeof(alpha)/sizeof(double));
	std::vector<double> beta_v(beta,beta+sizeof(beta)/sizeof(double));
	std::vector<double> gamma_v(GAMMA,GAMMA+sizeof(GAMMA)/sizeof(double));
	std::vector<double> epsilon_v(epsilon,epsilon+sizeof(epsilon)/sizeof(double));
	std::vector<double> a_v(a,a+sizeof(a)/sizeof(double));
	std::vector<double> b_v(b,b+sizeof(b)/sizeof(double));
	std::vector<double> A_v(A,A+sizeof(A)/sizeof(double));
	std::vector<double> B_v(B,B+sizeof(B)/sizeof(double));
	std::vector<double> C_v(C,C+sizeof(C)/sizeof(double));
	std::vector<double> D_v(D,D+sizeof(D)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));
	std::vector<double> gamma0_v(gamma0,gamma0+sizeof(gamma0)/sizeof(double));

	phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,51);
	phi_BC * phirg_ = new phir_gaussian(n_v,d_v,t_v,alpha_v,epsilon_v,beta_v,gamma_v,52,54);
	phi_BC * phirc_ = new phir_critical(n_v,d_v,t_v,a_v,b_v,beta_v,A_v,B_v,C_v,D_v,55,56);

	phirlist.push_back(phir_);
	phirlist.push_back(phirg_);
	phirlist.push_back(phirc_);

	phi_BC * phi0_lead_ = new phi0_lead(n0[1],n0[2]);
	phi_BC * phi0_logtau_ = new phi0_logtau(n0[3]);
	phi_BC * phi0_PE_ = new phi0_Planck_Einstein(n0_v,gamma0_v,4,8);
	
	phi0list.push_back(phi0_lead_);
	phi0list.push_back(phi0_logtau_);
	phi0list.push_back(phi0_PE_);

	// Critical parameters
	crit.rho = 322;
	crit.p = PressureUnit(22064, UNIT_KPA);
	crit.T = 647.096;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 18.015268;
	params.Ttriple = 273.16;
	params.ptriple = 0.611699223742;
	params.R_u = 8.314371357587;
	params.accentricfactor = 0.3442920843;

	// Fluid limits
	limits.Tmin = params.Ttriple;
	limits.Tmax = 2000;
	limits.pmax = 1000000.0;
	limits.rhomax = 73.96 * params.molemass;
	
	EOSReference.assign("\"The IAPWS Formulation 1995 for the Thermodynamic Properties"
						"of Ordinary Water Substance for General and Scientific Use\"," 
						"W. Wagner and A. Pruss, J. Phys. Chem. Ref. Data, v. 31, 2002");
		
	TransportReference.assign("Thermal Conductivity: Release on the IAPWS Formulation 2011 for the Thermal Conductivity of Ordinary Water Substance\n\n"
		"Viscosity: Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance\n\n"
		"Surface Tension: International Representation of the Surface Tension of Ordinary Water Substance 1994\n");

	name.assign("Water");
	aliases.push_back("water");
	aliases.push_back(std::string("WATER"));
	aliases.push_back("H2O");
	aliases.push_back("h2o");

	BibTeXKeys.EOS = "Wagner-JPCRD-2002";
	BibTeXKeys.VISCOSITY = "Huber-JPCRD-2009";
	BibTeXKeys.CONDUCTIVITY = "Huber-JPCRD-2012";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.ECS_LENNARD_JONES = "Poling-BOOK-2001";
}
static void visc_Helper(double Tbar, double rhobar, double *mubar_0, double *mubar_1)
{
	double H[6][7],sum;
	int i,j;

	// Dilute-gas component
	*mubar_0=100.0*sqrt(Tbar)/(1.67752+2.20462/Tbar+0.6366564/powInt(Tbar,2)-0.241605/powInt(Tbar,3));

	//Fill in zeros in H
	for (i=0;i<=5;i++)
	{
		for (j=0;j<=6;j++)
		{
			H[i][j]=0;
		}
	}

	//Set non-zero parameters of H
	H[0][0]=5.20094e-1;
	H[1][0]=8.50895e-2;
	H[2][0]=-1.08374;
	H[3][0]=-2.89555e-1;

	H[0][1]=2.22531e-1;
	H[1][1]=9.99115e-1;
	H[2][1]=1.88797;
	H[3][1]=1.26613;
	H[5][1]=1.20573e-1;

	H[0][2]=-2.81378e-1;
	H[1][2]=-9.06851e-1;
	H[2][2]=-7.72479e-1;
	H[3][2]=-4.89837e-1;
	H[4][2]=-2.57040e-1;

	H[0][3]=1.61913e-1;
	H[1][3]=2.57399e-1;

	H[0][4]=-3.25372e-2;
	H[3][4]=6.98452e-2;

	H[4][5]=8.72102e-3;

	H[3][6]=-4.35673e-3;
	H[5][6]=-5.93264e-4;

	// Finite density component
	sum=0;
	for (i=0;i<=5;i++)
	{
		for (j=0;j<=6;j++)
		{
			sum+=powInt(1/Tbar-1,i)*(H[i][j]*powInt(rhobar-1,j));
		}
	}
	*mubar_1=exp(rhobar*sum);
}
double WaterClass::conductivity_Trho(double T, double rho)
{
	// Many thanks to http://twt.mpei.ac.ru/mcs/worksheets/iapws/wspsTCRT.xmcd

	double L[5][6] = {{1.60397357,-0.646013523,0.111443906,0.102997357,-0.0504123634,0.00609859258},
				{2.33771842,-2.78843778,1.53616167,-0.463045512,0.0832827019,-0.00719201245},
				{2.19650529,-4.54580785,3.55777244,-1.40944978,0.275418278,-0.0205938816},
				{-1.21051378,1.60812989,-0.621178141,0.0716373224,0,0},
				{-2.7203370,4.57586331,-3.18369245,1.1168348,-0.19268305,0.012913842}};

	double lambdabar_0,lambdabar_1,lambdabar_2,rhobar,Tbar,sum,R_Water;
	double Tstar=647.096,rhostar=322,pstar=22064000,lambdastar=1e-3,mustar=1e-6;
	double tau,xi;
	int i,j;

	Tbar = T/Tstar;
	rhobar = rho/rhostar;
	R_Water = 8314.47215/params.molemass;

	// Dilute gas contribution
	lambdabar_0=sqrt(Tbar)/(2.443221e-3+1.323095e-2/Tbar+6.770357e-3/pow(Tbar,2)-3.454586e-3/pow(Tbar,3)+4.096266e-4/pow(Tbar,4));

	sum=0;
	for (i=0;i<=4;i++)
	{
		for (j=0;j<=5;j++)
		{
			//printf("L[%d][%d]: %0.12fg\n",i,j,L[i][j]);
			sum+=L[i][j]*powInt(1.0/Tbar-1.0,i)*powInt(rhobar-1,j);
		}
	}
	// Finite density contribution
	lambdabar_1=exp(rhobar*sum);

	double nu=0.630,GAMMA =177.8514,gamma=1.239,xi_0=0.13,Lambda_0=0.06,Tr_bar=1.5,qd_bar=1/0.4,pi=3.141592654,
		R=461.51805;//J/kg/K
	rhobar=rho/rhostar;
	Tbar=T/Tstar;
	tau=1/Tbar;
	
	double drhodp=1/(R*T*(1+2*rhobar*dphir_dDelta(1/Tbar,rhobar)+rhobar*rhobar*d2phir_dDelta2(1/Tbar,rhobar)));
	double drhobar_dpbar=pstar/rhostar*drhodp;
	double drhodp_Trbar=1/(R*Tr_bar*Tstar*(1+2*rhobar*dphir_dDelta(1/Tr_bar,rhobar)+rhobar*rhobar*d2phir_dDelta2(1/Tr_bar,rhobar)));
	double drhobar_dpbar_Trbar=pstar/rhostar*drhodp_Trbar;
	double cp=specific_heat_p_Trho(T,rho); // [J/kg/K]
	double cv=specific_heat_v_Trho(T,rho); // [J/kg/K]
	double cpbar=cp/R; //[-]
	double mubar = viscosity_Trho(T,rho)/mustar;
	double DELTAchibar_T=rhobar*(drhobar_dpbar-drhobar_dpbar_Trbar*Tr_bar/Tbar);
	if (DELTAchibar_T<0)
		xi=0;
	else
		xi=xi_0*pow(DELTAchibar_T/Lambda_0,nu/gamma);
	double y=qd_bar*xi;

	double Z;
	double kappa = cp/cv;
	if (y<1.2e-7)
		Z=0;
	else
		Z=2/(pi*y)*(((1-1/kappa)*atan(y)+y/kappa)-(1-exp(-1/(1/y+y*y/3/rhobar/rhobar))));

	lambdabar_2=GAMMA*rhobar*cpbar*Tbar/mubar*Z;
	return (lambdabar_0*lambdabar_1+lambdabar_2)*lambdastar;
}
double WaterClass::viscosity_Trho(double T, double rho)
{
	double x_mu=0.068,qc=1/1.9,qd=1/1.1,nu=0.630,gamma=1.239,zeta_0=0.13,LAMBDA_0=0.06,Tbar_R=1.5;
	double delta,tau,mubar_0,mubar_1,mubar_2,drhodp,drhodp_R,DeltaChibar,zeta,w,L,Y,psi_D,Tbar,rhobar;
	double drhobar_dpbar,drhobar_dpbar_R,R_Water;
	
	Tbar=T/reduce.T;
	rhobar=rho/reduce.rho;
	R_Water=R();

	// Dilute and finite gas portions
	visc_Helper(Tbar,rhobar,&mubar_0,&mubar_1);

	///************************ Critical Enhancement ************************
	delta=rhobar;
	// "Normal" calculation
	tau=1/Tbar;
	drhodp=1/(R_Water*T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta)));
	drhobar_dpbar=reduce.p.Pa/reduce.rho*drhodp;
	// "Reducing" calculation
	tau=1/Tbar_R;
	drhodp_R=1/(R_Water*Tbar_R*reduce.T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta)));
	drhobar_dpbar_R=reduce.p.Pa/reduce.rho*drhodp_R;
	
	DeltaChibar=rhobar*(drhobar_dpbar-drhobar_dpbar_R*Tbar_R/Tbar);
	if (DeltaChibar<0)
		DeltaChibar=0;
	zeta=zeta_0*pow(DeltaChibar/LAMBDA_0,nu/gamma);
	if (zeta<0.3817016416)
	{
		Y=1.0/5.0*qc*zeta*powInt(qd*zeta,5)*(1-qc*zeta+powInt(qc*zeta,2)-765.0/504.0*powInt(qd*zeta,2));
	}
	else
	{
		psi_D=acos(pow(1+powInt(qd*zeta,2),-1.0/2.0));
		w=sqrt(fabs((qc*zeta-1)/(qc*zeta+1)))*tan(psi_D/2.0);
		if (qc*zeta>1)
		{
			L=log((1+w)/(1-w));
		}
		else
		{
			L=2*atan(fabs(w));
		}
		Y=1.0/12.0*sin(3*psi_D)-1/(4*qc*zeta)*sin(2*psi_D)+1.0/powInt(qc*zeta,2)*(1-5.0/4.0*powInt(qc*zeta,2))*sin(psi_D)-1.0/powInt(qc*zeta,3)*((1-3.0/2.0*powInt(qc*zeta,2))*psi_D-pow(fabs(powInt(qc*zeta,2)-1),3.0/2.0)*L);
	}
	mubar_2=exp(x_mu*Y);

	return (mubar_0*mubar_1*mubar_2)/1e6;
}
double WaterClass::psat(double T)
{
	const double ti[]={0,1.0,1.5,3.0,3.5,4.0,7.5};
    const double ai[]={0,-7.85951783,1.84408259,-11.7866497,22.6807411,-15.9618719,1.80122502};
    double summer=0;
    int i;
    for (i=1;i<=6;i++)
    {
        summer=summer+ai[i]*pow(1-T/reduce.T,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double WaterClass::rhosatL(double T)
{
	const double ti[]={0,1.0/3.0,2.0/3.0,5.0/3.0,16.0/3.0,43.0/3.0,110.0/3.0};
    const double ai[]={0,1.99274064,1.09965342,-0.510839303,-1.75493479,-45.5170352,-1.75493479,-45.5170352,-6.74694450e5};
    double summer=1;
    int i;
    for (i=1;i<=6;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/reduce.T,ti[i]);
    }
    return reduce.rho*summer;
}
double WaterClass::rhosatV(double T)
{
	const double ti[]={0,2.0/6.0,4.0/6.0,8./6.0,18.0/6.0,37.0/6.0,71.0/6.0};
    const double ai[]={0,-2.03150240,-2.68302940,-5.38626492,-17.2991605,-44.7586581,-63.9201063};
    double summer=0;
    int i;
    for (i=1;i<=6;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/reduce.T,ti[i]);
    }
    return reduce.rho*exp(summer);
}

double WaterClass::surface_tension_T(double T)
{
	// From Mulero, 2012
	return -0.1306*pow(1-T/reduce.T,2.471)+0.2151*pow(1-T/reduce.T,1.233);
}

//double WaterClass::IsothermCompress_Water(double T, double p)
//{
//    // Isothermal compressibility with units of 1/Pa
//    double rho,delta,tau,dpdrho,R_Water;
//    rho=Props('D','T',T,'P',p,"Water");
//    delta=rho/rhoc;
//    tau=Tc/T;
//    R_Water=8.31447215/M_Water;
//    dpdrho=R_Water*T*(1.0+2.0*delta*dphir_dDelta(tau,delta)+powInt(delta,2)*d2phir_dDelta2(tau,delta));
//    return 1/(rho*dpdrho*1000.0);
//}
//
//double WaterClass::B_Water(double tau)
//{
//	// given by B*rhoc=lim(delta --> 0) [dphir_ddelta(tau)]
//	return 1.0/rhoc*dphir_dDelta(tau,1e-12);
//}
//
//double WaterClass::dBdT_Water(double tau)
//{
//	// given by B*rhoc^2=lim(delta --> 0) [dphir2_ddelta2(tau)]
//	return -1.0/rhoc*tau*tau/Tc*d2phir_dDelta_dTau(tau,1e-12);
//}
//
//double WaterClass::C_Water(double tau)
//{
//	// given by B*rhoc^2=lim(delta --> 0) [dphir2_ddelta2(tau)]
//	return 1.0/(rhoc*rhoc)*d2phir_dDelta2(tau,1e-12);
//}
//
//double WaterClass::dCdT_Water(double tau)
//{
//	// given by B*rhoc^2=lim(delta --> 0) [dphir2_ddelta2(tau)]
//	return -1.0/(rhoc*rhoc)*tau*tau/Tc*d3phir_dDelta2_dTau(tau,1e-12);
//}

/*
double dphir3_dDelta2_dTau_Water(double tau, double delta)
{ 
    // Derivation from Hermann 2009 (ASHRAE Report RP-1485)
    int i;
    double di, ci;
    double dphir3_dDelta2_dTau=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta,psi,dPSI2_dDelta2,dDELTAbi2_dDelta2,dDELTA2_dDelta2;
    double dPSI2_dDelta_dTau, dDELTAbi2_dDelta_dTau, dPSI_dTau, dDELTAbi_dTau;
    double Line1,Line2,Line3,dDELTA2_dDelta_dTau,dDELTA3_dDelta2_dTau,dDELTAbim1_dTau,dDELTAbim2_dTau;
    double dDELTA_dTau,dDELTAbi3_dDelta2_dTau;
    for (i=1;i<=7;i++)
    {
        di=(double)d[i];
        dphir3_dDelta2_dTau+=n[i]*di*t[i]*(di-1)*powInt(delta,d[i]-2)*pow(tau,t[i]-1.0);
    }
    for (i=8;i<=51;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir3_dDelta2_dTau+=n[i]*t[i]*exp(-powInt(delta,c[i]))*powInt(delta,d[i]-2)*pow(tau,t[i]-1.0)*(((di-ci*powInt(delta,c[i])))*(di-1-c[i]*powInt(delta,c[i]))-c[i]*c[i]*powInt(delta,c[i]));
    }
    for (i=52;i<=54;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir3_dDelta2_dTau+=n[i]*pow(tau,t[i])*psi*powInt(delta,d[i])*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]))*(-2.0*alpha[i]*(delta-epsilon[i])+4*alpha[i]*alpha[i]*pow(delta,d[i])*powInt(delta-epsilon[i],2)-4*d[i]*alpha[i]*pow(delta,d[i]-1)*(delta-epsilon[i])+d[i]*(d[i]-1)*powInt(delta,d[i]-2));
    }
    for (i=55;i<=56;i++)
    {
        theta=(1.0-tau)+A[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=powInt(theta,2)+B[i]*pow(powInt(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powInt(delta-1.0,2)-D[i]*powInt(tau-1.0,2));
        
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(powInt(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
        dPSI2_dDelta2=(2.0*C[i]*powInt(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
        dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+powInt(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(powInt(delta-1.0,2),a[i]-2.0)+2.0*powInt(A[i]/beta[i],2)*powInt(pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
        dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*powInt(dDELTA_dDelta,2));
        
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        
        dPSI2_dDelta_dTau=4.0*C[i]*D[i]*(delta-1.0)*(tau-1.0)*PSI;
        dDELTAbi2_dDelta_dTau=-A[i]*b[i]*2.0/beta[i]*pow(DELTA,b[i]-1.0)*(delta-1.0)*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0)-2.0*theta*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0)*dDELTA_dDelta;
	
		//Following Terms added for this derivative
		dDELTA_dTau=-2*((1-tau)+A[i]*pow(powInt(delta-1,2),1/(2*beta[i])-1)+2*B[i]*a[i]*pow(powInt(delta-1,2),a[i]-1));
		dDELTA2_dDelta_dTau=-(delta-1)*A[i]*(2/beta[i])*pow(powInt(delta-1,2),1/(2*beta[i])-1);
		dDELTA3_dDelta2_dTau=1/(delta-1)*dDELTA2_dDelta_dTau-powInt(delta-1,2)*A[i]*(4/beta[i])*(1/(2*beta[i])-1)*pow(powInt(delta-1,2),1/(2*beta[i])-2);
		
		dDELTAbim1_dTau=(b[i]-1)*pow(DELTA,b[i]-2)*dDELTA_dTau;
		dDELTAbim2_dTau=(b[i]-2)*pow(DELTA,b[i]-3)*dDELTA_dTau;
		Line1=dDELTAbim1_dTau*dDELTA2_dDelta2+pow(DELTA,b[i]-1)*dDELTA3_dDelta2_dTau;
		Line2=(b[i]-1)*(dDELTAbim2_dTau*powInt(dDELTA_dDelta,2)+pow(DELTA,b[i]-2)*2*dDELTA2_dDelta_dTau*dDELTA_dDelta);
		dDELTAbi3_dDelta2_dTau=b[i]*(Line1+Line2);
		
		Line1=pow(DELTA,b[i])*(2*delta*dPSI2_dDelta_dTau+delta*dDELTA3_dDelta2_dTau)+dDELTAbi_dTau*(2*dPSI_dDelta+delta*dPSI2_dDelta2);
		Line2=2*dDELTAbi2_dDelta_dTau*(PSI+delta*dPSI_dDelta)+2*dDELTAbi_dDelta*(dPSI_dTau+delta*dPSI2_dDelta_dTau);
		Line3=dDELTAbi3_dDelta2_dTau*delta*PSI+dDELTAbi2_dDelta2*delta*dPSI_dTau;
        dphir3_dDelta2_dTau+=n[i]*(Line1+Line2+Line3);
    }
    return dphir3_dDelta2_dTau;
}
*/
