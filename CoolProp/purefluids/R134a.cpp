/*
Properties for R134a.  
by Ian Bell

Thermo props from
"A International Standard Formulation for the Thermodynamic Properties of 1,1,1,2-Tetrafluoroethane 
(HFC-134a) for Temperatures from 170 K to 455 K and Pressures up to 70 MPa"
by Reiner Tillner-Roth and Hans Dieter Baehr, J. Phys. Chem. Ref. Data, v. 23, 1994, pp 657-729
*/

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include <cmath>
#include "string.h"
#include "stdio.h"
#include "CoolProp.h"
#include "FluidClass.h"
#include "R134a.h"

R134aClass::R134aClass()
{

	static const double n[]={
		 0.0,			//[0]
		 0.5586817e-1, 	//[1]
		 0.4982230e0,	//[2]
		 0.2458698e-1,	//[3]
		 0.8570145e-3,	//[4]
		 0.4788584e-3,	//[5]
		-0.1800808e1,	//[6]
		 0.2671641e0,	//[7]
		-0.4781652e-1,	//[8]
		 0.1423987e-1,	//[9]
		 0.3324062e0,	//[10]
		-0.7485907e-2,	//[11]
		 0.1017263e-3,	//[12]
		-0.5184567e+0,	//[13]
		-0.8692288e-1,	//[14]
 		 0.2057144e+0,	//[15]
		-0.5000457e-2,	//[16]
		 0.4603262e-3,	//[17]
		-0.3497836e-2,	//[18]
		 0.6995038e-2,	//[19]
		-0.1452184e-1,	//[20]
		-0.1285458e-3,	//[21]
	};

	static const double d[]={
		0,			//[0]
		2, 			//[1]
		1, 			//[2]
		3, 			//[3]
		6, 			//[4]
		6, 			//[5]
		1, 			//[6]
		1, 			//[7]
		2, 			//[8]
		5, 			//[9]
		2, 			//[10]
		2, 			//[11]
		4, 			//[12]
		1, 			//[13]
		4, 			//[14]
		1, 			//[15]
		2, 			//[16]
		4, 			//[17]
		1, 			//[18]
		5, 			//[19]
		3, 			//[20]
		10 			//[21]
	};

	static const double t[]={
		0.0,		//[0]
		-1.0/2.0,	//[1]
		0.0,		//[2]
		0.0,		//[3]
		0.0,		//[4]
		3.0/2.0,	//[5]
		3.0/2.0,	//[6]
		2.0,		//[7]
		2.0,		//[8]
		1.0,		//[9]
		3.0,		//[10]
		5.0,		//[11]
		1.0,		//[12]
		5.0, 		//[13]
		5.0,		//[14]
		6.0,		//[15]
		10.0,		//[16]
		10.0,		//[17]
		10.0,		//[18]
		18.0,		//[19]
		22.0,		//[20]
		50.0		//[21]
	};

	static const double c[]={
		0.0,		//[0]
		0.0,		//[1]
		0.0,		//[2]
		0.0,		//[3]
		0.0,		//[4]
		0.0,		//[5]
		0.0,		//[6]
		0.0,		//[7]
		0.0,		//[8]
		1.0,		//[9]
		1.0,		//[10]
		1.0,		//[11]
		2.0,		//[12]
		2.0, 		//[13]
		2.0,		//[14]
		2.0,		//[15]
		2.0,		//[16]
		2.0,		//[17]
		3.0,		//[18]
		3.0,		//[19]
		3.0,		//[20]
		4.0			//[21]
	};

	static const double a0[]={
		0.0,		//[0]
		-1.019535,	//[1]
		 9.047135,	//[2]
		-1.629789,	//[3]
		-9.723916,	//[4]
		-3.927170	//[5]
	};
	static const double t0[]={
		0.0,		//[0]
		0.0,		//[1]
		0.0,		//[2]
		0.0,		//[3]
		-1.0/2.0,	//[4]
		-3.0/4.0	//[5]
	};

	phirlist.push_back(new phir_power(n,d,t,c,1,21,22));

	// phi0=log(delta)+a0[1]+a0[2]*tau+a0[3]*log(tau)+a0[4]*pow(tau,-1.0/2.0)+a0[5]*pow(tau,-3.0/4.0);
	phi0list.push_back(new phi0_lead(a0[1],a0[2]));
	phi0list.push_back(new phi0_logtau(a0[3]));
	phi0list.push_back(new phi0_power(a0,t0,4,5,6));

	// Other fluid parameters
	params.molemass = 102.032;
	params.Ttriple = 169.85;
	params.ptriple = 0.3896;
	params.accentricfactor = 0.32684;
	params.R_u = 8.314471;

	// Critical parameters
	crit.rho = 5.017053*params.molemass;
	crit.p = PressureUnit(4059.28, UNIT_KPA);
	crit.T = 374.21;
	crit.v = 1.0/crit.rho;

	// Reducing parameters used in EOS
	reduce.p = PressureUnit(4059.28, UNIT_KPA);
	reduce.T = 374.18;
	reduce.rho = 4.978830171*params.molemass;
	reduce.v = 1.0/reduce.rho;

	preduce = &reduce;

	// Limits of EOS
	limits.Tmin = 169.85;
	limits.Tmax = 455.0;
	limits.pmax = 70000.0;
	limits.rhomax = 15.60*params.molemass;
	
	EOSReference.assign("\"A International Standard Formulation for the Thermodynamic Properties of 1,1,1,2-Tetrafluoroethane" 
					    "(HFC-134a) for Temperatures from 170 K to 455 K and Pressures up to 70 MPa\""
						"by Reiner Tillner-Roth and Hans Dieter Baehr, J. Phys. Chem. Ref. Data, v. 23, 1994, pp 657-729");
	TransportReference.assign("Viscosity: Marcia L. Huber, Arno Laesecke, and Richard A. Perkins, "
							"\"Model for the Viscosity and Thermal Conductivity of Refrigerants,"
							"Including a New Correlation for the Viscosity of R134a\"A Reference "
							"Ind. Eng. Chem. Res. 2003, 42, 3163-3178\n\n"
							"Conductivity: McLinden, M.O., and S.A. Klein and R.A. Perkins, \"An extended corresponding states model for the thermal conductivity of refrigerants and refrigerant mixtures\", International Journal of Refrigeration, 23 (2000) 43-63.\nSeveral typos: Table A1- a1 should be 8.00892e-5, R0 should be 1.03 . Eqn A5 - chi*(Tref,rho) should be multiplied by Tref/T. Eqn A4 - epsilon should be eta (viscosity)\n\n"
							"Surface Tension: Mulero, JPCRD, 2012");

	name.assign("R134a");
	aliases.push_back("R134A");

	BibTeXKeys.EOS = "TillnerRoth-JPCRD-1994";
	BibTeXKeys.VISCOSITY = "Huber-IECR-2003";
	BibTeXKeys.CONDUCTIVITY = "McLinden-IJR-2000";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R134aClass::conductivity_dilute(double T)
{
	// Typo in McLinden, 1999.  a1 should be x10-5
	double a0=-1.05248e-2, //[W/m/K]
		   a1=8.00982e-5;     //[W/m/K^2]
	double lambda = (a0+a1*T); //[W/m/K]
	return lambda; //[W/m/K]
}
double R134aClass::conductivity_residual(double T, double rho)
{
	double b1=1.836526,
		   b2=5.126143,
		   b3=-1.436883,
		   b4=0.626144,
		   lambda_reducing=2.055e-3; //[W/m/K]
	double delta = rho/(5.049886*params.molemass); // Does not use either the reduce.rho value or crit.rho value !!! Or the value listed in Huber, 2003 either
	double lambda_r = lambda_reducing*(b1*delta+b2*pow(delta,2)+b3*pow(delta,3)+b4*pow(delta,4)); //[W/m/K]
	return lambda_r; //[W/m/K]
}
double R134aClass::conductivity_background(double T, double rho) 
{
	return conductivity_residual(T,rho);
}
double R134aClass::conductivity_Trho(double T, double rho)
{
	return conductivity_dilute(T)+conductivity_residual(T,rho)+Fluid::conductivity_critical(T,rho,1.89202e9,0.0496,1.94e-10);
}
void R134aClass::ECSParams(double *e_k, double *sigma)
{
	// Huber, 2003
	*e_k = 299.363;  *sigma = 0.468932;
}
double R134aClass::viscosity_dilute(double T)
{	
	double eta_star, a[]={0.355404,-0.464337,0.257353e-1};
	double e_k = 299.363, sigma = 0.468932, theta_star, Tstar;
	Tstar = T/e_k;
	theta_star = exp(a[0]*pow(log(Tstar),0)+a[1]*pow(log(Tstar),1)+a[2]*pow(log(Tstar),2));
	eta_star = 0.021357*sqrt(params.molemass*T)/(pow(sigma,2)*theta_star)/1e6;
	return eta_star;
}
double R134aClass::viscosity_residual(double T, double rho)
{
	double sum=0,delta_0,DELTA_H_eta,B_eta_star,B_eta,N_A=6.02214129e23,tau,delta;
	double e_k = 299.363 /* K */, sigma = 0.468932/* nm */, Tstar, eta_r;
	double b[]={-19.572881,219.73999,-1015.3226,2471.0125,-3375.1717,2491.6597,-787.26086,14.085455,-0.34664158};
	double t[]={0,-0.25,-0.50,-0.75,-1.00,-1.25,-1.50,-2.50,-5.50};
	double c[]={0,
		-0.206900719e-1, //[1]
		 0.356029549e-3, //[2]
		 0.211101816e-2, //[3]
		 0.139601415e-1, //[4]
		-0.456435020e-2, //[5]
		-0.351593275e-2, //[6]
		 0.214763320,    //[7]
		-0.890173375e-1, //[8]
		 0.100035295,    //[9]
		 3.163695636     //[10]
	};
	tau = T / 374.21;
	/* From Huber (2003):
	The higher-density terms of eq 11, ∆Hη(F,T), were
	formulated in terms of the reduced density δ = F/Fc and
	the reduced temperature τ ) T/Tc with the critical
	parameters of R134a from the equation of state of
	Tillner-Roth and Baehr as reducing parameters, rhoc=
	511.9 kg/m3 */
	delta = rho / 511.9; 
	
	Tstar = T / e_k;
	for (unsigned int i=0;i<=8;i++){
		sum += b[i]*pow(Tstar,t[i]);
	}
	B_eta_star = sum; //[no units]
	B_eta = N_A*pow(sigma/1e9,3)*B_eta_star; //[m3/mol]

	delta_0=c[10]/(1+c[8]*tau+c[9]*pow(tau,2)); //[no units]
	DELTA_H_eta = c[1]*delta + (c[2]/pow(tau,6)+c[3]/pow(tau,2)+c[4]/sqrt(tau)+c[5]*pow(tau,2))*pow(delta,2)
		+c[6]*pow(delta,3)+c[7]/(delta_0-delta)-c[7]/delta_0; //[mPa-s]
	// B_eta*rho needs to be non-dimensional [m3/mol]*[kg/m3] so need divide by mole mass and multiply by 1000
	eta_r = viscosity_dilute(T)*B_eta*rho/params.molemass*1000+DELTA_H_eta/1e3;
	return eta_r;
}
double R134aClass::viscosity_background(double T, double rho)
{
	return viscosity_residual(T,rho);
}
double R134aClass::viscosity_Trho(double T, double rho)
{
	double v = viscosity_dilute(T)+viscosity_residual(T,rho);
	return v;
}
double R134aClass::psat(double T)
{
	const double ti[]={0,1,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.6896400207782598, 2.0859760566425463, -2.6755347075503888, 0.3010493765467589, -5.8583601582059233, 3.4788072104059631};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return crit.p.Pa*exp(summer*crit.T/T);
}
double R134aClass::rhosatL(double T)
{
	const double ti[]={0,0.30476365788328746, 1.3693916347775943, 1.47181218204009};
    const double Ni[]={0,1.4547755306897099, -1.8409578717703114, 1.7741041693195705};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=3;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double R134aClass::rhosatV(double T)
{
	const double ti[]={0,0.36671869805172808, 0.90909437225529799, 3.8975848316804713};
    const double Ni[]={0,-2.5530620653554106, -3.260103853948312, -5.5500764019235618};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=3;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer*crit.T/T);
}
double R134aClass::surface_tension_T(double T)
{
	return 0.05801*pow(1-T/crit.T,1.241);
}


//// STORAGE
///// Not used, but is in theory more accurate than current correlation
///* 
	//From "A Reference Multiparameter Viscosity Equation for R134a
	//with an Optimized Functional Form"
	//by G. Scalabrin and P. Marchi, R. Span
	//J. Phys. Chem. Ref. Data, Vol. 35, No. 2, 2006 
	//*/
	//double sum=0, Tr,rhor;
	//int i;
	//double g[]={0.0,0.0,0.0,1.0,1.0,2.0,2.0,4.0,0.0,2.0,5.0};
	//double h[]={0.0,2.0,20.0,0.0,3.0,0.0,4.0,14.0,1.0,1.0,3.0};
	//double n[]={0.0,0.6564868,0.6882417e-10,0.5668165,-0.2989820,-0.1061795,
	//	0.6245080e-1,0.2758366e-6,-0.1621088,0.1675102,-0.9224693e-1};

	//Tr=T/crit.T;
	//rhor=rho/crit.rho;

	//for (i=1;i<=7;i++)
	//{
	//	sum += n[i]*pow(Tr,g[i])*pow(rhor,h[i]);
	//}
	//for (i=8;i<=10;i++)
	//{
	//	sum += exp(-2*rhor*rhor)*n[i]*pow(Tr,g[i])*pow(rhor,h[i]);
	//}
	//return (exp(sum)-1.0)*25.17975/1e6;

///* 
//	From "A multiparameter thermal conductivity equation
//	for R134a with an optimized functional form"
//	by G. Scalabrin, P. Marchi, F. Finezzo, 
//	Fluid Phase Equilibria 245 (2006) 37-51 
//	*/
//	int i;
//	double sum=0, Tr,rhor,alpha,lambda_r_ce,lambda_r,num,den;
//	double g[]={0.0,0.0,0.5,1.0,1.5,4.0,5.5,6.0,0.0};
//	double h[]={0.0,1.0,1.0,6.0,0.0,3.0,0.0,0.0,1.0};
//	double n[]={0.0,23.504800,-15.261689,0.064403724,7.9735850,0.28675949,
//		8.1819842,-6.4852285,-4.8298888};
//	double nc=1.2478242;
//	double a[]={0.0,1.0,0.0,0.0,0.30,0.30,0.36525,
//		0.61221,0.94930,0.92162,0.15,0.08,0.14};
//
//	Tr=T/crit.T;
//	rhor=rho/crit.rho;
//	alpha=1.0-a[10]*acosh(1+a[11]*pow((1-Tr)*(1-Tr),a[12]));
//	num=rhor*exp(-pow(rhor,a[1])/a[1]-powInt(a[2]*(Tr-1.0),2)-powInt(a[3]*(rhor-1.0),2));
//	den=pow(pow(powInt((1.0-1.0/Tr)+a[4]*pow((rhor-1.0)*(rhor-1.0),1.0/(2.0*a[5])),2),a[6])+pow(a[7]*a[7]*(rhor-alpha)*(rhor-alpha),a[8]),a[9]);
//	lambda_r_ce=num/den;
//	for(i=1;i<=7;i++)
//	{
//		sum+=n[i]*pow(Tr,g[i])*pow(rhor,h[i]);
//	}
//	lambda_r=sum+n[8]*exp(-5.0*rhor*rhor)*pow(Tr,g[8])*pow(rhor,h[8])+nc*lambda_r_ce;
//	return 2.0547*lambda_r/1e6;
