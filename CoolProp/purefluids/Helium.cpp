/* Properties of Helium
by Ian Bell

Thermo properties from 
---------------------
Ortiz-Vega, D.O., Hall, K.R., Arp, V.D., and Lemmon, E.W.,
Interim equation,
to be published in Int. J. Thermophys., 2010.
#Using the EOS constants from REFPROP by permission while awaiting Ortiz-Vega publication in JPCRD

Transport properties from
-------------------------
Arp, V.D., McCarty, R.D., and Friend, D.G.,
"Thermophysical Properties of Helium-4 from 0.8 to 1500 K with
Pressures to 2000 MPa,"
NIST Technical Note 1334 (revised), 1998.

Hands, B.A. and Arp, V.D.,
"A Correlation of Thermal Conductivity Data for Helium,"
Cryogenics, 21(12):697-703, 1981.

*/

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
// The most important line
//#define new new(_NORMAL_BLOCK, __FILE__, __LINE__)
#else
#include <stdlib.h>
#endif

#include "math.h"
#include "stdio.h"
#include <string.h>
#include "CoolProp.h"
#include <vector>
#include <iostream>
#include <list>
#include "Helmholtz.h"
#include "FluidClass.h"
#include "Helium.h"

HeliumClass::HeliumClass()
{
	static const double n[]={0,
	0.009288766, //[1]
	0.9258069, //[2]
	-1.718156, //[3]
	0.7606137, //[4]
	-1.024864, //[5]
	0.1052455, //[6]
	-0.1875722, //[7]
	-0.1287812, //[8]
	-0.002227619, //[9]
	0.1823465, //[10]
	-0.04450014, //[11]
	-0.00008729033, //[12]
	0.0385432, //[13]
	-0.9585106, //[14]
	-0.0545401, //[15]
	-0.0368726, //[16]
	-0.001021851, //[17]
	0.06166348, //[18]
	0.02493437, //[19]
	-0.008127424, //[20]
	-0.008233032, //[21]
	};
	static const double d[]={0,
	4, //[1]
	1, //[2]
	1, //[3]
	2, //[4]
	2, //[5]
	3, //[6]
	1, //[7]
	1, //[8]
	3, //[9]
	2, //[10]
	2, //[11]
	8, //[12]
	1, //[13]
	1, //[14]
	1, //[15]
	2, //[16]
	2, //[17]
	2, //[18]
	3, //[19]
	3, //[20]
	2, //[21]
	};
	static const double t[]={0.00,
	1, //[1]
	0.28, //[2]
	0.735, //[3]
	0.64, //[4]
	0.82, //[5]
	1.16, //[6]
	1.28, //[7]
	2, //[8]
	0.41, //[9]
	1.33, //[10]
	4.2, //[11]
	0.6, //[12]
	3, //[13]
	1, //[14]
	8.2, //[15]
	1, //[16]
	2.71, //[17]
	1, //[18]
	1, //[19]
	2, //[20]
	1, //[21]
	};
	static const double cv[]={0,
	0, //[1]
	0, //[2]
	0, //[3]
	0, //[4]
	0, //[5]
	0, //[6]
	1, //[7]
	2, //[8]
	2, //[9]
	1, //[10]
	2, //[11]
	1, //[12]
	2, //[13]
	2, //[14]
	2, //[15]
	2, //[16]
	2, //[17]
	2, //[18]
	2, //[19]
	2, //[20]
	2, //[21]
	};
	// alpha is used here for consistency with the definitions in R744.c upon which Helium.c is based
	static const double alpha[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-12]
	1.0833, //[13]
	18.3824, //[14]
	5.0573, //[15]
	0.2832, //[16]
	6.0582, //[17]
	0.2444, //[18]
	0.0539, //[19]
	0.185, //[20]
	0.5941, //[21]
	};
	static const double beta[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-12]
	0.0385, //[13]
	19.8246, //[14]
	9.3799, //[15]
	0.8073, //[16]
	0.031, //[17]
	0.0061, //[18]
	0.3581, //[19]
	0.7518, //[20]
	7.4629, //[21]
	};
	static const double GAMMA[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-12]
	1.9776, //[13]
	1.6178, //[14]
	0.4371, //[15]
	0.5355, //[16]
	0.7777, //[17]
	0.4832, //[18]
	0.8162, //[19]
	1.2896, //[20]
	0.3577, //[21]
	};
	static const double epsilon[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0, // indices [0-12]
	0.6914, //[13]
	0.859, //[14]
	0.8787, //[15]
	2.7182, //[16]
	2.0301, //[17]
	0.89, //[18]
	1.179, //[19]
	0.568, //[20]
	1.6412, //[21]
	};

	phirlist.push_back(new phir_power(n,d,t,cv,1,12,13));
	phirlist.push_back(new phir_gaussian(n,d,t,alpha,epsilon,beta,GAMMA,13,21,22));

	// Critical parameters
	crit.rho = 4.002602 * 18.13;
	crit.p = PressureUnit(227.6, UNIT_KPA);
	crit.T = 5.1953;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 4.002602;
	params.Ttriple = 2.1768;
	params.ptriple = 5.05513477113;
	params.accentricfactor = -0.385 ;
	params.R_u = 8.314472;

	double T0 = 4.222,
		   rho0 = 124.95883288439697,
		   m,
		   c,
		   tau0 = crit.T/T0,
		   delta0 = rho0/crit.rho;
	
	// log(delta)+c+m*tau
	
	/// m multiplies the tau term in the leading term (slope)
	/// Constant determined by finding the value of dphi0_dTau that 
	/// yields the enthalpy of 0 at saturated liquid at NBP
	/// Necessary code:
	/// double phi0_dTau = 1/tau*(0.0-delta*DerivTerms("dphir_dDelta",T,rho,"Helium")-1.0)-DerivTerms("dphir_dTau",T,rho,"Helium");
	/// Then dphi0_dTau = m + 1.5/tau, or
	/// m = dphi0_dTau-1.5/tau0
	m = 1.7038767900158605-1.5/tau0;

	/// c is the constant term
	/// phi0 = log(delta) + c + m*tau + 1.5*log(tau)
	/// Necessary code: 
	/// double phi0 = tau*(DerivTerms("dphir_dTau",T,rho,"Helium")+DerivTerms("dphi0_dTau",T,rho,"Helium"))-DerivTerms("phir",T,rho,"Helium")-0;
	/// At reference state you know phi0 from entropy at reference state(0)
	/// phi0 = log(delta0) + c + m*tau0 + 1.5*log(tau0)
	c = 1.6384427034133615 - log(delta0)-m*tau0-1.5*log(tau0);

	phi_BC * phi0_lead_ = new phi0_lead(c,m);
	phi_BC * phi0_logtau_ = new phi0_logtau(1.5);

	phi0list.push_back(phi0_lead_);
	phi0list.push_back(phi0_logtau_);

	// Limits of EOS
	limits.Tmin = 2.1768;
	limits.Tmax = 2000.0;
	limits.pmax = 1000000.0;
	limits.rhomax = 141.22*params.molemass;
	
	EOSReference.assign("Ortiz-Vega, D.O., Hall, K.R., Arp, V.D., and Lemmon, E.W.,"
						"Interim equation to be published in Int. J. Thermophys., 2010."
						"\n\nNote: Using the EOS constants from REFPROP by permission while awaiting Ortiz-Vega publication in JPCRD");
		
	TransportReference.assign("Viscosity: Arp, V.D., McCarty, R.D., and Friend, D.G., "
							  "\"Thermophysical Properties of Helium-4 from 0.8 to 1500 K with Pressures to 2000 MPa\", "
							  "NIST Technical Note 1334 (revised), 1998.\n\n"
							  "Thermal Conductivity: Hands, B.A. and Arp, V.D., "
							  "\"A Correlation of Thermal Conductivity Data for Helium,\" "
							  "Cryogenics, 21(12):697-703, 1981. \n\n"
							  "Warning: Critical enhancement of conductivity not included");

	name.assign("Helium");
	aliases.push_back("helium");
	aliases.push_back("HELIUM");
	aliases.push_back("He");

	BibTeXKeys.EOS = "OrtizVega-2010";
	BibTeXKeys.VISCOSITY = "ARP-NIST-1998";
	BibTeXKeys.CONDUCTIVITY = "Hands-CRYO-1981";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double HeliumClass::viscosity_Trho(double T, double rho)
{
	double eta_0,eta_0_slash, eta_E_slash, B,C,D,ln_eta,x;
	//
	// Arp, V.D., McCarty, R.D., and Friend, D.G.,
	// "Thermophysical Properties of Helium-4 from 0.8 to 1500 K with Pressures to 2000 MPa", 
	// NIST Technical Note 1334 (revised), 1998.
	// 
	// Using Arp NIST report 
	// Report is not clear on viscosity, referring to REFPROP source code for clarity

	// Density in g/cm^3; kg/m^3 --> g/cm^3, divide by 1000
	// Yields viscosity in micro g/(cm-s); to get Pa-s, divide by 10 to get micro Pa-s, then another 1e6 to get Pa-s
	
	rho /= 1000.0;
	if (T <= 300){
		x = log(T);
	}
	else{
		x = log(300.0);
	}
	// Evaluate the terms B,C,D
	B = -47.5295259/x+87.6799309-42.0741589*x+8.33128289*x*x-0.589252385*x*x*x;
	C = 547.309267/x-904.870586+431.404928*x-81.4504854*x*x+5.37008433*x*x*x;
	D = -1684.39324/x+3331.08630-1632.19172*x+308.804413*x*x-20.2936367*x*x*x;
	eta_0_slash = -0.135311743/x+1.00347841+1.20654649*x-0.149564551*x*x+0.012520841*x*x*x;
	eta_E_slash = rho*B+rho*rho*C+rho*rho*rho*D;
	
	if (T<=100)
	{
		ln_eta = eta_0_slash + eta_E_slash;
		return exp(ln_eta)/10.0/1e6;
	}
	else
	{
		ln_eta = eta_0_slash + eta_E_slash;
		eta_0 = 196*pow(T,0.71938)*exp(12.451/T-295.67/T/T-4.1249);
		return (exp(ln_eta)+eta_0-exp(eta_0_slash))/10.0/1e6;
	}
}
double HeliumClass::conductivity_Trho(double T, double rho)
{
	/*
	What an incredibly annoying formulation!  Implied coefficients?? Not cool.
	*/
	double rhoc = 68.0, lambda_e, lambda_c;
	double summer = 3.739232544/T-2.620316969e1/T/T+5.982252246e1/T/T/T-4.926397634e1/T/T/T/T;
	double lambda_0 = 2.7870034e-3*pow(T,7.034007057e-1)*exp(summer);
	double c[]={ 1.862970530e-4,
				-7.275964435e-7,
				-1.427549651e-4,
				 3.290833592e-5,
				-5.213335363e-8,
				 4.492659933e-8,
				-5.924416513e-9,
				 7.087321137e-6,
				-6.013335678e-6,
				 8.067145814e-7,
				 3.995125013e-7};
	// Equation 17
	lambda_e = (c[0]+c[1]*T+c[2]*pow(T,1/3.0)+c[3]*pow(T,2.0/3.0))*rho
			   +(c[4]+c[5]*pow(T,1.0/3.0)+c[6]*pow(T,2.0/3.0))*rho*rho*rho
			   +(c[7]+c[8]*pow(T,1.0/3.0)+c[9]*pow(T,2.0/3.0)+c[10]/T)*rho*rho*log(rho/rhoc);
	
	lambda_c = 0.0;
	return lambda_0+lambda_e+lambda_c;
}
double HeliumClass::psat(double T)
{
	const double ti[]={0,1.0,1.5,1.85,2.7};
    const double ai[]={0,-0.399865e+01,0.870145e+00,0.171451e+00, 0.120927e+01};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1-T/reduce.T,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double HeliumClass::rhosatL(double T)
{
	const double ti[]={0, 1.17, 7.0, 15.0, 20.0};
    const double ai[]={0, 0.140808e+01, -0.543843e+00, 0.177220e+01, -0.344056e+01};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer += ai[i]*pow(pow(1.0-T/reduce.T,1.0/3.0),ti[i]);
    }
    return reduce.rho*(1.0+summer);
}
double HeliumClass::rhosatV(double T)
{
	const double ti[]={0,0.263, 1.04, 3.25, 8.5};
    const double ai[]={0,-0.126074e+01, -0.363425e+01, -0.487998e+01, -0.130581e+02};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/reduce.T,ti[i]);
    }
	return reduce.rho*exp(summer);
}
double HeliumClass::surface_tension_T(double T)
{
	// From Mulero, 2012, JPCRD
	return 0.0004656*pow(1-T/reduce.T,1.04)+0.001889*pow(1-T/reduce.T,2.468)+-0.002006*pow(1-T/reduce.T,2.661);
}
