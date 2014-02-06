/*
Properties for R32.  
by Ian Bell


Thermo props from
Tillner-Roth, R. and Yokozeki, A.,
"An international standard equation of state for difluoromethane (R-32)
for temperatures from the triple point at 136.34 K to 435 K and pressures
up to 70 MPa,"
J. Phys. Chem. Ref. Data, 26(6):1273-1328, 1997.

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
#include <stdlib.h>
#include "CoolProp.h"
#include "FluidClass.h"
#include "R32.h"
#include "R290.h"

R32Class::R32Class()
{
	static const double n[]={
		 0.0,			//[0]
		 0.1046634e+1, 	//[1]
		-0.5451165,		//[2]
		-0.2448595e-2,	//[3]
		-0.4877002e-1,	//[4]
		 0.3520158e-1,	//[5]
		 0.1622750e-2,	//[6]
		 0.2377225e-4,	//[7]
		 0.29149e-1,	//[8]
		 0.3386203e-2,	//[9]
		-0.4202444e-2,	//[10]
		 0.4782025e-3,	//[11]
		-0.5504323e-2,	//[12]
		-0.2418396e-1,	//[13]
		 0.4209034,		//[14]
 		-0.4616537,		//[15]
		-0.1200513e+1,	//[16]
		-0.2591550e+1,	//[17]
		-0.1400145e+1,	//[18]
		 0.8263017		//[19]
	};

	static const double d[]={
		0,			//[0]
		1, 			//[1]
		2, 			//[2]
		5, 			//[3]
		1, 			//[4]
		1, 			//[5]
		3, 			//[6]
		8, 			//[7]
		4, 			//[8]
		4, 			//[9]
		4, 			//[10]
		8, 			//[11]
		3, 			//[12]
		5, 			//[13]
		1, 			//[14]
		1, 			//[15]
		3, 			//[16]
		1, 			//[17]
		2, 			//[18]
		3 			//[19]
	};

	static const double t[]={
		0.0,		//[0]
		1.0/4.0,	//[1]
		1.0,		//[2]
		-1.0/4.0,	//[3]
		-1.0,		//[4]
		2.0,		//[5]
		2.0,		//[6]
		3.0/4.0,	//[7]
		1.0/4.0,	//[8]
		18.0,		//[9]
		26.0,		//[10]
		-1.0,		//[11]
		25.0, 		//[12]
		7.0/4.0,	//[13]
		4.0,		//[14]
		5.0,		//[15]
		1.0,		//[16]
		3.0/2.0,	//[17]
		1.0,		//[18]
		1.0/2.0		//[19]
	};

	static const double c[]={
		0,			//[0]
		0,			//[1]
		0,			//[2]
		0,			//[3]
		0,			//[4]
		0,			//[5]
		0,			//[6]
		0,			//[7]
		0,			//[8]
		4,			//[9]
		3,			//[10]
		1,			//[11]
		4,			//[12]
		1,			//[13]
		2,			//[14]
		2,			//[15]
		1,			//[16]
		1,			//[17]
		1,			//[18]
		1			//[19]
	};

	static const double a0[]={
		-8.258096,	//[0]
		6.353098,	//[1]
		3.004486,	//[2]
		1.160761,	//[3]
		2.645151,	//[4]
		5.794987,	//[5]
		1.129475	//[6]
	};
	static const double n0[]={
		0.0,		//[0]
		0.0,		//[1]
		0.0,		//[2]
		2.2718538,	//[3]
		11.9144210,	//[4]
		5.1415638,	//[5]
		32.7682170	//[6]
	};

	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> n0_v(n0,n0+sizeof(n0)/sizeof(double));

	phirlist.push_back(new phir_power(n,d,t,c,1,19,20));

	// return log(delta)+a0[0]+a0[1]*tau+a0[2]*log(tau)+a0[3]*log(1-exp(-n0[3]*tau))+a0[4]*log(1-exp(-n0[4]*tau))+a0[5]*log(1-exp(-n0[5]*tau))+a0[6]*log(1-exp(-n0[6]*tau));
	phi0list.push_back(new phi0_lead(a0[0],a0[1]));
	phi0list.push_back(new phi0_logtau(a0[2]));
	phi0list.push_back(new phi0_Planck_Einstein(a0_v,n0_v,3,6));

	// Critical parameters
	crit.rho = 8.1500846*52.024;
	crit.p = PressureUnit(5782, UNIT_KPA);
	crit.T = 351.255;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 52.024;
	params.Ttriple = 136.34;
	params.ptriple = 0.0480073825051;
	params.accentricfactor = 0.2769;
	params.R_u = 8.314471;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 435.0;
	limits.pmax = 70000.0;
	limits.rhomax = 27.4734*params.molemass;
	
	EOSReference.assign("Tillner-Roth, R. and Yokozeki, A.,"
						" \"An international standard equation of state for difluoromethane (R-32)"
						" for temperatures from the triple point at 136.34 K to 435 K and pressures"
						" up to 70 MPa,\""
						" J. Phys. Chem. Ref. Data, 26(6):1273-1328, 1997.");
	TransportReference.assign("Surface Tension: R. Heide, \"The surface tension of HFC refrigerants and mixtures\", Int J. Refrig. Vol. 20, No. 7, pp. 496-503, 1997");

	name.assign("R32");

	ECSReferenceFluid = "Propane";

	BibTeXKeys.EOS = "TillnerRoth-JPCRD-1997";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R32Class::psat(double T)
{
	const double ti[]={0,1,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.4655606703362523, 1.804181047910655, -1.5910694825428875, -0.72617761983564866, -3.0701624093185873, 1.4422289706087676};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return crit.p.Pa*exp(summer*crit.T/T);
}

double R32Class::rhosatL(double T)
{
	double theta, phi;
	phi=T/crit.T;
	theta=1-phi;

	return 424.0+434.55*pow(theta,1.0/4.0)+1296.53*pow(theta,2.0/3.0)-777.49*theta+366.84*pow(theta,5.0/3.0);
}
double R32Class::rhosatV(double T)
{
	const double ti[]={0,0.30280531791334198, 0.7420002747530019, 3.8727512808619684};
    const double Ni[]={0,-1.6156590624508722, -3.8438826469427365, -4.429417770439076};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=3;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer*crit.T/T);
}
double R32Class::ECS_psi_viscosity(double rhor)
{
	return 0.7954+5.426580e-2*rhor;
}
double R32Class::ECS_f_int(double T)
{
	return 0.000436654+0.00000178134*T;
}
double R32Class::ECS_chi_conductivity(double rhor)
{
	return 1.2942-0.0924549*rhor;
}
void R32Class::ECSParams(double *e_k, double *sigma)
{
	*e_k = 289.65;
	*sigma = 0.4098;
}
double R32Class::surface_tension_T(double T)
{
	// From Mulero, 2012, JPCRD
	return 0.07147*pow(1-T/reduce.T,1.246);
}
