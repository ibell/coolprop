/*
Properties for R22.  
by Ian Bell


Thermo props from
A. Kamei, S. W. Beyerlein, and R. T Jacobsen 
"Application of Nonlinear Regression in the 
Development of a Wide Range Formulation 
for HCFC-22"
International Journal of Thermophysics, Vol 16, No. 5, 1995 

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
#include "R22.h"

R22Class::R22Class()
{
	static const double n[]={0,
	0.0695645445236, //[1]
	25.2275419999, //[2]
	-202.351148311, //[3]
	350.063090302, //[4]
	-223.134648863, //[5]
	48.8345904592, //[6]
	0.0108874958556, //[7]
	0.590315073614, //[8]
	-0.689043767432, //[9]
	0.284224445844, //[10]
	0.125436457897, //[11]
	-0.0113338666416, //[12]
	-0.063138895917, //[13]
	0.00974021015232, //[14]
	-0.000408406844722, //[15]
	0.00074194877357, //[16]
	0.000315912525922, //[17]
	0.00000876009723338, //[18]
	-0.000110343340301, //[19]
	-0.0000705323356879, //[20]
	0.23585073151, //[21]
	-0.192640494729, //[22]
	0.00375218008557, //[23]
	-0.0000448926036678, //[24]
	0.0198120520635, //[25]
	-0.0356958425255, //[26]
	0.0319594161562, //[27]
	0.00000260284291078, //[28]
	-0.00897629021967, //[29]
	0.0345482791645, //[30]
	-0.00411831711251, //[31]
	0.00567428536529, //[32]
	-0.00563368989908, //[33]
	0.00191384919423, //[34]
	-0.00178930036389 //[35]
	};

	static const double d[]={0,
	1, //[1]
	1, //[2]
	1, //[3]
	1, //[4]
	1, //[5]
	1, //[6]
	1, //[7]
	2, //[8]
	2, //[9]
	2, //[10]
	3, //[11]
	3, //[12]
	4, //[13]
	5, //[14]
	6, //[15]
	7, //[16]
	7, //[17]
	7, //[18]
	8, //[19]
	8, //[20]
	2, //[21]
	2, //[22]
	2, //[23]
	2, //[24]
	3, //[25]
	4, //[26]
	4, //[27]
	4, //[28]
	4, //[29]
	6, //[30]
	6, //[31]
	6, //[32]
	8, //[33]
	8, //[34]
	8, //[35]
	};

	static const double t[]={0,
		-1, //[1]
		1.75, //[2]
		2.25, //[3]
		2.5, //[4]
		2.75, //[5]
		3, //[6]
		5.5, //[7]
		1.5, //[8]
		1.75, //[9]
		3.5, //[10]
		1, //[11]
		4.5, //[12]
		1.5, //[13]
		0.5, //[14]
		4.5, //[15]
		1, //[16]
		4, //[17]
		5, //[18]
		-0.5, //[19]
		3.5, //[20]
		5, //[21]
		7, //[22]
		12, //[23]
		15, //[24]
		3.5, //[25]
		3.5, //[26]
		8, //[27]
		15, //[28]
		25, //[29]
		3, //[30]
		9, //[31]
		19, //[32]
		2, //[33]
		7, //[34]
		13 //[35]
	};

	static const double l[]={0,
		0, //[1]
		0, //[2]
		0, //[3]
		0, //[4]
		0, //[5]
		0, //[6]
		0, //[7]
		0, //[8]
		0, //[9]
		0, //[10]
		0, //[11]
		0, //[12]
		0, //[13]
		0, //[14]
		0, //[15]
		0, //[16]
		0, //[17]
		0, //[18]
		0, //[19]
		0, //[20]
		2, //[21]
		2, //[22]
		2, //[23]
		2, //[24]
		3, //[25]
		2, //[26]
		2, //[27]
		2, //[28]
		4, //[29]
		2, //[30]
		2, //[31]
		4, //[32]
		2, //[33]
		2, //[34]
		4 //[35]
	};

	static const double B[]={0,
	4352.3095, //[1]
	1935.1591, //[2]
	1887.67936, //[3]
	1694.88284, //[4]
	1605.67848, //[5]
	1162.53424, //[6]
	857.51288, //[7]
	605.72638, //[8]
	530.90982, //[9]
	5.26140446e-3, //[10] 
	1.20662553e-4 //[11]
	};

	std::vector<double> B_v(B,B+sizeof(B)/sizeof(double));
	std::vector<double> b_v(10,1.0);

	phirlist.push_back(new phir_power(n,d,t,l,1,35,36));

	double T0=273.15, 
		   p0=1.0, 
		   R_=8.31451/86.549,
		   rho0=p0/(R_*T0),
		   m,
		   c,
		   R_u=8.31451,
		   Tc=369.295,
		   H0 = 35874.80, /// kJ/kmol
		   S0 = 205.2925, /// kJ/kmol/K
		   tau0=Tc/T0, 
		   delta0=rho0/523.84216696;
	
	// log(delta)+c+m*tau
	
	/// c is the constant term
	c=-S0/R_u-1+log(tau0/delta0);/*<< from the leading term*/

	/// m multiplies the tau term in the leading term (slope)
	m=H0/(R_u*Tc); /*<< from the leading term */

	for (unsigned int i=1; i<=9; i++){B_v[i]/=Tc;}

	phi_BC * phi0_lead_ = new phi0_lead(c,m);
	phi_BC * phi0_logtau_ = new phi0_logtau(-1.0);
	phi_BC * phi0_cp0_constant_ = new phi0_cp0_constant(B[10]+4,Tc,T0);/// checked - good
	phi_BC * phi0_cp0_poly_ = new phi0_cp0_poly(B[11],1,Tc,T0);/// checked - good
	//phi_BC * phi0_cp0_exponential_ = new phi0_cp0_exponential(b_v,B_v,Tc,T0,1,9);/// checked - good

	phi0list.push_back(phi0_lead_);
	phi0list.push_back(phi0_logtau_);
	phi0list.push_back(phi0_cp0_constant_);
	phi0list.push_back(phi0_cp0_poly_);
	//phi0list.push_back(phi0_cp0_exponential_);
	phi0list.push_back(new phi0_Planck_Einstein(b_v,B_v,1,9));
	
	
	// Other fluid parameters
	params.molemass = 86.468;
	params.Ttriple = 115.73;
	params.ptriple = 0.00037947;
	params.accentricfactor = 0.22082;
	params.R_u = 8.314510;

	// Critical parameters
	crit.rho = 6.05822*params.molemass;
	crit.p = PressureUnit(4990, UNIT_KPA);
	crit.T = 369.295;
	crit.v = 1.0/crit.rho;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 550.0;
	limits.pmax = 60000.0;
	limits.rhomax = 19.91*params.molemass;
	
	EOSReference.assign("Kamei, A., and S. W. Beyerlein, and R. T Jacobsen "
						"\"Application of Nonlinear Regression in the "
						"Development of a Wide Range Formulation "
						"for HCFC-22\""
						"International Journal of Thermophysics, Vol 16, No. 5, 1995 ");
	TransportReference.assign("Using ECS");
	name.assign("R22");

	BibTeXKeys.EOS = "Kamei-IJT-1995";
	BibTeXKeys.ECS_LENNARD_JONES = "McLinden-IJR-2000";
	BibTeXKeys.ECS_FITS = "McLinden-IJR-2000";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";

}
double R22Class::rhosatL(double T) 
{
    const double ti[]={0,0.31192915101182356, 2.8946435869795102, 0.67246871922761298};
    const double Ni[]={0,1.483021299160145, 0.080969023427476514, -0.20168767717342609};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=3;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double R22Class::rhosatV(double T) 
{
    const double ti[]={0,0.34946619849137917, 0.85911512967423576, 4.0559379494670686};
    const double Ni[]={0,-2.1778924368532797, -3.1490525869084243, -4.5796870359311495};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=3;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer*crit.T/T);
}
double R22Class::psat(double T) 
{
    const double ti[]={0,1,1.5,2.3,3.6,5.2,7.3};
    const double Ni[]={0,-7.0751284116491497, 1.6540574277960167, -1.6991342890871426, -0.28422811006931187, -3.5966160189328447, 1.217643375964375};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1;i<=6;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return crit.p.Pa*exp(summer*crit.T/T);
}
double R22Class::ECS_f_int(double T)
{
	// McLinden et al., 2000
	return 7.7817e-4+1.2636e-6*T;
}
double R22Class::ECS_chi_conductivity(double rhor)
{
	// McLinden et al, 2000
	return 1.0750-0.0385740*rhor;
}
double R22Class::ECS_psi_viscosity(double rhor)
{
	return 1.0272423-0.0198493*rhor;
}
void R22Class::ECSParams(double *e_k, double *sigma)
{
	*e_k = 284.7242;
	*sigma = 0.4666;
}
