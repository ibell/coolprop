#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "R125.h"
#include "REFPROP.h"

R125Class::R125Class()
{
	double n[] = {0.0, 5.280760, -8.676580, 0.7501127, 0.7590023, 0.01451899, 4.777189, -3.330988, 3.775673, -2.290919, 0.8888268, -0.6234864, -0.04127263, -0.08455389, -0.1308752, 0.008344962, -1.532005, -0.05883649, 0.02296658};
	double t[] = {0, 0.669, 1.05, 2.75, 0.956, 1.00, 2.00, 2.75, 2.38, 3.37, 3.47, 2.63, 3.45, 0.72, 4.23, 0.20, 4.5, 29.0, 24.0};
	double d[] = {0, 1, 1, 1, 2, 4, 1, 1, 2, 2, 3, 4, 5, 1, 5, 1, 2, 3, 5};
	double l[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 2, 3, 3};
	double m[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.7, 7.0, 6.0};

	//Critical parameters
	crit.rho = 4.779*120.0214; //[kg/m^3]
	crit.p = PressureUnit(3617.7, UNIT_KPA); //[kPa]
	crit.T = 339.173; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 120.0214;
	params.Ttriple = 172.52;
	params.accentricfactor = 0.3052;
	params.R_u = 8.314472;
	params.ptriple = 2.914;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_Lemmon2005( n,d,t,l,m,1,18,19));

	const double a1 = 37.2674, a2 = 8.88404, a3 = -49.8651;
	phi0list.push_back(new phi0_lead(a1,a2));
	phi0list.push_back(new phi0_logtau(-1));
	phi0list.push_back(new phi0_power(a3,-0.1));

	double a[] = {0,0,0,0, 2.303, 5.086, 7.3};
	double b[] = {0,0,0,0, 0.92578, 2.22895, 5.03283};

	std::vector<double> a_v(a,a+sizeof(a)/sizeof(double));
	std::vector<double> b_v(b,b+sizeof(b)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(a_v,b_v,4,6));

	EOSReference.assign("Lemmon-JPCRD-2005");
	TransportReference.assign("Viscosity: Huber IECR 2006\n\n Conductivity: Perkins JCED 2006");
	name.assign("R125");
	REFPROPname.assign("R125");

	BibTeXKeys.EOS = "Lemmon-JPCRD-2005";
	BibTeXKeys.VISCOSITY = "Huber-IECR-2006";
	BibTeXKeys.CONDUCTIVITY = "Perkins-JCED-2006";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2006";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double R125Class::psat(double T)
{
    const double ti[]={0,1.0,1.5,2.3,4.6};
    const double Ni[]={0, -7.5295, 1.9026, -2.2966, -3.4480};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}
double R125Class::rhosatL(double T)
{
    const double ti[]={0,1.0/3.0,0.6,2.9};
    const double Ni[]={0, 1.6684, 0.88415, 0.44383};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=3;i++)
    {
        summer+=Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*(summer+1);
}
double R125Class::rhosatV(double T)
{
    // Maximum absolute error is 0.161887 % between 87.800001 K and 419.289990 K
    const double ti[]={0,0.38,1.22,3.3,6.9};
    const double Ni[]={0, -2.8403, -7.2738, -21.890, -58.825};
    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1; i<=4; i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}
double R125Class::viscosity_Trho(double T, double rho)
{
	double b[] = {-19.572881, 219.73999, -1015.3226, 2471.01251, -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158};
	double N_A = 6.0221415e23;

	double e_k, sigma;
	this->ECSParams(&e_k,&sigma);
	double Tstar = T/e_k;

	double eta_0 = this->viscosity_dilute(T,e_k,sigma)*1e6; // uPa-s

	//Rainwater-Friend for Bstar
	double Bstar = b[0]*pow(Tstar,-0.25*0)+b[1]*pow(Tstar,-0.25*1)+b[2]*pow(Tstar,-0.25*2)+b[3]*pow(Tstar,-0.25*3)+b[4]*pow(Tstar,-0.25*4)+b[5]*pow(Tstar,-0.25*5)+b[6]*pow(Tstar,-0.25*6)+b[7]*pow(Tstar,-2.5)+b[8]*pow(Tstar,-5.5);
	double B = Bstar*N_A*sigma*sigma*sigma/1e27*1000;

	double e[4][4]; // init with zeros
	e[2][1] = 0; e[3][2] = 0;
	e[2][2] = 5.677448e-3; e[3][1] = -5.096662e-3;
	
	double c1 = 1.412564e-1, c2 = 3.033797, c3 = 2.992464e-1;

	double sumresid = 0;
	double tau = T/crit.T, delta = rho/crit.rho;
	for (int i = 2; i<=3; i++)
	{
		for (int j = 1; j <= 2; j++)
		{
			sumresid += e[i][j]*pow(delta,i)/pow(tau,j);
		}
	}
	double delta_0 = c2 + c3*sqrt(T/crit.T);
	double eta_r = (sumresid + c1*(delta/(delta_0-delta)-delta/delta_0))*1000; // uPa-s
	
	double rhobar = rho/params.molemass; // [mol/L]
	return (eta_0*(1+B*rhobar)+eta_r)/1e6;
}
void R125Class::ECSParams(double *e_k, double *sigma)
{
	// Huber IECR 2006
	*e_k = 237.077; *sigma = 0.5235;
}
double R125Class::conductivity_Trho(double T, double rho)
{
	// Perkins JCED 2006

	double sumresid = 0, tau = T/crit.T;
	double lambda_0 = (-4.6082e-3 + 1.68688e-2*tau + 4.88345e-3*tau*tau); //[W/m/K]

	double B1[] = {0, -7.29410e-3, 4.16339e-2, -3.11487e-2, 1.12682e-2, -1.38322e-3};
	double B2[] = {0, 1.10497e-2, -2.89236e-2, 2.78399e-2, -1.2110e-2, 2.11196e-3};
	for (int i = 1; i <= 5; i++)
	{
		sumresid += (B1[i]+B2[i]*T/crit.T)*pow(rho/crit.rho,i);
	}
	double lambda_r = sumresid; //[W/m/K]

	double lambda_c = this->conductivity_critical(T,rho,1.7139e9); //[W/m/K]

	return lambda_0 + lambda_r + lambda_c; //[W/m/K]
}

double R125Class::surface_tension_T(double T)
{
	return 0.05252*pow(1-T/crit.T, 1.237);
}
