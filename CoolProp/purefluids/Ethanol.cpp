#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "Ethanol.h"

EthanolClass::EthanolClass()
{
    double n[]= {0, 5.8200796E-02, 9.4391227E-01, -8.0941908E-01, 5.5359038E-01, -1.4269032E+00, 1.3448717E-01, 4.2671978E-01, -1.1700261E+00, -9.2405872E-01, 3.4891808E-01, -9.1327720E-01, 2.2629481E-02, -1.5513423E-01, 2.1055146E-01, -2.1997690E-01, -6.5857238E-03, 7.5564749E-01, 1.0694110E-01, -6.9533844E-02, -2.4947395E-01, 2.7177891E-02, -9.0539530E-04, -1.2310953E-01, -8.9779710E-02, -3.9512601E-01};
	double d[] = {0, 4, 1, 1, 2, 2, 3, 1, 1, 1, 3, 3, 2, 2, 6, 6, 8, 1, 1, 2, 3, 3, 2, 2, 2, 1};
	double t[] = {0, 1.00, 1.04, 2.72, 1.174, 1.329, 0.195, 2.43, 1.274, 4.16, 3.30, 4.177, 2.50, 0.81, 2.02, 1.606, 0.86, 2.50, 3.72, 1.19, 3.25, 3.00, 2.00, 2.00, 1.00, 1.00};
	double l[] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 2, 1, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double alpha[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.075, 0.463, 0.876, 1.108, 0.741, 4.032, 2.453, 2.300, 3.143};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.207, 0.0895, 0.581, 0.947, 2.356, 27.01, 4.542, 1.287, 3.090};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.194, 1.986, 1.583, 0.756, 0.495, 1.002, 1.077, 1.493, 1.542};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.779, 0.805, 1.869, 0.694, 1.312, 2.054, 0.441, 0.793, 0.313};

    // Critical parameters
    crit.rho = 5.93*46.06844;
    crit.p = 6268.0;
    crit.T = 514.71;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 46.06844;
    params.Ttriple = 159.1;
	params.ptriple = 7.2e-7;
    params.accentricfactor = 0.644;
    params.R_u = 8.314472;

    // Limits of EOS
	limits.Tmin = params.Ttriple;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

	// Residual part
    phirlist.push_back(new phir_power(n,d,t,l,1,16,26));
	phirlist.push_back(new phir_gaussian(n,d,t,alpha,epsilon,beta,gamma,17,25,26));

	// Ideal-gas part
	const double a[]={0.0, 3.43069, -12.7531, 9.39094, 2.14326, 5.09206, 6.60138, 5.70777};
    const double b[]={0.0, 0, 0, 0, 0.816771, 2.59175, 3.80408, 8.58736};
	std::vector<double> a_v(a, a+sizeof(a)/sizeof(double));
    std::vector<double> b_v(b, b+sizeof(b)/sizeof(double));
	//for (unsigned int i=1;i<v0_v.size();i++) { v0_v[i]/=crit.T; }

	phi0list.push_back(new phi0_lead(a[2], a[3]));
	phi0list.push_back(new phi0_logtau(a[1]));
	phi0list.push_back(new phi0_Planck_Einstein(a_v,b_v,4,7));

    EOSReference.assign("Schroeder Idaho Thesis, 2011");
    TransportReference.assign("Using ECS");

    name.assign("Ethanol");
    aliases.push_back(std::string("C2H6O")); 
    REFPROPname.assign("ETHANOL");

	BibTeXKeys.EOS = "Schroeder-MSTHESIS-2011"; 
	BibTeXKeys.VISCOSITY = "Kiselev-IECR-2005";
	BibTeXKeys.ECS_LENNARD_JONES = "Kiselev-IECR-2005";
	BibTeXKeys.CONDUCTIVITY = "__Kiselev-IECR-2005";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double EthanolClass::psat(double T)
{
    const double ti[]={0, 1.0, 1.5, 3.4, 3.7};
    const double Ni[]={0, -8.94161, 1.61761, -51.1428, 53.1360};

    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=4;i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.p*exp(reduce.T/T*summer);
}
double EthanolClass::rhosatL(double T)
{
    const double ti[] = {0, 0.5, 0.8, 1.1, 1.5, 3.3};
    const double Ni[] = {0, 9.00921, -23.1668, 30.9092, -16.5459, 3.64294};
    double summer=0;
    int i;
    double theta;
    theta=1-T/reduce.T;
    for (i=1;i<=5;i++)
    {
        summer += Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*(summer+1);
}
double EthanolClass::rhosatV(double T)
{
    const double ti[]={0, 0.21, 1.1, 3.4, 10.};
    const double Ni[]={0, -1.75362, -10.5323, -37.6407, -129.762};

    double summer=0,theta;
    int i;
    theta=1.0-T/reduce.T;
    for (i=1; i<=4; i++)
    {
        summer=summer+Ni[i]*pow(theta,ti[i]);
    }
    return reduce.rho*exp(summer);
}

double EthanolClass::viscosity_Trho(double T, double rho)
{
	double eta_0 = -1.03116 + 3.48379e-2*T - 6.50264e-6*T*T;

	// Rainwater-Friend for initial density dependence
	double e_k, sigma;
	this->ECSParams(&e_k,&sigma);
	double Tstar = T/e_k;
	double b[] = {-19.572881, 219.73999, -1015.3226, 2471.01251, -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158};
	double Bstar = b[0]*pow(Tstar,-0.25*0)+b[1]*pow(Tstar,-0.25*1)+b[2]*pow(Tstar,-0.25*2)+b[3]*pow(Tstar,-0.25*3)+b[4]*pow(Tstar,-0.25*4)+b[5]*pow(Tstar,-0.25*5)+b[6]*pow(Tstar,-0.25*6)+b[7]*pow(Tstar,-2.5)+b[8]*pow(Tstar,-5.5);
	double N_A = 6.0221415e23;

	double B = Bstar*N_A*sigma*sigma*sigma/1e27/1000;

	double e[4][3]; // init with zeros
	
	e[2][0] = 0.131194057;
	e[2][1] = -0.382240694; 
	e[2][2] = 0;
	e[3][0] = -0.0805700894;
	e[3][1] = 0.153811778;
	e[3][2] = -0.110578307;

	double c1 = 23.7222995, c2 = -3.38264465, c3 = 12.7568864;

	double sumresid = 0;
	double tau = T/513.9, delta = rho/5.991/46.06844;
	for (int j = 2; j<=3; j++)
	{
		for (int k = 0; k <= 2; k++)
		{
			sumresid += e[j][k]*pow(delta,j)/pow(tau,k);
		}
	}
	double delta_0 = c2 + c3*sqrt(T/crit.T);
	double eta_r = (sumresid + c1*(delta/(delta_0-delta)-delta/delta_0))*1000; // uPa-s

	double rhobar = rho/params.molemass; // [mol/L]
	return (eta_0*(1+B*rhobar)+eta_r)/1e6;
}

double EthanolClass::surface_tension_T(double T)
{
	return 0.05*pow(1-T/reduce.T,0.952);
}