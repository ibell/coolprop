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
	params.ptriple = 7.3504707747213536e-007;
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
	aliases.push_back(std::string("ethanol")); 
    REFPROPname.assign("ETHANOL");

	BibTeXKeys.EOS = "Schroeder-MSTHESIS-2011"; 
	BibTeXKeys.VISCOSITY = "Kiselev-IECR-2005";
	BibTeXKeys.ECS_LENNARD_JONES = "Kiselev-IECR-2005";
	BibTeXKeys.CONDUCTIVITY = "Assael-JPCRD-2013A";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double EthanolClass::psat(double T)
{
    // Maximum absolute error is 0.049039 % between 159.000001 K and 514.709990 K
    const double t[]={0, 1, 2, 3, 4, 5, 7, 10, 12, 20, 29};
    const double N[]={0, -0.075818069149825001, -7.7448823729682985, -4.5515296488692663, 19.789200573819908, -35.839811025062396, 37.318806908870094, -57.147969051116043, 44.527860526639977, -15.896759939203806, 10.79757895413991};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
    for (i=1;i<=10;i++)
    {
        summer += N[i]*pow(theta,t[i]/2);
    }
    return reduce.p*exp(reduce.T/T*summer);
}

double EthanolClass::rhosatL(double T)
{
    // Maximum absolute error is 0.121983 % between 159.000001 K and 514.709990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.5, 1.8333333333333333, 2.1666666666666665, 2.6666666666666665, 3.0};
    const double N[] = {0, 3.5573356580555995, -115.94068902907544, 1164.8482158801778, -5411.3846003602521, 12919.948271386753, -14718.337097208585, 4955.9389804136181, 6233.4802540917681, -10696.921112569838, 7841.2862300742336, -3088.2381629744345, 915.08823423558681};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
	for (i=1; i<=12; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double EthanolClass::rhosatV(double T)
{
    // Maximum absolute error is 0.249190 % between 159.000001 K and 514.709990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.5, 1.8333333333333333, 2.1666666666666665, 2.6666666666666665, 3.0, 3.5, 4.0, 4.833333333333333, 5.666666666666667, 6.5, 7.5};
    const double N[] = {0, -715.73679606003213, 27174.368520137847, -363760.44006316864, 2497623.4505360974, -9887215.6886598319, 22668282.538672678, -25633604.436509725, 26893871.018254004, -39217522.13402579, 42918084.402004614, -54504033.86915829, 54184715.632529818, -30631923.511888951, 13790730.526646532, -3790817.9701093137, 1376091.7531195069, -369423.82245897013, 42454.025141831691};
    double summer=0,theta;
    int i;
    theta=1-T/reduce.T;
	for (i=1; i<=18; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

double EthanolClass::viscosity_Trho(double T, double rho)
{
	double eta_0 = -1.03116 + 3.48379e-2*T - 6.50264e-6*T*T; //uPa-s

	// Rainwater-Friend for initial density dependence
	double e_k, sigma;
	this->ECSParams(&e_k,&sigma);
	double Tstar = T/e_k;
	double b[] = {-19.572881, 219.73999, -1015.3226, 2471.01251, -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158};
	double Bstar = b[0]*pow(Tstar,-0.25*0)+b[1]*pow(Tstar,-0.25*1)+b[2]*pow(Tstar,-0.25*2)+b[3]*pow(Tstar,-0.25*3)+b[4]*pow(Tstar,-0.25*4)+b[5]*pow(Tstar,-0.25*5)+b[6]*pow(Tstar,-0.25*6)+b[7]*pow(Tstar,-2.5)+b[8]*pow(Tstar,-5.5);
	double B = Bstar*0.60221415*sigma*sigma*sigma;

	double e[4][3]; // init with zeros
	
	e[2][0] = 0.131194057;
	e[2][1] = -0.382240694; 
	e[2][2] = 0;
	e[3][0] = -0.0805700894;
	e[3][1] = 0.153811778;
	e[3][2] = -0.110578307;

	double c1 = 23.7222995, c2 = -3.38264465, c3 = 12.7568864;

	double sumresid = 0;
	double tau = T/513.9, delta = rho/(5.991*46.06844);
	for (int j = 2; j<=3; j++)
	{
		for (int k = 0; k <= 2; k++)
		{
			sumresid += e[j][k]*pow(delta,j)/pow(tau,k);
		}
	}
	double delta_0 = c2 + c3*sqrt(T/513.9);
	double eta_r = (sumresid + c1*(delta/(delta_0-delta)-delta/delta_0))*1000; // uPa-s

	double rhobar = rho/params.molemass; // [mol/L]
	return (eta_0*(1+B*rhobar)+eta_r)/1e6;
}
double EthanolClass::conductivity_Trho(double T, double rho)
{
	double Tr = T/reduce.T;
	double lambda_0 = (-2.09575 + 19.9045*Tr-53.964*Tr*Tr+82.1223*Tr*Tr*Tr-1.98864*Tr*Tr*Tr*Tr-0.495513*Tr*Tr*Tr*Tr*Tr)/(0.17223-0.078273*Tr+Tr*Tr); // [mW/m/K]

	double sumresid = 0;
	double B1[] = {0, 2.67222e-2, 1.48279e-1, -1.30429e-1, 3.46232e-2, -2.44293e-3};
	double B2[] = {0, 1.77166e-2, -8.93088e-2, 6.84664e-2, -1.45702e-2, 8.09189e-4};

	for (int i = 1; i<= 5; i++){		
		sumresid += (B1[i]+B2[i]*(T/reduce.T))*pow(rho/reduce.rho,i);
	}

	double lambda_r = sumresid; // [W/m/K]

	double lambda_c = this->conductivity_critical(T,rho,1.0/(5.3e-10)); // [kW/m/K]

	return (lambda_0+lambda_r*1e3+lambda_c*1e6)/1e6;

}
double EthanolClass::surface_tension_T(double T)
{
	return 0.05*pow(1-T/reduce.T,0.952);
}