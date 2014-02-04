/* Properties of Carbon Dioxide (R744)
by Ian Bell

Themo properties from 
"A New Equation of State for Carbon Dioxide Covering the Fluid Region from the 
Triple Point Temperature to 1100 K at Pressures up to 800 MPa", 
R. Span and W. Wagner, J. Phys. Chem. Ref. Data, v. 25, 1996

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

#include "math.h"
#include "stdio.h"
#include <string.h>
#include "CoolProp.h"
#include <vector>
#include "FluidClass.h"
#include "R744.h"

R744Class::R744Class()
{

	static double alpha[40],beta[43],GAMMA[40],epsilon[40],a[43],b[43],A[43],B[43],C[43],D[43],a0[9],theta0[9];
	static double n[]={0,    
	 0.388568232032E+00,
	 0.293854759427E+01,
	-0.558671885350E+01,
	-0.767531995925E+00,
	 0.317290055804E+00,
	 0.548033158978E+00,
	 0.122794112203E+00,
	 
	 0.216589615432E+01,
	 0.158417351097E+01,
	-0.231327054055E+00,
	 0.581169164314E-01,
	-0.553691372054E-00,
	 0.489466159094E-00,
	-0.242757398435E-01,
	 0.624947905017E-01,
	-0.121758602252E+00,
	-0.370556852701E+00,
	-0.167758797004E-01,
	-0.119607366380E+00,
	-0.456193625088E-01,
	 0.356127892703E-01, 
	-0.744277271321E-02,
	-0.173957049024E-02,
	-0.218101212895E-01,
	 0.243321665592E-01,
	-0.374401334235E-01,
	 0.143387157569E-00,
	-0.134919690833E-00,
	-0.231512250535E-01,
	 0.123631254929E-01,
	 0.210583219729E-02,
	-0.339585190264E-03,
	 0.559936517716E-02,
	-0.303351180556E-03,

	-0.213654886883E+03,
	 0.266415691493E+05,
	-0.240272122046E+05,
	-0.283416034240E+03,
	 0.212472844002E+03,
	 
	-0.666422765408E+00,
	 0.726086323499E+00,
	 0.550686686128E-01};

	static double d[]={0,
	1,
	1,
	1,
	1,
	2,
	2,
	3,
	1,
	2,
	4,
	5,
	5,
	5,
	6,
	6,
	6,
	1,
	1,
	4,
	4,
	4,
	7,
	8,
	2,
	3,
	3,
	5,
	5,
	6,
	7,
	8,
	10,
	4,
	8,
	2,
	2,
	2,
	3,
	3};

	static double t[]={0.00,
	0.00,
	0.75,
	1.00,
	2.00,
	0.75,
	2.00,
	0.75,
	1.50,
	1.50,
	2.50,
	0.00,
	1.50,
	2.00,
	0.00,
	1.00,
	2.00,
	3.00,
	6.00,
	3.00,
	6.00,
	8.00,
	6.00,
	0.00,
	7.00,
	12.00,
	16.00,
	22.00,
	24.00,
	16.00,
	24.00,
	8.00,
	2.00,
	28.00,
	14.00,
	1.00,
	0.00,
	1.00,
	3.00,
	3.00};

	static double c[]={0,0,0,0,0,0,0,0,
	1,
	1,
	1,
	1,
	1,
	1,
	1,
	1,
	1,
	2,
	2,
	2,
	2,
	2,
	2,
	2,
	3,
	3,
	3,
	4,
	4,
	4,
	4,
	4,
	4,
	5,
	6};

	alpha[35]=25.0;
	alpha[36]=25.0;
	alpha[37]=25.0;
	alpha[38]=15.0;
	alpha[39]=20.0;

	beta[35]=325.0;
	beta[36]=300.0;
	beta[37]=300.0;
	beta[38]=275.0;
	beta[39]=275.0;

	GAMMA[35]=1.16;
	GAMMA[36]=1.19;
	GAMMA[37]=1.19;
	GAMMA[38]=1.25;
	GAMMA[39]=1.22;

	epsilon[35]=1.00;
	epsilon[36]=1.00;
	epsilon[37]=1.00;
	epsilon[38]=1.00;
	epsilon[39]=1.00;

	a[40]=3.5;
	a[41]=3.5;
	a[42]=3.0;

	b[40]=0.875;
	b[41]=0.925;
	b[42]=0.875;

	beta[40]=0.300;
	beta[41]=0.300;
	beta[42]=0.300;

	A[40]=0.700;
	A[41]=0.700;
	A[42]=0.700;

	B[40]=0.3;
	B[41]=0.3;
	B[42]=1.0;

	C[40]=10.0;
	C[41]=10.0;
	C[42]=12.5;

	D[40]=275.0;
	D[41]=275.0;
	D[42]=275.0;

	//Constants for ideal gas expression
	a0[1]=8.37304456;
	a0[2]=-3.70454304;
	a0[3]=2.500000;
	a0[4]=1.99427042;
	a0[5]=0.62105248;
	a0[6]=0.41195293;
	a0[7]=1.04028922;
	a0[8]=0.08327678;

	theta0[4]=3.15163;
	theta0[5]=6.11190;
	theta0[6]=6.77708;
	theta0[7]=11.32384;
	theta0[8]=27.08792;

	std::vector<double> a_v(a,a+sizeof(a)/sizeof(double));
	std::vector<double> b_v(b,b+sizeof(b)/sizeof(double));
	std::vector<double> A_v(A,A+sizeof(A)/sizeof(double));
	std::vector<double> B_v(B,B+sizeof(B)/sizeof(double));
	std::vector<double> C_v(C,C+sizeof(C)/sizeof(double));
	std::vector<double> D_v(D,D+sizeof(D)/sizeof(double));

	phirlist.push_back(new phir_power(n,d,t,c,1,34,35));
	phirlist.push_back(new phir_gaussian(n,d,t,alpha,epsilon,beta,GAMMA,35,39,40));
	phirlist.push_back(new phir_critical(n,d,t,a,b,beta,A,B,C,D,40,42,43));

	std::vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	std::vector<double> theta0_v(theta0,theta0+sizeof(theta0)/sizeof(double));
	
	phi0list.push_back(new phi0_lead(a0[1],a0[2]));
	phi0list.push_back(new phi0_logtau(a0[3]));
	phi0list.push_back(new phi0_Planck_Einstein(a0,theta0,4,8,9));

	// Critical parameters
	crit.rho = 44.0098*10.6249063;
	crit.p = PressureUnit(7377.3, UNIT_KPA);
	crit.T = 304.1282;
	crit.v = 1.0/crit.rho;

	// Other fluid parameters
	params.molemass = 44.0098;
	params.Ttriple = 216.592;
	params.ptriple = 517.996380545;
	params.R_u = 8.31451;
	params.accentricfactor = 0.22394;

	// Limit parameters
	limits.Tmin = 216.592;
	limits.Tmax = 2000.0;
	limits.pmax = 800000.0;
	limits.rhomax = 37.24 * params.molemass;
	
	EOSReference.assign("\"A New Equation of State for Carbon Dioxide Covering the Fluid Region from the "
						"Triple Point Temperature to 1100 K at Pressures up to 800 MPa\", "
						"R. Span and W. Wagner, J. Phys. Chem. Ref. Data, v. 25, 1996");
		
	TransportReference.assign("\"The Transport Properties of Carbon Dioxide\","
							  " V. Vesovic and W.A. Wakeham and G.A. Olchowy and "
							  "J.V. Sengers and J.T.R. Watson and J. Millat"
							  "J. Phys. Chem. Ref. Data, v. 19, 1990");

	name.assign("CarbonDioxide");
	aliases.push_back("R744");
	aliases.push_back("co2");
	aliases.push_back("CO2");
	aliases.push_back("carbondioxide");
	aliases.push_back(std::string("CARBONDIOXIDE"));
	REFPROPname.assign("CO2");

	// Adjust to the IIR reference state (h=200 kJ/kg, s = 1 kJ/kg for sat. liq at 0C)
    params.HSReferenceState = "IIR";

	BibTeXKeys.EOS = "Span-JPCRD-1996";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
	BibTeXKeys.VISCOSITY = "Vesovic-JPCRD-1990";
	BibTeXKeys.CONDUCTIVITY = "Vesovic-JPCRD-1990";
}
double R744Class::conductivity_Trho(double T, double rho)
{
	double e_k=251.196,Tstar;
	//double a[]={0.235156,-0.491266,5.211155e-2,5.347906e-2,-1.537102e-2};
	double b[]={0.4226159,0.6280115,-0.5387661,0.6735941,0,0,-0.4362677,0.2255388};
	double c[]={0,2.387869e-2,4.350794,-10.33404,7.981590,-1.940558};
	double d[]={0,2.447164e-2,8.705605e-5,-6.547950e-8,6.594919e-11};
	//double e[]={0,3.635074e-3,7.209997e-5, 0,0,0,0,3.00306e-20};

	//Vesovic Eq. 31 [no units]
	double summer = 0;
	for (int i=1; i<=5; i++)
		summer += c[i]*pow(T/100.0, 2-i);
	double cint_k = 1.0 + exp(-183.5/T)*summer;

	//Vesovic Eq. 12 [no units]
	double r = sqrt(2.0/5.0*cint_k);

	Tstar=T/e_k;
	//Vesovic Eq. 30 [no units]
	summer = 0;
	for (int i=0; i<=7; i++)
		summer += b[i]/pow(Tstar, i);
	double Gstar_lambda = summer;

	//Vesovic Eq. 29 [mW/m/K]
	double lambda_0 = 0.475598*sqrt(T)*(1+r*r)/Gstar_lambda;

	//Vesovic Eq. 63 [mW/m/K]
	summer = 0;
	for (int i=1; i<=4; i++)
		summer += d[i]*pow(rho, i);
	double delta_lambda = summer;

	// Use the simplified cross-over critical enhancement term
	double delta_c = conductivity_critical(T,rho,1/4e-10,0.052,1.50e-10)*1e6; //[mW/m/K]

	return (lambda_0+delta_lambda+delta_c)/1e3; //[W/m/K]
}
double R744Class::viscosity_Trho(double T, double rho)
{
	int i;
	double e_k=251.196,Tstar,sumGstar=0.0,Gstar,eta0;
	double a[]={0.235156,-0.491266,5.211155e-2,5.347906e-2,-1.537102e-2};
	double e[]={0,3.635074e-3,7.209997e-5, 0,0,0,0,3.00306e-20};

	//double d11=0.4071119e-2,d21=0.7198037e-4,d64=0.2411697e-16,d81=0.2971072e-22,d82=-0.1627888e-22;
	
	Tstar=T/e_k;
	for (i=0;i<=4;i++)
	{
		sumGstar += a[i]*pow(log(Tstar),(int)i);
	}
	Gstar=exp(sumGstar);
	eta0=1.00697*sqrt(T)/Gstar; //[mW/m/K]

	double summer=0;
	for (i=1;i<=7;i++)
		summer += e[i]*pow(rho,(int)i);
	double delta_eta_g = summer; // [mW/m/K]

	// No critical enhancement in viscosity
	//double delta_eta_c = 0.0;

	//double B = 18.56 + 0.014*T;
	//double V0 = 7.41e-4 - 3.3e-7*T;
	//double eta_l = delta_eta_c+1/(B*(1/rho-V0));

	return (eta0+delta_eta_g)/1e6;	
}
double R744Class::psat(double T)
{
	const double ti[]={0,1.0,1.5,2.0,4.0};
    const double ai[]={0,-7.0602087,1.9391218,-1.6463597,-3.2995634};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1-T/crit.T,ti[i]);
    }
    return crit.p.Pa*exp(crit.T/T*summer);
}
double R744Class::rhosatL(double T)
{
	const double ti[]={0,0.340,1.0/2.0,10.0/6.0,11.0/6.0};
    const double ai[]={0,1.9245108,-0.62385555,-0.32731127,0.39245142};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/crit.T,ti[i]);
    }
    return crit.rho*exp(summer);
}
double R744Class::rhosatV(double T)
{
	const double ti[]={0,0.340,1.0/2.0,1.0,7.0/3.0,14.0/3.0};
    const double ai[]={0,-1.7074879,-0.82274670,-4.6008549,-10.111178,-29.742252};
    double summer=0;
    int i;
    for (i=1;i<=5;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/crit.T,ti[i]);
    }
    return crit.rho*exp(summer);
}
