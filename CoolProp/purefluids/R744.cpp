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

using namespace std;

static double alpha[40],beta[43],GAMMA[40],epsilon[40],a[43],b[43],A[43],B[43],C[43],D[43],a0[9],theta0[9];
static const double Tc=304.128, M_R744=44.0098, rhoc=467.606, Pc=7377.3, _Ttriple=216.592;
             //    K               g/mol       kg/m^3          kPa              K
static const double n[]={0,    
 0.3885682320316100E+00,
 0.2938547594274000E+01,
-0.5586718853493400E+01,
-0.7675319959247700E+00,
 0.3172900558041600E+00,
 0.5480331589776700E+00,
 0.1227941122033500E+00,
 
 0.2165896154322000E+01,
 0.1584173510972400E+01,
-0.2313270540550300E+00,
 0.5811691643143600E-01,
-0.5536913720538200E-00,
 0.4894661590942200E-00,
-0.2427573984350100E-01,
 0.6249479050167800E-01,
-0.1217586022524600E+00,
-0.3705568527008600E+00,
-0.1677587970042600E-01,
-0.1196073663798700E+00,
-0.4561936250877800E-01,
 0.3561278927034600E-01, 
-0.7442772713205200E-02,
-0.1739570490243200E-02,
-0.2181012128952700E-01,
 0.2433216655923600E-01,
-0.3744013342346300E-01,
 0.1433871575687800E-00,
-0.1349196908328600E-00,
-0.2315122505348000E-01,
 0.1236312549290100E-01,
 0.2105832197294000E-02,
-0.3395851902636800E-03,
 0.5599365177159200E-02,
-0.3033511805564600E-03,

-0.2136548868832000E+03,
 0.2664156914927200E+05,
-0.2402721220455700E+05,
-0.2834160342399900E+03,
 0.2124728440017900E+03,
 
-0.6664227654075100E+00,
 0.7260863234989700E+00,
 0.5506866861284200E-01};

static const int d[]={0,
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

static const double t[]={0.00,
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

static const int c[]={0,0,0,0,0,0,0,0,
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

void setCoeffs(void)
{
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
}

R744Class::R744Class()
{
	setCoeffs();
	vector<double> n_v(n,n+sizeof(n)/sizeof(double));
	vector<double> d_v(d,d+sizeof(d)/sizeof(int));
	vector<double> t_v(t,t+sizeof(t)/sizeof(double));
	vector<double> l_v(c,c+sizeof(c)/sizeof(int));
	vector<double> alpha_v(alpha,alpha+sizeof(alpha)/sizeof(double));
	vector<double> beta_v(beta,beta+sizeof(beta)/sizeof(double));
	vector<double> gamma_v(GAMMA,GAMMA+sizeof(GAMMA)/sizeof(double));
	vector<double> epsilon_v(epsilon,epsilon+sizeof(epsilon)/sizeof(double));
	vector<double> a_v(a,a+sizeof(a)/sizeof(double));
	vector<double> b_v(b,b+sizeof(b)/sizeof(double));
	vector<double> A_v(A,A+sizeof(A)/sizeof(double));
	vector<double> B_v(B,B+sizeof(B)/sizeof(double));
	vector<double> C_v(C,C+sizeof(C)/sizeof(double));
	vector<double> D_v(D,D+sizeof(D)/sizeof(double));
	vector<double> a0_v(a0,a0+sizeof(a0)/sizeof(double));
	vector<double> theta0_v(theta0,theta0+sizeof(theta0)/sizeof(double));

	phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,34);
	phi_BC * phirg_ = new phir_gaussian(n_v,d_v,t_v,alpha_v,epsilon_v,beta_v,gamma_v,35,39);
	phi_BC * phirc_ = new phir_critical(n_v,d_v,t_v,a_v,b_v,beta_v,A_v,B_v,C_v,D_v,40,42);

	phirlist.push_back(phir_);
	phirlist.push_back(phirg_);
	phirlist.push_back(phirc_);

	phi_BC * phi0_lead_ = new phi0_lead(a0[1],a0[2]);
	phi_BC * phi0_logtau_ = new phi0_logtau(a0[3]);
	phi_BC * phi0_PE_ = new phi0_Planck_Einstein(a0_v,theta0_v,4,8);
	
	phi0list.push_back(phi0_lead_);
	phi0list.push_back(phi0_logtau_);
	phi0list.push_back(phi0_PE_);

	// Critical parameters
	crit.rho = 467.606;
	crit.p = 7377.3;
	crit.T = 304.128;
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
	REFPROPname.assign("CO2");

}
double R744Class::conductivity_critical(double T, double rho)
{
	double k=1.380658e-23, //[J/K]
		Tref = 1.5*reduce.T, //[K]
		Pcrit = reduce.p, //[kPa]

		//Critical exponents
		nu=0.63,
		gamma=1.2415,

		//Critical amplitudes
		R0=1.01,
		zeta0=1.50e-10, //[m]
		GAMMA = 0.052,
		
		qd = 1/(4.0e-10), //[1/m]
		cp,cv,delta,num,zeta,mu,
		OMEGA_tilde,OMEGA_tilde0,pi=M_PI,tau;

	delta = rho/reduce.rho;

	tau = reduce.T/T;
	double dp_drho=R()*T*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta));
	double X = Pcrit/pow(reduce.rho,2)*rho/dp_drho;
	tau = reduce.T/Tref;
	double dp_drho_ref=R()*Tref*(1+2*delta*dphir_dDelta(tau,delta)+delta*delta*d2phir_dDelta2(tau,delta));
	double Xref = Pcrit/pow(reduce.rho,2)*rho/dp_drho_ref*Tref/T;
	num=X-Xref;

	// no critical enhancement if numerator is negative
	if (num<0)
		return 0.0;
	else
		zeta=zeta0*pow(num/GAMMA,nu/gamma); //[m]

	cp=specific_heat_p_Trho(T,rho); //[kJ/kg/K]
	cv=specific_heat_v_Trho(T,rho); //[kJ/kg/K]
	mu=viscosity_Trho(T,rho)*1e6; //[uPa-s]

	OMEGA_tilde=2.0/pi*((cp-cv)/cp*atan(zeta*qd)+cv/cp*(zeta*qd)); //[-]
	OMEGA_tilde0=2.0/pi*(1.0-exp(-1.0/(1.0/(qd*zeta)+1.0/3.0*(zeta*qd)*(zeta*qd)/delta/delta))); //[-]

	double lambda=rho*cp*1e9*(R0*k*T)/(6*pi*mu*zeta)*(OMEGA_tilde-OMEGA_tilde0); //[W/m/K]
	return lambda*1e3; //[mW/m/K]
}
double R744Class::conductivity_Trho(double T, double rho)
{
	double e_k=251.196,Tstar;
	double a[]={0.235156,-0.491266,5.211155e-2,5.347906e-2,-1.537102e-2};
	double b[]={0.4226159,0.6280115,-0.5387661,0.6735941,0,0,-0.4362677,0.2255388};
	double c[]={0,2.387869e-2,4.350794,-10.33404,7.981590,-1.940558};
	double d[]={0,2.447164e-2,8.705605e-5,-6.547950e-8,6.594919e-11};
	double e[]={0,3.635074e-3,7.209997e-5, 0,0,0,0,3.00306e-20};

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
	double delta_c = conductivity_critical(T,rho); //[mW/m/K]

	return (lambda_0+delta_lambda+delta_c)/1e6; //[kW/m/K]
}
double R744Class::viscosity_Trho(double T, double rho)
{
	int i;
	double e_k=251.196,Tstar,sumGstar=0.0,Gstar,eta0;
	double a[]={0.235156,-0.491266,5.211155e-2,5.347906e-2,-1.537102e-2};
	double e[]={0,3.635074e-3,7.209997e-5, 0,0,0,0,3.00306e-20};

	double d11=0.4071119e-2,d21=0.7198037e-4,d64=0.2411697e-16,d81=0.2971072e-22,d82=-0.1627888e-22;
	double summer=0;
	Tstar=T/e_k;
	for (i=0;i<=4;i++)
	{
		sumGstar += a[i]*powInt(log(Tstar),i);
	}
	Gstar=exp(sumGstar);
	eta0=1.00697*sqrt(T)/Gstar;

	for (i=1;i<=7;i++)
		summer += e[i]*pow(rho,i);
	double delta_eta_g = summer;

	// No critical enhancement in viscosity
	double delta_eta_c = 0.0;

	double B = 18.56 + 0.014*T;
	double V0 = 7.41e-4 - 3.3e-7*T;
	double eta_l = delta_eta_c+1/(B*(1/rho-V0));

	return (eta0+delta_eta_g)/1e6;
}
double R744Class::psat(double T)
{
	const double ti[]={0,1.0,1.5,2.0,4.0};
    const double ai[]={0,-7.0602087,1.9391218,-1.6463597,-3.2995634};
    double summer=0;
    int i;
    setCoeffs();
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1-T/Tc,ti[i]);
    }
    return Pc*exp(Tc/T*summer);
}
double R744Class::rhosatL(double T)
{
	const double ti[]={0,0.340,1.0/2.0,10.0/6.0,11.0/6.0};
    const double ai[]={0,1.9245108,-0.62385555,-0.32731127,0.39245142};
    double summer=0;
    int i;
    for (i=1;i<=4;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(summer);
}
double R744Class::rhosatV(double T)
{
	const double ti[]={0,0.340,1.0/2.0,1.0,7.0/3.0,14.0/3.0};
    const double ai[]={0,-1.7074879,-0.82274670,-4.6008549,-10.111178,-29.742252};
    double summer=0;
    int i;
    for (i=1;i<=5;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(summer);
}
