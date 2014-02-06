#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "R236EA.h"

R236EAClass::R236EAClass()
{
	double n[] = {0.0, 0.051074, 2.5584, -2.9180, -0.71485, 0.15534, -1.5894, -0.784, 0.85767, -0.67235, -0.017953, 1.3165, -0.42023, -0.28053, -1.4134, -0.0000062617};
	double t[] = {0, 1.0, 0.264, 0.5638, 1.306, 0.2062, 2.207, 2.283, 1.373, 2.33, 0.6376, 1.08, 1.67, 3.502, 4.357, 0.6945};
	double d[] = {0, 4, 1, 1, 2, 3, 1, 3, 2, 2, 7, 1, 1, 3, 3, 2};
	double c[] = {0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 0, 0, 0};
	double eta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.019, 1.341, 1.034, 5.264, 24.44};
	double beta[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.30, 2.479, 1.068, 79.850, 49.06};
	double gamma[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.13, 0.6691, 0.465, 1.280, 0.8781};
	double epsilon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.7119, 0.9102, 0.678, 0.7091, 1.727};

	//Critical parameters
	crit.rho = 3.716*152.0384; //[kg/m^3]
	crit.p = PressureUnit(3420, UNIT_KPA); //[kPa]
	crit.T = 412.44; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 152.0384;
	params.Ttriple = 170;
	params.accentricfactor = 0.36878238039914035;
	params.R_u = 8.314472;
	params.ptriple = 17.525807103151166; // At Tmin of 243

	// Limits of EOS
	limits.Tmin = 243;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power( n,d,t,c,1,10,16));
	phirlist.push_back(new phir_gaussian( n,d,t, eta, epsilon, beta, gamma, 11,15,16));

	const double n5 = -14.1214241350, n6 = 10.2355589225, n0 = 3.762;
	phi0list.push_back(new phi0_lead(n5,n6));
	phi0list.push_back(new phi0_logtau(n0-1));

	const double u0[] = {0, 144/crit.T, 385/crit.T, 1536/crit.T, 7121/crit.T};
	const double v0[] = {0, 0.7762, 10.41, 12.18, 3.332};
	std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
	std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

	phi0list.push_back(new phi0_Planck_Einstein(v0_v,u0_v,1,4));

	EOSReference.assign("Xinfang Rui and Jiang Pan and Yugang Wang, \"An equation of state for the thermodynamic properties of 1,1,1,2,3,3-hexafluoropropane (R236ea)\", Fluid Phase Equilibria 341 (2013) 78-85\n\nNote: Erratum in paper: a1 should be -17.5983849 and a2 should be 8.87150449\n\nNote: REFPROP 9.0 does not use this EOS, they use thermodynamic corresponding states");
	TransportReference.assign("Using ECS in fully predictive mode.");

	name.assign("R236EA");
	aliases.push_back(std::string("R236ea"));
	REFPROPname.assign("R236EA");

	BibTeXKeys.EOS = "Rui-FPE-2013";
	BibTeXKeys.ECS_FITS = "Huber-IECR-2003";
	BibTeXKeys.ECS_LENNARD_JONES = "Huber-IECR-2003";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}

double R236EAClass::psat(double T)
{
	// Max error is 0.0529123403505 % between 243.0 and 412.439999 K
    const double t[]={0, 0.3525, 0.364, 0.38999999999999996, 1.0, 3.5, 8.5};
    const double N[]={0, -40.377867678370912, 62.473820408057939, -22.503110890281281, -6.5543045908240307, -4.1160029346162421, -19.382010565649018};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double R236EAClass::rhosatL(double T)
{
    // Maximum absolute error is 0.248712 % between 240.000000 K and 412.440000 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.6666666666666667, 1.8333333333333333, 2.1666666666666665, 2.5, 2.8333333333333335, 3.3333333333333335, 3.6666666666666665, 4.166666666666667, 4.666666666666667, 5.333333333333333};
    const double N[] = {0, -0.11265901539856658, 22.489063881386389, -1400.5354055892612, 39725.059739358134, -610848.71399103873, 5635725.3222955344, -33191999.940851696, 128885646.12191778, -330211446.70392752, 532015597.2670126, -446394365.84050828, 285861007.14091557, -276457308.02201492, 210222673.6578109, -161503261.17770496, 116030875.22871566, -40055688.151083663, 10870242.794050815, -1136826.2201711037};
    double summer=0,theta;
    theta=1-T/reduce.T;
    	
for (int i=1; i<=19; i++)
{
    summer += N[i]*pow(theta,t[i]);
}
return reduce.rho*(summer+1);

}

double R236EAClass::rhosatV(double T)
{
    // Maximum absolute error is 0.372083 % between 240.000000 K and 412.440000 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.6666666666666667, 1.8333333333333333, 2.1666666666666665, 2.5, 2.8333333333333335, 3.3333333333333335};
    const double N[] = {0, -0.44429053123604484, 47.74080392708337, -1685.097406350224, 28198.126487234345, -263633.51442363445, 1516482.7144380445, -5690546.1347721191, 14296761.827532202, -23912142.681091785, 25228142.657964848, -13830872.815913403, 3680297.3135379627, -1373672.8249656181, 351025.56021906296, -28436.490307792792};
    double summer=0,theta;
    theta=1-T/reduce.T;
    	
for (int i=1; i<=15; i++)
{
    summer += N[i]*pow(theta,t[i]);
}
return reduce.rho*exp(reduce.T/T*summer);

}
