#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "R12_R113.h"

R12Class::R12Class()
{
	double n[] = {0, 2.075343402E+00, -2.962525996E+00, 1.001589616E-02, 1.781347612E-02, 2.556929157E-02, 2.352142637E-03, -8.495553314E-05, -1.535945599E-02, -2.108816776E-01, -1.654228806E-02, -1.181316130E-02, -4.160295830E-05, 2.784861664E-05, 1.618686433E-06, -1.064614686E-01, 9.369665207E-04, 2.590095447E-02, -4.347025025E-02, 1.012308449E-01, -1.100003438E-01, -3.361012009E-03, 3.789190008E-04};
	double t[] = {0, 0.5, 1, 2, 2.5, -0.5, 0, 0, -0.5, 1.5, 2.5, -0.5, 0, 0.5, -0.5, 4, 4, 2, 4, 12, 14, 0, 14};
	double d[] = {0, 1, 1, 1, 2, 4, 6, 8, 1, 1, 5, 7, 12, 12, 14, 1, 9, 1, 1, 3, 3, 5, 9};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 4};

	//Critical parameters
	crit.rho = 565; //[kg/m^3]
	crit.p = PressureUnit(4136.1, UNIT_KPA); //[kPa]
	crit.T = 385.12; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 120.913;
	params.Ttriple = 116.099;
	params.accentricfactor = 0.17947831734355124;
	params.R_u = 8.314471;
	params.ptriple = 0.000242549651843;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power(n,d,t,c,1,22,23));

	const double a0[] = {0, 0.100100905e2, -0.466434985e1, 0.300361975e1, 0.316062357e1, 0.371258136, 0.356226039e1, 0.212152336e1};
	const double n0[] = {0, 0, 0, 0, 0.372204562e1, 0.630985083e1, 0.178037889e1, 0.107087607e1};
	phi0list.push_back(new phi0_lead(a0[1],a0[2]));
	phi0list.push_back(new phi0_logtau(a0[3]));
	phi0list.push_back(new phi0_Planck_Einstein(a0,n0,4,7,8));

	name.assign("R12");
	REFPROPname.assign("R12");
  
	ECSReferenceFluid = "R134a";

	BibTeXKeys.EOS = "Marx-BOOK-1992";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R12Class::psat(double T)
{
// Max error is  0.183805186373 % between 116.099 and 385.119999 K

    const double t[]={0, 0.061, 0.39599999999999996, 0.39949999999999997, 0.5, 1.1666666666666667, 4.5};
    const double N[]={0, 0.0045343799045095765, -176.10634960801124, 186.36514355948748, -12.261063583187864, -4.3088070968281045, -3.2069974876341432};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1;i<=6;i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double R12Class::rhosatL(double T)
{
    // Maximum absolute error is 0.161217 % between 116.099001 K and 385.119990 K
    const double t[] = {0, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.6666666666666667, 1.8333333333333333};
    const double N[] = {0, -15.83265376429673, 293.31579672214946, -1480.3857720260487, 3652.8805027085969, -4786.5713610565572, 2926.1361104385555, -923.27691484452203, 336.57663013621755};
    double summer=0,theta;
    theta=1-T/reduce.T;    	
	for (int i=1; i<=8; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double R12Class::rhosatV(double T)
{
    // Maximum absolute error is 0.273246 % between 116.099001 K and 385.119990 K
    const double t[] = {0, 0.16666666666666666, 0.3333333333333333, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.5, 1.8333333333333333, 2.1666666666666665, 2.8333333333333335};
    const double N[] = {0, 5.1535942287482994, -181.88689154055641, 2234.5003934456213, -14430.538256207248, 55957.145780671788, -138171.17272878395, 219037.45896919936, -212488.7912512334, 104039.14478204225, -19890.345389677819, 4094.3597114320387, -214.43675518312637};
    double summer=0,theta;
    theta=1-T/reduce.T;    	
	for (int i=1; i<=12; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

R113Class::R113Class()
{
	double n[] = {0, 8.432092286E-01, -2.019185967E+00, 2.920612996E-01, 5.323107661E-02, 3.214971931E-03, 4.667858574E-05, -1.227522799E-06, 8.167288718E-01, -1.340790803E+00, 4.065752705E-01, -1.534754634E-01, -2.414435149E-02, -2.113056197E-02, -3.565436205E-02, 1.364654968E-03, -1.251838755E-02, -1.385761351E-03, 7.206335486E-04};
	double t[] = {0, 0.5, 1.5, 1.5, -0.5, 2, 0, 3, -0.5, 0, 2, 1.5, 6, 2, 10, 6, 18, 15, 33};
	double d[] = {0, 1, 1, 2, 3, 4, 8, 8, 3, 3, 3, 5, 1, 2, 2, 9, 3, 7, 8};
	double c[] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4};

	//Critical parameters
	crit.rho = 560; //[kg/m^3]
	crit.p = PressureUnit(4988.5, UNIT_KPA); //[kPa]
	crit.T = 487.21; //[K]
	crit.v = 1/crit.rho; 

	// Other fluid parameters
	params.molemass = 187.375;
	params.Ttriple = 236.93;
	params.accentricfactor = 0.42002392161549462;
	params.R_u = 8.314471;
	params.ptriple = 1.87;

	// Limits of EOS
	limits.Tmin = params.Ttriple;
	limits.Tmax = 500.0;
	limits.pmax = 100000.0;
	limits.rhomax = 1000000.0*params.molemass;

	phirlist.push_back(new phir_power(n,d,t,c,1,18,19));

	const double a0[] = {0, 0.131479282e2, -0.540537150e1, 0.299999660e1, 0.124464495e2, 0.272181845e1, 0.692712415, 0.332248298e1};
	const double n0[] = {0, 0, 0, 0, 0.104971737e1, 0.329788641e1, 0.862650812e1, 0.329670446e1};
	phi0list.push_back(new phi0_lead(a0[1],a0[2]));
	phi0list.push_back(new phi0_logtau(a0[3]));
	phi0list.push_back(new phi0_Planck_Einstein(a0,n0,4,7,8));

	name.assign("R113");
	REFPROPname.assign("R113");
  
	ECSReferenceFluid = "R134a";

	BibTeXKeys.EOS = "Marx-BOOK-1992";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";
}
double R113Class::psat(double T)
{
	// Max error is  0.498029371079 % between 236.93 and 487.209999 K
    const double t[]={0, 0.07300000000000001, 0.139, 0.38799999999999996, 0.38849999999999996, 0.38999999999999996, 3.0};
    const double N[]={0, -5.327581315717687, 15.505549613757017, -1232825.8272321301, 1643735.9135705857, -410925.79234683927, -4.0193459454122111};
    double summer=0,theta;
    theta=1-T/reduce.T;
    for (int i=1; i<=6; i++)
    {
        summer += N[i]*pow(theta,t[i]);
    }
    return reduce.p.Pa*exp(reduce.T/T*summer);
}

double R113Class::rhosatL(double T)
{
    // Maximum absolute error is 0.043960 % between 236.930001 K and 487.209990 K
    const double t[] = {0, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.6666666666666667, 2.0};
    const double N[] = {0, 178.98095062724852, -1130.7872484737015, 3081.8200478483832, -4230.0445838570422, 2591.7469721955208, -602.56423497728849, 114.08084768423475};
    double summer=0,theta;
    theta=1-T/reduce.T;    	
	for (int i=1; i<=7; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*(summer+1);
}

double R113Class::rhosatV(double T)
{
    // Maximum absolute error is 0.228945 % between 236.930001 K and 487.209990 K
    const double t[] = {0, 0.5, 0.6666666666666666, 0.8333333333333334, 1.0, 1.1666666666666667, 1.3333333333333333, 1.6666666666666667, 2.0};
    const double N[] = {0, 24.762566760651598, -431.27133184173249, 2284.6772084467916, -6030.4478775727466, 8441.9421957649065, -5424.8443870386272, 1453.4553856467098, -328.27673845437351};
    double summer=0,theta;
    theta=1-T/reduce.T;    	
	for (int i=1; i<=8; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	return reduce.rho*exp(reduce.T/T*summer);
}

