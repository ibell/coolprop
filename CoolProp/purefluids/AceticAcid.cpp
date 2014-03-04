#include "CoolProp.h"
#include <vector>
#include "CPExceptions.h"
#include "FluidClass.h"
#include "AceticAcid.h"
#include "Solvers.h"

/// A stub class to do the density(T,p) calculations for near the critical point using Brent solver
class DensityTpResids : public FuncWrapper1D
{
private:
	double p,T;
	Fluid *pFluid;
public:
	DensityTpResids(Fluid *pFluid, double T, double p){this->pFluid = pFluid; this->p = p; this->T = T;};
	~DensityTpResids(){};
	
	double call(double rho)
	{
		return this->p - pFluid->pressure_Trho(T,rho);
	}
};

AceticAcidClass::AceticAcidClass()
{
	double n[] = {0, -0.15624834164583e1, -0.874703669570960e0, 0.46968858010355e1, 0.97367136204905e-2, -0.49055972708048e-2, 0.24499997808125e2, -0.31443235067567e2, -0.13768156877983e1, 0.14849435860881e1, 0.11374909453775e1, -0.26039791873344e1, -0.30484923493199e-1, 0.53316386834696e1, -0.56733952193640e1, -0.126785566440530e0};
	double d[] = {0, 1, 1, 2, 2, 6, 3, 3, 3, 4, 4, 5, 5, 5, 5, 2};
	double t[] = {0, -1.000, 1.375, 1.000, 1.375, 0.750, -0.250, 0.000, 2.250, 0.125, 2.125, 1.250, 2.250, 2.125, 2.375, 14.000};
	double l[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3};

    // Critical parameters
    crit.rho = 351;
    crit.T = 590.70;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 60.05196;
    params.Ttriple = 289.8;
	params.ptriple = _HUGE;
    params.accentricfactor = 0.48092;
    params.R_u = 8.314472;

    // Limits of EOS
	limits.Tmin = params.Ttriple;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;

	// Residual part
    phirlist.push_back(new phir_power(n,d,t,l,1,15,16));
	double m = 1.01871348, vbarn = 0.444215309e-1,kappabar = 0.109117041e-4 , epsilonbar = 12.2735737;
	phirlist.push_back(new phir_SAFT_associating_1(m, epsilonbar, vbarn, kappabar));

	// Ideal-gas part
	phi0list.push_back(new phi0_lead(-3.94616949, 5.48487930));
	phi0list.push_back(new phi0_logtau(3.66766530));
	phi0list.push_back(new phi0_power(-0.210687796, -1));
	phi0list.push_back(new phi0_power(-0.781330239, -2));
	phi0list.push_back(new phi0_power(0.130979005, -3));
	phi0list.push_back(new phi0_Planck_Einstein( 6.28891793, 2.09502491));

	// Set the critical pressure based on the evaluation of the EOS
	reduce = crit;
	crit.p = PressureUnit(pressure_Trho(crit.T,crit.rho+1e-10), UNIT_PA);

    name.assign("AceticAcid");
    aliases.push_back(std::string("ACETICACID"));
    REFPROPname.assign("N/A");

	BibTeXKeys.EOS = "Piazza-FPE-2011";
	BibTeXKeys.SURFACE_TENSION = "Mulero-JPCRD-2012";

	double Tt = params.Ttriple;
	double Tc = crit.T;
	double pc = crit.p.Pa;

	//double w = 6.67228479e-09*Tc*Tc*Tc-7.20464352e-06*Tc*Tc+3.16947758e-03*Tc-2.88760012e-01;
	//double q = -6.08930221451*w-5.42477887222;
	//double pt = exp(q*(Tc/Tt-1))*pc;

	/*if (1)
	{
		double rhoL,rhoV;

		std::string errstr;
		double T = params.Ttriple;
		double p0 = pt;
		double gibbsL,gibbsV;
		double rhoLstore = _HUGE;
		for (double rho = crit.rho; rho < crit.rho*7; rho += crit.rho/5)
		{
			DensityTpResids DTPR = DensityTpResids(this,T,p0);
			std::string errstr;
			try{
				rhoL = density_Tp(T,p0,rho);
				gibbsL = gibbs_Trho(T,rhoL);
				double pp = pressure_Trho(T,rhoL);
				if (dpdrho_Trho(T,rhoL)>0){
					rhoLstore = rhoL;
				};
				double rr = 0;
			}
			catch (std::exception &)
			{
			}
		}
		double rhoVstore = _HUGE;
		for (double rho = crit.rho; rho > 1e-16; rho /= 1.5)
		{
			try{
				rhoV = density_Tp(T,p0,rho);
				double pp = pressure_Trho(T,rhoV);
				gibbsV = gibbs_Trho(T,rhoV);
				if (dpdrho_Trho(T,rhoV)>0){
					rhoVstore = rhoV;
				};
				double rr = 0;
			}
			catch (std::exception &)
			{
			}
		}
		double p;
		rhoL = rhoLstore; rhoV = rhoVstore;
		FILE *fp;
		fp = fopen("aceticancillary.csv","w");
		fclose(fp);
		while (T < Tc)
		{
			if (T>486)
			{
				rhosatPure_Akasaka(T, &rhoL, &rhoV, &p, 0.1, true);
			}
			else
			{
				rhosatPure(T, &rhoL, &rhoV, &p, 1.0, true);
			}
			fp = fopen("aceticancillary.csv","a+");
			fprintf(fp,"%g,%g,%g,%g\n", T, p, rhoL, rhoV);
			fclose(fp);
			if (T < 580)
			{
				T += 1;
			}
			else
			{
				T += 0.001;
			}
		}
		fclose(fp);
		double rr = 0;
	}*/
}

double AceticAcidClass::psat(double T)
{
    // Max error is  0.104268020681 % between 290.0 and 588.0 K
    const double t[] = {0, 0.065, 0.367, 1.0, 4.0, 4.166666666666667};
    const double N[] = {0, 0.03788266794177342, -0.3014667718204638, -7.465403420229625, 51.68803799661326, -61.422675342107055};
    double summer = 0, theta;
    theta = 1 - T/crit.T;
    for (int i=1; i <= 5; i++)
    {
        summer += N[i]*pow(theta, t[i]);
    }
    return 5780000*exp(crit.T/T*summer);
}

double AceticAcidClass::rhosatL(double T)
{
    // Max error is  0.24407295586 % between 290.0 and 590.0 K
    const double t[] = {0, 0.13, 0.3585, 0.361, 0.365, 0.8333333333333334};
    const double N[] = {0, -8.839463207863323, 73076.09152165553, -118684.81658933076, 45623.75583577903, -3.5208689421322545};
    double summer=0,theta;
    theta = 1 - T/crit.T;
	for (int i=1; i <= 5; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	double val = crit.rho*(1+summer);
    return val;
}

double AceticAcidClass::rhosatV(double T)
{
    // Max error is  1.16918454075 % between 290.0 and 590.0 K
    const double t[] = {0, 0.14600000000000002, 0.38249999999999995, 0.39349999999999996, 0.39749999999999996, 5.333333333333333};
    const double N[] = {0, 7.738193007082631, -6937.911692305578, 25918.893876023878, -18994.78347845314, -14.42897045481017};
    double summer=0,theta;
    theta=1-T/reduce.T;	
	for (int i=1; i<=5; i++)
	{
		summer += N[i]*pow(theta,t[i]);
	}
	double val = crit.rho*exp(crit.T/T*summer);
    return val;
}

#ifndef DISABLE_CATCH
#include "Catch/catch.hpp"
TEST_CASE("Acetic acid validation from Piazza, FPE, 2011","[aceticacid],[validation]")
{
    // Temperature (K) Density (kg/m3) Pressure (MPa) Compress. factor (Z) Enthalpy (kJ/kg) Entropy (kJ/kgK) cv (kJ/kgK) cp (kJ/kgK) Sound speed (m/s)
    double data[10][9] = {
    {290,1e-16,1e-16,1.00000,669.46135,0,0.89579,1.03425,215.30848},
    {290,0.010,0.23586e-3,0.58743,238.11539,1.05462,4.64706,5.35799,158.69275},
    {290,0.025,0.55619e-3,0.55408,203.44932,0.86742,3.46806,3.91294,154.08504},
    {290,1060.,9.5458,0.22429,-218.59529,-0.66978,1.66179,1.99658,1221.38239},
    {290,1070.,22.580,0.52557,-209.68922,-0.68127,1.65651,1.98783,1278.95824},
    {450,1e-16,1e-16,1.00000,868.66312,0,1.31313,1.45159,262.43829},
    {450,1.0,0.55569e-1,0.89189,760.34138,1.94165,3.74492,4.41190,244.91797},
    {450,6.0,0.26514,0.70925,599.36342,1.40868,4.43708,5.53182,212.41678},
    {450,870.,2.7252,0.05028,152.00058,0.35592,2.48256,2.74525,536.47030},
    {450,880.,5.6865,0.10372,153.23341,0.35114,2.46536,2.73491,607.71381}};

    int N = 10;

    SECTION("P")
    {
        for (int i = 0; i < N; i++)
        {
            double T = data[i][0], rho = data[i][1], p = data[i][2]*1e6;
            double PCP = PropsSI("P","T",T,"D",rho,"AceticAcid");
            if (rho < 1e-14){p = 1; PCP = 1;}
            CAPTURE(rho);
            CAPTURE(p);
            CAPTURE(PCP);
            CHECK(fabs(PCP/p-1) < 1e-4);
        }
    }
    SECTION("H")
    {
        for (int i = 0; i < N; i++)
        {
            double T = data[i][0], rho = data[i][1], h = data[i][4]*1e3;
            double HCP = PropsSI("H","T",T,"D",rho,"AceticAcid");
            CAPTURE(h);
            CAPTURE(HCP);
            CHECK(fabs(HCP-h) < 1);
        }
    }
    SECTION("S")
    {
        for (int i = 0; i < N; i++)
        {
            double T = data[i][0], rho = data[i][1], s = data[i][5]*1e3;
            double SCP = PropsSI("S","T",T,"D",rho,"AceticAcid");
            if (rho < 1e-14){s = 1; SCP = 1;}
            CAPTURE(s);
            CAPTURE(SCP);
            CHECK(fabs(SCP-s) < 0.01);
        }
    }
    SECTION("O")
    {
        for (int i = 0; i < N; i++)
        {
            double T = data[i][0], rho = data[i][1], cv = data[i][6]*1e3;
            double OCP = PropsSI("O","T",T,"D",rho,"AceticAcid");
            CAPTURE(cv);
            CAPTURE(OCP);
            CHECK(fabs(OCP-cv) < 1e-2);
        }
    }
    SECTION("C")
    {
        for (int i = 0; i < N; i++)
        {
            double T = data[i][0], rho = data[i][1], cp = data[i][7]*1e3;
            double CCP = PropsSI("C","T",T,"D",rho,"AceticAcid");
            CAPTURE(cp);
            CAPTURE(CCP);
            CHECK(fabs(CCP-cp) < 1e-2);
        }
    }
    SECTION("A")
    {
        for (int i = 0; i < N; i++)
        {
            double T = data[i][0], rho = data[i][1], a = data[i][8];
            double ACP = PropsSI("A","T",T,"D",rho,"AceticAcid");
            CAPTURE(a);
            CAPTURE(ACP);
            CHECK(fabs(ACP-a) < 1e-3);
        }
    }
}
TEST_CASE("Check acetic derivatives residual","[aceticacid],[helmholtz]")
{
    AceticAcidClass Acetic = AceticAcidClass();
	double eps = sqrt(DBL_EPSILON);
	SECTION("dDelta")
	{
        double ANA = Acetic.dphir_dDelta(0.5, 0.5);
		double NUM = (Acetic.phir(0.5, 0.5+eps) - Acetic.phir(0.5,0.5-eps))/(2*eps);
		REQUIRE(fabs(NUM-ANA) < 1e-6);
	}
	SECTION("dTau")
	{
		double ANA = Acetic.dphir_dTau(0.5, 0.5);
		double NUM = (Acetic.phir(0.5+eps, 0.5) - Acetic.phir(0.5-eps,0.5))/(2*eps);
		REQUIRE(fabs(NUM-ANA) < 1e-6);
	}
	SECTION("dDelta2")
	{
		double ANA = Acetic.d2phir_dDelta2(0.5, 0.5);
		double NUM = (Acetic.dphir_dDelta(0.5, 0.5+eps) - Acetic.dphir_dDelta(0.5,0.5-eps))/(2*eps);
		REQUIRE(fabs(NUM-ANA) < 1e-6);
	}
	SECTION("dTau2")
	{
		double ANA = Acetic.d2phir_dTau2(0.5, 0.5);
		double NUM = (Acetic.dphir_dTau(0.5+eps, 0.5) -Acetic.dphir_dTau(0.5-eps,0.5))/(2*eps);
		REQUIRE(fabs(NUM-ANA) < 1e-6);
	}
	SECTION("dDeltadTau")
	{
		double ANA = Acetic.d2phir_dDelta_dTau(0.5, 0.5);
		double NUM = (Acetic.dphir_dTau(0.5, 0.5+eps) - Acetic.dphir_dTau(0.5,0.5-eps))/(2*eps);
		REQUIRE(fabs(NUM-ANA) < 1e-6);
	}
}

TEST_CASE("Check acetic derivatives ideal","[aceticacid],[helmholtz]")
{
    AceticAcidClass Acetic = AceticAcidClass();
	double eps = sqrt(DBL_EPSILON);

	SECTION("dDelta")
	{
        double ANA = Acetic.dphi0_dDelta(0.5, 0.5);
		double NUM = (Acetic.phi0(0.5, 0.5+eps) - Acetic.phi0(0.5,0.5-eps))/(2*eps);
		REQUIRE(fabs(NUM-ANA) < 1e-6);
	}
	SECTION("dTau")
	{
		double ANA = Acetic.dphi0_dTau(0.5, 0.5);
		double NUM = (Acetic.phi0(0.5+eps, 0.5) - Acetic.phi0(0.5-eps,0.5))/(2*eps);
		REQUIRE(fabs(NUM-ANA) < 1e-6);
	}
	SECTION("dDelta2")
	{
		double ANA = Acetic.d2phi0_dDelta2(0.5, 0.5);
		double NUM = (Acetic.dphi0_dDelta(0.5, 0.5+eps) - Acetic.dphi0_dDelta(0.5,0.5-eps))/(2*eps);
		REQUIRE(fabs(NUM-ANA) < 1e-6);
	}
	SECTION("dTau2")
	{
		double ANA = Acetic.d2phi0_dTau2(0.5, 0.5);
		double NUM = (Acetic.dphi0_dTau(0.5+eps, 0.5) -Acetic.dphi0_dTau(0.5-eps,0.5))/(2*eps);
		REQUIRE(fabs(NUM-ANA) < 1e-6);
	}
	SECTION("dDeltadTau")
	{
		double ANA = Acetic.d2phi0_dDelta_dTau(0.5, 0.5);
		double NUM = (Acetic.dphi0_dTau(0.5, 0.5+eps) - Acetic.dphi0_dTau(0.5,0.5-eps))/(2*eps);
		REQUIRE(fabs(NUM-ANA) < 1e-6);
	}
}

#endif