
#include "Mixtures.h"

Mixture::Mixture(std::vector<Fluid *> pFluids)
{
	this->pFluids = pFluids;
	
	std::vector<double> x(2, 0.5);

	STLMatrix beta_v, gamma_v, beta_T, gamma_T;

	/// Resize reducing parameter matrices to be the same size as x in both directions
	beta_v.resize(x.size(),std::vector<double>(x.size(),0));
	gamma_v.resize(x.size(),std::vector<double>(x.size(),0));
	beta_T.resize(x.size(),std::vector<double>(x.size(),0));
	gamma_T.resize(x.size(),std::vector<double>(x.size(),0));
	
	/// For now, only one combination is supported - binary Methane/Ethane
	for (int i = 0; i < x.size(); i++)
	{
		for (int j = 0; j < x.size(); j++)
		{
			if ((!pFluids[i]->get_name().compare("Methane") && !pFluids[j]->get_name().compare("Ethane")) || 
				(!pFluids[j]->get_name().compare("Methane") && !pFluids[i]->get_name().compare("Ethane"))
				)
			{
				beta_v[i][j] = 0.997547866;
				gamma_v[i][j] = 1.006617867;
				beta_T[i][j] = 0.996336508;
				gamma_T[i][j] = 1.049707697;
			}
		}
	}

	STLMatrix F;
	F.resize(x.size(),std::vector<double>(x.size(),1.0));
	/// Methane-Ethane
	GERGReducingFunction pRed = GERGReducingFunction(pFluids, beta_v, gamma_v, beta_T, gamma_T);
	GERGDepartureFunction Excess = GERGDepartureFunction(F);
	IdealMixture IdealMix = IdealMixture(pFluids);

	double Tr = pRed.Tr(x);
	double rhorbar = pRed.rhorbar(x);
	double dTr_dxi = pRed.dTr_dxi(x,0);

	double T = 300;
	double rhobar = 0.5;
	double tau = Tr/T;
	double Rbar = 8.314472;
	double delta = rhobar/rhorbar;

	double dphir_dDelta = Excess.dphir_dDelta(tau,delta,x) + IdealMix.dphir_dDelta(tau,delta,x);
	double p = Rbar*rhobar*T*(1 + delta*dphir_dDelta);

	double rr = 0;
}

double GERGReducingFunction::Tr(std::vector<double> x)
{
	double Tr = 0;
	for (int i = 0; i < x.size(); i++)
	{
		double xi = x[i];
		Tr += xi*xi*pFluids[i]->reduce.T;
	}
	for (int i = 0; i < x.size()-1; i++)
	{
		for (int j = i+1; j < x.size(); j++)
		{
			double xi = x[i], xj = x[j], beta_T_ij = beta_T[i][j];
			Tr += 2*xi*xj*beta_T_ij*gamma_T[i][j]*(xi+xj)/(beta_T_ij*beta_T_ij*xi+xj)*sqrt(pFluids[i]->reduce.T*pFluids[j]->reduce.T);
		}
	}
	return Tr;
}
double GERGReducingFunction::dTr_dxi(std::vector<double> x, int i)
{
	// See Table B9 from Kunz Wagner 2012 (GERG 2008)
	double xi = x[i];
	double dTr_dxi = 2*xi*pFluids[i]->reduce.T;
	for (int k = 0; k < i; k++)
	{
		double xk = x[k], beta_T_ki = beta_T[i][k], gamma_T_ki = gamma_T[i][k];
		double Tr_ki = sqrt(pFluids[i]->reduce.T*pFluids[k]->reduce.T);
		double term = xk*(xk+xi)/(beta_T_ki*beta_T_ki*xk+xi)+xk*xi/(beta_T_ki*beta_T_ki*xk+xi)*(1-(xk+xi)/(beta_T_ki*beta_T_ki*xk+xi));
		dTr_dxi += 2*beta_T_ki*gamma_T_ki*Tr_ki*term;
	}
	for (int k = i+1; k < x.size(); k++)
	{
		double xk = x[k], beta_T_ik = beta_T[i][k], gamma_T_ik = gamma_T[i][k];
		double Tr_ik = sqrt(pFluids[i]->reduce.T*pFluids[k]->reduce.T);
		double term = xk*(xi+xk)/(beta_T_ik*beta_T_ik*xi+xk)+xi*xk/(beta_T_ik*beta_T_ik*xi+xk)*(1-beta_T_ik*beta_T_ik*(xi+xk)/(beta_T_ik*beta_T_ik*xi+xk));
		dTr_dxi += 2*beta_T_ik*gamma_T_ik*Tr_ik*term;
	}
	return dTr_dxi;
}
double GERGReducingFunction::rhorbar(std::vector<double> x)
{
	double vrbar = 0;
	for (int i = 0; i < x.size(); i++)
	{
		double xi = x[i];
		vrbar += xi*xi/pFluids[i]->reduce.rhobar;
	}
	for (int i = 0; i < x.size()-1; i++)
	{
		for (int j = i+1; j < x.size(); j++)
		{
			double xi = x[i], xj = x[j], beta_v_ij = beta_v[i][j];
			vrbar += 2*xi*xj*beta_v_ij*gamma_v[i][j]*(xi+xj)/(beta_v_ij*beta_v_ij*xi+xj)/8.0*pow(pow(pFluids[i]->reduce.rhobar,-1.0/3.0)+pow(pFluids[j]->reduce.rhobar,-1.0/3.0),3.0);
		}
	}
	return 1/vrbar;
}

GERGDepartureFunction::GERGDepartureFunction(STLMatrix F)
{
	// Methane-Ethane
	double d[] = {0, 3, 4, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3};
	double t[] = {0, 0.65, 1.55, 3.1, 5.9, 7.05, 3.35, 1.2, 5.8, 2.7, 0.45, 0.55, 1.95};
	double n[] = {0, -8.0926050298746E-04, -7.5381925080059E-04, -4.1618768891219E-02, -2.3452173681569E-01, 1.4003840584586E-01, 6.3281744807738E-02, -3.4660425848809E-02, -2.3918747334251E-01, 1.9855255066891E-03, 6.1777746171555E+00, -6.9575358271105E+00, 1.0630185306388E+00};
	double eta[] = {0, 0, 0, 1, 1, 1, 0.875, 0.75, 0.5, 0, 0, 0, 0};
	double epsilon[] = {0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
	double beta[] = {0, 0, 0, 1, 1, 1, 1.25, 1.5, 2, 3, 3, 3, 3};
	double gamma[] = {0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

	phi1 = phir_power(n,d,t,1,2,13);
	phi2 = phir_GERG_gaussian(n,d,t,eta,epsilon,beta,gamma,3,12,13);

	this->F = F;
}
double GERGDepartureFunction::phir(double tau, double delta, std::vector<double> x)
{
	double term = phi1.base(tau, delta) + phi2.base(tau, delta);
	double summer = 0;
	for (int i = 0; i < x.size()-1; i++)
	{
		for (int j = i + 1; j < x.size(); j++)
		{	
			summer += x[i]*x[j]*F[i][j]*term;
		}
	}
	return summer;
}
double GERGDepartureFunction::dphir_dDelta(double tau, double delta, std::vector<double> x)
{
	double term = phi1.dDelta(tau, delta) + phi2.dDelta(tau, delta);
	double summer = 0;
	for (int i = 0; i < x.size()-1; i++)
	{
		for (int j = i + 1; j < x.size(); j++)
		{	
			summer += x[i]*x[j]*F[i][j]*term;
		}
	}
	return summer;
}

IdealMixture::IdealMixture(std::vector<Fluid*> pFluids)
{
	this->pFluids = pFluids;
}
double IdealMixture::dphir_dDelta(double tau, double delta, std::vector<double> x)
{
	double summer = 0;
	for (int i = 0; i< x.size(); i++)
	{
		summer += x[i]*pFluids[i]->dphir_dDelta(tau,delta);
	}
	return summer;
}