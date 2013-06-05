
#include "Mixtures.h"
#include "Solvers.h"

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
	for (unsigned int i = 0; i < x.size(); i++)
	{
		for (unsigned int j = 0; j < x.size(); j++)
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

	//// R32/R134a using form from Lemmon
	//double F_ij = 1.0; // -
	//double zeta_ij = 7.909; // K
	//double xi_ij = -0.002039; // [dm^3/mol] or [L/mol]
	//double N[] = {0,0.22909, 0.094074, 0.00039876, 0.021113};
	//double t[] = {0,1.9, 0.25, 0.07, 2.0};
	//double d[] = {0,1,3,8,1};
	//double l[] = {0,1,1,1,2};

	STLMatrix F;
	F.resize(x.size(),std::vector<double>(x.size(),1.0));
	/// Methane-Ethane
	pReducing = new GERGReducingFunction(pFluids, beta_v, gamma_v, beta_T, gamma_T);
	pExcess = new GERGDepartureFunction(F);
	pResidualIdealMix = new ResidualIdealMixture(pFluids);

	double Tr = pReducing->Tr(x);
	double rhorbar = pReducing->rhorbar(x);
	double dTr_dxi = pReducing->dTr_dxi(x,1);
	double rhorbar_dxi = pReducing->drhorbar_dxi(x,1);

	double T = 200;
	double rhobar = 0.5;
	double tau = Tr/T;
	double Rbar = 8.314472;
	double delta = rhobar/rhorbar;

	double _dphir_dDelta = dphir_dDelta(tau,delta,x);
	double p = Rbar*rhobar*T*(1 + delta*_dphir_dDelta);

	TpzFlash(T,p,x);

	double f0 = fugacity(tau, delta, x, 0);
	double f1 = fugacity(tau, delta, x, 1);

	double rr = 0;
}
Mixture::~Mixture()
{
	if (pReducing != NULL){
		delete pReducing; pReducing = NULL;
	}
	if (pExcess != NULL){
		delete pExcess; pExcess = NULL;
	}
	if (pResidualIdealMix != NULL){
		delete pResidualIdealMix; pResidualIdealMix = NULL;
	}
}
double Mixture::Wilson_lnK_factor(double T, double p, int i)
{
	return log(pFluids[i]->reduce.p/p)+5.373*(1+pFluids[i]->params.accentricfactor)*(1-pFluids[i]->reduce.T/T);
}
double Mixture::fugacity(double tau, double delta, std::vector<double> x, int i)
{
	double Tr = pReducing->Tr(x);
	double rhorbar = pReducing->rhorbar(x);
	double T = Tr/tau, rhobar = rhorbar*delta, Rbar = 8.314472;
	// Fugacity

	double summer_term1 = 0;
	for (unsigned int k = 0; k < x.size(); k++)
	{
		summer_term1 += x[k]*pReducing->drhorbar_dxi(x,k);
	}
	double term1 = delta*dphir_dDelta(tau,delta,x)*(1-1/rhorbar*(pReducing->drhorbar_dxi(x,i)-summer_term1));

	// The second line
	double summer_term2 = 0;
	for (unsigned int k = 0; k < x.size(); k++)
	{
		summer_term2 += x[k]*pReducing->dTr_dxi(x,k);
	}
	double term2 = tau*dphir_dTau(tau,delta,x)*(pReducing->dTr_dxi(x,i)-summer_term2)/Tr;

	// The third line
	double term3 = dphir_dxi(tau,delta,x,i);
	for (unsigned int k = 0; k < x.size(); k++)
	{
		term3 -= x[k]*dphir_dxi(tau,delta,x,k);
	}

	double dnphir_dni = phir(tau,delta,x) + term1 + term2 + term3;

	double f_i = x[i]*rhobar*Rbar*T*exp(dnphir_dni);
	return f_i;
}

double Mixture::phir(double tau, double delta, std::vector<double> x)
{
	return pResidualIdealMix->phir(tau,delta,x) + pExcess->phir(tau,delta,x);
}
double Mixture::dphir_dxi(double tau, double delta, std::vector<double> x, int i)
{	
	return pFluids[i]->phir(tau,delta) + pExcess->dphir_dxi(tau,delta,x,i);
}
double Mixture::dphir_dDelta(double tau, double delta, std::vector<double> x)
{
	return pResidualIdealMix->dphir_dDelta(tau,delta,x) + pExcess->dphir_dDelta(tau,delta,x);
}
double Mixture::dphir_dTau(double tau, double delta, std::vector<double> x)
{
	return pResidualIdealMix->dphir_dTau(tau,delta,x) + pExcess->dphir_dTau(tau,delta,x);
}

/// A wrapper function around the Rachford-Rice residual
class gRR_resid : public FuncWrapper1D
{
public:
	std::vector<double> z,lnK;
	Mixture *Mix;

	gRR_resid(Mixture *Mix, std::vector<double> z, std::vector<double> lnK){ this->z=z; this->lnK = lnK; this->Mix = Mix; };
	double call(double beta){return Mix->g_RachfordRice(z, lnK, beta); };
};

double Mixture::TpzFlash(double T, double p, std::vector<double> z)
{
	int N = z.size();
	double beta;
	std::vector<double> lnK(N), x(N), y(N);

	// Wilson k-factors for each component
	for (unsigned int i = 0; i < N; i++)
	{
		lnK[i] = Wilson_lnK_factor(T, p, i);
	}

	// Check which phase we are in using Wilson estimations
	double g_RR_0 = g_RachfordRice(z, lnK, 0);
	if (g_RR_0 < 0)
	{
		// Subcooled liquid - done
		return _HUGE;
	}
	else
	{
		double g_RR_1 = g_RachfordRice(z, lnK, 1);
		if (g_RR_1 > 0)
		{
			// Superheated vapor - done
			return _HUGE;
		}
	}
	// TODO: How do you know that you aren't in the two-phase region? Safety factor needed?
	
	// Now find the value of beta that satisfies Rachford-Rice

	gRR_resid Resid(this,z,lnK);
	std::string errstr;
	beta = Brent(&Resid,0,1,1e-16,1e-10,300,&errstr);

	// Evaluate mole fractions in liquid and vapor
	for (unsigned int i = 0; i < N; i++)
	{
		double Ki = exp(lnK[i]);
		double den = (1 - beta + beta*Ki); // Common denominator
		// Liquid mole fraction of component i
		x[i] = z[i]/den;
		// Vapor mole fraction of component i
		y[i] = Ki*z[i]/den;
	}

	// Evaluate the fugacities of each component
	double Tr = pReducing->Tr(x);
	double rhorbar = pReducing->rhorbar(x);
	double tau = Tr/T, delta = rhobar/rhorbar;

	return g_RR_0;
}
double Mixture::g_RachfordRice(std::vector<double> z, std::vector<double> lnK, double beta)
{
	// g function from Rashford-Rice
	double summer = 0;
	for (unsigned int i = 0; i < z.size(); i++)
	{
		double Ki = exp(lnK[i]);
		summer += z[i]*(Ki-1)/(1-beta+beta*Ki);
	}
	return summer;
}

double GERGReducingFunction::Tr(std::vector<double> x)
{
	double Tr = 0;
	for (unsigned int i = 0; i < x.size(); i++)
	{
		double xi = x[i];
		Tr += xi*xi*pFluids[i]->reduce.T;
	}
	for (unsigned int i = 0; i < x.size()-1; i++)
	{
		for (unsigned int j = i+1; j < x.size(); j++)
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
	for (unsigned int k = i+1; k < x.size(); k++)
	{
		double xk = x[k], beta_T_ik = beta_T[i][k], gamma_T_ik = gamma_T[i][k];
		double Tr_ik = sqrt(pFluids[i]->reduce.T*pFluids[k]->reduce.T);
		double term = xk*(xi+xk)/(beta_T_ik*beta_T_ik*xi+xk)+xi*xk/(beta_T_ik*beta_T_ik*xi+xk)*(1-beta_T_ik*beta_T_ik*(xi+xk)/(beta_T_ik*beta_T_ik*xi+xk));
		dTr_dxi += 2*beta_T_ik*gamma_T_ik*Tr_ik*term;
	}
	return dTr_dxi;
}
double GERGReducingFunction::drhorbar_dxi(std::vector<double> x, int i)
{
	// See Table B9 from Kunz Wagner 2012 (GERG 2008)
	double xi = x[i];
	double dvrbar_dxi = 2*xi/pFluids[i]->reduce.rhobar;
	for (int k = 0; k < i; k++)
	{
		double xk = x[k], beta_v_ki = beta_v[i][k], gamma_v_ki = gamma_v[i][k];
		double vrbar_ki = 1.0/8.0*pow(pow(pFluids[i]->reduce.rhobar,-1.0/3.0) + pow(pFluids[k]->reduce.rhobar,-1.0/3.0),3.0);
		double term = xk*(xk+xi)/(beta_v_ki*beta_v_ki*xk+xi)+xk*xi/(beta_v_ki*beta_v_ki*xk+xi)*(1-(xk+xi)/(beta_v_ki*beta_v_ki*xk+xi));
		dvrbar_dxi += 2*beta_v_ki*gamma_v_ki*vrbar_ki*term;
	}
	for (unsigned int k = i+1; k < x.size(); k++)
	{
		double xk = x[k], beta_v_ik = beta_v[i][k], gamma_v_ik = gamma_v[i][k];
		double vrbar_ik = 1.0/8.0*pow(pow(pFluids[i]->reduce.rhobar,-1.0/3.0) + pow(pFluids[k]->reduce.rhobar,-1.0/3.0),3.0);
		double term = xk*(xi+xk)/(beta_v_ik*beta_v_ik*xi+xk)+xi*xk/(beta_v_ik*beta_v_ik*xi+xk)*(1-beta_v_ik*beta_v_ik*(xi+xk)/(beta_v_ik*beta_v_ik*xi+xk));
		dvrbar_dxi += 2*beta_v_ik*gamma_v_ik*vrbar_ik*term;
	}
	return -pow(rhorbar(x),2)*dvrbar_dxi;
}

double GERGReducingFunction::rhorbar(std::vector<double> x)
{
	double vrbar = 0;
	for (unsigned int i = 0; i < x.size(); i++)
	{
		double xi = x[i];
		vrbar += xi*xi/pFluids[i]->reduce.rhobar;
	}
	for (unsigned int i = 0; i < x.size()-1; i++)
	{
		for (unsigned int j = i+1; j < x.size(); j++)
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
	for (unsigned int i = 0; i < x.size()-1; i++)
	{
		for (unsigned int j = i + 1; j < x.size(); j++)
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
	for (unsigned int i = 0; i < x.size()-1; i++)
	{
		for (unsigned int j = i + 1; j < x.size(); j++)
		{	
			summer += x[i]*x[j]*F[i][j]*term;
		}
	}
	return summer;
}
double GERGDepartureFunction::dphir_dTau(double tau, double delta, std::vector<double> x)
{
	double term = phi1.dTau(tau, delta) + phi2.dTau(tau, delta);
	double summer = 0;
	for (unsigned int i = 0; i < x.size()-1; i++)
	{
		for (unsigned int j = i + 1; j < x.size(); j++)
		{	
			summer += x[i]*x[j]*F[i][j]*term;
		}
	}
	return summer;
}
double GERGDepartureFunction::dphir_dxi(double tau, double delta, std::vector<double> x, int i)
{
	double summer = 0;
	double term = phi1.base(tau, delta) + phi2.base(tau, delta);
	for (unsigned int k = 0; k < x.size(); k++)
	{
		if (i != k)
		{
			summer += x[k]*F[i][k]*term;
		}
	}
	return summer;
}

ResidualIdealMixture::ResidualIdealMixture(std::vector<Fluid*> pFluids)
{
	this->pFluids = pFluids;
}
double ResidualIdealMixture::phir(double tau, double delta, std::vector<double> x)
{
	double summer = 0;
	for (unsigned int i = 0; i < x.size(); i++)
	{
		summer += x[i]*pFluids[i]->phir(tau,delta);
	}
	return summer;
}
double ResidualIdealMixture::dphir_dDelta(double tau, double delta, std::vector<double> x)
{
	double summer = 0;
	for (unsigned int i = 0; i < x.size(); i++)
	{
		summer += x[i]*pFluids[i]->dphir_dDelta(tau,delta);
	}
	return summer;
}
double ResidualIdealMixture::dphir_dTau(double tau, double delta, std::vector<double> x)
{
	double summer = 0;
	for (unsigned int i = 0; i < x.size(); i++)
	{
		summer += x[i]*pFluids[i]->dphir_dTau(tau,delta);
	}
	return summer;
}