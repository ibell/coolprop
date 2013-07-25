
#include "Mixtures.h"
#include "Solvers.h"
#include "CPExceptions.h"

enum PengRobinsonOptions{PR_SATL, PR_SATV};
Mixture::Mixture(std::vector<Fluid *> pFluids)
{
	Rbar = 8.314472;
	this->pFluids = pFluids;
	
	std::vector<double> z(2, 0.5);

	z[0] = 0.5;
	z[1] = 1-z[0];

	STLMatrix beta_v, gamma_v, beta_T, gamma_T;

	/// Resize reducing parameter matrices to be the same size as x in both directions
	beta_v.resize(z.size(),std::vector<double>(z.size(),0));
	gamma_v.resize(z.size(),std::vector<double>(z.size(),0));
	beta_T.resize(z.size(),std::vector<double>(z.size(),0));
	gamma_T.resize(z.size(),std::vector<double>(z.size(),0));
	
	/// For now, only one combination is supported - binary Methane/Ethane
	for (unsigned int i = 0; i < z.size(); i++)
	{
		for (unsigned int j = 0; j < z.size(); j++)
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
	F.resize(z.size(),std::vector<double>(z.size(),1.0));
	/// Methane-Ethane
	pReducing = new GERGReducingFunction(pFluids, beta_v, gamma_v, beta_T, gamma_T);
	pExcess = new GERGDepartureFunction(F);
	pResidualIdealMix = new ResidualIdealMixture(pFluids);

	double Tr = pReducing->Tr(&z);
	double rhorbar = pReducing->rhorbar(&z);
	double dTr_dxi = pReducing->dTr_dxi(&z,1);
	double rhorbar_dxi = pReducing->drhorbar_dxi(&z,1);

	double T = 200;
	double rhobar = 0.5;
	double tau = Tr/T;
	double Rbar = 8.314472;
	double delta = rhobar/rhorbar;

	double _dphir_dDelta = dphir_dDelta(tau,delta,&z);
	double p = Rbar*rhobar*T*(1 + delta*_dphir_dDelta);

	clock_t _t1,_t2;
	long N = 1000000;
	_t1 = clock();
	for (unsigned long long i = 1; i<N; i++)
	{
		z[0] = z[0]+1e-100*i;
		Tr = pReducing->Tr(&z);
	}
	_t2=clock();
	double telap = ((double)(_t2-_t1))/CLOCKS_PER_SEC/(double)N*1e6;
	std::cout << "elap" << telap << std::endl;


	double tau_liq = 0.5;
	double delta_liq  = 0.5;
	double dd1 = (phir(tau_liq+0.0001,delta_liq,&z)- phir(tau_liq-0.0001,delta_liq,&z) )/(2*0.0001);
	double dd2 = dphir_dTau(tau_liq+0.0001,delta_liq,&z);

	//p = 595.61824;
	std::vector<double> x,y;
	//TpzFlash(T, p, z, &rhobar, &x, &y);

	double Tsat;
	for (double x0 = 0; x0 <= 1.000000000000001; x0 += 0.01)
	{
		z[0] = x0; z[1] = 1-x0;
		Tsat = saturation_p(TYPE_BUBBLEPOINT, 1000, &z, &x, &y);
		std::cout << format("%g %g %g %g",x0,Tsat,y[0],y[1]);
		Tsat = saturation_p(TYPE_DEWPOINT, 1000, &z, &x, &y);
		std::cout << format(" %g %g %g",Tsat,x[0],x[1]);
			
		std::cout << std::endl;
	}

	for (double p = 100; p <= 1e9; p *= 1.1)
	{
		double x0 = 0.5;
		z[0] = x0; z[1] = 1-x0;
		Tsat = saturation_p(TYPE_BUBBLEPOINT, p, &z, &x, &y);
		if (!ValidNumber(Tsat)){break;}
		std::cout << format("%g %g %g %g\n",x0,Tsat,y[0],y[1]);
	}

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
double Mixture::fugacity(double tau, double delta, std::vector<double> *x, int i)
{
	double Tr = pReducing->Tr(x);
	double rhorbar = pReducing->rhorbar(x);
	double T = Tr/tau, rhobar = rhorbar*delta, Rbar = 8.314472;

	double dnphir_dni = phir(tau,delta,x) + ndphir_dni(tau,delta,x,i);

	double f_i = (*x)[i]*rhobar*Rbar*T*exp(dnphir_dni);
	return f_i;
}
double Mixture::ndphir_dni(double tau, double delta, std::vector<double> *x, int i)
{
	double Tr = pReducing->Tr(x);
	double rhorbar = pReducing->rhorbar(x);
	double T = Tr/tau, rhobar = rhorbar*delta, Rbar = 8.314472;

	double summer_term1 = 0;
	for (unsigned int k = 0; k < (*x).size(); k++)
	{
		summer_term1 += (*x)[k]*pReducing->drhorbar_dxi(x,k);
	}
	double term1 = delta*dphir_dDelta(tau,delta,x)*(1-1/rhorbar*(pReducing->drhorbar_dxi(x,i)-summer_term1));

	// The second line
	double summer_term2 = 0;
	for (unsigned int k = 0; k < (*x).size(); k++)
	{
		summer_term2 += (*x)[k]*pReducing->dTr_dxi(x,k);
	}
	double term2 = tau*dphir_dTau(tau,delta,x)*(pReducing->dTr_dxi(x,i)-summer_term2)/Tr;

	// The third line
	double term3 = dphir_dxi(tau,delta,x,i);
	for (unsigned int k = 0; k < (*x).size(); k++)
	{
		term3 -= (*x)[k]*dphir_dxi(tau,delta,x,k);
	}
	return term1 + term2 + term3;
}

double Mixture::dndphir_dni_dTau(double tau, double delta, std::vector<double> *x, int i)
{
	double Tr = pReducing->Tr(x);
	double rhorbar = pReducing->rhorbar(x);
	double T = Tr/tau, rhobar = rhorbar*delta, Rbar = 8.314472;

	double summer_term1 = 0;
	for (unsigned int k = 0; k < (*x).size(); k++)
	{
		summer_term1 += (*x)[k]*pReducing->drhorbar_dxi(x,k);
	}
	double term1 = delta*d2phir_dDelta_dTau(tau,delta,x)*(1-1/rhorbar*(pReducing->drhorbar_dxi(x,i)-summer_term1));

	// The second line
	double summer_term2 = 0;
	for (unsigned int k = 0; k < (*x).size(); k++)
	{
		summer_term2 += (*x)[k]*pReducing->dTr_dxi(x,k);
	}
	double term2 = (tau*d2phir_dTau2(tau,delta,x)+dphir_dTau(tau,delta,x))*(pReducing->dTr_dxi(x,i)-summer_term2)/Tr;

	// The third line
	double term3 = d2phir_dxi_dTau(tau,delta,x,i);
	for (unsigned int k = 0; k < (*x).size(); k++)
	{
		term3 -= (*x)[k]*d2phir_dxi_dTau(tau,delta,x,k);
	}
	return term1 + term2 + term3;
}

double Mixture::phir(double tau, double delta, std::vector<double> *x)
{
	return pResidualIdealMix->phir(tau,delta,x) + pExcess->phir(tau,delta,x);
}
double Mixture::dphir_dxi(double tau, double delta, std::vector<double> *x, int i)
{	
	return pFluids[i]->phir(tau,delta) + pExcess->dphir_dxi(tau,delta,x,i);
}
double Mixture::d2phir_dxi_dTau(double tau, double delta, std::vector<double> *x, int i)
{	
	return pFluids[i]->dphir_dTau(tau,delta) + pExcess->d2phir_dxi_dTau(tau,delta,x,i);
}
double Mixture::dphir_dDelta(double tau, double delta, std::vector<double> *x)
{
	return pResidualIdealMix->dphir_dDelta(tau,delta,x) + pExcess->dphir_dDelta(tau,delta,x);
}
double Mixture::d2phir_dDelta_dTau(double tau, double delta, std::vector<double> *x)
{
	return pResidualIdealMix->d2phir_dDelta_dTau(tau,delta,x) + pExcess->d2phir_dDelta_dTau(tau,delta,x);
}
double Mixture::d2phir_dTau2(double tau, double delta, std::vector<double> *x)
{
	return pResidualIdealMix->d2phir_dTau2(tau,delta,x) + pExcess->d2phir_dTau2(tau,delta,x);
}
double Mixture::dphir_dTau(double tau, double delta, std::vector<double> *x)
{
	return pResidualIdealMix->dphir_dTau(tau,delta,x) + pExcess->dphir_dTau(tau,delta,x);
}

/// A wrapper function around the Rachford-Rice residual
class gRR_resid : public FuncWrapper1D
{
public:
	std::vector<double> *z,*lnK;
	Mixture *Mix;

	gRR_resid(Mixture *Mix, std::vector<double> *z, std::vector<double> *lnK){ this->z=z; this->lnK = lnK; this->Mix = Mix; };
	double call(double beta){return Mix->g_RachfordRice(z, lnK, beta); };
};

/// A wrapper function around the density(T,p,x) residual
class rho_Tpz_resid : public FuncWrapper1D
{
protected:
	double T, p, Rbar, tau, Tr, rhorbar;
public:
	std::vector<double> *x;
	Mixture *Mix;

	rho_Tpz_resid(Mixture *Mix, double T, double p, std::vector<double> *x){ 
		this->x=x; this->T = T; this->p = p; this->Mix = Mix;
		Tr = Mix->pReducing->Tr(x);
		rhorbar = Mix->pReducing->rhorbar(x);
		tau = Tr/T;
		Rbar = 8.314472; // kJ/kmol/K
	};
	double call(double rhobar){	
		double delta = rhobar/rhorbar;
		double resid = Rbar*rhobar*T*(1 + delta*Mix->dphir_dDelta(tau, delta, x))-p;
		return resid;
	}
};
double Mixture::rhobar_Tpz(double T, double p, std::vector<double> *x, double rhobar0)
{
	rho_Tpz_resid Resid(this,T,p,x);
	std::string errstr;
	return Secant(&Resid, rhobar0, 0.00001, 1e-8, 100, &errstr);
}



//double Mixture::saturation_p_NewtonRaphson(int type, double T, double p, std::vector<double> *z, std::vector<double> *ln_phi_liq, std::vector<double> *ln_phi_vap, std::vector<double> *x, std::vector<double> *y)
//{
//
//	
//}


/*! A wrapper function around the residual to find the initial guess for the bubble point temperature
\f[
r = \sum_i \left[z_i K_i\right] - 1 
\f]
*/
class bubblepoint_WilsonK_resid : public FuncWrapper1D
{
public:
	double p;
	std::vector<double> *z;
	Mixture *Mix;

	bubblepoint_WilsonK_resid(Mixture *Mix, double p, std::vector<double> *z){ 
		this->z=z; this->p = p; this->Mix = Mix; 
	};
	double call(double T){
		double summer = 0;
		for (unsigned int i = 0; i< (*z).size(); i++) { summer += (*z)[i]*exp(Mix->Wilson_lnK_factor(T,p,i)); }
		return summer - 1; // 1 comes from the sum of the z_i which must sum to 1
	};
};
/*! A wrapper function around the residual to find the initial guess for the dew point temperature
\f[
r = 1- \sum_i \left[z_i/K_i\right]
\f]
*/
class dewpoint_WilsonK_resid : public FuncWrapper1D
{
public:
	double p;
	std::vector<double> *z;
	Mixture *Mix;

	dewpoint_WilsonK_resid(Mixture *Mix, double p, std::vector<double> *z){ 
		this->z=z; this->p = p; this->Mix = Mix; 
	};
	double call(double T){
		double summer = 0;
		for (unsigned int i = 0; i< (*z).size(); i++) { summer += (*z)[i]*(1-1/exp(Mix->Wilson_lnK_factor(T,p,i))); }
		return summer;
	};
};

double Mixture::saturation_p(int type, double p, std::vector<double> *z, std::vector<double> *x, std::vector<double> *y)
{
	int iter = 0;
	double change, T, f, dfdT;
	unsigned int N = (*z).size();
	std::vector<double> K(N), ln_phi_liq(N), ln_phi_vap(N);

	(*x).resize(N);
	(*y).resize(N);

	if (type == TYPE_BUBBLEPOINT)
	{
		// Liquid is at the bulk composition
		*x = (*z);
		// Find first guess for T using Wilson K-factors
		bubblepoint_WilsonK_resid Resid(this,p,z); //sum(z_i*K_i) - 1
		std::string errstr;
		double Tr = pReducing->Tr(z);
		// Try a range of different values for the temperature, hopefully one works
		for (double T_guess = Tr*0.9; T_guess > 0; T_guess -= Tr*0.1)
		{
			try{
				T = Secant(&Resid, T_guess, 0.001, 1e-10, 100, &errstr);
				if (!ValidNumber(T)){throw ValueError();}
				break;
			} catch (CoolPropBaseError) {}
		}
	}
	else if (type == TYPE_DEWPOINT)
	{
		// Vapor is at the bulk composition
		*y = (*z);
		// Find first guess for T using Wilson K-factors
		dewpoint_WilsonK_resid Resid(this,p,z); //1-sum(z_i/K_i)
		std::string errstr;
		double Tr = pReducing->Tr(z);
		// Try a range of different values for the temperature, hopefully one works
		for (double T_guess = Tr*0.9; T_guess > 0; T_guess -= Tr*0.1)
		{
			try{
				T = Secant(&Resid, T_guess, 0.001, 1e-10, 100, &errstr);
				if (!ValidNumber(T)){throw ValueError();}
				break;
			} catch (CoolPropBaseError) {}
		}
	}
	else
	{
		throw ValueError("Invalid type to saturation_p");
	}

	// Calculate the K factors for each component
	for (unsigned int i = 0; i < N; i++)
	{
		K[i] = exp(Wilson_lnK_factor(T,p,i));
	}	

	// Initial guess for mole fractions in the other phase
	if (type == TYPE_BUBBLEPOINT)
	{
		// Calculate the vapor molar fractions using the K factor and Rachford-Rice
		for (unsigned int i = 0; i < N; i++)
		{
			(*y)[i] = K[i]*(*z)[i];
		}
	}
	else
	{
		// Calculate the liquid molar fractions using the K factor and Rachford-Rice
		for (unsigned int i = 0; i < N; i++)
		{
			(*x)[i] = (*z)[i]/K[i];
		}
	}

	do
	{
		double rhobar_liq = rhobar_Tpz(T, p, x, rhobar_pengrobinson(T,p,x,PR_SATL));
		double rhobar_vap = rhobar_Tpz(T, p, y, rhobar_pengrobinson(T,p,y,PR_SATV));
		double Tr_liq = pReducing->Tr(x);
		double Tr_vap = pReducing->Tr(y); 
		double tau_liq = Tr_liq/T;
		double tau_vap = Tr_vap/T;
		double rhorbar_liq = pReducing->rhorbar(x);
		double rhorbar_vap = pReducing->rhorbar(y);
		double delta_liq = rhobar_liq/rhorbar_liq; 
		double delta_vap = rhobar_vap/rhorbar_vap; 
		double dtau_dT_liq = -Tr_liq/T/T;
		double dtau_dT_vap = -Tr_vap/T/T;

		double Z_liq = p/(Rbar*rhobar_liq*T);
		double Z_vap = p/(Rbar*rhobar_vap*T);

		double dZ_liq_dT = -p/(Rbar*rhobar_liq*T*T);
		double dZ_vap_dT = -p/(Rbar*rhobar_vap*T*T);

		double phir_liq = phir(tau_liq, delta_liq, x);
		double phir_vap = phir(tau_vap, delta_vap, y);

		double dphir_liq_dT = dphir_dTau(tau_liq, delta_liq, x)*dtau_dT_liq;
		double dphir_vap_dT = dphir_dTau(tau_vap, delta_vap, y)*dtau_dT_vap;

		f = 0;
		dfdT = 0;
		for (unsigned int i=0; i < N; i++)
		{
			ln_phi_liq[i] = phir_liq + ndphir_dni(tau_liq,delta_liq,x,i)-log(Z_liq);
			ln_phi_vap[i] = phir_vap + ndphir_dni(tau_vap,delta_vap,y,i)-log(Z_vap);

			double dln_phi_liq_dT = dphir_liq_dT + dndphir_dni_dTau(tau_liq, delta_liq, x, i)*dtau_dT_liq-1/Z_liq*dZ_liq_dT;
			double dln_phi_vap_dT = dphir_vap_dT + dndphir_dni_dTau(tau_vap, delta_vap, y, i)*dtau_dT_vap-1/Z_vap*dZ_vap_dT;

			double Ki = exp(ln_phi_liq[i] - ln_phi_vap[i]);
			
			if (type == TYPE_BUBBLEPOINT){
				f += (*z)[i]*(Ki-1);
				dfdT += (*z)[i]*Ki*(dln_phi_liq_dT-dln_phi_vap_dT);
			}
			else{
				f += (*z)[i]*(1-1/Ki);
				dfdT += (*z)[i]/Ki*(dln_phi_liq_dT-dln_phi_vap_dT);
			}
		}

		change = -f/dfdT;

		T += -f/dfdT;
		if (type == TYPE_BUBBLEPOINT)
		{
			// Calculate the vapor molar fractions using the K factor and Rachford-Rice
			double sumy = 0;
			for (unsigned int i = 0; i < N; i++)
			{
				(*y)[i] = (*z)[i]*exp(ln_phi_liq[i])/exp(ln_phi_vap[i]);
				sumy += (*y)[i];
			}
			// Normalize the components
			for (unsigned int i = 0; i < N; i++)
			{
				(*y)[i] /= sumy;
			}
		}
		else
		{
			double sumx = 0;
			for (unsigned int i = 0; i < N; i++)
			{
				(*x)[i] = (*z)[i]/exp(ln_phi_liq[i])*exp(ln_phi_vap[i]);
				sumx += (*x)[i];
			}
			// Normalize the components
			for (unsigned int i = 0; i < N; i++)
			{
				(*x)[i] /= sumx;
			}
		}

		iter += 1;
		if (iter > 50)
		{
			return _HUGE;
			//throw ValueError(format("saturation_p was unable to reach a solution within 50 iterations"));
		}
	}
	while(abs(f) > 1e-8);
	//std::cout << iter << std::endl;
	return T;
}
void Mixture::TpzFlash(double T, double p, std::vector<double> *z, double *rhobar, std::vector<double> *x, std::vector<double> *y)
{
	unsigned int N = (*z).size();
	double beta, change;
	std::vector<double> lnK(N);
	
	(*x).resize(N);
	(*y).resize(N);

	// Wilson k-factors for each component
	for (unsigned int i = 0; i < N; i++)
	{
		lnK[i] = Wilson_lnK_factor(T, p, i);
	}

	// Check which phase we are in using Wilson estimations
	double g_RR_0 = g_RachfordRice(z, &lnK, 0);
	if (g_RR_0 < 0)
	{
		// Subcooled liquid - done
		*rhobar = rhobar_Tpz(T,p,z,rhobar_pengrobinson(T,p,z,PR_SATL));
		return;
	}
	else
	{
		double g_RR_1 = g_RachfordRice(z, &lnK, 1);
		if (g_RR_1 > 0)
		{
			// Superheated vapor - done
			*rhobar = rhobar_Tpz(T,p,z,rhobar_pengrobinson(T,p,z,PR_SATV));
			return;
		}
	}
	// TODO: How can you be sure that you aren't in the two-phase region? Safety factor needed?
	// TODO: Calculate the dewpoint density
	do
	{
	// Now find the value of beta that satisfies Rachford-Rice using Brent's method

	gRR_resid Resid(this,z,&lnK);
	std::string errstr;
	beta = Brent(&Resid,0,1,1e-16,1e-10,300,&errstr);

	// Evaluate mole fractions in liquid and vapor
	for (unsigned int i = 0; i < N; i++)
	{
		double Ki = exp(lnK[i]);
		double den = (1 - beta + beta*Ki); // Common denominator
		// Liquid mole fraction of component i
		(*x)[i] = (*z)[i]/den;
		// Vapor mole fraction of component i
		(*y)[i] = Ki*(*z)[i]/den;
	}

	// Reducing parameters for each phase
	double tau_liq = pReducing->Tr(x)/T;
	double tau_vap = pReducing->Tr(y)/T;
	
	double rhobar_liq = rhobar_Tpz(T, p, x, rhobar_pengrobinson(T,p,x,PR_SATL));
	double rhobar_vap = rhobar_Tpz(T, p, y, rhobar_pengrobinson(T,p,y,PR_SATV));
	double rhorbar_liq = pReducing->rhorbar(x);
	double rhorbar_vap = pReducing->rhorbar(y);
	double delta_liq = rhobar_liq/rhorbar_liq; 
	double delta_vap = rhobar_vap/rhorbar_vap; 

	double Z_liq = p/(Rbar*rhobar_liq*T);
	double Z_vap = p/(Rbar*rhobar_vap*T);
	
	double phir_liq =this->phir(tau_liq, delta_liq, x); 
	double phir_vap =this->phir(tau_vap, delta_vap, y);
	
	// Evaluate fugacity coefficients in liquid and vapor
	for (unsigned int i = 0; i < N; i++)
	{
		double ln_phi_liq = phir_liq + this->ndphir_dni(tau_liq,delta_liq,x,i)-log(Z_liq);
		double ln_phi_vap = phir_vap + this->ndphir_dni(tau_vap,delta_vap,y,i)-log(Z_vap);
		
		double lnKold = lnK[i];
		// Recalculate the K-factor (log(exp(ln_phi_liq)/exp(ln_phi_vap)))
		lnK[i] = ln_phi_liq - ln_phi_vap;

		change = lnK[i] - lnKold;
	}
	
	}
	while( abs(change) > 1e-7);

	return;
}
double Mixture::rhobar_pengrobinson(double T, double p, std::vector<double> *x, int solution)
{
	double k_ij = 0, A  = 0, B = 0, m_i, m_j, a_i, a_j, b_i, a = 0, b = 0, Z, rhobar;

	for (unsigned int i = 0; i < (*x).size(); i++)
	{
		m_i = 0.37464 + 1.54226*pFluids[i]->params.accentricfactor-0.26992*pow(pFluids[i]->params.accentricfactor,2);
		b_i = 0.077796074*(Rbar*pFluids[i]->reduce.T)/(pFluids[i]->reduce.p);

		B += (*x)[i]*b_i*p/(Rbar*T);

		for (unsigned int j = 0; j < (*x).size(); j++)
		{
			
			m_j = 0.37464 + 1.54226*pFluids[j]->params.accentricfactor-0.26992*pow(pFluids[j]->params.accentricfactor,2);
			a_i = 0.45724*pow(Rbar*pFluids[i]->reduce.T,2)/pFluids[i]->reduce.p*pow(1+m_i*(1-sqrt(T/pFluids[i]->reduce.T)),2)*1000;
			a_j = 0.45724*pow(Rbar*pFluids[j]->reduce.T,2)/pFluids[j]->reduce.p*pow(1+m_j*(1-sqrt(T/pFluids[j]->reduce.T)),2)*1000;	

			A += (*x)[i]*(*x)[j]*sqrt(a_i*a_j)*p/(Rbar*Rbar*T*T)/1000;
		}
	}

	std::vector<double> solns = solve_cubic(1, -1+B, A-3*B*B-2*B, -A*B+B*B+B*B*B);

	// Erase negative solutions and unstable solutions
	// Stable solutions are those for which dpdrho is positive
	for (int i = (int)solns.size()-1; i >= 0; i--)
	{
		
		if (solns[i] < 0)
		{
			solns.erase(solns.begin()+i);
		}
		else
		{
			double v = (solns[i]*Rbar*T)/p; //[mol/L]
			double dpdrho = -v*v*(-Rbar*T/pow(v-b,2)+a*(2*v+2*b)/pow(v*v+2*b*v-b*b,2));
			if (dpdrho < 0)
			{
				solns.erase(solns.begin()+i);
			}
		}
	}

	if (solution == PR_SATL)
	{
		Z = *std::min_element(solns.begin(), solns.end());
	}
	else if (solution == PR_SATV)
	{
		Z = *std::max_element(solns.begin(), solns.end());
	}
	else 
	{
		throw ValueError();
	}

	rhobar = p/(Z*Rbar*T);

	return rhobar;
}
double Mixture::g_RachfordRice(std::vector<double> *z, std::vector<double> *lnK, double beta)
{
	// g function from Rashford-Rice
	double summer = 0;
	for (unsigned int i = 0; i < (*z).size(); i++)
	{
		double Ki = exp((*lnK)[i]);
		summer += (*z)[i]*(Ki-1)/(1-beta+beta*Ki);
	}
	return summer;
}

double GERGReducingFunction::Tr(std::vector<double> *x)
{
	double Tr = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		double xi = (*x)[i], Tci = pFluids[i]->reduce.T;
		Tr += xi*xi*Tci;
		
		// The last term is only used for the pure component
		if (i==N-1){
			break;
		}

		for (unsigned int j = i+1; j < N; j++)
		{
			double xj = (*x)[j], beta_T_ij = beta_T[i][j];
			Tr += 2*xi*xj*beta_T_ij*gamma_T[i][j]*(xi+xj)/(beta_T_ij*beta_T_ij*xi+xj)*sqrt(Tci*pFluids[j]->reduce.T);
		}
	}
	return Tr;
}
double GERGReducingFunction::dTr_dxi(std::vector<double> *x, int i)
{
	// See Table B9 from Kunz Wagner 2012 (GERG 2008)
	double xi = (*x)[i];
	double dTr_dxi = 2*xi*pFluids[i]->reduce.T;
	for (int k = 0; k < i; k++)
	{
		double xk = (*x)[k], beta_T_ki = beta_T[i][k], gamma_T_ki = gamma_T[i][k];
		double Tr_ki = sqrt(pFluids[i]->reduce.T*pFluids[k]->reduce.T);
		double term = xk*(xk+xi)/(beta_T_ki*beta_T_ki*xk+xi)+xk*xi/(beta_T_ki*beta_T_ki*xk+xi)*(1-(xk+xi)/(beta_T_ki*beta_T_ki*xk+xi));
		dTr_dxi += 2*beta_T_ki*gamma_T_ki*Tr_ki*term;
	}
	for (unsigned int k = i+1; k < N; k++)
	{
		double xk = (*x)[k], beta_T_ik = beta_T[i][k], gamma_T_ik = gamma_T[i][k];
		double Tr_ik = sqrt(pFluids[i]->reduce.T*pFluids[k]->reduce.T);
		double term = xk*(xi+xk)/(beta_T_ik*beta_T_ik*xi+xk)+xi*xk/(beta_T_ik*beta_T_ik*xi+xk)*(1-beta_T_ik*beta_T_ik*(xi+xk)/(beta_T_ik*beta_T_ik*xi+xk));
		dTr_dxi += 2*beta_T_ik*gamma_T_ik*Tr_ik*term;
	}
	return dTr_dxi;
}
double GERGReducingFunction::drhorbar_dxi(std::vector<double> *x, int i)
{
	// See Table B9 from Kunz Wagner 2012 (GERG 2008)
	double xi = (*x)[i];
	double dvrbar_dxi = 2*xi/pFluids[i]->reduce.rhobar;
	for (int k = 0; k < i; k++)
	{
		double xk = (*x)[k], beta_v_ki = beta_v[i][k], gamma_v_ki = gamma_v[i][k];
		double vrbar_ki = 1.0/8.0*pow(pow(pFluids[i]->reduce.rhobar,-1.0/3.0) + pow(pFluids[k]->reduce.rhobar,-1.0/3.0),3.0);
		double term = xk*(xk+xi)/(beta_v_ki*beta_v_ki*xk+xi)+xk*xi/(beta_v_ki*beta_v_ki*xk+xi)*(1-(xk+xi)/(beta_v_ki*beta_v_ki*xk+xi));
		dvrbar_dxi += 2*beta_v_ki*gamma_v_ki*vrbar_ki*term;
	}
	for (unsigned int k = i+1; k < N; k++)
	{
		double xk = (*x)[k], beta_v_ik = beta_v[i][k], gamma_v_ik = gamma_v[i][k];
		double vrbar_ik = 1.0/8.0*pow(pow(pFluids[i]->reduce.rhobar,-1.0/3.0) + pow(pFluids[k]->reduce.rhobar,-1.0/3.0),3.0);
		double term = xk*(xi+xk)/(beta_v_ik*beta_v_ik*xi+xk)+xi*xk/(beta_v_ik*beta_v_ik*xi+xk)*(1-beta_v_ik*beta_v_ik*(xi+xk)/(beta_v_ik*beta_v_ik*xi+xk));
		dvrbar_dxi += 2*beta_v_ik*gamma_v_ik*vrbar_ik*term;
	}
	return -pow(rhorbar(x),2)*dvrbar_dxi;
}

double GERGReducingFunction::rhorbar(std::vector<double> *x)
{
	double vrbar = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		double xi = (*x)[i];
		vrbar += xi*xi/pFluids[i]->reduce.rhobar;

		if (i == N-1){ break;}

		for (unsigned int j = i+1; j < N; j++)
		{
			double xj = (*x)[j], beta_v_ij = beta_v[i][j];
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
double GERGDepartureFunction::phir(double tau, double delta, std::vector<double> *x)
{
	double term = phi1.base(tau, delta) + phi2.base(tau, delta);
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size()-1; i++)
	{
		for (unsigned int j = i + 1; j < (*x).size(); j++)
		{	
			summer += (*x)[i]*(*x)[j]*F[i][j]*term;
		}
	}
	return summer;
}
double GERGDepartureFunction::dphir_dDelta(double tau, double delta, std::vector<double> *x)
{
	double term = phi1.dDelta(tau, delta) + phi2.dDelta(tau, delta);
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size()-1; i++)
	{
		for (unsigned int j = i + 1; j < (*x).size(); j++)
		{	
			summer += (*x)[i]*(*x)[j]*F[i][j]*term;
		}
	}
	return summer;
}
double GERGDepartureFunction::d2phir_dDelta_dTau(double tau, double delta, std::vector<double> *x)
{
	double term = phi1.dDelta_dTau(tau, delta) + phi2.dDelta_dTau(tau, delta);
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size()-1; i++)
	{
		for (unsigned int j = i + 1; j < (*x).size(); j++)
		{	
			summer += (*x)[i]*(*x)[j]*F[i][j]*term;
		}
	}
	return summer;
}
double GERGDepartureFunction::dphir_dTau(double tau, double delta, std::vector<double> *x)
{
	double term = phi1.dTau(tau, delta) + phi2.dTau(tau, delta);
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size()-1; i++)
	{
		for (unsigned int j = i + 1; j < (*x).size(); j++)
		{	
			summer += (*x)[i]*(*x)[j]*F[i][j]*term;
		}
	}
	return summer;
}
double GERGDepartureFunction::d2phir_dTau2(double tau, double delta, std::vector<double> *x)
{
	double term = phi1.dTau2(tau, delta) + phi2.dTau2(tau, delta);
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size()-1; i++)
	{
		for (unsigned int j = i + 1; j < (*x).size(); j++)
		{	
			summer += (*x)[i]*(*x)[j]*F[i][j]*term;
		}
	}
	return summer;
}
double GERGDepartureFunction::dphir_dxi(double tau, double delta, std::vector<double> *x, int i)
{
	double summer = 0;
	double term = phi1.base(tau, delta) + phi2.base(tau, delta);
	for (unsigned int k = 0; k < (*x).size(); k++)
	{
		if (i != k)
		{
			summer += (*x)[k]*F[i][k]*term;
		}
	}
	return summer;
}
double GERGDepartureFunction::d2phir_dxi_dTau(double tau, double delta, std::vector<double> *x, int i)
{
	double summer = 0;
	double term = phi1.dTau(tau, delta) + phi2.dTau(tau, delta);
	for (unsigned int k = 0; k < (*x).size(); k++)
	{
		if (i != k)
		{
			summer += (*x)[k]*F[i][k]*term;
		}
	}
	return summer;
}

ResidualIdealMixture::ResidualIdealMixture(std::vector<Fluid*> pFluids)
{
	this->pFluids = pFluids;
}
double ResidualIdealMixture::phir(double tau, double delta, std::vector<double> *x)
{
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size(); i++)
	{
		summer += (*x)[i]*pFluids[i]->phir(tau,delta);
	}
	return summer;
}
double ResidualIdealMixture::dphir_dDelta(double tau, double delta, std::vector<double> *x)
{
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size(); i++)
	{
		summer += (*x)[i]*pFluids[i]->dphir_dDelta(tau,delta);
	}
	return summer;
}
double ResidualIdealMixture::d2phir_dDelta_dTau(double tau, double delta, std::vector<double> *x)
{
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size(); i++)
	{
		summer += (*x)[i]*pFluids[i]->d2phir_dDelta_dTau(tau,delta);
	}
	return summer;
}
double ResidualIdealMixture::dphir_dTau(double tau, double delta, std::vector<double> *x)
{
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size(); i++)
	{
		summer += (*x)[i]*pFluids[i]->dphir_dTau(tau,delta);
	}
	return summer;
}
double ResidualIdealMixture::d2phir_dTau2(double tau, double delta, std::vector<double> *x)
{
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size(); i++)
	{
		summer += (*x)[i]*pFluids[i]->d2phir_dTau2(tau,delta);
	}
	return summer;
}