
#include "Mixtures.h"
#include "Solvers.h"
#include "CPExceptions.h"
#include "mixture_excess_JSON.h" // Loads the JSON code for the excess parameters, and makes a variable "std::string mixture_excess_JSON"
#include "mixture_reducing_JSON.h" // Loads the JSON code for the reducing parameters, and makes a variable "std::string mixture_reducing_JSON"

std::vector<double> JSON_double_array(const rapidjson::Value& a)
{
	std::vector<double> v;
	for (rapidjson::Value::ConstValueIterator itrC = a.Begin(); itrC != a.End(); ++itrC)
	{
		v.push_back(itrC->GetDouble());
	}
	return v;
}
std::map<std::string, double> Mixture::load_reducing_values(int i, int j)
{
	// Make an output map that maps the necessary keys to the arrays of values
	std::map<std::string, double > outputmap;

	std::string Model;
	rapidjson::Document JSON;

	JSON.Parse<0>(mixture_reducing_JSON.c_str());
	
	// Iterate over the entries for the excess term
	for (rapidjson::Value::ConstValueIterator itr = JSON.Begin(); itr != JSON.End(); ++itr)
	{
		// Get the Model
		if (itr->HasMember("Model") && (*itr)["Model"].IsString())
		{
			Model = (*itr)["Model"].GetString();
		}
		else
		{
			throw ValueError("Model not included for reducing term");
		}

		// Get the coefficients
		if (itr->HasMember("Coeffs") && (*itr)["Coeffs"].IsArray())
		{
			const rapidjson::Value& Coeffs = (*itr)["Coeffs"];
			
			// Get the coeffs
			for (rapidjson::Value::ConstValueIterator itrC = Coeffs.Begin(); itrC != Coeffs.End(); ++itrC)
			{
				std::string Name1,Name2,CAS1,CAS2;
				std::vector<std::string> Names1, Names2, CASs1, CASs2;

				if (itrC->HasMember("Name1") && (*itrC)["Name1"].IsString() && itrC->HasMember("Name2") && (*itrC)["Name2"].IsString())
				{
					Name1 = (*itrC)["Name1"].GetString();
					Name2 = (*itrC)["Name2"].GetString();
					CAS1 = (*itrC)["CAS1"].GetString();
					CAS2 = (*itrC)["CAS2"].GetString();
				}
				else if (itrC->HasMember("Names1") && (*itrC)["Names1"].IsArray() && itrC->HasMember("Names2") && (*itrC)["Names2"].IsArray())
				{
					std::cout << format("Not currently supporting lists of components in excess terms\n").c_str();
					continue;
				}				
				std::string FluidiCAS = pFluids[i]->params.CAS, FluidjCAS = pFluids[j]->params.CAS;

				// Check if CAS codes match with either the i,j or j,i fluids
				if (
					(!(FluidiCAS.compare(CAS1)) && !(FluidjCAS.compare(CAS2)))
					||
					(!(FluidjCAS.compare(CAS1)) && !(FluidiCAS.compare(CAS2)))
				   )
				{
					// See if it is the GERG-2008 formulation
					if (!Model.compare("Kunz-JCED-2012"))
					{
						outputmap.insert(std::pair<std::string, double >("betaT",(*itrC)["betaT"].GetDouble()));
						outputmap.insert(std::pair<std::string, double >("betaV",(*itrC)["betaV"].GetDouble()));
						outputmap.insert(std::pair<std::string, double >("gammaT",(*itrC)["gammaT"].GetDouble()));
						outputmap.insert(std::pair<std::string, double >("gammaV",(*itrC)["gammaV"].GetDouble()));
						outputmap.insert(std::pair<std::string, double >("F",(*itrC)["F"].GetDouble()));

						return outputmap;
					}
					else
					{
						throw ValueError(format("This model [%s] is not currently supported\n",Model.c_str()).c_str());
					}
				}
			}
		}
		else
		{
			throw ValueError("Coeffs not included for reducing term");
		}
	}
	return outputmap;
};

std::map<std::string,std::vector<double> > Mixture::load_excess_values(int i, int j)
{
	// Make an output map that maps the necessary keys to the arrays of values
	std::map<std::string,std::vector<double> > outputmap;

	std::string Model;
	rapidjson::Document JSON;

	JSON.Parse<0>(mixture_excess_JSON.c_str());
	
	// Iterate over the entries for the excess term
	for (rapidjson::Value::ConstValueIterator itr = JSON.Begin(); itr != JSON.End(); ++itr)
	{
		// Get the Model
		if (itr->HasMember("Model") && (*itr)["Model"].IsString())
		{
			Model = (*itr)["Model"].GetString();
		}
		else
		{
			throw ValueError("Model not included for excess term");
		}

		// Get the coefficients
		if (itr->HasMember("Coeffs") && (*itr)["Coeffs"].IsArray())
		{
			const rapidjson::Value& Coeffs = (*itr)["Coeffs"];
			
			// Get the coeffs
			for (rapidjson::Value::ConstValueIterator itrC = Coeffs.Begin(); itrC != Coeffs.End(); ++itrC)
			{
				std::string Name1,Name2,CAS1,CAS2;
				std::vector<std::string> Names1, Names2, CASs1, CASs2;

				if (itrC->HasMember("Name1") && (*itrC)["Name1"].IsString() && itrC->HasMember("Name2") && (*itrC)["Name2"].IsString())
				{
					Name1 = (*itrC)["Name1"].GetString();
					Name2 = (*itrC)["Name2"].GetString();
					CAS1 = (*itrC)["CAS1"].GetString();
					CAS2 = (*itrC)["CAS2"].GetString();
				}
				else if (itrC->HasMember("Names1") && (*itrC)["Names1"].IsArray() && itrC->HasMember("Names2") && (*itrC)["Names2"].IsArray())
				{
					std::cout << format("Not currently supporting lists of components in excess terms\n").c_str();
					continue;
				}				
				std::string FluidiCAS = pFluids[i]->params.CAS, FluidjCAS = pFluids[j]->params.CAS;

				// Check if CAS codes match with either the i,j or j,i fluids
				if (
					(!(FluidiCAS.compare(CAS1)) && !(FluidjCAS.compare(CAS2)))
					||
					(!(FluidjCAS.compare(CAS1)) && !(FluidiCAS.compare(CAS2)))
				   )
				{
					// See if it is the GERG-2008 formulation
					if (!Model.compare("Kunz-JCED-2012"))
					{
						const rapidjson::Value& _n = (*itrC)["n"];
						outputmap.insert(std::pair<std::string,std::vector<double> >("n",JSON_double_array(_n)));
						const rapidjson::Value& _t = (*itrC)["t"];
						outputmap.insert(std::pair<std::string,std::vector<double> >("t",JSON_double_array(_t)));
						const rapidjson::Value& _d = (*itrC)["d"];
						outputmap.insert(std::pair<std::string,std::vector<double> >("d",JSON_double_array(_d)));
						const rapidjson::Value& _eta = (*itrC)["eta"];
						outputmap.insert(std::pair<std::string,std::vector<double> >("eta",JSON_double_array(_eta)));
						const rapidjson::Value& _epsilon = (*itrC)["epsilon"];
						outputmap.insert(std::pair<std::string,std::vector<double> >("epsilon",JSON_double_array(_epsilon)));
						const rapidjson::Value& _beta = (*itrC)["beta"];
						outputmap.insert(std::pair<std::string,std::vector<double> >("beta",JSON_double_array(_beta)));
						const rapidjson::Value& _gamma = (*itrC)["gamma"];
						outputmap.insert(std::pair<std::string,std::vector<double> >("gamma",JSON_double_array(_gamma)));

						return outputmap;
					}
					else
					{
						throw ValueError(format("This model [%s] is not currently supported\n",Model.c_str()).c_str());
					}
				}
			}
		}
		else
		{
			throw ValueError("Coeffs not included for excess term");
		}
	}
	return outputmap;
};

enum PengRobinsonOptions{PR_SATL, PR_SATV};
Mixture::Mixture(std::vector<Fluid *> pFluids)
{
	// Make a copy of the list of pointers to fluids that compose this mixture
	this->pFluids = pFluids;

	STLMatrix F;
	F.resize(pFluids.size(),std::vector<double>(pFluids.size(),1.0));

	// One reducing function for the entire mixture
	pReducing = new GERG2008ReducingFunction(pFluids);

	// One Residual Ideal Mixture term
	pResidualIdealMix = new ResidualIdealMixture(pFluids);

	// A matrix of departure functions for each binary pair
	pExcess = new ExcessTerm(pFluids.size());

	for (unsigned int i = 0; i < pFluids.size(); i++)
	{
		for (unsigned int j = 0; j < pFluids.size(); j++)
		{
			if (i != j)
			{
				pExcess->DepartureFunctionMatrix[i][j] = new GERG2008DepartureFunction();
				std::map<std::string, double > reducing_map = load_reducing_values(i, j);
				std::map<std::string, std::vector<double> > excess_map = load_excess_values(i, j);

				pExcess->set_coeffs_from_map(i, j, excess_map);
				pReducing->set_coeffs_from_map(i, j, reducing_map);
				// TODO: Load F factors into excess term
			}
		}
	}

	///         END OF INITIALIZATION
	///         END OF INITIALIZATION
	///         END OF INITIALIZATION
	///         END OF INITIALIZATION
	
	Rbar = 8.314491; // As in REFRPOP TODO: weight based on mole fractions of each component
	
	std::vector<double> z(2, 0.5);

	z[0] = 0.0;
	z[1] = 1-z[0];

	double Tr = pReducing->Tr(&z); //[K]
	double rhorbar = pReducing->rhorbar(&z); //[mol/m^3]
	double dTr_dxi = pReducing->dTr_dxi(&z,1);
	double rhorbar_dxi = pReducing->drhorbar_dxi(&z,1);

	double T = 145; // [K]
	double rhobar = 4; // [mol/m^3]; to convert from mol/L, multiply by 1000
	double tau = Tr/T;
	double delta = rhobar/rhorbar;

	double _dphir_dDelta = dphir_dDelta(tau,delta,&z);
	double p = Rbar*rhobar*T*(1 + delta*_dphir_dDelta);

	double rhobar_pr = rhobar_pengrobinson(T, p, &z, PR_SATV);

	//clock_t _t1,_t2;
	//long N = 1000000;
	//_t1 = clock();
	//for (unsigned long long i = 1; i<N; i++)
	//{
	//	z[0] = z[0]+1e-16*i;
	//	Tr = pReducing->Tr(&z);
	//}
	//_t2=clock();
	//double telap = ((double)(_t2-_t1))/CLOCKS_PER_SEC/(double)N*1e6;
	//std::cout << "elap" << telap << std::endl;

	//double tau_liq = 0.5;
	//double delta_liq  = 0.5;
	//double dd1 = (phir(tau_liq+0.0001,delta_liq,&z)- phir(tau_liq-0.0001,delta_liq,&z) )/(2*0.0001);
	//double dd2 = dphir_dTau(tau_liq+0.0001,delta_liq,&z);

	//p = 595.61824;
	std::vector<double> x,y;
	//TpzFlash(T, p, z, &rhobar, &x, &y);

	
	double Tsat;
	double psat = 1e3;

	for (double x0 = 0; x0 <= 1.000000000000001; x0 += 0.01)
	{
		z[0] = x0; z[1] = 1-x0;
		Tsat = saturation_p(TYPE_BUBBLEPOINT, psat, &z, &x, &y);
		std::cout << format("%g %g %g %g",x0,Tsat,y[0],y[1]).c_str();
		Tsat = saturation_p(TYPE_DEWPOINT, psat, &z, &x, &y);
		std::cout << format(" %g %g %g",Tsat,x[0],x[1]).c_str()	;
		std::cout << std::endl;
	}

	for (double p = 100; p <= 1e9; p *= 1.5)
	{
		double x0 = 0.5;
		z[0] = x0; z[1] = 1-x0;
		Tsat = saturation_p(TYPE_BUBBLEPOINT, p, &z, &x, &y);
		if (!ValidNumber(Tsat)){break;}
		std::cout << format("%g %g %g %g %g\n",x0,Tsat,y[0],y[1],p).c_str();
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
	double pci = pFluids[i]->reduce.p.Pa;
	double wi = pFluids[i]->params.accentricfactor;
	double Tci = pFluids[i]->reduce.T;
	return log(pci/p)+5.373*(1 + wi)*(1-Tci/T);
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
double Mixture::d2phir_dDelta2(double tau, double delta, std::vector<double> *x)
{
	return pResidualIdealMix->d2phir_dDelta2(tau,delta,x) + pExcess->d2phir_dDelta2(tau,delta,x);
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
	double deriv(double beta){return Mix->dgdbeta_RachfordRice(z, lnK, beta); };
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
		Rbar = Mix->Rbar; // J/mol/K
	};
	double call(double rhobar){	
		double delta = rhobar/rhorbar;
		double resid = Rbar*rhobar*T*(1 + delta*Mix->dphir_dDelta(tau, delta, x))-p;
		return resid;
	}
	double deriv(double rhobar){
		double delta = rhobar/rhorbar;
		double dDelta_drhobar = 1/rhorbar;
		double val = Rbar*T*(1 + 2*delta*Mix->dphir_dDelta(tau, delta, x)+delta*delta*Mix->d2phir_dDelta2(tau, delta, x));
		return val;
	}
};
double Mixture::rhobar_Tpz(double T, double p, std::vector<double> *x, double rhobar0)
{
	rho_Tpz_resid Resid(this,T,p,x);
	std::string errstr;
	return Newton(&Resid, rhobar0, 1e-7, 100, &errstr);
}



//double Mixture::saturation_p_NewtonRaphson(int type, double T, double p, std::vector<double> *z, std::vector<double> *ln_phi_liq, std::vector<double> *ln_phi_vap, std::vector<double> *x, std::vector<double> *y)
//{
//
//	
//}


/*! A wrapper function around the residual to find the initial guess for the bubble point temperature
\f[
r = \sum_i \left[z_i\cdot K_i\right] - 1 
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
		double summer = 1;
		for (unsigned int i = 0; i < (*z).size(); i++) 
		{ 
			double K = exp(Mix->Wilson_lnK_factor(T,p,i));
			summer -= (*z)[i]/K; 
		}
		return summer;
	};
};

double Mixture::saturation_p(int type, double p, std::vector<double> *z, std::vector<double> *x, std::vector<double> *y)
{
	int iter = 0;
	double change, T, f, dfdT, rhobar_liq_new, rhobar_vap_new;
	unsigned int N = (*z).size();
	std::vector<double> K(N), ln_phi_liq(N), ln_phi_vap(N);

	double pseudo_ptriple = 0;
	double pseudo_pcrit = 0;
	double pseudo_Ttriple = 0;
	double pseudo_Tcrit = 0;

	// Check if a pure fluid, if so, mole fractions and saturation temperatures are equal to bulk composition
	for (unsigned int i = 0; i < N; i++)
	{
		pseudo_ptriple += pFluids[i]->params.ptriple*(*z)[i];
		pseudo_pcrit += pFluids[i]->crit.p.Pa*(*z)[i];
		pseudo_Ttriple += pFluids[i]->params.Ttriple*(*z)[i];
		pseudo_Tcrit += pFluids[i]->crit.T*(*z)[i];

		if (fabs((*z)[i]-1) < 10*DBL_EPSILON)
		{
			(*x) = (*z);
			(*y) = (*z);
			// Get temperature from pure-fluid saturation call
			double TL=0, TV=0, rhoL=0, rhoV=0;
			pFluids[i]->saturation_p(p,false,&TL,&TV,&rhoL,&rhoV);
			return TL;
		}
	}

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
			} 
			catch (CoolPropBaseError) 
			{
				double eee = 0;
			}
		}
	}
	else if (type == TYPE_DEWPOINT)
	{
		// Vapor is at the bulk composition
		*y = (*z);
		// Find first guess for T using Wilson K-factors
		dewpoint_WilsonK_resid Resid(this,p,z); //1-sum(z_i/K_i)
		std::string errstr;
		double T_guess = (pseudo_Tcrit-pseudo_Ttriple)/(pseudo_pcrit-pseudo_ptriple)*(p-pseudo_ptriple)+pseudo_Ttriple;
		T = Secant(&Resid, T_guess, 0.001, 1e-10, 100, &errstr);
		if (!ValidNumber(T)){throw ValueError();}
	}
	else
	{
		throw ValueError("Invalid type to saturation_p");
	}

	double res_dew = 0;
	double res_bub = 0;
	// Calculate the K factors for each component
	for (unsigned int i = 0; i < N; i++)
	{
		K[i] = exp(Wilson_lnK_factor(T,p,i));
		res_dew += (*z)[i]/K[i];
		res_bub += (*z)[i]*K[i];
	}	
	res_dew = 1-res_dew;
	res_bub -= 1;

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
		double sumx = 0;
		// Calculate the liquid molar fractions using the K factor and Rachford-Rice
		for (unsigned int i = 0; i < N; i++)
		{
			(*x)[i] = (*z)[i]/K[i];
			sumx += (*x)[i];
		}
		// Normalize the components
		for (unsigned int i = 0; i < N; i++)
		{
			(*x)[i] /= sumx;
		}
	}

	double rhobar_liq = rhobar_pengrobinson(T,p,x,PR_SATL); // [kg/m^3]
	double rhobar_vap = rhobar_pengrobinson(T,p,y,PR_SATV); // [kg/m^3]
	do
	{
		//// Recalculate using Peng-Robinson
		//rhobar_liq = rhobar_pengrobinson(T,p,x,PR_SATL); // [kg/m^3]
		//rhobar_vap = rhobar_pengrobinson(T,p,y,PR_SATV); // [kg/m^3]

		rhobar_liq_new = rhobar_Tpz(T, p, x, rhobar_liq); // [kg/m^3]
		rhobar_vap_new = rhobar_Tpz(T, p, y, rhobar_vap); // [kg/m^3]
		if (!ValidNumber(rhobar_vap_new))
		{
			rhobar_vap = rhobar_Tpz(T, p, y, rhobar_pengrobinson(T,p,y,PR_SATV)); // [kg/m^3]
		}
		rhobar_liq = rhobar_liq_new;
		rhobar_vap = rhobar_vap_new;
		double Tr_liq = pReducing->Tr(x); // [K]
		double Tr_vap = pReducing->Tr(y);  // [K]
		double tau_liq = Tr_liq/T; // [-]
		double tau_vap = Tr_vap/T; // [-]
		double rhorbar_liq = pReducing->rhorbar(x); //[kg/m^3]
		double rhorbar_vap = pReducing->rhorbar(y); //[kg/m^3]
		double delta_liq = rhobar_liq/rhorbar_liq;  //[-]
		double delta_vap = rhobar_vap/rhorbar_vap;  //[-] 
		double dtau_dT_liq = -Tr_liq/T/T;
		double dtau_dT_vap = -Tr_vap/T/T;

		double Z_liq = p/(Rbar*rhobar_liq*T); //[-]
		double Z_vap = p/(Rbar*rhobar_vap*T); //[-]

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

		if (!ValidNumber(change))
		{
			double rr = 0;
		}
		T += change;
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
			// Calculate the liquid molar fractions using the K factor and Rachford-Rice
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
	beta = Newton(&Resid,0,1e-10,300,&errstr);

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
	double A  = 0, B = 0, m_i, m_j, a_i, a_j, b_i, a = 0, b = 0, Z, rhobar;

	for (unsigned int i = 0; i < (*x).size(); i++)
	{
		m_i = 0.37464 + 1.54226*pFluids[i]->params.accentricfactor-0.26992*pow(pFluids[i]->params.accentricfactor,2);
		b_i = 0.077796074*(Rbar*pFluids[i]->reduce.T)/(pFluids[i]->reduce.p.Pa);

		B += (*x)[i]*b_i*p/(Rbar*T);

		for (unsigned int j = 0; j < (*x).size(); j++)
		{
			
			m_j = 0.37464 + 1.54226*pFluids[j]->params.accentricfactor-0.26992*pow(pFluids[j]->params.accentricfactor,2);
			a_i = 0.45724*pow(Rbar*pFluids[i]->reduce.T,2)/pFluids[i]->reduce.p.Pa*pow(1+m_i*(1-sqrt(T/pFluids[i]->reduce.T)),2)*1000;
			a_j = 0.45724*pow(Rbar*pFluids[j]->reduce.T,2)/pFluids[j]->reduce.p.Pa*pow(1+m_j*(1-sqrt(T/pFluids[j]->reduce.T)),2)*1000;	

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
	// g function from Rachford-Rice
	double summer = 0;
	for (unsigned int i = 0; i < (*z).size(); i++)
	{
		double Ki = exp((*lnK)[i]);
		summer += (*z)[i]*(Ki-1)/(1-beta+beta*Ki);
	}
	return summer;
}
double Mixture::dgdbeta_RachfordRice(std::vector<double> *z, std::vector<double> *lnK, double beta)
{
	// derivative of g function from Rachford-Rice with respect to beta
	double summer = 0;
	for (unsigned int i = 0; i < (*z).size(); i++)
	{
		double Ki = exp((*lnK)[i]);
		summer += -(*z)[i]*pow((Ki-1)/(1-beta+beta*Ki),2);
	}
	return summer;
}

double GERG2008ReducingFunction::Tr(std::vector<double> *x)
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
double GERG2008ReducingFunction::dTr_dxi(std::vector<double> *x, int i)
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
double GERG2008ReducingFunction::drhorbar_dxi(std::vector<double> *x, int i)
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

double GERG2008ReducingFunction::rhorbar(std::vector<double> *x)
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
void GERG2008ReducingFunction::set_coeffs_from_map(int i, int j, std::map<std::string,double > m)
{
	beta_v[i][j] = m.find("betaV")->second;
	beta_T[i][j] = m.find("betaT")->second;
	gamma_v[i][j] = m.find("gammaV")->second;
	gamma_T[i][j] = m.find("gammaT")->second;

	//std::map<std::string,std::vector<double> >::iterator it;
	//// Try to find using the ma
	//it = m.find("t")->second;
	//// If it is found the iterator will not be equal to end
	//if (it != param_map.end() )
	//{

	//}
	//// Try to find using the map
	//it = m.find("n");
	//// If it is found the iterator will not be equal to end
	//if (it != param_map.end() )
	//{

	//}
}

//GERG2008DepartureFunction::GERG2008DepartureFunction()
//{
//	if (phi1 != NULL){
//		delete phi1; phi1 = NULL;
//	}
//	if (phi2 != NULL){
//		delete phi2; phi2 = NULL;
//	}
//}
double GERG2008DepartureFunction::phir(double tau, double delta, std::vector<double> *x)
{
	return phi1.base(tau, delta) + phi2.base(tau, delta);
}
double GERG2008DepartureFunction::dphir_dDelta(double tau, double delta, std::vector<double> *x)
{
	return phi1.dDelta(tau, delta) + phi2.dDelta(tau, delta);
}
double GERG2008DepartureFunction::d2phir_dDelta2(double tau, double delta, std::vector<double> *x)
{
	return phi1.dDelta2(tau, delta) + phi2.dDelta2(tau, delta);
}
double GERG2008DepartureFunction::d2phir_dDelta_dTau(double tau, double delta, std::vector<double> *x)
{
	return phi1.dDelta_dTau(tau, delta) + phi2.dDelta_dTau(tau, delta);
}
double GERG2008DepartureFunction::dphir_dTau(double tau, double delta, std::vector<double> *x)
{
	return phi1.dTau(tau, delta) + phi2.dTau(tau, delta);
}
double GERG2008DepartureFunction::d2phir_dTau2(double tau, double delta, std::vector<double> *x)
{
	return phi1.dTau2(tau, delta) + phi2.dTau2(tau, delta);
}

void GERG2008DepartureFunction::set_coeffs_from_map(std::map<std::string,std::vector<double> > m)
{
	std::vector<double> n = m.find("n")->second;
	std::vector<double> t = m.find("t")->second;
	std::vector<double> d = m.find("d")->second;
	std::vector<double> eta = m.find("eta")->second;
	std::vector<double> epsilon = m.find("epsilon")->second;
	std::vector<double> beta = m.find("beta")->second;
	std::vector<double> gamma = m.find("gamma")->second;

	int iStart = 0;
	while (eta[iStart+1] < 0.5)
	{
		iStart++;
	}

	phi1 = phir_power(n,d,t,1,iStart);
	phi2 = phir_GERG2008_gaussian(n,d,t,eta,epsilon,beta,gamma,iStart+1,eta.size()-1);

	//std::map<std::string,std::vector<double> >::iterator it;
	//// Try to find using the ma
	//it = m.find("t")->second;
	//// If it is found the iterator will not be equal to end
	//if (it != param_map.end() )
	//{

	//}
	//// Try to find using the map
	//it = m.find("n");
	//// If it is found the iterator will not be equal to end
	//if (it != param_map.end() )
	//{

	//}
}

ExcessTerm::ExcessTerm(int N)
{
	F.resize(N,std::vector<double>(N,0));
	DepartureFunctionMatrix.resize(N,std::vector<DepartureFunction*>(N,NULL));
	this->N = N;
}
ExcessTerm::~ExcessTerm()
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{	
			if (DepartureFunctionMatrix[i][j] != NULL)
			{
				delete DepartureFunctionMatrix[i][j]; DepartureFunctionMatrix[i][j] = NULL;
			}
		}
	}
}

double ExcessTerm::phir(double tau, double delta, std::vector<double> *x)
{
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size()-1; i++)
	{
		for (unsigned int j = i + 1; j < (*x).size(); j++)
		{	
			summer += (*x)[i]*(*x)[j]*F[i][j]*DepartureFunctionMatrix[i][j]->phir(tau,delta,x);
		}
	}
	return summer;
}
double ExcessTerm::dphir_dTau(double tau, double delta, std::vector<double> *x)
{
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size()-1; i++)
	{
		for (unsigned int j = i + 1; j < (*x).size(); j++)
		{	
			summer += (*x)[i]*(*x)[j]*F[i][j]*DepartureFunctionMatrix[i][j]->dphir_dTau(tau,delta,x);
		}
	}
	return summer;
}
double ExcessTerm::dphir_dDelta(double tau, double delta, std::vector<double> *x)
{
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size()-1; i++)
	{
		for (unsigned int j = i + 1; j < (*x).size(); j++)
		{	
			summer += (*x)[i]*(*x)[j]*F[i][j]*DepartureFunctionMatrix[i][j]->dphir_dDelta(tau,delta,x);
		}
	}
	return summer;
}
double ExcessTerm::d2phir_dDelta2(double tau, double delta, std::vector<double> *x)
{
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size()-1; i++)
	{
		for (unsigned int j = i + 1; j < (*x).size(); j++)
		{	
			summer += (*x)[i]*(*x)[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2phir_dDelta2(tau,delta,x);
		}
	}
	return summer;
}
double ExcessTerm::d2phir_dTau2(double tau, double delta, std::vector<double> *x)
{
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size()-1; i++)
	{
		for (unsigned int j = i + 1; j < (*x).size(); j++)
		{	
			summer += (*x)[i]*(*x)[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2phir_dTau2(tau,delta,x);
		}
	}
	return summer;
}
double ExcessTerm::d2phir_dDelta_dTau(double tau, double delta, std::vector<double> *x)
{
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size()-1; i++)
	{
		for (unsigned int j = i + 1; j < (*x).size(); j++)
		{	
			summer += (*x)[i]*(*x)[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2phir_dDelta_dTau(tau,delta,x);
		}
	}
	return summer;
}
double ExcessTerm::dphir_dxi(double tau, double delta, std::vector<double> *x, unsigned int i)
{
	double summer = 0;
	for (unsigned int k = 0; k < (*x).size(); k++)
	{
		if (i != k)
		{
			summer += (*x)[k]*F[i][k]*DepartureFunctionMatrix[i][k]->phir(tau,delta,x);
		}
	}
	return summer;
}
double ExcessTerm::d2phir_dxi_dTau(double tau, double delta, std::vector<double> *x, unsigned int i)
{
	double summer = 0;
	for (unsigned int k = 0; k < (*x).size(); k++)
	{
		if (i != k)
		{
			summer += (*x)[k]*F[i][k]*DepartureFunctionMatrix[i][k]->dphir_dTau(tau,delta,x);
		}
	}
	return summer;
}
double ExcessTerm::d2phir_dxi_dDelta(double tau, double delta, std::vector<double> *x, unsigned int i)
{
	double summer = 0;
	for (unsigned int k = 0; k < (*x).size(); k++)
	{
		if (i != k)
		{
			summer += (*x)[k]*F[i][k]*DepartureFunctionMatrix[i][k]->dphir_dDelta(tau,delta,x);
		}
	}
	return summer;
}
void ExcessTerm::set_coeffs_from_map(int i, int j, std::map<std::string, std::vector<double> > m)
{
	DepartureFunctionMatrix[i][j]->set_coeffs_from_map(m);
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
double ResidualIdealMixture::d2phir_dDelta2(double tau, double delta, std::vector<double> *x)
{
	double summer = 0;
	for (unsigned int i = 0; i < (*x).size(); i++)
	{
		summer += (*x)[i]*pFluids[i]->d2phir_dDelta2(tau,delta);
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