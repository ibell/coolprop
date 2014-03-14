
#define _CRT_SECURE_NO_WARNINGS

#include "rapidjson_CoolProp.h"

#include "Mixtures.h"
#include "Solvers.h"
#include "CPExceptions.h"
#include "MatrixMath.h"
#include "mixture_excess_JSON.h" // Loads the JSON code for the excess parameters, and makes a variable "std::string mixture_excess_JSON"
#include "mixture_reducing_JSON.h" // Loads the JSON code for the reducing parameters, and makes a variable "std::string mixture_reducing_JSON"
#include <numeric>
#include "CoolProp.h"
#include "Spline.h"

#ifdef MPLSUPPORTED
#include "MPLPlot.h"
#endif

static const bool use_cache = true;

bool has_string_array_member(const rapidjson::Value& a, const char * member)
{
	if(a.HasMember(member) && a[member].IsArray())
	{
		for (rapidjson::Value::ConstValueIterator itrC = a[member].Begin(); itrC != a[member].End(); ++itrC)
		{	
			if (!itrC->IsString())
			{
				return false;
			}
		}
		return true;
	}
	else
	{
		return false;
	}
}

std::vector<double> JSON_double_array(const rapidjson::Value& a)
{
	std::vector<double> v;
	for (rapidjson::Value::ConstValueIterator itrC = a.Begin(); itrC != a.End(); ++itrC)
	{
		if (itrC->IsInt())
		{
			v.push_back((double)itrC->GetInt());
		}
		else
		{
			v.push_back(itrC->GetDouble());
		}
	}
	return v;
}
std::vector<std::string> JSON_string_array(const rapidjson::Value& a)
{
	std::vector<std::string> v;
	for (rapidjson::Value::ConstValueIterator itrC = a.Begin(); itrC != a.End(); ++itrC)
	{
		v.push_back(itrC->GetString());
	}
	return v;
}
bool match_CAS_entry(std::vector<std::string> CAS1, std::vector<std::string> CAS2, std::string FluidiCAS, std::string FluidjCAS)
{
	for (unsigned int k = 0; k < CAS1.size(); k++)
	{
		if (
			 (!FluidiCAS.compare(CAS1[k]) && !FluidjCAS.compare(CAS2[k]))
			 ||
			 (!FluidjCAS.compare(CAS1[k]) && !FluidiCAS.compare(CAS2[k]))
			){ 	return true; }
	}
	return false;
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
						if (!pReducing) 
						{ 
							// One reducing function for the entire mixture
							pReducing = new GERG2008ReducingFunction(pFluids);
						}
						if (!(FluidiCAS.compare(CAS1)) && !(FluidjCAS.compare(CAS2)))
						{
							outputmap.insert(std::pair<std::string, double >("betaT",(*itrC)["betaT"].GetDouble()));
							outputmap.insert(std::pair<std::string, double >("betaV",(*itrC)["betaV"].GetDouble()));
						}
						else
						{
							outputmap.insert(std::pair<std::string, double >("betaT",1/(*itrC)["betaT"].GetDouble()));
							outputmap.insert(std::pair<std::string, double >("betaV",1/(*itrC)["betaV"].GetDouble()));
						}
						outputmap.insert(std::pair<std::string, double >("gammaT",(*itrC)["gammaT"].GetDouble()));
						outputmap.insert(std::pair<std::string, double >("gammaV",(*itrC)["gammaV"].GetDouble()));
						outputmap.insert(std::pair<std::string, double >("F",(*itrC)["F"].GetDouble()));

						return outputmap;
					}
					else if (!Model.compare("Lemmon-JPCRD-2000") || !Model.compare("Lemmon-JPCRD-2004"))
					{
						if (!pReducing) 
						{ 
							// One reducing function for the entire mixture
							pReducing = new LemmonAirHFCReducingFunction(pFluids);
						}
						outputmap.insert(std::pair<std::string,double>("xi",(*itrC)["xi"].GetDouble()));
						outputmap.insert(std::pair<std::string,double>("zeta",(*itrC)["zeta"].GetDouble()));
						outputmap.insert(std::pair<std::string,double>("F",(*itrC)["F"].GetDouble()));		
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

void Mixture::load_excess_values(int i, int j)
{
	// Make an output map that maps the necessary keys to the arrays of values
	std::map<std::string,std::vector<double> > outputmap;

	std::string FluidiCAS = pFluids[i]->params.CAS, 
		        FluidjCAS = pFluids[j]->params.CAS;

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
			// It doesn't have the term Model
			throw ValueError("Model not included for excess term");
		}

		// Get the coefficients
		if (itr->HasMember("Coeffs") && (*itr)["Coeffs"].IsArray())
		{
			const rapidjson::Value& Coeffs = (*itr)["Coeffs"];
			
			// Get the coeffs
			for (rapidjson::Value::ConstValueIterator itrC = Coeffs.Begin(); itrC != Coeffs.End(); ++itrC)
			{
				std::vector<std::string> Names1, Names2, CASs1, CASs2;

				if (itrC->HasMember("Name1") && (*itrC)["Name1"].IsString() && itrC->HasMember("Name2") && (*itrC)["Name2"].IsString())
				{
					Names1 = std::vector<std::string>(1, (*itrC)["Name1"].GetString());
					Names2 = std::vector<std::string>(1, (*itrC)["Name2"].GetString());
					CASs1 = std::vector<std::string>(1, (*itrC)["CAS1"].GetString());
					CASs2 = std::vector<std::string>(1, (*itrC)["CAS2"].GetString());
				}
				else if (has_string_array_member((*itrC), "Name1") && has_string_array_member((*itrC), "Name2") 
					     && has_string_array_member((*itrC), "CAS1") && has_string_array_member((*itrC), "CAS2"))
				{
					Names1 = JSON_string_array((*itrC)["Name1"]);
					Names2 = JSON_string_array((*itrC)["Name2"]);
					CASs1 = JSON_string_array((*itrC)["CAS1"]);
					CASs2 = JSON_string_array((*itrC)["CAS2"]);
				}			

				// Check if CAS codes match with either the i,j or j,i fluids
				if (match_CAS_entry(CASs1,CASs2,FluidiCAS,FluidjCAS))
				{
					// See if it is the GERG-2008 formulation
					if (!Model.compare("Kunz-JCED-2012"))
					{
						// Create an instance for the departure function for GERG formulation
						pExcess->DepartureFunctionMatrix[i][j] = new GERG2008DepartureFunction();

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

						// Set the variables in the class
						pExcess->DepartureFunctionMatrix[i][j]->set_coeffs_from_map(outputmap);
						return;
					}
					else if (!Model.compare("Lemmon-JPCRD-2004"))
					{
						// Create an instance for the departure function for the HFC mixtures
						pExcess->DepartureFunctionMatrix[i][j] = new LemmonHFCDepartureFunction();

						const rapidjson::Value& _n = (*itrC)["n"];
						outputmap.insert(std::pair<std::string,std::vector<double> >("n",JSON_double_array(_n)));
						const rapidjson::Value& _t = (*itrC)["t"];
						outputmap.insert(std::pair<std::string,std::vector<double> >("t",JSON_double_array(_t)));
						const rapidjson::Value& _d = (*itrC)["d"];
						outputmap.insert(std::pair<std::string,std::vector<double> >("d",JSON_double_array(_d)));
						const rapidjson::Value& _l = (*itrC)["l"];
						outputmap.insert(std::pair<std::string,std::vector<double> >("l",JSON_double_array(_l)));

						// Set the variables in the class
						pExcess->DepartureFunctionMatrix[i][j]->set_coeffs_from_map(outputmap);
						return;
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
	throw ValueError("No excess parameters loaded for this binary pair");
};

void normalize_vector(std::vector<double> &x)
{
	double sumx = std::accumulate( x.begin(), x.end(), (double)0 );
	// Normalize the components
	for (unsigned int i = 0; i < x.size(); i++)
	{
		x[i] /= sumx;
	}
};

enum PengRobinsonOptions{PR_SATL, PR_SATV};

double Mixture::Rbar(const std::vector<double> &x)
{
	double R = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		R += x[i]*pFluids[i]->params.R_u;
	}
	return R;
}
Mixture::Mixture(std::vector<Fluid *> pFluids)
{
	// Make a copy of the list of pointers to fluids that compose this mixture
	this->pFluids = pFluids;
	this->N = pFluids.size();
	
	// Finish the initalization
	this->initialize();

}
Mixture::Mixture(std::string FluidsString)
{
	// Split at the '|'
	std::vector<std::string> fluids = strsplit(FluidsString,'|');

	// Number of components
	this->N = fluids.size();

	// Resize the pointers
	pFluids.resize(this->N);

	// Get all the pointers to the components
	for (unsigned int i = 0; i < fluids.size(); i++)
	{
		pFluids[i] = get_fluid(get_Fluid_index(fluids[i]));
	}

	// Finish the initalization
	this->initialize();
}

void Mixture::initialize()
{
	// Give child classes a pointer to this class
	SS.Mix = this;
	NRVLE.Mix = this;
	Envelope.Mix = this;

	STLMatrix F;
	F.resize(pFluids.size(),std::vector<double>(pFluids.size(),1.0));

	NRVLE.resize(this->N);

	// Reset the reducing parameter pointer
	pReducing = NULL;

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
				std::map<std::string, double > reducing_map = load_reducing_values(i, j);
				load_excess_values(i, j);
				pReducing->set_coeffs_from_map(i, j, reducing_map);
				pExcess->F[i][j] = reducing_map.find("F")->second;
			}
		}
	}
	
	///         END OF INITIALIZATION
	///         END OF INITIALIZATION
	///         END OF INITIALIZATION
	///         END OF INITIALIZATION
}

void Mixture::test()
{
	std::vector<double> x,y,z(2, 0.5);

	//Envelope.build(100000,z);

	z[0] = 0.3;
	z[1] = 1-z[0];

	double _Tsat = saturation_p(0, 440000, z, x, y);

	double Tr = pReducing->Tr(z); //[K]
	double rhorbar = pReducing->rhorbar(z); //[mol/m^3]

	double T = 145; // [K]
	double rhobar = 4; // [mol/m^3]; to convert from mol/L, multiply by 1000
	double tau = Tr/T;
	double delta = rhobar/rhorbar;
	//double dtau_dT = -Tr/T/T;

	//double _dphir_dDelta = dphir_dDelta(tau, delta, &z);
	//double p = Rbar(&z)*rhobar*T*(1 + delta*_dphir_dDelta)/1000; //[kPa]

	
	//check();
	//time_t t1,t2;
	//unsigned int N = 100;
	//double TTE = 0;
	//t1 = clock();
	//for (unsigned int i = 0; i < N; i++)
	//{
	//	TTE += saturation_p(TYPE_BUBBLEPOINT,10000,z,x,y);
	//	//TTE += Props("T",'P',10,'Q',0,"REFPROP-MIX:Methane[0.5]&Ethane[0.5]");
	//}
	//t2 = clock();
	//printf("val %g elapsed time/call %g us\n", TTE/(double)N, (double)(t2-t1)/(double)CLOCKS_PER_SEC/N*1e6);
//	return;
	//p = 595.61824;
	
	//TpzFlash(T, p, z, &rhobar, &x, &y);

	//double rhor = pReducing->rhorbar(&z);
	
	double Tsat;

	for (double psat = 1e6; psat < 6.9e6; psat += 1e6)
	{
		std::cout << "psat: " << psat << " Pa" << std::endl;
		for (double x0 = 0; x0 <= 1.000000000000001; x0 += 0.01)
		{
			z[0] = x0; z[1] = 1-x0;
			Tsat = saturation_p(0, psat, z, x, y);
			std::cout << format("%g %g %0.9g %0.9g",x0,Tsat,y[0],y[1]).c_str();
			Tsat = saturation_p(1, psat, z, x, y);
			std::cout << format(" %g %0.9g %0.9g",Tsat,x[0],x[1]).c_str()	;
			std::cout << std::endl;
		}
	}

	/*double x0 = 0.5;
	z[0] = x0; z[1] = 1-x0;
	double TL, TV;
	for (double p = 100; p <= 1e9; p *= 1.5)
	{
		TL = saturation_p(TYPE_BUBBLEPOINT, p, z, x, y);
		TV = saturation_p(TYPE_DEWPOINT, p, z, x, y);
		if (!ValidNumber(Tsat)){break;}
		std::cout << format("%g %g %g\n",TL,TV,p).c_str();
	}*/
}

void Mixture::check_WaterEthanol()
{
	std::vector<double> z(2, 0.5);

	double tol = 1e-10; // relative tolerance
	z[0] = 0.4;
	z[1] = 1-z[0];

	double Tr = pReducing->Tr(z);
	double Tr_RP91 = 574.10640748981632; // [K]
	if (fabs(Tr/Tr_RP91-1) > tol){throw ValueError();}

	double rhorbar = pReducing->rhorbar(z);
	double rhorbar_RP91 = 11090.627577544581; //[mol/m^3]
	if (fabs(rhorbar/rhorbar_RP91-1) > tol){throw ValueError();}

	double T = 500; // [K]
	double rhobar = 868.40163264165199; // [mol/m^3]; to convert from mol/L, multiply by 1000
	
	double tau = Tr/T;
	double delta = rhobar/rhorbar;

	double RT = Rbar(z)*T;
	double p = rhobar*RT*(1+delta*dphir_dDelta(tau, delta,z));
	double p_RP91 = 3011336.2097271665;
	if (fabs(p/p_RP91-1) > tol){throw ValueError();}

	double dtdn0 = pReducing->ndTrdni__constnj(z,0);
	double dtdn0_RP91 = -80.606224987885071;
	if (fabs(dtdn0/dtdn0_RP91-1) > tol){throw ValueError();}

	double dtdn1 = pReducing->ndTrdni__constnj(z,1);
	double dtdn1_RP91 = 53.737483325256562;
	if (fabs(dtdn1/dtdn1_RP91-1) > tol){throw ValueError();}

	double drhodn0 = pReducing->ndrhorbardni__constnj(z,0);
	double drhodn0_RP91 = -7259.8385028178765;
	if (fabs(drhodn0/drhodn0_RP91-1) > tol){throw ValueError();}

	double drhodn1 = pReducing->ndrhorbardni__constnj(z,1);
	double drhodn1_RP91 = 4839.8923352119168;
	if (fabs(drhodn1/drhodn1_RP91-1) > tol){throw ValueError();}

	double ddrdxn00 = pReducing->d_ndrhorbardni_dxj__constxi(z,0,0);
	double ddrdxn00_RP91 = 57050.801427407002;
	if (fabs(ddrdxn00/ddrdxn00_RP91-1) > tol){throw ValueError();}

	double ddrdxn01 = pReducing->d_ndrhorbardni_dxj__constxi(z,0,1);
	double ddrdxn01_RP91 = 35234.083487633299;
	if (fabs(ddrdxn01/ddrdxn01_RP91-1) > tol){throw ValueError();}

	double ddrdxn10 = pReducing->d_ndrhorbardni_dxj__constxi(z,1,0);
	double ddrdxn10_RP91 = 11034.621811573714;
	if (fabs(ddrdxn10/ddrdxn10_RP91-1) > tol){throw ValueError();}

	double ddrdxn11 = pReducing->d_ndrhorbardni_dxj__constxi(z,1,1);
	double ddrdxn11_RP91 = 5412.8823747065340;
	if (fabs(ddrdxn11/ddrdxn11_RP91-1) > tol){throw ValueError();}

	double dtrdxn00 = pReducing->d_ndTrdni_dxj__constxi(z,0,0);
	double dtrdxn00_RP91 = -1078.6599487910798;
	if (fabs(dtrdxn00/dtrdxn00_RP91-1) > tol){throw ValueError();}

	double dtrdxn01 = pReducing->d_ndTrdni_dxj__constxi(z,0,1);
	double dtrdxn01_RP91 = -1328.9251007518092;
	if (fabs(dtrdxn01/dtrdxn01_RP91-1) > tol){throw ValueError();}
	
	double dtrdxn10 = pReducing->d_ndTrdni_dxj__constxi(z,1,0);
	double dtrdxn10_RP91 = -1060.2376841255259;
	if (fabs(dtrdxn10/dtrdxn10_RP91-1) > tol){throw ValueError();}
	
	double dtrdxn11 = pReducing->d_ndTrdni_dxj__constxi(z,1,1);
	double dtrdxn11_RP91 = -1117.3004300069424;
	if (fabs(dtrdxn11/dtrdxn11_RP91-1) > tol){throw ValueError();}

	/*double dadxi0 = this->dphir_dxi(tau, delta, z, 0);
	double dadxi0_RP91 = -1.66733159546361936E-003;
	if (fabs(dadxi0/dadxi0_RP91-1) > tol){throw ValueError();}

	double dadxi1 = this->dphir_dxi(tau, delta, z, 1);
	double dadxi1_RP91 = -1.82138497681156902E-003;
	if (fabs(dadxi1/dadxi1_RP91-1) > tol){throw ValueError();}*/

	double daddx0 = this->d2phir_dxi_dDelta(tau,delta,z,0);
	double daddx0_RP91 = -1.9898117177651518;
	if (fabs(daddx0/daddx0_RP91-1) > tol){throw ValueError();}

	double daddx1 = this->d2phir_dxi_dDelta(tau,delta,z,1);
	double daddx1_RP91 = -2.1704856931980134;
	if (fabs(daddx1/daddx1_RP91-1) > tol){throw ValueError();}

	double dadtx0 = this->d2phir_dxi_dTau(tau,delta,z,0);
	double dadtx0_RP91 = -0.54706966169729132;
	if (fabs(dadtx0/dadtx0_RP91-1) > tol){throw ValueError();}

	double dadtx1 = this->d2phir_dxi_dTau(tau,delta,z,1);
	double dadtx1_RP91 = -0.45648092575441646;
	if (fabs(dadtx1/dadtx1_RP91-1) > tol){throw ValueError();}

	double dadn0 = this->ndphir_dni__constT_V_nj(tau, delta, z, 0);
	double dadn0_RP91 = -0.18827619536857490;
	if (fabs(dadn0/dadn0_RP91-1) > tol){throw ValueError();}

	double dadn1 = this->ndphir_dni__constT_V_nj(tau, delta, z, 1);
	double dadn1_RP91 = -0.15092181712014577;
	if (fabs(dadn1/dadn1_RP91-1) > tol){throw ValueError();}

	double dpdn0 = ndpdni__constT_V_nj(tau, delta, z, 0);
	double dpdn0_RP91 = 2359472.7345531628;
	if (fabs(dpdn0/dpdn0_RP91-1) > tol){throw ValueError();}

	double dpdn1 = ndpdni__constT_V_nj(tau, delta, z, 1);
	double dpdn1_RP91 = 2459536.5646813584;
	if (fabs(dpdn1/dpdn1_RP91-1) > tol){throw ValueError();}

	double vhat0 = partial_molar_volume(tau, delta, z, 0);
	double vhat0_RP91 = 1.1229663043954532/1000;
	if (fabs(vhat0/vhat0_RP91-1) > tol){throw ValueError();}
	
	double vhat1 = partial_molar_volume(tau, delta, z, 1);
	double vhat1_RP91 = 1.1705906349830217/1000;
	if (fabs(vhat1/vhat1_RP91-1) > tol){throw ValueError();}

	/*double d2adbn0 = d2nphir_dni_dT(tau, delta, z, 0);
	double d2adbn0_RP91 = 2.91400791470335643E-005;
	if (fabs(d2adbn0/d2adbn0_RP91-1) > tol){throw ValueError();}
	
	double d2adbn1 = d2nphir_dni_dT(tau, delta, z, 1);
	double d2adbn1_RP91 = 6.58714026367381391E-005;
	if (fabs(d2adbn1/d2adbn1_RP91-1) > tol){throw ValueError();}*/

	double dphidT0 = dln_fugacity_coefficient_dT__constp_n(tau,delta,z,0);
	double dphidT0_RP91 = 1.80275151505261211E-003;
	if (fabs(dphidT0/dphidT0_RP91-1) > tol){throw ValueError();}

	double dphidT1 = dln_fugacity_coefficient_dT__constp_n(tau,delta,z,1);
	double dphidT1_RP91 =  1.24216502708919410E-003;
	if (fabs(dphidT1/dphidT1_RP91-1) > tol){throw ValueError();}

	double dphidP0 = dln_fugacity_coefficient_dp__constT_n(tau,delta,z,0);
	double dphidP0_RP91 = -6.19532350116346726E-005/1000;
	if (fabs(dphidP0/dphidP0_RP91-1) > tol){throw ValueError();}

	double dphidP1 = dln_fugacity_coefficient_dp__constT_n(tau,delta,z,1);
	double dphidP1_RP91 = -5.04973839435411730E-005/1000;
	if (fabs(dphidP1/dphidP1_RP91-1) > tol){throw ValueError();}

	double dadxij00 = d2phirdxidxj(tau, delta, z, 0, 0);
	double dadxij00_RP91 = 0.0;
	if (fabs(dadxij00-dadxij00_RP91) > tol){throw ValueError();}

	double dadxij01 = d2phirdxidxj(tau, delta, z, 0, 1);
	double dadxij01_RP91 = 2.95156019175951646E-003;
	if (fabs(dadxij01/dadxij01_RP91-1) > tol){throw ValueError();}
	
	double dadxij10 = d2phirdxidxj(tau, delta, z, 1, 0);
	double dadxij10_RP91 = 2.95156019175951646E-003;
	if (fabs(dadxij10/dadxij10_RP91-1) > tol){throw ValueError();}

	double dadxij11 = d2phirdxidxj(tau, delta, z, 1, 1);
	double dadxij11_RP91 = 0.0;
	if (fabs(dadxij11-dadxij11_RP91) > tol){throw ValueError();}

	double d2adxn00 = d_ndphirdni_dxj__constdelta_tau_xi(tau,delta,z,0,0);
	double d2adxn00_RP91 = 1.4701046408163372;
	if (fabs(d2adxn00-d2adxn00_RP91) > tol){throw ValueError();}

	double d2adxn01 = d_ndphirdni_dxj__constdelta_tau_xi(tau,delta,z,0,1);
	double d2adxn01_RP91 = 1.4673267962070480;
	if (fabs(d2adxn01/d2adxn01_RP91-1) > tol){throw ValueError();}
	
	double d2adxn10 = d_ndphirdni_dxj__constdelta_tau_xi(tau,delta,z,1,0);
	double d2adxn10_RP91 = 1.5168952199179448;
	if (fabs(d2adxn10/d2adxn10_RP91-1) > tol){throw ValueError();}
	
	double d2adxn11 = d_ndphirdni_dxj__constdelta_tau_xi(tau,delta,z,1,1);
	double d2adxn11_RP91 = 1.4329117162994811;
	if (fabs(d2adxn11/d2adxn11_RP91-1) > tol){throw ValueError();}

	double d2addn0 = d_ndphirdni_dDelta(tau,delta,z,0);
	double d2addn0_RP91 = -2.3060567022876386;
	if (fabs(d2addn0/d2addn0_RP91-1) > tol){throw ValueError();}

	double d2addn1 = d_ndphirdni_dDelta(tau,delta,z,1);
	double d2addn1_RP91 = -1.9520671401909093;
	if (fabs(d2addn1/d2addn1_RP91-1) > tol){throw ValueError();}

	double d2adtn0 = d_ndphirdni_dTau(tau,delta,z,0);
	double d2adtn0_RP91 = -0.62562519155517959;
	if (fabs(d2adtn0/d2adtn0_RP91-1) > tol){throw ValueError();}
	
	double d2adtn1 = d_ndphirdni_dTau(tau,delta,z,1);
	double d2adtn1_RP91 = -0.43264230263327769;
	if (fabs(d2adtn1/d2adtn1_RP91-1) > tol){throw ValueError();}

	double d2ann00 = nd2nphirdnidnj__constT_V(tau, delta, z, 0, 0);
	double d2dnn00_RP91 = -0.38451299457845400;
	if (fabs(d2ann00/d2dnn00_RP91-1) > tol){throw ValueError();}

	double d2ann01 = nd2nphirdnidnj__constT_V(tau, delta, z, 0, 1);
	double d2dnn01_RP91 = -0.32103958765767326;
	if (fabs(d2ann01/d2dnn01_RP91-1) > tol){throw ValueError();}
	
	double d2ann10 = nd2nphirdnidnj__constT_V(tau, delta, z, 1, 0);
	double d2dnn10_RP91 = -0.32103958765767304;
	if (fabs(d2ann10/d2dnn10_RP91-1) > tol){throw ValueError();}
	
	double d2ann11 = nd2nphirdnidnj__constT_V(tau, delta, z, 1, 1);
	double d2dnn11_RP91 = -0.31715926219249480;
	if (fabs(d2ann11/d2dnn11_RP91-1) > tol){throw ValueError();}

	double dphidnj00 = ndln_fugacity_coefficient_dnj__constT_p(tau,delta,z,0,0);
	double dphidnj00_RP91 = -2.18661832047073457E-002;
	if (fabs(dphidnj00/dphidnj00_RP91-1) > tol){throw ValueError();}
	
	double dphidnj01 = ndln_fugacity_coefficient_dnj__constT_p(tau,delta,z,0,1);
	double dphidnj01_RP91 = 1.45774554698049341E-002;
	if (fabs(dphidnj01/dphidnj01_RP91-1) > tol){throw ValueError();}
	
	double dphidnj10 = ndln_fugacity_coefficient_dnj__constT_p(tau,delta,z,1,0);
	double dphidnj10_RP91 = 1.45774554698051562E-002;
	if (fabs(dphidnj10/dphidnj10_RP91-1) > tol){throw ValueError();}

	double dphidnj11 = ndln_fugacity_coefficient_dnj__constT_p(tau,delta,z,1,1);
	double dphidnj11_RP91 = -9.71830364653669676E-003;
	if (fabs(dphidnj11/dphidnj11_RP91-1) > tol){throw ValueError();}

}

void Mixture::check_MethaneEthane()
{
	std::vector<double> z(2, 0.5);

	double tol = 1e-10; // relative tolerance
	z[0] = 0.5;
	z[1] = 1-z[0];

	double Tr = pReducing->Tr(z);
	double Tr_RP91 = 250.57185990876351; // [K]
	if (fabs(Tr/Tr_RP91-1) > tol){throw ValueError();}

	double rhorbar = pReducing->rhorbar(z);
	double rhorbar_RP91 = 8205.7858837320694; //[mol/m^3]
	if (fabs(rhorbar/rhorbar_RP91-1) > tol){throw ValueError();}

	double T = 145; // [K]
	double rhobar = 4; // [mol/m^3]; to convert from mol/L, multiply by 1000
	
	double tau = Tr/T;
	double delta = rhobar/rhorbar;

	double RT = Rbar(z)*T;
	double p = rhobar*RT*(1+delta*dphir_dDelta(tau, delta,z));
	double p_RP91 = 4814.1558068139991;
	if (fabs(p/p_RP91-1) > tol){throw ValueError();}

	double dtdn0 = pReducing->ndTrdni__constnj(z,0);
	double dtdn0_RP91 = -56.914351037292874;
	if (fabs(dtdn0/dtdn0_RP91-1) > tol){throw ValueError();}

	double dtdn1 = pReducing->ndTrdni__constnj(z,1);
	double dtdn1_RP91 = 56.914351037292874;
	if (fabs(dtdn1/dtdn1_RP91-1) > tol){throw ValueError();}

	double drhodn0 = pReducing->ndrhorbardni__constnj(z,0);
	double drhodn0_RP91 = 1579.4307575322835;
	if (fabs(drhodn0/drhodn0_RP91-1) > tol){throw ValueError();}

	double drhodn1 = pReducing->ndrhorbardni__constnj(z,1);
	double drhodn1_RP91 = -1579.4307575322843;
	if (fabs(drhodn1/drhodn1_RP91-1) > tol){throw ValueError();}

	double ddrdxn00 = pReducing->d_ndrhorbardni_dxj__constxi(z,0,0);
	double ddrdxn00_RP91 = 10652.242638037194;
	if (fabs(ddrdxn00/ddrdxn00_RP91-1) > tol){throw ValueError();}

	double ddrdxn01 = pReducing->d_ndrhorbardni_dxj__constxi(z,0,1);
	double ddrdxn01_RP91 = 12694.316351697381;
	if (fabs(ddrdxn01/ddrdxn01_RP91-1) > tol){throw ValueError();}

	double ddrdxn10 = pReducing->d_ndrhorbardni_dxj__constxi(z,1,0);
	double ddrdxn10_RP91 = 19012.039381826526;
	if (fabs(ddrdxn10/ddrdxn10_RP91-1) > tol){throw ValueError();}

	double ddrdxn11 = pReducing->d_ndrhorbardni_dxj__constxi(z,1,1);
	double ddrdxn11_RP91 = 23287.688698295469;
	if (fabs(ddrdxn11/ddrdxn11_RP91-1) > tol){throw ValueError();}

	double dtrdxn00 = pReducing->d_ndTrdni_dxj__constxi(z,0,0);
	double dtrdxn00_RP91 = -506.39802892344619;
	if (fabs(dtrdxn00/dtrdxn00_RP91-1) > tol){throw ValueError();}

	double dtrdxn01 = pReducing->d_ndTrdni_dxj__constxi(z,0,1);
	double dtrdxn01_RP91 = -609.71811278619361;
	if (fabs(dtrdxn01/dtrdxn01_RP91-1) > tol){throw ValueError();}
	
	double dtrdxn10 = pReducing->d_ndTrdni_dxj__constxi(z,1,0);
	double dtrdxn10_RP91 = -382.06070863702217;
	if (fabs(dtrdxn10/dtrdxn10_RP91-1) > tol){throw ValueError();}
	
	double dtrdxn11 = pReducing->d_ndTrdni_dxj__constxi(z,1,1);
	double dtrdxn11_RP91 = -506.39802892344630;
	if (fabs(dtrdxn11/dtrdxn11_RP91-1) > tol){throw ValueError();}

	double dadxi0 = this->dphir_dxi(tau, delta, z, 0);
	double dadxi0_RP91 = -1.66733159546361936E-003;
	if (fabs(dadxi0/dadxi0_RP91-1) > tol){throw ValueError();}

	double dadxi1 = this->dphir_dxi(tau, delta, z, 1);
	double dadxi1_RP91 = -1.82138497681156902E-003;
	if (fabs(dadxi1/dadxi1_RP91-1) > tol){throw ValueError();}

	double daddx0 = this->d2phir_dxi_dDelta(tau,delta,z,0);
	double daddx0_RP91 = -3.4250061506157587;
	if (fabs(daddx0/daddx0_RP91-1) > tol){throw ValueError();}

	double dadtx0 = this->d2phir_dxi_dTau(tau,delta,z,0);
	double dadtx0_RP91 = -1.94597122635832131E-003;
	if (fabs(dadtx0/dadtx0_RP91-1) > tol){throw ValueError();}

	double daddx1 = this->d2phir_dxi_dDelta(tau,delta,z,1);
	double daddx1_RP91 = -3.7448162669516916;
	if (fabs(daddx1/daddx1_RP91-1) > tol){throw ValueError();}

	double dadtx1 = this->d2phir_dxi_dTau(tau,delta,z,1);
	double dadtx1_RP91 = -2.15974774776696646E-003;
	if (fabs(dadtx1/dadtx1_RP91-1) > tol){throw ValueError();}

	double dadn0 = this->ndphir_dni__constT_V_nj(tau, delta, z, 0);
	double dadn0_RP91 = -5.24342982584834368E-004;
	if (fabs(dadn0/dadn0_RP91-1) > tol){throw ValueError();}

	double dadn1 = this->ndphir_dni__constT_V_nj(tau, delta, z, 1);
	double dadn1_RP91 = -2.89676062124885614E-003;
	if (fabs(dadn1/dadn1_RP91-1) > tol){throw ValueError();}

	double dpdn0 = ndpdni__constT_V_nj(tau, delta, z, 0);
	double dpdn0_RP91 = 4811.6359520642318;
	if (fabs(dpdn0/dpdn0_RP91-1) > tol){throw ValueError();}

	double dpdn1 = ndpdni__constT_V_nj(tau, delta, z, 1);
	double dpdn1_RP91 = 4800.1319544625067;
	if (fabs(dpdn1/dpdn1_RP91-1) > tol){throw ValueError();}

	double vhat0 = partial_molar_volume(tau, delta, z, 0);
	double vhat0_RP91 = 0.25029921648425145;
	if (fabs(vhat0/vhat0_RP91-1) > tol){throw ValueError();}
	
	double vhat1 = partial_molar_volume(tau, delta, z, 1);
	double vhat1_RP91 = 0.24970078351574867;
	if (fabs(vhat1/vhat1_RP91-1) > tol){throw ValueError();}

	double d2adbn0 = d2nphir_dni_dT(tau, delta, z, 0);
	double d2adbn0_RP91 = 2.91400791470335643E-005;
	if (fabs(d2adbn0/d2adbn0_RP91-1) > tol){throw ValueError();}
	
	double d2adbn1 = d2nphir_dni_dT(tau, delta, z, 1);
	double d2adbn1_RP91 = 6.58714026367381391E-005;
	if (fabs(d2adbn1/d2adbn1_RP91-1) > tol){throw ValueError();}

	double dphidT0 = dln_fugacity_coefficient_dT__constp_n(tau,delta,z,0);
	double dphidT0_RP91 = 8.84377505112714235E-006;
	if (fabs(dphidT0/dphidT0_RP91-1) > tol){throw ValueError();}

	double dphidT1 = dln_fugacity_coefficient_dT__constp_n(tau,delta,z,1);
	double dphidT1_RP91 = 6.21123852185909847E-005;
	if (fabs(dphidT1/dphidT1_RP91-1) > tol){throw ValueError();}

	double dphidP0 = dln_fugacity_coefficient_dp__constT_n(tau,delta,z,0);
	double dphidP0_RP91 = -1.07128474189810419E-007;
	if (fabs(dphidP0/dphidP0_RP91-1) > tol){throw ValueError();}

	double dphidP1 = dln_fugacity_coefficient_dp__constT_n(tau,delta,z,1);
	double dphidP1_RP91 = -6.03505693277411881E-007;
	if (fabs(dphidP1/dphidP1_RP91-1) > tol){throw ValueError();}

	double dadxij00 = d2phirdxidxj(tau, delta, z, 0, 0);
	double dadxij00_RP91 = 0.0;
	if (fabs(dadxij00-dadxij00_RP91) > tol){throw ValueError();}

	double dadxij01 = d2phirdxidxj(tau, delta, z, 0, 1);
	double dadxij01_RP91 = -1.44900341276804124E-004;
	if (fabs(dadxij01/dadxij01_RP91-1) > tol){throw ValueError();}
	
	double dadxij10 = d2phirdxidxj(tau, delta, z, 1, 0);
	double dadxij10_RP91 = -1.44900341276804124E-004;
	if (fabs(dadxij10/dadxij10_RP91-1) > tol){throw ValueError();}

	double dadxij11 = d2phirdxidxj(tau, delta, z, 1, 1);
	double dadxij11_RP91 = 0.0;
	if (fabs(dadxij11-dadxij11_RP91) > tol){throw ValueError();}

	double d2adxn00 = d_ndphirdni_dxj__constdelta_tau_xi(tau,delta,z,0,0);
	double d2adxn00_RP91 = 9.52786141739760811E-003;
	if (fabs(d2adxn00-d2adxn00_RP91) > tol){throw ValueError();}

	double d2adxn01 = d_ndphirdni_dxj__constdelta_tau_xi(tau,delta,z,0,1);
	double d2adxn01_RP91 = 1.11090273394917772E-002;
	if (fabs(d2adxn01/d2adxn01_RP91-1) > tol){throw ValueError();}
	
	double d2adxn10 = d_ndphirdni_dxj__constdelta_tau_xi(tau,delta,z,1,0);
	double d2adxn10_RP91 = 8.82661064281056729E-003;
	if (fabs(d2adxn10/d2adxn10_RP91-1) > tol){throw ValueError();}
	
	double d2adxn11 = d_ndphirdni_dxj__constdelta_tau_xi(tau,delta,z,1,1);
	double d2adxn11_RP91 = 1.16784901260031035E-002;
	if (fabs(d2adxn11/d2adxn11_RP91-1) > tol){throw ValueError();}

	double d2addn0 = d_ndphirdni_dDelta(tau,delta,z,0);
	double d2addn0_RP91 = -1.0719438474166139;
	if (fabs(d2addn0/d2addn0_RP91-1) > tol){throw ValueError();}

	double d2addn1 = d_ndphirdni_dDelta(tau,delta,z,1);
	double d2addn1_RP91 = -5.9657336386754745;
	if (fabs(d2addn1/d2addn1_RP91-1) > tol){throw ValueError();}

	double d2adtn0 = d_ndphirdni_dTau(tau,delta,z,0);
	double d2adtn0_RP91 = -4.58046408525892678E-004;
	if (fabs(d2adtn0/d2adtn0_RP91-1) > tol){throw ValueError();}
	
	double d2adtn1 = d_ndphirdni_dTau(tau,delta,z,1);
	double d2adtn1_RP91 = -3.54010070086436396E-003;
	if (fabs(d2adtn1/d2adtn1_RP91-1) > tol){throw ValueError();}

	double d2ann00 = nd2nphirdnidnj__constT_V(tau, delta, z, 0, 0);
	double d2dnn00_RP91 = -1.55709211017489475E-003;
	if (fabs(d2ann00/d2dnn00_RP91-1) > tol){throw ValueError();}

	double d2ann01 = nd2nphirdnidnj__constT_V(tau, delta, z, 0, 1);
	double d2dnn01_RP91 = -2.90907297842576788E-003;
	if (fabs(d2ann01/d2dnn01_RP91-1) > tol){throw ValueError();}
	
	double d2ann10 = nd2nphirdnidnj__constT_V(tau, delta, z, 1, 0);
	double d2dnn10_RP91 = -2.90907297842577049E-003;
	if (fabs(d2ann10/d2dnn10_RP91-1) > tol){throw ValueError();}
	
	double d2ann11 = nd2nphirdnidnj__constT_V(tau, delta, z, 1, 1);
	double d2dnn11_RP91 = -6.32815473413224101E-003;
	if (fabs(d2ann11/d2dnn11_RP91-1) > tol){throw ValueError();}

	double dphidnj00 = ndln_fugacity_coefficient_dnj__constT_p(tau,delta,z,0,0);
	double dphidnj00_RP91 = -5.18202802448741728E-004;
	if (fabs(dphidnj00/dphidnj00_RP91-1) > tol){throw ValueError();}
	
	double dphidnj01 = ndln_fugacity_coefficient_dnj__constT_p(tau,delta,z,0,1);
	double dphidnj01_RP91 = 5.18202802448741728E-004;
	if (fabs(dphidnj01/dphidnj01_RP91-1) > tol){throw ValueError();}
	
	double dphidnj10 = ndln_fugacity_coefficient_dnj__constT_p(tau,delta,z,1,0);
	double dphidnj10_RP91 = 5.18202802448741728E-004;
	if (fabs(dphidnj10/dphidnj10_RP91-1) > tol){throw ValueError();}

	double dphidnj11 = ndln_fugacity_coefficient_dnj__constT_p(tau,delta,z,1,1);
	double dphidnj11_RP91 = -5.18202802448741728E-004;
	if (fabs(dphidnj11/dphidnj11_RP91-1) > tol){throw ValueError();}

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
double Mixture::fugacity(double tau, double delta, const std::vector<double> &x, int i)
{
	double Tr = pReducing->Tr(x);
	double rhorbar = pReducing->rhorbar(x);
	double T = Tr/tau, rhobar = rhorbar*delta, Rbar = 8.314472;

	double dnphir_dni = phir(tau,delta,x) + ndphir_dni__constT_V_nj(tau,delta,x,i);

	double f_i = x[i]*rhobar*Rbar*T*exp(dnphir_dni);
	return f_i;
}
double Mixture::ln_fugacity_coefficient(double tau, double delta, const std::vector<double> &x, int i)
{
	return phir(tau,delta,x) + ndphir_dni__constT_V_nj(tau,delta,x,i)-log(1+delta*dphir_dDelta(tau,delta,x));
}
double Mixture::dln_fugacity_coefficient_dT__constrho(double tau, double delta, const std::vector<double> &x, int i)
{
	double Tr = pReducing->Tr(x); // [K]
	double T = Tr/tau;
	double dtau_dT = -tau/T;
	return (dphir_dTau(tau,delta,x) + d_ndphirdni_dTau(tau,delta,x,i)-1/(1+delta*dphir_dDelta(tau,delta,x))*(delta*d2phir_dDelta_dTau(tau,delta,x)))*dtau_dT;
}
double Mixture::dnphir_dni__constT_V_nj(double tau, double delta, const std::vector<double> &x, int i)
{
	// GERG Equation 7.42
	return phir(tau,delta,x) + ndphir_dni__constT_V_nj(tau, delta, x, i);
}
double Mixture::d2nphir_dni_dT(double tau, double delta, const std::vector<double> &x, int i)
{
	double Tr = pReducing->Tr(x); // [K]
	double T = Tr/tau;
	return -tau/T*(dphir_dTau(tau,delta,x) + d_ndphirdni_dTau(tau,delta,x,i));
}
double Mixture::dln_fugacity_coefficient_dT__constp_n(double tau, double delta, const std::vector<double> &x, int i)
{
	double Tr = pReducing->Tr(x); // [K]
	double T = Tr/tau;
	return d2nphir_dni_dT(tau, delta, x, i) + 1/T-this->partial_molar_volume(tau,delta,x,i)/(Rbar(x)*T)*dpdT__constV_n(tau,delta,x,i);
}
double Mixture::partial_molar_volume(double tau, double delta, const std::vector<double> &x, int i)
{
	return -ndpdni__constT_V_nj(tau,delta,x,i)/ndpdV__constT_n(tau,delta,x,i);
}

double Mixture::dln_fugacity_coefficient_dp__constT_n(double tau, double delta, const std::vector<double> &x, int i)
{
	// GERG equation 7.30
	double Tr = pReducing->Tr(x); // [K]
	double T = Tr/tau;
	double rhorbar = pReducing->rhorbar(x);
	double rhobar = rhorbar*delta;
	double RT = Rbar(x)*T; // J/mol/K*K = J/mol = N*m/mol
	double p = rhobar*RT*(1.0+delta*dphir_dDelta(tau, delta, x)); // [Pa]
	double partial_molar_volume = this->partial_molar_volume(tau,delta,x,i); // [m^3/mol]
	double term1 = partial_molar_volume/RT; // m^3/mol/(N*m)*mol = m^2/N = 1/Pa
	double term2 = 1.0/p;
	return term1 - term2;
}



double Mixture::dln_fugacity_coefficient_dxj__constT_p_xi(double tau, double delta, const std::vector<double> &x, int i, int j)
{
	// Gernert 3.115
	double Tr = pReducing->Tr(x); //[K]
	double T = Tr/tau;
	double term1 = d2nphir_dni_dxj__constT_V(tau, delta, x, i, j);
	double term2 = - this->partial_molar_volume(tau,delta,x,i)/(Rbar(x)*T)*dpdxj__constT_V_xi(tau,delta,x,j);
	// partial molar volume is -dpdn/dpdV, so need to flip the sign here
	return d2nphir_dni_dxj__constT_V(tau, delta, x, i, j) - this->partial_molar_volume(tau,delta,x,i)/(Rbar(x)*T)*dpdxj__constT_V_xi(tau,delta,x,j);
}
double Mixture::dpdxj__constT_V_xi(double tau, double delta, const std::vector<double> &x, int j)
{
	// Gernert 3.130
	double rhorbar = pReducing->rhorbar(x);
	double rhobar = rhorbar*delta;
	double Tr = pReducing->Tr(x); //[K]
	double T = Tr/tau;

	return rhobar*Rbar(x)*T*(ddelta_dxj__constT_V_xi(tau,delta,x,j)*dphir_dDelta(tau, delta, x)+delta*d_dphirddelta_dxj__constT_V_xi(tau,delta,x,j));
}

double Mixture::d_dphirddelta_dxj__constT_V_xi(double tau, double delta, const std::vector<double> &x, int j)
{
	// Gernert Equation 3.134 (Catch test provided)
	return d2phir_dDelta2(tau, delta, x)*ddelta_dxj__constT_V_xi(tau, delta, x, j)
		 + d2phir_dDelta_dTau(tau, delta, x)*dtau_dxj__constT_V_xi(tau, delta, x, j)
		 + d2phir_dxi_dDelta(tau, delta, x, j);
}
double Mixture::d2nphir_dni_dxj__constT_V(double tau, double delta, const std::vector<double> &x, int i, int j)
{
	// Gernert 3.117
	return d_ndphirdni_dxj__constT_V_xi(tau, delta, x, i, j) + dphir_dxj__constT_V_xi(tau, delta, x, j);
}
double Mixture::dphir_dxj__constT_V_xi(double tau, double delta, const std::vector<double> &x, int j)
{
	//Gernert 3.119 (Catch test provided)
	return dphir_dDelta(tau, delta, x)*ddelta_dxj__constT_V_xi(tau, delta, x, j)+dphir_dTau(tau, delta, x)*dtau_dxj__constT_V_xi(tau, delta, x, j)+dphir_dxi(tau, delta, x, j);
}
double Mixture::d_ndphirdni_dxj__constT_V_xi(double tau, double delta, const std::vector<double> &x, int i, int j)
{
	// Gernert 3.118
	return d_ndphirdni_dxj__constdelta_tau_xi(tau,delta,x,i,j)
		  + ddelta_dxj__constT_V_xi(tau,delta,x,j)*d_ndphirdni_dDelta(tau,delta,x,i) 
		  + dtau_dxj__constT_V_xi(tau,delta,x,j)*d_ndphirdni_dTau(tau,delta,x,i);
}
double Mixture::ddelta_dxj__constT_V_xi(double tau, double delta, const std::vector<double> &x, int j)
{
	// Gernert 3.121 (Catch test provided)
	double rhorbar = pReducing->rhorbar(x);
	return -delta/rhorbar*pReducing->drhorbardxi__constxj(x,j);
}
double Mixture::dtau_dxj__constT_V_xi(double tau, double delta, const std::vector<double> &x, int j)
{
	// Gernert 3.122 (Catch test provided)
	double Tr = pReducing->Tr(x); //[K]
	double T = Tr/tau;
	return 1/T*pReducing->dTrdxi__constxj(x,j);
}



double Mixture::dpdT__constV_n(double tau, double delta, const std::vector<double> &x, int i)
{
	double rhorbar = pReducing->rhorbar(x);
	double rhobar = rhorbar*delta;
	return rhobar*Rbar(x)*(1+delta*dphir_dDelta(tau, delta, x)-delta*tau*d2phir_dDelta_dTau(tau, delta, x));
}
double Mixture::ndpdV__constT_n(double tau, double delta, const std::vector<double> &x, int i)
{
	double rhorbar = pReducing->rhorbar(x);
	double rhobar = rhorbar*delta;
	double Tr = pReducing->Tr(x); // [K]
	double T = Tr/tau;
	return -rhobar*rhobar*Rbar(x)*T*(1+2*delta*dphir_dDelta(tau,delta,x)+delta*delta*d2phir_dDelta2(tau, delta, x));
}
double Mixture::ndpdni__constT_V_nj(double tau, double delta, const std::vector<double> &x, int i)
{
	// Eqn 7.64 and 7.63
	double rhorbar = pReducing->rhorbar(x);
	double rhobar = rhorbar*delta;
	double Tr = pReducing->Tr(x); // [K]
	double T = Tr/tau;
	double ndrhorbar_dni__constnj = pReducing->ndrhorbardni__constnj(x,i);
	double ndTr_dni__constnj = pReducing->ndTrdni__constnj(x,i);
	double summer = 0;
	for (unsigned int k = 0; k < x.size(); k++)
	{
		summer += x[k]*d2phir_dxi_dDelta(tau,delta,x,k);
	}
	double nd2phir_dni_dDelta = delta*d2phir_dDelta2(tau,delta,x)*(1-1/rhorbar*ndrhorbar_dni__constnj)+tau*d2phir_dDelta_dTau(tau,delta,x)/Tr*ndTr_dni__constnj+d2phir_dxi_dDelta(tau,delta,x,i)-summer;
	return rhobar*Rbar(x)*T*(1+delta*dphir_dDelta(tau,delta,x)*(2-1/rhorbar*ndrhorbar_dni__constnj)+delta*nd2phir_dni_dDelta);
}

double Mixture::ndphir_dni__constT_V_nj(double tau, double delta, const std::vector<double> &x, int i)
{
	double Tr = pReducing->Tr(x);
	double rhorbar = pReducing->rhorbar(x);

	double term1 = delta*dphir_dDelta(tau,delta,x)*(1-1/rhorbar*pReducing->ndrhorbardni__constnj(x,i));
	double term2 = tau*dphir_dTau(tau,delta,x)*(1/Tr)*pReducing->ndTrdni__constnj(x,i);

	double s = 0;
	for (unsigned int k = 0; k < x.size(); k++)
	{
		s += x[k]*dphir_dxi(tau,delta,x,k);
	}
	double term3 = dphir_dxi(tau,delta,x,i);
	return term1 + term2 + term3 - s;
}
double Mixture::ndln_fugacity_coefficient_dnj__constT_p(double tau, double delta, const std::vector<double> &x, int i, int j)
{
	double Tr = pReducing->Tr(x);
	double T = Tr/tau;
	return nd2nphirdnidnj__constT_V(tau, delta, x, j, i) + 1 - partial_molar_volume(tau,delta,x,i)/(Rbar(x)*T)*ndpdni__constT_V_nj(tau,delta,x,j);
}
double Mixture::nddeltadni__constT_V_nj(double tau, double delta, const std::vector<double> &x, int i)
{
	double rhorbar = pReducing->rhorbar(x);
	return delta-delta/rhorbar*pReducing->ndrhorbardni__constnj(x,i);
}
double Mixture::ndtaudni__constT_V_nj(double tau, double delta, const std::vector<double> &x, int i)
{
	double Tr = pReducing->Tr(x);
	return tau/Tr*pReducing->ndTrdni__constnj(x,i);
}
double Mixture::d_ndphirdni_dxj__constdelta_tau_xi(double tau, double delta, const std::vector<double> &x, int i, int j)
{
	double rhorbar = pReducing->rhorbar(x);
	double Tr = pReducing->Tr(x); // [K]

	double line1 = delta*d2phir_dxi_dDelta(tau,delta,x,j)*(1-1/rhorbar*pReducing->ndrhorbardni__constnj(x,i));
	double line2 = -delta*dphir_dDelta(tau,delta,x)*(1/rhorbar)*(pReducing->d_ndrhorbardni_dxj__constxi(x,i,j)-1/rhorbar*pReducing->drhorbardxi__constxj(x,j)*pReducing->ndrhorbardni__constnj(x,i));
	double line3 = tau*d2phir_dxi_dTau(tau,delta,x,j)*(1/Tr)*pReducing->ndTrdni__constnj(x,i);
	double line4 = tau*dphir_dTau(tau,delta,x)*(1/Tr)*(pReducing->d_ndTrdni_dxj__constxi(x,i,j)-1/Tr*pReducing->dTrdxi__constxj(x,j)*pReducing->ndTrdni__constnj(x,i));
	double s = 0;
	for (unsigned int m = 0; m < x.size(); m++)
	{
		s += x[m]*d2phirdxidxj(tau,delta,x,j,m);
	}
	double line5 = d2phirdxidxj(tau,delta,x,i,j)-dphir_dxi(tau,delta,x,j)-s;
	return line1+line2+line3+line4+line5;
}
double Mixture::nd2nphirdnidnj__constT_V(double tau, double delta, const std::vector<double> &x, int i, int j)
{	
	double line0 = ndphir_dni__constT_V_nj(tau, delta, x, j); // First term from 7.46
	double line1 = this->d_ndphirdni_dDelta(tau, delta, x, i)*this->nddeltadni__constT_V_nj(tau,delta,x,j);
	double line2 = this->d_ndphirdni_dTau(tau, delta, x, i)*this->ndtaudni__constT_V_nj(tau,delta,x,j);
	double summer = 0;
	for (unsigned int k = 0; k < x.size(); k++)
	{
		summer += x[k]*this->d_ndphirdni_dxj__constdelta_tau_xi(tau, delta, x, i, k);
	}
	double line3 = this->d_ndphirdni_dxj__constdelta_tau_xi(tau, delta, x, i, j)-summer;
	return line0 + line1 + line2 + line3;
}
double Mixture::d_ndphirdni_dDelta(double tau, double delta, const std::vector<double> &x, int i)
{
	double Tr = pReducing->Tr(x);
	double rhorbar = pReducing->rhorbar(x);

	// The first line
	double term1 = (delta*d2phir_dDelta2(tau,delta,x)+dphir_dDelta(tau,delta,x))*(1-1/rhorbar*pReducing->ndrhorbardni__constnj(x,i));

	// The second line
	double term2 = tau*d2phir_dDelta_dTau(tau,delta,x)*(1/Tr)*pReducing->ndTrdni__constnj(x,i);

	// The third line
	double term3 = d2phir_dxi_dDelta(tau,delta,x,i);
	for (unsigned int k = 0; k < x.size(); k++)
	{
		term3 -= x[k]*d2phir_dxi_dDelta(tau,delta,x,k);
	}
	return term1 + term2 + term3;
}

double Mixture::d_ndphirdni_dTau(double tau, double delta, const std::vector<double> &x, int i)
{
	double Tr = pReducing->Tr(x);
	double rhorbar = pReducing->rhorbar(x);

	// The first line
	double term1 = delta*d2phir_dDelta_dTau(tau,delta,x)*(1-1/rhorbar*pReducing->ndrhorbardni__constnj(x,i));

	// The second line
	double term2 = (tau*d2phir_dTau2(tau,delta,x)+dphir_dTau(tau,delta,x))*(1/Tr)*pReducing->ndTrdni__constnj(x,i);

	// The third line
	double term3 = d2phir_dxi_dTau(tau,delta,x,i);
	for (unsigned int k = 0; k < x.size(); k++)
	{
		term3 -= x[k]*d2phir_dxi_dTau(tau,delta,x,k);
	}
	return term1 + term2 + term3;
}

double Mixture::phir(double tau, double delta, const std::vector<double> &x)
{
	return pResidualIdealMix->phir(tau,delta,x) + pExcess->phir(tau,delta,x);
}
double Mixture::dphir_dxi(double tau, double delta, const std::vector<double> &x, int i)
{	
	return pFluids[i]->phir(tau,delta) + pExcess->dphir_dxi(tau,delta,x,i);
}
double Mixture::d2phirdxidxj(double tau, double delta, const std::vector<double> &x, int i, int j)
{	
	return 0                           + pExcess->d2phirdxidxj(tau,delta,x,i,j);
}
double Mixture::d2phir_dxi_dTau(double tau, double delta, const std::vector<double> &x, int i)
{	
	return pFluids[i]->dphir_dTau(tau,delta) + pExcess->d2phir_dxi_dTau(tau,delta,x,i);
}
double Mixture::d2phir_dxi_dDelta(double tau, double delta, const std::vector<double> &x, int i)
{	
	return pFluids[i]->dphir_dDelta(tau,delta) + pExcess->d2phir_dxi_dDelta(tau,delta,x,i);
}
double Mixture::dphir_dDelta(double tau, double delta, const std::vector<double> &x)
{
	return pResidualIdealMix->dphir_dDelta(tau,delta,x) + pExcess->dphir_dDelta(tau,delta,x);
}
double Mixture::d2phir_dDelta2(double tau, double delta, const std::vector<double> &x)
{
	return pResidualIdealMix->d2phir_dDelta2(tau,delta,x) + pExcess->d2phir_dDelta2(tau,delta,x);
}
double Mixture::d2phir_dDelta_dTau(double tau, double delta, const std::vector<double> &x)
{
	return pResidualIdealMix->d2phir_dDelta_dTau(tau,delta,x) + pExcess->d2phir_dDelta_dTau(tau,delta,x);
}
double Mixture::d2phir_dTau2(double tau, double delta, const std::vector<double> &x)
{
	return pResidualIdealMix->d2phir_dTau2(tau,delta,x) + pExcess->d2phir_dTau2(tau,delta,x);
}
double Mixture::dphir_dTau(double tau, double delta, const std::vector<double> &x)
{
	return pResidualIdealMix->dphir_dTau(tau,delta,x) + pExcess->dphir_dTau(tau,delta,x);
}

/// A wrapper function around the Rachford-Rice residual
class gRR_resid : public FuncWrapper1D
{
public:
	const std::vector<double> * z;
	const std::vector<double> * lnK;
	Mixture *Mix;

	gRR_resid(Mixture *Mix, std::vector<double> const &z, std::vector<double> const& lnK){ this->z = &z; this->lnK = &lnK; this->Mix = Mix; };
	double call(double beta){return Mix->g_RachfordRice(*z, *lnK, beta); };
	double deriv(double beta){return Mix->dgdbeta_RachfordRice(*z, *lnK, beta); };
};

/// A wrapper function around the density(T,p,x) residual
class rho_Tpz_resid : public FuncWrapper1D
{
protected:
	double T, p, Rbar, tau, Tr, rhorbar,dphir_dDelta;
public:
	const std::vector<double> *x;
	Mixture *Mix;

	rho_Tpz_resid(Mixture *Mix, double T, double p, const std::vector<double> &x){ 
		this->x = &x; this->T = T; this->p = p; this->Mix = Mix;
		
		Tr = Mix->pReducing->Tr(x);
		rhorbar = Mix->pReducing->rhorbar(x);
		tau = Tr/T;
		Rbar = Mix->Rbar(x); // J/mol/K
	};
	double call(double rhobar){	
		double delta = rhobar/rhorbar;
		dphir_dDelta = Mix->dphir_dDelta(tau, delta, *x);
		double resid = Rbar*rhobar*T*(1 + delta*dphir_dDelta)-p;
		return resid;
	}
	double deriv(double rhobar){
		double delta = rhobar/rhorbar;
		double val = Rbar*T*(1 + 2*delta*Mix->dphir_dDelta(tau, delta, *x)+delta*delta*Mix->d2phir_dDelta2(tau, delta, *x));
		return val;
	}
};
double Mixture::rhobar_Tpz(double T, double p, const std::vector<double> &x, double rhobar0)
{
	rho_Tpz_resid Resid(this,T,p,x);
	double rhobar;
	std::string errstr;
	//FILE *fp;
	//fp = fopen("rhobar_resid.csv","w");
	//for (double rhobar = 0; rhobar < 20000; rhobar += 1)
	//{
	//	fprintf(fp,"%g, %g\n",rhobar, Resid.call(rhobar));
	//}
	//fclose(fp);

	rhobar = Newton(&Resid, rhobar0, 1e-16, 100, &errstr);
	if (!ValidNumber(rhobar)) 
	{ 
		rhobar = Secant(&Resid, rhobar0, 1.001*rhobar0, 1e-16, 100, &errstr);
		if (!ValidNumber(rhobar)) 
		{
			throw ValueError("Density solver failed"); 
		}
	}
	return rhobar;
}


/*! A wrapper function around the residual to find the initial guess for the bubble point temperature
\f[
r = \sum_i \frac{z_i(K_i-1)}{1-beta+beta*K_i}
\f]
*/
class WilsonK_resid : public FuncWrapper1D
{
public:
	double p, beta;
	const std::vector<double> *z;
	std::vector<double> K;
	Mixture *Mix;

	WilsonK_resid(Mixture *Mix, double beta, double p, std::vector<double> const& z){ 
		this->z = &z; this->p = p; this->Mix = Mix, this->beta = beta; 
		K = std::vector<double>(z.size(),_HUGE);
	};
	double call(double T){
		double summer = 0;
		for (unsigned int i = 0; i< (*z).size(); i++) {
			K[i] = exp(Mix->Wilson_lnK_factor(T,p,i));
			summer += (*z)[i]*(K[i]-1)/(1-beta+beta*K[i]);
		}
		return summer;
	};
};
double Mixture::saturation_p_preconditioner(double p, const std::vector<double> &z)
{
	double pseudo_ptriple = 0;
	double pseudo_pcrit = 0;
	double pseudo_Ttriple = 0;
	double pseudo_Tcrit = 0;

	// Check if a pure fluid, if so, mole fractions and saturation temperatures are equal to bulk composition
	for (unsigned int i = 0; i < N; i++)
	{
		pseudo_ptriple += pFluids[i]->params.ptriple*z[i];
		pseudo_pcrit += pFluids[i]->crit.p.Pa*z[i];
		pseudo_Ttriple += pFluids[i]->params.Ttriple*z[i];
		pseudo_Tcrit += pFluids[i]->crit.T*z[i];
	}

	// Preliminary guess based on interpolation of log(p) v. T
	return (pseudo_Tcrit-pseudo_Ttriple)/(pseudo_pcrit/pseudo_ptriple)*(p/pseudo_ptriple)+pseudo_Ttriple;
}
double Mixture::saturation_p_Wilson(double beta, double p, const std::vector<double> &z, double T_guess, std::vector<double> &K)
{
	double T;

	std::string errstr;

	// Find first guess for T using Wilson K-factors
	WilsonK_resid Resid(this, beta, p, z);
	T = Secant(&Resid, T_guess, 0.001, 1e-10, 100, &errstr);
	// Get the K factors from the residual wrapper
	K = Resid.K;
	
	if (!ValidNumber(T)){throw ValueError("saturation_p_Wilson failed to get good T");}
	return T;
}
double Mixture::saturation_p(double beta, double p, std::vector<double> const& z, std::vector<double> &x, std::vector<double> &y)
{
	double T;
	unsigned int N = z.size();
	std::vector<double> K(N), ln_phi_liq(N), ln_phi_vap(N);

	// Check if a pure fluid, if so, mole fractions and saturation temperatures are equal to bulk composition
	for (unsigned int i = 0; i < N; i++)
	{
		if (fabs(z[i]-1) < 10*DBL_EPSILON)
		{
			// Mole fractions are equal to bulk composition
			x = z;
			y = z;
			// Get temperature from pure-fluid saturation call
			double TL=0, TV=0, rhoL=0, rhoV=0;
			pFluids[i]->saturation_p(p,false, TL, TV, rhoL, rhoV);
			return TL; // Also equal to TV for pure fluid
		}
	}

	// Preliminary guess based on interpolation of log(p) v. T
	double T_guess = saturation_p_preconditioner(p, z);

	// Initialize the call using Wilson to get K-factors and temperature
	T = saturation_p_Wilson(beta, p, z, T_guess, K);

	// Call the successive substitution routine with the guess values provided above from Wilson
	// It will call the Newton-Raphson routine after a few steps of successive-substitution
	SS.useNR = true;
	T = SS.call(beta, T, p, z, K);

	x = SS.x;
	y = SS.y;

	return T;
}
void Mixture::TpzFlash(double T, double p, const std::vector<double> &z, double &rhobar, std::vector<double> &x, std::vector<double> &y)
{
	unsigned int N = z.size();
	double beta, change = 0;
	std::vector<double> lnK(N);
	
	x.resize(N);
	y.resize(N);

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
		rhobar = rhobar_Tpz(T,p,z,rhobar_pengrobinson(T,p,z,PR_SATL));
		return;
	}
	else
	{
		double g_RR_1 = g_RachfordRice(z, lnK, 1);
		if (g_RR_1 > 0)
		{
			// Superheated vapor - done
			rhobar = rhobar_Tpz(T,p,z,rhobar_pengrobinson(T,p,z,PR_SATV));
			return;
		}
	}
	// TODO: How can you be sure that you aren't in the two-phase region? Safety factor needed?
	// TODO: Calculate the dewpoint density
	do
	{
		// Now find the value of beta that satisfies Rachford-Rice using Brent's method

		gRR_resid Resid(this,z,lnK);
		std::string errstr;
		beta = Newton(&Resid,0,1e-10,300,&errstr);

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

		// Reducing parameters for each phase
		double tau_liq = pReducing->Tr(x)/T;
		double tau_vap = pReducing->Tr(y)/T;
		
		double rhobar_liq = rhobar_Tpz(T, p, x, rhobar_pengrobinson(T,p,x,PR_SATL));
		double rhobar_vap = rhobar_Tpz(T, p, y, rhobar_pengrobinson(T,p,y,PR_SATV));
		double rhorbar_liq = pReducing->rhorbar(x);
		double rhorbar_vap = pReducing->rhorbar(y);
		double delta_liq = rhobar_liq/rhorbar_liq; 
		double delta_vap = rhobar_vap/rhorbar_vap; 
		
		// Evaluate fugacity coefficients in liquid and vapor
		for (unsigned int i = 0; i < N; i++)
		{
			double ln_phi_liq = ln_fugacity_coefficient(tau_liq, delta_vap, x, i);
			double ln_phi_vap = ln_fugacity_coefficient(tau_vap, delta_liq, x, i);
			
			double lnKold = lnK[i];
			// Recalculate the K-factor (log(exp(ln_phi_liq)/exp(ln_phi_vap)))
			lnK[i] = ln_phi_liq - ln_phi_vap;

			change = lnK[i] - lnKold;
		}
	}
	while( fabs(change) > 1e-7);

	return;
}
void Mixture::x_and_y_from_K(double beta, const std::vector<double> &K, const std::vector<double> &z, std::vector<double> &x, std::vector<double> &y)
{
	for (unsigned int i=0; i < N; i++)
	{
		double denominator = (1-beta+beta*K[i]); // Common denominator
		x[i] = z[i]/denominator;
		y[i] = K[i]*z[i]/denominator;
	}
	normalize_vector(x);
	normalize_vector(y);
}
double Mixture::rhobar_pengrobinson(double T, double p, const std::vector<double> &x, int solution)
{ 
	double A  = 0, B = 0, a = 0, b = 0, m_i, m_j, a_i, a_j, b_i, R = Rbar(x), k_ij;
	
	for (unsigned int i = 0; i < N; i++)
	{
		
		b_i = 0.077796074*(R*pFluids[i]->reduce.T)/(pFluids[i]->reduce.p.Pa);
		b += x[i]*b_i;
		
		double omega_i = pFluids[i]->params.accentricfactor;
		//m_i = 0.480 + 1.574*omega_i-0.176*pow(omega_i,(int)2);
		m_i = 0.37464 + 1.54226*omega_i-0.26992*pow(omega_i,(int)2);

		for (unsigned int j = 0; j < N; j++)
		{
			double omega_j = pFluids[j]->params.accentricfactor;
			m_j = 0.37464 + 1.54226*omega_j-0.26992*pow(omega_j,(int)2);
			//m_j = 0.480 + 1.574*omega_j - 0.176*pow(omega_j,(int)2);
			a_i = 0.45724*pow(R*pFluids[i]->reduce.T,2)/pFluids[i]->reduce.p.Pa*pow(1+m_i*(1-sqrt(T/pFluids[i]->reduce.T)),2);
			
			
			a_j = 0.45724*pow(R*pFluids[j]->reduce.T,2)/pFluids[j]->reduce.p.Pa*pow(1+m_j*(1-sqrt(T/pFluids[j]->reduce.T)),2);
			if (i!=j)
			{
				k_ij = -0.3; // Hard coded for Ethanol-Water
			}
			else
			{
				k_ij = 0;
			}
			a += x[i]*x[j]*sqrt(a_i*a_j)*(1-k_ij);
		}
	}
	A = a*p/(R*R*T*T);
	B = b*p/(R*T);
	
	double Z0, Z1, Z2;
	solve_cubic(1, -1+B, A-3*B*B-2*B, -A*B+B*B+B*B*B, &Z0, &Z1, &Z2);

	if (solution == PR_SATL)
	{
		return p/(min3(Z0,Z1,Z2)*R*T);
	}
	else if (solution == PR_SATV)
	{
		return p/(max3(Z0,Z1,Z2)*R*T);
	}
	else 
	{
		throw ValueError("Solution should be one of PR_SATL, PR_SATV");
	}
}

double Mixture::g_RachfordRice(const std::vector<double> &z, const std::vector<double> &lnK, double beta)
{
	// g function from Rachford-Rice
	double summer = 0;
	for (unsigned int i = 0; i < z.size(); i++)
	{
		double Ki = exp(lnK[i]);
		summer += z[i]*(Ki-1)/(1-beta+beta*Ki);
	}
	return summer;
}
double Mixture::dgdbeta_RachfordRice(const std::vector<double> &z, const std::vector<double> &lnK, double beta)
{
	// derivative of g function from Rachford-Rice with respect to beta
	double summer = 0;
	for (unsigned int i = 0; i < z.size(); i++)
	{
		double Ki = exp(lnK[i]);
		summer += -z[i]*pow((Ki-1)/(1-beta+beta*Ki),2);
	}
	return summer;
}

double GERG2008ReducingFunction::Tr(const std::vector<double> &x)
{
	double Tr = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		double xi = x[i], Tci = pFluids[i]->reduce.T;
		Tr += xi*xi*Tci;
		
		// The last term is only used for the pure component, as it is sum_{i=1}^{N-1}sum_{j=1}^{N}
		if (i==N-1){ break; }

		for (unsigned int j = i+1; j < N; j++)
		{
			Tr += c_Y_ij(i, j, &beta_T, &gamma_T, &T_c)*f_Y_ij(x, i, j, &beta_T);
		}
	}
	return Tr;
}
double GERG2008ReducingFunction::dTrdxi__constxj(const std::vector<double> &x, int i)
{
	// See Table B9 from Kunz Wagner 2012 (GERG 2008)
	double xi = x[i];
	double dTr_dxi = 2*xi*pFluids[i]->reduce.T;
	for (int k = 0; k < i; k++)
	{
		dTr_dxi += c_Y_ji(k,i,&beta_T,&gamma_T,&T_c)*dfYkidxi__constxk(x,k,i,&beta_T);
	}
	for (unsigned int k = i+1; k < N; k++)
	{
		dTr_dxi += c_Y_ij(i,k,&beta_T,&gamma_T,&T_c)*dfYikdxi__constxk(x,i,k,&beta_T);
	}
	return dTr_dxi;
}
double GERG2008ReducingFunction::d2Trdxi2__constxj(const std::vector<double> &x, int i)
{
	// See Table B9 from Kunz Wagner 2012 (GERG 2008)
	double d2Tr_dxi2 = 2*pFluids[i]->reduce.T;
	for (int k = 0; k < i; k++)
	{
		d2Tr_dxi2 += c_Y_ij(k,i,&beta_T,&gamma_T,&T_c)*d2fYkidxi2__constxk(x,k,i,&beta_T);
	}
	for (unsigned int k = i+1; k < N; k++)
	{
		d2Tr_dxi2 += c_Y_ij(i,k,&beta_T,&gamma_T,&T_c)*d2fYikdxi2__constxk(x,i,k,&beta_T);
	}
	return d2Tr_dxi2;
}
double GERG2008ReducingFunction::d2Trdxidxj(const std::vector<double> &x, int i, int j)
{
	if (i == j)
	{
		return d2Trdxi2__constxj(x,i);
	}
	else
	{
		// See Table B9 from Kunz Wagner 2012 (GERG 2008)
		return c_Y_ij(i,j,&beta_T,&gamma_T,&T_c)*d2fYijdxidxj(x,i,j,&beta_T);
	}
}
double GERG2008ReducingFunction::dvrbardxi__constxj(const std::vector<double> &x, int i)
{
	// See Table B9 from Kunz Wagner 2012 (GERG 2008)
	double xi = x[i];
	double dvrbar_dxi = 2*xi/pFluids[i]->reduce.rhobar;

	for (int k = 0; k < i; k++)
	{
		dvrbar_dxi += c_Y_ij(k, i, &beta_v, &gamma_v, &v_c)*dfYkidxi__constxk(x, k, i, &beta_v);
	}
	for (unsigned int k = i+1; k < N; k++)
	{
		dvrbar_dxi += c_Y_ij(i, k, &beta_v, &gamma_v, &v_c)*dfYikdxi__constxk(x, i, k, &beta_v);	
	}
	return dvrbar_dxi;
}
double GERG2008ReducingFunction::d2vrbardxidxj(const std::vector<double> &x, int i, int j)
{
	if (i == j)
	{
		return d2vrbardxi2__constxj(x, i);
	}
	else
	{
		return c_Y_ij(i, j, &beta_v, &gamma_v, &v_c)*d2fYijdxidxj(x, i, j, &beta_v);
	}
}
double GERG2008ReducingFunction::drhorbardxi__constxj(const std::vector<double> &x, int i)
{
	return -pow(rhorbar(x),2)*dvrbardxi__constxj(x,i);
}
double GERG2008ReducingFunction::d2vrbardxi2__constxj(const std::vector<double> &x, int i)
{
	// See Table B9 from Kunz Wagner 2012 (GERG 2008)
	double d2vrbardxi2 = 2/pFluids[i]->reduce.rhobar;

	for (int k = 0; k < i; k++)
	{
		d2vrbardxi2 += c_Y_ij(k, i, &beta_v, &gamma_v, &v_c)*d2fYkidxi2__constxk(x, k, i, &beta_v);
	}
	for (unsigned int k = i+1; k < N; k++)
	{
		d2vrbardxi2 += c_Y_ij(i, k, &beta_v, &gamma_v, &v_c)*d2fYikdxi2__constxk(x, i, k, &beta_v);	
	}
	return d2vrbardxi2;
}
double GERG2008ReducingFunction::d2rhorbardxi2__constxj(const std::vector<double> &x, int i)
{
	double rhor = this->rhorbar(x);
	double dvrbardxi = this->dvrbardxi__constxj(x,i);
	return 2*pow(rhor,(int)3)*pow(dvrbardxi,(int)2)-pow(rhor,(int)2)*this->d2vrbardxi2__constxj(x,i);
}
double GERG2008ReducingFunction::d2rhorbardxidxj(const std::vector<double> &x, int i, int j)
{
	double rhor = this->rhorbar(x);
	double dvrbardxi = this->dvrbardxi__constxj(x,i);
	double dvrbardxj = this->dvrbardxi__constxj(x,j);
	return 2*pow(rhor,(int)3)*dvrbardxi*dvrbardxj-pow(rhor,(int)2)*this->d2vrbardxidxj(x,i,j);
}

double GERG2008ReducingFunction::rhorbar(const std::vector<double> &x)
{
	double vrbar = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		double xi = x[i];
		vrbar += xi*xi/pFluids[i]->reduce.rhobar;

		if (i == N-1){ break; }

		for (unsigned int j = i+1; j < N; j++)
		{	
			vrbar += c_Y_ij(i, j, &beta_v, &gamma_v, &v_c)*f_Y_ij(x,i,j,&beta_v);
		}
	}
	return 1/vrbar;
}
double GERG2008ReducingFunction::dfYkidxi__constxk(const std::vector<double> &x, int k, int i, std::vector< std::vector< double> > * beta)
{
	double xk = x[k], xi = x[i], beta_Y = (*beta)[k][i];
	return xk*(xk+xi)/(beta_Y*beta_Y*xk+xi)+xk*xi/(beta_Y*beta_Y*xk+xi)*(1-(xk+xi)/(beta_Y*beta_Y*xk+xi));
}
double GERG2008ReducingFunction::dfYikdxi__constxk(const std::vector<double> &x, int i, int k, std::vector< std::vector< double> > * beta)
{
	double xk = x[k], xi = x[i], beta_Y = (*beta)[i][k];
	return xk*(xi+xk)/(beta_Y*beta_Y*xi+xk)+xi*xk/(beta_Y*beta_Y*xi+xk)*(1-beta_Y*beta_Y*(xi+xk)/(beta_Y*beta_Y*xi+xk));
}
double GERG2008ReducingFunction::c_Y_ij(int i, int j, std::vector< std::vector< double> > * beta, std::vector< std::vector< double> > *gamma, std::vector< std::vector< double> > *Y_c)
{
	return 2*((*beta)[i][j])*((*gamma)[i][j])*((*Y_c)[i][j]);
}
double GERG2008ReducingFunction::c_Y_ji(int j, int i, std::vector< std::vector< double> > * beta, std::vector< std::vector< double> > *gamma, std::vector< std::vector< double> > *Y_c)
{
	return 2/((*beta)[i][j])*((*gamma)[i][j])*((*Y_c)[i][j]);
}
double GERG2008ReducingFunction::f_Y_ij(const std::vector<double> &x, int i, int j, std::vector< std::vector< double> > * beta)
{
	double xi = x[i], xj = x[j], beta_Y = (*beta)[i][j];
	return xi*xj*(xi+xj)/(beta_Y*beta_Y*xi+xj);
}
double GERG2008ReducingFunction::d2fYikdxi2__constxk(const std::vector<double> &x, int i, int k, std::vector< std::vector< double> > * beta)
{
	double xi = x[i], xk = x[k], beta_Y = (*beta)[i][k];
	return 1/(beta_Y*beta_Y*xi+xk)*(1-beta_Y*beta_Y*(xi+xk)/(beta_Y*beta_Y*xi+xk))*(2*xk-xi*xk*2*beta_Y*beta_Y/(beta_Y*beta_Y*xi+xk));
}
double GERG2008ReducingFunction::d2fYkidxi2__constxk(const std::vector<double> &x, int k, int i, std::vector< std::vector< double> > * beta)
{
	double xi = x[i], xk = x[k], beta_Y = (*beta)[k][i];
	return 1/(beta_Y*beta_Y*xk+xi)*(1-(xk+xi)/(beta_Y*beta_Y*xk+xi))*(2*xk-xk*xi*2/(beta_Y*beta_Y*xk+xi));
}
double GERG2008ReducingFunction::d2fYijdxidxj(const std::vector<double> &x, int i, int j, std::vector< std::vector< double> > * beta)
{
	double xi = x[i], xj = x[j], beta_Y = (*beta)[i][j], beta_Y2 = beta_Y*beta_Y;
	return (xi+xj)/(beta_Y2*xi+xj) + xj/(beta_Y2*xi+xj)*(1-(xi+xj)/(beta_Y2*xi+xj))
		+xi/(beta_Y2*xi+xj)*(1-beta_Y2*(xi+xj)/(beta_Y2*xi+xj))
		-xi*xj/pow(beta_Y2*xi+xj,(int)2)*(1+beta_Y2-2*beta_Y2*(xi+xj)/(beta_Y2*xi+xj));
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
double GERG2008DepartureFunction::phir(double tau, double delta)
{
	if (double_equal(tau,cache.phir.tau) && double_equal(delta,cache.phir.delta))
	{
		return cache.phir.cached_val;
	}
	else
	{
		double sum;
		if (using_gaussian){
			sum = phi1.base(tau, delta) + phi2.base(tau, delta);
		}
		else{
			sum = phi1.base(tau, delta);
		}
		cache.phir.tau = tau;
		cache.phir.delta = delta;
		cache.phir.cached_val = sum;
		return sum;
	}
}
double GERG2008DepartureFunction::dphir_dDelta(double tau, double delta)
{
	if (double_equal(tau,cache.dphir_dDelta.tau) && double_equal(delta,cache.dphir_dDelta.delta))
	{
		return cache.dphir_dDelta.cached_val;
	}
	else
	{
		double sum;
		if (using_gaussian){
			sum = phi1.dDelta(tau, delta) + phi2.dDelta(tau, delta);
		}
		else{
			sum = phi1.dDelta(tau, delta);
		}
		cache.dphir_dDelta.tau = tau;
		cache.dphir_dDelta.delta = delta;
		cache.dphir_dDelta.cached_val = sum;
		return sum;
	}
}
double GERG2008DepartureFunction::d2phir_dDelta2(double tau, double delta)
{
	if (double_equal(tau,cache.d2phir_dDelta2.tau) && double_equal(delta,cache.d2phir_dDelta2.delta))
	{
		return cache.d2phir_dDelta2.cached_val;
	}
	else
	{
		double sum;
		if (using_gaussian){
			sum = phi1.dDelta2(tau, delta) + phi2.dDelta2(tau, delta);
		}
		else{
			sum = phi1.dDelta2(tau, delta);
		}
		cache.d2phir_dDelta2.tau = tau;
		cache.d2phir_dDelta2.delta = delta;
		cache.d2phir_dDelta2.cached_val = sum;
		return sum;
	}
}
double GERG2008DepartureFunction::d2phir_dDelta_dTau(double tau, double delta)
{
	
	if (double_equal(tau,cache.d2phir_dDelta_dTau.tau) && double_equal(delta,cache.d2phir_dDelta_dTau.delta))
	{
		return cache.d2phir_dDelta_dTau.cached_val;
	}
	else
	{
		double sum;
		if (using_gaussian){
			sum = phi1.dDelta_dTau(tau, delta) + phi2.dDelta_dTau(tau, delta);
		}
		else{
			sum = phi1.dDelta_dTau(tau, delta);
		}
		cache.d2phir_dDelta_dTau.tau = tau;
		cache.d2phir_dDelta_dTau.delta = delta;
		cache.d2phir_dDelta_dTau.cached_val = sum;
		return sum;
	}
}
double GERG2008DepartureFunction::dphir_dTau(double tau, double delta)
{
	
	if (double_equal(tau,cache.dphir_dTau.tau) && double_equal(delta,cache.dphir_dTau.delta))
	{
		return cache.dphir_dTau.cached_val;
	}
	else
	{
		double sum;
		if (using_gaussian){
			sum = phi1.dTau(tau, delta) + phi2.dTau(tau, delta);
		}
		else{
			sum = phi1.dTau(tau, delta);
		}
		cache.dphir_dTau.tau = tau;
		cache.dphir_dTau.delta = delta;
		cache.dphir_dTau.cached_val = sum;
		return sum;
	}
}
double GERG2008DepartureFunction::d2phir_dTau2(double tau, double delta)
{
	
	if (double_equal(tau,cache.d2phir_dTau2.tau) && double_equal(delta,cache.d2phir_dTau2.delta))
	{
		return cache.d2phir_dTau2.cached_val;
	}
	else
	{
		double sum;
		if (using_gaussian){
			sum = phi1.dTau2(tau, delta) + phi2.dTau2(tau, delta);
		}
		else{
			sum = phi1.dTau2(tau, delta);
		}
		cache.d2phir_dTau2.tau = tau;
		cache.d2phir_dTau2.delta = delta;
		cache.d2phir_dTau2.cached_val = sum;
		return sum;
	}
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

	unsigned int iStart;
	// If sum of entries in eta are zero, there is no gaussian term
	// All terms are power-like
	if (double_equal(std::accumulate(eta.begin(), eta.end(), 0.0),0.0))
	{
		using_gaussian = false;
		iStart = eta.size()-1;
	}
	else
	{
		using_gaussian = true;
		iStart = 0;
		while (eta[iStart+1] < 0.25)
		{
			iStart++;
		}
	}

	phi1 = phir_power(n,d,t,1,iStart);
	if (using_gaussian) { phi2 = phir_GERG2008_gaussian(n,d,t,eta,epsilon,beta,gamma,iStart+1,eta.size()-1); }

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

double LemmonHFCDepartureFunction::phir(double tau, double delta)
{
	if (use_cache && double_equal(tau,cache.phir.tau) && double_equal(delta,cache.phir.delta))
	{
		return cache.phir.cached_val;
	}
	else
	{
		double sum = phi1.base(tau, delta);
		cache.phir.tau = tau;
		cache.phir.delta = delta;
		cache.phir.cached_val = sum;
		return sum;
	}
}
double LemmonHFCDepartureFunction::dphir_dDelta(double tau, double delta)
{
	if (use_cache && double_equal(tau,cache.dphir_dDelta.tau) && double_equal(delta,cache.dphir_dDelta.delta))
	{
		return cache.dphir_dDelta.cached_val;
	}
	else
	{
		double sum = phi1.dDelta(tau, delta);
		cache.dphir_dDelta.tau = tau;
		cache.dphir_dDelta.delta = delta;
		cache.dphir_dDelta.cached_val = sum;
		return sum;
	}
}
double LemmonHFCDepartureFunction::d2phir_dDelta2(double tau, double delta)
{
	if (use_cache && double_equal(tau,cache.d2phir_dDelta2.tau) && double_equal(delta,cache.d2phir_dDelta2.delta))
	{
		return cache.d2phir_dDelta2.cached_val;
	}
	else
	{
		double sum = phi1.dDelta2(tau, delta);
		cache.d2phir_dDelta2.tau = tau;
		cache.d2phir_dDelta2.delta = delta;
		cache.d2phir_dDelta2.cached_val = sum;
		return sum;
	}
}
double LemmonHFCDepartureFunction::d2phir_dDelta_dTau(double tau, double delta)
{
	
	if (use_cache && double_equal(tau,cache.d2phir_dDelta_dTau.tau) && double_equal(delta,cache.d2phir_dDelta_dTau.delta))
	{
		return cache.d2phir_dDelta_dTau.cached_val;
	}
	else
	{
		double sum = phi1.dDelta_dTau(tau, delta);
		cache.d2phir_dDelta_dTau.tau = tau;
		cache.d2phir_dDelta_dTau.delta = delta;
		cache.d2phir_dDelta_dTau.cached_val = sum;
		return sum;
	}
}
double LemmonHFCDepartureFunction::dphir_dTau(double tau, double delta)
{
	if (use_cache && double_equal(tau,cache.dphir_dTau.tau) && double_equal(delta,cache.dphir_dTau.delta))
	{
		return cache.dphir_dTau.cached_val;
	}
	else
	{
		double sum = phi1.dTau(tau, delta);
		cache.dphir_dTau.tau = tau;
		cache.dphir_dTau.delta = delta;
		cache.dphir_dTau.cached_val = sum;
		return sum;
	}
}
double LemmonHFCDepartureFunction::d2phir_dTau2(double tau, double delta)
{
	
	if (use_cache && double_equal(tau,cache.d2phir_dTau2.tau) && double_equal(delta,cache.d2phir_dTau2.delta))
	{
		return cache.d2phir_dTau2.cached_val;
	}
	else
	{
		double sum = phi1.dTau2(tau, delta);
		cache.d2phir_dTau2.tau = tau;
		cache.d2phir_dTau2.delta = delta;
		cache.d2phir_dTau2.cached_val = sum;
		return sum;
	}
}

void LemmonHFCDepartureFunction::set_coeffs_from_map(std::map<std::string,std::vector<double> > m)
{
	std::vector<double> n = m.find("n")->second;
	std::vector<double> t = m.find("t")->second;
	std::vector<double> d = m.find("d")->second;
	std::vector<double> l = m.find("l")->second;

	// All terms are power-like
	phi1 = phir_power(n,d,t,l,1,n.size()-1);
}

void LemmonAirHFCReducingFunction::set_coeffs_from_map(int i, int j, std::map<std::string,double > m)
{
	double xi_ij = m.find("xi")->second;
	double zeta_ij = m.find("zeta")->second;
	beta_T[i][j] = 1;
	beta_v[i][j] = 1;
	gamma_T[i][j] = (pFluids[i]->reduce.T + pFluids[j]->reduce.T + xi_ij)/(2*sqrt(pFluids[i]->reduce.T*pFluids[j]->reduce.T));
	double v_i = 1/pFluids[i]->reduce.rhobar;
	double v_j = 1/pFluids[j]->reduce.rhobar;
	double one_third = 0.33333333333333333;
	gamma_v[i][j] = (v_i + v_j + zeta_ij)/(0.25*pow(pow(v_i, one_third)+pow(v_j, one_third),(int)3));
}

ExcessTerm::ExcessTerm(int N)
{
	F.resize(N,std::vector<double>(N,1.0));
	DepartureFunctionMatrix.resize(N,std::vector<DepartureFunction*>(N));
	this->N = N;
}
ExcessTerm::~ExcessTerm()
{
	for (unsigned int i = 0; i < N; i++)
	{
		for (unsigned int j = 0; j < N; j++)
		{	
			if (DepartureFunctionMatrix[i][j] != NULL)
			{
				delete DepartureFunctionMatrix[i][j]; DepartureFunctionMatrix[i][j] = NULL;
			}
		}
	}
}

double ExcessTerm::phir(double tau, double delta, const std::vector<double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{	
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->phir(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::dphir_dTau(double tau, double delta, const std::vector<double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{	
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->dphir_dTau(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::dphir_dDelta(double tau, double delta, const std::vector<double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{	
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->dphir_dDelta(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2phir_dDelta2(double tau, double delta, const std::vector<double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{	
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2phir_dDelta2(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2phir_dTau2(double tau, double delta, const std::vector<double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{	
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2phir_dTau2(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2phir_dDelta_dTau(double tau, double delta, const std::vector<double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{	
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2phir_dDelta_dTau(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::dphir_dxi(double tau, double delta, const std::vector<double> &x, unsigned int i)
{
	double summer = 0;
	for (unsigned int k = 0; k < N; k++)
	{
		if (i != k)
		{
			summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->phir(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2phirdxidxj(double tau, double delta, const std::vector<double> &x, unsigned int i, unsigned int j)
{
	if (i != j)
	{
		return F[i][j]*DepartureFunctionMatrix[i][j]->phir(tau,delta);
	}
	else
	{
		return 0;
	}
}
double ExcessTerm::d2phir_dxi_dTau(double tau, double delta, const std::vector<double> &x, unsigned int i)
{
	double summer = 0;
	for (unsigned int k = 0; k < N; k++)
	{
		if (i != k)
		{
			summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->dphir_dTau(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2phir_dxi_dDelta(double tau, double delta, const std::vector<double> &x, unsigned int i)
{
	double summer = 0;
	for (unsigned int k = 0; k < N; k++)
	{
		if (i != k)
		{
			summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->dphir_dDelta(tau,delta);
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
	this->N = pFluids.size();
}
double ResidualIdealMixture::phir(double tau, double delta, const std::vector<double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		summer += x[i]*pFluids[i]->phir(tau,delta);
	}
	return summer;
}
double ResidualIdealMixture::dphir_dDelta(double tau, double delta, const std::vector<double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		summer += x[i]*pFluids[i]->dphir_dDelta(tau,delta);
	}
	return summer;
}
double ResidualIdealMixture::d2phir_dDelta2(double tau, double delta, const std::vector<double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		summer += x[i]*pFluids[i]->d2phir_dDelta2(tau,delta);
	}
	return summer;
}
double ResidualIdealMixture::d2phir_dDelta_dTau(double tau, double delta, const std::vector<double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		summer += x[i]*pFluids[i]->d2phir_dDelta_dTau(tau,delta);
	}
	return summer;
}
double ResidualIdealMixture::dphir_dTau(double tau, double delta, const std::vector<double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		summer += x[i]*pFluids[i]->dphir_dTau(tau,delta);
	}
	return summer;
}
double ResidualIdealMixture::d2phir_dTau2(double tau, double delta, const std::vector<double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		summer += x[i]*pFluids[i]->d2phir_dTau2(tau,delta);
	}
	return summer;
}
double ReducingFunction::d_ndTrdni_dxj__constxi(const std::vector<double> &x, int i, int j)
{
	double s = 0;
	for (unsigned int k = 0; k < N; k++)
	{
		s += x[k]*d2Trdxidxj(x,j,k);
	}
	return d2Trdxidxj(x,i,j)-dTrdxi__constxj(x,j)-s;
}
double ReducingFunction::d_ndrhorbardni_dxj__constxi(const std::vector<double> &x, int i, int j)
{
	double s = 0;
	for (unsigned int k = 0; k < N; k++)
	{
		s += x[k]*d2rhorbardxidxj(x,j,k);
	}
	return d2rhorbardxidxj(x,j,i)-drhorbardxi__constxj(x,j)-s;
}
double ReducingFunction::ndrhorbardni__constnj(const std::vector<double> &x, int i)
{
	double summer_term1 = 0;
	for (unsigned int j = 0; j < N; j++)
	{
		summer_term1 += x[j]*drhorbardxi__constxj(x,j);
	}
	return drhorbardxi__constxj(x,i)-summer_term1;
}
double ReducingFunction::ndTrdni__constnj(const std::vector<double> &x, int i)
{
	// GERG Equation 7.54
	double summer_term1 = 0;
	for (unsigned int j = 0; j < N; j++)
	{
		summer_term1 += x[j]*dTrdxi__constxj(x,j);
	}
	return dTrdxi__constxj(x,i)-summer_term1;
}

double SuccessiveSubstitutionVLE::call(double beta, double T, double p, const std::vector<double> &z, std::vector<double> &K)
{
	int iter = 1;
	double change, f, dfdT, rhobar_liq_new, rhobar_vap_new;
	unsigned int N = z.size();
	K.resize(N); ln_phi_liq.resize(N); ln_phi_vap.resize(N); x.resize(N); y.resize(N);

	Mix->x_and_y_from_K(beta, K, z, x, y);
	
	rhobar_liq = 1.3 * Mix->rhobar_pengrobinson(T, p, x, PR_SATL); // [kg/m^3]
	rhobar_vap = Mix->rhobar_pengrobinson(T, p, y, PR_SATV); // [kg/m^3]

	do
	{
		rhobar_liq_new = Mix->rhobar_Tpz(T, p, x, rhobar_liq); // [kg/m^3]
		rhobar_vap_new = Mix->rhobar_Tpz(T, p, y, rhobar_vap); // [kg/m^3]
		rhobar_liq = rhobar_liq_new;
		rhobar_vap = rhobar_vap_new;

		double Tr_liq = Mix->pReducing->Tr(x); // [K]
		double Tr_vap = Mix->pReducing->Tr(y);  // [K]
		double tau_liq = Tr_liq/T; // [-]
		double tau_vap = Tr_vap/T; // [-]

		double rhorbar_liq = Mix->pReducing->rhorbar(x); //[kg/m^3]
		double rhorbar_vap = Mix->pReducing->rhorbar(y); //[kg/m^3]
		double delta_liq = rhobar_liq/rhorbar_liq;  //[-]
		double delta_vap = rhobar_vap/rhorbar_vap;  //[-] 

		f = 0;
		dfdT = 0;

		for (unsigned int i=0; i < N; i++)
		{
			// Loop over the liquids to take advantage of caching since cached derivatives 
			// will be used as long as tau and delta don't change
			ln_phi_liq[i] = Mix->ln_fugacity_coefficient(tau_liq, delta_liq, x, i);
			double dln_phi_liq_dT = Mix->dln_fugacity_coefficient_dT__constp_n(tau_liq, delta_liq, x, i);

			// And then loop over the vapors
			ln_phi_vap[i] = Mix->ln_fugacity_coefficient(tau_vap, delta_vap, y, i);
			double dln_phi_vap_dT = Mix->dln_fugacity_coefficient_dT__constp_n(tau_vap, delta_vap, y, i);

			K[i] = exp(ln_phi_liq[i]-ln_phi_vap[i]);
			
			f += z[i]*(K[i]-1)/(1-beta+beta*K[i]);

			double dfdK = K[i]*z[i]/pow(1-beta+beta*K[i],(int)2);
			dfdT += dfdK*(dln_phi_liq_dT-dln_phi_vap_dT);
		}
		
		change = -f/dfdT;

		//double Tnew_over_Told = (T + change)/T;
		T += change;

		Mix->x_and_y_from_K(beta,K,z,x,y);

		iter += 1;
		if (iter > 50)
		{
			return _HUGE;
			//throw ValueError(format("saturation_p was unable to reach a solution within 50 iterations"));
		}
		// If SS takes a large step, help by re-guessing using Peng-Robinson
		if (fabs(change) > 3)
		{
			rhobar_liq = 1.3 * Mix->rhobar_pengrobinson(T, p, x, PR_SATL); // [kg/m^3]
			rhobar_vap = Mix->rhobar_pengrobinson(T, p, y, PR_SATV); // [kg/m^3]
		}
	}
	while(fabs(f) > 1e-3 && iter < Nstep_max);
	
	if (!useNR)
	{
		this->rhobar_liq = rhobar_liq;
		this->rhobar_vap = rhobar_vap;
		return T;
	}
	else
	{
		// Pass off to Newton-Raphson to polish the solution
		double Treturn = Mix->NRVLE.call(beta, T, p, rhobar_liq, rhobar_vap, z, K, N+1, log(p));
		this->rhobar_liq = Mix->NRVLE.rhobar_liq;
		this->rhobar_vap = Mix->NRVLE.rhobar_vap;
		return Treturn;
	}
}

void NewtonRaphsonVLE::resize(unsigned int N)
{
	this->N = N;
	x.resize(N); 
	y.resize(N); 
	phi_ij_liq.resize(N); 
	phi_ij_vap.resize(N);

	r.resize(N+2);
	J.resize(N+2, std::vector<double>(N+2, 0));

	neg_dFdS.resize(N+2);
	dXdS.resize(N+2);

	// Fill the vector -dFdS with zeros (Gerg Eqn. 7.132)
	std::fill(neg_dFdS.begin(), neg_dFdS.end(), (double)0.0);
	// Last entry is 1
	neg_dFdS[N+1] = 1.0;
}
double NewtonRaphsonVLE::call(double beta, double T, double p, double rhobar_liq, double rhobar_vap, const std::vector<double> &z, std::vector<double> & K, int spec_index, double spec_value)
{
	int iter = 0;

	// Reset all the variables and resize
	pre_call();
	resize(K.size());

	//std::cout << format("NRVLE beta %g, T %g, p %g, rhobar_liq %g, rhobar_vap %g , z %s, K %s, index %d, val %g \n",beta,T,p,rhobar_liq,rhobar_vap,vec_to_string(z,"%4.3g").c_str(),vec_to_string(K,"%4.3g").c_str(),spec_index,spec_value);

	//~ double old_rhobar_liq = rhobar_liq, 
		   //~ old_rhobar_vap = rhobar_vap,
		   //~ old_T = T,
		   //~ old_p = p;
	std::vector<double> old_K = K;

	this->rhobar_liq = rhobar_liq;
	this->rhobar_vap = rhobar_vap;
	do
	{
		// Build the Jacobian and residual vectors for given inputs of K_i,T,p
		build_arrays(beta,T,p,z,K,spec_index,spec_value);

		// Solve for the step; v is the step with the contents 
		// [delta(lnK0), delta(lnK1), ..., delta(lnT), delta(lnp)]
		std::vector<double> v = linsolve(J, r);

		// Set the variables again, the same structure independent of the specified variable
		for (unsigned int i = 0; i < N; i++)
		{
			K[i] = exp(log(K[i]) + v[i]);
			if (!ValidNumber(K[i]))
			{
				throw ValueError(format("K[i] (%g) is invalid",K[i]).c_str());
			}
		}
		T = exp(log(T) + v[N]);
		p = exp(log(p) + v[N+1]);

		if (fabs(T) > 1e6 || fabs(p) > 1e10)
		{
			/*std::cout << "J = " << vec_to_string(J,"%16.15g");
			std::cout << "nr = " << vec_to_string(r,"%16.15g");*/
			throw ValueError("Temperature or p has bad value");
		}
		
		//std::cout << iter << " " << T << " " << p << " " << error_rms << std::endl;
		iter++;
	}
	while(this->error_rms > 1e-11 && iter < Nsteps_max);
	// Store new values since they were passed by value
	this->T = T;
	this->p = p;
	this->K = K;
	this->Nsteps = iter;

	//std::cout << format("NRVLE beta %g, T %g, p %g, rhobar_liq %g, rhobar_vap %g , z %s, K %s, index %d, val %g \n",beta,T,p,rhobar_liq,rhobar_vap,vec_to_string(z,"%4.3g").c_str(),vec_to_string(K,"%4.3g").c_str(),spec_index,exp(spec_value));
	//std::cout << format("liq_change %g vap_change %g\n", this->rhobar_liq/old_rhobar_liq-1, this->rhobar_vap/old_rhobar_vap-1);
	return T;
}

void NewtonRaphsonVLE::build_arrays(double beta, double T, double p, const std::vector<double> &z, std::vector<double> &K, int spec_index, double spec_value)
{
	// Step 0:
	// --------
	// Calculate the mole fractions in liquid and vapor phases
	Mix->x_and_y_from_K(beta, K, z, x, y);

	// Step 1:
	// -------
	// Calculate the new reducing and reduced parameters for each phase
	// based on the current values of the molar fractions
	double old_rhobar_liq = this->rhobar_liq;
	double old_rhobar_vap = this->rhobar_vap;
	
	this->rhobar_liq = Mix->rhobar_Tpz(T, p, x, this->rhobar_liq); // [kg/m^3] (Not exact due to solver convergence)
	this->rhobar_vap = Mix->rhobar_Tpz(T, p, y, this->rhobar_vap); // [kg/m^3] (Not exact due to solver convergence)

	if (this->rhobar_vap > this->rhobar_liq)
	{
		// Phases are inverted, swap x&y, recalculate densities
		std::swap(x, y);
		std::cout << format("Phase inversion");
		this->rhobar_liq = Mix->rhobar_Tpz(T, p, x, old_rhobar_liq); // [kg/m^3] (Not exact due to solver convergence)
		this->rhobar_vap = Mix->rhobar_Tpz(T, p, y, old_rhobar_vap); // [kg/m^3] (Not exact due to solver convergence)
	}

	//std::cout << format("rho liq: %g rho vap: %g k: %s T: %g P: %g x: %s\n", this->rhobar_liq, this->rhobar_vap, vec_to_string(K,"%6.5g").c_str(),T,p, vec_to_string(x,"%6.5g").c_str()).c_str();

	if (!ValidNumber(this->rhobar_liq)){ throw ValueError(format("Liquid density solver has failed with guess value %g",old_rhobar_liq).c_str()); }
	if (!ValidNumber(this->rhobar_vap)){ throw ValueError(format("Vapor density solver has failed with guess value %g",old_rhobar_vap).c_str()); }
	
	double Tr_liq = Mix->pReducing->Tr(x); // [K]
	double Tr_vap = Mix->pReducing->Tr(y);  // [K]
	double tau_liq = Tr_liq/T; // [-]
	double tau_vap = Tr_vap/T; // [-]

	double rhorbar_liq = Mix->pReducing->rhorbar(x); //[kg/m^3]
	double rhorbar_vap = Mix->pReducing->rhorbar(y); //[kg/m^3]
	double delta_liq = this->rhobar_liq/rhorbar_liq;  //[-]
	double delta_vap = this->rhobar_vap/rhorbar_vap;  //[-]

	double p_liq = this->rhobar_liq*Mix->Rbar(x)*T*(1+delta_liq*Mix->dphir_dDelta(tau_liq,delta_liq,x));
	double p_vap = this->rhobar_vap*Mix->Rbar(y)*T*(1+delta_vap*Mix->dphir_dDelta(tau_vap,delta_vap,y));

	// Step 2:
	// -------
	// Build the residual vector and the Jacobian matrix

	// For the residuals F_i
	for (unsigned int i = 0; i < N; i++)
	{
		double ln_phi_liq = Mix->ln_fugacity_coefficient(tau_liq, delta_liq, x, i);
		double phi_iT_liq = Mix->dln_fugacity_coefficient_dT__constp_n(tau_liq, delta_liq, x, i);
		double phi_ip_liq = Mix->dln_fugacity_coefficient_dp__constT_n(tau_liq, delta_liq, x, i);
		for (unsigned int j = 0; j < N; j++)
		{
			phi_ij_liq[j] = Mix->ndln_fugacity_coefficient_dnj__constT_p(tau_liq,delta_liq,x,i,j); // 7.126 from GERG monograph
		}

		double ln_phi_vap = Mix->ln_fugacity_coefficient(tau_vap, delta_vap, y, i);
		double phi_iT_vap = Mix->dln_fugacity_coefficient_dT__constp_n(tau_vap, delta_vap, y, i);
		double phi_ip_vap = Mix->dln_fugacity_coefficient_dp__constT_n(tau_vap, delta_vap, y, i);
		for (unsigned int j = 0; j < N; j++)
		{
			phi_ij_vap[j] = Mix->ndln_fugacity_coefficient_dnj__constT_p(tau_vap,delta_vap,y,i,j); // 7.126 from GERG monograph
		}
		
		r[i] = ln_phi_vap - ln_phi_liq + log(K[i]);
		for (unsigned int j = 0; j < N; j++)
		{	
			J[i][j] = K[j]*z[j]/pow(1-beta+beta*K[j],(int)2)*((1-beta)*phi_ij_vap[j]+beta*phi_ij_liq[j])+Kronecker_delta(i,j);
		}
		// dF_{i}/d(ln(T))
		J[i][N] = T*(phi_iT_vap-phi_iT_liq);
		// dF_{i}/d(ln(p))
		J[i][N+1] = p*(phi_ip_vap-phi_ip_liq);
	}

	double summer1 = 0;
	for (unsigned int i = 0; i < N; i++)
	{
		// Although the definition of this term is given by 
		// y[i]-x[i], when x and y are normalized, you get 
		// the wrong values.  Why? No idea.
		summer1 += z[i]*(K[i]-1)/(1-beta+beta*K[i]); 
	}
	r[N] = summer1;

	// For the residual term F_{N+1}, only non-zero derivatives are with respect
	// to ln(K[i])
	for (unsigned int j = 0; j < N; j++)
	{
		J[N][j] = K[j]*z[j]/pow(1-beta+beta*K[j],(int)2);
	}
	
	// For the specification term F_{N+2}, with index of N+1
	if (spec_index == N)
		{r[N+1] = log(T) - spec_value;}
	else if (spec_index == N+1)
		{r[N+1] = log(p) - spec_value;}
	else
		{r[N+1] = log(K[spec_index]) - spec_value;}

	// The row in the Jacobian for the specification
	for (unsigned int i = 0; i < N+2; i++)
	{
		J[N+1][i] = Kronecker_delta(i,spec_index);
	}

	// Flip all the signs of the entries in the residual vector since we are solving Jv = -r, not Jv=r
	// Also calculate the rms error of the residual vector at this step
	error_rms = 0;
	for (unsigned int i = 0; i < N+2; i++)
	{
		r[i] *= -1;
		error_rms += r[i]*r[i]; // Sum the squares
	}
	error_rms = sqrt(error_rms); // Square-root (The R in RMS)

	//std::cout << format("r: %s\n",vec_to_string(r,"%6.5g").c_str());
	
	// Step 3:
	// =======
	// Calculate the sensitivity vector dXdS
	dXdS = linsolve(J,neg_dFdS);
}


void PhaseEnvelope::build(double p0, const std::vector<double> &z, double beta_envelope)
{
	double S, Sold, DELTAS, abs_Kdeviation;
	K.resize(Mix->N);
	bubble.K.resize(Mix->N);
	bubble.lnK.resize(Mix->N);
	bubble.x.resize(Mix->N);
	bubble.y.resize(Mix->N);
	dew.K.resize(Mix->N);
	dew.lnK.resize(Mix->N);
	dew.x.resize(Mix->N);
	dew.y.resize(Mix->N);

	for (unsigned int betaI = 0; betaI < 2; betaI++)
	{
		// If betaI is 0, use beta_envelope, otherwise use 1-beta_envelope
		double beta = (betaI==0) ? beta_envelope : 1-beta_envelope;
		
		// Reference for data points to either bubble or dew
		PhaseEnvelopeLog &data = (betaI==0) ? bubble : dew;

		// Use the preconditioner to get the very rough guess for saturation temperature 
		// (interpolation based on pseudo-critical pressure)
		double T_guess = Mix->saturation_p_preconditioner(p0, z);

		// Initialize the call using Wilson to get K-factors and temperature
		double T = Mix->saturation_p_Wilson(beta, p0, z, T_guess, K);

		// Call successive substitution and Newton-Raphson to get updated guess for K-factors using imposed pressure
		T = Mix->SS.call(beta, T, p0, z, K);
		rhobar_liq = Mix->SS.rhobar_liq;
		rhobar_vap = Mix->SS.rhobar_vap;

		double p = p0;
		double lnT = log(T);
		double lnP = log(p0);

		// Find the most sensitive input (the one with the largest absolute value in dX/dS
		double max_abs_val = -1;
		int i_S = -1;
		for (unsigned int i = 0; i < Mix->N+1; i++) // N+1 because specification not allowed to be selected
		{
			if (fabs(Mix->NRVLE.dXdS[i])>max_abs_val)
			{
				max_abs_val = fabs(Mix->NRVLE.dXdS[i]);
				i_S = i;
			}
		}
		// Determine the value of the specified variable based on the index of the specified variable
		if (i_S >= (int)Mix->N)
		{
			if (i_S == Mix->N) 
			{ 
				Sold = lnT; 
				DELTAS = log(1.01); // Temperatures increase as we approach the critical point for both branches
			}
			else
			{ 
				Sold = lnP; 
				DELTAS = log(1.01); // Pressures increase as we approach the critical point for both branches
			}
		}
		else
		{
			// K will go towards 1 as we move up on the envelope
			Sold = log(K[i_S]);
			if (Sold < 0) // K_i < 1
			{
				DELTAS = log(1.1); // K_i will increase
			}
			else
			{
				DELTAS = log(0.9); // K_i will decrease
			}
		}

		// Run once with the specified variable set
		Mix->NRVLE.call(beta, T, p, rhobar_liq, rhobar_vap, z, K, i_S, Sold);

		// Store the variables in the log if the step worked ok
		data.store_variables(T,p,Mix->NRVLE.rhobar_liq,Mix->NRVLE.rhobar_vap,K,i_S, Mix->NRVLE.x, Mix->NRVLE.y, Mix->N);

		int iter = 1;
		do
		{
			bool step_accepted = false;
			std::vector<double> dXdSold = Mix->NRVLE.dXdS; //Copy from the last good run
			double Told = T; // Copy from the last good run
			double pold = p; // Copy from the last good run
			std::vector<double> Kold = K; // Copy from the last good run

			// Loop while the step size isn't small enough - usually only requires one downsizing of the step
			while (step_accepted == false)
			{
				// Make a copy of the stored values from the last iteration
				std::vector<double> DELTAX(dXdSold.size()), dXdS = dXdSold;
				T = Told;
				p = pold;
				K = Kold;
				
				for (unsigned int i = 0; i < Mix->N+2; i++)
				{
					DELTAX[i] = dXdS[i]*DELTAS;
				}

				// Update the temperature
				double DELTALNT = DELTAX[Mix->N];
				T = exp(log(Told)+DELTALNT);
				lnT = log(T);

				// Update the pressure
				double DELTALNP = DELTAX[Mix->N+1];
				p = exp(log(pold)+DELTALNP);
				lnP = log(p);

				// Update the K-factors
				for (unsigned int i = 0; i < Mix->N; i++)
				{
					K[i] = exp(log(Kold[i])+DELTAX[i]);
				}

				// Update the specified variable
				S = Sold + DELTAS;
				
				//std::cout << format("S: %g Sold %g DELTAS %g exp(S) %g exp(Sold) %g K %s\n",S,Sold,DELTAS, exp(S), exp(Sold),vec_to_string(K,"%0.5g").c_str());

				// UPDATE THE GUESSES
				// The specified variable is known directly. The others must be obtained either through the use of dXdS and/or extrapolation
				if (data.T.size() > 2)
				{
					int M = data.T.size();
					// Update the densities using quadratic extrapolation
					rhobar_liq = QuadInterp(data.lnK[1][M-3],data.lnK[1][M-2],data.lnK[1][M-1],data.rhobar_liq[M-3],data.rhobar_liq[M-2],data.rhobar_liq[M-1],log(K[1]));
					rhobar_vap = exp(QuadInterp(data.lnK[1][M-3],data.lnK[1][M-2],data.lnK[1][M-1],log(data.rhobar_vap[M-3]),log(data.rhobar_vap[M-2]),log(data.rhobar_vap[M-1]),log(K[1])));
					T = exp(QuadInterp(data.lnK[1][M-3],data.lnK[1][M-2],data.lnK[1][M-1],data.lnT[M-3],data.lnT[M-2],data.lnT[M-1],log(K[1])));
					p = exp(QuadInterp(data.lnK[1][M-3],data.lnK[1][M-2],data.lnK[1][M-1],data.lnp[M-3],data.lnp[M-2],data.lnp[M-1],log(K[1])));
				}
				else if (data.T.size() > 3)
				{
					int M = data.T.size();
					// First we use a cubic spline to find the next value for the vapor molar density
					SplineClass SC1;
					SC1.add_4value_constraints(data.lnK[0][M-4],data.lnK[0][M-3],data.lnK[0][M-2],data.lnK[0][M-1],
											   data.lnrhobar_vap[M-4],data.lnrhobar_vap[M-3],data.lnrhobar_vap[M-2],data.lnrhobar_vap[M-1]);
					SC1.build();
					rhobar_vap = exp(SC1.evaluate(log(K[0])));

					SplineClass SC2;
					SC2.add_4value_constraints(log(data.rhobar_vap[M-4]), log(data.rhobar_vap[M-3]), log(data.rhobar_vap[M-2]), log(data.rhobar_vap[M-1]),
											   data.lnT[M-4], data.lnT[M-3], data.lnT[M-2], data.lnT[M-1]);
					SC2.build();
					T = exp(SC2.evaluate(log(rhobar_vap)));

					SplineClass SC3;
					SC3.add_4value_constraints(log(data.rhobar_vap[M-4]), log(data.rhobar_vap[M-3]), log(data.rhobar_vap[M-2]), log(data.rhobar_vap[M-1]),
											   data.lnp[M-4], data.lnp[M-3], data.lnp[M-2], data.lnp[M-1]);
					SC3.build();
					p = exp(SC3.evaluate(log(rhobar_vap)));

					SplineClass SC4;
					/*SC4.add_4value_constraints(log(data.rhobar_vap[M-4]), log(data.rhobar_vap[M-3]), log(data.rhobar_vap[M-2]), log(data.rhobar_vap[M-1]),
											   log(data.rhobar_liq[M-4]), log(data.rhobar_liq[M-3]), log(data.rhobar_liq[M-2]), log(data.rhobar_liq[M-1]));
					SC4.build();
					rhobar_liq = exp(SC4.evaluate(log(rhobar_vap)));*/
					SC4.add_4value_constraints(data.rhobar_vap[M-4], data.rhobar_vap[M-3], data.rhobar_vap[M-2], data.rhobar_vap[M-1],
											   data.rhobar_liq[M-4], data.rhobar_liq[M-3], data.rhobar_liq[M-2], data.rhobar_liq[M-1]);
					SC4.build();
					rhobar_liq = SC4.evaluate(rhobar_vap);
				}
				else
				{
					// Treat the liquid as being incompressible, don't update the guess
					// Treat the vapor as being ideal, use pressure ratio to update pressure

					// Start with T&p from the last specified state
					T = Mix->NRVLE.T;
					p = Mix->NRVLE.p;

					// Start with the densities from the last specified state
					rhobar_liq = Mix->NRVLE.rhobar_liq;
					rhobar_vap = Mix->NRVLE.rhobar_vap;
				}

				// Run with the selected specified variable
				try
				{
					
					Mix->NRVLE.call(beta, T, p, rhobar_liq, rhobar_vap, z, K, i_S, S);
					// Throw an exception if an invalid value returned but no exception was thrown
					if (!ValidNumber(Mix->NRVLE.T))
					{
						throw ValueError(format("T [%g] is not valid", Mix->NRVLE.T).c_str());
					}
					
					T = Mix->NRVLE.T;
					p = Mix->NRVLE.p;
					rhobar_liq = Mix->NRVLE.rhobar_liq;
					rhobar_vap = Mix->NRVLE.rhobar_vap;

					//std::vector<double> x(K.size()), y(K.size());
					//double T2 = Mix->NRVLE.T;
					//double T3 = Mix->saturation_p(beta_envelope,p,z,x,y);
					//std::cout << T2 << " " << T3;
				}
				catch (CoolPropBaseError &e)
				{
					printf("error: %s\n",e.what());
					// Decrease the step size by a factor of 10
					printf("Failure; step downsize\n");

					DELTAS *= 0.1;
					
					continue;
				}

				#ifdef MPLSUPPORTED
				if (K[0] < 1.01)
				{
				Dictionary d;
				d.add("lw", 1);
				d.add("color", "r");
				d.add("linestyle", "-");
				d.add("marker", "o");

				PyPlotter plt;
				plt.plot(data.lnK[0],data.lnrhobar_vap, &d);
				plt.show();
				}
				#endif

				// Store the variables in the log if the step worked ok
				data.store_variables(T, p, Mix->NRVLE.rhobar_liq, Mix->NRVLE.rhobar_vap, K, i_S, Mix->NRVLE.x, Mix->NRVLE.y, Mix->N);

				step_accepted = true;

				if (Mix->NRVLE.Nsteps > 4)
				{
					std::cout << "Downsizing, too many iterations\n";
					DELTAS *= 0.25;
				}
				else if (Mix->NRVLE.Nsteps < 4)
				{
					DELTAS *= 2;
				}

				double _max_abs_val = -1;
				int iii_S = -1;
				for (unsigned int i = 0; i < Mix->N+1; i++)
				{
					if (fabs(Mix->NRVLE.dXdS[i]) > _max_abs_val)
					{
						_max_abs_val = fabs(Mix->NRVLE.dXdS[i]);
						iii_S = i;
					}
				}
				//std::cout << format("T,P,Nstep,K : %g %g %d %g %g %d %d %s %s\n",T,p,Mix->NRVLE.Nsteps, rhobar_liq, rhobar_vap, iii_S, i_S, vec_to_string(Mix->NRVLE.x,"%6.5g").c_str(), vec_to_string(K,"%6.5g").c_str());
			}
			
			// Update step counter
			iter++;

			// Reset last specified value
			Sold = S;
			
			if (!ValidNumber(T))
			{
				double rr = 0;
			}

			std::vector<double> Kdeviation(K.size(),0);
			for (unsigned int i = 0; i < Mix->N; i++)
			{
				Kdeviation[i] = abs(K[i]-1);
			}
			abs_Kdeviation = *std::max_element(Kdeviation.begin(), Kdeviation.end());
		}
		while ((p > p0 && abs_Kdeviation > 0.01 && iter < 10000) || iter < 3);
	}
}

void PhaseEnvelope::to_python_files(std::string base_fname)
{
	FILE *fp;
	for (unsigned int i = 0; i < 2; i++)
	{	
		// reference for data points to bubble if i==0, otherwise dew
		PhaseEnvelopeLog &data = (i==0) ? bubble : dew;
		std::string fname;
		if (i==0){ 
			fname = base_fname+std::string("_bubble.py");
		}
		else{
			fname = base_fname+std::string("_dew.py");
		}
		fp = fopen(fname.c_str(), "w");
		
		fprintf(fp, "import matplotlib.pyplot as plt\n");
		fprintf(fp, "T = %s\n",vec_to_string(data.T,"%10.9g").c_str());
		fprintf(fp, "p = %s\n",vec_to_string(data.p,"%10.9g").c_str());
		fprintf(fp, "rhobar_liq = %s\n", vec_to_string(data.rhobar_liq,"%10.9g").c_str());
		fprintf(fp, "rhobar_vap = %s\n", vec_to_string(data.rhobar_vap,"%10.9g").c_str());
		fprintf(fp, "K0 = %s\n", vec_to_string(data.K[0],"%10.9g").c_str());
		fprintf(fp, "K1 = %s\n", vec_to_string(data.K[1],"%10.9g").c_str());
		fprintf(fp, "lnK0 = %s\n", vec_to_string(data.lnK[0],"%10.9g").c_str());
		fprintf(fp, "lnK1 = %s\n", vec_to_string(data.lnK[1],"%10.9g").c_str());
		fprintf(fp, "if __name__=='__main__':\n\tplt.plot(T,p,'o-')\n\tplt.show()");
		fprintf(fp, "\n\tplt.plot(K0,p,'o-')\n\tplt.show()");
		fclose(fp);
	}
}



/////    #####################   TESTS #######################

#ifndef DISABLE_CATCH
#include "Catch/catch.hpp"
TEST_CASE("Mixture derivative checks", "")
	{
		Mixture Mix("Methane|Ethane");
		std::vector<double> z(2,0.5);
		double rhorbar = Mix.pReducing->rhorbar(z);
		double Tr = Mix.pReducing->Tr(z);
		int i = 0;
		double epsT = sqrt(DBL_EPSILON), epsp = 0.1, epsz = 1e-5;
		double T0 = 300, p0 = 100000;
		double rho0 = Mix.rhobar_Tpz(T0,p0, z, 5);
		double delta0 = rho0/rhorbar;
		double tau0 = Tr/T0;

		SECTION("check against RP9.1")
		{
			Mix.check_MethaneEthane();
		}
		SECTION("check p")
		{
			double p1 = rho0*Mix.Rbar(z)*T0*(1+delta0*Mix.dphir_dDelta(tau0, delta0, z));
			CAPTURE(p0);
			CAPTURE(p1);
			REQUIRE(fabs(p0/p1-1) < 1e-8);
		}
		SECTION("dlnphi_dT")
		{
			double ANA = Mix.dln_fugacity_coefficient_dT__constp_n(tau0, delta0, z, i);

			double T1 = T0+epsT;
			double delta1 = Mix.rhobar_Tpz(T1,p0,z, 5)/rhorbar;
			double tau1 = Tr/T1;
			double f1 = Mix.ln_fugacity_coefficient(tau1, delta1, z, i);

			double T2 = T0-epsT;
			double delta2 = Mix.rhobar_Tpz(T2,p0,z, 5)/rhorbar;
			double tau2 = Tr/T2;
			double f2 = Mix.ln_fugacity_coefficient(tau2, delta2, z, i);
			
			double NUM = (f1-f2)/(2*epsT);
			
			CAPTURE(NUM);
			CAPTURE(ANA);
			REQUIRE(fabs(NUM/ANA-1) < 1e-8);
		}
		SECTION("dlnphi_dp")
		{
			double ANA = Mix.dln_fugacity_coefficient_dp__constT_n(tau0, delta0, z, i);

			double p1 = p0+epsp;
			double delta1 = Mix.rhobar_Tpz(T0,p1,z, 5)/rhorbar;
			double tau1 = Tr/T0;
			double f1 = Mix.ln_fugacity_coefficient(tau1, delta1, z, i);

			double p2 = p0-epsp;
			double delta2 = Mix.rhobar_Tpz(T0,p2,z, 5)/rhorbar;
			double tau2 = Tr/T0;
			double f2 = Mix.ln_fugacity_coefficient(tau2, delta2, z, i);
			
			double NUM = (f1-f2)/(2*epsp);
			
			CAPTURE(NUM);
			CAPTURE(ANA);
			REQUIRE(fabs(NUM/ANA-1) < 1e-7);
		}
		SECTION("dln_fugacity_coefficient_dxj__constT_p_xi")
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					double ANA = Mix.dln_fugacity_coefficient_dxj__constT_p_xi(tau0, delta0, z, i, j);

					std::vector<double> z1 = z;
					z1[j] += epsz;
					double delta1 = Mix.rhobar_Tpz(T0, p0, z1, 5)/Mix.pReducing->rhorbar(z1);
					double tau1 = Mix.pReducing->Tr(z1)/T0;
					double f1 = Mix.ln_fugacity_coefficient(tau1, delta1, z1, i);

					std::vector<double> z2 = z;
					z2[j] -= epsz;
					double delta2 = Mix.rhobar_Tpz(T0, p0, z2, 5)/Mix.pReducing->rhorbar(z2);
					double tau2 = Mix.pReducing->Tr(z2)/T0;
					double f2 = Mix.ln_fugacity_coefficient(tau2, delta2, z2, i);
					
					double NUM = (f1-f2)/(2*epsz);
					
					CAPTURE(NUM);
					CAPTURE(ANA);
					CAPTURE(i);
					CAPTURE(j);
					CHECK(fabs(NUM-ANA) < 1e-7);
				}
			}
		}
		SECTION("d_ndphirdni_dxj__constdelta_tau_xi")
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					double ANA = Mix.d_ndphirdni_dxj__constdelta_tau_xi(tau0, delta0, z, i, j);
					std::vector<double> z1 = z, z2 = z;
					
					z1[j] += epsz;
					double f1 = Mix.ndphir_dni__constT_V_nj(tau0, delta0, z1, i);

					z2[j] -= epsz;
					double f2 = Mix.ndphir_dni__constT_V_nj(tau0, delta0, z2, i);

					double NUM = (f1-f2)/(2*epsz);
					CAPTURE(NUM);
					CAPTURE(ANA);
					CAPTURE(i);
					CAPTURE(j);
					CHECK(fabs(NUM/ANA-1) < 1e-10);
				}
			}
		}
		SECTION("dphir_dxj__constT_V_xi")
		{
			for (int j = 0; j < 2; j++)
			{
				double ANA = Mix.dphir_dxj__constT_V_xi(tau0, delta0, z, j);
				
				std::vector<double> z1 = z;
				z1[j] += epsz;
				double delta1 = rho0/Mix.pReducing->rhorbar(z1);
				double tau1 = Mix.pReducing->Tr(z1)/T0;
				double f1 = Mix.phir(tau1, delta1, z1);

				std::vector<double> z2 = z;
				z2[j] -= epsz;
				double delta2 = rho0/Mix.pReducing->rhorbar(z2);
				double tau2 = Mix.pReducing->Tr(z2)/T0;
				double f2 = Mix.phir(tau2, delta2, z2);

				double NUM = (f1-f2)/(2*epsz);
				CAPTURE(NUM);
				CAPTURE(ANA);
				CAPTURE(j);
				CHECK(fabs(NUM/ANA-1) < 1e-9);
			}
		}
		SECTION("d_dphirddelta_dxj__constT_V_xi")
		{
			for (int j = 0; j < 2; j++)
			{
				double ANA = Mix.d_dphirddelta_dxj__constT_V_xi(tau0, delta0, z, j);
				
				std::vector<double> z1 = z;
				z1[j] += epsz;
				double delta1 = rho0/Mix.pReducing->rhorbar(z1);
				double tau1 = Mix.pReducing->Tr(z1)/T0;
				double f1 = Mix.dphir_dDelta(tau1, delta1, z1);

				std::vector<double> z2 = z;
				z2[j] -= epsz;
				double delta2 = rho0/Mix.pReducing->rhorbar(z2);
				double tau2 = Mix.pReducing->Tr(z2)/T0;
				double f2 = Mix.dphir_dDelta(tau2, delta2, z2);

				double NUM = (f1-f2)/(2*epsz);
				CAPTURE(NUM);
				CAPTURE(ANA);
				CAPTURE(j);
				CHECK(fabs(NUM/ANA-1) < 1e-9);
			}
		}
		SECTION("d2nphir_dni_dxj__constT_V")
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					double ANA = Mix.d2nphir_dni_dxj__constT_V(tau0, delta0, z, i, j);
					
					std::vector<double> z1 = z;
					z1[j] += epsz;
					double delta1 = rho0/Mix.pReducing->rhorbar(z1);
					double tau1 = Mix.pReducing->Tr(z1)/T0;
					double f1 = Mix.dnphir_dni__constT_V_nj(tau1, delta1, z1, i);

					std::vector<double> z2 = z;
					z2[j] -= epsz;
					double delta2 = rho0/Mix.pReducing->rhorbar(z2);
					double tau2 = Mix.pReducing->Tr(z2)/T0;
					double f2 = Mix.dnphir_dni__constT_V_nj(tau2, delta2, z2, i);

					double NUM = (f1-f2)/(2*epsz);
					CAPTURE(NUM);
					CAPTURE(ANA);
					CAPTURE(j);
					CHECK(fabs(NUM/ANA-1) < 1e-8);
				}
			}
		}
		SECTION("dpdxj__constT_V_xi")
		{
			for (int j = 0; j < 2; j++)
			{
				double ANA = Mix.dpdxj__constT_V_xi(tau0, delta0, z, j);
				
				std::vector<double> z1 = z;
				z1[j] += epsz;
				double delta1 = rho0/Mix.pReducing->rhorbar(z1);
				double tau1 = Mix.pReducing->Tr(z1)/T0;
				double p1 = rho0*Mix.Rbar(z1)*T0*(1+delta1*Mix.dphir_dDelta(tau1,delta1,z1));

				std::vector<double> z2 = z;
				z2[j] -= epsz;
				double delta2 = rho0/Mix.pReducing->rhorbar(z2);
				double tau2 = Mix.pReducing->Tr(z2)/T0;
				double p2 = rho0*Mix.Rbar(z2)*T0*(1+delta2*Mix.dphir_dDelta(tau2,delta2,z2));

				double NUM = (p1-p2)/(2*epsz);
				CAPTURE(NUM);
				CAPTURE(ANA);
				CAPTURE(j);
				CHECK(fabs(NUM/ANA-1) < 1e-9);
			}
		}
		SECTION("ddelta_dxj__constT_V_xi")
		{
			for (int j = 0; j < 2; j++)
			{
				double ANA = Mix.ddelta_dxj__constT_V_xi(tau0, delta0, z, j);
				
				std::vector<double> z1 = z;
				z1[j] += epsz;
				double f1 = rho0/Mix.pReducing->rhorbar(z1);

				std::vector<double> z2 = z;
				z2[j] -= epsz;
				double f2 = rho0/Mix.pReducing->rhorbar(z2);

				double NUM = (f1-f2)/(2*epsz);
				CAPTURE(NUM);
				CAPTURE(ANA);
				CAPTURE(j);
				CHECK(fabs(NUM/ANA-1) < 1e-10);
			}
		}
		SECTION("ddtau_dxj__constT_V_xi")
		{
			for (int j = 0; j < 2; j++)
			{
				double ANA = Mix.dtau_dxj__constT_V_xi(tau0, delta0, z, j);
				
				std::vector<double> z1 = z;
				z1[j] += epsz;
				double f1 = Mix.pReducing->Tr(z1)/T0;

				std::vector<double> z2 = z;
				z2[j] -= epsz;
				double f2 = Mix.pReducing->Tr(z2)/T0;

				double NUM = (f1-f2)/(2*epsz);
				CAPTURE(NUM);
				CAPTURE(ANA);
				CAPTURE(j);
				CHECK(fabs(NUM/ANA-1) < 1e-10);
			}
		}
	}
#endif