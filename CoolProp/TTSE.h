#ifndef TTSE_H
#define TTSE_H

class TTSESinglePhaseTable
{
protected:
	unsigned int Nrow, Ncol;
	Fluid *pFluid;
	CoolPropStateClass *pCPS;
	double dh,dp;

public:
	
	TTSESinglePhaseTable(Fluid *pFluid, int Nrow=100, int Ncol=100);
	~TTSESinglePhaseTable(){delete pCPS;};

	double hmin,hmax,pmin,pmax;

	// Variables with h, p
	std::vector<std::vector<double>> T,dTdh,dTdp,d2Tdh2,d2Tdp2,d2Tdhdp;
	std::vector<std::vector<double>> rho,drhodh,drhodp,d2rhodh2,d2rhodp2,d2rhodhdp;
	std::vector<std::vector<double>> s,dsdh,dsdp,d2sdh2,d2sdp2,d2sdhdp;
	std::vector<double> h,p;

	/// Build the tables
	/// @param hmin Minimum enthalpy [kJ/kg]
	/// @param hmax Maximum enthalpy [kJ/kg]
	/// @param pmin Minimum pressure [kJ/kg]
	/// @param pmax Maximum pressure [kJ/kg]
	/// @param logp Flag whether to use log(p) rather than p [true/false]
	double build(double hmin, double hmax, double pmin, double pmax, bool logp = false);
	
	/// Evaluate a property in the single-phase region
	/// @param iParam Index of desired output
	/// @param h Enthalpy [kJ/kg]
	/// @param p Pressure (absolute) [kPa]
	double evaluate(long iParam, double h, double p);

	/// Randomly select a point within the range, and evaluate the property using TTSE and the EOS
	/// @param iParam Index of desired output
	/// @param N Number of runs to do
	/// @param h std::vector of enthalpies
	/// @param p std::vector of enthalpies
	/// @param EOS std::vector of values from Equation of State
	/// @param TTSE std::vector of values from TTSE
	///
	/// Note: h,p,EOS, TTSE should be empty std::vector passed by reference
	double check_randomly(long iParam, unsigned int N, std::vector<double> *h, std::vector<double> *p, std::vector<double> *EOS, std::vector<double> *TTSE);
};

#endif