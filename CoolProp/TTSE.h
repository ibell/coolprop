


#ifndef TTSE_H
#define TTSE_H

class TTSETwoPhaseTableClass
{
protected:
	unsigned int N;
	Fluid *pFluid;
	double dh,dp;

public:
	/// Default Instantiator
	TTSETwoPhaseTableClass(){};
	/// Instantiator
	/// @param pFluid Pointer to an instance of a Fluid class
	/// @param N Number of elements in arrays
	/// @param Q Quality [kg/kg], in [0,1]
	TTSETwoPhaseTableClass(Fluid *pFluid, double Q);
	~TTSETwoPhaseTableClass(){};

	void set_size(unsigned int N);

	double pmin,pmax,Q,logpmin,logpmax;

	// Variables with h, p
	std::vector<double> T,dTdp,d2Tdp2;
	std::vector<double> rho,drhodp,d2rhodp2,logrho;
	std::vector<double> s,dsdp,d2sdp2;
	std::vector<double> h,dhdp,d2hdp2;
	std::vector<double> p,logp;

	/// Build the tables along the saturation curves
	/// @param pmin Minimum pressure [kJ/kg]
	/// @param pmax Maximum pressure [kJ/kg]
	double build(double pmin, double pmax);
	
	/// Evaluate a property in the two-phase region using the TTSE method and
	/// @param iParam Index of desired output
	/// @param p Pressure (absolute) [kPa]
 	double evaluate(long iParam, double p);

	/// Randomly evaluate a property in the two-phase region using the TTSE method
	/// @param iParam Index of desired output
	/// @param N Number of runs to do
 	double evaluate_randomly(long iParam, unsigned int N);

	/// Randomly select points within the range, and evaluate the property using TTSE and the EOS
	/// @param iParam Index of desired output
	/// @param N Number of runs to do
	/// @param p std::vector of pressures
	/// @param EOS std::vector of values from Equation of State
	/// @param TTSE std::vector of values from TTSE
	///
	/// Note: p,EOS, TTSE should be empty std::vector passed by reference
	double check_randomly(long iParam, unsigned int N, std::vector<double> *p, std::vector<double> *EOS, std::vector<double> *TTSE);
};

class TTSESinglePhaseTableClass
{
protected:
	unsigned int Nrow, Ncol;
	Fluid *pFluid;
	double dh,dp;

public:
	TTSESinglePhaseTableClass(){};
	TTSESinglePhaseTableClass(Fluid *pFluid);
	void set_size(unsigned int Nrow=100, unsigned int Ncol=100);

	double hmin,hmax,pmin,pmax;

	// Variables with h, p
	std::vector<std::vector<double> > T,dTdh,dTdp,d2Tdh2,d2Tdp2,d2Tdhdp;
	std::vector<std::vector<double> > rho,drhodh,drhodp,d2rhodh2,d2rhodp2,d2rhodhdp;
	std::vector<std::vector<double> > s,dsdh,dsdp,d2sdh2,d2sdp2,d2sdhdp;
	std::vector<double> h,p;

	/// Build the tables
	/// @param hmin Minimum enthalpy [kJ/kg]
	/// @param hmax Maximum enthalpy [kJ/kg]
	/// @param pmin Minimum pressure [kJ/kg]
	/// @param pmax Maximum pressure [kJ/kg]
	/// @param logp Flag whether to use log(p) rather than p [true/false]
	double build(double hmin, double hmax, double pmin, double pmax, bool logp = false);

	/// Build the tables, but using the saturation TTSE LUT to determine where you are
	/// @param hmin Minimum enthalpy [kJ/kg]
	/// @param hmax Maximum enthalpy [kJ/kg]
	/// @param pmin Minimum pressure [kJ/kg]
	/// @param pmax Maximum pressure [kJ/kg]
	/// @param TTSESatL Saturated liquid TTSE LUT
	/// @param TTSESatV Saturated vapor TTSE LUT
	double build_TTSESat(double hmin, double hmax, double pmin, double pmax, TTSETwoPhaseTableClass *SatL, TTSETwoPhaseTableClass *SatV);
	
	/// Evaluate a property in the single-phase region
	/// @param iParam Index of desired output
	/// @param h Enthalpy [kJ/kg]
	/// @param p Pressure (absolute) [kPa]
	double evaluate(long iParam, double h, double p);

	/// Randomly evaluate a property in the single phase region using the TTSE method
	/// @param iParam Index of desired output
	/// @param N Number of runs to do
 	double evaluate_randomly(long iParam, unsigned int N);

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

	/// Write a representation of the ph surface to file with O in each "good" spot and "X" in each "bad" one or two-phase
	void write_dotdrawing_tofile(char fName[]);
};



#endif