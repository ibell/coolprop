#ifndef TTSE_H
#define TTSE_H

#include "FluidClass.h"

class Fluid;

class TTSETwoPhaseTableClass
{
protected:
	unsigned int N;
	Fluid *pFluid;
	double dh,dp;
	double pratio,logpratio;
public:
	/// Default Instantiator
	TTSETwoPhaseTableClass(){};
	/// Instantiator
	/// @param pFluid Pointer to an instance of a Fluid class
	/// @param Q Quality [kg/kg], in [0,1]
	TTSETwoPhaseTableClass(Fluid *pFluid, double Q);
	/// Destructor
	~TTSETwoPhaseTableClass(){};

	/// Set the size of the Two-Phase table
	/// @param N Number of elements in arrays
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
	/// @param other TTSETwoPhaseTableClass for the other phase boundary (liquid for the vapor, or vice versa)
	double build(double pmin, double pmax, TTSETwoPhaseTableClass *other = NULL);
	
	/// Evaluate a property in the two-phase region using the TTSE method with p as input
	/// @param iParam Index of desired output
	/// @param p Pressure (absolute) [kPa]
 	double evaluate(long iParam, double p);

	/// Evaluate a property in the two-phase region using the TTSE method with p as input
	/// @param T Temperature [K]
 	double evaluate_T(double T);

	/// Evaluate the derivative of a property along the saturation curve using the TTSE method
	/// @param iParam Index of desired output
	/// @param p Pressure (absolute) [kPa]
	double evaluate_sat_derivative(long iParam, double p);

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
	unsigned int Nh, Np;
	Fluid *pFluid;
	double pratio,logpratio,logpmin;
	int jpcrit_floor, jpcrit_ceil;
	std::vector<double> T_Trho, rho_Trho;
	std::vector<std::vector<int> > i_Trho, j_Trho;

	void ph_Trho(int i, int j, double T, double rho, double *pout, double *hout);
public:
	TTSESinglePhaseTableClass();
	TTSESinglePhaseTableClass(Fluid *pFluid);
	TTSETwoPhaseTableClass *SatL, *SatV;
	void set_size(unsigned int Nh = 100, unsigned int Np = 100);
	
	std::string root_path;
	std::vector<double> TL,SL,DL,TV,SV,DV;

	double hmin,hmax,pmin,pmax;
	bool enable_writing_tables_to_files;

	// Variables with h, p
	std::vector<std::vector<double> > T,dTdh,dTdp,d2Tdh2,d2Tdp2,d2Tdhdp;
	std::vector<std::vector<double> > rho,drhodh,drhodp,d2rhodh2,d2rhodp2,d2rhodhdp;
	std::vector<std::vector<double> > s,dsdh,dsdp,d2sdh2,d2sdp2,d2sdhdp;
	std::vector<std::vector<double> > mu,dmudh,dmudp,d2mudh2,d2mudp2,d2mudhdp;
	std::vector<std::vector<double> > k,dkdh,dkdp,d2kdh2,d2kdp2,d2kdhdp;
	std::vector<double> h, p;

	/// Indices of the last indices within the two-phase region for pressures 
	/// from pmin to pmax.  If supercritical, iL and iV entries both -1
	std::vector<int> iL, iV;

	bool read_all_from_file(std::string root_path);
	void write_all_to_file(std::string root_path = std::string());
	void matrix_to_file(std::string fName, std::vector< std::vector<double> > *A);
	void vector_to_file(std::string fName, std::vector<double> *A);
	void vector_from_file(std::string fName, int N, std::vector<double> *A);
	void matrix_from_file(std::string fName, std::vector< std::vector<double> > *A);

	void update_saturation_boundary_indices();
	
	/// Update the map from the values of T,rho to the indices i,j
	void update_Trho_map();



	/// Build the tables
	/// @param hmin Minimum enthalpy [kJ/kg]
	/// @param hmax Maximum enthalpy [kJ/kg]
	/// @param pmin Minimum pressure [kJ/kg]
	/// @param pmax Maximum pressure [kJ/kg]
	/// @param TTSESatL Saturated liquid TTSE LUT
	/// @param TTSESatV Saturated vapor TTSE LUT
	double build(double hmin, double hmax, double pmin, double pmax, TTSETwoPhaseTableClass *SatL = NULL, TTSETwoPhaseTableClass *SatV = NULL);
	
	/// Evaluate a property in the single-phase region
	/// @param iParam Index of desired output
	/// @param p Pressure (absolute) [kPa]
	/// @param h Enthalpy [kJ/kg]
	double evaluate(long iParam, double p, double h);

	/// Evaluate the TTSE using P,S or P,T
	double evaluate_one_other_input(long iInput1, double Param1, long iOther, double Other);
	
	/// Evaluate the TTSE using T,D or XXX to give P,H
	void evaluate_two_other_inputs(long iInput1, double Input1, long iInput2, double Input2, double *pout, double *hout);

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

	/// Find the nearest neighbor density and temperature if they exist to speed up calcs for the calculation of T(p,h) and rho(p,h)
	/// @param i Index in h
	/// @param j Index in p
	/// @param T0 Temperature
	/// @param rho0 Density
	void nearest_neighbor(int i, int j, double *T0, double *rho0);

	/// Find the nearest neighbor indices that have good values if i,j are not good
	/// @param i Index in h
	/// @param j Index in p
	void nearest_good_neighbor(int *i, int *j);

	void get_partner(int i, int j, double p, double h, int *ipartner, int *jpartner);

	/// Evaluate the first partial derivative
	/// @param iOF Index in numerator
	/// @param iWRT Index of denominator
	/// @param iCONSTANT Index of property held constant
	/// @param p Pressure [kPa]
	/// @param h Enthalpy [kJ/kg]
	double evaluate_first_derivative(long iOF, long iWRT, long iCONSTANT, double p, double h);

};



#endif