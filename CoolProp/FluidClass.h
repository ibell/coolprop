#include "Helmholtz.h"
#include <list>
#include <string>
#include <exception>
#include "CPExceptions.h"

#ifndef FLUIDCLASS_H
#define FLUIDCLASS_H

struct OtherParameters
{
	double molemass, Ttriple, accentricfactor, R_u;
};
struct CriticalStruct
{
	double rho, T, p, v;
};
struct FluidLimits
{
	double Tmin, Tmax, pmax, rhomax;	
};

struct SatLUTStruct
{
	std::vector<double> T,rhoL,rhoV,p,hL,hV,tau,logp;
	int N;
	bool built;
};

// Fluid is the abstract base class that is employed by all the other fluids
class Fluid
{
	protected:
		std::string name;
		std::list <std::string> aliases;
		std::list <phi_BC*> phirlist;
		std::list <phi_BC*> phi0list;
		std::string EOSReference;
		std::string TransportReference;
		bool isPure;
		SatLUTStruct SatLUT;
		double _Dekker_Tsat(double x_min, double x_max, double eps, double p, double Q);
    public:
		// Constructor/destructor - this gets called before the derived type constructor
		// This allows to set default values that can be overwritten by the derived class
		Fluid(){
			params.R_u=8.314472; 
			isPure=true; 
			SatLUT.N = 200;
			SatLUT.built = false;
			preduce = &crit;
		};
		virtual ~Fluid();

		// Fluid-specific parameters
		struct CriticalStruct crit;
		struct FluidLimits limits;
		struct OtherParameters params;
		struct CriticalStruct * preduce; // The point that is used to reduce the T and rho for EOS
		struct CriticalStruct reduce; // The point that is used to reduce the T and rho for EOS

		// Member Access functions
		std::string get_name(){return name;};
		bool pure(){return isPure;};
		double R(){return params.R_u/params.molemass;};

		// These MUST be implemented by derived class
        virtual double conductivity_Trho(double tau, double delta) = 0;
        virtual double viscosity_Trho(double tau, double delta) = 0;
		
		// These Helmholtz energy terms are provided by the base class
		double phir(double tau, double delta);
		double dphir_dDelta(double tau, double delta);
		double d2phir_dDelta2(double tau, double delta);
		double d2phir_dDelta_dTau(double tau, double delta);
		double d3phir_dDelta2_dTau(double tau, double delta);
		double dphir_dTau(double tau, double delta);
		double d2phir_dTau2(double tau, double delta);
		double phi0(double tau, double delta);
		double dphi0_dDelta(double tau, double delta);
		double d2phi0_dDelta2(double tau, double delta);
		double d2phi0_dDelta_dTau(double tau, double delta);
		double dphi0_dTau(double tau, double delta);
		double d2phi0_dTau2(double tau, double delta);

		// These thermodynamic properties as a function of temperature and density are
		// provided by the base class
		double pressure_Trho(double T, double rho);
		double enthalpy_Trho(double T, double rho);
		double entropy_Trho(double T, double rho);
		double internal_energy_Trho(double T, double rho);
		double speed_sound_Trho(double T, double rho);
		double specific_heat_p_Trho(double T, double rho);
		double specific_heat_v_Trho(double T, double rho);
		double gibbs_Trho(double T, double rho);
		double density_Tp(double T, double p, double rho_guess);

		// Optional ancillary functions can be overloaded, will throw a NotImplementedError
		// to be caught by calling function if not implemented
		virtual double psat(double T){throw NotImplementedError(); return 1e99;};
		virtual double psatL(double T){throw NotImplementedError(); return 1e99;};
		virtual double psatV(double T){throw NotImplementedError(); return 1e99;};
		virtual double rhosatL(double T){throw NotImplementedError(); return 1e99;};
		virtual double rhosatV(double T){throw NotImplementedError(); return 1e99;};

		// Saturation properties
		void saturation(double T, bool UseLUT, double *psatLout, double *psatVout, double *rhoLout, double *rhoVout);
		void rhosatPure(double T, double *rhoLout, double *rhoVout, double *pout);
		void BuildSaturationLUT();
		double ApplySaturationLUT(char *OutPropName,char *InPropName,double InPropVal);
		
		// Tsat functions - without UseLUT provided, assumed to be false
		double Tsat(double p, double Q, double T_guess);
		double Tsat(double p, double Q, double T_guess, bool UseLUT);
		
		bool isAlias(std::string name);
};

#endif