#include "Helmholtz.h"
#include <list>
#include <map>
#include <string>
#include <exception>
#include "CPExceptions.h"
#include "PropMacros.h"
#include <Eigen/Dense>

#ifndef FLUIDCLASS_H
#define FLUIDCLASS_H

struct OtherParameters
{
	double molemass, Ttriple, ptriple, accentricfactor, R_u;
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

struct OnePhaseLUTStruct
{
	std::vector< std::vector<double> > hmat,rhomat,cpmat,cp0mat, smat,cvmat,umat,viscmat,kmat,pmat,dpdTmat;
	std::vector<double> Tvec,pvec;
	double Tmin,Tmax,pmin,pmax;
	int nT,np;
	bool built,forcebuild;
};

/// Fluid is the abstract base class that is employed by all the other fluids
class Fluid
{
	protected:
		std::string name; /// The name of the fluid
		std::string REFPROPname; /// The REFPROP-compliant name if REFPROP-"name" is not a compatible fluid name.  If not included, "name" is assumed to be a valid name for REFPROP
		std::list <std::string> aliases; /// A list of aliases of names for the Fluid, each element is a std::string instance
		
		std::string EOSReference; /// A std::string that contains a reference for thermo properties for the fluid
		std::string TransportReference; /// A std::string that contains a reference for the transport properties of the fluid
		bool isPure; /// True if it is a pure fluid, false otherqwise
		SatLUTStruct SatLUT; /// The private Saturation lookup structure
		OnePhaseLUTStruct LUT; /// The private single-phase lookup structure
		
		/// Obtain a guess value for the density of the fluid for a given set of temperature and pressure
		/// @param T Temperature [K]
		/// @param p Pressure [kPa(abs)]
		double _get_rho_guess(double T, double p); 

		std::vector<std::vector<double> >* _get_LUT_ptr(std::string Prop);
    public:

		std::list <phi_BC*> phirlist; /// A list of instances of the phi_BC classes for the residual Helmholtz energy contribution
		std::list <phi_BC*> phi0list; /// A list of instances of the phi_BC classes for the ideal-gas Helmholtz energy contribution

		/// Constructor for the Fluid class.  This is an abstract base class that 
		/// is not meant to be instantiated directly.  Rather it should be subclassed
		/// to implement new fluids.  Examples of this behavior are the R134aClass or the 
		/// R290Class classes.
		Fluid(){
			params.R_u=8.314472; /// Default molar gas constant [kJ/kmol/K]
			isPure=true;  /// true if fluid is one-component pure fluid, false otherwise
			SatLUT.N = 200; /// Number of points in saturation lookup table [-]
			SatLUT.built = false; /// reset the saturation LUT built flag
			LUT.built = false; /// reset the LUT built flag
				
			preduce = &crit; /// pointer to the reducing parameters

			// Default limits for the Single-phase LUT
			// Can be overwritten by calling set_1phase_LUT_params()
			LUT.nT = 200;
			LUT.np = 200;
			LUT.Tmax = limits.Tmax-0.1;
			LUT.Tmin = limits.Tmin+0.1;
			LUT.pmin = 0.0001;
			LUT.pmax = limits.pmax-0.1;
		};
		virtual ~Fluid();

		// Some post-loading things happen here
		void post_load(void);

		// Fluid-specific parameters
		struct CriticalStruct crit;
		struct FluidLimits limits;
		struct OtherParameters params;
		struct CriticalStruct * preduce; /// A pointer to the point that is used to reduce the T and rho for EOS
		struct CriticalStruct reduce; /// The point that is used to reduce the T and rho for EOS

		// Member Access functions
		/// Returns a std::string with the name of the fluid
		std::string get_name(){return name;};
		/// Returns a char* with the name of the fluid
		char * get_namec(){return (char *)name.c_str();};
		std::string get_REFPROPname(){return REFPROPname;};
		std::string get_EOSReference(){return EOSReference;};
		std::string get_TransportReference(){return TransportReference;};
		/// Returns true if the fluid is pure, false if pseudo-pure or a mixture
		bool pure(){return isPure;};
		/// Returns the mass-specific gas constant for the fluid [kJ/kg/K]
		double R(){return params.R_u/params.molemass;};

		// These MUST be implemented by derived class
        virtual double conductivity_Trho(double T, double rho);
		virtual double viscosity_Trho(double T, double rho);
		
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
		double specific_heat_p_ideal_Trho(double T);
		double specific_heat_v_Trho(double T, double rho);
		double gibbs_Trho(double T, double rho);
		double dpdT_Trho(double T,double rho);
		double density_Tp(double T, double p);
		double density_Tp(double T, double p, double rho_guess);
		std::string phase_Tp(double T, double p);

		// Optional ancillary functions can be overloaded, will throw a NotImplementedError
		// to be caught by calling function if not implemented
		virtual double psat(double T){
			throw NotImplementedError(std::string("psat not implemented for this fluid")); return _HUGE;
		};
		virtual double psatL(double T){
			throw NotImplementedError(std::string("psatL not implemented for this fluid")); return _HUGE;
		};
		virtual double psatV(double T){
			throw NotImplementedError(std::string("psatV not implemented for this fluid")); return _HUGE;
		};
		virtual double rhosatL(double T){
			throw NotImplementedError(std::string("rhosatL not implemented for this fluid")); return _HUGE;
		};
		virtual double rhosatV(double T){
			throw NotImplementedError(std::string("rhosatV not implemented for this fluid")); return _HUGE;
		};

		Eigen::Vector2d ConformalTemperature(Fluid *InterestFluid, Fluid *ReferenceFluid,double T, double rho, std::string *errstring);
		// Extended corresponding states functions for fluids that do not have their own high-accuracy
		// transport property implementation
		virtual void ECSParams(double *e_k, double *sigma){
			throw NotImplementedError("ECSParams not implemented for this fluid");
		};
		/// This function is optional, and returns a NotImplementedError if the derived class does not implement it.
		/// Hopefully the calling function catches the error
		///@param rhor The reduced density where rhor = rho/rhoc
		virtual double ECS_psi_viscosity(double rhor){
			throw NotImplementedError("ECS_psi_viscosity not implemented for this fluid");
		};
		/// This function is optional, and returns a NotImplementedError if the derived class does not implement it.
		/// Hopefully the calling function catches the error
		///@param rhor The reduced density where rhor = rho/rhoc
		virtual double ECS_chi_conductivity(double rhor){
			throw NotImplementedError("ECS_chi_conductivity not implemented for this fluid");
		};
		/// This function is optional, and returns a NotImplementedError if the derived class does not implement it.
		/// Hopefully the calling function catches the error
		///@param T The temperature
		virtual double ECS_f_int(double T){
			throw NotImplementedError("ECS_f_int not implemented for this fluid");
		};
		/// This function is optional, and returns a NotImplementedError if the derived class does not implement it.
		/// Hopefully the calling function catches the error
		///@param T Temperature [K]
		///@param rho Density [kg/m3]
		///@param ReferenceFluid A pointer to an instance of a Fluid class that is used as the reference fluid for the ECS model
		double viscosity_ECS_Trho(double T, double rho, Fluid * ReferenceFluid);

		/// This function is optional, and returns a NotImplementedError if the derived class does not implement it.
		/// Hopefully the calling function catches the error
		///@param T Temperature [K]
		///@param rho Density [kg/m3]
		///@param ReferenceFluid A pointer to an instance of a Fluid class that is used as the reference fluid for the ECS model
		double conductivity_ECS_Trho(double T, double rho, Fluid * ReferenceFluid);

		/// The background viscosity for the fluid.  The viscosity minus the critical contribution minus the dilute-gas contribution
		/// This function is optional, and returns a NotImplementedError if the derived class does not implement it.
		/// Hopefully the calling function catches the error
		/// @see viscosity_Trho
		/// @param T Temperature [K]
		/// @param rho Density [kg/m3]
		virtual double conductivity_background(double T, double rho){
			throw NotImplementedError(std::string("conductivity_background not implemented for this fluid")); return _HUGE;
		};

		virtual double conductivity_critical(double T, double rho);

		/// This function returns the dilute portion of the viscosity
		/// @param T Temperature [K]
		/// @param e_k epsilon/k_B for the fluid [K]
		/// @param sigma [nm]
		/// @returns mu Dilute-gas viscosity in the limit of zero density [Pa-s]
		double viscosity_dilute(double T, double e_k, double sigma);

		/// The residual viscosity for the fluid.  The sum of the background and critical contributions. 
		/// This function is optional, and returns a NotImplementedError if the derived class does not implement it.
		/// Hopefully the calling function catches the error
		/// @see viscosity_Trho
		/// @param T Temperature [K]
		/// @param rho Density [kg/m3]
		virtual double viscosity_residual(double T, double rho){
			throw NotImplementedError(std::string("viscosity_residual not implemented for this fluid")); return _HUGE;
		};
		/// The background viscosity for the fluid.  The viscosity minus the critical contribution minus the dilute-gas contribution
		/// This function is optional, and returns a NotImplementedError if the derived class does not implement it.
		/// Hopefully the calling function catches the error
		/// @see viscosity_Trho
		/// @param T Temperature [K]
		/// @param rho Density [kg/m3]
		virtual double viscosity_background(double T, double rho){
			throw NotImplementedError(std::string("viscosity_background not implemented for this fluid")); return _HUGE;
		};

		virtual double surface_tension_T(double T);

		/// Saturation pressure and saturated liquid and vapor densities as a function of the temperature.
		/// @param UseLUT If True, use the Saturation Lookup tables, otherwise use the EOS and the equal gibbs function and equal pressure criterion to determine saturation state
		/// @param T Temperature [K]
		/// @param psatLout Saturated liquid pressure [kPa(abs)]
		/// @param psatVout Saturated vapor pressure [kPa(abs)]
		/// @param rhoLout Saturated liquid pressure [kg/m3]
		/// @param rhoVout Saturated vapor pressure [kg/m3]
		void saturation(double T, bool UseLUT, double *psatLout, double *psatVout, double *rhoLout, double *rhoVout);
		
		/// NB: Only valid for pure fluids - no pseudo-pure or mixtures.
		/// Get the saturated liquid, vapor densities and the saturated pressure
		/// @param T Temperature [K]]
		/// @param pout Saturated pressure [kPa(abs)]
		/// @param rhoLout Saturated liquid pressure [kg/m3]
		/// @param rhoVout Saturated vapor pressure [kg/m3]
		void rhosatPure(double T, double *rhoLout, double *rhoVout, double *pout);

		/// NB: Only valid for pure fluids - no pseudo-pure or mixtures.
		/// Get the saturated liquid, vapor densities and the saturated pressure using the method from Akasaka given by
		/// Ryo Akasaka, "A reliable and useful Method to Determine the Saturation State from Helmholtz Energy Equations of State", Journal of Thermal Science and Technology, v.3 n3 2008, 442-451
		/// @param T Temperature [K]]
		/// @param pout Saturated pressure [kPa(abs)]
		/// @param rhoLout Saturated liquid pressure [kg/m3]
		/// @param rhoVout Saturated vapor pressure [kg/m3]
		void rhosatPure_Akasaka(double T, double *rhoLout, double *rhoVout, double *pout);

		/// Tries to build the Saturation Lookup tables for Temperature, pressure, densities and enthalpies.  
		///  If the LUT are already built, don't rebuild them, but if the parameter SatLUT.force=True, rebuild.
		void BuildSaturationLUT(void);

		/// Apply the Saturation lookup tables to a given set of inputs and output
		/// @param OutPropName Single-char corresponding to the output of interest. @see Props
		/// @param InPropName Single-char corresponding to the input of interest. @see Props
		/// @param InPropVal Value of the input of interest. @see Props
		double ApplySaturationLUT(char *OutPropName,char *InPropName,double InPropVal);
		
		/// The saturation temperature of the fluid for a given pressure and quality. If the 
		/// fluid is pure, the saturated vapor and saturated liquid temperatures are the same.  
		/// If not, give the "correct" saturation temperature.
		/// @param p Saturation pressure [kPa(abs)]
		/// @param Q Vapor quality in the range [0,1]
		/// @param T_guess  Not currently used due to the use of the Brent solver. Used to be the guess temperature
		double Tsat(double p, double Q, double T_guess);

		/// The saturation temperature of the fluid for a given pressure and quality. If the 
		/// fluid is pure, the saturated vapor and saturated liquid temperatures are the same.  
		/// If not, give the "correct" saturation temperature.
		/// @param p Saturation pressure [kPa(abs)]
		/// @param Q Vapor quality in the range [0,1]
		/// @param T_guess  Not currently used due to the use of the Brent solver. Used to be the guess temperature
		/// @param UseLUT Use the Lookup table (True), or not (False)
		double Tsat(double p, double Q, double T_guess, bool UseLUT);

		/// Build the single-phase LUT in the range[Tmin,Tmax]x[pmin,pmax] if not already built
		void BuildLookupTable();

		/// Set the single-phase LUT range, including the number of steps and the range for each variable
		/// @param nT The number of points to use for the temperature
		/// @param np The number of points to use for the pressure
		/// @param Tmin The minimum temperature [K]
		/// @param Tmax The maximum temperature [K]
		/// @param pmin The minimum pressure [kPa(abs)]
		/// @param pmax The maximum pressure [kPa(abs)]
		void set_1phase_LUT_params(int nT, int np, double Tmin, double Tmax, double pmin, double pmax);

		/// Get the single-phase LUT range, including the number of steps and the range for each variable
		/// @param nT The number of points to use for the temperature
		/// @param np The number of points to use for the pressure
		/// @param Tmin The minimum temperature [K]
		/// @param Tmax The maximum temperature [K]
		/// @param pmin The minimum pressure [kPa(abs)]
		/// @param pmax The maximum pressure [kPa(abs)]
		void get_1phase_LUT_params(int *nT, int *np, double *Tmin, double *Tmax, double *pmin, double *pmax);
		
		/// Use the single-phase lookup table to determine	the value for the desired property as a function of temperature and pressure.
		/// The definition of the keys for Prop are given by Props
		/// @see Props
		/// @param Prop Single character corresponding to the property of interest
		/// @param T Temperature [K]
		/// @param p Temperature [kPa(abs)]
		/// @returns Value corresponding to the key Prop [varies]
		double LookupValue_TP(std::string Prop, double T, double p);

		/// Use the single-phase lookup table to determine	the value for the desired property as a function of temperature and density.
		/// The definition of the keys for Prop are given by
		/// @see Props
		/// @param Prop Single character corresponding to the property of interest
		/// @param T Temperature [K]
		/// @param rho Density [kg/m3]
		double LookupValue_Trho(std::string Prop, double T, double rho);
		
		/// Returns true if the given name is an alias of the Fluid name.  (Case-sensitive!!)
		/// @param name The given name
		bool isAlias(std::string name);
};

	/// This class contains pointers to all the classes of the fluids that can be used in addition  
	/// to some convenience functions
	class FluidsContainer
	{
	private:
		std::map<std::string,Fluid*> fluid_name_map; ///< maps fluid names to pointers to the fluid
		std::list <Fluid*> FluidsList; ///< A list of pointers to the instances of the fluids
	public:
		/// Constructor for the FluidsContainer class
		/// @see FluidsContainer
		FluidsContainer();

		/// Destructor for the FluidsContainer class.  Deletes each fluid in tern
		/// @see FluidsContainer
		~FluidsContainer();

		/// Accessor.  Throws a NotImplementedError if the fluid given by name is not found.  Also searches aliases
		/// @param name Fluid to be searched for
		Fluid * get_fluid(std::string name);

		/// Returns a std::string of a comma-separated list of the CoolProp names of all the fluids that are loaded.
		std::string FluidList();
	};

#endif
