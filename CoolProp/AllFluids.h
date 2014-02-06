
#ifndef ALLFLUIDS_H
#define ALLFLUIDS_H

#include "rapidjson_CoolProp.h"

#include "FluidClass.h"

/// This class contains pointers to all the classes of the fluids that can be used in addition  
/// to some convenience functions
class FluidsContainer
{
private:
	rapidjson::Document JSON, JSON_CAS;
	std::map<std::string,long> fluid_index_map; ///< maps fluid names to 0-based index
	std::map<std::string,Fluid*> fluid_name_map; ///< maps fluid names to pointers to the fluid
	
public:

	std::vector <Fluid*> FluidsList; ///< A list of pointers to the instances of the fluids
	
	/// Constructor for the FluidsContainer class
	/// @see FluidsContainer
	FluidsContainer();

	/// Destructor for the FluidsContainer class.  Deletes each fluid in tern
	/// @see FluidsContainer
	~FluidsContainer();

	/// Accessor.  Throws a NotImplementedError if the fluid given by name is not found.  Also searches aliases
	/// @param name Fluid to be searched for
	Fluid * get_fluid(std::string name);

	/// Accessor.  Throws a NotImplementedError if the fluid index is invalid
	/// @param name Fluid index to be used
	Fluid * get_fluid(long iFluid);

	/// Add a REFPROP fluid to the container
	/// @param FluidName 
	/// @param xmol std::vector of mole fractions - must add to 1
	bool add_REFPROP_fluid(std::string FluidName, std::vector<double> xmol);
	
	/// Get the index of the fluid name
	/// @param pFluid pointer to Fluid instance
	long get_fluid_index(Fluid* pFluid);

	/// Returns a std::string of a comma-separated list of the CoolProp names of all the fluids that are loaded.
	std::string FluidList();
};

#endif
