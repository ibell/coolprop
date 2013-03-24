/*
The name of this file is a slight misnomer.  It includes the analysis
aqueous solutions as well as single-phase liquids.
*/


#ifndef BRINE_H
#define BRINE_H

#include <string>
//#include <boost/units/physical_dimensions/thermal_conductivity.hpp>
//#include <boost/units/systems/si/prefixes.hpp>
//#include <boost/units/systems/si/pressure.hpp>
//#include <boost/units/systems/si/energy.hpp>
//#include <boost/units/systems/si/time.hpp>
//#include <boost/units/systems/si/temperature.hpp>
//#include <boost/units/physical_dimensions/specific_energy.hpp>
//#include <boost/units/physical_dimensions/specific_heat_capacity.hpp>
//#include <boost/units/systems/si/io.hpp>
//
//using namespace boost::units;
//using namespace boost::units::si;
//
//namespace boost {
//  namespace units {
//    namespace si {
//
//		typedef make_system< pascal_base_unit, second_base_unit >::type system;
//
//      typedef unit< thermal_conductivity_dimension, si::system > thermal_conductivity;
//      static const thermal_conductivity watt_per_meter_kelvin;
//	  static const thermal_conductivity W_per_m_K;
//
//	  typedef unit< specific_energy_dimension, si::system > specific_energy;
//      static const specific_energy joule_per_kilogram;
//	  static const specific_energy J_per_kg;
//
//	  typedef unit< specific_heat_capacity_dimension, si::system > entropy;
//	  typedef unit< specific_heat_capacity_dimension, si::system > specific_heat_capacity;
//      static const specific_heat_capacity joule_per_kilogram_kelvin;
//	  static const specific_heat_capacity J_per_kg_K;
//	  
//    }
//  }
//}

double SecFluids(char Output, double T, double p, char * Ref);
int Brine(char * Mix, double T, double C, /*in --- out */double *Tfreeze, double *Tmax, double *rho, double *cp, double *k, double *visc, double *h, double *s);

//int Brine(char * Mix, 
//		  quantity<temperature> T, 
//		  double C, 
//		  /*in --- out */ 
//		  quantity<temperature> *Tfreeze, 
//		  quantity<temperature> *Tmax, 
//		  quantity<mass_density> *rho, 
//		  quantity<specific_heat_capacity> *cp, 
//		  quantity<thermal_conductivity> *k, 
//		  quantity<dynamic_viscosity> *visc, 
//		  quantity<specific_energy> *h, 
//		  quantity<entropy> *s);


#endif