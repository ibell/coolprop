
#include "CoolPropC.h"

#include <stdbool.h>
#include <string>
#include <CoolProp/CoolPropDLL.h>

double CONVENTION props1_(char * FluidName, char * Output)
{
  //return Props1(std::string(FluidName), std::string(Output));
  return Props1(FluidName, Output);
}

double CONVENTION props_(char * Output,char Name1, double Prop1, char Name2, double Prop2, char * FluidName)
{
  //return Props(std::string(Output), Name1, Prop1, Name2, Prop2, std::string(FluidName));
  return Props(Output, Name1, Prop1, Name2, Prop2, FluidName);
}

double CONVENTION derivterms_(char * Term, double T, double rho, char * FluidName)
{
  //return DerivTerms(std::string(Term), T, rho, std::string(FluidName));
  return DerivTerms(Term, T, rho, FluidName);
}

int CONVENTION set_reference_states_(char * FluidName, char * reference_state)
{
  //return set_reference_stateS(std::string(FluidName), std::string(reference_state));
  return set_reference_stateS(FluidName, reference_state);
}
//int set_reference_stateS(std::string FluidName, std::string reference_state);

/// Enable the TTSE
bool CONVENTION enable_ttse_lut_(char *FluidName)
{
  return enable_TTSE_LUT(FluidName);
}
//bool CONVENTION enable_TTSE_LUT(char *FluidName);

/// Check if TTSE is enabled
bool CONVENTION isenabled_ttse_lut_(char *FluidName)
{
  return isenabled_TTSE_LUT(FluidName);
}
//bool CONVENTION isenabled_TTSE_LUT(char *FluidName);

/// Disable the TTSE
bool CONVENTION disable_ttse_lut_(char *FluidName)
{
  return disable_TTSE_LUT(FluidName);
}
//bool CONVENTION disable_TTSE_LUT(char *FluidName);

/// Set the TTSE mode (normal or bicubic)
int CONVENTION set_ttse_mode_(char *FluidName, char * Value)
{
  return set_TTSE_mode(FluidName, Value);
}
//int CONVENTION set_TTSE_mode(char *FluidName, char * Value);