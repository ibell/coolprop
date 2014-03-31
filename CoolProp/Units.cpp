#include "Units.h"
#include "GlobalConstants.h"
#include "CoolProp.h"

double convert_from_unit_system_to_SI(long iInput, double value, int old_system)
{
	if (get_debug_level() > 8)
	{
		std::cout << format("%s:%d: convert_from_unit_system_to_SI(i=%d,value=%g,old_system=%d)\n",__FILE__,__LINE__,iInput,value,old_system).c_str();
	}
	if (old_system == UNIT_SYSTEM_SI)
	{
		return value;
	}

	switch (iInput)
	{
	case iP:
	case iPtriple:
	case iPcrit:
	case iPsat:
	case iC:
	case iC0:
	case iS:
	case iG:
	case iO:
	case iH:
	case iU:
	case iL:
	case iDERdp_dT__rho:
	case iDERdp_drho__T:
	case iDERdh_dT__rho:
	case iDERdh_drho__T:
		switch (old_system)
		{
		case UNIT_SYSTEM_KSI:
			return value*1000.0;
			break;
		default:
			throw ValueError(format("Unit System [%d] is undefined",old_system).c_str());
			break;
		}
		break;
	case iDERdrho_dh__p:
	case iDERdrho_dp__h:
	case iDERdrho_smoothed_dh:
	case iDERdrho_smoothed_dp:
	case iDERdrhodh_constp_smoothed:
	case iDERdrhodp_consth_smoothed:
	case iDERIsothermalCompressibility:
		switch (old_system)
		{
		case UNIT_SYSTEM_KSI:
			return value/1000.0;
			break;
		default:
			throw ValueError(format("Unit System [%d] is undefined",old_system).c_str());
			break;
		}
		break;
	case iD:
	case iA:
	case iQ:
	case iV:
	case iT:
	case iTcrit:
	case iTtriple: 
	case iTreduce:
	case iPrandtl:
	case iTmin:
	case iTmax:
	case iTfreeze:
	case iPHASE_LIQUID:
	case iPHASE_GAS:
	case iPHASE_SUPERCRITICAL:
	case iPHASE_TWOPHASE:
	case iODP:
	case iGWP100:
	case iMM:
	case iRhocrit:
	case iRhoreduce: 
	case iAccentric:
	case iI:
	case iCritSplineT:
	case iDERdh_dp__rho:
	case iDERdh_dp__v:
	case iDERZ:
	case iDERdZ_dDelta:
	case iDERdZ_dTau:
	case iDERB:
	case iDERdB_dT:
	case iDERC:
	case iDERdC_dT:
	case iDERphir:
	case iDERdphir_dTau:
	case iDERdphir_dDelta:
	case iDERd2phir_dTau2:
	case iDERd2phir_dDelta2:
	case iDERd2phir_dDelta_dTau:
	case iDERd3phir_dDelta3:
	case iDERd3phir_dDelta2_dTau:
	case iDERd3phir_dDelta_dTau2:
	case iDERd3phir_dTau3:
	case iDERphi0:
	case iDERdphi0_dTau:
	case iDERd2phi0_dTau2:
	case iDERdphi0_dDelta:
	case iDERd2phi0_dDelta2:
	case iDERd2phi0_dDelta_dTau:
	case iDERd3phi0_dTau3:
	case iDERdrho_dT__p:
	case iDERrho_smoothed:
		return value;
		break;
	default:
		throw ValueError(format("index [%d] is invalid in convert_from_unit_system_to_SI",iInput).c_str());
		break;
	}
	return _HUGE;
}

double convert_from_SI_to_unit_system(long iInput, double value, int new_system)
{
	if (get_debug_level() > 8)
	{
		std::cout << format("%s:%d: convert_from_SI_to_unit_system(%d,%g,%d)\n",__FILE__,__LINE__,iInput,value,new_system).c_str();
	}

	if (new_system == UNIT_SYSTEM_SI)
	{
		return value;
	}

	switch (iInput)
	{
	case iP:
	case iPtriple:
	case iPcrit:
	case iPsat:
	case iC:
	case iC0:
	case iS:
	case iG:
	case iO:
	case iH:
	case iU:
	case iL:
	case iDERdp_dT__rho:
	case iDERdp_drho__T:
	case iDERdh_dT__rho:
	case iDERdh_drho__T:
		switch (new_system)
		{
		case UNIT_SYSTEM_KSI:
			return value/1000.0;
			break;
		default:
			throw ValueError(format("Unit System [%d] is undefined",new_system).c_str());
			break;
		}
		break;
	case iDERdrho_dh__p:
	case iDERdrho_dp__h:
	case iDERdrho_smoothed_dh:
	case iDERdrho_smoothed_dp:
	case iDERdrhodh_constp_smoothed:
	case iDERdrhodp_consth_smoothed:
	case iDERIsothermalCompressibility:
		switch (new_system)
		{
		case UNIT_SYSTEM_KSI:
			return value*1000.0;
			break;
		default:
			throw ValueError(format("Unit System [%d] is undefined",new_system).c_str());
			break;
		}
		break;
	case iD:
	case iQ:
	case iA:
	case iV:
	case iT:
	case iTmin:
	case iTmax:
	case iTfreeze:
	case iTcrit:
	case iTtriple: 
	case iTreduce:
	case iPHASE_LIQUID:
	case iPHASE_GAS:
	case iPHASE_SUPERCRITICAL:
	case iPHASE_TWOPHASE:
	case iODP:
	case iGWP100:
	case iMM:
	case iRhocrit:
	case iRhoreduce: 
	case iAccentric:
	case iI:
	case iCritSplineT:
	case iDERdh_dp__rho:
	case iDERdh_dp__v:
	case iDERZ:
	case iDERdZ_dDelta:
	case iDERdZ_dTau:
	case iDERB:
	case iDERdB_dT:
	case iDERC:
	case iDERdC_dT:
	case iDERphir:
	case iDERdphir_dTau:
	case iDERdphir_dDelta:
	case iDERd2phir_dTau2:
	case iDERd2phir_dDelta2:
	case iDERd2phir_dDelta_dTau:
	case iDERd3phir_dDelta3:
	case iDERd3phir_dDelta2_dTau:
	case iDERd3phir_dDelta_dTau2:
	case iDERd3phir_dTau3:
	case iDERphi0:
	case iDERdphi0_dTau:
	case iDERd2phi0_dTau2:
	case iDERdphi0_dDelta:
	case iDERd2phi0_dDelta2:
	case iDERd2phi0_dDelta_dTau:
	case iDERd3phi0_dTau3:
	case iDERdrho_dT__p:
	case iDERrho_smoothed:
    case iDERdh_dp__T:
		return value;
		break;
	default:
		throw ValueError(format("index [%d] is invalid in convert_from_SI_to_unit_system",iInput).c_str());
		break;
	}
	return _HUGE;
}

double convert_from_SI_to_unit_system(std::string input, double value, std::string new_system)
{
	int iSys = -1;
	if (!new_system.compare("kSI")) iSys = UNIT_SYSTEM_KSI;
	if (!new_system.compare( "SI")) iSys = UNIT_SYSTEM_SI;
	double val = convert_from_SI_to_unit_system(get_param_index(input), value, iSys);
	return val;
}
double convert_from_unit_system_to_SI(std::string input, double value, std::string old_system)
{
	int iSys = -1;
	if (!old_system.compare("kSI")) iSys = UNIT_SYSTEM_KSI;
	if (!old_system.compare( "SI")) iSys = UNIT_SYSTEM_SI;
	double val = convert_from_unit_system_to_SI(get_param_index(input), value, iSys);
	return val;
}
double convert_from_SI_to_unit_system(char *input, double value, char *new_system)
{
	return convert_from_SI_to_unit_system(std::string(input), value, std::string(new_system));
}
double convert_from_unit_system_to_SI(char *input, double value, char *old_system)
{
	return convert_from_unit_system_to_SI(std::string(input), value, std::string(old_system));
}

double conversion_factor(std::string num)
{
	double factor = 1.0;
	
	// Parse the first entry
	if (num[0] != '1')
	{
		std::string key;
		key.resize(1);
		key[0] = num[0];
		factor *= convert_from_SI_to_unit_system(get_param_index(key),1,get_standard_unit_system());
	}
	
	for (int i = 1; i < (int)num.size(); i += 2)
	{
		char op = num[i];

		std::string key;
		key.resize(1);
		key[0] = num[i+1];
		double term = convert_from_SI_to_unit_system(get_param_index(key),1,get_standard_unit_system());
		if (op == '*')
		{
			factor *= term;
		}
		else
		{
			factor /= term;
		}

	}
	return factor;
}
