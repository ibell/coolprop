// In this file, the code for all the fluids is included verbatim.  This file is in turn included into FluidClass.cpp

#include <map>
#include <string>
#include <exception>
#include <vector>

#include "rapidjson_CoolProp.h"

#include "CPExceptions.h"
#include "CoolPropTools.h"
#include "CoolPropDLL.h"
#include "Helmholtz.h"
#include "Units.h"
#include "AllFluids.h"
#include "CriticalSplineConstants.h"

#include "pseudopurefluids/Air.cpp"
#include "pseudopurefluids/R507A.cpp"
#include "pseudopurefluids/R407C.cpp"
#include "pseudopurefluids/R410A.cpp"
#include "pseudopurefluids/R404A.cpp"
#include "pseudopurefluids/SES36.cpp"
#include "pseudopurefluids/R407F.cpp"

#include "purefluids/AceticAcid.cpp"
#include "purefluids/Argon.cpp"
#include "purefluids/Alkanes.cpp"
#include "purefluids/Benzene.cpp"
#include "purefluids/Butenes.cpp"
#include "purefluids/Cyclopentane.cpp"
#include "purefluids/Cyclopropane_Propyne.cpp"
#include "purefluids/Deuterium.cpp"
#include "purefluids/DMC.cpp"
#include "purefluids/Ethanol.cpp"
#include "purefluids/Ether.cpp"
#include "purefluids/Ethylene.cpp"
#include "purefluids/FAME.cpp"
#include "purefluids/Fluorine.cpp"
#include "purefluids/Helium.cpp"
#include "purefluids/HFE143m.cpp"
#include "purefluids/Hydrogen.cpp"
#include "purefluids/IndustrialFluids.cpp"
#include "purefluids/Methanol.cpp"
#include "purefluids/Neon.cpp"
#include "purefluids/Nitrogen.cpp"
#include "purefluids/Oxygen.cpp"
#include "purefluids/Propylene.cpp"
#include "purefluids/R12_R113.cpp"
#include "purefluids/R124.cpp"
#include "purefluids/R125.cpp"
#include "purefluids/R1234yf.cpp"
#include "purefluids/R1234ze.cpp"
#include "purefluids/R1233zd(E).cpp"
#include "purefluids/R134a.cpp"
#include "purefluids/R143A.cpp"
#include "purefluids/R161_Fluoroethane.cpp"
#include "purefluids/R22.cpp"
#include "purefluids/R227EA_R365MFC.cpp"
#include "purefluids/R23.cpp"
#include "purefluids/R236FA.cpp"
#include "purefluids/R236EA.cpp"
#include "purefluids/R290.cpp"
#include "purefluids/R32.cpp"
#include "purefluids/R717.cpp"
#include "purefluids/R744.cpp"
#include "purefluids/RC318_R21_R114_R13_R14.cpp"
#include "purefluids/Siloxanes.cpp"
#include "purefluids/Span_Short.cpp"
#include "purefluids/SulfurHexafluoride.cpp"
#include "purefluids/Undecane.cpp"
#include "purefluids/Water.cpp"
#include "purefluids/Xylene_EthylBenzene.cpp"

// ------------------------------
// FluidsContainer Implementation
// ------------------------------

FluidsContainer::FluidsContainer()
{
	
	// The pure fluids
	FluidsList.push_back(new WaterClass());
	FluidsList.push_back(new R134aClass());
	FluidsList.push_back(new HeliumClass());
	FluidsList.push_back(new OxygenClass());
	FluidsList.push_back(new HydrogenClass());
	FluidsList.push_back(new ParaHydrogenClass());
	FluidsList.push_back(new OrthoHydrogenClass());
	FluidsList.push_back(new ArgonClass());
	FluidsList.push_back(new R744Class());
	FluidsList.push_back(new NitrogenClass());
	FluidsList.push_back(new R290Class());
	FluidsList.push_back(new R717Class());
	FluidsList.push_back(new R1234yfClass());
	FluidsList.push_back(new R1234zeClass());
	FluidsList.push_back(new R32Class());
	FluidsList.push_back(new R22Class());
	FluidsList.push_back(new SES36Class());
	FluidsList.push_back(new EthyleneClass());
	FluidsList.push_back(new SulfurHexafluorideClass());
	FluidsList.push_back(new EthanolClass());
	FluidsList.push_back(new DimethylEtherClass());	
	FluidsList.push_back(new DimethylCarbonateClass());	
	FluidsList.push_back(new R143AClass());	
	FluidsList.push_back(new R23Class());
	FluidsList.push_back(new nDodecaneClass());	
	FluidsList.push_back(new PropyleneClass());	
	FluidsList.push_back(new CyclopentaneClass());	
	FluidsList.push_back(new R236FAClass());
	FluidsList.push_back(new R236EAClass());
	FluidsList.push_back(new R227EAClass());
	FluidsList.push_back(new R365MFCClass());
	FluidsList.push_back(new R161Class());
	FluidsList.push_back(new HFE143mClass());
	FluidsList.push_back(new BenzeneClass());
	FluidsList.push_back(new UndecaneClass());
	FluidsList.push_back(new R125Class());
	FluidsList.push_back(new CycloPropaneClass());
	FluidsList.push_back(new NeonClass());
	FluidsList.push_back(new R124Class());
	FluidsList.push_back(new PropyneClass());
	FluidsList.push_back(new FluorineClass());
	FluidsList.push_back(new MethanolClass());
	FluidsList.push_back(new RC318Class());
	FluidsList.push_back(new R21Class());
	FluidsList.push_back(new R114Class());
	FluidsList.push_back(new R13Class());
	FluidsList.push_back(new R14Class());
	FluidsList.push_back(new R12Class());
	FluidsList.push_back(new R113Class());
	FluidsList.push_back(new R1234zeZClass());
	FluidsList.push_back(new R1233zdEClass());
	FluidsList.push_back(new AceticAcidClass());

	// The industrial fluids
	FluidsList.push_back(new R245faClass());
	FluidsList.push_back(new R41Class());
	FluidsList.push_back(new CarbonMonoxideClass());
	FluidsList.push_back(new CarbonylSulfideClass());
	FluidsList.push_back(new DecaneClass());
	FluidsList.push_back(new HydrogenSulfideClass());
	FluidsList.push_back(new IsopentaneClass());
	FluidsList.push_back(new NeopentaneClass());
	FluidsList.push_back(new IsohexaneClass());
	FluidsList.push_back(new KryptonClass());
	FluidsList.push_back(new NonaneClass());
	FluidsList.push_back(new TolueneClass());
	FluidsList.push_back(new XenonClass());
	FluidsList.push_back(new R116Class());
	FluidsList.push_back(new AcetoneClass());
	FluidsList.push_back(new NitrousOxideClass());
	FluidsList.push_back(new SulfurDioxideClass());
	FluidsList.push_back(new R141bClass());
	FluidsList.push_back(new R142bClass());
	FluidsList.push_back(new R218Class());
	
	FluidsList.push_back(new MethaneClass());
	FluidsList.push_back(new EthaneClass());
	FluidsList.push_back(new nButaneClass());
	FluidsList.push_back(new IsoButaneClass());
	// Span Non-Polar
	FluidsList.push_back(new nPentaneClass());
	FluidsList.push_back(new nHexaneClass());
	FluidsList.push_back(new nHeptaneClass());
	FluidsList.push_back(new nOctaneClass());
	FluidsList.push_back(new CyclohexaneClass());
	// Span Polar
	FluidsList.push_back(new R152AClass());
	FluidsList.push_back(new R123Class());
	FluidsList.push_back(new R11Class());

	// The Siloxanes
	FluidsList.push_back(new OctamethyltrisiloxaneClass()); //MDM
	FluidsList.push_back(new DecamethyltetrasiloxaneClass()); //MD2M
	FluidsList.push_back(new DodecamethylpentasiloxaneClass()); //MD3M
	FluidsList.push_back(new DodecamethylcyclohexasiloxaneClass()); //D6
	FluidsList.push_back(new HexamethyldisiloxaneClass());//MM
	FluidsList.push_back(new TetradecamethylhexasiloxaneClass()); //MD4M
	FluidsList.push_back(new OctamethylcyclotetrasiloxaneClass()); //D4
	FluidsList.push_back(new DecamethylcyclopentasiloxaneClass()); //D5

	// The butenes
	FluidsList.push_back(new OneButeneClass());
	FluidsList.push_back(new IsoButeneClass());
	FluidsList.push_back(new Cis2ButeneClass());
	FluidsList.push_back(new Trans2ButeneClass());
	
	// The methyl ester components of biodiesel
	FluidsList.push_back(new MethylPalmitateClass());
	FluidsList.push_back(new MethylStearateClass());
	FluidsList.push_back(new MethylOleateClass());
	FluidsList.push_back(new MethylLinoleateClass());
	FluidsList.push_back(new MethylLinolenateClass());

	// Xylene isomers and EthylBenzene
	FluidsList.push_back(new oXyleneClass());
	FluidsList.push_back(new mXyleneClass());
	FluidsList.push_back(new pXyleneClass());
	FluidsList.push_back(new EthylBenzeneClass());

	// Deuterium and spin isomers
	FluidsList.push_back(new DeuteriumClass());
	FluidsList.push_back(new ParaDeuteriumClass());
	FluidsList.push_back(new OrthoDeuteriumClass());

	// The pseudo-pure fluids
	FluidsList.push_back(new AirClass());
	FluidsList.push_back(new R404AClass());
	FluidsList.push_back(new R410AClass());
	FluidsList.push_back(new R407CClass());
	FluidsList.push_back(new R507AClass());
	FluidsList.push_back(new R407FClass());

	// Includes the C++ JSON code for all the fluids as the variable JSON_code
	#include "JSON_code.h"

	// Includes the C++ JSON code for the CAS number lookup as the variable JSON_cas
	#include "JSON_CAS.h"

	JSON.Parse<0>(JSON_code);
	JSON_CAS.Parse<0>(JSON_cas);

	// Build the map of fluid names mapping to pointers to the Fluid class instances
	for (std::vector<Fluid*>::iterator it = FluidsList.begin(); it != FluidsList.end(); it++)
	{
		//// Load all the parameters related to the Critical point spline
		set_critical_spline_constants((*it));
		// Call the post_load routine
		(*it)->post_load(JSON, JSON_CAS);
		// Load up entry in map
		fluid_name_map[(*it)->get_name()] = *it;

		std::string name = (*it)->get_name();
		std::string ucasename = upper(name);
		
		if (!((*it)->isAlias(ucasename)) && ucasename.compare(name))
		{
			(*it)->add_alias(ucasename);
		}
	}
}

// Destructor
FluidsContainer::~FluidsContainer()
{
	while (!FluidsList.empty())
	{
		delete FluidsList.back();
		FluidsList.pop_back();
	}
}

bool FluidsContainer::add_REFPROP_fluid(std::string FluidName, std::vector<double> xmol)
{
	// Some fluids from REFPROP that are not included in CoolProp due to not having a Helmholtz energy EOS
	Fluid * pREFPROPFluid = new REFPROPFluidClass(FluidName,xmol);
	FluidsList.push_back(pREFPROPFluid);
	// Add entry to the fluid name map
	fluid_name_map.insert(std::pair<std::string,Fluid*>(FluidName,pREFPROPFluid));
	return true;
}

Fluid * FluidsContainer::get_fluid(long iFluid)
{
	if (iFluid > -1)
	{
		return FluidsList[iFluid];
	}
	else
	{
		return NULL;
	}
}
Fluid * FluidsContainer::get_fluid(std::string name)
{
	std::map<std::string,Fluid*>::iterator it;
	// Try to find using the map if Fluid name is provided
	it = fluid_name_map.find(name);
	// If it is found the iterator will not be equal to end
	if (it != fluid_name_map.end() )
	{
		// Return a pointer to the class
		return (*it).second;
	}

	// Wasn't found, now we need to check for an alias
	for (std::vector<Fluid*>::iterator it = FluidsList.begin(); it != FluidsList.end(); it++)
	{
		if ( (*it)->isAlias(name) )
		{
			return *it;
		}
	}
	return NULL;
}

long FluidsContainer::get_fluid_index(Fluid* pFluid)
{
	int i = 0;
	// Iterate to find the 0-based index of the fluid
	for (std::vector<Fluid*>::iterator it = FluidsList.begin(); it != FluidsList.end(); it++)
	{
		if ((*it) == pFluid)
		{
			return i;
		}
		i++;
	}
	return -1;
}

std::string FluidsContainer::FluidList()
{
	// Return a std::string with the list of fluids
	std::string FL;
	for (std::vector<Fluid*>::iterator it = FluidsList.begin(); it != FluidsList.end(); it++)
	{
		FL+=(*it)->get_name();
		FL+=",";
	}
	//Remove the tailing comma
	FL = FL.substr (0,FL.length()-1);
	return FL;
}

#ifndef DISABLE_CATCH
#include "Catch/catch.hpp"

TEST_CASE((char*)"Check ancillary curves for pure and pseudo-pure fluids","[slow],[ancillary]")
{
	FluidsContainer Fluids = FluidsContainer();

	SECTION((char*)"Saturated Liquid Pressure Ancillary")
	{
		for (std::vector<Fluid*>::const_iterator it = Fluids.FluidsList.begin(); it != Fluids.FluidsList.end(); it++)
		{
			std::string name = (*it)->get_name();
			double Tmin = (*it)->limits.Tmin;
			double Tmax = (*it)->crit.T-1;
			int N = 5;
			for (double T = Tmin; T <= Tmax; T += (Tmax-Tmin)/(N-1))
			{
				CAPTURE(name);
				double pL,pV,rhoL,rhoV;
				CHECK_NOTHROW((*it)->saturation_T(T,false,pL,pV,rhoL,rhoV));
				double p_EOS = pL;
				double p_ancillary = (*it)->pure() ? (*it)->psat(T): (*it)->psatL(T);
				CAPTURE(p_EOS);
				CAPTURE(p_ancillary);
				CHECK(fabs(p_EOS/p_ancillary -1) <= 2e-2);
			}
		}
	}
	SECTION((char*)"Saturated Vapor Pressure Ancillary")
	{
		for (std::vector<Fluid*>::const_iterator it = Fluids.FluidsList.begin(); it != Fluids.FluidsList.end(); it++)
		{
			std::string name = (*it)->get_name();
			double Tmin = (*it)->limits.Tmin;
			double Tmax = (*it)->crit.T-1;
			int N = 5;
			for (double T = Tmin; T <= Tmax; T += (Tmax-Tmin)/(N-1))
			{
				double pL,pV,rhoL,rhoV;
				(*it)->saturation_T(T,false,pL,pV,rhoL,rhoV);
				double p_EOS = pV;
				double p_ancillary = (*it)->pure() ? (*it)->psat(T): (*it)->psatV(T);
				
				CAPTURE(name);
				CAPTURE(p_EOS);
				CAPTURE(p_ancillary);
				CHECK(fabs(p_EOS/p_ancillary -1) <= 2e-2);
			}
		}
	}

	SECTION((char*)"Saturated Liquid Density Ancillary")
	{
		for (std::vector<Fluid*>::const_iterator it = Fluids.FluidsList.begin(); it != Fluids.FluidsList.end(); it++)
		{
			std::string name = (*it)->get_name();
			double Tmin = (*it)->limits.Tmin;
			double Tmax = (*it)->crit.T-1;
			int N = 5;
			for (double T = Tmin; T <= Tmax; T += (Tmax-Tmin)/(N-1))
			{
				double pL,pV,rhoL,rhoV;
				(*it)->saturation_T(T,false,pL,pV,rhoL,rhoV);
				double rho_EOS = rhoL;
				double rho_ancillary = (*it)->rhosatL(T);
				CAPTURE(name);
				CAPTURE(rho_EOS);
				CAPTURE(rho_ancillary);
				CHECK(fabs(rho_EOS/rho_ancillary -1) <= 2e-2);
			}
		}
	}

	SECTION((char*)"Saturated Vapor Density Ancillary")
	{
		for (std::vector<Fluid*>::const_iterator it = Fluids.FluidsList.begin(); it != Fluids.FluidsList.end(); it++)
		{
			std::string name = (*it)->get_name();
			double Tmin = (*it)->limits.Tmin;
			double Tmax = (*it)->crit.T-1;
			int N = 5;
			for (double T = Tmin; T <= Tmax; T += (Tmax-Tmin)/(N-1))
			{
				double pL,pV,rhoL,rhoV;
				(*it)->saturation_T(T,false,pL,pV,rhoL,rhoV);
				double rho_EOS = rhoV;
				double rho_ancillary = (*it)->rhosatV(T);
				CAPTURE(name);
				CAPTURE(rho_EOS);
				CAPTURE(rho_ancillary);
				CHECK(fabs(rho_EOS/rho_ancillary -1) <= 2e-2);
			}
		}
	}
}

TEST_CASE((char*)"Fluid parameter checks not requiring saturation","[fast]")
{
	FluidsContainer Fluids = FluidsContainer();

	SECTION((char*)"Check Tmin > Ttriple")
	{
		for (std::vector<Fluid*>::const_iterator it = Fluids.FluidsList.begin(); it != Fluids.FluidsList.end(); it++)
		{
			REQUIRE((*it)->params.Ttriple <= (*it)->limits.Tmin);
		}
	}
	if (REFPROPFluidClass::refpropSupported())
	{
		SECTION((char*)"Check pcrit matches REFPROP")
		{
			for (std::vector<Fluid*>::const_iterator it = Fluids.FluidsList.begin(); it != Fluids.FluidsList.end(); it++)
			{
				std::string RPName = std::string("REFPROP-")+get_fluid_param_string((*it)->get_name(), "REFPROPName");
				if (RPName.compare("REFPROP-N/A"))
				{
					// Skip if not in REFPROP
					double pcrit_CP = (*it)->crit.p.Pa;
					double pcrit_RP = PropsSI((char*)"pcrit",(char*)"T",300,(char*)"D",1e-10,(char*)RPName.c_str());
					CAPTURE(pcrit_CP);
					CAPTURE(pcrit_RP);
					CAPTURE(RPName);
					CHECK(fabs(pcrit_RP/pcrit_CP-1) < 0.01);
				}
			}
		}
	}
}

TEST_CASE((char*)"Fluid parameter checks requiring saturation","[slow]")
{
	FluidsContainer Fluids = FluidsContainer();

	SECTION((char*)"Check ptriple")
	{
		for (std::vector<Fluid*>::const_iterator it = Fluids.FluidsList.begin(); it != Fluids.FluidsList.end(); it++)
		{
			std::string name = (*it)->get_name();
			double pL,pV,rhoL,rhoV;
			(*it)->saturation_T((*it)->limits.Tmin,false,pL,pV,rhoL,rhoV);
			double ptriple_EOS = pV;
			CAPTURE(name);
			CAPTURE(ptriple_EOS);
			INFO(name);
			REQUIRE((*it)->params.ptriple == ptriple_EOS);
		}
	}
	SECTION((char*)"Check accentric factor")
	{
		for (std::vector<Fluid*>::const_iterator it = Fluids.FluidsList.begin(); it != Fluids.FluidsList.end(); it++)
		{
			if ((*it)->limits.Tmin < 0.7*(*it)->crit.T)
			{
				std::string name = (*it)->get_name();
				double pL,pV,rhoL,rhoV;
				(*it)->saturation_T((*it)->crit.T*0.7,false,pL,pV,rhoL,rhoV);
				double accentric_EOS = -log10(pV/(*it)->crit.p.Pa)-1;
				double accentric_Fluid = (*it)->params.accentricfactor;
				CAPTURE(name);
				CAPTURE(accentric_EOS);
				CAPTURE(accentric_Fluid);
				INFO(format("accentric factor should be %0.7g",accentric_EOS));
				REQUIRE(fabs(accentric_Fluid/accentric_EOS-1) < 1e-2);
			}
		}
	}
}
TEST_CASE((char*)"Saturation consistency checks", (char*)"[slow],[consistency]" )
{
	FluidsContainer Fluids = FluidsContainer();
	SECTION((char*)"saturation_T")
	{
		for (std::vector<Fluid*>::const_iterator it = Fluids.FluidsList.begin(); it != Fluids.FluidsList.end(); it++)
		{
			double Tt = (*it)->limits.Tmin;
			double Tc = (*it)->crit.T;
			double N = 30;
			for (double T = Tt; T<Tc; T+=(Tc-Tt)/(N-1))
			{
				double p, psatV, rhoL, rhoV;
				std::string name = (*it)->get_name();
				CAPTURE(name);
				CAPTURE(T);
				CAPTURE(Tc);
				CHECK_NOTHROW((*it)->saturation_T(T, false, p, psatV, rhoL, rhoV));
			}
		}
	}

	SECTION((char*)"TL->pL->TL")
	{
		for (std::vector<Fluid*>::const_iterator it = Fluids.FluidsList.begin(); it != Fluids.FluidsList.end(); it++)
		{
			if (!(*it)->pure()){continue;}
			double Tt = (*it)->limits.Tmin;
			double Tc = (*it)->crit.T;
			double pc = (*it)->crit.p.Pa;
			double N = 50;
			for (double T = Tt; T<Tc; T+=(Tc-0.0001-Tt)/(N-1))
			{
				double psatL, psatV, TsatL, TsatV,rhoL,rhoV;
				std::string name = (*it)->get_name();
				CAPTURE(name);
				CAPTURE(T);
				CAPTURE(Tc);
				CHECK_NOTHROW((*it)->saturation_T(T, false, psatL, psatV, rhoL, rhoV));
				CAPTURE(pc);
				CHECK_NOTHROW((*it)->saturation_p(psatL, false, TsatL, TsatV, rhoL, rhoV));
				CAPTURE(TsatL);
				INFO(name);
				CHECK(fabs(T/TsatL-1) < 1e-3);
			}
		}
	}
	SECTION((char*)"TV->pV->TV")
	{
		for (std::vector<Fluid*>::const_iterator it = Fluids.FluidsList.begin(); it != Fluids.FluidsList.end(); it++)
		{
			if (!(*it)->pure()){continue;}
			double Tt = (*it)->limits.Tmin;
			double Tc = (*it)->crit.T;
			double pc = (*it)->crit.p.Pa;
			double N = 50;
			for (double T = Tt; T<Tc; T+=(Tc-0.0001-Tt)/(N-1))
			{
				double psatL, psatV, TsatL, TsatV, rhoL, rhoV;
				std::string name = (*it)->get_name();
				CAPTURE(name);
				CAPTURE(T);
				CAPTURE(Tc);
				CHECK_NOTHROW((*it)->saturation_T(T, false, psatL, psatV, rhoL, rhoV));
				CAPTURE(psatV);
				CAPTURE(pc);
				CHECK_NOTHROW((*it)->saturation_p(psatV, false, TsatL, TsatV, rhoL, rhoV));
				CAPTURE(TsatV);
				INFO(name);
				CHECK(fabs(T/TsatV-1) < 1e-3);
			}
		}
	}
}
#endif