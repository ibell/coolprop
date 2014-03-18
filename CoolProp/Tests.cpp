


#include "CoolPropDLL.h"
#include "CoolProp.h"

#ifndef DISABLE_CATCH
	#include "Tests.h"
	#define CATCH_CONFIG_RUNNER
	#include "Catch/catch.hpp"
	TEST_CASE( "Check reference state",  "[reference_state]" ) 
	{
		SECTION( "IIR")
		{
			set_reference_stateS("Propane","IIR");
			double h_EOS = PropsSI("H","T",273.15,"Q",0,"Propane");
			double s_EOS = PropsSI("S","T",273.15,"Q",0,"Propane");
			double h_target = 200000;
			double s_target = 1000;

			REQUIRE(abs(h_target-h_EOS) < 1e-6);
			REQUIRE(abs(s_target-s_EOS) < 1e-6);
		}
		SECTION( "NBP")
		{
			set_reference_stateS( "Propane", "NBP");
			double h_EOS = PropsSI("H","P",101325,"Q",0,"Propane");
			double s_EOS = PropsSI("S","P",101325,"Q",0,"Propane");
			double h_target = 0;
			double s_target = 0;

			REQUIRE(abs(h_target-h_EOS) < 1e-6);
			REQUIRE(abs(s_target-s_EOS) < 1e-6);
		}
		SECTION("ASHRAE")
		{
			set_reference_stateS( "Propane", "ASHRAE");
			double h_EOS = PropsSI( "H", "T",233.15, "Q",0, "Propane");
			double s_EOS = PropsSI( "S", "T",233.15, "Q",0, "Propane");
			double h_target = 0;
			double s_target = 0;

			REQUIRE(abs(h_target-h_EOS) < 1e-6);
			REQUIRE(abs(s_target-s_EOS) < 1e-6);
		}
	}

	TEST_CASE("Check units of fluid constants",  "[fast]" ) 
	{
		SECTION( "kSI")
		{
			set_standard_unit_system(UNIT_SYSTEM_KSI);

			SECTION( "pcrit")
			{
				double p_Props1SI = Props1SI( "R134a", "pcrit");
				double p_Props1 = Props1( "R134a", "pcrit")*1000;
				double p_Props = Props( "pcrit",'T',300,'Q',0, "R134a")*1000;
				double p_PropsSI = PropsSI( "pcrit", "T",300, "Q",0, "R134a");
                double p_IProps = IProps(iPcrit,iT,0,iP,0,get_Fluid_index("R134a"))*1000;
                double p_IPropsSI = IPropsSI(iPcrit,iT,0,iP,0,get_Fluid_index("R134a"));
				REQUIRE(fabs(p_Props1/p_Props-1) < 1e-6);
				REQUIRE(fabs(p_PropsSI/p_Props-1) < 1e-6);
				REQUIRE(fabs(p_PropsSI/p_Props1SI-1) < 1e-6);
                REQUIRE(fabs(p_IProps/p_Props1SI-1) < 1e-6);
                REQUIRE(fabs(p_IPropsSI/p_Props1SI-1) < 1e-6);
			}
		}
		SECTION( "SI")
		{
			set_standard_unit_system(UNIT_SYSTEM_SI);

			SECTION( "pcrit")
			{
				double p_Props1SI = Props1( "R134a", "pcrit");
				double p_Props1 = Props1( "R134a", "pcrit");
				double p_Props = Props( "pcrit",'T',300,'Q',0, "R134a");
				double p_PropsSI = PropsSI( "pcrit", "T",300, "Q",0, "R134a");
                double p_IProps = IProps(iPcrit,iT,0,iP,0,get_Fluid_index("R134a"));
                double p_IPropsSI = IPropsSI(iPcrit,iT,0,iP,0,get_Fluid_index("R134a"));
				REQUIRE(fabs(p_Props1/p_Props-1) < 1e-6);
				REQUIRE(fabs(p_PropsSI/p_Props-1) < 1e-6);
				REQUIRE(fabs(p_PropsSI/p_Props1SI-1) < 1e-6);
                REQUIRE(fabs(p_IProps/p_Props1SI-1) < 1e-6);
                REQUIRE(fabs(p_IProps/p_Props1SI-1) < 1e-6);
			}
            set_standard_unit_system(UNIT_SYSTEM_KSI);
		}
	}

    TEST_CASE("Normal boiling point of water",  "[fast]" )
	{
        SECTION( "kSI")
		{
            int us = get_standard_unit_system();
			set_standard_unit_system(UNIT_SYSTEM_KSI);
			SECTION( "pcrit")
			{
				double T_Props = Props("T","P",101.325,"Q",0,"Water");
                double T_PropsSI = PropsSI("T","P",101325,"Q",0,"Water");
                double T_IProps = IProps(iT,iP,101.325,iQ,0,get_Fluid_index("Water"));
                double T_IPropsSI = IPropsSI(iT,iP,101325,iQ,0,get_Fluid_index("Water"));
				double T_nbp = 99.98+273.15;
                REQUIRE(fabs(T_Props-T_nbp) < 0.01);
				REQUIRE(fabs(T_PropsSI-T_Props) < 0.01);
				REQUIRE(fabs(T_IProps-T_PropsSI) < 0.01);
                REQUIRE(fabs(T_IPropsSI-T_IProps) < 0.01);
			}
            set_standard_unit_system(us);
		}
		SECTION( "SI")
		{
            int us = get_standard_unit_system();
			set_standard_unit_system(UNIT_SYSTEM_SI);
			SECTION( "pcrit")
			{
				double T_Props = Props("T","P",101325,"Q",0,"Water");
                double T_PropsSI = PropsSI("T","P",101325,"Q",0,"Water");
                double T_IProps = IProps(iT,iP,101325,iQ,0,get_Fluid_index("Water"));
                double T_IPropsSI = IPropsSI(iT,iP,101325,iQ,0,get_Fluid_index("Water"));
				double T_nbp = 99.98+273.15;
                REQUIRE(fabs(T_Props-T_nbp) < 0.01);
				REQUIRE(fabs(T_PropsSI-T_Props) < 0.01);
				REQUIRE(fabs(T_IProps-T_PropsSI) < 0.01);
                REQUIRE(fabs(T_IPropsSI-T_IProps) < 0.01);
			}
            set_standard_unit_system(us);
		}
	}

	static Catch::Session session; // There must be exactly once instance

	int run_fast_tests()
	{
		Catch::ConfigData &config = session.configData();
		config.testsOrTags.clear();
		config.testsOrTags.push_back("[fast]");
		session.useConfigData(config);
		return session.run();
	}

	int run_not_slow_tests()
	{
		Catch::ConfigData &config = session.configData();
		config.testsOrTags.clear();
		config.testsOrTags.push_back("~[slow]");
		session.useConfigData(config);

		time_t t1, t2;
		t1 = clock();
		session.run();
		t2 = clock();
		printf("Elapsed time for not slow tests: %g s",(double)(t2-t1)/CLOCKS_PER_SEC);

		return 1;
	}

	int run_user_defined_tests(const std::vector<std::string> & tests_or_tags)
	{
		Catch::ConfigData &config = session.configData();
		config.testsOrTags.clear();
		for (unsigned int i = 0; i < tests_or_tags.size(); i++)
		{
			config.testsOrTags.push_back(tests_or_tags[i]);
		}
		session.useConfigData(config);

		time_t t1, t2;
		t1 = clock();
		session.run();
		t2 = clock();
		printf("Elapsed time for user defined tests: %g s",(double)(t2-t1)/CLOCKS_PER_SEC);

		return 1;
	}

	void run_tests()
	{
        Catch::ConfigData &config = session.configData();
		config.testsOrTags.clear();
        session.useConfigData(config);
		session.run();
	}


#endif