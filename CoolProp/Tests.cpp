


#include "CoolPropDLL.h"
#include "CoolProp.h"

#ifndef DISABLE_CATCH
	#include "Tests.h"
	#define CATCH_CONFIG_RUNNER
	#include "Catch/catch.hpp"
	TEST_CASE( "Check reference state", "[reference_state]" ) 
	{
		SECTION("IIR")
		{
			set_reference_stateS("Propane","IIR");
			double h_EOS = PropsSI("H","T",273.15,"Q",0,"Propane");
			double s_EOS = PropsSI("S","T",273.15,"Q",0,"Propane");
			double h_target = 200000;
			double s_target = 1000;

			REQUIRE(abs(h_target-h_EOS) < 1e-6);
			REQUIRE(abs(s_target-s_EOS) < 1e-6);
		}
		SECTION("NBP")
		{
			set_reference_stateS("Propane","NBP");
			double h_EOS = PropsSI("H","P",101325,"Q",0,"Propane");
			double s_EOS = PropsSI("S","P",101325,"Q",0,"Propane");
			double h_target = 0;
			double s_target = 0;

			REQUIRE(abs(h_target-h_EOS) < 1e-6);
			REQUIRE(abs(s_target-s_EOS) < 1e-6);
		}
		SECTION("ASHRAE")
		{
			set_reference_stateS("Propane","ASHRAE");
			double h_EOS = PropsSI("H","T",233.15,"Q",0,"Propane");
			double s_EOS = PropsSI("S","T",233.15,"Q",0,"Propane");
			double h_target = 0;
			double s_target = 0;

			REQUIRE(abs(h_target-h_EOS) < 1e-6);
			REQUIRE(abs(s_target-s_EOS) < 1e-6);
		}
	}

	TEST_CASE( "Check units of fluid constants", "[fast]" ) 
	{
		SECTION("kSI")
		{
			set_standard_unit_system(UNIT_SYSTEM_KSI);

			SECTION("pcrit")
			{
				double p_Props1SI = Props1SI("R134a","pcrit");
				double p_Props1 = Props1("R134a","pcrit")*1000;
				double p_Props = Props("pcrit",'T',300,'Q',0,"R134a")*1000;
				double p_PropsSI = PropsSI("pcrit","T",300,"Q",0,"R134a");
				REQUIRE(fabs(p_Props1/p_Props-1) < 1e-6);
				REQUIRE(fabs(p_PropsSI/p_Props-1) < 1e-6);
				REQUIRE(fabs(p_PropsSI/p_Props1SI-1) < 1e-6);
			}
		}
		SECTION("SI")
		{
			set_standard_unit_system(UNIT_SYSTEM_SI);

			SECTION("pcrit")
			{
				double p_Props1SI = Props1("R134a","pcrit");
				double p_Props1 = Props1("R134a","pcrit");
				double p_Props = Props("pcrit",'T',300,'Q',0,"R134a");
				double p_PropsSI = PropsSI("pcrit","T",300,"Q",0,"R134a");
				REQUIRE(fabs(p_Props1/p_Props-1) < 1e-6);
				REQUIRE(fabs(p_PropsSI/p_Props-1) < 1e-6);
				REQUIRE(fabs(p_PropsSI/p_Props1SI-1) < 1e-6);
			}
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
		session.run();
	}


#endif