


#include "CoolPropDLL.h"
#include "CoolProp.h"

#ifndef DISABLE_CATCH
	#include "Tests.h"
	#define CATCH_CONFIG_RUNNER
	#include "Catch/catch.hpp"
	TEST_CASE((char*)"Check reference state", (char*)"[reference_state]" ) 
	{
		SECTION((char*)"IIR")
		{
			set_reference_stateS("Propane","IIR");
			double h_EOS = PropsSI("H","T",273.15,"Q",0,"Propane");
			double s_EOS = PropsSI("S","T",273.15,"Q",0,"Propane");
			double h_target = 200000;
			double s_target = 1000;

			REQUIRE(abs(h_target-h_EOS) < 1e-6);
			REQUIRE(abs(s_target-s_EOS) < 1e-6);
		}
		SECTION((char*)"NBP")
		{
			set_reference_stateS((char*)"Propane",(char*)"NBP");
			double h_EOS = PropsSI("H","P",101325,"Q",0,"Propane");
			double s_EOS = PropsSI("S","P",101325,"Q",0,"Propane");
			double h_target = 0;
			double s_target = 0;

			REQUIRE(abs(h_target-h_EOS) < 1e-6);
			REQUIRE(abs(s_target-s_EOS) < 1e-6);
		}
		SECTION("ASHRAE")
		{
			set_reference_stateS((char*)"Propane",(char*)"ASHRAE");
			double h_EOS = PropsSI((char*)"H",(char*)"T",233.15,(char*)"Q",0,(char*)"Propane");
			double s_EOS = PropsSI((char*)"S",(char*)"T",233.15,(char*)"Q",0,(char*)"Propane");
			double h_target = 0;
			double s_target = 0;

			REQUIRE(abs(h_target-h_EOS) < 1e-6);
			REQUIRE(abs(s_target-s_EOS) < 1e-6);
		}
	}

	TEST_CASE( (char*)"Check units of fluid constants", (char*)"[fast]" ) 
	{
		SECTION((char*)"kSI")
		{
			set_standard_unit_system(UNIT_SYSTEM_KSI);

			SECTION((char*)"pcrit")
			{
				double p_Props1SI = Props1SI((char*)"R134a",(char*)"pcrit");
				double p_Props1 = Props1((char*)"R134a",(char*)"pcrit")*1000;
				double p_Props = Props((char*)"pcrit",'T',300,'Q',0,(char*)"R134a")*1000;
				double p_PropsSI = PropsSI((char*)"pcrit",(char*)"T",300,(char*)"Q",0,(char*)"R134a");
				REQUIRE(fabs(p_Props1/p_Props-1) < 1e-6);
				REQUIRE(fabs(p_PropsSI/p_Props-1) < 1e-6);
				REQUIRE(fabs(p_PropsSI/p_Props1SI-1) < 1e-6);
			}
		}
		SECTION((char*)"SI")
		{
			set_standard_unit_system(UNIT_SYSTEM_SI);

			SECTION((char*)"pcrit")
			{
				double p_Props1SI = Props1((char*)"R134a",(char*)"pcrit");
				double p_Props1 = Props1((char*)"R134a",(char*)"pcrit");
				double p_Props = Props((char*)"pcrit",'T',300,'Q',0,(char*)"R134a");
				double p_PropsSI = PropsSI((char*)"pcrit",(char*)"T",300,(char*)"Q",0,(char*)"R134a");
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