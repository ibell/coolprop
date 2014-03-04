


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
				REQUIRE(fabs(p_Props1/p_Props-1) < 1e-6);
				REQUIRE(fabs(p_PropsSI/p_Props-1) < 1e-6);
				REQUIRE(fabs(p_PropsSI/p_Props1SI-1) < 1e-6);
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
				REQUIRE(fabs(p_Props1/p_Props-1) < 1e-6);
				REQUIRE(fabs(p_PropsSI/p_Props-1) < 1e-6);
				REQUIRE(fabs(p_PropsSI/p_Props1SI-1) < 1e-6);
			}
		}
	}

    TEST_CASE("Check ancillaries and surface tension", "[fast]" ) 
	{
		SECTION("surface tension")
		{
            REQUIRE_THROWS(_Props("I","P",101325,"D",1,"Propane"));
            REQUIRE_NOTHROW(_Props("I","T",300,"D",1,"Propane"));
		}
        SECTION("rhosatLanc")
		{
            REQUIRE_THROWS(_Props("rhosatLanc","P",101325,"D",1,"Propane"));
            REQUIRE_NOTHROW(_Props("rhosatLanc","T",300,"D",1,"Propane"));
		}
        SECTION("rhosatVanc")
		{
            REQUIRE_THROWS(_Props("rhosatVanc","P",101325,"D",1,"Propane"));
            REQUIRE_NOTHROW(_Props("rhosatVanc","T",300,"D",1,"Propane"));
		}
        SECTION("psatLanc")
		{
            REQUIRE_THROWS(_Props("psatLanc","P",101325,"D",1,"Propane"));
            REQUIRE_NOTHROW(_Props("psatLanc","T",300,"D",1,"Propane"));
		}
        SECTION("psatVanc")
		{
            REQUIRE_THROWS(_Props("psatVanc","P",101325,"D",1,"Propane"));
            REQUIRE_NOTHROW(_Props("psatVanc","T",300,"D",1,"Propane"));
		}
	}
	
	class validator_element
	{
	public:
		std::string in1, in3, in5;
		double in2, in4, in6;
		validator_element(std::string in1, double in2, std::string in3, double in4, std::string in5, double in6) : in1(in1), in2(in2), in3(in3), in4(in4), in5(in5), in6(in6) {};
	};

	TEST_CASE(  "Cyclohexane",  "[cyclohexane],[validation]" ) 
	{
		Fluid *CHEX = get_fluid(get_Fluid_index("Cyclohexane"));
		double mm = CHEX->params.molemass;
		validator_element data[] = {
			validator_element("T",300.0,"D",9.4*mm,"P",24.173705*1e6),
			validator_element("T",300.0,"D",9.4*mm,"O",115.28612/mm*1000),
			validator_element("T",300.0,"D",9.4*mm,"C",154.76968/mm*1000),
			validator_element("T",300.0,"D",9.4*mm,"A",1383.3876),
			validator_element("T",500.0,"D",6.5*mm,"P",3.9246630*1e6),
			validator_element("T",500.0,"D",6.5*mm,"O",192.52079/mm*1000),
			validator_element("T",500.0,"D",6.5*mm,"C",255.57110/mm*1000),
			validator_element("T",500.0,"D",6.5*mm,"A",434.13058),
			validator_element("T",500.0,"D",0.7*mm,"P",1.9981172*1e6),
			validator_element("T",500.0,"D",0.7*mm,"O",191.96468/mm*1000),
			validator_element("T",500.0,"D",0.7*mm,"C",235.52304/mm*1000),
			validator_element("T",500.0,"D",0.7*mm,"A",155.34798),
			validator_element("T",600.0,"D",3.5*mm,"P",6.8225506*1e6),  
			validator_element("T",600.0,"D",3.5*mm,"O",232.79249/mm*1000),
			validator_element("T",600.0,"D",3.5*mm,"C",388.55212/mm*1000),
			validator_element("T",600.0,"D",3.5*mm,"A",150.53314),
			validator_element("T",553.6,"D",3.3*mm,"P",4.0805433*1e6),
			validator_element("T",553.6,"D",3.3*mm,"O",224.19580/mm*1000), 
			validator_element("T",553.6,"D",3.3*mm,"C",199224.62/mm*1000), 
			validator_element("T",553.6,"D",3.3*mm,"A",87.913862)
		};

		//Now actually construct the vector
		std::vector<validator_element> elements(data, data + sizeof(data) / sizeof(validator_element));		
		
		for (std::vector<validator_element>::iterator it = elements.begin(); it != elements.end(); it++)
		{
			validator_element &el = *it;
			double eos = PropsSI( el.in5.c_str(),  el.in1.c_str(), el.in2,  el.in3.c_str(), el.in4, "Cyclohexane");
			double valid = el.in6;
            CAPTURE(eos);
            CAPTURE(valid);
            CAPTURE(el.in1);
            CAPTURE(el.in3);
            CAPTURE(el.in5);
			CHECK(fabs(valid/eos-1) < 1e-8);
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