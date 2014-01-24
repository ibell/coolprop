#define CATCH_CONFIG_RUNNER
#include "Catch/catch.hpp"
#include "Tests.h"
#include "CoolPropDLL.h"
#include "CoolProp.h"

TEST_CASE( "Factorials are computed", "[factorial]" ) {
		REQUIRE( 1 == 1 );
	}

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




void run_tests()
{
	Catch::Session session; // There must be exactly once instance

	session.run();
}