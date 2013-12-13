#ifdef WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#include <stdint.h>
#endif

#include <vector>
#include <math.h>
#include <string.h>
#include <iostream>
#include <map>
#include <utility>

#include "../../CoolProp/CoolPropTools.h"
#include "../../CoolProp/IncompBase.h"
#include "../../CoolProp/IncompLiquid.h"
#include "../../CoolProp/IncompSolution.h"

/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
 * windows and linux.
 * Taken from http://stackoverflow.com/questions/1861294/how-to-calculate-execution-time-of-a-code-snippet-in-c
 * */
uint64_t getTimeValue() {
#ifdef WIN32
	/* Windows */
	FILETIME ft;
	LARGE_INTEGER li;

	/* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
	 * to a LARGE_INTEGER structure. */
	GetSystemTimeAsFileTime(&ft);
	li.LowPart = ft.dwLowDateTime;
	li.HighPart = ft.dwHighDateTime;

	uint64 ret = li.QuadPart;
	ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
	ret /= 10; /* From 100 nano seconds (10^-7) to 1 microsecond (10^-6) intervals */

	return ret;
#else
	/* Linux */
	struct timeval tv;

	gettimeofday(&tv, NULL);

	uint64_t ret = tv.tv_usec;
	/* In micro seconds (10^-6), add the seconds (10^0)
	 * after converting them to microseconds (10^-6) */
	ret += (tv.tv_sec * 1000000);

	return ret;
#endif
}


double testObject(IncompressibleFluid* fluid, std::vector<double> params, int exponent, bool print){


	uint64_t runs = 10;
	for (int i=1; i<=exponent; i++){
		runs *= 10;
	}


	fluid->setDebug(false);

	uint64_t start = getTimeValue();
	for (int i=0; i<runs; i++){
		switch (params.size()) {
			case 0: break;
		  	case 1: break;
		  	case 2: fluid->cp(params[0], params[1]*(1+i/runs)); break;
		  	case 2: fluid->cp(params[0], params[1]*(1+i/runs), params[2]); break;
		  	default: break;
		}
	}
	uint64_t end = getTimeValue();

	fluid->setDebug(print);
	switch (params.size()) {
		case 0: break;
		case 1: break;
		case 2: fluid->cp(params[0], params[1]); break;
		case 2: fluid->cp(params[0], params[1], params[2]); break;
		default: break;
	}

	return (end-start)/runs;
}


int main(int argc, const char* argv[]) {



	std::vector<double> params;

	SimpleIncompressible* liquid = new DowthermQClass();
	params.push_back(50.0 + 273.15);
	params.push_back(3e5);
	double time_liq = testObject(liquid, params, 6, false);

	IncompressibleSolution* solution = new EGSolution();
	params.push_back(0.25);
	double time_sol = testObject(solution, params, 6, false);


	std::cout << "Time consumption liquid: "   << time_liq << " µs per call from " << runs << " calls." << std::endl;
	std::cout << "Time consumption solution: " << time_sol << " µs per call from " << runs << " calls." << std::endl;


//	SecCoolSolution* obj = new MethanolSolution();
//	double x = 0.25;
//	double T = 5.0 + 273.15;
//	double p = 3e5;
//
//	obj->testInputs(T + 00, p, x);
//	obj->testInputs(T + 05, p, x);
//	obj->testInputs(T + 10, p, x);
//	obj->testInputs(T + 15, p, x);

}

