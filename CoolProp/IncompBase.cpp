
#include <string>
#include <vector>
#include <math.h>
#include "CoolPropTools.h"
#include "CPExceptions.h"
#include "IncompBase.h"
#include "IncompLiquid.h"
#include "IncompSolution.h"

#include <stdio.h>

int main() {

	SimpleIncompressible* liquid = new DowthermQClass();
	double AT      =  150.0 + 273.15;
	double Ap      =  3e5;
    liquid->testInputs(AT,Ap);


	SecCoolSolutionClass* obj = new MethanolSolutionClass();
    double x      =   0.25;
    double T      =   5.0 + 273.15;
    double p      =   3e5;

	obj->testInputs(T+00,p,x);
	obj->testInputs(T+05,p,x);
	obj->testInputs(T+10,p,x);
	obj->testInputs(T+15,p,x);


}
