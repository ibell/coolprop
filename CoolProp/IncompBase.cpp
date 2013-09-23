
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


//	SimpleIncompressible* liquid = new DowthermQClass();
//	double Aresult = 0.;
//	double AT =  150.0 + 273.15;
//	double Ap =  300.0;
//    Aresult = liquid->rho(AT,Ap);
//    printf("From object:       %f \n\n",Aresult);

	SecCoolSolutionClass* obj = new MethanolSolutionClass();

	double result = 0.;
    double x =   0.25;
    double T =   5.0 + 273.15;
    double p =  300.0;

	obj->testInputs(T,p,x);
	obj->testInputs(T-15,p,x);


//	printf("The parameters %s: \n","are");
//	int counter = 0;
//	for(int i=0; i<obj->cRho.size(); i++) {
//		for(int j=0; j<obj->cRho[i].size(); j++) {
//			counter += 1;
//			printf("Parameter %d is: %f \n",counter,obj->cRho[i][j]);
//		}
//	}


}
