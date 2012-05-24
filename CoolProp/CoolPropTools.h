#include <vector>
#include <exception>
#include <iostream>
#include <list>
#include "FluidClass.h"

#if defined(_WIN32) || defined(__WIN32__) || defined(_WIN64) || defined(__WIN64__)
#define __ISWINDOWS__
#endif

#ifndef COOLPROPTOOLS_H
#define COOLPROPTOOLS_H

#include <string>
#include <cstdarg>

//missing string printf
//this is safe and convenient but not exactly efficient
inline std::string format(const char* fmt, ...){
    int size = 512;
    char* buffer = 0;
    buffer = new char[size];
    va_list vl;
    va_start(vl,fmt);
    int nsize = vsnprintf(buffer,size,fmt,vl);
    if(size<=nsize){//fail delete buffer and try again
        delete buffer; buffer = 0;
        buffer = new char[nsize+1];//+1 for /0
        nsize = vsnprintf(buffer,size,fmt,vl);
    }
    std::string ret(buffer);
    va_end(vl);
    delete buffer;
    return ret;
}

	#define OK 1
	#define FAIL 0

	struct fcnPointers{
		// A structure to contain points to residual Helmholtz Energy functions for the fluid
		double(*phir)(double,double);
		double(*dphir_dDelta)(double,double);
		double(*dphir2_dDelta2)(double,double);
		double(*dphir2_dDelta_dTau)(double,double);
		double(*dphir_dTau)(double,double);
		double(*dphir2_dTau2)(double,double);
		double(*phi0)(double,double);
		double(*dphi0_dDelta)(double,double);
		double(*dphi02_dDelta2)(double,double);
		double(*dphi0_dTau)(double,double);
		double(*dphi02_dTau2)(double,double);

		// Saturation curve functions
		double (*rhosatV)(double);
		double (*rhosatL)(double);
		//For the pure fluids
		double (*psat)(double);
		// For the psedo-pure fluids
		double (*p_dp)(double);
		double (*p_bp)(double);

		//Transport properties
		double (*visc)(double,double);
		double (*cond)(double,double);
	};

	struct LUTVals{
		double Tmin;
		double Tmax;
		double pmin;
		double pmax;
	};

	struct fluidParamsVals{
		struct LUTVals LUT;
		struct fcnPointers funcs;
		int Type;
		double Tc;
		double pc;
		double rhoc;
		double MM;
		double Tt;
		char Reference[1000];
	};

	double powInt(double x, int y);
	double QuadInterp(double x0, double x1, double x2, double f0, double f1, double f2, double x);
	double CubicInterp( double x0, double x1, double x2, double x3, double f0, double f1, double f2, double f3, double x);
	int ValidNumber(double x);

	int BuildLookupTable(char *Ref, struct fluidParamsVals *Fluid);
	double LookupValue_TP(char Prop, double T, double p, char *Ref, struct fluidParamsVals *Fluid);
    double LookupValue_Trho(char Prop, double T, double rho, char *Ref, struct fluidParamsVals *Fluid);
	int WriteLookup2File(int ILUT);

	int ValidateFluid(void);

	//void Append2ErrorString(char *string);
	//void PrintError(void);

	// This class contains all the classes for the fluids that can be used
	class FluidsContainer
	{
	private:
		std::list <Fluid*> FluidsList;
	public:
		// Constructor
		FluidsContainer();

		// Destructor
		~FluidsContainer();

		// Accessor
		Fluid * get_fluid(std::string name);

		//
		std::string FluidList();
	};

#endif
