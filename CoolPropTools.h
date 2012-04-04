#ifndef COOLPROPTOOLS_H
#define COOLPROPTOOLS_H

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
	int ValidNumber(double x);

	int BuildLookupTable(char *Ref, struct fluidParamsVals *Fluid);
	double LookupValue(char Prop, double T, double p, char *Ref, struct fluidParamsVals *Fluid);
	int WriteLookup2File(int ILUT);

	int ValidateFluid(void);

	// Error handling things
	char CP_errString[5000];
	int ErrorFlag;
	void Append2ErrorString(char *string);
	void PrintError(void);
#endif
