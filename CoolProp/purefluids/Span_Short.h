#ifndef SPAN_SHORT_H
#define SPAN_SHORT_H

class nPentaneClass : public Fluid {

public:
    nPentaneClass();
    ~nPentaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.08015*pow(1-T/reduce.T,1.408)+0.004384*pow(1-T/reduce.T,1.031)-0.03437*pow(1-T/reduce.T,1.818);
	}
};

class nHexaneClass : public Fluid {

public:
    nHexaneClass();
    ~nHexaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double viscosity_Trho(double, double);
	double conductivity_Trho(double, double);
	void ECSParams(double *e_k, double *sigma){*e_k = 378.4; *sigma = 0.6334;};
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.210952*pow(1-T/reduce.T,1.0962)+-0.158485*pow(1-T/reduce.T,1.05893);
	}
};

class nHeptaneClass : public Fluid {

public:
    nHeptaneClass();
    ~nHeptaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double conductivity_Trho(double, double);
	void ECSParams(double *e_k, double *sigma){
		// Chichester, 2008
		*e_k = 400; *sigma = 0.64947;
	};
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.07765*pow(1-T/reduce.T,1.319)+-0.02599*pow(1-T/reduce.T,1.6);
	}
};

class nOctaneClass : public Fluid {

public:
    nOctaneClass();
    ~nOctaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma);
	double viscosity_Trho(double, double);
	double conductivity_Trho(double, double);
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.34338*pow(1-T/reduce.T,1.6607)-0.50634*pow(1-T/reduce.T,1.9632)+0.2238*pow(1-T/reduce.T,2.3547);
	}
};

class nDodecaneClass : public Fluid {

public:
    nDodecaneClass();
    ~nDodecaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

	double viscosity_Trho(double T, double rho);
	double conductivity_Trho(double T, double rho);
	double viscosity_dilute(double T);
	double viscosity_background(double T, double rho);
	double conductivity_dilute(double T);
	double conductivity_background(double T, double rho);

	void ECSParams(double *e_k, double *sigma);
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.0154*pow(1-T/reduce.T,4.18)+0.048*pow(1-T/reduce.T,1.17);
	}
};

class CyclohexaneClass : public Fluid {

public:
    CyclohexaneClass();
    ~CyclohexaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){
		// Poling
		*e_k = 297.1; *sigma = 0.6182;
	};
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.06485*pow(1-T/reduce.T,1.263);
	}
};

class R152AClass : public Fluid {

public:
    R152AClass();
    ~R152AClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma);
	double viscosity_Trho(double, double);
	double conductivity_Trho(double, double);
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.05808*pow(1-T/reduce.T,1.2115);
	}
};

class R123Class : public Fluid {

public:
    R123Class();
    ~R123Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma);
	double viscosity_Trho(double T, double rho);
	double conductivity_Trho(double T, double rho);
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.056151*pow(1-T/reduce.T,1.2367);
	}
};

class R11Class : public Fluid {

public:
    R11Class();
    ~R11Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){
		// McLinden 2000
		*e_k = 363.61; *sigma = 0.5447;
	};
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.06212*pow(1-T/reduce.T, 1.247);
	}

};


#endif
