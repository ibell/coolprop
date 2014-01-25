#ifndef R227EA_R365MFC_H
#define R227EA_R365MFC_H

class R227EAClass : public Fluid {

public:
    R227EAClass();
    ~R227EAClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T){ 
		// From Mulero, 2012, JPCRD
		return 0.06127*pow(1-T/reduce.T,1.192)-0.009516*pow(1-T/reduce.T,0.9795)-0.00192*pow(1-T/reduce.T,1.421);
	};
	void ECSParams(double *e_k, double *sigma)
	{
		// From Huber (2003)
		*e_k = 289.34; *sigma = 0.5746;
	}
	double ECS_f_int(double T)
	{
		// From Huber (2003)
		return 1.42313e-3+8.31496e-9*T;
	}
	double ECS_psi_viscosity(double rhor)
	{
		// From Huber (2003)
		return 0.76758+0.254482*rhor-5.33748e-2*rhor*rhor;
	}
	double ECS_chi_conductivity(double rhor)
	{
		// From Huber (2003)
		return 1.3122-8.74448e-2*rhor;
	}
};

class R365MFCClass : public Fluid {

public:
    R365MFCClass();
    ~R365MFCClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double surface_tension_T(double T){ 
		// From Mulero, 2012, JPCRD
		return 0.0534*pow(1-T/reduce.T,1.210);
	};

};

#endif
