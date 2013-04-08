#ifndef R227EA_R365MFC_H
#define R227EA_R365MFC_H

class R227EAClass : public Fluid {

public:
    R227EAClass();
    ~R227EAClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){
		// From Chichester NISTIR report 6650
		*e_k = 289.34; *sigma = 0.5746;};
	double surface_tension_T(double T){ 
		// From Mulero, 2012, JPCRD
		return 0.06127*pow(1-T/reduce.T,1.192)-0.009516*pow(1-T/reduce.T,0.9795)-0.00192*pow(1-T/reduce.T,1.421);
	};
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
