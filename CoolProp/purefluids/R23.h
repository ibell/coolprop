#ifndef R23_H
#define R23_H

class R23Class : public Fluid {

public:
    R23Class();
    ~R23Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	//void ECSParams(double *e_k, double *sigma){*e_k = 341.1; *sigma = 0.5784;};
};

#endif
