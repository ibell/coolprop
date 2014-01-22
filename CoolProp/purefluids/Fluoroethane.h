#ifndef FLUOROETHANE_H
#define FLUOROETHANE_H

	class FluoroethaneClass : public Fluid{

	public:
		FluoroethaneClass();
		~FluoroethaneClass(){};
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		//void ECSParams(double *e_k, double *sigma);
		//double surface_tension_T(double T);
	};

#endif
