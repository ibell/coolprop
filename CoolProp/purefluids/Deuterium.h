#ifndef DEUTERIUM_H
#define DEUTERIUM_H

	class DeuteriumClass : public Fluid{

	public:
		DeuteriumClass();
		~DeuteriumClass(){};
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		//void ECSParams(double *e_k, double *sigma)
		//{
		//	// Poling
		//	*e_k = 59.7;
		//	*sigma = 0.2827;
		//}
		double surface_tension_T(double);
	};

	class ParaDeuteriumClass : public Fluid{

	public:
		ParaDeuteriumClass();
		~ParaDeuteriumClass(){};
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
	};

	class OrthoDeuteriumClass : public Fluid{

	public:
		OrthoDeuteriumClass();
		~OrthoDeuteriumClass(){};
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
	};


#endif
