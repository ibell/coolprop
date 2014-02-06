#ifndef Hydrogen_H
#define Hydrogen_H

	class HydrogenClass : public Fluid{

	public:
		HydrogenClass();
		~HydrogenClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		void ECSParams(double *e_k, double *sigma)
		{
			// Poling
			*e_k = 59.7;
			*sigma = 0.2827;
		}
		double surface_tension_T(double);
	};

	class ParaHydrogenClass : public Fluid{

	public:
		ParaHydrogenClass();
		~ParaHydrogenClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
		double surface_tension_T(double);
	};

	class OrthoHydrogenClass : public Fluid{

	public:
		OrthoHydrogenClass();
		~OrthoHydrogenClass(){};
		double psat(double);
		double rhosatL(double);
		double rhosatV(double);
	};


#endif
