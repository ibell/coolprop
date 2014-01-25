#ifndef AIR_H
#define AIR_H

	class AirClass : public Fluid{

	public:
		AirClass();
		~AirClass(){};
		virtual double conductivity_Trho(double, double);
		virtual double viscosity_Trho(double, double);
		double psatL(double);
		double psatV(double);
		double rhosatL(double);
		double rhosatV(double);

		struct CriticalStruct maxcondT;

		double X_tilde(double T,double tau,double delta);
	};

#endif
