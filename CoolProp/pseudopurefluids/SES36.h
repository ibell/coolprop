#ifndef SES36_H
#define SES36_H

	class SES36Class : public Fluid{

	public:
		SES36Class();
		~SES36Class(){};
		double psatL(double);
		double psatV(double);
		double rhosatL(double);
		double rhosatV(double);
		double surface_tension_T(double);
	};

#endif
