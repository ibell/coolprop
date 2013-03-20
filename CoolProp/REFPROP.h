
#ifndef REFPROP_H
#define REFPROP_H

	// Only add REFPROP if build on Windows or Linux platform
	#if defined(__ISWINDOWS__) || defined(__ISLINUX__)
	
	class REFPROPFluidClass: public Fluid
	{
	private:
		std::vector<double> xmol;
	public:
		REFPROPFluidClass();
		REFPROPFluidClass(std::string FluidName, std::vector<double> xmol);

		double phir(double tau, double delta);
		double dphir_dDelta(double tau, double delta);
		double dphir_dTau(double tau, double delta);
		double d2phir_dTau2(double tau, double delta);

		double phi0(double tau, double delta);
		double dphi0_dTau(double tau, double delta);

		double viscosity_Trho(double T, double rho);
		double conductivity_Trho(double T, double rho);

		void saturation_T(double T, bool UseLUT, double *psatLout, double *psatVout, double *rhosatLout, double *rhosatVout);
		void saturation_p(double p, bool UseLUT, double *TsatLout, double *TsatVout, double *rhosatLout, double *rhosatVout);

		// Flash routines
		void temperature_ph(double p, double h, double *Tout, double *rhoout, double *rhoLout, double *rhoVout, double *TsatLout, double *TsatVout, double T0, double rho0);
		double density_Tp(double T, double p);
		double density_Tp(double T, double p, double rho_guess);

		double psat(double T);

	};
	bool set_REFPROP_fluid(std::string Ref);
	std::string get_REFPROP_fluid_path();

	double REFPROP(char Output,        char Name1,        double Prop1, char Name2,        double Prop2, char * Ref);
	double REFPROP(std::string Output, std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref);
	#endif

#endif
