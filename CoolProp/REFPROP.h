
#ifndef REFPROP_H
#define REFPROP_H
	
	class REFPROPFluidClass: public Fluid
	{
	private:
		std::vector<double> xmol;
		static bool supported;

	public:
		REFPROPFluidClass();
		REFPROPFluidClass(std::string FluidName, std::vector<double> xmol);

		double R(void){return params.R_u/params.molemass;};
		static bool refpropSupported();

		double phir(double tau, double delta);
		double dphir_dDelta(double tau, double delta);
		double dphir_dTau(double tau, double delta);
		double d2phir_dDelta_dTau(double tau, double delta);
		double d2phir_dDelta2(double tau, double delta);
		double d2phir_dTau2(double tau, double delta);

		double phi0(double tau, double delta);
		double dphi0_dTau(double tau, double delta);
		double d2phi0_dTau2(double tau, double delta);

		double viscosity_Trho(double T, double rho);
		double conductivity_Trho(double T, double rho);
		double surface_tension_T(double T);

		void saturation_T(double T, bool UseLUT, double &psatLout, double &psatVout, double &rhosatLout, double &rhosatVout);
		void saturation_p(double p, bool UseLUT, double &TsatLout, double &TsatVout, double &rhosatLout, double &rhosatVout);

		// Flash routines
		void temperature_ph(double p, double h, double &Tout, double &rhoout, double &rhoLout, double &rhoVout, double &TsatLout, double &TsatVout, double T0, double rho0);
		void temperature_ps(double p, double s, double &Tout, double &rhoout, double &rhoLout, double &rhoVout, double &TsatLout, double &TsatVout);
		void temperature_hs(double h, double s, double &Tout, double &rhoout, double &rhoLout, double &rhoVout, double &TsatLout, double &TsatVout);
		double density_Tp(double T, double p);
		double density_Tp(double T, double p, double rho_guess);

		double psat(double T);
		double rhosatL(double T);
		double rhosatV(double T);

	};
	bool set_REFPROP_fluid(std::string Ref, std::vector<double> &x);
	std::string get_REFPROP_fluid_path();

    double REFPROPSI(long iOutput,     long iName1,       double iProp1, long iName2, double iProp2, std::string Ref);
	double REFPROP(char Output,        char Name1,        double Prop1, char Name2,        double Prop2, char * Ref);
	double REFPROP(std::string Output, std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref);

#endif
