#ifndef INDUSTRIALFLUIDS_H
#define INDUSTRIALFLUIDS_H


class CarbonMonoxideClass : public Fluid {

public:
    CarbonMonoxideClass();
    ~CarbonMonoxideClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){ 
		// Poling
		*e_k = 91.7; *sigma = 0.369; 
	};
	double surface_tension_T(double T){
		// From Mulero, 2012 JPCRD
		return 0.02843*pow(1-T/reduce.T,1.148);
	};

};

class CarbonylSulfideClass : public Fluid {

public:
    CarbonylSulfideClass();
    ~CarbonylSulfideClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){
		// Poling
		*e_k = 336; *sigma = 0.413; 
	};
	double surface_tension_T(double T)
	{
		// From Mulero, 2012 JPCRD
		return 0.07246*pow(1-T/reduce.T,1.407);
	}

};

class DecaneClass : public Fluid {

public:
    DecaneClass();
    ~DecaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma);
	double viscosity_Trho(double, double);
	double conductivity_Trho(double, double);
	double surface_tension_T(double T)
	{
		// From Mulero, 2012 JPCRD
		return 0.05473*pow(1-T/reduce.T,1.29);
	}
};

class HydrogenSulfideClass : public Fluid {

public:
    HydrogenSulfideClass();
    ~HydrogenSulfideClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma);
	double viscosity_Trho(double T, double rho);
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.078557*pow(1-T/reduce.T,1.2074);
	}

};

class IsopentaneClass : public Fluid {

public:
    IsopentaneClass();
    ~IsopentaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){
		//Chichester
		*e_k = 341.06; *sigma = 0.56232; 
	};
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.051*pow(1-T/reduce.T,1.209);
	}

};

class NeopentaneClass : public Fluid {

public:
    NeopentaneClass();
    ~NeopentaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){ 
		// Chichester
		*e_k = 191; *sigma = 0.644; 
	};

};

class IsohexaneClass : public Fluid {

public:
    IsohexaneClass();
    ~IsohexaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){
		// Chichester
		*e_k = 395.2; *sigma = 0.5799; 
	};
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.05024*pow(1-T/reduce.T,1.194);
	}
};

class KryptonClass : public Fluid {

public:
    KryptonClass();
    ~KryptonClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){ 
		// Poling
		*e_k = 178.9; *sigma = 0.3655; 
	};
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.0447*pow(1-T/reduce.T,1.245);
	}

};

class NonaneClass : public Fluid {

public:
    NonaneClass();
    ~NonaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma);
	double viscosity_Trho(double, double);
	double conductivity_Trho(double, double);
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.05388*pow(1-T/reduce.T,1.262);
	}
};

class TolueneClass : public Fluid {

public:
    TolueneClass();
    ~TolueneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double conductivity_Trho(double T, double rho);
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.06897*pow(1-T/reduce.T,1.291);
	}
};

class XenonClass : public Fluid {

public:
    XenonClass();
    ~XenonClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){
		// Poling
		*e_k = 231; *sigma = 0.4047; 
	};
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return -0.11538*pow(1-T/reduce.T,1.0512)+0.16598*pow(1-T/reduce.T,1.098);
	}
};

class R116Class : public Fluid {

public:
    R116Class();
    ~R116Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
    void ECSParams(double *e_k, double *sigma);
    double ECS_psi_viscosity(double rhor);
    double ECS_chi_conductivity(double rhor);
    double ECS_f_int(double T);
	double surface_tension_T(double T)
	{
		// From Mulero, 2012, JPCRD
		return 0.047593*pow(1-T/reduce.T,1.2666)-0.0073402*pow(1-T/reduce.T,1.9892);
	}
};

class AcetoneClass : public Fluid {

public:
    AcetoneClass();
    ~AcetoneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){
		// Poling
		*e_k = 560.2; *sigma = 0.46; 
	};
	double surface_tension_T(double T){
		// From Mulero, 2012, JPCRD
		return 0.0633*pow(1-T/reduce.T,1.16);
	};
};

class R245faClass : public Fluid {

public:
    R245faClass();
    ~R245faClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
    void ECSParams(double *e_k, double *sigma);
    double ECS_psi_viscosity(double rhor);
    double ECS_chi_conductivity(double rhor);
    double ECS_f_int(double T);
	double surface_tension_T(double T);
	//double viscosity_Trho(double T, double rho);
};

class R41Class : public Fluid {

public:
    R41Class();
    ~R41Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){
		// Chichester
		*e_k = 244.88; *sigma = 0.4123; 
	};
	double surface_tension_T(double T){
		// From Mulero, 2012, JPCRD
		return 0.05049*pow(1-T/reduce.T,1.242);
	};
};

class NitrousOxideClass : public Fluid {

public:
    NitrousOxideClass();
    ~NitrousOxideClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){
		// Poling
		*e_k = 232.4; *sigma = 0.3828; 
	};
	double surface_tension_T(double T){
		// From Mulero, 2012, JPCRD
		return 0.07087*pow(1-T/reduce.T,1.204);
	};
};

class SulfurDioxideClass : public Fluid {

public:
    SulfurDioxideClass();
    ~SulfurDioxideClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	void ECSParams(double *e_k, double *sigma){
		// Poling
		*e_k = 335.4; *sigma = 0.4112; 
	};
	double surface_tension_T(double T){
		// From Mulero, 2012, JPCRD
		return 0.0803*pow(1-T/reduce.T,0.928)+0.0139*pow(1-T/reduce.T,1.57)-0.0114*pow(1-T/reduce.T,0.364);
	};
};

class R141bClass : public Fluid {

public:
    R141bClass();
    ~R141bClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
    void ECSParams(double *e_k, double *sigma);
    double ECS_psi_viscosity(double rhor);
    double ECS_chi_conductivity(double rhor);
    double ECS_f_int(double T);
	double surface_tension_T(double T){
		// From Mulero, 2012, JPCRD
		return 0.000073958*pow(1-T/reduce.T,0.066331)+0.059941*pow(1-T/reduce.T,1.2214);
	};
};

class R142bClass : public Fluid {

public:
    R142bClass();
    ~R142bClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
    void ECSParams(double *e_k, double *sigma);
    double ECS_psi_viscosity(double rhor);
    double ECS_chi_conductivity(double rhor);
    double ECS_f_int(double T);
	double surface_tension_T(double T){
		// From Mulero, 2012, JPCRD
		return 0.05685*pow(1-T/reduce.T,1.237);
	};
};

class R218Class : public Fluid {

public:
    R218Class();
    ~R218Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
    void ECSParams(double *e_k, double *sigma);
    double ECS_psi_viscosity(double rhor);
    double ECS_chi_conductivity(double rhor);
    double ECS_f_int(double T);
	double surface_tension_T(double T){
		// From Mulero, 2012, JPCRD
		return 0.04322*pow(1-T/reduce.T,1.224);
	};
};


#endif
