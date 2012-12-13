#ifndef INDUSTRIALFLUIDS_H
#define INDUSTRIALFLUIDS_H


class CarbonMonoxideClass : public Fluid {

public:
    CarbonMonoxideClass();
    ~CarbonMonoxideClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

};

class CarbonylSulfideClass : public Fluid {

public:
    CarbonylSulfideClass();
    ~CarbonylSulfideClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

};

class DecaneClass : public Fluid {

public:
    DecaneClass();
    ~DecaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

};

class HydrogenSulfideClass : public Fluid {

public:
    HydrogenSulfideClass();
    ~HydrogenSulfideClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

};

class IsopentaneClass : public Fluid {

public:
    IsopentaneClass();
    ~IsopentaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

};

class NeopentaneClass : public Fluid {

public:
    NeopentaneClass();
    ~NeopentaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

};

class IsohexaneClass : public Fluid {

public:
    IsohexaneClass();
    ~IsohexaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

};

class KryptonClass : public Fluid {

public:
    KryptonClass();
    ~KryptonClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

};

class NonaneClass : public Fluid {

public:
    NonaneClass();
    ~NonaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

};

class TolueneClass : public Fluid {

public:
    TolueneClass();
    ~TolueneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

};

class XenonClass : public Fluid {

public:
    XenonClass();
    ~XenonClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

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
};

class AcetoneClass : public Fluid {

public:
    AcetoneClass();
    ~AcetoneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

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
};

class R41Class : public Fluid {

public:
    R41Class();
    ~R41Class(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

};

class NitrousOxideClass : public Fluid {

public:
    NitrousOxideClass();
    ~NitrousOxideClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

};

class SulfurDioxideClass : public Fluid {

public:
    SulfurDioxideClass();
    ~SulfurDioxideClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);

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
};


#endif
