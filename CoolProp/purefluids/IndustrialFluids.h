#ifndef INDUSTRIALFLUIDS_H
#define INDUSTRIALFLUIDS_H


class CarbonMonoxideClass : public Fluid {

public:
    CarbonMonoxideClass();
    ~CarbonMonoxideClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class CarbonylSulfideClass : public Fluid {

public:
    CarbonylSulfideClass();
    ~CarbonylSulfideClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class DecaneClass : public Fluid {

public:
    DecaneClass();
    ~DecaneClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class HydrogenSulfideClass : public Fluid {

public:
    HydrogenSulfideClass();
    ~HydrogenSulfideClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class IsopentaneClass : public Fluid {

public:
    IsopentaneClass();
    ~IsopentaneClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class NeopentaneClass : public Fluid {

public:
    NeopentaneClass();
    ~NeopentaneClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class IsohexaneClass : public Fluid {

public:
    IsohexaneClass();
    ~IsohexaneClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class KryptonClass : public Fluid {

public:
    KryptonClass();
    ~KryptonClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class NonaneClass : public Fluid {

public:
    NonaneClass();
    ~NonaneClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class TolueneClass : public Fluid {

public:
    TolueneClass();
    ~TolueneClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class XenonClass : public Fluid {

public:
    XenonClass();
    ~XenonClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class R116Class : public Fluid {

public:
    R116Class();
    ~R116Class(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class AcetoneClass : public Fluid {

public:
    AcetoneClass();
    ~AcetoneClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class NitrousOxideClass : public Fluid {

public:
    NitrousOxideClass();
    ~NitrousOxideClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class SulfurDioxideClass : public Fluid {

public:
    SulfurDioxideClass();
    ~SulfurDioxideClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class R141bClass : public Fluid {

public:
    R141bClass();
    ~R141bClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class R142bClass : public Fluid {

public:
    R142bClass();
    ~R142bClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class R218Class : public Fluid {

public:
    R218Class();
    ~R218Class(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class R245faClass : public Fluid {

public:
    R245faClass();
    ~R245faClass(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

class R41Class : public Fluid {

public:
    R41Class();
    ~R41Class(){};
    virtual double conductivity_Trho(double, double);
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};
#endif
