#ifndef SILOXANES_H
#define SILOXANES_H

//MDM
class OctamethyltrisiloxaneClass : public Fluid {

public:
    OctamethyltrisiloxaneClass();
    ~OctamethyltrisiloxaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double ECS_psi_viscosity(double rhor)
	{
		// Fit of saturated liquid data from DIPPR 801
		return 1.42859-0.144404*rhor;
	};
};

//MD2M
class DecamethyltetrasiloxaneClass : public Fluid {

public:
    DecamethyltetrasiloxaneClass();
    ~DecamethyltetrasiloxaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

//MD3M
class DodecamethylpentasiloxaneClass : public Fluid {

public:
    DodecamethylpentasiloxaneClass();
    ~DodecamethylpentasiloxaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

//D6
class DodecamethylcyclohexasiloxaneClass : public Fluid {

public:
    DodecamethylcyclohexasiloxaneClass();
    ~DodecamethylcyclohexasiloxaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

//MM
class HexamethyldisiloxaneClass : public Fluid {

public:
    HexamethyldisiloxaneClass();
    ~HexamethyldisiloxaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
	double ECS_psi_viscosity(double rhor)
	{
		// Fit of saturated liquid datasets from DIPPR 801
		return 1.35719-0.127828*rhor;
	};
};

//MD4M
class TetradecamethylhexasiloxaneClass : public Fluid {

public:
    TetradecamethylhexasiloxaneClass();
    ~TetradecamethylhexasiloxaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

//D4
class OctamethylcyclotetrasiloxaneClass : public Fluid {

public:
    OctamethylcyclotetrasiloxaneClass();
    ~OctamethylcyclotetrasiloxaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

//D5
class DecamethylcyclopentasiloxaneClass : public Fluid {

public:
    DecamethylcyclopentasiloxaneClass();
    ~DecamethylcyclopentasiloxaneClass(){};
    double psat(double);
    double rhosatL(double);
    double rhosatV(double);
};

#endif
