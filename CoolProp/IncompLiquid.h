
#include <string>
#include "CPExceptions.h"


#ifndef INCOMPRESSIBLE_LIQUID_H
#define INCOMPRESSIBLE_LIQUID_H

/**
Notes for developers:

If you want to add a fluid, add its definition to the header 
and then add it to the map in the constructor in LiquidsContainer 
in IncompLiquid.cpp
**/

/// The abstract base class for the liquids
class IncompressibleLiquid{
	
protected:
	std::string name,description,reference;
public:
	// Constructor
	IncompressibleLiquid(){};

	///< Destructor.  No implementation
	~IncompressibleLiquid(){};

	virtual double rho(double T_K){return 0;};
	virtual double cp(double T_K){return 0;};
	virtual double s(double T_K){return 0;};
	virtual double u(double T_K){return 0;};
	virtual double visc(double T_K){return 0;};
	virtual double cond(double T_K){return 0;};
	
	double h(double T_K, double p)
	{
		return u(T_K)+p/rho(T_K);
	};

	std::string get_name(){
		return name;
	}
};

bool IsIncompressibleLiquid(std::string name);
double IncompLiquid(long iOutput, double T, double p, long iFluid);
double IncompLiquid(long iOutput, double T, double p, std::string name);

class DEBLiquidClass : public IncompressibleLiquid{

public:

    // Constructor
    DEBLiquidClass(){
        name = std::string("DEB");
        description = std::string("Diethylbenzene mixture");
        reference = std::string("Dowtherm J Dow Chemical Co. - from Ake Melinder, 2010, \"Properties of Secondary Working Fluids for Indirect Systems\", IIR");
    };

    ///< Destructor.  No implementation
    ~DEBLiquidClass(){};

    double rho(double T_K){
        return -0.731182*T_K+1076.5;
    }
    double cp(double T_K){
        return 0.00287576*T_K+0.999729;
    }
    double u(double T_K){
        return 0.00287576*(T_K*T_K-298*298)/2.0+0.999729*(T_K-298);
    }
    double s(double T_K){
        return 0.00287576*(T_K-298)/2.0+0.999729*log(T_K/298);
    }
    double visc(double T_K){
        return exp(7.03331e-05*T_K*T_K+-0.0566396*T_K+3.5503);
    }
    double cond(double T_K){
        return -2.06364e-07*T_K+0.000189132;
    }
};


class HCMLiquidClass : public IncompressibleLiquid{

public:

    // Constructor
    HCMLiquidClass(){
        name = std::string("HCM");
        description = std::string("Hydrocarbon mixture (synthetic)");
        reference = std::string("Therminol D12 (Gilotherm D12) Solutia - from Ake Melinder, 2010, \"Properties of Secondary Working Fluids for Indirect Systems\", IIR");
    };

    ///< Destructor.  No implementation
    ~HCMLiquidClass(){};

    double rho(double T_K){
        return -0.718788*T_K+971.725;
    }
    double cp(double T_K){
        return 0.00431212*T_K+0.844023;
    }
    double u(double T_K){
        return 0.00431212*(T_K*T_K-298*298)/2.0+0.844023*(T_K-298);
    }
    double s(double T_K){
        return 0.00431212*(T_K-298)/2.0+0.844023*log(T_K/298);
    }
    double visc(double T_K){
        return exp(0.000209096*T_K*T_K+-0.14706*T_K+18.3237);
    }
    double cond(double T_K){
        return -1.51212e-07*T_K+0.000153716;
    }
};


class HFELiquidClass : public IncompressibleLiquid{

public:

    // Constructor
    HFELiquidClass(){
        name = std::string("HFE");
        description = std::string("Hydrofluoroether");
        reference = std::string("HFE-7100 3M Novec - from Ake Melinder, 2010, \"Properties of Secondary Working Fluids for Indirect Systems\", IIR");
    };

    ///< Destructor.  No implementation
    ~HFELiquidClass(){};

    double rho(double T_K){
        return -0.918485*T_K+1822.37;
    }
    double cp(double T_K){
        return 0.000858788*T_K+0.871834;
    }
    double u(double T_K){
        return 0.000858788*(T_K*T_K-298*298)/2.0+0.871834*(T_K-298);
    }
    double s(double T_K){
        return 0.000858788*(T_K-298)/2.0+0.871834*log(T_K/298);
    }
    double visc(double T_K){
        return exp(7.39823e-06*T_K*T_K+-0.0114765*T_K+-4.22878);
    }
    double cond(double T_K){
        return -8.33333e-08*T_K+9.92958e-05;
    }
};


class PMS1LiquidClass : public IncompressibleLiquid{

public:

    // Constructor
    PMS1LiquidClass(){
        name = std::string("PMS1");
        description = std::string("Polydimethylsiloxan 1.");
        reference = std::string("Baysilone KT3 - from Ake Melinder, 2010, \"Properties of Secondary Working Fluids for Indirect Systems\", IIR");
    };

    ///< Destructor.  No implementation
    ~PMS1LiquidClass(){};

    double rho(double T_K){
        return -0.9025*T_K+1172.35;
    }
    double cp(double T_K){
        return 0.00148417*T_K+1.22369;
    }
    double u(double T_K){
        return 0.00148417*(T_K*T_K-298*298)/2.0+1.22369*(T_K-298);
    }
    double s(double T_K){
        return 0.00148417*(T_K-298)/2.0+1.22369*log(T_K/298);
    }
    double visc(double T_K){
        return exp(7.51428e-05*T_K*T_K+-0.0636352*T_K+6.36183);
    }
    double cond(double T_K){
        return -2.84167e-07*T_K+0.000207526;
    }
};


class PMS2LiquidClass : public IncompressibleLiquid{

public:

    // Constructor
    PMS2LiquidClass(){
        name = std::string("PMS2");
        description = std::string("Polydimethylsiloxan 2.");
        reference = std::string("Syltherm XLT Dow Corning Co. - from Ake Melinder, 2010, \"Properties of Secondary Working Fluids for Indirect Systems\", IIR");
    };

    ///< Destructor.  No implementation
    ~PMS2LiquidClass(){};

    double rho(double T_K){
        return -1.02576*T_K+1155.94;
    }
    double cp(double T_K){
        return 0.00210788*T_K+1.15355;
    }
    double u(double T_K){
        return 0.00210788*(T_K*T_K-298*298)/2.0+1.15355*(T_K-298);
    }
    double s(double T_K){
        return 0.00210788*(T_K-298)/2.0+1.15355*log(T_K/298);
    }
    double visc(double T_K){
        return exp(8.09988e-05*T_K*T_K+-0.065582*T_K+5.66926);
    }
    double cond(double T_K){
        return -2.11212e-07*T_K+0.000172305;
    }
};


class SABLiquidClass : public IncompressibleLiquid{

public:

    // Constructor
    SABLiquidClass(){
        name = std::string("SAB");
        description = std::string("Synthetic alkyl benzene");
        reference = std::string("Marlotherm X - from Ake Melinder, 2010, \"Properties of Secondary Working Fluids for Indirect Systems\", IIR");
    };

    ///< Destructor.  No implementation
    ~SABLiquidClass(){};

    double rho(double T_K){
        return -0.801667*T_K+1102.34;
    }
    double cp(double T_K){
        return 0.00151667*T_K+1.36094;
    }
    double u(double T_K){
        return 0.00151667*(T_K*T_K-298*298)/2.0+1.36094*(T_K-298);
    }
    double s(double T_K){
        return 0.00151667*(T_K-298)/2.0+1.36094*log(T_K/298);
    }
    double visc(double T_K){
        return exp(8.5066e-05*T_K*T_K+-0.0665792*T_K+5.21288);
    }
    double cond(double T_K){
        return -2.61667e-07*T_K+0.000208374;
    }
};


class HCBLiquidClass : public IncompressibleLiquid{

public:

    // Constructor
    HCBLiquidClass(){
        name = std::string("HCB");
        description = std::string("Hydrocarbon blend");
        reference = std::string("Dynalene MV - from Ake Melinder, 2010, \"Properties of Secondary Working Fluids for Indirect Systems\", IIR");
    };

    ///< Destructor.  No implementation
    ~HCBLiquidClass(){};

    double rho(double T_K){
        return -0.772024*T_K+1071.78;
    }
    double cp(double T_K){
        return 0.00352976*T_K+0.761393;
    }
    double u(double T_K){
        return 0.00352976*(T_K*T_K-298*298)/2.0+0.761393*(T_K-298);
    }
    double s(double T_K){
        return 0.00352976*(T_K-298)/2.0+0.761393*log(T_K/298);
    }
    double visc(double T_K){
        return exp(0.000130604*T_K*T_K+-0.0863212*T_K+7.16819);
    }
    double cond(double T_K){
        return -2.3869e-07*T_K+0.000203186;
    }
};


class TCOLiquidClass : public IncompressibleLiquid{

public:

    // Constructor
    TCOLiquidClass(){
        name = std::string("TCO");
        description = std::string("Terpene from citrus oils");
        reference = std::string("d-Limonene - from Ake Melinder, 2010, \"Properties of Secondary Working Fluids for Indirect Systems\", IIR");
    };

    ///< Destructor.  No implementation
    ~TCOLiquidClass(){};

    double rho(double T_K){
        return -0.778166*T_K+1071.02;
    }
    double cp(double T_K){
        return 0.0052159*T_K+0.223775;
    }
    double u(double T_K){
        return 0.0052159*(T_K*T_K-298*298)/2.0+0.223775*(T_K-298);
    }
    double s(double T_K){
        return 0.0052159*(T_K-298)/2.0+0.223775*log(T_K/298);
    }
    double visc(double T_K){
        return exp(1.14086e-06*T_K*T_K+-0.0107031*T_K+-3.47971);
    }
    double cond(double T_K){
        return -1.85052e-07*T_K+0.000174156;
    }
};


#endif