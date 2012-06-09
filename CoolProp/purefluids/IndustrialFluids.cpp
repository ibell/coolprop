

// **** WARNING ******
// **** WARNING ******
// **** WARNING ******

// Do NOT modify this file.  It is created by a script in the industrialfluidsbuilder folder within the source

#include "CoolProp.h"
#include "IndustrialFluids.h"
#include <vector>
#include "CPExceptions.h"

static const double d_nonpolar[] =
{
0,
1.0, //[1]
1.0, //[2]
1.0, //[3]
2.0, //[4]
3.0, //[5]
7.0, //[6]
2.0, //[7]
5.0, //[8]
1.0, //[9]
4.0, //[10]
3.0, //[11]
4.0, //[12]
};

static const double t_nonpolar[] =
{
0,
0.25,  //[1]
1.125, //[2]
1.5,   //[3]
1.375, //[4]
0.25,  //[5]
0.875, //[6]
0.625, //[7]
1.75,  //[8]
3.625, //[9]
3.625, //[10]
14.5,  //[11]
12.0,  //[12]
};

static const double c_nonpolar[] =
{
0,
0.0, //[1]
0.0, //[2]
0.0, //[3]
0.0, //[4]
0.0, //[5]
0.0, //[6]
1.0, //[7]
1.0, //[8]
2.0, //[9]
2.0, //[10]
3.0, //[11]
3.0, //[12]
};

static const double d_polar[] =
{
0,
1.0, //[1]
1.0, //[2]
1.0, //[3]
3.0, //[4]
7.0, //[5]
1.0, //[6]
2.0, //[7]
5.0, //[8]
1.0, //[9]
1.0, //[10]
4.0, //[11]
2.0, //[12]
};

static const double t_polar[] =
{
0,
0.25,  //[1]
1.25,  //[2]
1.5,   //[3]
0.25,  //[4]
0.875, //[5]
2.375, //[6]
2.0,   //[7]
2.125, //[8]
3.5,   //[9]
6.5,   //[10]
4.75,  //[11]
12.5,  //[12]
};

static const double c_polar[] =
{
0,
0.0, //[1]
0.0, //[2]
0.0, //[3]
0.0, //[4]
0.0, //[5]
1.0, //[6]
1.0, //[7]
1.0, //[8]
2.0, //[9]
2.0, //[10]
2.0, //[11]
3.0, //[12]
};



CarbonMonoxideClass::CarbonMonoxideClass()
{
    const double n[]={0.0,0.905540000,-2.451500000,0.531490000,0.024173000,0.072156000,0.000188180,0.194050000,-0.043268000,-0.127780000,-0.027896000,-0.034154000,0.016329000};
    const double u0[]={0.0,3089.0};
    const double v0[]={0.0,1.0128};

    // Critical parameters
    crit.rho = 303.909585;
    crit.p = 3494.0;
    crit.T = 132.86;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 28.0101;
    params.Ttriple = 68.16;
    params.accentricfactor = 0.0497;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 68.16;
    limits.Tmax = 500.0;
    limits.pmax = 100000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-3.3728318564,3.3683460039);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(3.5-1);
    phi0list.push_back(phi0_logtau_);

    phi_BC * phi0_power_ = new phi0_power(-2.2311e-07*pow(crit.T,1.5)/(1.5*(1.5+1)),-1.5);
	phi0list.push_back(phi0_power_);
    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("CarbonMonoxide");
    aliases.push_back(std::string("CO")); 
    REFPROPname.assign("CO");
}
double CarbonMonoxideClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 314.41+562.237*pow(THETA,0.370304)+201.906*pow(THETA,0.921208);
}
double CarbonMonoxideClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.119344-3.94598*pow(THETA,0.529036)-8.93952*pow(THETA,1.85865)-33.3358*pow(THETA,5.11753);
    return exp(RHS)*crit.rho;
}
double CarbonMonoxideClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = 0.00381083-5.39772*pow(THETA,0.966354)-2.44727*pow(THETA,4.60196);
    return exp(crit.T/T*RHS)*crit.p;
}






CarbonylSulfideClass::CarbonylSulfideClass()
{
    const double n[]={0.0,0.943740000,-2.534800000,0.590580000,-0.021488000,0.082083000,0.000246890,0.212260000,-0.041251000,-0.223330000,-0.050828000,-0.028333000,0.016983000};
    const double u0[]={0.0,768.0,1363.0,3175.0,12829.0};
    const double v0[]={0.0,2.1651,0.93456,1.0623,0.34269};

    // Critical parameters
    crit.rho = 445.156491;
    crit.p = 6370.0;
    crit.T = 378.77;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 60.0751;
    params.Ttriple = 134.3;
    params.accentricfactor = 0.0978;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 134.3;
    limits.Tmax = 650.0;
    limits.pmax = 50000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-3.6587449805,3.7349245016);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(3.5-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("CarbonylSulfide");
    aliases.push_back(std::string("COS")); 
    REFPROPname.assign("COS");
}
double CarbonylSulfideClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 475.543+1033.64*pow(THETA,0.451523)+119.239*pow(THETA,3.25939);
}
double CarbonylSulfideClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.23617-14.9863*pow(THETA,2.50487)-5.12115*pow(THETA,0.604117)-57.4114*pow(THETA,7.01427);
    return exp(RHS)*crit.rho;
}
double CarbonylSulfideClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = 0.000586188-5.62785*pow(THETA,0.96735)-2.63574*pow(THETA,4.34562);
    return exp(crit.T/T*RHS)*crit.p;
}






DecaneClass::DecaneClass()
{
    const double n[]={0.0,1.046100000,-2.480700000,0.743720000,-0.525790000,0.153150000,0.000328650,0.841780000,0.055424000,-0.735550000,-0.185070000,-0.020775000,0.012335000};
    const double u0[]={0.0,1193.0,2140.0,4763.0,10862.0};
    const double v0[]={0.0,25.685,28.233,12.417,10.035};

    // Critical parameters
    crit.rho = 233.3419552;
    crit.p = 2103.0;
    crit.T = 617.7;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 142.28168;
    params.Ttriple = 243.5;
    params.accentricfactor = 0.4884;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 243.5;
    limits.Tmax = 675.0;
    limits.pmax = 800000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(13.9361966549,-10.5265128286);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(19.109-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Decane");
    aliases.push_back(std::string("decane")); 
    REFPROPname.assign("decane");
}
double DecaneClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 221.841+614.913*pow(THETA,0.380059)+141.621*pow(THETA,2.57228);
}
double DecaneClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.205737-22.9208*pow(THETA,2.51421)-6.45011*pow(THETA,0.587562)-84.889*pow(THETA,6.49204);
    return exp(RHS)*crit.rho;
}
double DecaneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.000777868-7.47275*pow(THETA,0.973801)-6.00257*pow(THETA,3.54267);
    return exp(crit.T/T*RHS)*crit.p;
}






HydrogenSulfideClass::HydrogenSulfideClass()
{
    const double n[]={0.0,0.876410000,-2.036700000,0.216340000,-0.050199000,0.066994000,0.000190760,0.202270000,-0.004534800,-0.222300000,-0.034714000,-0.014885000,0.007415400};
    const double u0[]={0.0,1823.0,3965.0};
    const double v0[]={0.0,1.1364,1.9721};

    // Critical parameters
    crit.rho = 347.2841672;
    crit.p = 9000.0;
    crit.T = 373.1;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 34.08088;
    params.Ttriple = 187.7;
    params.accentricfactor = 0.1005;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 187.7;
    limits.Tmax = 760.0;
    limits.pmax = 170000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-4.0740770957,3.7632137341);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);

    phi_BC * phi0_power_ = new phi0_power(-1.4327e-06*pow(crit.T,1.5)/(1.5*(1.5+1)),-1.5);
	phi0list.push_back(phi0_power_);
    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("HydrogenSulfide");
    aliases.push_back(std::string("H2S")); 
    REFPROPname.assign("H2S");
}
double HydrogenSulfideClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 341.396+807.261*pow(THETA,0.403989)+113.571*pow(THETA,1.42227);
}
double HydrogenSulfideClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.0867325-4.31985*pow(THETA,0.523612)-10.4004*pow(THETA,1.98277)-37.8195*pow(THETA,5.34652);
    return exp(RHS)*crit.rho;
}
double HydrogenSulfideClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = 0.00182417-2.48038*pow(THETA,4.1794)-5.63481*pow(THETA,0.966803);
    return exp(crit.T/T*RHS)*crit.p;
}






IsopentaneClass::IsopentaneClass()
{
    const double n[]={0.0,1.096300000,-3.040200000,1.031700000,-0.154100000,0.115350000,0.000298090,0.395710000,-0.045881000,-0.358040000,-0.101070000,-0.035484000,0.018156000};
    const double u0[]={0.0,442.0,1109.0,2069.0,4193.0};
    const double v0[]={0.0,7.4056,9.5772,15.765,12.119};

    // Critical parameters
    crit.rho = 235.99865938;
    crit.p = 3378.0;
    crit.T = 460.35;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 72.14878;
    params.Ttriple = 112.65;
    params.accentricfactor = 0.2274;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 112.65;
    limits.Tmax = 500.0;
    limits.pmax = 1000000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(2.5822330405,1.1609103419);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Isopentane");
    aliases.push_back(std::string("ipentane")); 
    REFPROPname.assign("ipentane");
}
double IsopentaneClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 263.636+558.226*pow(THETA,0.449723)+90.2293*pow(THETA,3.59523);
}
double IsopentaneClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.307787-7.37525*pow(THETA,0.764929)-104.82*pow(THETA,9.92148)-30.5914*pow(THETA,3.48774);
    return exp(RHS)*crit.rho;
}
double IsopentaneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.0029609-6.35049*pow(THETA,0.979701)-3.60812*pow(THETA,4.15062);
    return exp(crit.T/T*RHS)*crit.p;
}






NeopentaneClass::NeopentaneClass()
{
    const double n[]={0.0,1.113600000,-3.179200000,1.141100000,-0.104670000,0.117540000,0.000340580,0.295530000,-0.074765000,-0.314740000,-0.099401000,-0.039569000,0.023177000};
    const double u0[]={0.0,710.0,1725.0,3280.0,7787.0};
    const double v0[]={0.0,14.422,12.868,17.247,12.663};

    // Critical parameters
    crit.rho = 235.9265106;
    crit.p = 3196.0;
    crit.T = 433.74;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 72.14878;
    params.Ttriple = 256.6;
    params.accentricfactor = 0.1961;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 256.6;
    limits.Tmax = 550.0;
    limits.pmax = 200000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(0.8702452614,1.6071746358);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Neopentane");
    aliases.push_back(std::string("neopentn")); 
    REFPROPname.assign("neopentn");
}
double NeopentaneClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 259.311+759.07*pow(THETA,0.457452)-205.291*pow(THETA,0.457453);
}
double NeopentaneClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.105179-4.16981*pow(THETA,0.513219)-8.9438*pow(THETA,1.72519)-33.8828*pow(THETA,4.4748);
    return exp(RHS)*crit.rho;
}
double NeopentaneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = 0.0015322-3.14966*pow(THETA,3.75951)-6.09657*pow(THETA,0.970415);
    return exp(crit.T/T*RHS)*crit.p;
}






IsohexaneClass::IsohexaneClass()
{
    const double n[]={0.0,1.102700000,-2.969900000,1.029500000,-0.212380000,0.118970000,0.000277380,0.401030000,-0.034238000,-0.435840000,-0.116930000,-0.019262000,0.008078300};
    const double u0[]={0.0,325.0,1150.0,2397.0,5893.0};
    const double v0[]={0.0,7.9127,16.871,19.257,14.075};

    // Critical parameters
    crit.rho = 233.9661024;
    crit.p = 3040.0;
    crit.T = 497.7;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 86.17536;
    params.Ttriple = 119.6;
    params.accentricfactor = 0.2797;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 119.6;
    limits.Tmax = 550.0;
    limits.pmax = 1000000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(6.9259123919,-0.3128629679);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Isohexane");
    aliases.push_back(std::string("ihexane")); 
    REFPROPname.assign("ihexane");
}
double IsohexaneClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 238.626+578.165*pow(THETA,0.398393)+109.839*pow(THETA,2.87839);
}
double IsohexaneClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.241105-7.59813*pow(THETA,0.746707)-33.6159*pow(THETA,3.5051)-115.543*pow(THETA,10.0502);
    return exp(RHS)*crit.rho;
}
double IsohexaneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.00331591-4.1175*pow(THETA,4.06252)-6.62477*pow(THETA,0.982416);
    return exp(crit.T/T*RHS)*crit.p;
}






KryptonClass::KryptonClass()
{
    const double n[]={0.0,0.835610000,-2.372500000,0.545670000,0.014361000,0.066502000,0.000193100,0.168180000,-0.033133000,-0.150080000,-0.022897000,-0.021454000,0.006939700};
    const double u0[]={0.0,0};
    const double v0[]={0.0,0};

    // Critical parameters
    crit.rho = 909.2083;
    crit.p = 5525.0;
    crit.T = 209.48;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 83.798;
    params.Ttriple = 115.77;
    params.accentricfactor = -0.00089;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 115.77;
    limits.Tmax = 750.0;
    limits.pmax = 200000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-3.7506412806,3.7798018435);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(2.5-1);
    phi0list.push_back(phi0_logtau_);




    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Krypton");
    aliases.push_back(std::string("krypton")); 
    REFPROPname.assign("krypton");
}
double KryptonClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 809.494+1434.1*pow(THETA,0.270621)+916.309*pow(THETA,0.795785);
}
double KryptonClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.0392213-3.75692*pow(THETA,0.509614)-7.68651*pow(THETA,1.77205)-25.2873*pow(THETA,4.75346);
    return exp(RHS)*crit.rho;
}
double KryptonClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = 0.00424691-5.11498*pow(THETA,0.958996)-2.92151*pow(THETA,5.82349);
    return exp(crit.T/T*RHS)*crit.p;
}






NonaneClass::NonaneClass()
{
    const double n[]={0.0,1.115100000,-2.702000000,0.834160000,-0.388280000,0.137600000,0.000281850,0.620370000,0.015847000,-0.617260000,-0.150430000,-0.012982000,0.004432500};
    const double u0[]={0.0,1221.0,2244.0,5008.0,11724.0};
    const double v0[]={0.0,24.926,24.842,11.188,17.483};

    // Critical parameters
    crit.rho = 232.141731;
    crit.p = 2281.0;
    crit.T = 594.55;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 128.2551;
    params.Ttriple = 219.7;
    params.accentricfactor = 0.4433;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 219.7;
    limits.Tmax = 600.0;
    limits.pmax = 800000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(10.7927224829,-8.2418318753);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(17.349-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Nonane");
    aliases.push_back(std::string("nonane")); 
    REFPROPname.assign("nonane");
}
double NonaneClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 226.867+582.22*pow(THETA,0.361657)+143.829*pow(THETA,2.04318);
}
double NonaneClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.0977264-23.6568*pow(THETA,2.59633)-6.38283*pow(THETA,0.593543)-87.0597*pow(THETA,6.83533);
    return exp(RHS)*crit.rho;
}
double NonaneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.000227121-5.62753*pow(THETA,3.58749)-7.25067*pow(THETA,0.971535);
    return exp(crit.T/T*RHS)*crit.p;
}






TolueneClass::TolueneClass()
{
    const double n[]={0.0,0.964640000,-2.785500000,0.867120000,-0.188600000,0.118040000,0.000251810,0.571960000,-0.029287000,-0.433510000,-0.125400000,-0.028207000,0.014076000};
    const double u0[]={0.0,190.0,797.0,1619.0,3072.0,7915.0};
    const double v0[]={0.0,1.6994,8.0577,17.059,8.4567,8.6423};

    // Critical parameters
    crit.rho = 291.98665298;
    crit.p = 4126.0;
    crit.T = 591.75;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 92.13842;
    params.Ttriple = 178.0;
    params.accentricfactor = 0.2657;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 178.0;
    limits.Tmax = 700.0;
    limits.pmax = 500000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(3.5241174832,1.1360823464);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Toluene");
    aliases.push_back(std::string("toluene")); 
    REFPROPname.assign("toluene");
}
double TolueneClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 306.67+722.044*pow(THETA,0.418475)+130.709*pow(THETA,2.9141);
}
double TolueneClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.191773-6.8994*pow(THETA,0.691253)-25.5456*pow(THETA,3.08489)-89.434*pow(THETA,8.33106);
    return exp(RHS)*crit.rho;
}
double TolueneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = 0.00129047-4.14206*pow(THETA,3.99689)-6.4543*pow(THETA,0.969688);
    return exp(crit.T/T*RHS)*crit.p;
}






XenonClass::XenonClass()
{
    const double n[]={0.0,0.831150000,-2.355300000,0.539040000,0.014382000,0.066309000,0.000196490,0.149960000,-0.035319000,-0.159290000,-0.027521000,-0.023305000,0.008694100};
    const double u0[]={0.0,0};
    const double v0[]={0.0,0};

    // Critical parameters
    crit.rho = 1102.8612;
    crit.p = 5842.0;
    crit.T = 289.733;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 131.293;
    params.Ttriple = 161.4;
    params.accentricfactor = 0.00363;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 161.4;
    limits.Tmax = 750.0;
    limits.pmax = 700000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-3.8227178129,3.8416395351);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(2.5-1);
    phi0list.push_back(phi0_logtau_);




    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Xenon");
    aliases.push_back(std::string("Xe")); aliases.push_back(std::string("xenon")); 
    REFPROPname.assign("xenon");
}
double XenonClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 1032.78+1710.94*pow(THETA,0.284822)+1086.69*pow(THETA,0.779685);
}
double XenonClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.0761819-3.68385*pow(THETA,0.50101)-7.66799*pow(THETA,1.74507)-25.3412*pow(THETA,4.69948);
    return exp(RHS)*crit.rho;
}
double XenonClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = 0.00359722-5.14301*pow(THETA,0.960574)-2.64349*pow(THETA,5.51079);
    return exp(crit.T/T*RHS)*crit.p;
}






R116Class::R116Class()
{
    const double n[]={0.0,1.163200000,-2.812300000,0.772020000,-0.143310000,0.102270000,0.000246290,0.308930000,-0.028499000,-0.303430000,-0.068793000,-0.027218000,0.010665000};
    const double u0[]={0.0,190.0,622.0,1470.0};
    const double v0[]={0.0,2.4818,7.0622,7.9951};

    // Critical parameters
    crit.rho = 613.32452808;
    crit.p = 3048.0;
    crit.T = 293.03;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 138.01182;
    params.Ttriple = 173.1;
    params.accentricfactor = 0.2566;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 173.1;
    limits.Tmax = 425.0;
    limits.pmax = 50000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_nonpolar,d_nonpolar+sizeof(d_nonpolar)/sizeof(double));
    std::vector<double> t_v(t_nonpolar,t_nonpolar+sizeof(t_nonpolar)/sizeof(double));
    std::vector<double> l_v(c_nonpolar,c_nonpolar+sizeof(c_nonpolar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-10.7088650331,8.9148979056);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("R116");

    REFPROPname.assign("R116");
}
double R116Class::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 642.924+1289.38*pow(THETA,0.373632)+338.774*pow(THETA,1.06117);
}
double R116Class::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.161563-4.35824*pow(THETA,0.519469)-9.81421*pow(THETA,1.76751)-36.972*pow(THETA,4.5405);
    return exp(RHS)*crit.rho;
}
double R116Class::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = 0.00172169-3.43163*pow(THETA,3.63334)-6.36934*pow(THETA,0.968831);
    return exp(crit.T/T*RHS)*crit.p;
}
void R116Class::ECSParams(double *e_k, double *sigma)
{
    *e_k = 226.16;
    *sigma = 0.5249;
}
double R116Class::ECS_f_int(double T)
{
    return 0.00132+0*T;
}
double R116Class::ECS_psi_viscosity(double rhor)
{
    return 1.21996-0.0647835*rhor+0*rhor*rhor;
}
double R116Class::ECS_chi_conductivity(double rhor)
{
    return 1.1804-0.0539975*rhor;
}


AcetoneClass::AcetoneClass()
{
    const double n[]={0.0,0.90041000000,-2.12670000000,-0.08340900000,0.06568300000,0.00016527000,-0.03966300000,0.72085000000,0.00923180000,-0.17217000000,-0.14961000000,-0.07612400000,-0.01816600000};
    const double u0[]={0.0,310.0,3480.0,1576.0};
    const double v0[]={0.0,3.7072,7.0675,11.012};

    // Critical parameters
    crit.rho = 272.971958;
    crit.p = 4700.0;
    crit.T = 508.1;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 58.07914;
    params.Ttriple = 178.5;
    params.accentricfactor = 0.3071;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 178.5;
    limits.Tmax = 550.0;
    limits.pmax = 700000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-9.4883659997,7.1422719708);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("Acetone");
    aliases.push_back(std::string("acetone")); 
    REFPROPname.assign("acetone");
}
double AcetoneClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 312.663+703.191*pow(THETA,0.464385)+108.895*pow(THETA,3.3934);
}
double AcetoneClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.41049-21.187*pow(THETA,2.75107)-6.65354*pow(THETA,0.661221)-74.2311*pow(THETA,7.23056);
    return exp(RHS)*crit.rho;
}
double AcetoneClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = 0.000795221-3.84184*pow(THETA,3.82733)-6.6419*pow(THETA,0.968566);
    return exp(crit.T/T*RHS)*crit.p;
}






NitrousOxideClass::NitrousOxideClass()
{
    const double n[]={0.0,0.88045000000,-2.42350000000,0.38237000000,0.06891700000,0.00020367000,0.13122000000,0.46032000000,-0.00369850000,-0.23263000000,-0.00042859000,-0.04281000000,-0.02303800000};
    const double u0[]={0.0,879.0,2372.0,5447.0};
    const double v0[]={0.0,2.1769,1.6145,0.48393};

    // Critical parameters
    crit.rho = 452.011456;
    crit.p = 7245.0;
    crit.T = 309.52;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 44.0128;
    params.Ttriple = 182.33;
    params.accentricfactor = 0.1613;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 182.33;
    limits.Tmax = 525.0;
    limits.pmax = 50000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-4.4262736272,4.3120475243);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(3.5-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("NitrousOxide");
    aliases.push_back(std::string("N2O")); 
    REFPROPname.assign("N2O");
}
double NitrousOxideClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 403.808+595.753*pow(THETA,0.227187)+669.834*pow(THETA,0.741133);
}
double NitrousOxideClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.123997-3.75546*pow(THETA,0.47858)-8.99126*pow(THETA,1.66409)-32.9961*pow(THETA,4.52787);
    return exp(RHS)*crit.rho;
}
double NitrousOxideClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = 0.00124529-5.94036*pow(THETA,0.970591)-2.93925*pow(THETA,3.86231);
    return exp(crit.T/T*RHS)*crit.p;
}






SulfurDioxideClass::SulfurDioxideClass()
{
    const double n[]={0.0,0.93061000000,-1.95280000000,-0.17467000000,0.06152400000,0.00017711000,0.21615000000,0.51353000000,0.01041900000,-0.25286000000,-0.05472000000,-0.05985600000,-0.01652300000};
    const double u0[]={0.0,775.0,1851.0};
    const double v0[]={0.0,1.0620,1.9401};

    // Critical parameters
    crit.rho = 525.002841;
    crit.p = 7884.0;
    crit.T = 430.64;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 64.0638;
    params.Ttriple = 197.7;
    params.accentricfactor = 0.2557;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 197.7;
    limits.Tmax = 525.0;
    limits.pmax = 35000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-4.5328346436,4.4777967379);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);

    phi_BC * phi0_power_ = new phi0_power(-7.2453e-05*pow(crit.T,1.0)/(1.0*(1.0+1)),-1.0);
	phi0list.push_back(phi0_power_);
    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("SulfurDioxide");
    aliases.push_back(std::string("SO2")); 
    REFPROPname.assign("SO2");
}
double SulfurDioxideClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 556.41+1266.62*pow(THETA,0.413726)+210.182*pow(THETA,1.5635);
}
double SulfurDioxideClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.24882-13.6192*pow(THETA,2.10106)-4.90322*pow(THETA,0.553285)-53.2983*pow(THETA,5.65125);
    return exp(RHS)*crit.rho;
}
double SulfurDioxideClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.000667106-3.64759*pow(THETA,3.75619)-6.41012*pow(THETA,0.974263);
    return exp(crit.T/T*RHS)*crit.p;
}






R141bClass::R141bClass()
{
    const double n[]={0.0,1.14690000000,-3.67990000000,1.34690000000,0.08332900000,0.00025137000,0.32720000000,0.46946000000,-0.02982900000,-0.31621000000,-0.02621900000,-0.07804300000,-0.02049800000};
    const double u0[]={0.0,502.0,1571.0,4603.0};
    const double v0[]={0.0,6.8978,7.8157,3.2039};

    // Critical parameters
    crit.rho = 458.55946002;
    crit.p = 4212.0;
    crit.T = 477.5;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 116.94962;
    params.Ttriple = 169.68;
    params.accentricfactor = 0.2195;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 169.68;
    limits.Tmax = 500.0;
    limits.pmax = 400000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-15.5074814985,9.1871858933);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("R141b");

    REFPROPname.assign("R141b");
}
double R141bClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 509.965+1091.12*pow(THETA,0.439413)+185.865*pow(THETA,2.64809);
}
double R141bClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.267267-18.7606*pow(THETA,2.58933)-5.57637*pow(THETA,0.619539)-73.7862*pow(THETA,7.15925);
    return exp(RHS)*crit.rho;
}
double R141bClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.00308472-3.97593*pow(THETA,4.09844)-6.29057*pow(THETA,0.97953);
    return exp(crit.T/T*RHS)*crit.p;
}
void R141bClass::ECSParams(double *e_k, double *sigma)
{
    *e_k = 370.44;
    *sigma = 0.5493;
}
double R141bClass::ECS_f_int(double T)
{
    return 0.000521722+0.00000292456*T;
}
double R141bClass::ECS_psi_viscosity(double rhor)
{
    return 0.92135+0.041091*rhor+0*rhor*rhor;
}
double R141bClass::ECS_chi_conductivity(double rhor)
{
    return 1.0867-0.0216469*rhor;
}


R142bClass::R142bClass()
{
    const double n[]={0.0,1.00380000000,-2.76620000000,0.42921000000,0.08136300000,0.00024174000,0.48246000000,0.75542000000,-0.00743000000,-0.41460000000,-0.01655800000,-0.10644000000,-0.02170400000};
    const double u0[]={0.0,473.0,1256.0,2497.0,6840.0};
    const double v0[]={0.0,5.0385,6.8356,4.0591,2.8136};

    // Critical parameters
    crit.rho = 445.99694314;
    crit.p = 4055.0;
    crit.T = 410.26;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 100.49503;
    params.Ttriple = 142.72;
    params.accentricfactor = 0.2321;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 142.72;
    limits.Tmax = 470.0;
    limits.pmax = 60000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-12.6016527149,8.3160183265);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("R142b");

    REFPROPname.assign("R142b");
}
double R142bClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 500.671+1074.63*pow(THETA,0.442711)+200.08*pow(THETA,2.81729);
}
double R142bClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.252933-19.0218*pow(THETA,2.62874)-5.8241*pow(THETA,0.632573)-71.2366*pow(THETA,7.13992);
    return exp(RHS)*crit.rho;
}
double R142bClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.000556636-3.66973*pow(THETA,3.92032)-6.29906*pow(THETA,0.97278);
    return exp(crit.T/T*RHS)*crit.p;
}
void R142bClass::ECSParams(double *e_k, double *sigma)
{
    *e_k = 278.2;
    *sigma = 0.5362;
}
double R142bClass::ECS_f_int(double T)
{
    return 0.000940725+0.000000988196*T;
}
double R142bClass::ECS_psi_viscosity(double rhor)
{
    return 0.9716+0.019181*rhor+0*rhor*rhor;
}
double R142bClass::ECS_chi_conductivity(double rhor)
{
    return 1.0749-0.0177916*rhor;
}


R218Class::R218Class()
{
    const double n[]={0.0,1.32700000000,-3.84330000000,0.92200000000,0.11360000000,0.00036195000,1.10010000000,1.18960000000,-0.02514700000,-0.65923000000,-0.02796900000,-0.18330000000,-0.02163000000};
    const double u0[]={0.0,326.0,595.0,1489.0};
    const double v0[]={0.0,7.2198,7.2692,11.599};

    // Critical parameters
    crit.rho = 627.9845622;
    crit.p = 2640.0;
    crit.T = 345.02;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 188.01933;
    params.Ttriple = 125.45;
    params.accentricfactor = 0.3172;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 125.45;
    limits.Tmax = 440.0;
    limits.pmax = 20000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-15.6587335175,11.4531412796);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("R218");

    REFPROPname.assign("R218");
}
double R218Class::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 649.764+1507.68*pow(THETA,0.404379)+331.916*pow(THETA,2.60611);
}
double R218Class::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.12954-21.0428*pow(THETA,2.61018)-5.96899*pow(THETA,0.611065)-76.2274*pow(THETA,6.98835);
    return exp(RHS)*crit.rho;
}
double R218Class::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.00117049-4.27976*pow(THETA,3.67253)-6.68013*pow(THETA,0.972918);
    return exp(crit.T/T*RHS)*crit.p;
}
void R218Class::ECSParams(double *e_k, double *sigma)
{
    *e_k = 266.35;
    *sigma = 0.58;
}
double R218Class::ECS_f_int(double T)
{
    return 0.000892659+0.00000114912*T;
}
double R218Class::ECS_psi_viscosity(double rhor)
{
    return 1.10225-0.00550442*rhor+0*rhor*rhor;
}
double R218Class::ECS_chi_conductivity(double rhor)
{
    return 1.2877-0.0758811*rhor;
}


R245faClass::R245faClass()
{
    const double n[]={0.0,1.29000000000,-3.21540000000,0.50693000000,0.09314800000,0.00027638000,0.71458000000,0.87252000000,-0.01507700000,-0.40645000000,-0.11701000000,-0.13062000000,-0.02295200000};
    const double u0[]={0.0,222.0,1010.0,2450.0};
    const double v0[]={0.0,5.5728,10.385,12.554};

    // Critical parameters
    crit.rho = 516.084569;
    crit.p = 3651.0;
    crit.T = 427.16;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 134.04794;
    params.Ttriple = 171.05;
    params.accentricfactor = 0.3776;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 171.05;
    limits.Tmax = 440.0;
    limits.pmax = 200000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-13.4283638514,9.87236538);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);


    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("R245fa");

    REFPROPname.assign("R245fa");
}
double R245faClass::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 580.728+1264.86*pow(THETA,0.437088)+230.04*pow(THETA,2.8002);
}
double R245faClass::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.204607-20.291*pow(THETA,2.4712)-5.94045*pow(THETA,0.618455)-76.7531*pow(THETA,6.51668);
    return exp(RHS)*crit.rho;
}
double R245faClass::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.00254547-5.01251*pow(THETA,3.66623)-7.05974*pow(THETA,0.984797);
    return exp(crit.T/T*RHS)*crit.p;
}
void R245faClass::ECSParams(double *e_k, double *sigma)
{
    *e_k = 329.72;
    *sigma = 0.5529;
}
double R245faClass::ECS_f_int(double T)
{
    return 0.00164999-0.000000328868*T;
}
double R245faClass::ECS_psi_viscosity(double rhor)
{
    return 1.13848-0.0332328*rhor+0*rhor*rhor;
}
double R245faClass::ECS_chi_conductivity(double rhor)
{
    return 1.1627-0.0473491*rhor;
}

R41Class::R41Class()
{
    const double n[]={0.0,0.85316000000,-2.63660000000,0.69129000000,0.05468100000,0.00012796000,-0.37093000000,0.33920000000,-0.00174130000,-0.09541700000,-0.07885200000,-0.03072900000,-0.01149700000};
    const double u0[]={0.0,1841.0,4232.0};
    const double v0[]={0.0,5.6936,2.9351};

    // Critical parameters
    crit.rho = 316.506156;
    crit.p = 5897.0;
    crit.T = 317.28;
    crit.v = 1.0/crit.rho;

    // Other fluid parameters
    params.molemass = 34.03292;
    params.Ttriple = 129.82;
    params.accentricfactor = 0.2004;
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = 129.82;
    limits.Tmax = 425.0;
    limits.pmax = 70000.0;
    limits.rhomax = 1000000.0*params.molemass;    

    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v(d_polar,d_polar+sizeof(d_polar)/sizeof(double));
    std::vector<double> t_v(t_polar,t_polar+sizeof(t_polar)/sizeof(double));
    std::vector<double> l_v(c_polar,c_polar+sizeof(c_polar)/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));

    for (unsigned int i=0;i<u0_v.size();i++) { u0_v[i]/=crit.T; }

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead(-4.867644116,4.2527951258);
    phi0list.push_back(phi0_lead_);

    phi_BC * phi0_logtau_ = new phi0_logtau(4.0-1);
    phi0list.push_back(phi0_logtau_);

    phi_BC * phi0_power_ = new phi0_power(-0.00016937*pow(crit.T,1.0)/(1.0*(1.0+1)),-1.0);
	phi0list.push_back(phi0_power_);
    phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);
	phi0list.push_back(phi0_Planck_Einstein_);

    EOSReference.assign("Lemmon, E.W., and R. Span, \"Short Fundamental Equations of State for 20 Industrial Fluids,\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign("Using ECS");

    name.assign("R41");

    REFPROPname.assign("R41");
}
double R41Class::rhosatL(double T) 
{
    double THETA = 1-T/crit.T;
    return 369.324+1847.91*pow(THETA,0.501338)-1021.15*pow(THETA,0.501344);
}
double R41Class::rhosatV(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = -0.16056-13.6245*pow(THETA,2.23602)-5.29706*pow(THETA,0.574996)-52.0276*pow(THETA,6.16363);
    return exp(RHS)*crit.rho;
}
double R41Class::psat(double T) 
{
    double THETA = 1-T/crit.T;
    double RHS = 0.0041837-6.13932*pow(THETA,0.964557)-2.73068*pow(THETA,4.21954);
    return exp(crit.T/T*RHS)*crit.p;
}




