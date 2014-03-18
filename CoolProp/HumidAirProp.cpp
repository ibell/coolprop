#if defined(_MSC_VER)
#define _CRT_SECURE_NO_WARNINGS
#endif



#include <stdlib.h>
#include "math.h"
#include "time.h"
#include "stdio.h"
#include <string.h>
#include <iostream>

#include "CoolProp.h"
#include "Ice.h"
#include "HumidAirProp.h"
#include "Solvers.h"
#include "CoolPropTools.h"

#include "purefluids/Water.h"
#include "pseudopurefluids/Air.h"

static WaterClass Water = WaterClass();
static AirClass Air = AirClass();

enum givens{GIVEN_TDP,GIVEN_HUMRAT,GIVEN_V,GIVEN_TWB,GIVEN_RH,GIVEN_ENTHALPY,GIVEN_T,GIVEN_P,GIVEN_VISC,GIVEN_COND};

static const char ITc='B';
static double epsilon=0.621945,R_bar=8.314472;
static int FlagUseVirialCorrelations=0,FlagUseIsothermCompressCorrelation=0,FlagUseIdealGasEnthalpyCorrelations=0;
double f_factor(double T, double p);

// A couple of convenience functions that are needed quite a lot
static double MM_Air(void)
{
    return Air.params.molemass;
}
static double MM_Water(void)
{
    return Water.params.molemass;
}
static double B_Air(double T)
{
    return DerivTerms((char *)"B",T,1e-12,(char *)"Air");
}
static double dBdT_Air(double T)
{
    return DerivTerms((char *)"dBdT",T,1e-12,(char *)"Air");
}
static double B_Water(double T)
{
    return DerivTerms((char *)"B",T,1e-12,(char *)"Water");
}
static double dBdT_Water(double T)
{
    return DerivTerms((char *)"dBdT",T,1e-12,(char *)"Water");
}
static double C_Air(double T)
{
    return DerivTerms((char *)"C",T,1e-12,(char *)"Air");
}
static double dCdT_Air(double T)
{
    return DerivTerms((char *)"dCdT",T,1e-12,(char *)"Air");
}
static double C_Water(double T)
{
    return DerivTerms((char *)"C",T,1e-12,(char *)"Water");
}
static double dCdT_Water(double T)
{
    return DerivTerms((char *)"dCdT",T,1e-12,(char *)"Water");
}
void UseVirialCorrelations(int flag)
{
    if (flag==0 || flag==1)
    {
        FlagUseVirialCorrelations=flag;
    }        
    else
    {
        printf("UseVirialCorrelations takes an integer, either 0 (no) or 1 (yes)\n");
    }
    
}
void UseIsothermCompressCorrelation(int flag)
{
    if (flag==0 || flag==1)
    {
        FlagUseIsothermCompressCorrelation=flag;
    }        
    else
    {
        printf("UseIsothermCompressCorrelation takes an integer, either 0 (no) or 1 (yes)\n");
    }
}
void UseIdealGasEnthalpyCorrelations(int flag)
{
    if (flag==0 || flag==1)
    {
        FlagUseIdealGasEnthalpyCorrelations=flag;
    }        
    else
    {
        printf("UseIdealGasEnthalpyCorrelations takes an integer, either 0 (no) or 1 (yes)\n");
    }
}
static double Brent_HAProps_T(char *OutputName, char *Input1Name, double Input1, char *Input2Name, double Input2, double TargetVal, double T_min, double T_max)
{
    double T;
    class BrentSolverResids : public FuncWrapper1D
    {
    private:
        double Input1,Input2,TargetVal;
        char *OutputName,*Input1Name,*Input2Name;
    public:
        double rhoL, rhoV;
        BrentSolverResids(char *OutputName, char *Input1Name, double Input1, char *Input2Name, double Input2, double TargetVal)
        {
            this->OutputName = OutputName;
            this->Input1Name = Input1Name;
            this->Input2Name = Input2Name;
            this->Input1 = Input1;
            this->Input2 = Input2;
            this->TargetVal = TargetVal;
        };
        ~BrentSolverResids(){};
        
        double call(double T){
            return HAProps(OutputName,(char *)"T",T,Input1Name,Input1,Input2Name,Input2)-TargetVal;
        }
    };

    BrentSolverResids *BSR = new BrentSolverResids(OutputName, Input1Name, Input1, Input2Name, Input2, TargetVal);

    std::string errstr;
    T = Brent(BSR,T_min,T_max,1e-7,1e-4,50,&errstr);

    delete BSR;
    return T;
}

static double Secant_HAProps_T(char *OutputName, char *Input1Name, double Input1, char *Input2Name, double Input2, double TargetVal, double T_guess)
{
    // Use a secant solve in order to yield a target output value for HAProps by altering T
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=5e-7,f=999,T=300,change;
    int iter=1;

    while ((iter<=3 || (fabs(f)>eps && fabs(change)>1e-10)) && iter<100)
    {
        if (iter==1){x1=T_guess; T=x1;}
        if (iter==2){x2=T_guess+0.001; T=x2;}
        if (iter>2) {T=x2;}
            f=HAProps(OutputName,(char *)"T",T,Input1Name,Input1,Input2Name,Input2)-TargetVal;
        if (iter==1){y1=f;}
        if (iter>1)
        {
            y2=f;
            x3=x2-y2/(y2-y1)*(x2-x1);
            change = y2/(y2-y1)*(x2-x1);
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
    }
    return T;
}

static double Secant_HAProps_W(const char *OutputName, const char *Input1Name, double Input1, const char *Input2Name, double Input2, double TargetVal, double W_guess)
{
    // Use a secant solve in order to yield a target output value for HAProps by altering humidity ratio
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,f=999,W=0.0001;
    int iter=1;

    while ((iter<=3 || fabs(f)>eps) && iter<100)
    {
        if (iter == 1){x1 = W_guess; W = x1;}
        if (iter == 2){x2 = W_guess+0.001; W = x2;}
        if (iter > 2) {W = x2;}
            f = HAProps(OutputName,(char *)"W",W,Input1Name,Input1,Input2Name,Input2)-TargetVal;
        if (iter == 1){y1 = f;}
        if (iter > 1)
        {
            y2=f;
            x3=x2-y2/(y2-y1)*(x2-x1);
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
    }
    return W;
}

// Mixed virial components
static double _B_aw(double T)
{
    // Function return has units of dm^3/mol, to convert to m^3/mol, divide by 1000
    double a[]={0,0.665687e2,-0.238834e3,-0.176755e3};
    double b[]={0,-0.237,-1.048,-3.183};
    double rhobarstar=1000,Tstar=100;
    return 1/rhobarstar*(a[1]*pow(T/Tstar,b[1])+a[2]*pow(T/Tstar,b[2])+a[3]*pow(T/Tstar,b[3]));
}

static double _dB_aw_dT(double T)
{
    // Function return has units of dm^3/mol/K, to convert to m^3/mol/K, divide by 1000
    double a[]={0,0.665687e2,-0.238834e3,-0.176755e3};
    double b[]={0,-0.237,-1.048,-3.183};
    double rhobarstar=1000,Tstar=100;
    return 1/rhobarstar/Tstar*(a[1]*b[1]*pow(T/Tstar,b[1]-1)+a[2]*b[2]*pow(T/Tstar,b[2]-1)+a[3]*b[3]*pow(T/Tstar,b[3]-1));
}

static double _C_aaw(double T)
{
    // Function return has units of dm^6/mol^2
    double c[]={0,0.482737e3,0.105678e6,-0.656394e8,0.294442e11,-0.319317e13};
    double rhobarstar=1000,Tstar=1,summer=0; int i;
    for (i=1;i<=5;i++)
    {
        summer+=c[i]*pow(T/Tstar,1-i);
    }
    return 1.0/rhobarstar/rhobarstar*summer;
}

static double _dC_aaw_dT(double T)
{
    // Function return has units of dm^6/mol^2/K
    double c[]={0,0.482737e3,0.105678e6,-0.656394e8,0.294442e11,-0.319317e13};
    double rhobarstar=1000,Tstar=1,summer=0; int i;
    for (i=2;i<=5;i++)
    {
        summer+=c[i]*(1-i)*pow(T/Tstar,-i);
    }
    return 1.0/rhobarstar/rhobarstar/Tstar*summer;
}

static double _C_aww(double T)
{
    // Function return has units of dm^6/mol^2
    double d[]={0,-0.1072887e2,0.347804e4,-0.383383e6,0.334060e8};
    double rhobarstar=1,Tstar=1,summer=0; int i;
    for (i=1;i<=4;i++)
    {
        summer+=d[i]*pow(T/Tstar,1-i);
    }
    return -1.0/rhobarstar/rhobarstar*exp(summer);
}

static double _dC_aww_dT(double T)
{
    // Function return has units of dm^6/mol^2/K
    double d[]={0,-0.1072887e2,0.347804e4,-0.383383e6,0.334060e8};
    double rhobarstar=1,Tstar=1,summer1=0,summer2=0; int i;
    for (i=1;i<=4;i++)
    {
        summer1+=d[i]*pow(T/Tstar,1-i);
    }
    for (i=2;i<=4;i++)
    {
        summer2+=d[i]*(1-i)*pow(T/Tstar,-i);
    }
    return -1.0/rhobarstar/rhobarstar/Tstar*exp(summer1)*summer2;
}


static double B_m(double T, double psi_w)
{
    // Bm has units of m^3/mol
    double Tj,tau_Air,tau_Water,B_aa,B_ww,B_aw;
    // NDG for fluid EOS for virial terms
    Tj=132.6312;
    tau_Air=Tj/T;
    tau_Water=Water.reduce.T/T;
    if (FlagUseVirialCorrelations==1)
    {
        B_aa=-0.000721183853646 +1.142682674467e-05*T -8.838228412173e-08*pow(T,2) 
        +4.104150642775e-10*pow(T,3) -1.192780880645e-12*pow(T,4) +2.134201312070e-15*pow(T,5) 
        -2.157430412913e-18*pow(T,6) +9.453830907795e-22*pow(T,7);
        B_ww=-10.8963128394 +2.439761625859e-01*T -2.353884845100e-03*pow(T,2) 
        +1.265864734412e-05*pow(T,3) -4.092175700300e-08*pow(T,4) +7.943925411344e-11*pow(T,5) 
        -8.567808759123e-14*pow(T,6) +3.958203548563e-17*pow(T,7);
    }
    else
    {
        B_aa=B_Air(T)*MM_Air()/1e3; //[m^3/kg] to [m^3/mol]
        B_ww=B_Water(T)*MM_Water()/1e3; //[m^3/kg] to [m^3/mol]
    }
    
    B_aw=_B_aw(T)/1e3; //[dm^3/mol] to [m^3/mol]
    return pow(1-psi_w,2)*B_aa+2*(1-psi_w)*psi_w*B_aw+psi_w*psi_w*B_ww;
}

static double dB_m_dT(double T, double psi_w)
{
    //dBm_dT has units of m^3/mol/K
    double Tj,tau_Air,tau_Water,dB_dT_aa,dB_dT_ww,dB_dT_aw;
    // NDG for fluid EOS for virial terms
    Tj=132.6312;
    tau_Air=Tj/T;
    tau_Water=Water.reduce.T/T;
    if (FlagUseVirialCorrelations)
    {
        dB_dT_aa=1.65159324353e-05 -3.026130954749e-07*T +2.558323847166e-09*pow(T,2) -1.250695660784e-11*pow(T,3) +3.759401946106e-14*pow(T,4) -6.889086380822e-17*pow(T,5) +7.089457032972e-20*pow(T,6) -3.149942145971e-23*pow(T,7);
        dB_dT_ww=0.65615868848 -1.487953162679e-02*T +1.450134660689e-04*pow(T,2) -7.863187630094e-07*pow(T,3) +2.559556607010e-09*pow(T,4) -4.997942221914e-12*pow(T,5) +5.417678681513e-15*pow(T,6) -2.513856275241e-18*pow(T,7);
    }
    else
    {
        dB_dT_aa=dBdT_Air(T)*MM_Air()/1e3; //[m^3/kg] to [m^3/mol]
        dB_dT_ww=dBdT_Water(T)*MM_Water()/1e3; //[m^3/kg] to [m^3/mol]
    }
        dB_dT_aw=_dB_aw_dT(T)/1e3; //[dm^3/mol] to [m^3/mol]
    return pow(1-psi_w,2)*dB_dT_aa+2*(1-psi_w)*psi_w*dB_dT_aw+psi_w*psi_w*dB_dT_ww;
}

static double C_m(double T, double psi_w)
{
    // Cm has units of m^6/mol^2
    double Tj,tau_Air,tau_Water,C_aaa,C_www,C_aww,C_aaw;
    // NDG for fluid EOS for virial terms
    Tj=132.6312;
    tau_Air=Tj/T;
    tau_Water=Water.reduce.T/T;
    if (FlagUseVirialCorrelations)
    {
        C_aaa=1.29192158975e-08 -1.776054020409e-10*T +1.359641176409e-12*pow(T,2) 
        -6.234878717893e-15*pow(T,3) +1.791668730770e-17*pow(T,4) -3.175283581294e-20*pow(T,5) 
        +3.184306136120e-23*pow(T,6) -1.386043640106e-26*pow(T,7);
        C_www=-0.580595811134 +1.365952762696e-02*T -1.375986293288e-04*pow(T,2) 
        +7.687692259692e-07*pow(T,3) -2.571440816920e-09*pow(T,4) +5.147432221082e-12*pow(T,5) 
        -5.708156494894e-15*pow(T,6) +2.704605721778e-18*pow(T,7);
    }
    else
    {
        C_aaa=C_Air(T)*MM_Air()*MM_Air()/1e6; //[m^6/kg^2] to [m^6/mol^2]
        C_www=C_Water(T)*MM_Water()*MM_Water()/1e6; //[m^6/kg^2] to [m^6/mol^2]
    }
    C_aaw=_C_aaw(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
    C_aww=_C_aww(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
    return pow(1-psi_w,3)*C_aaa+3*pow(1-psi_w,2)*psi_w*C_aaw+3*(1-psi_w)*psi_w*psi_w*C_aww+pow(psi_w,3)*C_www;
}

static double dC_m_dT(double T, double psi_w)
{
    // dCm_dT has units of m^6/mol^2/K
    
    double Tj,tau_Air,tau_Water,dC_dT_aaa,dC_dT_www,dC_dT_aww,dC_dT_aaw;
    // NDG for fluid EOS for virial terms
    Tj=132.6312;
    tau_Air=Tj/T;
    tau_Water=Water.reduce.T/T;
    if (FlagUseVirialCorrelations)
    {
        dC_dT_aaa=-2.46582342273e-10 +4.425401935447e-12*T -3.669987371644e-14*pow(T,2) +1.765891183964e-16*pow(T,3) -5.240097805744e-19*pow(T,4) +9.502177003614e-22*pow(T,5) -9.694252610339e-25*pow(T,6) +4.276261986741e-28*pow(T,7);
        dC_dT_www=0.0984601196142 -2.356713397262e-03*T +2.409113323685e-05*pow(T,2) -1.363083778715e-07*pow(T,3) +4.609623799524e-10*pow(T,4) -9.316416405390e-13*pow(T,5) +1.041909136255e-15*pow(T,6) -4.973918480607e-19*pow(T,7);
    }
    else
    {
        dC_dT_aaa=dCdT_Air(T)*MM_Air()*MM_Air()/1e6; //[m^6/kg^2] to [m^6/mol^2]
        dC_dT_www=dCdT_Water(T)*MM_Water()*MM_Water()/1e6; //[m^6/kg^2] to [m^6/mol^2]
    }
    dC_dT_aaw=_dC_aaw_dT(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
    dC_dT_aww=_dC_aww_dT(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
    return pow(1-psi_w,3)*dC_dT_aaa+3*pow(1-psi_w,2)*psi_w*dC_dT_aaw+3*(1-psi_w)*psi_w*psi_w*dC_dT_aww+pow(psi_w,3)*dC_dT_www;
}
double HumidityRatio(double psi_w)
{
    return psi_w*epsilon/(1-psi_w);
}

static double HenryConstant(double T)
{
    // Result has units of 1/Pa
    double p_ws,beta_N2,beta_O2,beta_Ar,beta_a,tau,Tr,Tc=647.096;
    Tr=T/Tc; 
    tau=1-Tr;
    p_ws=Props("P","T",T,"Q",1.0,"Water"); //[kPa]
    beta_N2=p_ws*exp(-9.67578/Tr+4.72162*pow(tau,0.355)/Tr+11.70585*pow(Tr,-0.41)*exp(tau));
    beta_O2=p_ws*exp(-9.44833/Tr+4.43822*pow(tau,0.355)/Tr+11.42005*pow(Tr,-0.41)*exp(tau));
    beta_Ar=p_ws*exp(-8.40954/Tr+4.29587*pow(tau,0.355)/Tr+10.52779*pow(Tr,-0.41)*exp(tau));
    beta_a=1/(0.7812/beta_N2+0.2095/beta_O2+0.0093/beta_Ar);
    return 1/(1.01325*beta_a)/1000.0;
}

double f_factor(double T, double p)
{
    double f=0,Rbar=8.314371,eps=1e-8,Tj;
    double x1=0,x2=0,x3,y1=0,y2,change=_HUGE;
    int iter=1;
    double p_ws,tau_Air,tau_Water,B_aa,B_aw,B_ww,C_aaa,C_aaw,C_aww,C_www,
        line1,line2,line3,line4,line5,line6,line7,line8,k_T,beta_H,LHS,RHS,psi_ws,
        vbar_ws;

    // Get total pressure in Pa from kPa
    p*=1000;
    
    // Saturation pressure [Pa]
    if (T>273.16)
    {
        // It is liquid water
        p_ws=PropsSI("P","T",T,"Q",0,"Water");
        if (FlagUseIsothermCompressCorrelation)
        {
            k_T = 1.6261876614E-22*pow(T,6) - 3.3016385196E-19*pow(T,5) + 2.7978984577E-16*pow(T,4)
                - 1.2672392901E-13*pow(T,3) + 3.2382864853E-11*pow(T,2) - 4.4318979503E-09*T + 2.5455947289E-07;
        }
        else
        {
            double rho = PropsSI("D","T",T,"P",p,"Water");
            k_T=DerivTerms((char *)"IsothermalCompressibility",T,rho,(char *)"Water")/1000; //[1/Pa]
        }
        beta_H=HenryConstant(T); //[1/Pa]
        vbar_ws=1.0/Props("D","T",T,"Q",0,"Water")*MM_Water()/1000; //[m^3/mol]
    }
    else
    {
        // It is ice
        p_ws=psub_Ice(T)*1000; // [Pa]
        k_T=IsothermCompress_Ice(T,p/1000); //[1/Pa]
        beta_H=0;
        vbar_ws=dg_dp_Ice(T,p/1000)*MM_Water()/1000/1000; //[m^3/mol]
    }
    
    // Hermann: In the iteration process of the enhancement factor in Eq. (3.25), k_T is set to zero for pw,s (T) > p.
    if (p_ws>p)
    {
        k_T=0;
        beta_H=0;
    }
    
    // NDG for fluid EOS for virial terms
    Tj=132.6312;
    tau_Air=Tj/T;
    tau_Water=Water.reduce.T/T;
    if (FlagUseVirialCorrelations)
    {
        B_aa=-0.000721183853646 +1.142682674467e-05*T -8.838228412173e-08*pow(T,2) 
        +4.104150642775e-10*pow(T,3) -1.192780880645e-12*pow(T,4) +2.134201312070e-15*pow(T,5) 
        -2.157430412913e-18*pow(T,6) +9.453830907795e-22*pow(T,7);
        B_ww=-10.8963128394 +2.439761625859e-01*T -2.353884845100e-03*pow(T,2) 
        +1.265864734412e-05*pow(T,3) -4.092175700300e-08*pow(T,4) +7.943925411344e-11*pow(T,5) 
        -8.567808759123e-14*pow(T,6) +3.958203548563e-17*pow(T,7);
        C_aaa=1.29192158975e-08 -1.776054020409e-10*T +1.359641176409e-12*pow(T,2) 
        -6.234878717893e-15*pow(T,3) +1.791668730770e-17*pow(T,4) -3.175283581294e-20*pow(T,5) 
        +3.184306136120e-23*pow(T,6) -1.386043640106e-26*pow(T,7);
        C_www=-0.580595811134 +1.365952762696e-02*T -1.375986293288e-04*pow(T,2) 
        +7.687692259692e-07*pow(T,3) -2.571440816920e-09*pow(T,4) +5.147432221082e-12*pow(T,5) 
        -5.708156494894e-15*pow(T,6) +2.704605721778e-18*pow(T,7);
    }
    else
    {
        B_aa=B_Air(T)*MM_Air()/1e3; //[m^3/kg] to [m^3/mol]
        C_aaa=C_Air(T)*MM_Air()*MM_Air()/1e6; //[m^6/kg^2] to [m^6/mol^2]
        B_ww=B_Water(T)*MM_Water()/1e3; //[m^3/kg] to [m^3/mol]
        C_www=C_Water(T)*MM_Water()*MM_Water()/1e6; //[m^6/kg^2] to [m^6/mol^2]
    }
    B_aw=_B_aw(T)/1e3; //[dm^3/mol] to [m^3/mol]
    C_aaw=_C_aaw(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
    C_aww=_C_aww(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
    
    // Use a little secant loop to find f iteratively
    // Start out with a guess value of 1 for f
    while ((iter<=3 || change>eps) && iter<100)
    {
        if (iter==1){x1=1.00; f=x1;}
        if (iter==2){x2=1.00+0.000001; f=x2;}
        if (iter>2) {f=x2;}
        
            // Left-hand-side of Equation 3.25
            LHS=log(f);
            // Eqn 3.24
            psi_ws=f*p_ws/p;
            
            // All the terms forming the RHS of Eqn 3.25
            line1=((1+k_T*p_ws)*(p-p_ws)-k_T*(p*p-p_ws*p_ws)/2.0)/(Rbar*T)*vbar_ws+log(1-beta_H*(1-psi_ws)*p);
            line2=pow(1-psi_ws,2)*p/(Rbar*T)*B_aa-2*pow(1-psi_ws,2)*p/(Rbar*T)*B_aw-(p-p_ws-pow(1-psi_ws,2)*p)/(Rbar*T)*B_ww;
            line3=pow(1-psi_ws,3)*p*p/pow(Rbar*T,2)*C_aaa+(3*pow(1-psi_ws,2)*(1-2*(1-psi_ws))*p*p)/(2*pow(Rbar*T,2))*C_aaw;
            line4=-3*pow(1-psi_ws,2)*psi_ws*p*p/pow(Rbar*T,2)*C_aww-((3-2*psi_ws)*psi_ws*psi_ws*p*p-p_ws*p_ws)/(2*pow(Rbar*T,2))*C_www;
            line5=-(pow(1-psi_ws,2)*(-2+3*psi_ws)*psi_ws*p*p)/pow(Rbar*T,2)*B_aa*B_ww;
            line6=-(2*pow(1-psi_ws,3)*(-1+3*psi_ws)*p*p)/pow(Rbar*T,2)*B_aa*B_aw;
            line7=(6*pow(1-psi_ws,2)*psi_ws*psi_ws*p*p)/pow(Rbar*T,2)*B_ww*B_aw-(3*pow(1-psi_ws,4)*p*p)/(2*pow(Rbar*T,2))*B_aa*B_aa;
            line8=-(2*pow(1-psi_ws,2)*psi_ws*(-2+3*psi_ws)*p*p)/pow(Rbar*T,2)*B_aw*B_aw-(p_ws*p_ws-(4-3*psi_ws)*pow(psi_ws,3)*p*p)/(2*pow(Rbar*T,2))*B_ww*B_ww;
            RHS=line1+line2+line3+line4+line5+line6+line7+line8;
        
        if (iter==1){y1=LHS-RHS;}
        if (iter>1)
        {
            y2=LHS-RHS;
            x3=x2-y2/(y2-y1)*(x2-x1);
            change=fabs(y2/(y2-y1)*(x2-x1));
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
    }
    if (f>=1.0)
        return f;
    else
        return 1.0;
}
void HAHelp(void)
{
    printf("Sorry, Need to update!");
}
int returnHumAirCode(const char * Code)
{
    if (!strcmp(Code,"GIVEN_TDP"))
        return GIVEN_TDP;
    else if (!strcmp(Code,"GIVEN_HUMRAT"))
        return GIVEN_HUMRAT;
    else if (!strcmp(Code,"GIVEN_TWB"))
        return GIVEN_TWB;
    else if (!strcmp(Code,"GIVEN_RH"))
        return GIVEN_RH;
    else if (!strcmp(Code,"GIVEN_ENTHALPY"))
        return GIVEN_ENTHALPY;
    else
    {
        fprintf(stderr,"Code to returnHumAirCode in HumAir.c [%s] not understood",Code);
        return -1;
    }
}
double Viscosity(double T, double p, double psi_w)
{
    /*
    Using the method of:
    
    P.T. Tsilingiris, 2009, Thermophysical and transport properties of humid air at temperature range between 0 and 100 oC, Energy Conversion and Management, 49, 1098-1010
    
    but using the detailed measurements for pure fluid from IAPWS formulations
    */
    double mu_a,mu_w,Phi_av,Phi_va,Ma,Mw;
    Mw=MM_Water();
    Ma=MM_Air();
    // Viscosity of dry air at dry-bulb temp and total pressure
    mu_a=Props("V","T",T,"P",p,"Air");
    // Viscosity of pure saturated water at dry-bulb temperature
    mu_w=Props("V","P",p,"Q",1,"Water");
    Phi_av=sqrt(2.0)/4.0*pow(1+Ma/Mw,-0.5)*pow(1+sqrt(mu_a/mu_w)*pow(Mw/Ma,0.25),2); //[-]
    Phi_va=sqrt(2.0)/4.0*pow(1+Mw/Ma,-0.5)*pow(1+sqrt(mu_w/mu_a)*pow(Ma/Mw,0.25),2); //[-]
    return (1-psi_w)*mu_a/((1-psi_w)+psi_w*Phi_av)+psi_w*mu_w/(psi_w+(1-psi_w)*Phi_va);
}
double Conductivity(double T, double p, double psi_w)
{
    /*
    Using the method of:
    
    P.T. Tsilingiris, 2009, Thermophysical and transport properties of humid air at temperature range between 0 and 100 oC, Energy Conversion and Management, 49, 1098-1010
    
    but using the detailed measurements for pure fluid from IAPWS formulations
    */
    double mu_a,mu_w,k_a,k_w,Phi_av,Phi_va,Ma,Mw;
    Mw=MM_Water();
    Ma=MM_Air();
    // Viscosity of dry air at dry-bulb temp and total pressure
    k_a=Props("L","T",T,"P",p,"Air");
    mu_a=Props("V","T",T,"P",p,"Air");
    // Viscosity of pure saturated water vapor at dry-bulb temperature
    k_w=Props("L","P",p,"Q",1,"Water");
    mu_w=Props("V","P",p,"Q",1,"Water");
    Phi_av=sqrt(2.0)/4.0*pow(1+Ma/Mw,-0.5)*pow(1+sqrt(mu_a/mu_w)*pow(Mw/Ma,0.25),2); //[-]
    Phi_va=sqrt(2.0)/4.0*pow(1+Mw/Ma,-0.5)*pow(1+sqrt(mu_w/mu_a)*pow(Ma/Mw,0.25),2); //[-]
    return (1-psi_w)*k_a/((1-psi_w)+psi_w*Phi_av)+psi_w*k_w/(psi_w+(1-psi_w)*Phi_va);
}
double MolarVolume(double T, double p, double psi_w)
{
    // Output in m^3/mol
    int iter;
    double v_bar0, v_bar=0, R_bar=8.314472,x1=0,x2=0,x3,y1=0,y2,resid,eps,Bm,Cm;
    
    // -----------------------------
    // Iteratively find molar volume
    // -----------------------------
    
    // Start by assuming it is an ideal gas to get initial guess
    v_bar0=R_bar*T/p/1000;
    
    //Bring outside the loop since not a function of v_bar
    Bm=B_m(T,psi_w);
    Cm=C_m(T,psi_w);

    iter=1; eps=1e-8; resid=999;
    while ((iter<=3 || fabs(resid)>eps) && iter<100)
    {
        if (iter==1){x1=v_bar0; v_bar=x1;}
        if (iter==2){x2=v_bar0+0.000001; v_bar=x2;}
        if (iter>2) {v_bar=x2;}
        
            // factor of 1000 is to deal with kmol and mol conversion
            // want v_bar in m^3/mol and R_bar in kJ/mol-K
            resid=p-(R_bar/1000)*T/v_bar*(1+Bm/v_bar+Cm/(v_bar*v_bar));
        
        if (iter==1){y1=resid;}
        if (iter>1)
        {
            y2=resid;
            x3=x2-y2/(y2-y1)*(x2-x1);
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
    }
    return v_bar;
}
double IdealGasMolarEnthalpy_Water(double T, double v_bar)
{
    double hbar_w_0,tau,rhobar,hbar_w,rho;
    // Ideal-Gas contribution to enthalpy of water
    hbar_w_0=-0.01102303806;//[kJ/kmol]
    tau = Water.crit.T/T;
    rhobar = 1/v_bar; //[kmol/m^3]
    rho = rhobar * Water.params.molemass;
    hbar_w = hbar_w_0+R_bar*T*(1+tau*Water.dphi0_dTau(tau,rho/Water.crit.rho));
    return hbar_w;
}
double IdealGasMolarEntropy_Water(double T, double p)
{
    double sbar_w,tau,R_bar,rho;
    R_bar = 8.314371; //[kJ/kmol/K]
    tau = Water.crit.T/T;
    rho = p/(R_bar/MM_Water()*T); //[kg/m^3]
    sbar_w = R_bar*(tau*Water.dphi0_dTau(tau,rho/Water.crit.rho)-Water.phi0(tau,rho/Water.crit.rho)); //[kJ/kmol/K]
    return sbar_w; 
}
double IdealGasMolarEnthalpy_Air(double T, double v_bar)
{
    double hbar_a_0,tau,rhobar,hbar_a,R_bar_Lemmon, rho;
    // Ideal-Gas contribution to enthalpy of air
    hbar_a_0=-7914.149298; //[kJ/kmol]
    //Tj and rhoj are given by 132.6312 and 302.5507652 respectively
    tau=132.6312/T;
    rhobar=1/v_bar; //[kmol/m^3]
    rho = rhobar * Props1("Air","molemass");
    R_bar_Lemmon=8.314510; //[kJ/kmol/K]
    hbar_a=hbar_a_0+R_bar_Lemmon*T*(1+tau*DerivTerms((char *)"dphi0_dTau",T,rho,(char *)"Air")); //[kJ/kmol]
    return hbar_a;
}
double IdealGasMolarEntropy_Air(double T, double v_bar_a)
{
    double sbar_0_Lem,tau,sbar_a,R_bar_Lemmon,T0=273.15,p0=101.325,v_0,v_bar_0, rho_a,rho_bar_a, rho_bar_0,rho_0;
    R_bar_Lemmon=8.314510; //[kJ/kmol/K]
    // Ideal-Gas contribution to entropy of air
    sbar_0_Lem=-196.1375815; //[kJ/kmol/K]
    //Tj and rhoj are given by 132.6312 and 302.5507652 respectively
    tau=132.6312/T; //[no units]
    v_0 = R_bar_Lemmon/MM_Air()*T0/p0; //[m^3/kg]
    rho_bar_a = 1/v_bar_a;
    rho_a = rho_bar_a * Props1("Air","molemass");
    v_bar_0 = R_bar_Lemmon*T0/p0; //[m^3/kmol]
    rho_bar_0 = 1/v_bar_0;
    rho_0 = rho_bar_0 * Props1("Air","molemass");
    sbar_a=sbar_0_Lem+R_bar_Lemmon*(tau*DerivTerms((char *)"dphi0_dTau",T,rho_0,(char *)"Air")-DerivTerms((char *)"phi0",T,rho_0,(char *)"Air"))+R_bar_Lemmon*log(v_bar_a/v_bar_0); //[kJ/kmol/K]
    return sbar_a; //[kJ/kmol/K]
}

double MolarEnthalpy(double T, double p, double psi_w, double v_bar)
{
    // In units of kJ/kmol
    
    // vbar (molar volume) in m^3/kg
    
    double hbar_0,hbar_a,hbar_w,hbar,R_bar=8.314472;
    // ----------------------------------------
    //      Enthalpy
    // ----------------------------------------
    // Constant for enthalpy
    // Not clear why getting rid of this term yields the correct values in the table, but enthalpies are equal to an additive constant, so not a big deal
    hbar_0=0.0;//2.924425468; //[kJ/kmol]
    
    if (FlagUseIdealGasEnthalpyCorrelations)
    {
    hbar_w=2.7030251618E-03*T*T + 3.1994361015E+01*T + 3.6123174929E+04;
    hbar_a=9.2486716590E-04*T*T + 2.8557221776E+01*T - 7.8616129429E+03;
    }
    else
    {
    hbar_w=IdealGasMolarEnthalpy_Water(T,v_bar);
    hbar_a=IdealGasMolarEnthalpy_Air(T,v_bar);
    }
    
    hbar=hbar_0+(1-psi_w)*hbar_a+psi_w*hbar_w+R_bar*T*((B_m(T,psi_w)-T*dB_m_dT(T,psi_w))/v_bar+(C_m(T,psi_w)-T/2.0*dC_m_dT(T,psi_w))/(v_bar*v_bar));
    return hbar; //[kJ/kmol]
}
double MassEnthalpy(double T, double p, double psi_w)
{
    double v_bar = MolarVolume(T, p, psi_w); //[m^3/mol_ha]
    double h_bar = MolarEnthalpy(T, p, psi_w, v_bar); //[kJ/kmol_ha]
    double W = HumidityRatio(psi_w); //[kg_w/kg_da]
    double M_ha = MM_Water()*psi_w+(1-psi_w)*28.966;
    return h_bar*(1+W)/M_ha; //[kJ/kg_da]
}

double MolarEntropy(double T, double p, double psi_w, double v_bar)
{
    // In units of kJ/kmol/K

    // Serious typo in RP-1485 - should use total pressure rather than
    // reference pressure in density calculation for water vapor molar entropy
    
    // vbar (molar volume) in m^3/kmol
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,f=999,R_bar_Lem=8.314510;
    int iter=1;
    double sbar_0,sbar_a=0,sbar_w=0,sbar,R_bar=8.314472,vbar_a_guess, Baa, Caaa,vbar_a=0;
    double B,dBdT,C,dCdT;
    // Constant for entropy
    sbar_0=0.02366427495;  //[kJ/kmol/K]

    //Calculate vbar_a, the molar volume of dry air
    // B_m, C_m, etc. functions take care of the units
    Baa = B_m(T,0); 
    B = B_m(T,psi_w);
    dBdT = dB_m_dT(T,psi_w);
    Caaa = C_m(T,0);
    C = C_m(T,psi_w);
    dCdT = dC_m_dT(T,psi_w);

    vbar_a_guess = R_bar_Lem*T/p; //[m^3/mol] since p in [kPa]
    
    while ((iter<=3 || fabs(f)>eps) && iter<100)
    {
        if (iter==1){x1=vbar_a_guess; vbar_a=x1;}
        if (iter==2){x2=vbar_a_guess+0.001; vbar_a=x2;}
        if (iter>2) {vbar_a=x2;}
            f=R_bar_Lem*T/vbar_a*(1+Baa/vbar_a+Caaa/pow(vbar_a,2))-p;
        if (iter==1){y1=f;}
        if (iter>1)
        {
            y2=f;
            x3=x2-y2/(y2-y1)*(x2-x1);
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
        if (iter>100){ return _HUGE; }
    }
    
    if (FlagUseIdealGasEnthalpyCorrelations)
    {
        std::cout << "Not implemented" << std::endl;
    }
    else
    {
        sbar_w=IdealGasMolarEntropy_Water(T,p);
        sbar_a=IdealGasMolarEntropy_Air(T,vbar_a);
    }
    if (psi_w!=0)
    {
        sbar=sbar_0+(1-psi_w)*sbar_a+psi_w*sbar_w-R_bar*(
        (B+T*dBdT)/v_bar+(C+T*dCdT)/(2*pow(v_bar,2))+
        (1-psi_w)*log(1-psi_w)+psi_w*log(psi_w));
    }
    else{
        sbar=sbar_0+sbar_a;
    }
    return sbar; //[kJ/kmol/K]
}

double DewpointTemperature(double T, double p, double psi_w)
{
    int iter;
    double p_w,eps,resid,Tdp=0,x1=0,x2=0,x3,y1=0,y2,T0;
    double p_ws_dp,f_dp;
    
    // Make sure it isn't dry air, return an impossible temperature otherwise
    if ((1-psi_w)<1e-16)
    {
        return -1;
    }
    // ------------------------------------------
    // Iteratively find the dewpoint temperature
    // ------------------------------------------

    // The highest dewpoint temperature possible is the dry-bulb temperature.  
    // When they are equal, the air is saturated (R=1)
    
    p_w = psi_w*p;

    // 0.61165... is the triple point pressure of water in kPa
    if (p_w > 0.6116547241637944){
        T0 = Props("T","P",p_w,"Q",1.0,"Water");
    }
    else{
        T0 = 268;
    }
    // A good guess for Tdp is that enhancement factor is unity, which yields
    // p_w_s = p_w, and get guess for T from saturation temperature
    
    iter=1; eps=1e-8; resid=999;
    while ((iter<=3 || fabs(resid)>eps) && iter<100)
    {
        if (iter==1){x1 = T0; Tdp=x1;}
        if (iter==2){x2 = x1 + 0.1; Tdp=x2;}
        if (iter>2) {Tdp=x2;}
        
            if (Tdp >= 273.16)
            {
                // Saturation pressure at dewpoint [kPa]
                p_ws_dp=Props("P","T",Tdp,"Q",0,"Water");
            }
            else
            {
                // Sublimation pressure at icepoint [kPa]
                p_ws_dp=psub_Ice(Tdp);
            }
            // Enhancement Factor at dewpoint temperature [-]
            f_dp=f_factor(Tdp,p);
            // Error between target and actual pressure [kPa]
            resid=p_w-p_ws_dp*f_dp;
        
        if (iter==1){y1=resid;}
        if (iter>1)
        {
            y2=resid;
            x3=x2-y2/(y2-y1)*(x2-x1);
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
    }
    return Tdp;
}

class WetBulbSolver : public FuncWrapper1D
{
private:
    double _T,_p,_W,LHS,v_bar_w,M_ha;
public:
    WetBulbSolver(double T, double p, double psi_w){
        _T = T;
        _p = p;
        _W = epsilon*psi_w/(1-psi_w);

        //These things are all not a function of Twb
        v_bar_w = MolarVolume(T,p,psi_w);
        M_ha = MM_Water()*psi_w+(1-psi_w)*28.966;
        LHS = MolarEnthalpy(T,p,psi_w,v_bar_w)*(1+_W)/M_ha;
    };
    ~WetBulbSolver(){};
    double call(double Twb)
    {
        double epsilon=0.621945;
        double f_wb,p_ws_wb,p_s_wb,W_s_wb,h_w,M_ha_wb,psi_wb,v_bar_wb;

        // Enhancement Factor at wetbulb temperature [-]
        f_wb=f_factor(Twb,_p);
        if (Twb > 273.16)
        {
            // Saturation pressure at wetbulb temperature [kPa]
            p_ws_wb=Props("P","T",Twb,"Q",0,"Water");
        }
        else
        {
            // Sublimation pressure at wetbulb temperature [kPa]
            p_ws_wb=psub_Ice(Twb);
        }
            
        // Vapor pressure
        p_s_wb = f_wb*p_ws_wb;
        // wetbulb humidity ratio
        W_s_wb = epsilon*p_s_wb/(_p-p_s_wb);
        // wetbulb water mole fraction
        psi_wb = W_s_wb/(epsilon+W_s_wb);
        if (Twb > 273.16)
        {
            // Enthalpy of water [kJ/kg_water]
            h_w=Props("H","T",Twb,"P",_p,"Water");
        }
        else
        {
            // Enthalpy of ice [kJ/kg_water]
            h_w=h_Ice(Twb,_p)/1000;
        }
        // Mole masses of wetbulb and humid air
        
        M_ha_wb=MM_Water()*psi_wb+(1-psi_wb)*28.966;
        v_bar_wb=MolarVolume(Twb,_p,psi_wb);
        double RHS = (MolarEnthalpy(Twb,_p,psi_wb,v_bar_wb)*(1+W_s_wb)/M_ha_wb+(_W-W_s_wb)*h_w);
        if (!ValidNumber(LHS-RHS)){throw ValueError();}
        return LHS - RHS;
    }
};

class WetBulbTminSolver : public FuncWrapper1D
{
public:
    double p,hair_dry,r, RHS;
    WetBulbTminSolver(double p, double hair_dry){
        this->p = p;
        this->hair_dry = hair_dry;
    };
    ~WetBulbTminSolver(){};
    double call(double Ts)
    {
        RHS = HAProps("H","T",Ts,"P",p,"R",1);
        if (!ValidNumber(RHS)){throw ValueError();}
        r = RHS - this->hair_dry;
        return r;
    }
};

double WetbulbTemperature(double T, double p, double psi_w)
{
    // ------------------------------------------
    // Iteratively find the wetbulb temperature
    // ------------------------------------------
    //  
    // If the temperature is less than the saturation temperature of water
    // for the given atmospheric pressure, the highest wetbulb temperature that is possible is the dry bulb
    // temperature
    //
    // If the temperature is above the saturation temperature corresponding to the atmospheric pressure,
    // then the maximum value for the wetbulb temperature is the saturation temperature
    double Tmax = T;
    double Tsat = Props("T","P",p,"Q",1.0,"Water");
    if (T >= Tsat)
    {
        Tmax = Tsat;
    }

    // Instantiate the solver container class
    WetBulbSolver WBS(T,p,psi_w);

    std::string errstr;
    
    double return_val;
    try{
        return_val = Secant(&WBS,Tmax,0.0001,1e-8,50,&errstr);
        
        // Solution obtained is out of range (T>Tmax)
        if (return_val > Tmax) {throw ValueError();}
    }
    catch(std::exception &)
    {
        // The lowest wetbulb temperature that is possible for a given dry bulb temperature 
        // is the saturated air temperature which yields the enthalpy of dry air at dry bulb temperature

        try{
            double hair_dry = MassEnthalpy(T,p,0);

            // Directly solve for the saturated temperature that yields the enthalpy desired
            WetBulbTminSolver WBTS(p,hair_dry);
            double Tmin = Brent(&WBTS,210,Tsat-1,1e-12,1e-12,50,&errstr);

            return_val = Brent(&WBS,Tmin-30,Tmax-1,1e-12,1e-12,50,&errstr);
        }
        catch(std::exception)
        {
            return_val = _HUGE;
        }
    }
    return return_val;	
}
static int Name2Type(const char *Name)
{
    if (!strcmp(Name,"Omega") || !strcmp(Name,"HumRat") || !strcmp(Name,"W"))
        return GIVEN_HUMRAT;
    else if (!strcmp(Name,"Tdp") || !strcmp(Name,"T_dp") || !strcmp(Name,"DewPoint") || !strcmp(Name,"D"))
        return GIVEN_TDP;
    else if (!strcmp(Name,"Twb") || !strcmp(Name,"T_wb") || !strcmp(Name,"WetBulb") || !strcmp(Name,"B"))
        return GIVEN_TWB;
    else if (!strcmp(Name,"Enthalpy") || !strcmp(Name,"H"))
        return GIVEN_ENTHALPY;
    else if (!strcmp(Name,"RH") || !strcmp(Name,"RelHum") || !strcmp(Name,"R"))
        return GIVEN_RH;
    else if (!strcmp(Name,"Tdb") || !strcmp(Name,"T_db") || !strcmp(Name,"T"))
        return GIVEN_T;
    else if (!strcmp(Name,"P"))
        return GIVEN_P;
    else if (!strcmp(Name,"V") || !strcmp(Name,"Vda"))
        return GIVEN_V;
    else if (!strcmp(Name,"mu") || !strcmp(Name,"Visc") || !strcmp(Name,"M"))
        return GIVEN_VISC;
    else if (!strcmp(Name,"k") || !strcmp(Name,"Conductivity") || !strcmp(Name,"K"))
        return GIVEN_COND;
    else
        printf("Sorry, your input [%s] was not understood to Name2Type in HumAir.c. Acceptable values are T,P,R,W,D,B,H,M,K and aliases thereof\n",Name);
        return -1;
}
int TypeMatch(int TypeCode, const char *Input1Name, const char *Input2Name, const char *Input3Name)
{
    // Return the index of the input variable that matches the input, otherwise return -1 for failure
    if (TypeCode==Name2Type(Input1Name))
        return 1;
    if (TypeCode==Name2Type(Input2Name))
        return 2;
    if (TypeCode==Name2Type(Input3Name))
        return 3;
    else
        return -1;
}
double MoleFractionWater(double T, double p, int HumInput, double InVal)
{
    double p_ws,f,W,epsilon=0.621945,Tdp,p_ws_dp,f_dp,p_w_dp,p_s,RH;
    
    if (HumInput==GIVEN_HUMRAT) //(2)
    {
        W=InVal;
        return W/(epsilon+W);
    }
    else if (HumInput==GIVEN_RH)
    {
        if (T>=273.16)
        {
            // Saturation pressure [kPa]
            p_ws=Props("P","T",T,"Q",0,"Water");
        }
        else
        {
            // Sublimation pressure [kPa]
            p_ws=psub_Ice(T);
            
        }
        // Enhancement Factor [-]
        f=f_factor(T,p);

        // Saturation pressure [kPa]
        p_s=f*p_ws;
        RH=InVal;
        W=epsilon*RH*p_s/(p-RH*p_s);
        return W/(epsilon+W);
    }
    else if (HumInput==GIVEN_TDP)
    {
        Tdp=InVal;
        // Saturation pressure at dewpoint [kPa]
        if (Tdp>=273.16)
        {
            p_ws_dp=Props("P","T",Tdp,"Q",0,"Water");
        }
        else{
            // Sublimation pressure [kPa]
            p_ws_dp=psub_Ice(Tdp);
        }

        // Enhancement Factor at dewpoint temperature [-]
        f_dp=f_factor(Tdp,p);
        // Water vapor pressure at dewpoint [kPa]
        p_w_dp=f_dp*p_ws_dp;
        // Water mole fraction [-]
        return p_w_dp/p;
    }
    else
    {
        return -1000000;
    }
}

double RelativeHumidity(double T, double p, double psi_w)
{
    double p_ws,f,p_s,W;
    if (T>=273.16)
    {
        // Saturation pressure [kPa]
        p_ws=Props("P","T",T,"Q",0,"Water");
    }
    else
    {
        // sublimation pressure [kPa]
        p_ws=psub_Ice(T);
        
    }
    // Enhancement Factor [-]
    f=f_factor(T,p);

    // Saturation pressure [kPa]
    p_s=f*p_ws;
    // Find humidity ratio
    W=HumidityRatio(psi_w);
    // Find relative humidity using W/e=phi*p_s/(p-phi*p_s)
    return W/epsilon*p/(p_s*(1+W/epsilon));
}
EXPORT_CODE double CONVENTION HAProps(const char *OutputName, const char *Input1Name, double Input1, const char *Input2Name, double Input2, const char *Input3Name, double Input3)
{
    try
    {
        int In1Type, In2Type, In3Type,iT,iW,iTdp,iRH,ip,Type1,Type2;
        double vals[3],p,T,RH,W,Tdp,psi_w,M_ha,v_bar,h_bar,s_bar,MainInputValue,SecondaryInputValue,T_guess;
        double Value1,Value2,W_guess;
        char MainInputName[100], SecondaryInputName[100],Name1[100],Name2[100];
        
        vals[0]=Input1;
        vals[1]=Input2;
        vals[2]=Input3;
        
        // First figure out what kind of inputs you have, convert names to Macro expansions
        In1Type=Name2Type(Input1Name);
        In2Type=Name2Type(Input2Name);
        In3Type=Name2Type(Input3Name);
        
        // Pressure must be included
        ip=TypeMatch(GIVEN_P,Input1Name,Input2Name,Input3Name);
        if (ip>0)
            p=vals[ip-1];
        else
            return -1000;

        // -----------------------------------------------------------------------------------------------------
        // Check whether the remaining values give explicit solution for mole fraction of water - nice and fast
        // -----------------------------------------------------------------------------------------------------
        
        // Find the codes if they are there
        iT=       TypeMatch(GIVEN_T,Input1Name,Input2Name,Input3Name);
        iRH=     TypeMatch(GIVEN_RH,Input1Name,Input2Name,Input3Name);
        iW=  TypeMatch(GIVEN_HUMRAT,Input1Name,Input2Name,Input3Name);
        iTdp=   TypeMatch(GIVEN_TDP,Input1Name,Input2Name,Input3Name);
        
        if (iT>0) // Found T (or alias) as an input
        {
            T=vals[iT-1];
            if (iRH>0) //Relative Humidity is provided
            {
                RH=vals[iRH-1];
                psi_w=MoleFractionWater(T,p,GIVEN_RH,RH);
            }
            else if (iW>0)
            {
                W=vals[iW-1];
                psi_w=MoleFractionWater(T,p,GIVEN_HUMRAT,W);
            }
            else if (iTdp>0)
            {
                Tdp=vals[iTdp-1];
                psi_w=MoleFractionWater(T,p,GIVEN_TDP,Tdp);
            }
            else
            {
                // Temperature and pressure are known, figure out which variable holds the other value
                if (In1Type!=GIVEN_T && In1Type!=GIVEN_P)
                {
                    strcpy(SecondaryInputName,Input1Name);
                    SecondaryInputValue=Input1;
                }
                else if (In2Type!=GIVEN_T && In2Type!=GIVEN_P)
                {
                    strcpy(SecondaryInputName,Input2Name);
                    SecondaryInputValue=Input2;
                }
                else if (In3Type!=GIVEN_T && In3Type!=GIVEN_P)
                {
                    strcpy(SecondaryInputName,Input3Name);
                    SecondaryInputValue=Input3;
                }
                else{
                    return _HUGE;
                }
                // Find the value for W
                W_guess=0.0001;
                W=Secant_HAProps_W(SecondaryInputName,"P",p,"T",T,SecondaryInputValue,W_guess);
                // Mole fraction of water
                psi_w=MoleFractionWater(T,p,GIVEN_HUMRAT,W);
                // And on to output...
            }
        }
        else
        {
            // Need to iterate to find dry bulb temperature since temperature is not provided
            
            // Pick one input, and alter T to match the other input
            
            // Get the variables and their values that are NOT pressure for simplicity
            // because you know you need pressure as an input and you already have
            // its value in variable p
            if (ip==1) // Pressure is in slot 1
            {
                strcpy(Name1,Input2Name);
                Value1=Input2;
                strcpy(Name2,Input3Name);
                Value2=Input3;
            }
            else if (ip==2) // Pressure is in slot 2
            {
                strcpy(Name1,Input1Name);
                Value1=Input1;
                strcpy(Name2,Input3Name);
                Value2=Input3;
            }
            else if (ip==3) // Pressure is in slot 3
            {
                strcpy(Name1,Input1Name);
                Value1=Input1;
                strcpy(Name2,Input2Name);
                Value2=Input2;
            }
            else{
            return _HUGE;
            }
                
            // Get the integer type codes
            Type1=Name2Type(Name1);
            Type2=Name2Type(Name2);
            
            // First, if one of the inputs is something that can potentially yield 
            // an explicit solution at a given iteration of the solver, use it
            if (Type1==GIVEN_RH || Type1==GIVEN_HUMRAT || Type1==GIVEN_TDP)
            {
                // First input variable is a "nice" one
                
                // MainInput is the one that you are using in the call to HAProps
                MainInputValue=Value1;
                strcpy(MainInputName,Name1);
                // SecondaryInput is the one that you are trying to match
                SecondaryInputValue=Value2;
                strcpy(SecondaryInputName,Name2);
            }
            else if (Type2==GIVEN_RH || Type2==GIVEN_HUMRAT || Type2==GIVEN_TDP)
            {
                // Second input variable is a "nice" one
                
                // MainInput is the one that you are using in the call to HAProps
                MainInputValue=Value2;
                strcpy(MainInputName,Name2);
                // SecondaryInput is the one that you are trying to match
                SecondaryInputValue=Value1;
                strcpy(SecondaryInputName,Name1);
            }
            else
            {
                printf("Sorry, but currently at least one of the variables as an input to HAProps() must be temperature, relative humidity, humidity ratio, or dewpoint\n  Eventually will add a 2-D NR solver to find T and psi_w simultaneously, but not included now\n");
                return -1000;
            }

            double T_min = 210;
            double T_max = 450;

            T = -1;

            // First try to use the secant solver to find T at a few different temperatures
            for (T_guess = 210; T_guess < 450; T_guess += 60)
            {
                try{
                    T = Secant_HAProps_T(SecondaryInputName,(char *)"P",p,MainInputName,MainInputValue,SecondaryInputValue,T_guess);
                    double val = HAProps(SecondaryInputName,(char *)"T",T,(char *)"P",p,MainInputName,MainInputValue);
                    if (!ValidNumber(T) || !ValidNumber(val) || !(T_min < T && T < T_max) || fabs(val-SecondaryInputValue)>1e-6)
                    { 
                        throw ValueError(); 
                    }
                    else
                    {
                        break;
                    }
                }
                catch (std::exception &){};
            }
            
            if (T < 0) // No solution found using secant
            {
                // Use the Brent's method solver to find T
                T = Brent_HAProps_T(SecondaryInputName,(char *)"P",p,MainInputName,MainInputValue,SecondaryInputValue,T_min,T_max);
            }
            
            // If you want the temperature, return it
            if (Name2Type(OutputName)==GIVEN_T)
                return T;
            else
            {
                // Otherwise, find psi_w for further calculations in the following section
                W=HAProps((char *)"W",(char *)"T",T,(char *)"P",p,MainInputName,MainInputValue);
                psi_w=MoleFractionWater(T,p,GIVEN_HUMRAT,W);
            }
        }
        M_ha=(1-psi_w)*28.966+MM_Water()*psi_w; //[kg_ha/kmol_ha]
        
        // -----------------------------------------------------------------
        // Calculate and return the desired value for known set of T,p,psi_w
        // -----------------------------------------------------------------
        if (!strcmp(OutputName,"Vda") || !strcmp(OutputName,"V"))
        {
            v_bar=MolarVolume(T,p,psi_w); //[m^3/mol_ha]
            W=HumidityRatio(psi_w);
            return v_bar*(1+W)/M_ha*1000; //[m^3/kg_da]
        }
        else if (!strcmp(OutputName,"Vha"))
        {
            v_bar=MolarVolume(T,p,psi_w); //[m^3/mol_ha]
            return v_bar/M_ha*1000; //[m^3/kg_ha]
        }
        else if (!strcmp(OutputName,"Y"))
        {
            return psi_w; //[mol_w/mol]
        }
        else if (!strcmp(OutputName,"Hda") || !strcmp(OutputName,"H"))
        {
            return MassEnthalpy(T,p,psi_w);
        }
        else if (!strcmp(OutputName,"Hha"))
        {
            v_bar=MolarVolume(T,p,psi_w); //[m^3/mol_ha]
            h_bar=MolarEnthalpy(T,p,psi_w,v_bar); //[kJ/kmol_ha]
            return h_bar/M_ha; //[kJ/kg_ha]
        }
        else if (!strcmp(OutputName,"S") || !strcmp(OutputName,"Entropy"))
        {
            v_bar=MolarVolume(T,p,psi_w); //[m^3/mol_ha]
            s_bar=MolarEntropy(T,p,psi_w,v_bar); //[kJ/kmol_ha]
            W=HumidityRatio(psi_w); //[kg_w/kg_da]
            return s_bar*(1+W)/M_ha; //[kJ/kg_da]
        }
        else if (!strcmp(OutputName,"C") || !strcmp(OutputName,"cp"))
        {
            double v_bar1,v_bar2,h_bar1,h_bar2, cp_bar, dT = 1e-3;
            v_bar1=MolarVolume(T-dT,p,psi_w); //[m^3/mol_ha]
            h_bar1=MolarEnthalpy(T-dT,p,psi_w,v_bar1); //[kJ/kmol_ha]
            v_bar2=MolarVolume(T+dT,p,psi_w); //[m^3/mol_ha]
            h_bar2=MolarEnthalpy(T+dT,p,psi_w,v_bar2); //[kJ/kmol_ha]
            W=HumidityRatio(psi_w); //[kg_w/kg_da]
            cp_bar = (h_bar2-h_bar1)/(2*dT);
            return cp_bar*(1+W)/M_ha; //[kJ/kg_da]
        }
        else if (!strcmp(OutputName,"Cha") || !strcmp(OutputName,"cp_ha"))
        {
            double v_bar1,v_bar2,h_bar1,h_bar2, cp_bar, dT = 1e-3;
            v_bar1=MolarVolume(T-dT,p,psi_w); //[m^3/mol_ha]
            h_bar1=MolarEnthalpy(T-dT,p,psi_w,v_bar1); //[kJ/kmol_ha]
            v_bar2=MolarVolume(T+dT,p,psi_w); //[m^3/mol_ha]
            h_bar2=MolarEnthalpy(T+dT,p,psi_w,v_bar2); //[kJ/kmol_ha]
            W=HumidityRatio(psi_w); //[kg_w/kg_da]
            cp_bar = (h_bar2-h_bar1)/(2*dT);
            return cp_bar/M_ha; //[kJ/kg_da]
        }
        else if (!strcmp(OutputName,"Tdp") || !strcmp(OutputName,"D"))
        {
            return DewpointTemperature(T,p,psi_w); //[K]
        }
        else if (!strcmp(OutputName,"Twb") || !strcmp(OutputName,"T_wb") || !strcmp(OutputName,"WetBulb") || !strcmp(OutputName,"B"))
        {
            return WetbulbTemperature(T,p,psi_w); //[K]
        }
        else if (!strcmp(OutputName,"Omega") || !strcmp(OutputName,"HumRat") || !strcmp(OutputName,"W"))
        {
            return HumidityRatio(psi_w);
        }
        else if (!strcmp(OutputName,"RH") || !strcmp(OutputName,"RelHum") || !strcmp(OutputName,"R"))
        {
            return RelativeHumidity(T,p,psi_w);
        }
        else if (!strcmp(OutputName,"mu") || !strcmp(OutputName,"Visc") || !strcmp(OutputName,"M"))
        {
            return Viscosity(T,p,psi_w);
        }
        else if (!strcmp(OutputName,"k") || !strcmp(OutputName,"Conductivity") || !strcmp(OutputName,"K"))
        {
            return Conductivity(T,p,psi_w);
        }
        else
        {
            return -1000;
        }
    }
    catch (std::exception &e)
    {
        set_err_string(e.what());
        return _HUGE;
    }
    catch (...)
    {
        return _HUGE;
    }
}

EXPORT_CODE double CONVENTION HAProps_Aux(const char* Name,double T, double p, double W, char *units)
{
    // This function provides some things that are not usually needed, but could be interesting for debug purposes.
    
    // Requires W since it is nice and fast and always defined.  Put a dummy value if you want something that doesn't use humidity
    
    // Takes temperature, pressure, and humidity ratio W as inputs;
    double psi_w,Tj,tau_Water,tau_Air,B_aa,C_aaa,B_ww,C_www,B_aw,C_aaw,C_aww,v_bar,delta, tau;
    
    Tj=132.6312;
    tau_Air=Tj/T;
    tau_Water=Water.reduce.T/T;
    
    try{
    if (!strcmp(Name,"Baa"))
    {
        B_aa=B_Air(T)*MM_Air()/1e3; //[m^3/kg] to [m^3/mol]
        strcpy(units,"m^3/mol");
        return B_aa;
    }
    else if (!strcmp(Name,"Caaa"))
    {
        C_aaa=C_Air(T)*MM_Air()*MM_Air()/1e6; //[m^6/kg^2] to [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_aaa;
    }
    else if (!strcmp(Name,"Bww"))
    {
        B_ww=B_Water(T)*MM_Water()/1e3; //[m^3/kg] to [m^3/mol]
        strcpy(units,"m^3/mol");
        return B_ww;
    }
    else if (!strcmp(Name,"Cwww"))
    {
        C_www=C_Water(T)*MM_Water()*MM_Water()/1e6; //[m^6/kg^2] to [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_www;
    }
    else if (!strcmp(Name,"dBaa"))
    {
        B_aa=dBdT_Air(T)*MM_Air()/1e3; //[m^3/kg] to [m^3/mol]
        strcpy(units,"m^3/mol");
        return B_aa;
    }
    else if (!strcmp(Name,"dCaaa"))
    {
        C_aaa=dCdT_Air(T)*MM_Air()*MM_Air()/1e6; //[m^6/kg^2] to [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_aaa;
    }
    else if (!strcmp(Name,"dBww"))
    {
        B_ww=dBdT_Water(T)*MM_Water()/1e3; //[m^3/kg] to [m^3/mol]
        strcpy(units,"m^3/mol");
        return B_ww;
    }
    else if (!strcmp(Name,"dCwww"))
    {
        C_www=dCdT_Water(T)*MM_Water()*MM_Water()/1e6; //[m^6/kg^2] to [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_www;
    }
    else if (!strcmp(Name,"Baw"))
    {
        B_aw=_B_aw(T)/1e3; //[dm^3/mol] to [m^3/mol]
        strcpy(units,"m^3/mol");
        return B_aw;
    }
    else if (!strcmp(Name,"Caww"))
    {
        C_aww=_C_aww(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_aww;
    }
    else if (!strcmp(Name,"Caaw"))
    {
        C_aaw=_C_aaw(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_aaw;
    }
    else if (!strcmp(Name,"dBaw"))
    {
        double dB_aw=_dB_aw_dT(T)/1e3; //[dm^3/mol] to [m^3/mol]
        strcpy(units,"m^3/mol");
        return dB_aw;
    }
    else if (!strcmp(Name,"dCaww"))
    {
        double dC_aww=_dC_aww_dT(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return dC_aww;
    }
    else if (!strcmp(Name,"dCaaw"))
    {
        double dC_aaw=_dC_aaw_dT(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return dC_aaw;
    }
    else if (!strcmp(Name,"beta_H"))
    {
        strcpy(units,"1/Pa");
        return HenryConstant(T);
    }
    else if (!strcmp(Name,"kT"))
    {
        strcpy(units,"1/Pa");
        if (T>273.16)
        {
            double rho = Props("D","T",T,"P",p,"Water");
            return DerivTerms((char *)"IsothermalCompressibility",T,rho,(char *)"Water")/1000; //[1/Pa]
        }
        else
            return IsothermCompress_Ice(T,p); //[1/Pa]
    }
    else if (!strcmp(Name,"p_ws"))
    {
        strcpy(units,"kPa");
        if (T>273.16)
            return Props("P","T",T,"Q",0,"Water");
        else
            return psub_Ice(T);
    }
    else if (!strcmp(Name,"vbar_ws"))
    {
        strcpy(units,"m^3/mol");
        if (T>273.16)
        {
            // It is liquid water
            return 1.0/Props("D","T",T,"Q",0,"Water")*MM_Water()/1000; //[m^3/mol]
        }
        else
        {
            // It is ice
            return dg_dp_Ice(T,p/1000)*MM_Water()/1000/1000; //[m^3/mol]
        }
    }
    else if (!strcmp(Name,"f"))
    {
        strcpy(units,"-");
        return f_factor(T,p);
    }
    // Get psi_w since everything else wants it
    psi_w=MoleFractionWater(T,p,GIVEN_HUMRAT,W);
    if (!strcmp(Name,"Bm"))
    {
        strcpy(units,"m^3/mol");
        return B_m(T,psi_w);
    }
    else if (!strcmp(Name,"Cm"))
    {
        strcpy(units,"m^6/mol^2");
        return C_m(T,psi_w);
    }
    else if (!strcmp(Name,"hvirial"))
    {
        v_bar=MolarVolume(T,p,psi_w);
        return 8.3145*T*((B_m(T,psi_w)-T*dB_m_dT(T,psi_w))/v_bar+(C_m(T,psi_w)-T/2.0*dC_m_dT(T,psi_w))/(v_bar*v_bar));
    }
    else if (!strcmp(Name,"ha"))
    {
        delta=1.1/322; tau=132/T;
        return 1+tau*DerivTerms((char *)"dphi0_dTau",tau,delta,(char *)"Water");
    }
    else if (!strcmp(Name,"hw"))
    {
        //~ return Props('D','T',T,'P',p,"Water")/322; tau=647/T;
        delta=1000/322; tau=647/T;
        //~ delta=rho_Water(T,p,TYPE_TP);tau=647/T;
        return 1+tau*DerivTerms((char *)"dphi0_dTau",tau,delta,(char *)"Water");
    }
    else if (!strcmp(Name,"hbaro_w"))
    {
        v_bar=MolarVolume(T,p,psi_w);
        return IdealGasMolarEnthalpy_Water(T,v_bar);
    }
    else if (!strcmp(Name,"hbaro_a"))
    {
        v_bar=MolarVolume(T,p,psi_w);
        return IdealGasMolarEnthalpy_Air(T,v_bar);
    }
    else
    {
        printf("Sorry I didn't understand your input [%s] to HAProps_Aux\n",Name);
        return -1;
    }
    }
    catch(std::exception &)
    {
        return _HUGE;
    }
    return _HUGE;
}
double cair_sat(double T)
{
    // Air saturation specific heat 
    // Based on a correlation from EES, good from 250K to 300K.  
    // No error bound checking is carried out
    // T: [K]
    // cair_s: [kJ/kg-K]
    return 2.14627073E+03-3.28917768E+01*T+1.89471075E-01*T*T-4.86290986E-04*T*T*T+4.69540143E-07*T*T*T*T;
}

double IceProps(const char* Name, double T, double p)
{
    if (!strcmp(Name,"s"))
    {
        return s_Ice(T,p*1000.0);
    }
    else if (!strcmp(Name,"rho"))
    {
        return rho_Ice(T,p*1000.0);
    }
    else if (!strcmp(Name,"h"))
    {
        return h_Ice(T,p*1000.0);
    }
    else
    {
        return 1e99;
    }
}


#ifndef CATCH_DISABLED
#include <math.h>
#include "Catch/catch.hpp"
TEST_CASE((char*)"Tests from ASHRAE RP-1485",(char*)"[RP1485]")
{
    SECTION((char*)"Table A.8.1")
    {
        double p1 = 101.325, T= 473.15;
                              // W      B    V      H      S
        std::string rows1[] ={ "0.00 45.07 1.341 202.52 0.5558",
                               "0.05 55.38 1.448 346.49 1.0299",
                               "0.10 61.85 1.556 490.43 1.4736",
                               "0.20 69.95 1.771 778.24 2.3336",
                               "0.30 75.00 1.986 1066.00 3.1751",
                               "0.40 78.51 2.201 1353.71 4.0059",
                               "0.50 81.12 2.416 1641.40 4.8295",
                               "0.60 83.14 2.630 1929.06 5.6479",
                               "0.70 84.76 2.845 2216.70 6.4623",
                               "0.80 86.09 3.060 2504.32 7.2736",
                               "0.90 87.20 3.274 2791.94 8.0824",
                               "1.00 88.15 3.489 3079.55 8.8890"};
        for (int i = 0; i < 12; i++)
        {
            std::vector<std::string> elements = strsplit(rows1[i],' ');
            double W = strtod(elements[0].c_str(),NULL);
            SECTION((char*)"B")
            {
                double B = strtod(elements[1].c_str(),NULL);
                double BCP = HAProps("B","T",T,"W",W,"P",p1) - 273.15;
                CAPTURE(B);
                CAPTURE(BCP);
                CHECK(fabs(BCP-B) < 0.01);
            }
            SECTION((char*)"V")
            {
                double V = strtod(elements[2].c_str(),NULL);
                double VCP = HAProps("V","T",T,"W",W,"P",p1);
                CAPTURE(V);
                CAPTURE(VCP);
                CHECK(fabs(VCP-V) < 0.01);
            }
            SECTION((char*)"H")
            {
                double H = strtod(elements[3].c_str(),NULL);
                double HCP = HAProps("H","T",T,"W",W,"P",p1);
                CAPTURE(H);
                CAPTURE(HCP);
                CHECK(fabs(HCP-H) < 0.01);
            }
            SECTION((char*)"S")
            {
                double S = strtod(elements[4].c_str(),NULL);
                double SCP = HAProps("S","T",T,"W",W,"P",p1);
                CAPTURE(S);
                CAPTURE(SCP);
                CHECK(fabs(SCP-S) < 0.01);
            }
        }
    }
}
#endif