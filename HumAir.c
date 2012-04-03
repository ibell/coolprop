#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include "math.h"
#include "time.h"
#include "stdio.h"
#include <string.h>

#include "PropMacros.h"
#include "CoolProp.h"
#include "Ice.h"
#include "SolverFunctions.h"
#include "HumAir.h"

static double epsilon=0.621945,R_bar=8.314472;
static int FlagUseVirialCorrelations=0,FlagUseIsothermCompressCorrelation=0,FlagUseIdealGasEnthalpyCorrelations=0;
double f_factor(double T, double p);

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
static double Secant_HAProps_T(char *OutputName, char *Input1Name, double Input1, char *Input2Name, double Input2, double TargetVal, double T_guess)
{
    // Use a secant solve in order to yield a target output value for HAProps by altering T
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f=999,T=300;
	int iter=1;

	while ((iter<=3 || fabs(f)>eps) && iter<100)
	{
		if (iter==1){x1=T_guess; T=x1;}
		if (iter==2){x2=T_guess+0.001; T=x2;}
		if (iter>2) {T=x2;}
			f=HAProps(OutputName,"T",T,Input1Name,Input1,Input2Name,Input2)-TargetVal;
		if (iter==1){y1=f;}
		if (iter>1)
		{
			y2=f;
			x3=x2-y2/(y2-y1)*(x2-x1);
			y1=y2; x1=x2; x2=x3;
		}
		iter=iter+1;
	}
	return T;
}

static double Secant_HAProps_W(char *OutputName, char *Input1Name, double Input1, char *Input2Name, double Input2, double TargetVal, double W_guess)
{
    // Use a secant solve in order to yield a target output value for HAProps by altering humidity ratio
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f=999,W=0.0001;
	int iter=1;

	while ((iter<=3 || fabs(f)>eps) && iter<100)
	{
		if (iter==1){x1=W_guess; W=x1;}
		if (iter==2){x2=W_guess+0.001; W=x2;}
		if (iter>2) {W=x2;}
			f=HAProps(OutputName,"W",W,Input1Name,Input1,Input2Name,Input2)-TargetVal;
		if (iter==1){y1=f;}
		if (iter>1)
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
    tau_Water=Tcrit_Water()/T;
    if (FlagUseVirialCorrelations==1)
    {
        B_aa=-0.000721183853646 +1.142682674467e-05*T -8.838228412173e-08*powI(T,2) 
        +4.104150642775e-10*powI(T,3) -1.192780880645e-12*powI(T,4) +2.134201312070e-15*powI(T,5) 
        -2.157430412913e-18*powI(T,6) +9.453830907795e-22*powI(T,7);
        B_ww=-10.8963128394 +2.439761625859e-01*T -2.353884845100e-03*powI(T,2) 
        +1.265864734412e-05*powI(T,3) -4.092175700300e-08*powI(T,4) +7.943925411344e-11*powI(T,5) 
        -8.567808759123e-14*powI(T,6) +3.958203548563e-17*powI(T,7);
    }
    else
    {
        B_aa=B_Air(tau_Air)*MM_Air()/1e3; //[m^3/kg] to [m^3/mol]
        B_ww=B_Water(tau_Water)*MM_Water()/1e3; //[m^3/kg] to [m^3/mol]
    }
    
    B_aw=_B_aw(T)/1e3; //[dm^3/mol] to [m^3/mol]
    return powI(1-psi_w,2)*B_aa+2*(1-psi_w)*psi_w*B_aw+psi_w*psi_w*B_ww;
}

static double dB_m_dT(double T, double psi_w)
{
    //dBm_dT has units of m^3/mol/K
    double Tj,tau_Air,tau_Water,dB_dT_aa,dB_dT_ww,dB_dT_aw;
    // NDG for fluid EOS for virial terms
    Tj=132.6312;
    tau_Air=Tj/T;
    tau_Water=Tcrit_Water()/T;
    if (FlagUseVirialCorrelations)
    {
        dB_dT_aa=1.65159324353e-05 -3.026130954749e-07*T +2.558323847166e-09*powI(T,2) -1.250695660784e-11*powI(T,3) +3.759401946106e-14*powI(T,4) -6.889086380822e-17*powI(T,5) +7.089457032972e-20*powI(T,6) -3.149942145971e-23*powI(T,7);
        dB_dT_ww=0.65615868848 -1.487953162679e-02*T +1.450134660689e-04*powI(T,2) -7.863187630094e-07*powI(T,3) +2.559556607010e-09*powI(T,4) -4.997942221914e-12*powI(T,5) +5.417678681513e-15*powI(T,6) -2.513856275241e-18*powI(T,7);
    }
    else
    {
        dB_dT_aa=dBdT_Air(tau_Air)*MM_Air()/1e3; //[m^3/kg] to [m^3/mol]
        dB_dT_ww=dBdT_Water(tau_Water)*MM_Water()/1e3; //[m^3/kg] to [m^3/mol]
    }
        dB_dT_aw=_dB_aw_dT(T)/1e3; //[dm^3/mol] to [m^3/mol]
    return powI(1-psi_w,2)*dB_dT_aa+2*(1-psi_w)*psi_w*dB_dT_aw+psi_w*psi_w*dB_dT_ww;
}

static double C_m(double T, double psi_w)
{
    // Cm has units of m^6/mol^2
    double Tj,tau_Air,tau_Water,C_aaa,C_www,C_aww,C_aaw;
    // NDG for fluid EOS for virial terms
    Tj=132.6312;
    tau_Air=Tj/T;
    tau_Water=Tcrit_Water()/T;
    if (FlagUseVirialCorrelations)
    {
        C_aaa=1.29192158975e-08 -1.776054020409e-10*T +1.359641176409e-12*powI(T,2) 
        -6.234878717893e-15*powI(T,3) +1.791668730770e-17*powI(T,4) -3.175283581294e-20*powI(T,5) 
        +3.184306136120e-23*powI(T,6) -1.386043640106e-26*powI(T,7);
        C_www=-0.580595811134 +1.365952762696e-02*T -1.375986293288e-04*powI(T,2) 
        +7.687692259692e-07*powI(T,3) -2.571440816920e-09*powI(T,4) +5.147432221082e-12*powI(T,5) 
        -5.708156494894e-15*powI(T,6) +2.704605721778e-18*powI(T,7);
    }
    else
    {
        C_aaa=C_Air(tau_Air)*MM_Air()*MM_Air()/1e6; //[m^6/kg^2] to [m^6/mol^2]
        C_www=C_Water(tau_Water)*MM_Water()*MM_Water()/1e6; //[m^6/kg^2] to [m^6/mol^2]
    }
    C_aaw=_C_aaw(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
    C_aww=_C_aww(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
    return powI(1-psi_w,3)*C_aaa+3*powI(1-psi_w,2)*psi_w*C_aaw+3*(1-psi_w)*psi_w*psi_w*C_aww+powI(psi_w,3)*C_www;
}

static double dC_m_dT(double T, double psi_w)
{
    // dCm_dT has units of m^6/mol^2/K
    
    double Tj,tau_Air,tau_Water,dC_dT_aaa,dC_dT_www,dC_dT_aww,dC_dT_aaw;
    // NDG for fluid EOS for virial terms
    Tj=132.6312;
    tau_Air=Tj/T;
    tau_Water=Tcrit_Water()/T;
    if (FlagUseVirialCorrelations)
    {
        dC_dT_aaa=-2.46582342273e-10 +4.425401935447e-12*T -3.669987371644e-14*powI(T,2) +1.765891183964e-16*powI(T,3) -5.240097805744e-19*powI(T,4) +9.502177003614e-22*powI(T,5) -9.694252610339e-25*powI(T,6) +4.276261986741e-28*powI(T,7);
        dC_dT_www=0.0984601196142 -2.356713397262e-03*T +2.409113323685e-05*powI(T,2) -1.363083778715e-07*powI(T,3) +4.609623799524e-10*powI(T,4) -9.316416405390e-13*powI(T,5) +1.041909136255e-15*powI(T,6) -4.973918480607e-19*powI(T,7);
    }
    else
    {
        dC_dT_aaa=dCdT_Air(tau_Air)*MM_Air()*MM_Air()/1e6; //[m^6/kg^2] to [m^6/mol^2]
        dC_dT_www=dCdT_Water(tau_Water)*MM_Water()*MM_Water()/1e6; //[m^6/kg^2] to [m^6/mol^2]
    }
    dC_dT_aaw=_dC_aaw_dT(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
    dC_dT_aww=_dC_aww_dT(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
    return powI(1-psi_w,3)*dC_dT_aaa+3*powI(1-psi_w,2)*psi_w*dC_dT_aaw+3*(1-psi_w)*psi_w*psi_w*dC_dT_aww+powI(psi_w,3)*dC_dT_www;
}

static double HenryConstant(double T)
{
    // Result has units of 1/Pa
    double p_ws,beta_N2,beta_O2,beta_Ar,beta_a,tau,Tr,Tc=647.096;
    Tr=T/Tc; 
    tau=1-Tr;
    p_ws=Props('P','T',T,'Q',1.0,"Water"); //[kPa]
    beta_N2=p_ws*exp(-9.67578/Tr+4.72162*pow(tau,0.355)/Tr+11.70585*pow(Tr,-0.41)*exp(tau));
    beta_O2=p_ws*exp(-9.44833/Tr+4.43822*pow(tau,0.355)/Tr+11.42005*pow(Tr,-0.41)*exp(tau));
    beta_Ar=p_ws*exp(-8.40954/Tr+4.29587*pow(tau,0.355)/Tr+10.52779*pow(Tr,-0.41)*exp(tau));
    beta_a=1/(0.7812/beta_N2+0.2095/beta_O2+0.0093/beta_Ar);
    return 1/(1.01325*beta_a)/1000.0;
}

double f_factor(double T, double p)
{
    double f,Rbar=8.314371,eps=1e-8,Tj;
    double x1,x2,x3,y1,y2,change;
    int iter=1;
    double p_ws,tau_Air,tau_Water,B_aa,B_aw,B_ww,C_aaa,C_aaw,C_aww,C_www,
        line1,line2,line3,line4,line5,line6,line7,line8,k_T,beta_H,LHS,RHS,psi_ws,
        vbar_ws;
    
    //// Use correlation for f at atmospheric pressure
    //if (p<1.001*101.325 && p>0.999*101.325 && T>213.15 && T<373.15)
    //{
    //    return 7.1285352695 -1.580759239351e-01*T +1.758777627765e-03*powI(T,2) -1.091297459853e-05*powI(T,3) +4.075109681028e-08*powI(T,4) -9.156215217527e-11*powI(T,5) +1.146098035510e-13*powI(T,6) -6.163466181798e-17*powI(T,7);
    //}

    // Get total pressure in Pa from kPa
    p*=1000;
    
    // Saturation pressure [Pa]
    if (T>=273.15)
    {
        // It is liquid water
        p_ws=Props('P','T',T,'Q',0,"Water")*1000;
		if (FlagUseIsothermCompressCorrelation)
		{
			k_T = 1.6261876614E-22*powI(T,6) - 3.3016385196E-19*powI(T,5) + 2.7978984577E-16*powI(T,4)
				- 1.2672392901E-13*powI(T,3) + 3.2382864853E-11*powI(T,2) - 4.4318979503E-09*T + 2.5455947289E-07;
		}
        else
		{
			k_T=IsothermCompress_Water(T,p/1000); //[1/Pa]
		}
        beta_H=HenryConstant(T); //[1/Pa]
        vbar_ws=1.0/Props('D','T',T,'Q',0,"Water")*MM_Water()/1000; //[m^3/mol]
    }
    else
    {
        // It is ice
        p_ws=psub_Ice(T)*1000;
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
    tau_Water=Tcrit_Water()/T;
	if (FlagUseVirialCorrelations)
	{
		B_aa=-0.000721183853646 +1.142682674467e-05*T -8.838228412173e-08*powI(T,2) 
        +4.104150642775e-10*powI(T,3) -1.192780880645e-12*powI(T,4) +2.134201312070e-15*powI(T,5) 
        -2.157430412913e-18*powI(T,6) +9.453830907795e-22*powI(T,7);
        B_ww=-10.8963128394 +2.439761625859e-01*T -2.353884845100e-03*powI(T,2) 
        +1.265864734412e-05*powI(T,3) -4.092175700300e-08*powI(T,4) +7.943925411344e-11*powI(T,5) 
        -8.567808759123e-14*powI(T,6) +3.958203548563e-17*powI(T,7);
		C_aaa=1.29192158975e-08 -1.776054020409e-10*T +1.359641176409e-12*powI(T,2) 
        -6.234878717893e-15*powI(T,3) +1.791668730770e-17*powI(T,4) -3.175283581294e-20*powI(T,5) 
        +3.184306136120e-23*powI(T,6) -1.386043640106e-26*powI(T,7);
        C_www=-0.580595811134 +1.365952762696e-02*T -1.375986293288e-04*powI(T,2) 
        +7.687692259692e-07*powI(T,3) -2.571440816920e-09*powI(T,4) +5.147432221082e-12*powI(T,5) 
        -5.708156494894e-15*powI(T,6) +2.704605721778e-18*powI(T,7);
	}
	else
	{
		B_aa=B_Air(tau_Air)*MM_Air()/1e3; //[m^3/kg] to [m^3/mol]
		C_aaa=C_Air(tau_Air)*MM_Air()*MM_Air()/1e6; //[m^6/kg^2] to [m^6/mol^2]
		B_ww=B_Water(tau_Water)*MM_Water()/1e3; //[m^3/kg] to [m^3/mol]
		C_www=C_Water(tau_Water)*MM_Water()*MM_Water()/1e6; //[m^6/kg^2] to [m^6/mol^2]
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
            line2=powI(1-psi_ws,2)*p/(Rbar*T)*B_aa-2*powI(1-psi_ws,2)*p/(Rbar*T)*B_aw-(p-p_ws-powI(1-psi_ws,2)*p)/(Rbar*T)*B_ww;
            line3=powI(1-psi_ws,3)*p*p/powI(Rbar*T,2)*C_aaa+(3*powI(1-psi_ws,2)*(1-2*(1-psi_ws))*p*p)/(2*powI(Rbar*T,2))*C_aaw;
            line4=-3*powI(1-psi_ws,2)*psi_ws*p*p/powI(Rbar*T,2)*C_aww-((3-2*psi_ws)*psi_ws*psi_ws*p*p-p_ws*p_ws)/(2*powI(Rbar*T,2))*C_www;
            line5=-(powI(1-psi_ws,2)*(-2+3*psi_ws)*psi_ws*p*p)/powI(Rbar*T,2)*B_aa*B_ww;
            line6=-(2*powI(1-psi_ws,3)*(-1+3*psi_ws)*p*p)/powI(Rbar*T,2)*B_aa*B_aw;
            line7=(6*powI(1-psi_ws,2)*psi_ws*psi_ws*p*p)/powI(Rbar*T,2)*B_ww*B_aw-(3*powI(1-psi_ws,4)*p*p)/(2*powI(Rbar*T,2))*B_aa*B_aa;
            line8=-(2*powI(1-psi_ws,2)*psi_ws*(-2+3*psi_ws)*p*p)/powI(Rbar*T,2)*B_aw*B_aw-(p_ws*p_ws-(4-3*psi_ws)*powI(psi_ws,3)*p*p)/(2*powI(Rbar*T,2))*B_ww*B_ww;
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
    return f;
}
void HAHelp()
{
    double Tj=132.6312,M_Air=28.9586,M_Water=18.015;
	printf("Sorry, Need to update!");
}
int returnHumAirCode(char * Code)
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
    mu_a=Props('V','T',T,'P',p,"Air");
    // Viscosity of pure saturated water at dry-bulb temperature
    mu_w=Props('V','P',p,'Q',1,"Water");
    Phi_av=sqrt(2)/4.0*pow(1+Ma/Mw,-0.5)*powI(1+sqrt(mu_a/mu_w)*pow(Mw/Ma,0.25),2); //[-]
    Phi_va=sqrt(2)/4.0*pow(1+Mw/Ma,-0.5)*powI(1+sqrt(mu_w/mu_a)*pow(Ma/Mw,0.25),2); //[-]
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
	UseSaturationLUT(1); // Use the lookup table
    // Viscosity of dry air at dry-bulb temp and total pressure
    k_a=Props('L','T',T,'P',p,"Air");
    mu_a=Props('V','T',T,'P',p,"Air");
    // Viscosity of pure saturated water at dry-bulb temperature
    k_w=Props('L','P',p,'Q',1,"Water");
    mu_w=Props('V','P',p,'Q',1,"Water");
    Phi_av=sqrt(2)/4.0*pow(1+Ma/Mw,-0.5)*powI(1+sqrt(mu_a/mu_w)*pow(Mw/Ma,0.25),2); //[-]
    Phi_va=sqrt(2)/4.0*pow(1+Mw/Ma,-0.5)*powI(1+sqrt(mu_w/mu_a)*pow(Ma/Mw,0.25),2); //[-]
    return (1-psi_w)*k_a/((1-psi_w)+psi_w*Phi_av)+psi_w*k_w/(psi_w+(1-psi_w)*Phi_va);
}
double MolarVolume(double T, double p, double psi_w)
{
    // Output in m^3/mol
    int iter;
    double v_bar0, v_bar, R_bar=8.314472,x1,x2,x3,y1,y2,resid,change,eps,Bm,Cm;
    
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
			change=fabs(y2/(y2-y1)*(x2-x1));
			y1=y2; x1=x2; x2=x3;
		}
		iter=iter+1;
	}
    return v_bar;
}
double IdealGasMolarEnthalpy_Water(double T, double v_bar)
{
	double hbar_w_0,tau,rhobar,delta,hbar_w;
	// Ideal-Gas contribution to enthalpy of water
    hbar_w_0=-0.01102303806;//[kJ/kmol]
    tau=Tcrit_Water()/T; 
    rhobar=322/MM_Water()*1000;
	delta=1/(v_bar*rhobar);
	hbar_w=hbar_w_0+R_bar*T*(1+tau*dphi0_dTau_Water(tau,delta));
	return hbar_w;
}
double IdealGasMolarEnthalpy_Air(double T, double v_bar)
{
	double hbar_a_0,tau,rhobar,delta,hbar_a,R_bar_Lemmon;
	// Ideal-Gas contribution to enthalpy of air
    hbar_a_0=-7914.149298; //[kJ/kmol]
    //Tj and rhoj are given by 132.6312 and 302.5507652 respectively
    tau=132.6312/T;
    rhobar=302.5507652/MM_Air()*1000;
    delta=1/(v_bar*rhobar);
    R_bar_Lemmon=8.314510; //[kJ/kmol/K]
    hbar_a=hbar_a_0+R_bar_Lemmon*T*(1+tau*dphi0_dTau_Air(tau,delta)); //[kJ/kmol]
	return hbar_a;
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
    hbar_0=0;//2.924425468; //[kJ/kmol]
    
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
double DewpointTemperature(double T, double p, double psi_w)
{
    int iter;
    double p_w,epsilon=0.621945,eps,resid,Tdp,x1,x2,x3,y1,y2,change;
    double p_ws_dp,f_dp;
    
    // Make sure it isn't dry air, return an impossible temperature otherwise
    if ((1-psi_w)<1e-16)
    {
        return -1;
    }
    // ------------------------------------------
    // Iteratively find the dewpoint temperature
    // ------------------------------------------
    
    p_w=psi_w*p;
    
    iter=1; eps=1e-8; resid=999;
    while ((iter<=3 || fabs(resid)>eps) && iter<100)
	{
		if (iter==1){x1=T; Tdp=x1;}
		if (iter==2){x2=T+0.000001; Tdp=x2;}
		if (iter>2) {Tdp=x2;}
        
            if (T>=273.15)
            {
                // Saturation pressure at dewpoint [kPa]
                UseSaturationLUT(1);
                p_ws_dp=Props('P','T',Tdp,'Q',0,"Water");
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
			change=fabs(y2/(y2-y1)*(x2-x1));
			y1=y2; x1=x2; x2=x3;
		}
		iter=iter+1;
	}
    return Tdp;
}

double WetbulbTemperature(double T, double p, double psi_w)
{
    int iter;
    double epsilon=0.621945,eps,resid,Twb,x1,x2,x3,y1,y2,change;
    double p_ws_wb,f_wb,W,W_s_wb,p_s_wb,h_w,psi_wb,M_ha_wb,M_ha;
    double v_bar_w,v_bar_wb;
    
    // ------------------------------------------
    // Iteratively find the wetbulb temperature
    // ------------------------------------------
    
    W=epsilon*psi_w/(1-psi_w);
    
    iter=1; eps=1e-8; resid=999;
    while ((iter<=3 || fabs(resid)>eps) && iter<100)
	{
		if (iter==1){x1=T; Twb=x1;}
		if (iter==2){x2=T+0.000001; Twb=x2;}
		if (iter>2) {Twb=x2;}
        
            // Enhancement Factor at wetbulb temperature [-]
            f_wb=f_factor(Twb,p);
            if (T>273.15)
            {
                // Saturation pressure at wetbulb temperature [kPa]
                p_ws_wb=Props('P','T',Twb,'Q',0,"Water");
            }
            else
            {
                // Sublimation pressure at wetbulb temperature [kPa]
                p_ws_wb=psub_Ice(Twb);
            }
                
            // Vapor pressure
            p_s_wb=f_wb*p_ws_wb;
            // wetbulb humidity ratio
            W_s_wb=epsilon*p_s_wb/(p-p_s_wb);
            // wetbulb water mole fraction
            psi_wb=W_s_wb/(epsilon+W_s_wb);
            if (T>=273.15)
            {
                // Enthalpy of water [kJ/kg_water]
                h_w=Props('H','T',Twb,'P',p,"Water");
            }
            else
            {
                // Enthalpy of ice [kJ/kg_water]
                h_w=h_Ice(T,p)/1000;
            }
            // Mole masses of wetbulb and humid air
            M_ha=MM_Water()*psi_w+(1-psi_w)*28.966;
            M_ha_wb=MM_Water()*psi_wb+(1-psi_wb)*28.966;
            v_bar_w=MolarVolume(T,p,psi_w);
            v_bar_wb=MolarVolume(T,p,psi_wb);
            // Error between target and actual pressure [kJ/kg_da]
            resid=MolarEnthalpy(T,p,psi_w,v_bar_w)/M_ha*(1+W)-(MolarEnthalpy(Twb,p,psi_wb,v_bar_wb)/M_ha_wb*(1+W_s_wb)+(W-W_s_wb)*h_w);
        
        if (iter==1){y1=resid;}
		if (iter>1)
		{
			y2=resid;
			x3=x2-y2/(y2-y1)*(x2-x1);
			change=fabs(y2/(y2-y1)*(x2-x1));
			y1=y2; x1=x2; x2=x3;
		}
		iter=iter+1;
	}
    return Twb;
}
static int Name2Type(char *Name)
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
    else if (!strcmp(Name,"mu") || !strcmp(Name,"Visc") || !strcmp(Name,"M"))
        return GIVEN_VISC;
    else if (!strcmp(Name,"k") || !strcmp(Name,"Conductivity") || !strcmp(Name,"K"))
        return GIVEN_COND;
    else
        printf("Sorry, your input [%s] was not understood to Name2Type in HumAir.c. Acceptable values are T,P,R,W,D,B,H,M,K and aliases thereof\n");
        return -1;
}
int TypeMatch(int TypeCode,char *Input1Name, char *Input2Name, char *Input3Name)
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
    UseSaturationLUT(1); // Enable the use of lookup tables for saturation properties for speed
    
    if (HumInput==GIVEN_HUMRAT) //(2)
    {
        W=InVal;
        return W/(epsilon+W);
    }
    else if (HumInput==GIVEN_RH)
    {
        if (T>=273.15)
        {
            // Saturation pressure [kPa]
            p_ws=Props('P','T',T,'Q',0,"Water");
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
        p_ws_dp=Props('P','T',Tdp,'Q',0,"Water");
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
double HumidityRatio(double psi_w)
{
    return psi_w*epsilon/(1-psi_w);
}
double RelativeHumidity(double T, double p, double psi_w)
{
    double p_ws,f,p_s,W;
    if (T>=273.15)
    {
        // Saturation pressure [kPa]
        p_ws=Props('P','T',T,'Q',0,"Water");
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
double HAProps(char *OutputName, char *Input1Name, double Input1, char *Input2Name, double Input2, char *Input3Name, double Input3)
{
    int In1Type, In2Type, In3Type,iT,iW,iTdp,iRH,ip,Type1,Type2;
    double vals[3],p,T,RH,W,Tdp,psi_w,M_ha,v_bar,h_bar,MainInputValue,SecondaryInputValue,T_guess;
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
            // Find the value for W
            W_guess=0.001;
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
        T_guess=278.15;
        
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

		//if (!strcmp(SecondaryInputName,"H"))
		//	h_star=log(Value2+33);
		//else
		//	h_star=log(Value1+33);
		//T_guess= -7.4251055543E-02*powI(h_star,6) + 1.0661647745E-01*powI(h_star,5) + 8.5881364720E+00*powI(h_star,4) - 7.2409797021E+01*powI(h_star,3) + 2.3508812707E+02*powI(h_star,2) - 2.8041007078E+02*h_star + 3.5922309997E+02;

        // Use the secant solver to find T
        T=Secant_HAProps_T(SecondaryInputName,"P",p,MainInputName,MainInputValue,SecondaryInputValue,T_guess);
        
        // If you want the temperature, return it
        if (Name2Type(OutputName)==GIVEN_T)
            return T;
        else
        {
            // Otherwise, find psi_w for further calculations in the following section
            W=HAProps("W","T",T,"P",p,MainInputName,MainInputValue);
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
        v_bar=MolarVolume(T,p,psi_w); //[m^3/mol_ha]
        h_bar=MolarEnthalpy(T,p,psi_w,v_bar); //[kJ/kmol_ha]
        W=HumidityRatio(psi_w); //[kg_w/kg_da]
        return h_bar*(1+W)/M_ha; //[kJ/kg_da]
    }
    else if (!strcmp(OutputName,"Hha"))
    {
        v_bar=MolarVolume(T,p,psi_w); //[m^3/mol_ha]
        h_bar=MolarEnthalpy(T,p,psi_w,v_bar); //[kJ/kmol_ha]
        return h_bar/M_ha; //[kJ/kg_ha]
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

double HAProps_Aux(char* Name,double T, double p, double W, char *units)
{
    // This function provides some things that are not usually needed, but could be interesting for debug purposes.
    
    // Requires W since it is nice and fast and always defined.  Put a dummy value if you want something that doesn't use humidity
    
    // Takes temperature, pressure, and humidity ratio W as inputs;
    double psi_w,Tj,tau_Water,tau_Air,B_aa,C_aaa,B_ww,C_www,B_aw,C_aaw,C_aww,p_ws,v_bar,delta, tau;
    
    Tj=132.6312;
    tau_Air=Tj/T;
    tau_Water=Tcrit_Water()/T;
    
    if (!strcmp(Name,"Baa"))
    {
        B_aa=B_Air(tau_Air)*MM_Air()/1e3; //[m^3/kg] to [m^3/mol]
        strcpy(units,"m^3/mol");
        return B_aa;
    }
    else if (!strcmp(Name,"Caaa"))
    {
        C_aaa=C_Air(tau_Air)*MM_Air()*MM_Air()/1e6; //[m^6/kg^2] to [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_aaa;
    }
    else if (!strcmp(Name,"Bww"))
    {
        B_ww=B_Water(tau_Water)*MM_Water()/1e3; //[m^3/kg] to [m^3/mol]
        strcpy(units,"m^3/mol");
        return B_ww;
    }
    else if (!strcmp(Name,"Cwww"))
    {
        C_www=C_Water(tau_Water)*MM_Water()*MM_Water()/1e6; //[m^6/kg^2] to [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_www;
    }
    else if (!strcmp(Name,"dBaa"))
    {
        B_aa=dBdT_Air(tau_Air)*MM_Air()/1e3; //[m^3/kg] to [m^3/mol]
        strcpy(units,"m^3/mol");
        return B_aa;
    }
    else if (!strcmp(Name,"dCaaa"))
    {
        C_aaa=dCdT_Air(tau_Air)*MM_Air()*MM_Air()/1e6; //[m^6/kg^2] to [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_aaa;
    }
    else if (!strcmp(Name,"dBww"))
    {
        B_ww=dBdT_Water(tau_Water)*MM_Water()/1e3; //[m^3/kg] to [m^3/mol]
        strcpy(units,"m^3/mol");
        return B_ww;
    }
    else if (!strcmp(Name,"dCwww"))
    {
        C_www=dCdT_Water(tau_Water)*MM_Water()*MM_Water()/1e6; //[m^6/kg^2] to [m^6/mol^2]
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
    else if (!strcmp(Name,"beta_H"))
    {
        strcpy(units,"1/Pa");
        return HenryConstant(T);
    }
    else if (!strcmp(Name,"kT"))
    {
        strcpy(units,"1/Pa");
        if (T>273.15)
            return IsothermCompress_Water(T,p); //[1/Pa]
        else
            return IsothermCompress_Ice(T,p); //[1/Pa]
    }
    else if (!strcmp(Name,"p_ws"))
    {
        strcpy(units,"kPa");
        if (T>=273.15)
            return Props('P','T',T,'Q',0,"Water");
        else
            return psub_Ice(T);
    }
    else if (!strcmp(Name,"vbar_ws"))
    {
        strcpy(units,"m^3/mol");
        if (T>=273.15)
        {
            // It is liquid water
            return 1.0/Props('D','T',T,'Q',0,"Water")*MM_Water()/1000; //[m^3/mol]
        }
        else
        {
            // It is ice
            p_ws=psub_Ice(T)*1000;
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
        return 1+tau*dphi0_dTau_Air(tau,delta);
    }
    else if (!strcmp(Name,"hw"))
    {
        //~ return Props('D','T',T,'P',p,"Water")/322; tau=647/T;
        delta=1000/322; tau=647/T;
        //~ delta=rho_Water(T,p,TYPE_TP);tau=647/T;
        return 1+tau*dphi0_dTau_Water(tau,delta);
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
    }
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
