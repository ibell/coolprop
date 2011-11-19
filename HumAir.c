#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include "math.h"
#include "stdio.h"
#include <string.h>
#include "HumAir.h"
#include "PropMacros.h"
#include "Water.h"
#include "Air.h"
#include "CoolProp.h"
#include "Ice.h"

static double powI(double x, int y);
double f_factor(double T, double p);
static int HumAir_func(double T, double p, int HumInput, double xSI, /* in --- out */ double *Tdp_out, double *W_out, double *h_out, double *RH_out, double *v_out);

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
    B_aa=B_Air(tau_Air)*MM_Air()/1e3; //[m^3/kg] to [m^3/mol]
    B_ww=B_Water(tau_Water)*MM_Water()/1e3; //[m^3/kg] to [m^3/mol]
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
    dB_dT_aa=dBdT_Air(tau_Air)*MM_Air()/1e3; //[m^3/kg] to [m^3/mol]
    dB_dT_ww=dBdT_Water(tau_Water)*MM_Water()/1e3; //[m^3/kg] to [m^3/mol]
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
    C_aaa=C_Air(tau_Air)*MM_Air()*MM_Air()/1e6; //[m^6/kg^2] to [m^6/mol^2]
    C_www=C_Water(tau_Water)*MM_Water()*MM_Water()/1e6; //[m^6/kg^2] to [m^6/mol^2]
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
    dC_dT_aaa=dCdT_Air(tau_Air)*MM_Air()*MM_Air()/1e6; //[m^6/kg^2] to [m^6/mol^2]
    dC_dT_www=dCdT_Water(tau_Water)*MM_Water()*MM_Water()/1e6; //[m^6/kg^2] to [m^6/mol^2]
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
    
    // Get total pressure in Pa from kPa
    p*=1000;
    
    // Saturation pressure [Pa]
    if (T>=273.15)
    {
        // It is liquid water
        p_ws=Props('P','T',T,'Q',0,"Water")*1000;
        k_T=IsothermCompress_Water(T,p/1000); //[1/Pa]
        beta_H=HenryConstant(T); //[1/Pa]
        vbar_ws=1.0/Props('D','T',T,'P',p_ws/1000,"Water")*MM_Water()/1000; //[m^3/mol]
    }
    else
    {
        // It is ice
        p_ws=psub_Ice(T)*1000;
        k_T=IsothermCompress_Ice(T,p/1000); //[1/Pa]
        beta_H=0;
        vbar_ws=dg_dp_Ice(T,p/1000)*MM_Water()/1000/1000; //[m^3/mol]
    }
    
    // NDG for fluid EOS for virial terms
    Tj=132.6312;
    tau_Air=Tj/T;
    tau_Water=Tcrit_Water()/T;
    B_aa=B_Air(tau_Air)*MM_Air()/1e3; //[m^3/kg] to [m^3/mol]
    C_aaa=C_Air(tau_Air)*MM_Air()*MM_Air()/1e6; //[m^6/kg^2] to [m^6/mol^2]
    B_ww=B_Water(tau_Water)*MM_Water()/1e3; //[m^3/kg] to [m^3/mol]
    C_www=C_Water(tau_Water)*MM_Water()*MM_Water()/1e6; //[m^6/kg^2] to [m^6/mol^2]
    B_aw=_B_aw(T)/1e3; //[dm^3/mol] to [m^3/mol]
    C_aaw=_C_aaw(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
    C_aww=_C_aww(T)/1e6; //[dm^6/mol] to [m^6/mol^2]
    
    //~ printf("B_aa %g C_aaa %g B_ww %g C_www %g B_aw %g C_aaw %g C_aww %g\n",B_aa*1e6, C_aaa*1e12, B_ww*1e6, C_www*1e12, B_aw*1e6,C_aaw*1e12,C_aww*1e12);
    
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
    //~ double T,p,Tdp_out,W_out,h_out,RH_out,v_out;
	printf(
		"First two inputs are always T [K], p [kPa absolute]\n"
		"Next input tells what the next input is. Codes available:\n"
		"%d : Dewpoint temperature [K]\n"
		"%d : Humidity ratio [kg water / kg dry air]\n"
		"%d : Given wet bulb temperature [K] (not currently coded)\n"
		"%d : Relative humidity in range [0-->1]\n"
		"%d : Enthalpy in [kJ/kg-K]\n"
		,GIVEN_TDP,GIVEN_HUMRAT,GIVEN_TWB,GIVEN_RH,GIVEN_ENTHALPY);
    
    //~ T=80+273.15;
    //~ p=101.325;
    //~ HumAir3(T,p,4,0.0001, /* in --- out */ &Tdp_out, &W_out, &h_out, &RH_out, &v_out);
    //~ printf("%g C: DP %0.8g W %0.8f h %0.8f RH %0.8f v %0.8f\n",T-273.15,Tdp_out,W_out,h_out,RH_out,v_out);
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

double MolarVolume(double T, double p, double psi_w)
{
    // Output in m^3/mol
    int iter;
    double v_bar0, v_bar, R_bar=8.314472,x1,x2,x3,y1,y2,resid,change,eps;
    
    // -----------------------------
    // Iteratively find molar volume
    // -----------------------------
    
    // Start by assuming it is an ideal gas to get initial guess
    v_bar0=R_bar*T/p/1000;
    
    iter=1; eps=1e-8; resid=999;
    while ((iter<=3 || fabs(resid)>eps) && iter<100)
	{
		if (iter==1){x1=v_bar0; v_bar=x1;}
		if (iter==2){x2=v_bar0+0.000001; v_bar=x2;}
		if (iter>2) {v_bar=x2;}
        
            // factor of 1000 is to deal with kmol and mol conversion
            // want v_bar in m^3/mol and R_bar in kJ/mol-K
            resid=p-(R_bar/1000)*T/v_bar*(1+B_m(T,psi_w)/v_bar+C_m(T,psi_w)/(v_bar*v_bar)); 
        
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
double MolarEnthalpy(double T, double p, double psi_w)
{
    // In units of kJ/kmol
    
    double hbar_0,hbar_a_0,tau,R_bar_Lemmon,hbar_a,hbar_w_0,hbar_w,v_bar,delta,hbar,R_bar=8.314472;
    // ----------------------------------------
    //      Enthalpy
    // ----------------------------------------
    
    //Get the molar volume
    v_bar=MolarVolume(T,p,psi_w);
    
    // Constant for enthalpy
    // Not clear why getting rid of this term yields the correct values in the table, but enthalpies are equal to an additive constant, so not a big deal
    hbar_0=0;//2.924425468; //[kJ/kmol]
    
    // Ideal-Gas contribution to enthalpy of air
    hbar_a_0=-7914.149298; //[kJ/kmol]
    //Tj and rhoj are given by 132.6312 and 302.5507652 respectively
    tau=132.6312/T; delta=Props('D','T',T,'P',p,"Air")/302.5507652;
    R_bar_Lemmon=8.314510; //[kJ/kmol/K]
    hbar_a=hbar_a_0+R_bar_Lemmon*T*(1+tau*dphi0_dTau_Air(tau,delta)); //[kJ/kmol]
    
    // Ideal-Gas contribution to enthalpy of water
    hbar_w_0=-0.01102303806;//[kJ/kmol]
    tau=Tcrit_Water()/T; delta=Props('D','T',T,'P',p,"Water")/322.0;
    hbar_w=hbar_w_0+R_bar*T*(1+tau*dphi0_dTau_Water(tau,delta));
    
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
                
            // 
            p_s_wb=f_wb*p_ws_wb;
            // wetbulb humidity ratio
            W_s_wb=epsilon*p_s_wb/(p-p_s_wb);
            // wetbulb water mole fraction
            psi_wb=W_s_wb/(epsilon+W_s_wb);
            if (T>=273.15)
            {
                // Enthalpy of water [kJ/kg]
                h_w=Props('H','T',Twb,'P',p,"Water");
            }
            else
            {
                // Enthalpy of ice
                h_w=h_Ice(T,p);
            }
            // Mole masses of wetbulb and humid air
            M_ha=MM_Water()*psi_w+(1-psi_w)*28.966;
            M_ha_wb=MM_Water()*psi_wb+(1-psi_wb)*28.966;
            // Error between target and actual pressure [kJ/kg_da]
            resid=MolarEnthalpy(T,p,psi_w)/M_ha*(1+W)-(MolarEnthalpy(Twb,p,psi_wb)/M_ha_wb*(1+W_s_wb)+(W-W_s_wb)*h_w);
        
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

static int HumAir_func(double T, double p, int HumInput, double xSI, /* in --- out */ double *Tdp_out, double *W_out, double *h_out, double *RH_out, double *v_out)
{
    /*
    Based (mostly) on the analysis of Herrmann, Kretzschmar, and Gatley in
    
    S. Herrmann, H.-J. Kretzschmar and D.P. Gatley, "ASHRAE RP-1485: Thermodynamic Properties of Real Moist Air, Dry Air, Steam, Water and Ice"
    
    There is a typo in the December 11, 2009 version, the molar rho_star value for C_aww should be 1 rather than 1000 so far as I (Ian Bell) can tell
    */
    double R_bar=8.314472,epsilon=0.621945,p_ws,p_s,f,RH,Tdp,W,psi_w,p_ws_dp,f_dp,p_w_dp;
    double M_ha,v_bar,v_da,hbar,Twb;
    
    UseSaturationLUT(1); // Enable the use of lookup tables for saturation properties for speed
    if (T>=273.15)
    {
        // Saturation pressure [kPa]
        p_ws=Props('P','T',T,'Q',0,"Water");
    }
    else
    {
        // Saturation pressure [kPa]
        p_ws=psub_Ice(T);
        
    }
    // Enhancement Factor [-]
    f=f_factor(T,p);
    
    // Saturation pressure [kPa]
    p_s=f*p_ws;
    
    if (HumInput==GIVEN_RH)
    {
        RH=xSI;
        W=epsilon*RH*p_s/(p-RH*p_s);
        psi_w=W/(epsilon+W);
    }
    else if (HumInput==GIVEN_TDP)
    {
        Tdp=xSI;
        // Saturation pressure at dewpoint [kPa]
        p_ws_dp=Props('P','T',Tdp,'Q',0,"Water");
        // Enhancement Factor at dewpoint temperature [-]
        f_dp=f_factor(Tdp,p);
        // Water vapor pressure at dewpoint [kPa]
        p_w_dp=f_dp*p_ws_dp;
        // Water mole fraction [-]
        psi_w=p_w_dp/p;
        // Humidity ratio [-]
        W=psi_w*epsilon/(1-psi_w);
    }
    else if (HumInput==GIVEN_TWB)
    {
        // Not implemented yet
        printf("Wet-bulb input not implemented yet");
    }
    else if (HumInput==GIVEN_HUMRAT) //(2)
    {
        W=xSI;
        psi_w=W/(epsilon+W);
    }
    
    M_ha=(1-psi_w)*28.966+MM_Water()*psi_w;

    // Find relative humidity using W/e=phi*p_s/(p-phi*p_s)
    RH=W/epsilon*p/(p_s*(1+W/epsilon));
    
    // Get molar volume
    v_bar=MolarVolume(T,p,psi_w);
    
    // Get molar enthalpy
    hbar=MolarEnthalpy(T,p,psi_w);
    
    // Get dewpoint temperature
    Tdp=DewpointTemperature(T,p,psi_w);
    
    //~ //Get wet-bulb temperature
    //~ Twb=WetbulbTemperature(T,p,psi_w);
    
    // Specific volume v_da in m^3/kg_da, v_bar has units of m^3/mol
    v_da=v_bar*(1+W)/M_ha*1000;
    
    *v_out=v_da;
    *h_out=hbar*(1+W)/M_ha;
    *W_out=W;
    *RH_out=RH;
    *Tdp_out=Tdp;
    //~ *Twb_out=Twb;

	return 1;
}

double HumAir_Single(double T, double pSI, char *HumInputStr, double xSI, char *OutputStr)
{
	int HumInput;
	double Tdp_out, w_out, h_out, RH_out, v_out;

	if (!strcmp(HumInputStr,"RH"))
		HumInput=GIVEN_RH;
	else if (!strcmp(HumInputStr,"Omega"))
		HumInput=GIVEN_HUMRAT;
	else if (!strcmp(HumInputStr,"DewPoint"))
		HumInput=GIVEN_TDP;
	else if (!strcmp(HumInputStr,"WetBulb"))
		HumInput=GIVEN_TWB;
	else if (!strcmp(HumInputStr,"Enthalpy"))
		HumInput=GIVEN_ENTHALPY;
	else
	{
		fprintf(stderr,"Type of input [%s] not acceptable in HumAir_Single. Valid values are Omega,RH,WetBulb,DewPoint,Enthalpy.",HumInputStr);
		return -1;
	}
	HumAir(T,pSI,HumInput,xSI,&Tdp_out,&w_out,&h_out,&RH_out,&v_out);
	if (!strcmp(OutputStr,"RH"))
		return RH_out;
	else if (!strcmp(OutputStr,"Omega"))
		return w_out;
	else if (!strcmp(OutputStr,"DewPoint"))
		return Tdp_out;
	else if (!strcmp(OutputStr,"Enthalpy"))
		return h_out;
	else if (!strcmp(OutputStr,"Volume"))
		return v_out;
	else
	{
		fprintf(stderr,"Type of output [%s] not acceptable in HumAir_Single. Valid values are RH,Omega,DewPoint,Enthalpy,Volume.",OutputStr);
		return -1;
	}

}
void HumAir(double T, double pSI, int HumInput, double xSI, /* in --- out */ double *Tdp_out, double *w_out, double *h_out, double *RH_out, double *v_out)
{
	double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f,RH;

	int iter=1;

	//HumAir1(T,pSI,HumInput,xSI,Tdp_out,w_out,h_out,RH_out,v_out);

	// The humid air input is something that is not enthalpy
	if (HumInput==GIVEN_TDP || HumInput==GIVEN_HUMRAT || HumInput==GIVEN_TWB
		|| HumInput==GIVEN_RH )
	{
		HumAir_func(T,pSI,HumInput,xSI,Tdp_out,w_out,h_out,RH_out,v_out);
	}
	else
	{
		// Use a secant solve to find the relative humidity which gives you the enthalpy you desire
		while ((iter<=3 || change>eps) && iter<100)
		{
			if (iter==1){x1=0.2; RH=x1;}
			if (iter==2){x2=0.6; RH=x2;}
			if (iter>2) {RH=x2;}
				HumAir_func(T,pSI,GIVEN_RH,RH,Tdp_out,w_out,h_out,RH_out,v_out);
				f=*h_out-xSI;
			if (iter==1){y1=f;}
			if (iter>1)
			{
				y2=f;
				x3=x2-y2/(y2-y1)*(x2-x1);
				change=fabs(y2/(y2-y1)*(x2-x1));
				y1=y2; x1=x2; x2=x3;
			}
			iter=iter+1;
		}
		HumAir_func(T,pSI,GIVEN_RH,RH,Tdp_out,w_out,h_out,RH_out,v_out);
	}
	return;
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

//~ double hair_sat(double T)
//~ {
	//~ // Saturated air enthalpy
	//~ // Based on a correlation from EES, good from 250K to 300K.  
	//~ // No error bound checking is carried out
	//~ // T: [K]
	//~ // hair_sat: [kJ/kg]
	//~ return 3.70618205E+04-5.69352831E+02*T+3.27662170E+00*T*T-8.39389630E-03*T*T*T+8.09476615E-06*T*T*T*T;
//~ }

double T_hss(double h, double p, double T_guess)
{
	// Function finds the temperature that yields 
	// the saturated surface enthalpy using the secant method
	// T[K],h[kJ/kg],p[kPa]
	double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f,T=300;
	double Tdp_dummy, omega_dummy, h_s, RH_dummy, v_dummy;

	int iter=1;

	while ((iter<=3 || change>eps) && iter<100)
	{
		if (iter==1){x1=T_guess; T=x1;}
		if (iter==2){x2=T_guess+0.1; T=x2;}
		if (iter>2) {T=x2;}
			HumAir(T, p, GIVEN_RH, 1.0, &Tdp_dummy, &omega_dummy, &h_s, &RH_dummy, &v_dummy);
			f=h_s-h;
		if (iter==1){y1=f;}
		if (iter>1)
		{
			y2=f;
			x3=x2-y2/(y2-y1)*(x2-x1);
			change=fabs(y2/(y2-y1)*(x2-x1));
			y1=y2; x1=x2; x2=x3;
		}
		iter=iter+1;
	}
	return T;
}

double T_homega(double h, double omega, double P, double T_guess, double deltaT){

	// To compute the temperature based on enthalpy and humidity ratio
	double T_dummy, rh_dummy, D_dummy, h_dummy, v_dummy, omega_dummy, h_dummy2;
	double error = 1;

	while(error>0.0005){
		HumAir(T_guess, P, GIVEN_HUMRAT, omega, &D_dummy, &omega_dummy, &h_dummy, &rh_dummy, &v_dummy);
		HumAir(T_guess+0.001, P, GIVEN_HUMRAT, omega, &D_dummy, &omega_dummy, &h_dummy2, &rh_dummy, &v_dummy);
		T_dummy = T_guess - (h_dummy - h)/(h_dummy2-h_dummy)*0.001;
		error = fabs(T_dummy-T_guess);
		T_guess = T_dummy;
	}
	return T_dummy;
}

static double powI(double x, int y)
{
    int k;
    double product=1.0;
    double x_in;
    int y_in;
    
    if (y==0)
    {
        return 1.0;
    }
    
    if (y<0)
    {
        x_in=1/x;
        y_in=-y;
    }
	else
	{
		x_in=x;
		y_in=y;
	}

    if (y_in==1)
    {
        return x_in;
    }    
    
    product=x_in;
    for (k=1;k<y_in;k++)
    {
        product=product*x_in;
    }
    
    return product;
}