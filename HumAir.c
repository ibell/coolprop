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

 

void Help()
{
	printf(
		"First two inputs are always T [K], p [kPa absolute]\n"
		"Next input tells what the next input is. Codes available:\n"
		"%d : Dewpoint temperature [K]\n"
		"%d : Humidity ratio [kg water / kg dry air]\n"
		"%d : Given wet bulb temperature [K] \n"
		"%d : Relative humidity in range [0-->1]\n"
		"%d : Enthalpy in [kJ/kg-K]"
		,GIVEN_TDP,GIVEN_HUMRAT,GIVEN_TWB,GIVEN_RH,GIVEN_ENTHALPY);

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

static int HumAir1(double tSI, double pSI, int HumInput, double xSI, /* in --- out */ double *Tdp_out, double *W_out, double *h_out, double *RH_out, double *v_out)
{
	double t,p,T,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18;
	double W,td,Td,pw,pws,rh,a,tstar,Tstar,pwsstar,Wsstar,hSI,v;
	tSI = tSI-273.15;       // Conversion from K to °C
	t = tSI*1.8+32;         // Conversion from °C to °F
	p = pSI/6.894757;       // Conversion from kPa to psi

	T = t + 459.67;         // Conversion from Fahrenheit to Rankine scales

	// Constants for vapor pressure and dewpoint calculations (#5, #6 and #40)
	c1 = -1.0214165e04;
	c2 = -4.8932428;
	c3 = -5.3765794e-03;
	c4 = 1.9202377e-07;
	c5 = 3.5575832e-10;
	c6 = -9.0344688e-14;
	c7 = 4.1635019;
	c8 = -1.0440397e04;
	c9 = -1.1294650e01;
	c10= -2.7022355e-02;
	c11= 1.2890360e-05;
	c12= -2.4780681e-09;
	c13= 6.5459673;
	c14= 100.45;
	c15= 33.193;
	c16= 2.319;
	c17= 0.17074;
	c18= 1.2063;

	//========================================================================
	// GIVEN DEWPOINT
	if (HumInput == 1)
	{
	    td = xSI * 1.8 + 32; // Conversion of dewpoint from °C to °F
	    Td = td + 459.67;   // Conversion of dewpoint from Fahrenheit to Rankine scales
	    // Some relationships have different formulations for above or below water
	    // freezing temperature
	    if (t < 32)
		{
	        pws = exp(c1/T + c2 + c3*T + c4* T*T + c5 * T*T*T + c6 * T*T*T*T*T + c7*log(T));
	        pw = exp(c1/Td + c2 + c3*Td + c4* Td*Td + c5 * Td*Td*Td + c6 * Td*Td*Td*Td*Td + c7*log(Td));
		}
	    else
		{
	        pws = exp(c8/T + c9 + c10*T + c11* T*T + c12 * T*T*T + c13*log(T));
	        pw = exp(c8/Td + c9 + c10*Td + c11* Td*T + c12 * Td*Td*Td + c13*log(Td));
		}
	    W = 0.621945*(pw/(p - pw));     // humidity ratio (#22)
	    rh = pw/pws;                    // relative humidity calculation (#24)
	}
    //========================================================================
    // GIVEN HUMIDITY RATIO
	if (HumInput == 2)
	{
	    W = xSI;
	    pw = p*W/(0.621945 + W);         // Water vapor partial pressure (derived from #22)
	    // Some relationships have formulations for above or below water
	    // freezing temperature
	    a = log(pw);                     // Variable alpha from ASHRAE Equation #40
	    if (t < 32)
		{
	        pws = exp(c1/T + c2 + c3*T + c4* T*T + c5 * T*T*T + c6 * T*T*T*T + c7*log(T));
	        td = 90.12 + 26.142*a + 0.8927*a*a;  // Dewpoint (#40)
		}
	    else
		{
	        pws = exp(c8/T + c9 + c10*T + c11* T*T + c12 * T*T*T + c13*log(T));
	        td = c14 + c15*a + c16*a*a + c17*a*a*a + c18*pow(pw,0.1984); // Dewpoint (#39)
		}
	    rh = pw/pws;                     // relative humidity calculation (#24)
	}
    //========================================================================
    // GIVEN WET BULB TEMPERATURE
	if (HumInput == 3)
	{
	    xSI = xSI - 273.15;               // Conversion from K to °C
	    tstar = xSI*1.8+32;                // Conversion from °C to °F   
	    Tstar = tstar + 459.67;          // Express wetbulb temperature in Rankine
	    // Some relationships have formulations for above or below water
	    // freezing temperature
	    if (t < 32)
		{
	        pws = exp(c1/T + c2 + c3*T + c4* T*T + c5 * T*T*T + c6 * T*T*T*T + c7*log(T));
	        pwsstar = exp(c1/Tstar + c2 + c3*Tstar + c4* Tstar*Tstar + c5 * Tstar*Tstar*Tstar + c6 * Tstar*Tstar*Tstar*Tstar + c7*log(Tstar));
	        Wsstar = 0.621945*(pwsstar/(p - pwsstar));  // humidity ratio at saturation (#23)
	        W = ((1220 - 0.04*tstar) * Wsstar - 0.24 * (t - tstar))/(1220 + 0.444*t - 0.48*tstar);  // (#35)
	        pw = p*W/(0.621945 + W);     // Water vapor partial pressure (derived from #22)
	        a = log(pw);                 // Variable alpha from ASHRAE Equation #40
	        td = 90.12 + 26.142*a + 0.8927*a*a;  // Dewpoint (#40)
		}
	    else
		{
	        pws = exp(c8/T + c9 + c10*T + c11* T*T + c12 * T*T*T + c13*log(T));
	        pwsstar = exp(c8/Tstar + c9 + c10*Tstar + c11* Tstar*Tstar + c12 * Tstar*Tstar*Tstar + c13*log(Tstar));
	        Wsstar = 0.621945*(pwsstar/(p - pwsstar));  // humidity ratio at saturation (#23)
	        W = ((1093 - 0.556*tstar) * Wsstar - 0.24 * (t - tstar))/(1093 + 0.444*t - tstar);  // (#35)
	        pw = p*W/(0.621945 + W);     // Water vapor partial pressure (derived from #22)
	        a = log(pw);                 // Variable alpha from ASHRAE Equation #40
	        td = c14 + c15*a + c16*a*a + c17*a*a*a + c18*pow(pw,0.1984); // Dewpoint (#39)
		}
	    rh = pw/pws;                     // relative humidity calculation (#24)
	}
    //========================================================================
    // GIVEN RELATIVE HUMIDITY
	if (HumInput == 4)
	{
	    rh = xSI;
	    // Some relationships have formulations for above or below water
	    // freezing temperature
	    if (t < 32)
		{
	        pws = exp(c1/T + c2 + c3*T + c4* T*T + c5 * T*T*T + c6 * T*T*T*T + c7*log(T));
	        pw = rh*pws;    // Water vapor partial pressure (derived from #24)
	        a = log(pw);    // Variable alpha from ASHRAE Equation #40
	        td = 90.12 + 26.142*a + 0.8927*a*a;  // Dewpoint (#40)
		}
	    else
		{
	        pws = exp(c8/T + c9 + c10*T + c11* T*T + c12 * T*T*T + c13*log(T));
	        pw = rh*pws;   // Water vapor partial pressure (derived from #24)
	        a = log(pw);   // Variable alpha from ASHRAE Equation #40
	        td = c14 + c15*a + c16*a*a + c17*a*a*a + c18*pow(pw,0.1984); // Dewpoint (#39)
		}
	    W = 0.621945*(pw/(p - pw));      // humidity ratio (#22)
	}
	// GIVEN ENTHALPY [kJ/kg]
	if (HumInput == 5)
	{
	    hSI = xSI;
	    // Some relationships have formulations for above or below water
	    // freezing temperature
	    if (t < 32)
		{
	        pws = exp(c1/T + c2 + c3*T + c4* T*T + c5 * T*T*T + c6 * T*T*T*T + c7*log(T));
	        W = (hSI-1.006*(tSI+273.15))/(2501 + 1.86*(tSI+273.15)); // humidity ratio (derived from #32)
	        pw = W*p/(.621945+W); // Water vapor partial pressure (derived from #22)
	        a = log(pw);    // Variable alpha from ASHRAE Equation #40
	        td = 90.12 + 26.142*a + 0.8927*a*a;  // Dewpoint (#40)
		}
	    else
		{
	        pws = exp(c8/T + c9 + c10*T + c11* T*T + c12 * T*T*T + c13*log(T));
	        W = (hSI-1.006*(tSI+273.15))/(2501 + 1.86*(tSI+273.15)); // humidity ratio (derived from #32)
	        pw = W*p/(.621945+W); // Water vapor partial pressure (derived from #22)
	        a = log(pw);   // Variable alpha from ASHRAE Equation #40
	        td = c14 + c15*a + c16*a*a + c17*a*a*a + c18*pow(pw,0.1984); // Dewpoint (#39)
		}
	    rh = pw/pws;    // Water vapor partial pressure (derived from #24)
	}

	if (HumInput>5 || HumInput<1)
	{
		*h_out = -1;	
    	*Tdp_out = -1;			
    	*v_out = -1;
		*RH_out=-1;
		*W_out=-1;
		printf("Incorrect HumInput index");
		return -1;

	}

	tSI = tSI + 273.15; // Convert temperature from °C to K
    v = 0.370486 * T * (1 + 1.607858*W)/p; // Specific volume [ft^3]
	
	//Prepare outputs
	*h_out = 1.006*tSI + W*(2501 + 1.86*tSI);	// Specific enthalpy [kJ/kg] (#32)
    *Tdp_out = (td-32)/1.8 + 273.15;			// Convert dewpoint from °F to K
    *v_out = v*0.062428;						// Convert ft^3/lbmda to m^3/kgda
	*RH_out=rh;
	*W_out=W;

	return 0;
}

double calcpw(double T)
{
	double ln_pw;
	//T-=273.15;
	if (T<273.15)
	{
		ln_pw=-5674.5359/T+6.3925247-0.009677843*T+6.2215701e-7*T*T+2.0747825e-9*T*T*T-9.484024e-13*T*T*T*T+4.1635019*log(T);
	}
	else
	{
		// Fit from EES.  Between 0 and 200C, max relative error 0.05% on pressure
		ln_pw=-7.22649229E+01+7.16353354E-01*T-2.68371078E-03*T*T+5.51697367E-06*T*T*T-5.97500136E-09*T*T*T*T+2.67600257E-12*T*T*T*T*T;
	}
	return exp(ln_pw);
}


static int HumAir2(double T, double p, int HumInput, double xSI, /* in --- out */ double *Tdp_out, double *W_out, double *h_out, double *RH_out, double *v_out)
{
	double R_a,f,Z_ms,ln_pv,p_w,v_ms,x_as,x_ws,W_s,v_hs,h_ms,v_ss,s_ms,
		v_a,Z_a,v_ha,h_a,v_sa,s_a,RH,x_wm,x_am,W_m,mu,v_m,h_m,
		s_m,Tdp,x_as_dp,x_ws_dp,ln_pg,pv;

	// ------------------------
	// Saturated Air properties
	// ------------------------

	// Convert pressure in kPa to pressure in Pa
		p*=1000.0; //[kPa --> Pa]

	// Ideal gas constant for air
		R_a=0.28702899; //[kJ/kg-K]

	// Enhancement factor
		f=2.2770286-0.02406584*T+1.8213945e-4*T*T-6.8894708e-7*T*T*T+1.297668e-9*T*T*T*T-9.7078508e-13*T*T*T*T*T+3.9945654e-8*p;
	// Saturated air compressibility factor
		Z_ms=1.9208388-0.018226313*T+1.4278754e-4*T*T-5.5526317e-7*T*T*T+1.073933e-9*T*T*T*T-8.2740746e-13*T*T*T*T*T+3.9755302e-9*p;

		p_w=calcpw(T);

		v_ms=R_a*T/(p-f*p_w)*Z_ms*1000.0; // 1000.0 from the unit conversion from kJ to J in numerator
		x_as=(p-f*p_w)/p;
		x_ws=f*p_w/p;
		W_s=0.62198*x_ws/x_as;

		v_hs=-4.1877292e-18*p*p*p+1.8317367e-12*p*p-3.0980589e-6*p+0.28801940;
		h_ms=1.0041923*(T-273.15)+W_s*(1.8642569*(T-273.15)+2500.7876)+v_hs;
	
		v_ss=1.3900054e-23*p*p*p*p-6.5157721e-18*p*p*p+1.1762918e-12*p*p-1.0177460e-7*p+3.6530569e-3;
		s_ms=(1.0041923*log(T)-5.63354)+W_s*(1.8642569*log(T)-3.661889)-R_a/x_as*log(p/101325)+R_a*log(Z_ms/x_as)+x_ws/x_as*R_a*log(Z_ms/x_ws)+v_ss;


	// ------------------------
	// Dry Air Properties
	// ------------------------

		Z_a=0.97901976+1.8181227e-4*T-5.2676556e-7*T*T+5.2528153e-10*T*T*T-6.6519834e-9*p;
		v_a=R_a*T/p*Z_a*1000.0; // 1000.0 from the unit conversion from kJ to J in numerator
		
		v_ha=3.1327784-3.1755794e-2*T+1.2006509e-4*T*T-1.5693167e-7*T*T*T-2.7847739e-6*p/1000+6.951284e-14*p/1000/1000*p;
		h_a=1.0041923*(T-273.15)+v_ha;

		v_sa=5.990033e-3-5.0193469e-5*T+2.2352256e-7*T*T-3.2817759e-10*T*T*T-7.4664626e-8*p+9.7299467e-13*p*p-7.3212392e-18*p*p*p+2.9856674e-23*p*p*p*p-5.1169449e-18*p*p*p*p*p;
		s_a=1.0041923*log(T)-5.63354-R_a*log(p/101325)+R_a*log(Z_a)+v_sa;

	//--------------------------
	// Mixture Properties
	// -------------------------

		if (HumInput==GIVEN_RH) //(4)
		{
			RH=xSI;
			x_wm=RH*x_ws;
			x_am=1-x_wm;
			W_m=0.62198*x_wm/x_am;
			mu=W_m/W_s;
		}
		else if (HumInput==GIVEN_TDP)
		{
			Tdp=xSI;
			// Enhancement factor
				f=2.2770286-0.02406584*T+1.8213945e-4*T*T-6.8894708e-7*T*T*T+1.297668e-9*T*T*T*T-9.7078508e-13*T*T*T*T*T+3.9945654e-8*p;
			p_w=calcpw(Tdp);
			x_as_dp=(p-f*p_w)/p;
			x_ws_dp=f*p_w/p;
			W_m=0.62198*x_ws_dp/x_as_dp;
			mu=W_m/W_s;
			x_wm=W_m/(0.62198+W_m);
		}
		else if (HumInput==GIVEN_TWB)
		{
			// Not implemented yet
			printf("Wet-bulb input not implemented yet");
		}
		else if (HumInput==GIVEN_HUMRAT) //(2)
		{
			W_m=xSI;
			mu=W_m/W_s;
			x_wm=W_m/(0.62198+W_m);
		}


		RH=x_wm/x_ws;
		if (T>273.15)
		{
			ln_pg=-7.44159509E+01
				+7.59005181E-01*T
				-3.02612526E-03*T*T
				+6.94585158E-06*T*T*T
				-9.25199787E-09*T*T*T*T
				+6.60133477E-12*T*T*T*T*T
				-1.92266328E-15*T*T*T*T*T*T;
		}
		else
		{
			ln_pg=1.18874020E+03
				-3.05184253E+01*T
				+3.17994863E-01*T*T
				-1.74454352E-03*T*T*T
				+5.35287821E-06*T*T*T*T
				-8.73024361E-09*T*T*T*T*T
				+5.91855083E-12*T*T*T*T*T*T;
		}

		pv=RH*exp(ln_pg);
		ln_pv=log(pv);
		if (p_w<603.7)
		{
			// solid-vapor equilibrium
			// Refprop sublimation pressure of solid-gas taken from 215 K to 273.15 K, then fit in EES. 
			// Below 220K, max error in dewpoint pressure is 0.14%, above 220K, max error 0.04%
			Tdp=2.12555164E+02
				+7.42023501E+00*ln_pv
				+1.84598611E-01*ln_pv*ln_pv
				+3.98898424E-02*ln_pv*ln_pv*ln_pv
				-6.96777141E-03*ln_pv*ln_pv*ln_pv*ln_pv
				+8.56960521E-04*ln_pv*ln_pv*ln_pv*ln_pv*ln_pv
				-3.80010431E-05*ln_pv*ln_pv*ln_pv*ln_pv*ln_pv*ln_pv;
		}
		else
		{
			// vapor-mixture equilibrium
			// Data fit from EES
			// Max error is 0.0138 K in the range 273.15 K to 200 K
			Tdp=2.38239452E+02
				-1.33711768E+01*ln_pv
				+6.25997592E+00*ln_pv*ln_pv
				-8.88780833E-01*ln_pv*ln_pv*ln_pv
				+7.71822220E-02*ln_pv*ln_pv*ln_pv*ln_pv
				-3.49670180E-03*ln_pv*ln_pv*ln_pv*ln_pv*ln_pv
				+7.14014963E-05*ln_pv*ln_pv*ln_pv*ln_pv*ln_pv*ln_pv;
		}

		v_m=v_a+mu*(v_ms-v_a);
		h_m=h_a+mu*(h_ms-h_a);
		s_m=s_a+mu*(s_ms-s_a);

		*v_out=v_m;
		*h_out=h_m;
		*W_out=W_m;
		*RH_out=x_wm/x_ws;
		*Tdp_out=Tdp;

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
int HumAir(double T, double pSI, int HumInput, double xSI, /* in --- out */ double *Tdp_out, double *w_out, double *h_out, double *RH_out, double *v_out)
{
	double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f,RH;

	int iter=1;

	//HumAir1(T,pSI,HumInput,xSI,Tdp_out,w_out,h_out,RH_out,v_out);

	// The humid air input is something that is not enthalpy
	if (HumInput==GIVEN_TDP || HumInput==GIVEN_HUMRAT || HumInput==GIVEN_TWB
		|| HumInput==GIVEN_RH )
	{
		HumAir2(T,pSI,HumInput,xSI,Tdp_out,w_out,h_out,RH_out,v_out);
	}
	else
	{
		// Use a secant solve to find the relative humidity which gives you the enthalpy you desire
		while ((iter<=3 || change>eps) && iter<100)
		{
			if (iter==1){x1=0.2; RH=x1;}
			if (iter==2){x2=0.6; RH=x2;}
			if (iter>2) {RH=x2;}
				HumAir2(T,pSI,GIVEN_RH,RH,Tdp_out,w_out,h_out,RH_out,v_out);
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
		HumAir2(T,pSI,GIVEN_RH,RH,Tdp_out,w_out,h_out,RH_out,v_out);
	}
	return 1;
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

double hair_sat(double T)
{
	// Saturated air enthalpy
	// Based on a correlation from EES, good from 250K to 300K.  
	// No error bound checking is carried out
	// T: [K]
	// hair_sat: [kJ/kg]
	return 3.70618205E+04-5.69352831E+02*T+3.27662170E+00*T*T-8.39389630E-03*T*T*T+8.09476615E-06*T*T*T*T;
}

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
	int flag;
	double error = 1;

	while(error>0.0005){
		flag = HumAir(T_guess, P, GIVEN_HUMRAT, omega, &D_dummy, &omega_dummy, &h_dummy, &rh_dummy, &v_dummy);
		flag = HumAir(T_guess+0.001, P, GIVEN_HUMRAT, omega, &D_dummy, &omega_dummy, &h_dummy2, &rh_dummy, &v_dummy);
		T_dummy = T_guess - (h_dummy - h)/(h_dummy2-h_dummy)*0.001;
		error = fabs(T_dummy-T_guess);
		T_guess = T_dummy;
	}
	return T_dummy;
}