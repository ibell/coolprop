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


int HumAir(double tSI, double pSI, int HumInput, double xSI, /* in --- out */ double *Tdp_out, double *W_out, double *h_out, double *RH_out, double *v_out)
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

	if (HumInput == 5)
	{
	    hSI = xSI;
	    // Some relationships have formulations for above or below water
	    // freezing temperature
	    if (t < 32)
		{
	        pws = exp(c1/T + c2 + c3*T + c4* T*T + c5 * T*T*T + c6 * T*T*T*T + c7*log(T));
	        W = (hSI-1.006*tSI)/(2501 + 1.86*tSI); // humidity ratio (derived from #32)
	        pw = W*p/(.621945+W); // Water vapor partial pressure (derived from #22)
	        a = log(pw);    // Variable alpha from ASHRAE Equation #40
	        td = 90.12 + 26.142*a + 0.8927*a*a;  // Dewpoint (#40)
		}
	    else
		{
	        pws = exp(c8/T + c9 + c10*T + c11* T*T + c12 * T*T*T + c13*log(T));
	        W = (hSI-1.006*tSI)/(2501 + 1.86*tSI); // humidity ratio (derived from #32)
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

	*h_out = 1.006*tSI + W*(2501 + 1.86*tSI);	// Specific enthlpay [kJ/kg] (#32)
    *Tdp_out = (td-32)/1.8 + 273.15;			// Convert dewpoint from °F to K
    *v_out = v*0.062428;						// Convert ft^3/lbmda to m^3/kgda
	*RH_out=rh;
	*W_out=W;

	return 0;
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
	// the saturated surface temperature using the secant method
	// T[K],h[kJ/kg],p[kPa]
	double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f,T=300;
	double Tdp_dummy, omega_dummy, h_s, RH_dummy, v_dummy;

	int iter=1;

	while ((iter<=3 || change>eps) && iter<100)
	{
		if (iter==1){x1=T_guess; T=x1;}
		if (iter==2){x2=T_guess+0.1; T=x2;}
		if (iter>2) {T=x2;}
			HumAir(T, p, 4, 1.0, &Tdp_dummy, &omega_dummy, &h_s, &RH_dummy, &v_dummy);
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