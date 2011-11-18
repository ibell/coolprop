#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#endif

#if defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#include "REFPROP.h"
#endif

#include <stdlib.h>
#include "string.h"
#include <stdio.h>
#include "CoolProp.h"

// Some constants for REFPROP... defined by macros for ease of use 
#define refpropcharlength 255
#define filepathlength 255
#define lengthofreference 3
#define errormessagelength 255
#define ncmax 20		// Note: ncmax is the max number of components
#define numparams 72 
#define maxcoefs 50

char LoadedREFPROPRef[255];
int FluidType;

#define FLUIDTYPE_REFPROP 0
#define FLUIDTYPE_BRINE 1
#define FLUIDTYPE_REFRIGERANT_PURE 2
#define FLUIDTYPE_REFRIGERANT_PSEUDOPURE 3

double R_u=8.314472; //From Lemmon et al 2009 (propane)
// Function prototypes
double _Dekker_Tsat(double x_min, double x_max, double eps, double p, double Q,char *Ref);

// Global residual Helmholtz derivatives function pointers
double (*p_func)(double,double);
double (*h_func)(double,double,int);
double (*s_func)(double,double,int);
double (*u_func)(double,double,int);
double (*rho_func)(double,double,int);
double (*cp_func)(double,double,int);
double (*cv_func)(double,double,int);
double (*visc_func)(double,double,int);
double (*k_func)(double,double,int);
double (*w_func)(double,double,int);
double (*MM_func)(void);
double (*rhocrit_func)(void);
double (*Tcrit_func)(void);
double (*pcrit_func)(void);
double (*Ttriple_func)(void);

//For the pure fluids
double (*psat_func)(double); 
// For the psedo-pure fluids
double (*pdp_func)(double);
double (*pbp_func)(double);

// Global residual Helmholtz derivatives function pointers
double (*phir_func)(double,double);
double (*dphir_dDelta_func)(double,double);
double (*dphir2_dDelta2_func)(double,double);
double (*dphir2_dDelta_dTau_func)(double,double);
double (*dphir_dTau_func)(double,double);
double (*dphir2_dTau2_func)(double,double);
double (*phi0_func)(double,double);
double (*dphi0_dDelta_func)(double,double);
double (*dphi02_dDelta2_func)(double,double);
double (*dphi0_dTau_func)(double,double);
double (*dphi02_dTau2_func)(double,double);

double (*rhosatV_func)(double);
double (*rhosatL_func)(double);


static double QuadInterpolate(double x0, double x1, double x2, double f0, double f1, double f2, double x)
{
    double L0, L1, L2;
    L0=((x-x1)*(x-x2))/((x0-x1)*(x0-x2));
    L1=((x-x0)*(x-x2))/((x1-x0)*(x1-x2));
    L2=((x-x0)*(x-x1))/((x2-x0)*(x2-x1));
    return L0*f0+L1*f1+L2*f2;
}

void MatInv_2(double A[2][2] , double B[2][2])
{
	double Det;
	//Using Cramer's Rule to solve

	Det=A[0][0]*A[1][1]-A[1][0]*A[0][1];
	B[0][0]=1.0/Det*A[1][1];
	B[1][1]=1.0/Det*A[0][0];
	B[1][0]=-1.0/Det*A[1][0];
	B[0][1]=-1.0/Det*A[0][1];
}

double Pressure_Trho(double T, double rho)
{
	double delta,tau,R;
	R=R_u/MM_func();
	tau=Tcrit_func()/T;
	delta=rho/rhocrit_func();
    return R*T*rho*(1.0+delta*dphir_dDelta_func(tau,delta));
}
double IntEnergy_Trho(double T, double rho)
{
    double delta,tau,R;
    R=R_u/MM_func();
	tau=Tcrit_func()/T;
	delta=rho/rhocrit_func();
    return R*T*tau*(dphi0_dTau_func(tau,delta)+dphir_dTau_func(tau,delta));
}
double Enthalpy_Trho(double T, double rho)
{
    double delta,tau,R;
    R=R_u/MM_func();
	tau=Tcrit_func()/T;
	delta=rho/rhocrit_func();
    return R*T*(1.0+tau*(dphi0_dTau_func(tau,delta)+dphir_dTau_func(tau,delta))+delta*dphir_dDelta_func(tau,delta));
}
double Entropy_Trho(double T, double rho)
{
    double delta,tau,R;
    R=R_u/MM_func();
	tau=Tcrit_func()/T;
	delta=rho/rhocrit_func();
    return R*(tau*(dphi0_dTau_func(tau,delta)+dphir_dTau_func(tau,delta))-phi0_func(tau,delta)-phir_func(tau,delta));
}
//static double SpecHeatV_Trho(double T, double rho)
//{
//    double delta,tau;
//    delta=rho/rhoc;
//    tau=Tc/T;
//    
//    return -R_R290*powI(tau,2)*(dphi02_dTau2_R290(tau,delta)+dphir2_dTau2_R290(tau,delta));
//}
//static double SpecHeatP_Trho(double T, double rho)
//{
//    double delta,tau,c1,c2;
//    delta=rho/rhoc;
//    tau=Tc/T;
//
//    c1=powI(1.0+delta*dphir_dDelta_R290(tau,delta)-delta*tau*dphir2_dDelta_dTau_R290(tau,delta),2);
//    c2=(1.0+2.0*delta*dphir_dDelta_R290(tau,delta)+powI(delta,2)*dphir2_dDelta2_R290(tau,delta));
//    return R_R290*(-powI(tau,2)*(dphi02_dTau2_R290(tau,delta)+dphir2_dTau2_R290(tau,delta))+c1/c2);
//}
//
//static double SpeedSound_Trho(double T, double rho)
//{
//    double delta,tau,c1,c2;
//    delta=rho/rhoc;
//    tau=Tc/T;
//
//    c1=-SpecHeatV_Trho(T,rho)/R_R290;
//    c2=(1.0+2.0*delta*dphir_dDelta_R290(tau,delta)+powI(delta,2)*dphir2_dDelta2_R290(tau,delta));
//    return sqrt(-c2*T*SpecHeatP_Trho(T,rho)*1000/c1);
//}

double Gibbs_Trho(double T,double rho)
{
	return Enthalpy_Trho(T,rho)-T*Entropy_Trho(T,rho);
}

double Density_Tp(double T, double p, double rho)
{
	double delta,tau,dpdrho,error=999,R,rhoc;
	R=R_u/MM_func();
	tau=Tcrit_func()/T;
	rhoc=rhocrit_func();
	while (fabs(error)>1e-7)
	{
		delta=rho/rhoc;
		// Use Newton's method to find the saturation density since the derivative of pressure w.r.t. density is known from EOS
		dpdrho=R*T*(1+2*delta*dphir_dDelta_func(tau,delta)+delta*delta*dphir2_dDelta2_func(tau,delta));
		// Update the step using Newton's method
		rho=rho-(Pressure_Trho(T,rho)-p)/dpdrho;
		// Residual
		error=Pressure_Trho(T,rho)-p;
	}		
	return rho;
}

void rhosatPure(char *Ref, double T, double *rhoLout, double *rhoVout, double *pout)
{
	// Only works for pure fluids (no blends)
	// At equilibrium, saturated vapor and saturated liquid are at the same pressure and the same Gibbs energy
	double rhoL,rhoV,p,error=999,x1,x2,x3,y1,y2,f,p_guess;
	int iter;

	// Call a property function to make sure that the pointers to Helmholtz functions are updated
	Props('B','T',0,'P',0,Ref); //Critical temperature

	// Use the density ancillary function as the starting point for the secant solver
	rhoL=rhosatL_func(T);
	rhoV=rhosatV_func(T);
	p_guess=Pressure_Trho(T,rhoV);

	iter=1;
	// Use a secant method to obtain pressure
	while ((iter<=3 || fabs(error)>1e-7) && iter<100)
	{
		if (iter==1){x1=p_guess; p=x1;}
		if (iter==2){x2=1.0001*p_guess; p=x2;}
		if (iter>2) {p=x2;}
			//Recalculate the densities based on the current pressure
			rhoL=Density_Tp(T,p,rhoL);
			rhoV=Density_Tp(T,p,rhoV);
			// Residual between saturated liquid and saturated vapor gibbs function
			f=Gibbs_Trho(T,rhoL)-Gibbs_Trho(T,rhoV);
		if (iter==1){y1=f;}
		if (iter>1)
		{
			y2=f;
			x3=x2-y2/(y2-y1)*(x2-x1);
			error=f;
			y1=y2; x1=x2; x2=x3;
		}
		iter=iter+1;
		if (iter>100)
		{
			printf("rhosatPure failed, current values of rhoL and rhoV are %g,%g\n",rhoL,rhoV);
			return;
		}
	}
	*rhoLout=rhoL;
	*rhoVout=rhoV;
	*pout=p;
	return;
}

double SecFluids(char Output, double T, double p,char * Ref)
{
	double Tfreeze,Tmax,rho,cp,k,mu,u,s,TC,C_gly;
	// Temperature and Pressure are the inputs

	if (!strcmp(Ref,"HC-10"))
	{
		// Curve fits generated from Microsoft Excel fit of manufacturer data
		// Input temperature is in deg C
		TC=T-273.15;
		switch(Output)
		{
			case 'D':
				return -4.52609118E-01*TC + 1.19919457E+03;
			case 'C':
				return 2.4797494E-03*TC + 3.2708330E+00;
			case 'L':
				return 1.00000E-06*TC + 5.04400E-04;
			case 'V':
				return 2.937072591E-12*TC*TC*TC*TC - 1.713938938E-09*TC*TC*TC + 3.826311605E-07*TC*TC - 4.253611683E-05*TC + 2.509839772E-03;
			case 'F':
				return 263.15;
			case 'M':
				return 390;
			case 'H':
				return (2.4797494E-03*TC*TC/2.0 + 3.2708330E+00*TC)/1000+p/(-4.52609118E-01*TC + 1.19919457E+03);
			case 'S':
				return 2.4797494E-03*(T-298.15) + 3.2708330E+00*log(T/298.15);
			default:
				return _HUGE;
		}
	}


	// Ethylene glycol blends and pure water
	// "EG-10%" is 10% by mass ethylene glycol
	else if (Ref[0]=='E' && Ref[1]=='G')
	{
		switch (Ref[3])
		{
			case '1':
				C_gly=10; break;
			case '2':
				C_gly=20; break;
			case '3':
				C_gly=30; break;
			case '4':
				C_gly=40; break;
			case '5':
				C_gly=50; break;
			default:
				return _HUGE;
		}

		/*
		Brine() takes inputs of temperature [C]  and secondary fluid concentration in water [mass %]
		Outputs are freezing temperature [C], density [kg/m^3], Specific Heat [J/kg-K], Conductivity [W/m-K], Viscosity [mPa-s]
		*/
		
		// Set the temperature to within the band to avoid zero
		if (Output=='M' || Output=='F'){ T=300;}
		Brine("EG", T - 273.15, C_gly, &Tfreeze, &Tmax, &rho, &cp, &k, &mu,&u,&s);
		switch(Output)
		{
			case 'D':
				return rho;
			case 'C':
				return cp/1000;
			case 'L':
				return k/1000;
			case 'V':
				return mu / 1000.0;
			case 'F':
				return Tfreeze+273.15;
			case 'M':
				return Tmax+273.15;
			case 'H':
				return u/1000+p/rho;
			case 'S':
				return s/1000;
			default:
				return _HUGE;
		}
	}

	// Propylene glycol blends and pure water
	// "PG-10%" is 10% by mass propylene glycol
	else if (Ref[0]=='P' && Ref[1]=='G')
	{
		switch (Ref[3])
		{
			case '1':
				C_gly=10; break;
			case '2':
				C_gly=20; break;
			case '3':
				C_gly=30; break;
			case '4':
				C_gly=40; break;
			case '5':
				C_gly=50; break;
			default:
				return _HUGE;
		}

		/*
		Brine() takes inputs of temperature [C]  and secondary fluid concentration in water [mass %]
		Outputs are freezing temperature [C], density [kg/m^3], Specific Heat [J/kg-K], Conductivity [W/m-K], Viscosity [mPa-s]
		*/
		// Set the temperature to within the band to avoid zero
		if (Output=='M' || Output=='F'){ T=300;}
		Brine("PG", T - 273.15, C_gly, &Tfreeze, &Tmax, &rho, &cp, &k, &mu,&u,&s);
		switch(Output)
		{
			case 'D':
				return rho;
			case 'C':
				return cp/1000;
			case 'L':
				return k/1000;
			case 'V':
				return mu / 1000.0;
			case 'F':
				return Tfreeze+273.15;
			case 'M':
				return Tmax+273.15;
			case 'H':
				return u/1000+p/rho;
			case 'S':
				return s/1000;
			default:
				return _HUGE;
		}
	}

	else if (strncmp(Ref,"Methanol",8)==0)
	{
		switch (Ref[9])
		{
			case '1':
				C_gly=10; break;
			case '2':
				C_gly=20; break;
			case '3':
				C_gly=30; break;
			case '4':
				C_gly=40; break;
			default:
				return _HUGE;
		}

		/*
		Brine() takes inputs of temperature [C]  and secondary fluid concentration in water [mass %]
		Outputs are freezing temperature [C], density [kg/m^3], Specific Heat [J/kg-K], Conductivity [W/m-K], Viscosity [mPa-s]
		*/
		// Set the temperature to within the band to avoid zero
		if (Output=='M' || Output=='F'){ T=300;}
		Brine("Methanol", T - 273.15, C_gly, &Tfreeze, &Tmax, &rho, &cp, &k, &mu,&u,&s);
		switch(Output)
		{
			case 'D':
				return rho;
			case 'C':
				return cp/1000;
			case 'L':
				return k/1000;
			case 'V':
				return mu / 1000.0;
			case 'F':
				return Tfreeze+273.15;
			case 'M':
				return Tmax+273.15;
			case 'H':
				return u/1000+p/rho;
			case 'S':
				return s/1000;
			default:
				return _HUGE;
		}
	}

	else if (strncmp(Ref,"NH3/H2O",7)==0)
	{
		switch (Ref[8])
		{
			case '1':
				C_gly=10; break;
			case '2':
				C_gly=20; break;
			default:
				return _HUGE;
		}

		/*
		Brine() takes inputs of temperature [C]  and secondary fluid concentration in water [mass %]
		Outputs are freezing temperature [C], density [kg/m^3], Specific Heat [J/kg-K], Conductivity [W/m-K], Viscosity [mPa-s]
		*/
		// Set the temperature to within the band to avoid zero
		if (Output=='M' || Output=='F'){ T=300;}
		Brine("NH3", T - 273.15, C_gly, &Tfreeze, &Tmax, &rho, &cp, &k, &mu,&u,&s);
		switch(Output)
		{
			case 'D':
				return rho;
			case 'C':
				return cp/1000;
			case 'L':
				return k/1000;
			case 'V':
				return mu / 1000.0;
			case 'F':
				return Tfreeze+273.15;
			case 'M':
				return Tmax+273.15;
			case 'H':
				return u/1000+p/rho;
			case 'S':
				return s/1000;
			default:
				return _HUGE;
		}
	}
	else
	{
		return _HUGE;
	}
}
#if defined(_WIN32) || defined(__WIN32__) //Check if it is a windows machine, if not, hide this function
double REFPROP(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
	int j;
	long i,ierr=0;
	char hf[refpropcharlength*ncmax], hrf[lengthofreference+1],
	herr[errormessagelength+1],hfmix[refpropcharlength+1];
	
	double x[ncmax],xliq[ncmax],xvap[ncmax];
	char RefString[255];
	double T,p=0,d,dl,dv,q,e,h,s,cv,cp,w,MW,hl,hv,sl,sv,ul,uv,pl,pv,eta,tcx,Q,Tcrit,pcrit,dcrit,rho;

	// First create a pointer to an instance of the library
	// Then have windows load the library.
		HINSTANCE RefpropdllInstance;
		RefpropdllInstance = LoadLibrary("c:/Program Files/REFPROP/refprop.dll");

	// Then get pointers into the dll to the actual functions.
		ABFL1dll = (fp_ABFL1dllTYPE) GetProcAddress(RefpropdllInstance,"ABFL1dll");
		ABFL2dll = (fp_ABFL2dllTYPE) GetProcAddress(RefpropdllInstance,"ABFL2dll");
		ACTVYdll = (fp_ACTVYdllTYPE) GetProcAddress(RefpropdllInstance,"ACTVYdll");
		AGdll = (fp_AGdllTYPE) GetProcAddress(RefpropdllInstance,"AGdll");
		CCRITdll = (fp_CCRITdllTYPE) GetProcAddress(RefpropdllInstance,"CCRITdll");
		CP0dll = (fp_CP0dllTYPE) GetProcAddress(RefpropdllInstance,"CP0dll");
		CRITPdll = (fp_CRITPdllTYPE) GetProcAddress(RefpropdllInstance,"CRITPdll");
		CSATKdll = (fp_CSATKdllTYPE) GetProcAddress(RefpropdllInstance,"CSATKdll");
		CV2PKdll = (fp_CV2PKdllTYPE) GetProcAddress(RefpropdllInstance,"CV2PKdll");
		CVCPKdll = (fp_CVCPKdllTYPE) GetProcAddress(RefpropdllInstance,"CVCPKdll");
		CVCPdll = (fp_CVCPdllTYPE) GetProcAddress(RefpropdllInstance,"CVCPdll"); 
		DBDTdll = (fp_DBDTdllTYPE) GetProcAddress(RefpropdllInstance,"DBDTdll");
		DBFL1dll = (fp_DBFL1dllTYPE) GetProcAddress(RefpropdllInstance,"DBFL1dll");
		DBFL2dll = (fp_DBFL2dllTYPE) GetProcAddress(RefpropdllInstance,"DBFL2dll");
		DDDPdll = (fp_DDDPdllTYPE) GetProcAddress(RefpropdllInstance,"DDDPdll");
		DDDTdll = (fp_DDDTdllTYPE) GetProcAddress(RefpropdllInstance,"DDDTdll");
		DEFLSHdll = (fp_DEFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"DEFLSHdll");
		DHD1dll = (fp_DHD1dllTYPE) GetProcAddress(RefpropdllInstance,"DHD1dll");
		DHFLSHdll = (fp_DHFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"DHFLSHdll");
		DIELECdll = (fp_DIELECdllTYPE) GetProcAddress(RefpropdllInstance,"DIELECdll");
		DOTFILLdll = (fp_DOTFILLdllTYPE) GetProcAddress(RefpropdllInstance,"DOTFILLdll");
		DPDD2dll = (fp_DPDD2dllTYPE) GetProcAddress(RefpropdllInstance,"DPDD2dll");
		DPDDKdll = (fp_DPDDKdllTYPE) GetProcAddress(RefpropdllInstance,"DPDDKdll");
		DPDDdll = (fp_DPDDdllTYPE) GetProcAddress(RefpropdllInstance,"DPDDdll");
		DPDTKdll = (fp_DPDTKdllTYPE) GetProcAddress(RefpropdllInstance,"DPDTKdll");
		DPDTdll = (fp_DPDTdllTYPE) GetProcAddress(RefpropdllInstance,"DPDTdll");
		DPTSATKdll = (fp_DPTSATKdllTYPE) GetProcAddress(RefpropdllInstance,"DPTSATKdll");
		DSFLSHdll = (fp_DSFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"DSFLSHdll");
		ENTHALdll = (fp_ENTHALdllTYPE) GetProcAddress(RefpropdllInstance,"ENTHALdll"); //**
		ENTROdll = (fp_ENTROdllTYPE) GetProcAddress(RefpropdllInstance,"ENTROdll");
		ESFLSHdll = (fp_ESFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"ESFLSHdll");
		FGCTYdll = (fp_FGCTYdllTYPE) GetProcAddress(RefpropdllInstance,"FGCTYdll");
		FPVdll = (fp_FPVdllTYPE) GetProcAddress(RefpropdllInstance,"FPVdll");
		GERG04dll = (fp_GERG04dllTYPE) GetProcAddress(RefpropdllInstance,"GERG04dll");
		GETFIJdll = (fp_GETFIJdllTYPE) GetProcAddress(RefpropdllInstance,"GETFIJdll");
		GETKTVdll = (fp_GETKTVdllTYPE) GetProcAddress(RefpropdllInstance,"GETKTVdll");
		GIBBSdll = (fp_GIBBSdllTYPE) GetProcAddress(RefpropdllInstance,"GIBBSdll");
		HSFLSHdll = (fp_HSFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"HSFLSHdll");
		INFOdll = (fp_INFOdllTYPE) GetProcAddress(RefpropdllInstance,"INFOdll");
		LIMITKdll = (fp_LIMITKdllTYPE) GetProcAddress(RefpropdllInstance,"LIMITKdll");
		LIMITSdll = (fp_LIMITSdllTYPE) GetProcAddress(RefpropdllInstance,"LIMITSdll");
		LIMITXdll = (fp_LIMITXdllTYPE) GetProcAddress(RefpropdllInstance,"LIMITXdll");
		MELTPdll = (fp_MELTPdllTYPE) GetProcAddress(RefpropdllInstance,"MELTPdll");
		MELTTdll = (fp_MELTTdllTYPE) GetProcAddress(RefpropdllInstance,"MELTTdll");
		MLTH2Odll = (fp_MLTH2OdllTYPE) GetProcAddress(RefpropdllInstance,"MLTH2Odll");
		NAMEdll = (fp_NAMEdllTYPE) GetProcAddress(RefpropdllInstance,"NAMEdll");
		PDFL1dll = (fp_PDFL1dllTYPE) GetProcAddress(RefpropdllInstance,"PDFL1dll");
		PDFLSHdll = (fp_PDFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"PDFLSHdll");
		PEFLSHdll = (fp_PEFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"PEFLSHdll");
		PHFL1dll = (fp_PHFL1dllTYPE) GetProcAddress(RefpropdllInstance,"PHFL1dll");
		PHFLSHdll = (fp_PHFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"PHFLSHdll");
		PQFLSHdll = (fp_PQFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"PQFLSHdll");
		PREOSdll = (fp_PREOSdllTYPE) GetProcAddress(RefpropdllInstance,"PREOSdll");
		PRESSdll = (fp_PRESSdllTYPE) GetProcAddress(RefpropdllInstance,"PRESSdll");
		PSFL1dll = (fp_PSFL1dllTYPE) GetProcAddress(RefpropdllInstance,"PSFL1dll");
		PSFLSHdll = (fp_PSFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"PSFLSHdll");
		PUREFLDdll = (fp_PUREFLDdllTYPE) GetProcAddress(RefpropdllInstance,"PUREFLDdll");
		QMASSdll = (fp_QMASSdllTYPE) GetProcAddress(RefpropdllInstance,"QMASSdll");
		QMOLEdll = (fp_QMOLEdllTYPE) GetProcAddress(RefpropdllInstance,"QMOLEdll");
		SATDdll = (fp_SATDdllTYPE) GetProcAddress(RefpropdllInstance,"SATDdll");
		SATEdll = (fp_SATEdllTYPE) GetProcAddress(RefpropdllInstance,"SATEdll");
		SATHdll = (fp_SATHdllTYPE) GetProcAddress(RefpropdllInstance,"SATHdll");
		SATPdll = (fp_SATPdllTYPE) GetProcAddress(RefpropdllInstance,"SATPdll");
		SATSdll = (fp_SATSdllTYPE) GetProcAddress(RefpropdllInstance,"SATSdll");
		SATTdll = (fp_SATTdllTYPE) GetProcAddress(RefpropdllInstance,"SATTdll");
		SETAGAdll = (fp_SETAGAdllTYPE) GetProcAddress(RefpropdllInstance,"SETAGAdll");
		SETKTVdll = (fp_SETKTVdllTYPE) GetProcAddress(RefpropdllInstance,"SETKTVdll");
		SETMIXdll = (fp_SETMIXdllTYPE) GetProcAddress(RefpropdllInstance,"SETMIXdll");
		SETMODdll = (fp_SETMODdllTYPE) GetProcAddress(RefpropdllInstance,"SETMODdll");
		SETREFdll = (fp_SETREFdllTYPE) GetProcAddress(RefpropdllInstance,"SETREFdll");
		SETUPdll = (fp_SETUPdllTYPE) GetProcAddress(RefpropdllInstance,"SETUPdll");
		SPECGRdll = (fp_SPECGRdllTYPE) GetProcAddress(RefpropdllInstance,"SPECGRdll");
		SUBLPdll = (fp_SUBLPdllTYPE) GetProcAddress(RefpropdllInstance,"SUBLPdll");
		SUBLTdll = (fp_SUBLTdllTYPE) GetProcAddress(RefpropdllInstance,"SUBLTdll");
		SURFTdll = (fp_SURFTdllTYPE) GetProcAddress(RefpropdllInstance,"SURFTdll");
		SURTENdll = (fp_SURTENdllTYPE) GetProcAddress(RefpropdllInstance,"SURTENdll");
		TDFLSHdll = (fp_TDFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"TDFLSHdll");
		TEFLSHdll = (fp_TEFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"TEFLSHdll");
		THERM0dll = (fp_THERM0dllTYPE) GetProcAddress(RefpropdllInstance,"THERM0dll");
		THERM2dll = (fp_THERM2dllTYPE) GetProcAddress(RefpropdllInstance,"THERM2dll");
		THERM3dll = (fp_THERM3dllTYPE) GetProcAddress(RefpropdllInstance,"THERM3dll");
		THERMdll = (fp_THERMdllTYPE) GetProcAddress(RefpropdllInstance,"THERMdll");
		THFLSHdll = (fp_THFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"THFLSHdll");
		TPFLSHdll = (fp_TPFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"TPFLSHdll");
		TPRHOdll = (fp_TPRHOdllTYPE) GetProcAddress(RefpropdllInstance,"TPRHOdll");
		TQFLSHdll = (fp_TQFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"TQFLSHdll");
		TRNPRPdll = (fp_TRNPRPdllTYPE) GetProcAddress(RefpropdllInstance,"TRNPRPdll");
		TSFLSHdll = (fp_TSFLSHdllTYPE) GetProcAddress(RefpropdllInstance,"TSFLSHdll");
		VIRBdll = (fp_VIRBdllTYPE) GetProcAddress(RefpropdllInstance,"VIRBdll");
		VIRCdll = (fp_VIRCdllTYPE) GetProcAddress(RefpropdllInstance,"VIRCdll");
		WMOLdll = (fp_WMOLdllTYPE) GetProcAddress(RefpropdllInstance,"WMOLdll");
		XMASSdll = (fp_XMASSdllTYPE) GetProcAddress(RefpropdllInstance,"XMASSdll");
		XMOLEdll = (fp_XMOLEdllTYPE) GetProcAddress(RefpropdllInstance,"XMOLEdll");
		
		// If the fluid name starts with the string "REFPROP-", chop off the "REFPROP-"
		if (!strncmp(Ref,"REFPROP-",8))
		{
		char *REFPROPRef=NULL,*RefCopy=NULL;
		double prop;
			
		// Allocate space for refrigerant name
			RefCopy=malloc(strlen(Ref)+1);
		// Make a backup copy
			strcpy(RefCopy,Ref);
		// Chop off the "REFPROP-"
			REFPROPRef = strtok(RefCopy,"-");
			REFPROPRef = strtok(NULL,"-");
		// Run with the stripped Refrigerant name
			prop=REFPROP(Output,Name1,Prop1,Name2,Prop2,REFPROPRef);
		// Free allocated memory
			free(RefCopy);
		// Return the new value
			return prop;
		}
		
		if (!strncmp(Ref,"MIX",3))
		{
			// Sample is "REFPROP-MIX:R32[0.697615]&R125[0.302385]"
			char *REFPROPRef=NULL,*RefCopy=NULL,RefString[255],*Refs[20],*Refrigerant;
			double molefraction;

			//Set global fluid type flag
			FluidType=FLUIDTYPE_REFPROP;
			// Allocate space for refrigerant name
			RefCopy=malloc(strlen(Ref)+1);
			// Make a backup copy
			strcpy(RefCopy,Ref);
			// Chop off the "MIX"
			REFPROPRef = strtok(RefCopy,":");
			i=1;
			while (REFPROPRef!=NULL)
			{
				Refs[i-1]=strtok(NULL,"&");
				if (Refs[i-1]==NULL)
				{
					i--;
					break;
				}
				else
					i++;
			}
			//Flush out RefString
			sprintf(RefString,"");
			for (j=0;j<i;j++)
			{	
				//Get component and its mole fraction
				Refrigerant=strtok(Refs[j],"[]");
				molefraction=strtod(strtok(NULL,"[]"),NULL);
				x[j]=molefraction;
				if (j==0)
					sprintf(RefString,"%s%s.fld",RefString,Refs[j]);
				else
					sprintf(RefString,"%s|%s.fld",RefString,Refs[j]);
			}
			// Free allocated memory
			free(RefCopy);
		}
		else if (!strcmp(Ref,"R508B"))
		{
			i=2;
			strcpy(RefString,"R23.fld|R116.fld");
			x[0]=0.62675;
			x[1]=0.37325;
		}
		else if (!strcmp(Ref,"R410A"))
		{
			i=2;
			strcpy(RefString,"R32.fld|R125.fld");
			x[0]=0.697615;
			x[1]=0.302385;
		}
		else if (!strcmp(Ref,"R404A"))
		{
			i=3;
			strcpy(RefString,"R125.fld|R134a.fld|R143a.fld");
			x[0]=0.35782;
			x[1]=0.038264;
			x[2]=0.60392;
		}
        else if (!strcmp(Ref,"Air"))
		{
			i=3;
			strcpy(RefString,"Nitrogen.fld|Oxygen.fld|Argon.fld");
			x[0]=0.7812;
			x[1]=0.2096;
			x[2]=0.0092;
		}
		else
		{
			i=1;
			strcpy(RefString,"");
			strcat(RefString,Ref);
			strcat(RefString,".fld");
			x[0]=1.0;     //Pure fluid
		}

		strcpy(hf,RefString);
		strcpy(hfmix,"hmx.bnc");
		strcpy(hrf,"DEF");
		strcpy(herr,"Ok");
		
		// If the name of the refrigerant doesn't match 
		// that of the currently loaded refrigerant
		if (strcmp(LoadedREFPROPRef,Ref))
		{
			//...Call SETUP to initialize the program
			SETUPdll(&i, hf, hfmix, hrf, &ierr, herr,
				refpropcharlength*ncmax,refpropcharlength,
				lengthofreference,errormessagelength);
			if (ierr != 0) printf("REFPROP setup gives this error during SETUP: %s\n",herr);
			//Copy the name of the loaded refrigerant back into the temporary holder
			strcpy(LoadedREFPROPRef,Ref);
		}

		// Get the molar mass of the fluid
		WMOLdll(x,&MW);
		if (Output=='B')
		{
			// Critical temperature
			CRITPdll(x,&Tcrit,&pcrit,&dcrit,&ierr,herr,255);
			return Tcrit;
		}
		else if (Output=='E')
		{
			// Critical pressure
			CRITPdll(x,&Tcrit,&pcrit,&dcrit,&ierr,herr,255);
			return pcrit;
		}
		else if (Output=='R')
		{
			long icomp;
			double wmm,Ttriple,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas;
			// Triple point temperature
			icomp=1;
			if (i>1)
			{
				fprintf(stderr,"Error: Triple point temperature only defined for pure fluids\n");
				return 200;
			}
			INFOdll(&icomp,&wmm,&Ttriple,&tnbpt,&tc,&pc,&Dc,&Zc,&acf,&dip,&Rgas);
			return Ttriple;
		}
		else if (Output=='M')
		{
			// mole mass
			return MW;
		}
		else if (Name1=='T' && Name2=='P')
		{
			// T in K, P in kPa

			// Use flash routine to find properties
			T=Prop1;
			p=Prop2;  
			TPFLSHdll(&T,&p,x,&d,&dl,&dv,xliq,xvap,&q,&e,&h,&s,&cv,&cp,&w,&ierr,herr,errormessagelength);
			if (Output=='H') return h/MW;
			else if (Output=='D') return d*MW;
			else if (Output=='S') return s/MW;
			else if (Output=='U') return e/MW;
			else if (Output=='C') return cp/MW;
			else if (Output=='O') return cv/MW;
			else if (Output=='V') 
			{
				TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return eta/1.0e6; //uPa-s to Pa-s
			} 
			else if (Output=='L')
			{
				TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return tcx/1000.0; //W/m-K to kW/m-K
			}
			else
				return _HUGE;
		}
		else if (Name1=='T' && Name2=='D')
		{
			// T in K, D in kg/m^3
			// This is the explicit formulation of the EOS
			T=Prop1;
			rho=Prop2/MW;
			
			TDFLSHdll(&T,&rho,x,&p,&dl,&dv,xliq,xvap,&q,&e,&h,&s,&cv,&cp,&w,&ierr,herr,errormessagelength);

			if (Output=='P')
			{
				return p;
			}
			if (Output=='H')
			{
				return h/MW;
			}
			else if (Output=='S')
			{
				return s/MW;
			}
			else if (Output=='U')
			{
				return (h-p/rho)/MW;
			}
			else if (Output=='C')
			{
				return cp/MW;
			}
			else if (Output=='O')
			{
				return cv/MW;
			}
			else if (Output=='V') 
			{
				TRNPRPdll(&T,&rho,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return eta/1.0e6; //uPa-s to Pa-s
			} 
			else if (Output=='L')
			{
				TRNPRPdll(&T,&rho,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return tcx/1000.0; //W/m-K to kW/m-K
			}
			else
				return _HUGE;
		}
		else if (Name1=='T' && Name2=='Q')
		{
			T=Prop1;
			Q=Prop2;
			// Saturation Density
			SATTdll(&T,x,&i,&p,&dl,&dv,xliq,xvap,&ierr,herr,errormessagelength);
			if (Output=='D') 
			{
				return 1/(Q/dv+(1-Q)/dl)*MW;
			}
			else if (Output=='P') 
			{
				PRESSdll(&T,&dl,xliq,&pl);
				PRESSdll(&T,&dv,xvap,&pv);
				return (pv*Q+pl*(1-Q));
			}
			else if (Output=='H') 
			{
				ENTHALdll(&T,&dl,xliq,&hl);
				ENTHALdll(&T,&dv,xvap,&hv);
				return (hv*Q+hl*(1-Q))/MW; // J/kg to kJ/kg
			}
			else if (Output=='S') 
			{
				ENTROdll(&T,&dl,xliq,&sl);
				ENTROdll(&T,&dv,xvap,&sv);
				return (sv*Q+sl*(1-Q))/MW; // J/kg-K to kJ/kg-K
			}
			else if (Output=='U') 
			{
				ENTHALdll(&T,&dl,xliq,&hl);
				ENTHALdll(&T,&dv,xvap,&hv);
				ul=hl-p/dl;
				uv=hv-p/dv;
				return (uv*Q+ul*(1-Q))/MW; // J/kg to kJ/kg
			}
			else if (Output=='C') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				CVCPdll(&T,&d,x,&cv,&cp);
				return cp/MW; // J/kg-K to kJ/kg-K
			}
			else if (Output=='O') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				CVCPdll(&T,&d,x,&cv,&cp);
				return cv/MW; // J/kg-K to kJ/kg-K
			}
			else if (Output=='V') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return eta/1.0e6; //uPa-s to Pa-s
			}
			else if (Output=='L') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return tcx/1000.0; //W/m-K to kW/m-K
			}
			else
				return _HUGE;
		}
		else if (Name1=='P' && Name2=='Q')
		{
			p=Prop1;
			Q=Prop2;
			// Saturation Density
			SATPdll(&p,x,&i,&T,&dl,&dv,xliq,xvap,&ierr,herr,errormessagelength);
			if (Output=='T')
			{
				return T;
			}
			else if (Output=='D') 
			{
				return 1/(Q/dv+(1-Q)/dl)*MW;
			}
			else if (Output=='P') 
			{
				PRESSdll(&T,&dl,xliq,&pl);
				PRESSdll(&T,&dv,xvap,&pv);
				return (pv*Q+pl*(1-Q));
			}
			else if (Output=='H') 
			{
				ENTHALdll(&T,&dl,xliq,&hl);
				ENTHALdll(&T,&dv,xvap,&hv);
				return (hv*Q+hl*(1-Q))/MW; // J/kg to kJ/kg
			}
			else if (Output=='S') 
			{
				ENTROdll(&T,&dl,xliq,&sl);
				ENTROdll(&T,&dv,xvap,&sv);
				return (sv*Q+sl*(1-Q))/MW; // J/kg-K to kJ/kg-K
			}
			else if (Output=='U') 
			{
				ENTHALdll(&T,&dl,xliq,&hl);
				ENTHALdll(&T,&dv,xvap,&hv);
				ul=hl-p/dl;
				uv=hv-p/dv;
				return (uv*Q+ul*(1-Q))/MW; // J/kg to kJ/kg
			}
			else if (Output=='C') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				CVCPdll(&T,&d,x,&cv,&cp);
				return cp/MW; // J/kg-K to kJ/kg-K
			}
			else if (Output=='O') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				CVCPdll(&T,&d,x,&cv,&cp);
				return cv/MW; // J/kg-K to kJ/kg-K
			}
			else if (Output=='V') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return eta/1.0e6; //uPa-s to Pa-s
			}
			else if (Output=='L') 
			{
				d=1/(Q/dv+(1-Q)/dl);
				TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
				return tcx/1000.0; //W/m-K to kW/m-K
			}
			else
				return _HUGE;
		}
		else
			return _HUGE;
}
#endif

void Help()
{
	printf("CoolProp Help\n");
	printf("CoolProp is written by Ian Bell (ihb2@cornell.edu)\n");
	printf("\n");
	printf("Following the naming conventions of MATLAB linked with REFPROP,\n");
	printf("each output property is represented by one character:\n");
	printf("\n");
	printf("P   Pressure [kPa]\n");
	printf("T   Temperature [K]\n");
	printf("D   Density [kg/m3]\n");
	printf("H   Enthalpy [kJ/kg]\n");
	printf("S   Entropy [kJ/(kg/K)]\n");
	printf("U   Internal energy [kJ/kg]\n");
	printf("C   Cp [kJ/(kg K)]\n");
	printf("O   Cv [kJ/(kg K)]\n");
	printf("K   Ratio of specific heats (Cp/Cv) [-]\n");
	printf("A   Speed of sound [m/s]\n");
	printf("X   liquid phase and gas phase composition (mass fractions)\n");
	printf("V   Dynamic viscosity [Pa*s]\n");
	printf("L   Thermal conductivity [kW/(m K)]\n");
	printf("Q   Quality (vapor fraction) (kg/kg)\n");
	printf("I   Surface tension [N/m]\n");
	printf("F   Freezing point of secondary fluid [K] **NOT IN MATLAB-REFPROP **\n");
	printf("M   Maximum temperature for secondary fluid [K] **NOT IN MATLAB-REFPROP **\n");
	printf("M   Molar mass for non-secondary fluid [g/mol] **NOT IN MATLAB-REFPROP **\n");
	printf("B   Critical Temperature [K] **NOT IN MATLAB-REFPROP **\n");
	printf("E   Critical Pressure [kPa] **NOT IN MATLAB-REFPROP **\n");
	printf("R   Triple point temperature [K] **NOT IN MATLAB-REFPROP **\n");
	printf("\n");
	printf("******** To call **************\n");
	printf("To call the function Props, for instance for R410A at 300K, 400 kPa, you would do:\n");
	printf("Props(\"H\",\"T\",300,\"P\",400,\"R410A\")\n");
	printf("\n");
	printf("Or to call a pure fluid from REFPROP (for instance Propane).  \n");
	printf("The name of the refrigerant is \"REPFROP-\" plus the REFPROP defined name of the fluid, for instance\n");
	printf("\"Propane\" for propane (R290)\n");
	printf("\n");
	printf("See the folder C:\\Program Files\\REFPROP\\fluids for the names of the fluids\n");
	printf("\n");
	printf("To call Propane from REFPROP:\n");
	printf("Props(\"H\",\"T\",300,\"P\",400,\"REFPROP-Propane\")\n");
	printf("\n");
	printf("**************** Inputs ***************\n");
	printf("The limited list of inputs that are allowed are:\n");
	printf("\n");
	printf("Prop1    ||    Prop2\n");
	printf("--------------------\n");
	printf("  T      ||      P\n");
	printf("  T      ||      Q\n");
	printf("  T      ||      D\n");

}

int IsFluidType(char *Ref, char *Type)
{
	// Call a function to set the fluid-specific properties
	Props('M','T',0,'P',0,Ref);
	
	if (((FluidType==FLUIDTYPE_REFRIGERANT_PURE) || (FluidType==FLUIDTYPE_REFPROP)) && !strcmp(Type,"PureFluid"))
	{
		return 1;
	}
	else if (FluidType==FLUIDTYPE_BRINE && !strcmp(Type,"Brine"))
	{
		return 1;
	}
	else if (FluidType==FLUIDTYPE_REFRIGERANT_PSEUDOPURE && !strcmp(Type,"PseudoPure"))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
int Phase(double T, double rho, char * Ref)
{
	double rhosatL,rhosatV,p,Tbubble,Tdew;
	if (fabs(SecFluids('D',0,0,Ref))<1e10)
	{
		//It's a secondary fluid, always subcooled
		return PHASE_SUBCOOLED;
	}
	p=Props('P','T',T,'D',rho,Ref);
	if (p>pcrit(Ref))
	{
		return PHASE_SUPERCRITICAL;
	}
	else
	{
		Tbubble=Tsat(Ref, p, 0.0,T);
		Tdew=Tsat(Ref, p, 1.0,T);
		rhosatV=Props('D','T',Tdew,'Q',1,Ref);
		rhosatL=Props('D','T',Tbubble,'Q',0,Ref);

		if (rho<rhosatV)
			return PHASE_SUPERHEATED;
		else if (rho>rhosatL)
			return PHASE_SUBCOOLED;
		else
			return PHASE_TWOPHASE;
	}
}

double Props(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
	double T,p,Q,rhoV,rhoL,Value,rho,pdp,pbp;
	int errCode,isTwoPhase;
	char errString[ERRSTRLENGTH];
	
	/*
	Following the naming conventions of MATLAB linked with REFPROP,
	each output property is represented by one character:

	P   Pressure [kPa]
	T   Temperature [K]
	D   Density [kg/m3]
	H   Enthalpy [kJ/kg]
	S   Entropy [kJ/(kg/K)]
	U   Internal energy [kJ/kg]
	C   Cp [kJ/(kg K)]
	O   Cv [kJ/(kg K)]
	K   Ratio of specific heats (Cp/Cv) [-]
	A   Speed of sound [m/s]
	X   liquid phase and gas phase composition (mass fractions)
	V   Dynamic viscosity [Pa*s]
	L   Thermal conductivity [kW/(m K)]
	Q   Quality (vapor fraction) (kg/kg)
	I   Surface tension [N/m]
	F	Freezing point of secondary fluid [K] **NOT IN MATLAB-REFPROP **
	M	Maximum temperature for secondary fluid [K] **NOT IN MATLAB-REFPROP **
	B	Critical Temperature [K] **NOT IN MATLAB-REFPROP **
	E	Critical Pressure [K] **NOT IN MATLAB-REFPROP **
	R   Triple point temperature [K]
	*/

	errCode=0;
	
	


	// **********************************************************************************
	// **********************************************************************************
	//                                   REFPROP
	// **********************************************************************************
	// **********************************************************************************
	
	/* 
	If the fluid name is not actually a refrigerant name, but a string beginning with "REFPROP-",
	then REFPROP is used to calculate the desired property.
	*/
	#if defined(_WIN32) || defined(__WIN32__)
	if (strncmp(Ref,"REFPROP-",8)==0)  // First eight characters match "REFPROP-"
	{
	#else
	if (0) // Automatically skip it because REFPROP is not supported on this platform
	{
	#endif
		FluidType=FLUIDTYPE_REFPROP;
		return REFPROP(Output,Name1,Prop1,Name2,Prop2,Ref);
	}

	// **********************************************************************************
	// **********************************************************************************
	//                                Normal Property evaluation
	// **********************************************************************************
	// **********************************************************************************

	// It's a secondary fluid
	else if (!strcmp(Ref,"HC-10") || (Ref[0]=='E' && Ref[1]=='G') || 
		(Ref[0]=='P' && Ref[1]=='G') || strncmp(Ref,"Methanol",8)==0 || strncmp(Ref,"NH3/H2O",7)==0)
	{
		if (Name1!='T' || Name2!='P')
		{
			printf("Warning: For brine, Name1 must be 'T' and Name2 must be 'P'\n");
		}
		FluidType=FLUIDTYPE_BRINE;
		return SecFluids(Output,Prop1,Prop2,Ref);
	}
	else 
	{
		T=Prop1; 
		
		// Wire up the function pointers for the given refrigerant
		if (!strcmp(Ref,"Argon"))
		{
			FluidType=FLUIDTYPE_REFRIGERANT_PURE;
			p_func=p_Argon;
			h_func=h_Argon;
			s_func=s_Argon;
			u_func=u_Argon;
			rho_func=rho_Argon;
			cp_func=cp_Argon;
			cv_func=cv_Argon;
			visc_func=visc_Argon;
			k_func=k_Argon;
			w_func=w_Argon;
			Ttriple_func=Ttriple_Argon;
			Tcrit_func=Tcrit_Argon;
			pcrit_func=pcrit_Argon;
			rhocrit_func=rhocrit_Argon;
			MM_func=MM_Argon;
			rhosatV_func=rhosatV_Argon;
			rhosatL_func=rhosatL_Argon;
			psat_func=psat_Argon;

			phir_func=phir_Argon;
			dphir_dDelta_func=dphir_dDelta_Argon;
			dphir2_dDelta2_func=dphir2_dDelta2_Argon;
			dphir2_dDelta_dTau_func=dphir2_dDelta_dTau_Argon;
			dphir_dTau_func=dphir_dTau_Argon;
			dphir2_dTau2_func=dphir2_dTau2_Argon;
			phi0_func=phi0_Argon;
			dphi0_dDelta_func=dphi0_dDelta_Argon;
			dphi02_dDelta2_func=dphi02_dDelta2_Argon;
			dphi0_dTau_func=dphi0_dTau_Argon;
			dphi02_dTau2_func=dphi02_dTau2_Argon;

		}
		else if (!strcmp(Ref,"Nitrogen") || !strcmp(Ref,"N2"))
		{
			FluidType=FLUIDTYPE_REFRIGERANT_PURE;
			p_func=p_Nitrogen;
			h_func=h_Nitrogen;
			s_func=s_Nitrogen;
			u_func=u_Nitrogen;
			rho_func=rho_Nitrogen;
			cp_func=cp_Nitrogen;
			cv_func=cv_Nitrogen;
			visc_func=visc_Nitrogen;
			k_func=k_Nitrogen;
			w_func=w_Nitrogen;
			Ttriple_func=Ttriple_Nitrogen;
			Tcrit_func=Tcrit_Nitrogen;
			pcrit_func=pcrit_Nitrogen;
			rhocrit_func=rhocrit_Nitrogen;
			MM_func=MM_Nitrogen;
			rhosatV_func=rhosatV_Nitrogen;
			rhosatL_func=rhosatL_Nitrogen;
			psat_func=psat_Nitrogen;

			phir_func=phir_Nitrogen;
			dphir_dDelta_func=dphir_dDelta_Nitrogen;
			dphir2_dDelta2_func=dphir2_dDelta2_Nitrogen;
			dphir2_dDelta_dTau_func=dphir2_dDelta_dTau_Nitrogen;
			dphir_dTau_func=dphir_dTau_Nitrogen;
			dphir2_dTau2_func=dphir2_dTau2_Nitrogen;
			phi0_func=phi0_Nitrogen;
			dphi0_dDelta_func=dphi0_dDelta_Nitrogen;
			dphi02_dDelta2_func=dphi02_dDelta2_Nitrogen;
			dphi0_dTau_func=dphi0_dTau_Nitrogen;
			dphi02_dTau2_func=dphi02_dTau2_Nitrogen;
		}
		else if (!strcmp(Ref,"R744") || !strcmp(Ref,"CO2"))
		{
			FluidType=FLUIDTYPE_REFRIGERANT_PURE;
			p_func=p_R744;
			h_func=h_R744;
			s_func=s_R744;
			u_func=u_R744;
			rho_func=rho_R744;
			cp_func=cp_R744;
			cv_func=cv_R744;
			visc_func=visc_R744;
			k_func=k_R744;
			w_func=w_R744;
			Ttriple_func=Ttriple_R744;
			Tcrit_func=Tcrit_R744;
			pcrit_func=pcrit_R744;
			rhocrit_func=rhocrit_R744;
			MM_func=MM_R744;
			rhosatV_func=rhosatV_R744;
			rhosatL_func=rhosatL_R744;
			psat_func=psat_R744;

			phir_func=phir_R744;
			dphir_dDelta_func=dphir_dDelta_R744;
			dphir2_dDelta2_func=dphir2_dDelta2_R744;
			dphir2_dDelta_dTau_func=dphir2_dDelta_dTau_R744;
			dphir_dTau_func=dphir_dTau_R744;
			dphir2_dTau2_func=dphir2_dTau2_R744;
			phi0_func=phi0_R744;
			dphi0_dDelta_func=dphi0_dDelta_R744;
			dphi02_dDelta2_func=dphi02_dDelta2_R744;
			dphi0_dTau_func=dphi0_dTau_R744;
			dphi02_dTau2_func=dphi02_dTau2_R744;
		}
		else if (!strcmp(Ref,"R718") || !strcmp(Ref,"Water") || !strcmp(Ref,"H2O"))
		{
			FluidType=FLUIDTYPE_REFRIGERANT_PURE;
			p_func=p_Water;
			h_func=h_Water;
			s_func=s_Water;
			u_func=u_Water;
			rho_func=rho_Water;
			cp_func=cp_Water;
			cv_func=cv_Water;
			visc_func=visc_Water;
			k_func=k_Water;
			w_func=w_Water;
			Ttriple_func=Ttriple_Water;
			Tcrit_func=Tcrit_Water;
			pcrit_func=pcrit_Water;
			rhocrit_func=rhocrit_Water;
			MM_func=MM_Water;
			rhosatV_func=rhosatV_Water;
			rhosatL_func=rhosatL_Water;
			psat_func=psat_Water;

			phir_func=phir_Water;
			dphir_dDelta_func=dphir_dDelta_Water;
			dphir2_dDelta2_func=dphir2_dDelta2_Water;
			dphir2_dDelta_dTau_func=dphir2_dDelta_dTau_Water;
			dphir_dTau_func=dphir_dTau_Water;
			dphir2_dTau2_func=dphir2_dTau2_Water;
			phi0_func=phi0_Water;
			dphi0_dDelta_func=dphi0_dDelta_Water;
			dphi02_dDelta2_func=dphi02_dDelta2_Water;
			dphi0_dTau_func=dphi0_dTau_Water;
			dphi02_dTau2_func=dphi02_dTau2_Water;
		}
		else if (!strcmp(Ref,"R134a"))
		{
			FluidType=FLUIDTYPE_REFRIGERANT_PURE;
			p_func=p_R134a;
			h_func=h_R134a;
			s_func=s_R134a;
			u_func=u_R134a;
			rho_func=rho_R134a;
			cp_func=cp_R134a;
			cv_func=cv_R134a;
			visc_func=visc_R134a;
			k_func=k_R134a;
			w_func=w_R134a;
			Ttriple_func=Ttriple_R134a;
			Tcrit_func=Tcrit_R134a;
			pcrit_func=pcrit_R134a;
			rhocrit_func=rhocrit_R134a;
			MM_func=MM_R134a;
			rhosatV_func=rhosatV_R134a;
			rhosatL_func=rhosatL_R134a;
			psat_func=psat_R134a;
			
			phir_func=phir_R134a;
			dphir_dDelta_func=dphir_dDelta_R134a;
			dphir2_dDelta2_func=dphir2_dDelta2_R134a;
			dphir2_dDelta_dTau_func=dphir2_dDelta_dTau_R134a;
			dphir_dTau_func=dphir_dTau_R134a;
			dphir2_dTau2_func=dphir2_dTau2_R134a;
			phi0_func=phi0_R134a;
			dphi0_dDelta_func=dphi0_dDelta_R134a;
			dphi02_dDelta2_func=dphi02_dDelta2_R134a;
			dphi0_dTau_func=dphi0_dTau_R134a;
			dphi02_dTau2_func=dphi02_dTau2_R134a;
		}

		else if (!strcmp(Ref,"R290"))
		{
			FluidType=FLUIDTYPE_REFRIGERANT_PURE;
			p_func=p_R290;
			h_func=h_R290;
			s_func=s_R290;
			u_func=u_R290;
			rho_func=rho_R290;
			cp_func=cp_R290;
			cv_func=cv_R290;
			visc_func=visc_R290;
			k_func=k_R290;
			w_func=w_R290;
			Ttriple_func=Ttriple_R290;
			Tcrit_func=Tcrit_R290;
			pcrit_func=pcrit_R290;
			rhocrit_func=rhocrit_R290;
			MM_func=MM_R290;
			rhosatV_func=rhosatV_R290;
			rhosatL_func=rhosatL_R290;
			psat_func=psat_R290;

			phir_func=phir_R290;
			dphir_dDelta_func=dphir_dDelta_R290;
			dphir2_dDelta2_func=dphir2_dDelta2_R290;
			dphir2_dDelta_dTau_func=dphir2_dDelta_dTau_R290;
			dphir_dTau_func=dphir_dTau_R290;
			dphir2_dTau2_func=dphir2_dTau2_R290;
			phi0_func=phi0_R290;
			dphi0_dDelta_func=dphi0_dDelta_R290;
			dphi02_dDelta2_func=dphi02_dDelta2_R290;
			dphi0_dTau_func=dphi0_dTau_R290;
			dphi02_dTau2_func=dphi02_dTau2_R290;
		}
		else if (!strcmp(Ref,"R717") || !strcmp(Ref,"NH3") || !strcmp(Ref,"Ammonia") || !strcmp(Ref,"ammonia"))
		{
			FluidType=FLUIDTYPE_REFRIGERANT_PURE;
			p_func=p_R717;
			h_func=h_R717;
			s_func=s_R717;
			u_func=u_R717;
			rho_func=rho_R717;
			cp_func=cp_R717;
			cv_func=cv_R717;
			visc_func=visc_R717;
			k_func=k_R717;
			w_func=w_R717;
			Ttriple_func=Ttriple_R717;
			Tcrit_func=Tcrit_R717;
			pcrit_func=pcrit_R717;
			rhocrit_func=rhocrit_R717;
			MM_func=MM_R717;
			rhosatV_func=rhosatV_R717;
			rhosatL_func=rhosatL_R717;
			psat_func=psat_R717;

			phir_func=phir_R717;
			dphir_dDelta_func=dphir_dDelta_R717;
			dphir2_dDelta2_func=dphir2_dDelta2_R717;
			dphir2_dDelta_dTau_func=dphir2_dDelta_dTau_R717;
			dphir_dTau_func=dphir_dTau_R717;
			dphir2_dTau2_func=dphir2_dTau2_R717;
			phi0_func=phi0_R717;
			dphi0_dDelta_func=dphi0_dDelta_R717;
			dphi02_dDelta2_func=dphi02_dDelta2_R717;
			dphi0_dTau_func=dphi0_dTau_R717;
			dphi02_dTau2_func=dphi02_dTau2_R717;
		}
		else if (!strcmp(Ref,"R32"))
		{
			FluidType=FLUIDTYPE_REFRIGERANT_PURE;
			p_func=p_R32;
			h_func=h_R32;
			s_func=s_R32;
			u_func=u_R32;
			rho_func=rho_R32;
			cp_func=cp_R32;
			cv_func=cv_R32;
			visc_func=visc_R32;
			k_func=k_R32;
			w_func=w_R32;
			Ttriple_func=Ttriple_R32;
			Tcrit_func=Tcrit_R32;
			pcrit_func=pcrit_R32;
			MM_func=MM_R32;
			rhosatV_func=rhosatV_R32;
			rhosatL_func=rhosatL_R32;
			psat_func=psat_R32;
		}
        else if (!strcmp(Ref,"Air"))
		{
			FluidType=FLUIDTYPE_REFRIGERANT_PSEUDOPURE;
			p_func=p_Air;
			h_func=h_Air;
			s_func=s_Air;
			u_func=u_Air;
			rho_func=rho_Air;
			cp_func=cp_Air;
			cv_func=cv_Air;
			visc_func=visc_Air;
			k_func=k_Air;
			w_func=w_Air;
			Ttriple_func=Ttriple_Air;
			Tcrit_func=Tcrit_Air;
			pcrit_func=pcrit_Air;
			MM_func=MM_Air;
			rhosatV_func=rhosatV_Air;
			rhosatL_func=rhosatL_Air;
			pdp_func=pdp_Air;
			pbp_func=pbp_Air;
		}
		else if (!strcmp(Ref,"R410A"))
		{
			FluidType=FLUIDTYPE_REFRIGERANT_PSEUDOPURE;
			p_func=p_R410A;
			h_func=h_R410A;
			s_func=s_R410A;
			u_func=u_R410A;
			rho_func=rho_R410A;
			cp_func=cp_R410A;
			cv_func=cv_R410A;
			visc_func=visc_R410A;
			k_func=k_R410A;
			w_func=w_R410A;
			Ttriple_func=Ttriple_R410A;
			Tcrit_func=Tcrit_R410A;
			pcrit_func=pcrit_R410A;
			MM_func=MM_R410A;
			rhosatV_func=rhosatV_R410A;
			rhosatL_func=rhosatL_R410A;
			pdp_func=p_dp_R410A;
			pbp_func=p_bp_R410A;
		}
		else if (!strcmp(Ref,"R404A"))
		{
			FluidType=FLUIDTYPE_REFRIGERANT_PSEUDOPURE;
			p_func=p_R404A;
			h_func=h_R404A;
			s_func=s_R404A;
			u_func=u_R404A;
			rho_func=rho_R404A;
			cp_func=cp_R404A;
			cv_func=cv_R404A;
			visc_func=visc_R404A;
			k_func=k_R404A;
			w_func=w_R404A;
			Ttriple_func=Ttriple_R404A;
			Tcrit_func=Tcrit_R404A;
			pcrit_func=pcrit_R404A;
			MM_func=MM_R404A;
			rhosatV_func=rhosatV_R404A;
			rhosatL_func=rhosatL_R404A;
			pdp_func=p_dp_R404A;
			pbp_func=p_bp_R404A;
		}
		else if (!strcmp(Ref,"R407C"))
		{
			FluidType=FLUIDTYPE_REFRIGERANT_PSEUDOPURE;
			p_func=p_R407C;
			h_func=h_R407C;
			s_func=s_R407C;
			u_func=u_R407C;
			rho_func=rho_R407C;
			cp_func=cp_R407C;
			cv_func=cv_R407C;
			visc_func=visc_R407C;
			k_func=k_R407C;
			w_func=w_R407C;
			Ttriple_func=Ttriple_R407C;
			Tcrit_func=Tcrit_R407C;
			pcrit_func=pcrit_R407C;
			MM_func=MM_R407C;
			rhosatV_func=rhosatV_R407C;
			rhosatL_func=rhosatL_R407C;
			pdp_func=p_dp_R407C;
			pbp_func=p_bp_R407C;
		}
		else if (!strcmp(Ref,"R507A"))
		{
			FluidType=FLUIDTYPE_REFRIGERANT_PSEUDOPURE;
			p_func=p_R507A;
			h_func=h_R507A;
			s_func=s_R507A;
			u_func=u_R507A;
			rho_func=rho_R507A;
			cp_func=cp_R507A;
			cv_func=cv_R507A;
			visc_func=visc_R507A;
			k_func=k_R507A;
			w_func=w_R507A;
			Ttriple_func=Ttriple_R507A;
			Tcrit_func=Tcrit_R507A;
			pcrit_func=pcrit_R507A;
			MM_func=MM_R507A;
			rhosatV_func=rhosatV_R507A;
			rhosatL_func=rhosatL_R507A;
			pdp_func=p_dp_R507A;
			pbp_func=p_bp_R507A;
		}
		else
		{
			fprintf(stderr,"Refrigerant %s not allowed\n",Ref);
			return _HUGE;
		}

		if (Name1=='T' && (Name2=='P' || Name2=='D' || Name2=='Q'))
		{
			// Temperature and Pressure are the inputs
			if (Name2=='P')
			{
				p=Prop2;
				switch (Output)
				{
					case 'H':
						Value=h_func(T,p,TYPE_TPNoLookup); break;
					case 'D':
						Value=rho_func(T,p,TYPE_TPNoLookup); break;
					case 'S':
						Value=s_func(T,p,TYPE_TPNoLookup); break;
					case 'U':
						Value=u_func(T,p,TYPE_TPNoLookup); break;
					case 'C':
						Value=cp_func(T,p,TYPE_TPNoLookup); break;
					case 'O':
						Value=cv_func(T,p,TYPE_TPNoLookup); break;
					case 'V':
						Value=visc_func(T,p,TYPE_TPNoLookup); break;
					case 'L':
						Value=k_func(T,p,TYPE_TPNoLookup); break;
					case 'A':
						Value=w_func(T,p,TYPE_TPNoLookup); break;
					case 'M':
						Value=MM_func(); break;
					case 'E':
						Value=pcrit_func(); break;
					case 'B':
						Value=Tcrit_func(); break;
					case 'R':
						Value=Ttriple_func(); break;
					default:
						strcpy(errString,"Invalid Output Name");
						return -100;
				}
				return Value;
			}
			// Temperature and density are the inputs
			else if (Name2=='D')
			{
				rho=Prop2;
				
				// Check if it is one of the outputs that do not require any state information
				switch (Output)
				{
					case 'M':
						return MM_func(); 
					case 'E':
						return pcrit_func();
					case 'B':
						return Tcrit_func();
					case 'R':
						return Ttriple_func();
				}
				
				// Determine if it is somewhere in the two-phase region.

				// First check if the condition is well away from saturation 
				// by using the ancillary equations because this is much faster
				if (rho>1.02*rhosatL_func(T) || rho<0.98*rhosatV_func(T))
				{
					isTwoPhase=0;
				}
				else
				{
					if (!strcmp(Ref,"R404A") || !strcmp(Ref,"R410A") || !strcmp(Ref,"R407C") || !strcmp(Ref,"R507A"))
					{
						rhoV=rhosatV_func(T);
						rhoL=rhosatL_func(T);
					}
					else
					{
						rhosatPure(Ref,T,&rhoL,&rhoV,&p);
					}
					if (rho<=rhoL && rho>=rhoV)
						isTwoPhase=1;
					else
						isTwoPhase=0;
				}
				
				if (isTwoPhase==0)
				{
					// It is not two-phase, and use EOS
					switch (Output)
					{
						case 'P':
							Value=p_func(T,rho); break;
						case 'H':
							Value=h_func(T,rho,TYPE_Trho); break;
						case 'S':
							Value=s_func(T,rho,TYPE_Trho); break;
						case 'U':
							Value=u_func(T,rho,TYPE_Trho); break;
						case 'C':
							Value=cp_func(T,rho,TYPE_Trho); break;
						case 'O':
							Value=cv_func(T,rho,TYPE_Trho); break;
						case 'V':
							Value=visc_func(T,rho,TYPE_Trho); break;
						case 'L':
							Value=k_func(T,rho,TYPE_Trho); break;
						case 'A':
							Value=w_func(T,rho,TYPE_Trho); break;
						default:
							strcpy(errString,"Invalid Output Name");
							return -100;
					}
					return Value;
				}
				else
				{
					// It is two-phase. Find the quality and call Props again with the quality
					Q=(1/rho-1/rhoL)/(1/rhoV-1/rhoL);
					return Props(Output,'T',T,'Q',Q,Ref);
				}
			}

			// Temperature and quality are the inputs
			else if (Name2=='Q')
			{
				Q=Prop2;
				if (!strcmp(Ref,"R290") || !strcmp(Ref,"Argon") ||!strcmp(Ref,"Nitrogen")  ||!strcmp(Ref,"R134a") ||!strcmp(Ref,"R717") || !strcmp(Ref,"CO2") || !strcmp(Ref,"R744") || !strcmp(Ref,"Water") )
				{
					rhosatPure(Ref,T,&rhoL,&rhoV,&p);
				}
				else
				{
					pbp=pbp_func(T);
					pdp=pdp_func(T);
					rhoV=rhosatV_func(T);
					rhoL=rhosatL_func(T);
					p=Q*pdp+(1-Q)*pbp;
				}
				rho=1/(Q/rhoV+(1-Q)/rhoL);

				switch (Output)
				{
					case 'P':
						Value=p; 
						break;
					case 'H':
						Value=Q*h_func(T,rhoV,TYPE_Trho)+(1-Q)*h_func(T,rhoL,TYPE_Trho); 
						break;
					case 'D':
						Value=rho; break;
					case 'S':
						Value=Q*s_func(T,rhoV,TYPE_Trho)+(1-Q)*s_func(T,rhoL,TYPE_Trho); 
						break;
					case 'U':
						Value=Q*u_func(T,rhoV,TYPE_Trho)+(1-Q)*u_func(T,rhoL,TYPE_Trho); 
						break;
					case 'C':
						Value=Q*cp_func(T,rhoV,TYPE_Trho)+(1-Q)*cp_func(T,rhoL,TYPE_Trho); 
						break;
					case 'O':
						Value=Q*cv_func(T,rhoV,TYPE_Trho)+(1-Q)*cv_func(T,rhoL,TYPE_Trho); 
						break;
					case 'V':
						Value=Q*visc_func(T,rhoV,TYPE_Trho)+(1-Q)*visc_func(T,rhoL,TYPE_Trho); 
						break;
					case 'L':
						Value=Q*k_func(T,rhoV,TYPE_Trho)+(1-Q)*k_func(T,rhoL,TYPE_Trho); 
						break;
					case 'A':
						Value=Q*w_func(T,rhoV,TYPE_Trho)+(1-Q)*w_func(T,rhoL,TYPE_Trho); 
						break;
					case 'M':
						Value=MM_func(); break;
					case 'E':
						Value=pcrit_func(); break;
					case 'B':
						Value=Tcrit_func(); break;
					case 'R':
						Value=Ttriple_func(); break;
					default:
						strcpy(errString,"Invalid Output Name");
						return -100;
				}
				return Value;
			}
		}
		else if (Name1=='P' && Name2=='Q')
		{
			T=Tsat(Ref,Prop1,Prop2,0);
			return Props(Output,'T',T,'Q',Prop2,Ref);
		}
		else if (Output=='T' && Name1=='P' && Name2=='Q')
		{
			return Tsat(Ref,Prop1,Prop2,0);
		}
		else
		{
			fprintf(stderr,"Names of input properties invalid (%c,%c) with refrigerant %s.  Valid choices are T,P or T,Q or T,D or P,Q",Name1,Name2,Ref);
			return _HUGE;
		}
	}
}

double pcrit(char *Ref)
{
	// Call a function to set the global constants
	Props('M','T',0,'P',0,Ref);
	
	// Brines do not have critical pressures, set it to a big number
	if (IsFluidType(Ref,"Brine"))
	{
		return 100000;
	}
	#if defined(_WIN32) || defined(__WIN32__) 
	if (strncmp(Ref,"REFPROP-",8)==0)  // First eight characters match "REFPROP-"
	{
		char *REFPROPRef=NULL,*RefCopy=NULL;
		double pcrit;
		RefCopy=malloc(strlen(Ref)+1);
		strcpy(RefCopy,Ref);
		REFPROPRef = strtok(RefCopy,"-");
		REFPROPRef = strtok(NULL,"-");
		// 'E' is the code for the critical pressure (ran out of sensible characters).  
		pcrit=REFPROP('E','T',0.0,'Q',0.0,REFPROPRef);
		free(RefCopy);
		return pcrit;
	}
	#else
	if (0){} // Automatically skip it because REFPROP is not supported on this platform
	#endif
	else{
		return Props('E','T',0.0,'Q',0.0,Ref);
	}
}

double Tcrit(char *Ref)
{	
	// Call a function to set the global constants
	Props('M','T',0,'P',0,Ref);
	
	// Brines do not have critical temperatures, set it to a big number
	if (IsFluidType(Ref,"Brine"))
	{
		return 100000;
	}
	#if defined(_WIN32) || defined(__WIN32__) 
	if (strncmp(Ref,"REFPROP-",8)==0)  // First eight characters match "REFPROP-"
	{
		char *REFPROPRef=NULL,*RefCopy=NULL;
		double Tc;
		RefCopy=malloc(strlen(Ref)+1);
		strcpy(RefCopy,Ref);
		REFPROPRef = strtok(RefCopy,"-");
		REFPROPRef = strtok(NULL,"-");
		// 'B' is the code for the critical pressure (ran out of sensible characters).
		Tc=REFPROP('B','T',0.0,'Q',0.0,REFPROPRef);
		free(RefCopy);
		return Tc;
	}
	#else
	if (0){} // Automatically skip it because REFPROP is not supported on this platform
	#endif
	else
	{
		return Props('B','T',0.0,'P',0.0,Ref);
	}
}
double Ttriple(char *Ref)
{
	// Call a function to set the global constants
	Props('M','T',0,'P',0,Ref);
	
	// Brines do not have triple point temperatures, set it to a big number
	if (IsFluidType(Ref,"Brine"))
	{
		return 100000;
	}
	#if defined(_WIN32) || defined(__WIN32__) 
	if (strncmp(Ref,"REFPROP-",8)==0)  // First eight characters match "REFPROP-"
	{
		char *REFPROPRef=NULL,*RefCopy=NULL;
		double Ttriple;
		RefCopy=malloc(strlen(Ref)+1);
		strcpy(RefCopy,Ref);
		REFPROPRef = strtok(RefCopy,"-");
		REFPROPRef = strtok(NULL,"-");
		// 'R' is the code for the triple point temperature (ran out of sensible characters).
		Ttriple=REFPROP('R','T',0.0,'Q',0.0,REFPROPRef);
		free(RefCopy);
		return Ttriple;
	}
	#else
	if (0){} // Automatically skip it because REFPROP is not supported on this platform
	#endif
	else
	{
		return Props('R','T',0.0,'Q',0.0,Ref);
	}
}

int errCode(char * Ref)
{
	if (!strcmp(Ref,"R134a"))
		return errCode_R134a();
	if (!strcmp(Ref,"R290"))
		return errCode_R290();
	if (!strcmp(Ref,"R744"))
		return errCode_R744();
	if (!strcmp(Ref,"R717"))
		return errCode_R717();
	if (!strcmp(Ref,"R32"))
		return errCode_R32();
	if (!strcmp(Ref,"R404A"))
		return errCode_R404A();
	if (!strcmp(Ref,"R410A"))
		return errCode_R410A();
	return -1;
}


double T_hp(char *Ref, double h, double p, double T_guess)
{
	double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f,T=300;
	int iter=1;
	while ((iter<=3 || change>eps) && iter<100)
	{
		if (iter==1){x1=T_guess; T=x1;}
		if (iter==2){x2=T_guess+0.1; T=x2;}
		if (iter>2) {T=x2;}
			f=Props('H','T',T,'P',p,Ref)-h;
		if (iter==1){y1=f;}
		if (iter>1)
		{
			y2=f;
			x3=x2-y2/(y2-y1)*(x2-x1);
			change=fabs(y2/(y2-y1)*(x2-x1));
			y1=y2; x1=x2; x2=x3;
		}
		iter=iter+1;
		if (iter>60)
		{
			//printf("%d: T_hp not converging with inputs(%s,%g,%g,%g) value: %0.12g\n",iter,Ref,h,p,T_guess,f);
		}
	}
	return T;
}

double h_sp(char *Ref, double s, double p, double T_guess)
{
	double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f=999,T=300;
	int iter=1;

	
	while ((iter<=3 || change>eps) && iter<100)
	{
		if (iter==1){x1=T_guess; T=x1;}
		if (iter==2){x2=T_guess+1.0; T=x2;}
		if (iter>2) {T=x2;}

			// Find the temperature which gives the same entropy
			f=Props('S','T',T,'P',p,Ref)-s;

		if (iter==1){y1=f;}
		if (iter>1)
		{
			y2=f;
			x3=x2-y2/(y2-y1)*(x2-x1);
			change=fabs(y2/(y2-y1)*(x2-x1));
			y1=y2; x1=x2; x2=x3;
		}
		iter=iter+1;
		if (iter>50)
		{
			printf("h_sp not converging with inputs(%s,%g,%g,%g)\n",Ref,s,p,T_guess);
		}
	}
	return Props('H','T',T,'P',p,Ref);
}

double Tsat(char *Ref, double p, double Q, double T_guess)
{
	double x1=0,x2=0,x3=0,y1=0,y2=0,y3=0,eps=1e-8,change=999,f=999,T=300,tau,tau1,tau3,tau2,logp1,logp2,logp3,T1,T2,T3;
	int iter=1;
	double Tc,Tmax,Tmin;

	// Brines do not have saturation temperatures, set it to a big number
	if (IsFluidType(Ref,"Brine"))
	{
		return 100000;
	}
	#if defined(_WIN32) || defined(__WIN32__)
		// It's a REFPROP fluid - use REFPROP to do all the calculations
		if (!strncmp(Ref,"REFPROP-",8))
		{
			return REFPROP('T','P',p,'Q',Q,Ref);
		}
	#endif
	
	Tc=Tcrit(Ref);
	Tmax=Tc-1;
	Tmin=Ttriple(Ref)+1;
	
	// Plotting Tc/T versus log(p) tends to give very close to straight line
	// Use this fact to figure out a reasonable guess temperature
	
	if (FluidType==FLUIDTYPE_REFRIGERANT_PURE)
	{
		T1=Ttriple(Ref)+1;
		T3=Tcrit(Ref)-1;
		T2=(T1+T3)/2;
		tau1=Tc/T1;
		tau2=Tc/T2;
		tau3=Tc/T3;
		logp1=log(psat_func(T1));
		logp2=log(psat_func(T2));
		logp3=log(psat_func(T3));

		T_guess=Tc/QuadInterpolate(logp1,logp2,logp3,tau1,tau2,tau3,log(p));
		if (T_guess+5<Tmax)
			Tmax=T_guess+5;
		if (T_guess-5>Tmin)
			Tmin=T_guess-5;
	}

	return _Dekker_Tsat(Tmin,Tmax,0.001,p,Q,Ref);
	
	while ((iter<=3 || exp(f)-1>eps) && iter<100)
	{
		if (iter==1){x1=Tc/Tmin; tau=x1;}
		if (iter==2){x2=Tc/Tmax; tau=x2;}
		if (iter>2) {tau=x2;}

			T=Tc/tau;
			f=log(Props('P','T',T,'Q',Q,Ref))-log(p);

			//printf("T: %g f %g\n",T,f);

		if (iter==1){y1=f;}
		if (iter>1)
		{
			y2=f;
			x3=x2-y2/(y2-y1)*(x2-x1);
			change=fabs(y2/(y2-y1)*(x2-x1));
			y1=y2; x1=x2; x2=x3;
		}
		iter=iter+1;
		if (iter>50)
		{
			printf("Tsat not converging with inputs(%s,%g,%g,%g)\n",Ref,p,Q,T_guess);
		}
	}
	return Tc/x2;
}

double K2F(double T)
{
	return T * 9 / 5 - 459.67;
}

double F2K(double T_F)
{
	return (T_F + 459.67) * 5 / 9;
}

void PrintSaturationTable(char *FileName, char * Ref,double Tmin, double Tmax)
{

	double T;
	FILE *f;
	f=fopen(FileName,"w");
	fprintf(f,"T,pL,pV,rhoL,rhoV,uL,uV,hL,hV,sL,sV,cvL,cvV,cpL,cpV,kL,kV,muL,muV\n");
	fprintf(f,"[K],[kPa],[kPa],[kg/m^3],[kg/m^3],[kJ/kg],[kJ/kg],[kJ/kg],[kJ/kg],[kJ/kg-K],[kJ/kg-K],[kJ/kg-K],[kJ/kg-K],[kJ/kg-K],[kJ/kg-K],[kW/m-K],[kW/m-K],[Pa-s],[Pa-s]\n");
	fprintf(f,"--------------------------------------------------------------------------\n");

	for (T=Tmin;T<Tmax;T+=1.0)
	{
	fprintf(f,"%0.3f,",T);
	fprintf(f,"%0.6f,",Props('P','T',T,'Q',0,Ref));
	fprintf(f,"%0.6f,",Props('P','T',T,'Q',1,Ref));
	fprintf(f,"%0.6f,",Props('D','T',T,'Q',0,Ref));
	fprintf(f,"%0.6f,",Props('D','T',T,'Q',1,Ref));
	fprintf(f,"%0.6f,",Props('U','T',T,'Q',0,Ref));
	fprintf(f,"%0.6f,",Props('U','T',T,'Q',1,Ref));
	fprintf(f,"%0.6f,",Props('H','T',T,'Q',0,Ref));
	fprintf(f,"%0.6f,",Props('H','T',T,'Q',1,Ref));
	fprintf(f,"%0.6f,",Props('S','T',T,'Q',0,Ref));
	fprintf(f,"%0.6f,",Props('S','T',T,'Q',1,Ref));
	fprintf(f,"%0.6f,",Props('O','T',T,'Q',0,Ref));
	fprintf(f,"%0.6f,",Props('O','T',T,'Q',1,Ref));
	fprintf(f,"%0.6f,",Props('C','T',T,'Q',0,Ref));
	fprintf(f,"%0.6f,",Props('C','T',T,'Q',1,Ref));
	fprintf(f,"%0.9f,",Props('L','T',T,'Q',0,Ref));
	fprintf(f,"%0.9f,",Props('L','T',T,'Q',1,Ref));
	fprintf(f,"%0.6g,",Props('V','T',T,'Q',0,Ref));
	fprintf(f,"%0.6g,",Props('V','T',T,'Q',1,Ref));
	fprintf(f,"\n");
	}
	fclose(f);
}

double DerivTerms(char *Term,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
	// Pointers to the functions
	double (*dhdT)(double,double,int);
	double (*dpdT)(double,double,int);
	double (*dhdrho)(double,double,int);
	double (*dpdrho)(double,double,int);
	
	// Wire up the pointers for the given refrigerant
	if (!strcmp(Ref,"Argon"))
	{
		dhdT=dhdT_Argon;
		dpdT=dpdT_Argon;
		dhdrho=dhdrho_Argon;
		dpdrho=dpdrho_Argon;
	}
	else if (!strcmp(Ref,"Nitrogen") || !strcmp(Ref,"N2"))
	{
		dhdT=dhdT_Nitrogen;
		dpdT=dpdT_Nitrogen;
		dhdrho=dhdrho_Nitrogen;
		dpdrho=dpdrho_Nitrogen;
	}
	else if (!strcmp(Ref,"R744") || !strcmp(Ref,"CO2"))
	{
		dhdT=dhdT_R744;
		dpdT=dpdT_R744;
		dhdrho=dhdrho_R744;
		dpdrho=dpdrho_R744;
	}
	else if (!strcmp(Ref,"R410A"))
	{
		dhdT=dhdT_R410A;
		dpdT=dpdT_R410A;
		dhdrho=dhdrho_R410A;
		dpdrho=dpdrho_R410A;
	}
	else
	{
		printf("Bad Refrigerant Name in DerivTerms [%s]\n",Ref);
	}
	
	if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdT"))
		return dhdT(Prop1,Prop2,TYPE_TPNoLookup);
	else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdrho"))
		return dhdrho(Prop1,Prop2,TYPE_TPNoLookup);
	else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdT"))
		return dpdT(Prop1,Prop2,TYPE_TPNoLookup);
	else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdrho"))
		return dpdrho(Prop1,Prop2,TYPE_TPNoLookup);
	
	else if (Name1=='T' && Name2=='D' && !strcmp(Term,"dhdT"))
		return dhdT(Prop1,Prop2,TYPE_Trho);
	else if (Name1=='T' && Name2=='D' && !strcmp(Term,"dhdrho"))
		return dhdrho(Prop1,Prop2,TYPE_Trho);
	else if (Name1=='T' && Name2=='D' && !strcmp(Term,"dpdT"))
		return dpdT(Prop1,Prop2,TYPE_Trho);
	else if (Name1=='T' && Name2=='D' && !strcmp(Term,"dpdrho"))
		return dpdrho(Prop1,Prop2,TYPE_Trho);
	
	else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdTnum"))
		return dhdT(Prop1,Prop2,99);
	else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdrhonum"))
		return dhdrho(Prop1,Prop2,99);
	else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdTnum"))
		return dpdT(Prop1,Prop2,99);
	else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdrhonum"))
		return dpdrho(Prop1,Prop2,99);

	else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdTnum"))
		return dhdT(Prop1,Prop2,99);
	else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dhdrhonum"))
		return dhdrho(Prop1,Prop2,99);
	else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdTnum"))
		return dpdT(Prop1,Prop2,99);
	else if (Name1=='T' && Name2=='P' && !strcmp(Term,"dpdrhonum"))
		return dpdrho(Prop1,Prop2,99);
	else
	{
		printf("Bad pair of properties[%c,%c] and derivative name [%s]\n",Name1,Name2,Term);
		return _HUGE;
	}	
}

static void swap(double *x, double *y)
{
    double temp;
    temp = *x;
    *x = *y;
    *y = temp;
}

double _Dekker_Tsat(double x_min, double x_max, double eps, double p, double Q,char *Ref)
{
	double a_k,b_k,f_ak,f_bk,f_bkn1,error, 
		b_kn1,b_kp1,s,m,a_kp1,f_akp1,f_bkp1,x,Tc;

	int iter=1;

	Tc=Tcrit_func();
	x_min=Tc/x_min;
	x_max=Tc/x_max;

	// Loop for the solver method
    while ((iter <= 1 || fabs(error) > eps) && iter < 100)
	{
		// Start with the maximum value
        if (iter == 1)
		{
            a_k = x_max;
            x = a_k;
		}
		// End with the minimum value
        if (iter == 2)
		{
            b_k = x_min;
            x = b_k;
		}
        if (iter > 2)
            x = b_k;

		// Evaluate residual
        error = log(Props('P','T',Tc/x,'Q',Q,Ref))-log(p);

		// First time through, store the outputs
        if(iter == 1)
		{
            f_ak = error;
            b_kn1 = a_k;
            f_bkn1 = error;
		}
        if (iter > 1)
		{
            f_bk = error;
			//Secant solution
            s = b_k - (b_k - b_kn1) / (f_bk - f_bkn1) * f_bk;
			//Midpoint solution
            m = (a_k + b_k) / 2.0;

            if (s > b_k && s < m)
			{
                //Use the secant solution
                b_kp1 = s;
			}
            else
			{
				//Use the midpoint solution
                b_kp1 = m;
			}

            //See if the signs of iterate and contrapoint are the same
            f_bkp1 = log(Props('P','T',Tc/b_kp1,'Q',Q,Ref))-log(p);

            if (f_ak / fabs(f_ak) != f_bkp1 / fabs(f_bkp1))
			{
                // If a and b have opposite signs, 
				//  keep the same contrapoint
                a_kp1 = a_k;
                f_akp1 = f_ak;
			}
            else
			{
                //Otherwise, keep the iterate
                a_kp1 = b_k;
                f_akp1 = f_bk;
			}

            if (fabs(f_akp1) < fabs(f_bkp1))
			{
                //a_k+1 is a better guess than b_k+1, so swap a and b values
				swap(&a_kp1, &b_kp1);
                swap(&f_akp1, &f_bkp1);
			}

            //Update variables
            //Old values
            b_kn1 = b_k;
            f_bkn1 = f_bk;
            //values at this iterate
            b_k = b_kp1;
            a_k = a_kp1;
            f_ak = f_akp1;
            f_bk = f_bkp1;
		}
        iter++;
		if (iter>90 && fabs(error)>eps)
		{
			printf("Dekker for Tsat has failed with inputs (%g,%g,%g,%g,%g,%s); value: %g\n",x_min,x_max,eps,p,Q,Ref,error);
		}
	}
	return Tc/b_k;
}
