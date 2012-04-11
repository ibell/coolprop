/*

Properties of Water (R718)
by Ian Bell

Themo properties from 
"The IAPWS Formulation 1995 for the Thermodynamic Properties
of Ordinary Water Substance for General and Scientific Use", 
W. Wagner and A. Pruss, J. Phys. Chem. Ref. Data, v. 31, 2000

*/

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
#include "CoolProp.h"


static const double Tc=647.096, M_Water=18.015268,rhoc=322,     Pc=22064, _Ttriple=273.16;
             //    K                    g/mol          kg/m^3      kPa             K
static const double n[]={0,    
	 0.125335479355233e-1, //[1]
	 0.78957634722828e1, //[2]
	-0.87803203303561e1, //[3]
	 0.31802509345418, //[4]
	-0.26145533859358, //[5]
	-0.78199751687981e-2, //[6]
	 0.88089493102134e-2, //[7]
	-0.66856572307965, //[8]
	 0.20433810950965, //[9]
	-0.66212605039687e-4, //[10]
	-0.19232721156002, //[11]
	-0.25709043003438, //[12]
	 0.16074868486251, //[13]
	-0.40092828925807e-1, //[14]
	 0.39343422603254e-6, //[15]
	-0.75941377088144e-5, //[16]
	 0.56250979351888e-3, //[17]
	-0.15608652257135e-4, //[18]
	 0.11537996422951e-8, //[19]
	 0.36582165144204e-6, //[20]
	-0.13251180074668e-11, //[21]
	-0.62639586912454e-9, //[22]
	-0.10793600908932, //[23]
	 0.17611491008752e-1, //[24]
	 0.22132295167546, //[25]
	-0.40247669763528, //[26]
	 0.58083399985759, //[27]
	 0.49969146990806e-2, //[28]
	-0.31358700712549e-1, //[29]
	-0.74315929710341, //[30]
	 0.47807329915480, //[31]
	 0.20527940895948e-1, //[32]
	-0.13636435110343, //[33]
	 0.14180634400617e-1, //[34]
	 0.83326504880713e-2, //[35]
	-0.29052336009585e-1, //[36]
	 0.38615085574206e-1, //[37]
	-0.20393486513704e-1, //[38]
	-0.16554050063734e-2, //[39]
	 0.19955571979541e-2, //[40]
	 0.15870308324157e-3, //[41]
	-0.16388568342530e-4, //[42]
	 0.43613615723811e-1, //[43]
	 0.34994005463765e-1, //[44]
	-0.76788197844621e-1, //[45]
	 0.22446277332006e-1, //[46]
	-0.62689710414685e-4, //[47]
	-0.55711118565645e-9, //[48]
	-0.19905718354408, //[49]
	 0.31777497330738, //[50]
	-0.11841182425981, //[51]
	-0.31306260323435e2, //[52]
	 0.31546140237781e2, //[53]
	-0.25213154341695e4, //[54]
	-0.14874640856724, //[55]
	 0.31806110878444, //[56]
};

static const int d[]={0,
	1, //[1]
	1, //[2]
	1, //[3]
	2, //[4]
	2, //[5]
	3, //[6]
	4, //[7]
	1, //[8]
	1, //[9]
	1, //[10]
	2, //[11]
	2, //[12]
	3, //[13]
	4, //[14]
	4, //[15]
	5, //[16]
	7, //[17]
	9, //[18]
	10, //[19]
	11, //[20]
	13, //[21]
	15, //[22]
	1, //[23]
	2, //[24]
	2, //[25]
	2, //[26]
	3, //[27]
	4, //[28]
	4, //[29]
	4, //[30]
	5, //[31]
	6, //[32]
	6, //[33]
	7, //[34]
	9, //[35]
	9, //[36]
	9, //[37]
	9, //[38]
	9, //[39]
	10, //[40]
	10, //[41]
	12, //[42]
	3, //[43]
	4, //[44]
	4, //[45]
	5, //[46]
	14, //[47]
	3, //[48]
	6, //[49]
	6, //[50]
	6, //[51]
	3, //[52]
	3, //[53]
	3, //[54]
};

static const double t[]={0.00,
	-0.5, //[1]
	0.875, //[2]
	1, //[3]
	0.5, //[4]
	0.75, //[5]
	0.375, //[6]
	1, //[7]
	4, //[8]
	6, //[9]
	12, //[10]
	1, //[11]
	5, //[12]
	4, //[13]
	2, //[14]
	13, //[15]
	9, //[16]
	3, //[17]
	4, //[18]
	11, //[19]
	4, //[20]
	13, //[21]
	1, //[22]
	7, //[23]
	1, //[24]
	9, //[25]
	10, //[26]
	10, //[27]
	3, //[28]
	7, //[29]
	10, //[30]
	10, //[31]
	6, //[32]
	10, //[33]
	10, //[34]
	1, //[35]
	2, //[36]
	3, //[37]
	4, //[38]
	8, //[39]
	6, //[40]
	9, //[41]
	8, //[42]
	16, //[43]
	22, //[44]
	23, //[45]
	23, //[46]
	10, //[47]
	50, //[48]
	44, //[49]
	46, //[50]
	50, //[51]
	0, //[52]
	1, //[53]
	4, //[54]
};

static const int c[]={0,
	0, //[1]
	0, //[2]
	0, //[3]
	0, //[4]
	0, //[5]
	0, //[6]
	0, //[7]
	1, //[8]
	1, //[9]
	1, //[10]
	1, //[11]
	1, //[12]
	1, //[13]
	1, //[14]
	1, //[15]
	1, //[16]
	1, //[17]
	1, //[18]
	1, //[19]
	1, //[20]
	1, //[21]
	1, //[22]
	2, //[23]
	2, //[24]
	2, //[25]
	2, //[26]
	2, //[27]
	2, //[28]
	2, //[29]
	2, //[30]
	2, //[31]
	2, //[32]
	2, //[33]
	2, //[34]
	2, //[35]
	2, //[36]
	2, //[37]
	2, //[38]
	2, //[39]
	2, //[40]
	2, //[41]
	2, //[42]
	3, //[43]
	3, //[44]
	3, //[45]
	3, //[46]
	4, //[47]
	6, //[48]
	6, //[49]
	6, //[50]
	6, //[51]
	0, //[52]
	0, //[53]
	0, //[54]
};

static const double alpha[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 51
	20, //[52]
	20, //[53]
	20, //[54]
};

static const double beta[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 51
	150, //[52]
	150, //[53]
	250, //[54]
	0.3, //[55]
	0.3, //[56]
};

static const double GAMMA[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 51
	1.21, //[52]
	1.21, //[53]
	1.25, //[54]
};

static const double epsilon[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 51
	1, //[52]
	1, //[53]
	1, //[54]
};

static const double a[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 54
	3.5, //[55]
	3.5, //[56]
};

static const double b[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 54
	0.85, //[55]
	0.95, //[56]
};

static const double A[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 54
	0.32, //[55]
	0.32, //[56]
};

static const double B[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 54
	0.2, //[55]
	0.2, //[56]
};

static const double C[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 54
	28., //[55]
	32., //[56]
};

static const double D[]={
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //0 to 54
	700.0, //[55]
	800.0, //[56]
};

// Constant for ideal gas part
static const double n0[]={0,
	-8.3204464837497, //[1] //From Herrmann 2099 ASHRAE RP-1485
	6.6832105275932, //[2] //From Herrmann 2099 ASHRAE RP-1485
	3.00632, //[3]
	0.012436, //[4]
	0.97315, //[5]
	1.27950, //[6]
	0.96956, //[7]
	0.24873, //[8]
};
static const double gamma0[]={0,
	0,0,0, // 1 to 3
	1.28728967, //[4]
	3.53734222, //[5]
	7.74073708, //[6]
	9.24437796, //[7]
	27.5075105 //[8]
};

int Load_Water(struct fluidParamsVals *Fluid)
{
    // Function pointers
    Fluid->funcs.phir=phir_Water;
    Fluid->funcs.dphir_dDelta=dphir_dDelta_Water;
    Fluid->funcs.dphir2_dDelta2=dphir2_dDelta2_Water;
    Fluid->funcs.dphir2_dDelta_dTau=dphir2_dDelta_dTau_Water;
    Fluid->funcs.dphir_dTau=dphir_dTau_Water;
    Fluid->funcs.dphir2_dTau2=dphir2_dTau2_Water;
    Fluid->funcs.phi0=phi0_Water;
    Fluid->funcs.dphi0_dDelta=dphi0_dDelta_Water;
    Fluid->funcs.dphi02_dDelta2=dphi02_dDelta2_Water;
    Fluid->funcs.dphi0_dTau=dphi0_dTau_Water;
    Fluid->funcs.dphi02_dTau2=dphi02_dTau2_Water;
    Fluid->funcs.rhosatL=rhosatL_Water;
    Fluid->funcs.rhosatV=rhosatV_Water;
    Fluid->funcs.psat=psat_Water;

    Fluid->funcs.visc=Viscosity_Trho_Water;
    Fluid->funcs.cond=Conductivity_Trho_Water;

    //Lookup table parameters
    Fluid->LUT.Tmin=200.0;
    Fluid->LUT.Tmax=800.0;
    Fluid->LUT.pmin=500;
    Fluid->LUT.pmax=16000;

    //Fluid parameters
    Fluid->Type=FLUIDTYPE_REFRIGERANT_PURE;
    Fluid->Tc=Tc;
    Fluid->rhoc=rhoc;
    Fluid->MM=M_Water;
    Fluid->pc=Pc;
    Fluid->Tt=_Ttriple;
    return 1;
}


double rhosatV_Water(double T)
{
    const double ti[]={0,2.0/6.0,4.0/6.0,8./6.0,18.0/6.0,37.0/6.0,71.0/6.0};
    const double ai[]={0,-2.03150240,-2.68302940,-5.38626492,-17.2991605,-44.7586581,-63.9201063};
    double summer=0;
    int i;
    for (i=1;i<=6;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*exp(summer);
}

double rhosatL_Water(double T)
{
    const double ti[]={0,1.0/3.0,2.0/3.0,5.0/3.0,16.0/3.0,43.0/3.0,110.0/3.0};
    const double ai[]={0,1.99274064,1.09965342,-0.510839303,-1.75493479,-45.5170352,-1.75493479,-45.5170352,-6.74694450e5};
    double summer=1;
    int i;
    for (i=1;i<=6;i++)
    {
        summer=summer+ai[i]*pow(1.0-T/Tc,ti[i]);
    }
    return rhoc*summer;
}

double psat_Water(double T)
{
    const double ti[]={0,1.0,1.5,3.0,3.5,4.0,7.5};
    const double ai[]={0,-7.85951783,1.84408259,-11.7866497,22.6807411,-15.9618719,1.80122502};
    double summer=0;
    int i;
    for (i=1;i<=6;i++)
    {
        summer=summer+ai[i]*pow(1-T/Tc,ti[i]);
    }
    return Pc*exp(Tc/T*summer);
}

double IsothermCompress_Water(double T, double p)
{
    // Isothermal compressibility with units of 1/Pa
    double rho,delta,tau,dpdrho,R_Water;
    rho=Props('D','T',T,'P',p,"Water");
    delta=rho/rhoc;
    tau=Tc/T;
    R_Water=8.31447215/M_Water;
    dpdrho=R_Water*T*(1.0+2.0*delta*dphir_dDelta_Water(tau,delta)+powInt(delta,2)*dphir2_dDelta2_Water(tau,delta));
    return 1/(rho*dpdrho*1000.0);
}

static void visc_Helper(double Tbar, double rhobar, double *mubar_0, double *mubar_1)
{
	double H[6][7],sum;
	int i,j;

	// Dilute-gas component
	*mubar_0=100.0*sqrt(Tbar)/(1.67752+2.20462/Tbar+0.6366564/powInt(Tbar,2)-0.241605/powInt(Tbar,3));

	//Fill in zeros in H
	for (i=0;i<=5;i++)
	{
		for (j=0;j<=6;j++)
		{
			H[i][j]=0;
		}
	}

	//Set non-zero parameters of H
	H[0][0]=5.20094e-1;
	H[1][0]=8.50895e-2;
	H[2][0]=-1.08374;
	H[3][0]=-2.89555e-1;

	H[0][1]=2.22531e-1;
	H[1][1]=9.99115e-1;
	H[2][1]=1.88797;
	H[3][1]=1.26613;
	H[5][1]=1.20573e-1;

	H[0][2]=-2.81378e-1;
	H[1][2]=-9.06851e-1;
	H[2][2]=-7.72479e-1;
	H[3][2]=-4.89837e-1;
	H[4][2]=-2.57040e-1;

	H[0][3]=1.61913e-1;
	H[1][3]=2.57399e-1;

	H[0][4]=-3.25372e-2;
	H[3][4]=6.98452e-2;

	H[4][5]=8.72102e-3;

	H[3][6]=-4.35673e-3;
	H[5][6]=-5.93264e-4;

	// Finite density component
	sum=0;
	for (i=0;i<=5;i++)
	{
		for (j=0;j<=6;j++)
		{
			sum+=powInt(1/Tbar-1,i)*(H[i][j]*powInt(rhobar-1,j));
		}
	}
	*mubar_1=exp(rhobar*sum);
}

double Viscosity_Trho_Water(double T,double rho)
{
	double x_mu=0.068,qc=1/1.9,qd=1/1.1,nu=0.630,gamma=1.239,zeta_0=0.13,LAMBDA_0=0.06,Tbar_R=1.5;
	double delta,tau,mubar_0,mubar_1,mubar_2,drhodp,drhodp_R,DeltaChibar,zeta,w,L,Y,psi_D,Tbar,rhobar;
	double drhobar_dpbar,drhobar_dpbar_R,R_Water;
	
	Tbar=T/Tc;
	rhobar=rho/rhoc;
	R_Water=8.31447215/M_Water;

	// Dilute and finite gas portions
	visc_Helper(Tbar,rhobar,&mubar_0,&mubar_1);

	///************************ Critical Enhancement ************************
	delta=rhobar;
	// "Normal" calculation
	tau=1/Tbar;
	drhodp=1/(R_Water*T*(1+2*delta*dphir_dDelta_Water(tau,delta)+delta*delta*dphir2_dDelta2_Water(tau,delta)));
	drhobar_dpbar=Pc/rhoc*drhodp;
	// "Reducing" calculation
	tau=1/Tbar_R;
	drhodp_R=1/(R_Water*Tbar_R*Tc*(1+2*delta*dphir_dDelta_Water(tau,delta)+delta*delta*dphir2_dDelta2_Water(tau,delta)));
	drhobar_dpbar_R=Pc/rhoc*drhodp_R;
	
	DeltaChibar=rhobar*(drhobar_dpbar-drhobar_dpbar_R*Tbar_R/Tbar);
	if (DeltaChibar<0)
		DeltaChibar=0;
	zeta=zeta_0*pow(DeltaChibar/LAMBDA_0,nu/gamma);
	if (zeta<0.3817016416)
	{
		Y=1.0/5.0*qc*zeta*powInt(qd*zeta,5)*(1-qc*zeta+powInt(qc*zeta,2)-765.0/504.0*powInt(qd*zeta,2));
	}
	else
	{
		psi_D=acos(pow(1+powInt(qd*zeta,2),-1.0/2.0));
		w=sqrt(fabs((qc*zeta-1)/(qc*zeta+1)))*tan(psi_D/2.0);
		if (qc*zeta>1)
		{
			L=log((1+w)/(1-w));
		}
		else
		{
			L=2*atan(fabs(w));
		}
		Y=1.0/12.0*sin(3*psi_D)-1/(4*qc*zeta)*sin(2*psi_D)+1.0/powInt(qc*zeta,2)*(1-5.0/4.0*powInt(qc*zeta,2))*sin(psi_D)-1.0/powInt(qc*zeta,3)*((1-3.0/2.0*powInt(qc*zeta,2))*psi_D-pow(fabs(powInt(qc*zeta,2)-1),3.0/2.0)*L);
	}
	mubar_2=exp(x_mu*Y);

	return (mubar_0*mubar_1*mubar_2)/1e6;
}

double Conductivity_Trho_Water(double T,double rho)
{
	double L[5][6]={{1.3293046,-0.40452437,0.2440949,0.018660751,-0.12961068,0.044809953},
	{1.7018363,-2.2156845,1.6511057,-0.76736002,0.37283344,-0.1120316},
	{5.2246158,-10.124111,4.9874687,-0.27297694,-0.43083393,0.13333849},
	{8.7127675,-9.5000611,4.3786606,-0.91783782,0,0},
	{-1.8525999,0.9340469,0,0,0,0}};
	double lambdabar_0,lambdabar_1,lambdabar_2,rhobar,Tbar,sum,R_Water;
	double Tstar=647.226,rhostar=317.763,pstar=22115,lambdastar=0.0004945;
	double drhodp,drhobar_dpbar,CHIbar_T,dpbar_dTbar,mubar_0,mubar_1,delta,tau;
	int i,j;

	Tbar=T/Tstar;
	rhobar=rho/rhostar;
	R_Water=8.31447215/M_Water;

	// Dilute gas contribution
	lambdabar_0=sqrt(Tbar)/(1.000+6.978267/Tbar+2.599096/powInt(Tbar,2)-0.998254/powInt(Tbar,3));

	sum=0;
	for (i=0;i<=4;i++)
	{
		for (j=0;j<=5;j++)
		{
			//printf("L[%d][%d]: %0.12fg\n",i,j,L[i][j]);
			sum+=L[i][j]*powInt(1.0/Tbar-1.0,i)*powInt(rhobar-1,j);
		}
	}
	// Finite density contribution
	lambdabar_1=exp(rhobar*sum);

	delta=rho/rhoc;
	tau=T/Tc;
	drhodp=1/(R_Water*T*(1+2*delta*dphir_dDelta_Water(tau,delta)+delta*delta*dphir2_dDelta2_Water(tau,delta)));
	drhobar_dpbar=pstar/rhostar*drhodp;
	CHIbar_T=rhobar*drhobar_dpbar;
	dpbar_dTbar=Tstar/pstar*rho*R_Water*(1+delta*dphir_dDelta_Water(tau,delta)-delta*tau*dphir2_dDelta_dTau_Water(tau,delta));

	// Dilute and finite gas portions
	visc_Helper(Tbar,rhobar,&mubar_0,&mubar_1);

	lambdabar_2=55.071*0.0013848/(mubar_0*mubar_1)*powInt(Tbar/rhobar,2)*powInt(dpbar_dTbar,2)*pow(CHIbar_T,0.4678)*sqrt(rhobar)*exp(-18.66*powInt(Tbar-1,2)-powInt(rhobar-1,4));
	return (lambdabar_0*lambdabar_1+lambdabar_2)*lambdastar;
}

double B_Water(double tau)
{
	// given by B*rhoc=lim(delta --> 0) [dphir_ddelta(tau)]
	return 1.0/rhoc*dphir_dDelta_Water(tau,1e-12);
}

double dBdT_Water(double tau)
{
	// given by B*rhoc^2=lim(delta --> 0) [dphir2_ddelta2(tau)]
	return -1.0/rhoc*tau*tau/Tc*dphir2_dDelta_dTau_Water(tau,1e-12);
}

double C_Water(double tau)
{
	// given by B*rhoc^2=lim(delta --> 0) [dphir2_ddelta2(tau)]
	return 1.0/(rhoc*rhoc)*dphir2_dDelta2_Water(tau,1e-12);
}

double dCdT_Water(double tau)
{
	// given by B*rhoc^2=lim(delta --> 0) [dphir2_ddelta2(tau)]
	return -1.0/(rhoc*rhoc)*tau*tau/Tc*dphir3_dDelta2_dTau_Water(tau,1e-12);
}

double phir_Water(double tau, double delta)
{ 
    
    int i;
    double phir=0,theta,DELTA,PSI,psi;
    
    for (i=1;i<=7;i++)
    {
        phir+=n[i]*powInt(delta,d[i])*pow(tau,t[i]);
    }
    
    for (i=8;i<=51;i++)
    {
        phir+=n[i]*powInt(delta,d[i])*pow(tau,t[i])*exp(-powInt(delta,c[i]));
    }
    
    for (i=52;i<=54;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        phir+=n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi;
    }
    
    for (i=55;i<=56;i++)
    {
        theta=(1.0-tau)+A[i]*pow(powInt(delta-1.0,2),1/(2*beta[i]));
        DELTA=powInt(theta,2)+B[i]*pow(powInt(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powInt(delta-1,2)-D[i]*powInt(tau-1,2));
        phir+=n[i]*pow(DELTA,b[i])*delta*PSI;
    }
    return phir;
}

double dphir_dDelta_Water(double tau, double delta)
{ 
    int i;
    double dphir_dDelta=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta,psi;
    for (i=1;i<=7;i++)
    {
        dphir_dDelta=dphir_dDelta+n[i]*d[i]*powInt(delta,d[i]-1)*pow(tau,t[i]);
    }
    for (i=8;i<=51;i++)
    {
        dphir_dDelta=dphir_dDelta+n[i]*exp(-powInt(delta,c[i]))*(powInt(delta,d[i]-1)*pow(tau,t[i])*(d[i]-c[i]*powInt(delta,c[i])));
    }
    for (i=52;i<=54;i++)
    {        
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir_dDelta=dphir_dDelta+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(d[i]/delta-2.0*alpha[i]*(delta-epsilon[i]));
    }
    for (i=55;i<=56;i++)
    {
        theta=(1.0-tau)+A[i]*pow(powInt(delta-1.0,2),1/(2*beta[i]));
        DELTA=powInt(theta,2)+B[i]*pow(powInt(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powInt(delta-1,2)-D[i]*powInt(tau-1,2));

        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(powInt(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        dphir_dDelta=dphir_dDelta+n[i]*(pow(DELTA,b[i])*(PSI+delta*dPSI_dDelta)+dDELTAbi_dDelta*delta*PSI);
    }
    return dphir_dDelta;
}

double dphir2_dDelta2_Water(double tau, double delta)
{ 
    
    int i;
    double di,ci;
    
    double dphir2_dDelta2=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta,psi,dPSI2_dDelta2,dDELTAbi2_dDelta2,dDELTA2_dDelta2;
    for (i=1;i<=7;i++)
    {
        di=(double)d[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*di*(di-1.0)*powInt(delta,d[i]-2)*pow(tau,t[i]);
    }
    for (i=8;i<=51;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta2=dphir2_dDelta2+n[i]*exp(-powInt(delta,c[i]))*(powInt(delta,d[i]-2)*pow(tau,t[i])*( (di-ci*powInt(delta,c[i]))*(di-1.0-ci*powInt(delta,c[i])) - ci*ci*powInt(delta,c[i])));
    }
    for (i=52;i<=54;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dDelta2=dphir2_dDelta2+n[i]*pow(tau,t[i])*psi*(-2.0*alpha[i]*powInt(delta,d[i])+4.0*powInt(alpha[i],2)*powInt(delta,d[i])*powInt(delta-epsilon[i],2)-4.0*di*alpha[i]*powInt(delta,d[i]-1)*(delta-epsilon[i])+di*(di-1.0)*powInt(delta,d[i]-2));
    }
    for (i=55;i<=56;i++)
    {
               
        theta=(1.0-tau)+A[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=powInt(theta,2)+B[i]*pow(powInt(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powInt(delta-1.0,2)-D[i]*powInt(tau-1.0,2));
        
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(powInt(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
        dPSI2_dDelta2=(2.0*C[i]*powInt(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
        dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+powInt(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(powInt(delta-1.0,2),a[i]-2.0)+2.0*powInt(A[i]/beta[i],2)*powInt(pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
        dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*powInt(dDELTA_dDelta,2));
        
        dphir2_dDelta2=dphir2_dDelta2+n[i]*(pow(DELTA,b[i])*(2.0*dPSI_dDelta+delta*dPSI2_dDelta2)+2.0*dDELTAbi_dDelta*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta2*delta*PSI);
    }
    return dphir2_dDelta2;
}

    
double dphir2_dDelta_dTau_Water(double tau, double delta)
{
    int i;
    double di, ci;
    double dphir2_dDelta_dTau=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta,psi;
    double dPSI2_dDelta_dTau, dDELTAbi2_dDelta_dTau, dPSI_dTau, dDELTAbi_dTau;
    for (i=1;i<=7;i++)
    {
        di=(double)d[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*di*t[i]*powInt(delta,d[i]-1)*pow(tau,t[i]-1.0);
    }
    for (i=8;i<=51;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*exp(-powInt(delta,c[i]))*powInt(delta,d[i]-1)*t[i]*pow(tau,t[i]-1.0)*(di-ci*powInt(delta,c[i]));
    }
    for (i=52;i<=54;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(di/delta-2.0*alpha[i]*(delta-epsilon[i]))*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    for (i=55;i<=56;i++)
    {
        
        theta=(1.0-tau)+A[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=powInt(theta,2)+B[i]*pow(powInt(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powInt(delta-1.0,2)-D[i]*powInt(tau-1.0,2));
        
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(powInt(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        
        dPSI2_dDelta_dTau=4.0*C[i]*D[i]*(delta-1.0)*(tau-1.0)*PSI;
        dDELTAbi2_dDelta_dTau=-A[i]*b[i]*2.0/beta[i]*pow(DELTA,b[i]-1.0)*(delta-1.0)*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0)-2.0*theta*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0)*dDELTA_dDelta;
        
        dphir2_dDelta_dTau=dphir2_dDelta_dTau+n[i]*(pow(DELTA,b[i])*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+delta*dDELTAbi_dDelta*dPSI_dTau+ dDELTAbi_dTau*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta_dTau*delta*PSI);
    }
    return dphir2_dDelta_dTau;
}

double dphir_dTau_Water(double tau, double delta)
{ 
    
    int i;
    double dphir_dTau=0,theta,DELTA,PSI,dPSI_dTau,dDELTAbi_dTau,psi;
    
    for (i=1;i<=7;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powInt(delta,d[i])*pow(tau,t[i]-1.0);
    }
    
    for (i=8;i<=51;i++)
    {
        dphir_dTau=dphir_dTau+n[i]*t[i]*powInt(delta,d[i])*pow(tau,t[i]-1.0)*exp(-powInt(delta,c[i]));
    }
    
    for (i=52;i<=54;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir_dTau=dphir_dTau+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]));
    }
    
    for (i=55;i<=56;i++)
    {
        theta=(1.0-tau)+A[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=powInt(theta,2)+B[i]*pow(powInt(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powInt(delta-1.0,2)-D[i]*powInt(tau-1.0,2));
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        dphir_dTau=dphir_dTau+n[i]*delta*(dDELTAbi_dTau*PSI+pow(DELTA,b[i])*dPSI_dTau);
    }
    return dphir_dTau;
}


double dphir2_dTau2_Water(double tau, double delta)
{ 
    int i;
    double dphir2_dTau2=0,theta,DELTA,PSI,dPSI_dTau,dDELTAbi_dTau,psi,dPSI2_dTau2,dDELTAbi2_dTau2;
    
    for (i=1;i<=7;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powInt(delta,d[i])*pow(tau,t[i]-2.0);
    }
    
    for (i=8;i<=51;i++)
    {
        dphir2_dTau2=dphir2_dTau2+n[i]*t[i]*(t[i]-1.0)*powInt(delta,d[i])*pow(tau,t[i]-2.0)*exp(-powInt(delta,c[i]));
    }
    
    for (i=52;i<=54;i++)
    {
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir2_dTau2=dphir2_dTau2+n[i]*powInt(delta,d[i])*pow(tau,t[i])*psi*(powInt(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]),2)-t[i]/powInt(tau,2)-2.0*beta[i]);
    }
    
    for (i=55;i<=56;i++)
    {
        theta=(1.0-tau)+A[i]*pow(powInt(delta-1.0,2),1/(2*beta[i]));
        DELTA=powInt(theta,2)+B[i]*pow(powInt(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powInt(delta-1.0,2)-D[i]*powInt(tau-1.0,2));
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        dPSI2_dTau2=(2.0*D[i]*powInt(tau-1.0,2)-1.0)*2.0*D[i]*PSI;
        dDELTAbi2_dTau2=2.0*b[i]*pow(DELTA,b[i]-1.0)+4.0*powInt(theta,2)*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0);
        dphir2_dTau2=dphir2_dTau2+n[i]*delta*(dDELTAbi2_dTau2*PSI+2.0*dDELTAbi_dTau*dPSI_dTau+pow(DELTA,b[i])*dPSI2_dTau2);
    }
    
    return dphir2_dTau2;
}

double dphir3_dDelta2_dTau_Water(double tau, double delta)
{ 
    // Derivation from Hermann 2009 (ASHRAE Report RP-1485)
    int i;
    double di, ci;
    double dphir3_dDelta2_dTau=0,theta,DELTA,PSI,dPSI_dDelta,dDELTA_dDelta,dDELTAbi_dDelta,psi,dPSI2_dDelta2,dDELTAbi2_dDelta2,dDELTA2_dDelta2;
    double dPSI2_dDelta_dTau, dDELTAbi2_dDelta_dTau, dPSI_dTau, dDELTAbi_dTau;
    double Line1,Line2,Line3,dDELTA2_dDelta_dTau,dDELTA3_dDelta2_dTau,dDELTAbim1_dTau,dDELTAbim2_dTau;
    double dDELTA_dTau,dDELTAbi3_dDelta2_dTau;
    for (i=1;i<=7;i++)
    {
        di=(double)d[i];
        dphir3_dDelta2_dTau+=n[i]*di*t[i]*(di-1)*powInt(delta,d[i]-2)*pow(tau,t[i]-1.0);
    }
    for (i=8;i<=51;i++)
    {
        di=(double)d[i];
        ci=(double)c[i];
        dphir3_dDelta2_dTau+=n[i]*t[i]*exp(-powInt(delta,c[i]))*powInt(delta,d[i]-2)*pow(tau,t[i]-1.0)*(((di-ci*powInt(delta,c[i])))*(di-1-c[i]*powInt(delta,c[i]))-c[i]*c[i]*powInt(delta,c[i]));
    }
    for (i=52;i<=54;i++)
    {
        di=(double)d[i];
        psi=exp(-alpha[i]*powInt(delta-epsilon[i],2)-beta[i]*powInt(tau-GAMMA[i],2));
        dphir3_dDelta2_dTau+=n[i]*pow(tau,t[i])*psi*powInt(delta,d[i])*(t[i]/tau-2.0*beta[i]*(tau-GAMMA[i]))*(-2.0*alpha[i]*(delta-epsilon[i])+4*alpha[i]*alpha[i]*pow(delta,d[i])*powInt(delta-epsilon[i],2)-4*d[i]*alpha[i]*pow(delta,d[i]-1)*(delta-epsilon[i])+d[i]*(d[i]-1)*powInt(delta,d[i]-2));
    }
    for (i=55;i<=56;i++)
    {
        
        theta=(1.0-tau)+A[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i]));
        DELTA=powInt(theta,2)+B[i]*pow(powInt(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*powInt(delta-1.0,2)-D[i]*powInt(tau-1.0,2));
        
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(powInt(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
        dPSI2_dDelta2=(2.0*C[i]*powInt(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
        dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+powInt(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(powInt(delta-1.0,2),a[i]-2.0)+2.0*powInt(A[i]/beta[i],2)*powInt(pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
        dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*powInt(dDELTA_dDelta,2));
        
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        
        dPSI2_dDelta_dTau=4.0*C[i]*D[i]*(delta-1.0)*(tau-1.0)*PSI;
        dDELTAbi2_dDelta_dTau=-A[i]*b[i]*2.0/beta[i]*pow(DELTA,b[i]-1.0)*(delta-1.0)*pow(powInt(delta-1.0,2),1.0/(2.0*beta[i])-1.0)-2.0*theta*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0)*dDELTA_dDelta;
	
	//Following Terms added for this derivative
	dDELTA_dTau=-2*((1-tau)+A[i]*pow(powInt(delta-1,2),1/(2*beta[i])-1)+2*B[i]*a[i]*pow(powInt(delta-1,2),a[i]-1));
	dDELTA2_dDelta_dTau=-(delta-1)*A[i]*(2/beta[i])*pow(powInt(delta-1,2),1/(2*beta[i])-1);
	dDELTA3_dDelta2_dTau=1/(delta-1)*dDELTA2_dDelta_dTau-powInt(delta-1,2)*A[i]*(4/beta[i])*(1/(2*beta[i])-1)*pow(powInt(delta-1,2),1/(2*beta[i])-2);
	
	dDELTAbim1_dTau=(b[i]-1)*pow(DELTA,b[i]-2)*dDELTA_dTau;
	dDELTAbim2_dTau=(b[i]-2)*pow(DELTA,b[i]-3)*dDELTA_dTau;
	Line1=dDELTAbim1_dTau*dDELTA2_dDelta2+pow(DELTA,b[i]-1)*dDELTA3_dDelta2_dTau;
	Line2=(b[i]-1)*(dDELTAbim2_dTau*powInt(dDELTA_dDelta,2)+pow(DELTA,b[i]-2)*2*dDELTA2_dDelta_dTau*dDELTA_dDelta);
	dDELTAbi3_dDelta2_dTau=b[i]*(Line1+Line2);
	
	Line1=pow(DELTA,b[i])*(2*delta*dPSI2_dDelta_dTau+delta*dDELTA3_dDelta2_dTau)+dDELTAbi_dTau*(2*dPSI_dDelta+delta*dPSI2_dDelta2);
	Line2=2*dDELTAbi2_dDelta_dTau*(PSI+delta*dPSI_dDelta)+2*dDELTAbi_dDelta*(dPSI_dTau+delta*dPSI2_dDelta_dTau);
	Line3=dDELTAbi3_dDelta2_dTau*delta*PSI+dDELTAbi2_dDelta2*delta*dPSI_dTau;
        dphir3_dDelta2_dTau+=n[i]*(Line1+Line2+Line3);
    }
    return dphir3_dDelta2_dTau;
}

double phi0_Water(double tau, double delta)
{
    double phi0=0;
    int i;
    
    phi0=log(delta)+n0[1]+n0[2]*tau+n0[3]*log(tau);
    for (i=4;i<=8;i++)
    {
        phi0=phi0+n0[i]*log(1.0-exp(-gamma0[i]*tau));
    }
    return phi0;
}

double dphi0_dDelta_Water(double tau, double delta)
{
    return 1/delta;
}

double dphi02_dDelta2_Water(double tau, double delta)
{
    return -1.0/powInt(delta,2);
}

double dphi0_dTau_Water(double tau, double delta)
{
    double dphi0_dTau=0;
    int i;
    dphi0_dTau=n0[2]+n0[3]/tau;
    for (i=4;i<=8;i++)
    {
        dphi0_dTau=dphi0_dTau+n0[i]*gamma0[i]*(1.0/(1.0-exp(-gamma0[i]*tau))-1.0);
    }
    return dphi0_dTau;
}

double dphi02_dTau2_Water(double tau, double delta)
{
    double dphi02_dTau2=0;
    int i;
    
    dphi02_dTau2=-n0[3]/powInt(tau,2);
    for (i=4;i<=8;i++)
    {
        dphi02_dTau2=dphi02_dTau2-n0[i]*powInt(gamma0[i],2)*exp(-gamma0[i]*tau)/powInt(1.0-exp(-gamma0[i]*tau),2);
    }
    return dphi02_dTau2;
}
