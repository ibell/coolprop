#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#endif
#include <stdlib.h>

#include "math.h"
#include "stdio.h"
#include <string.h>
#include "PropMacros.h"
#include "CoolPropTools.h"
#include "CoolProp.h"
#include "R410A.h"

#define nT 200
#define nP 200
#define nLUT 5

char RefLUT[255][nLUT]; // The names of the fluids that are loaded in LUT
static double hmat[nT][nP][nLUT];
static double rhomat[nT][nP][nLUT];
static double cpmat[nT][nP][nLUT];
static double smat[nT][nP][nLUT];
static double cvmat[nT][nP][nLUT];
static double umat[nT][nP][nLUT];
static double viscmat[nT][nP][nLUT];
static double kmat[nT][nP][nLUT];
static double Tvec[nT][nLUT];
static double pvec[nP][nLUT];

static int isNAN(double x)
{
	// recommendation from http://www.devx.com/tips/Tip/42853
	return x != x;
}

static int isINFINITY(double x)
{
	// recommendation from http://www.devx.com/tips/Tip/42853
	if ((x == x) && ((x - x) != 0.0))
		return 1;//return (x < 0.0 ? -1 : 1); // This will tell you whether positive or negative infinity
	else
		return 0;
}

int ValidNumber(double x)
{
	if (!isNAN(x) && !isINFINITY(x))
		return 1;
	else
		return 0;
}
double powInt(double x, int y)
{
    // Raise a double to an integer power
    // Overload not provided in math.h
    int i;
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
    for (i=1;i<y_in;i++)
    {
        product=product*x_in;
    }
    return product;
}

double QuadInterp(double x0, double x1, double x2, double f0, double f1, double f2, double x)
{
    /* Quadratic interpolation.  
    Based on method from Kreyszig, 
    Advanced Engineering Mathematics, 9th Edition 
    */
    double L0, L1, L2;
    L0=((x-x1)*(x-x2))/((x0-x1)*(x0-x2));
    L1=((x-x0)*(x-x2))/((x1-x0)*(x1-x2));
    L2=((x-x0)*(x-x1))/((x2-x0)*(x2-x1));
    return L0*f0+L1*f1+L2*f2;
}

int WriteLookup2File(int ILUT)
{
    int i,j;
    FILE *fp_h,*fp_s,*fp_rho,*fp_u,*fp_cp,*fp_cv,*fp_visc;
    fp_h=fopen("h.csv","w");
    fp_s=fopen("s.csv","w");
    fp_u=fopen("u.csv","w");
    fp_cp=fopen("cp.csv","w");
    fp_cv=fopen("cv.csv","w");
    fp_rho=fopen("rho.csv","w");
    fp_visc=fopen("visc.csv","w");

    // Write the pressure header row
    for (j=0;j<nP;j++)
    {
        fprintf(fp_h,",%0.12f",pvec[j][ILUT]);
        fprintf(fp_s,",%0.12f",pvec[j][ILUT]);
        fprintf(fp_rho,",%0.12f",pvec[j][ILUT]);
        fprintf(fp_u,",%0.12f",pvec[j][ILUT]);
        fprintf(fp_cp,",%0.12f",pvec[j][ILUT]);
        fprintf(fp_cv,",%0.12f",pvec[j][ILUT]);
        fprintf(fp_visc,",%0.12f",pvec[j][ILUT]);
    }
    fprintf(fp_h,"\n");
    fprintf(fp_s,"\n");
    fprintf(fp_rho,"\n");
    fprintf(fp_u,"\n");
    fprintf(fp_cp,"\n");
    fprintf(fp_cv,"\n");
    fprintf(fp_visc,"\n");
    
    for (i=1;i<nT;i++)
    {
        fprintf(fp_h,"%0.12f",Tvec[i][ILUT]);
        fprintf(fp_s,"%0.12f",Tvec[i][ILUT]);
        fprintf(fp_rho,"%0.12f",Tvec[i][ILUT]);
        fprintf(fp_u,"%0.12f",Tvec[i][ILUT]);
        fprintf(fp_cp,"%0.12f",Tvec[i][ILUT]);
        fprintf(fp_cv,"%0.12f",Tvec[i][ILUT]);
        fprintf(fp_visc,"%0.12f",Tvec[i][ILUT]);
        for (j=0;j<nP;j++)
        {
            fprintf(fp_h,",%0.12f",hmat[i][j][ILUT]);
            fprintf(fp_s,",%0.12f",smat[i][j][ILUT]);
            fprintf(fp_rho,",%0.12f",rhomat[i][j][ILUT]);
            fprintf(fp_u,",%0.12f",umat[i][j][ILUT]);
            fprintf(fp_cp,",%0.12f",cpmat[i][j][ILUT]);
            fprintf(fp_cv,",%0.12f",cvmat[i][j][ILUT]);
            fprintf(fp_visc,",%0.12f",viscmat[i][j][ILUT]);
        }
        fprintf(fp_h,"\n");
        fprintf(fp_s,"\n");
        fprintf(fp_rho,"\n");
        fprintf(fp_u,"\n");
        fprintf(fp_cp,"\n");
        fprintf(fp_cv,"\n");
        fprintf(fp_visc,"\n");
    }
    return 1;
}

int BuildLookupTable(char *Ref, struct fluidParamsVals *Fluid)
{
    int i,j,k,Tmin,Tmax,pmin,pmax,OldUseLUT,test1,test2,test3;
    double Tc;
    double (*p_dp)(double);
    double (*p_bp)(double);
    
    Tc=Fluid->Tc;
    OldUseLUT=SinglePhaseLUTStatus();
    UseSinglePhaseLUT(0);
    
    if (Fluid->Type==FLUIDTYPE_REFRIGERANT_PURE)
    {
        p_dp=Fluid->funcs.psat;
        p_bp=Fluid->funcs.psat;
    }
    else if (Fluid->Type==FLUIDTYPE_REFRIGERANT_PSEUDOPURE)
    {
        p_dp=Fluid->funcs.p_dp;
        p_bp=Fluid->funcs.p_bp;
    }
    else
    {
        printf("Invalid fluid type for building LUT");
    }
    
    // Properties evaluated at all points with X in the 
    // following p-h plot:
    
    /*      Supercritical
    ||X X X X X X X X X X X X X X X X X X X X X X
    ||X X X X X X X X X X X X X X X X X X X X X X
    ||	--------  X X X X X X X X X X X X
    ||    /		 \  X X X X X X X X X X X
p	||   /		  | X X X X Superheated Gas
    ||  /	Two   / X X X X X X X X X X X
    || /  Phase	 /X X X X X X X X X X X X 
    ||/		      / X X X X X X X X X X X X 
    ||=	= = = = = = = = = = = = = = = = = = = = =
                        Enthalpy										
    */
    
    // Check if the LUT is already in memory
    for (k=0;k<nLUT;k++)
    {
        // If it is found, break out of loop and stop; no LUT needed
        if (!strcmp(RefLUT[k],Ref)) return k;
        // You made it to an empty string, need a LUT, jump out of for loop
        if (!strcmp(RefLUT[k],"")) break;
    }
    if (k==nLUT)
    {
        // Ran out of fluids
        printf("Sorry, ran out of saturation LUT slots, increase nLUT in CoolPropTools.h and recompile\n");
        return -1;
    }
        
    // Therefore it hasn't been found yet, assign the refrigerant name
    strcpy(RefLUT[k],Ref);
    
    Tmin=Fluid->LUT.Tmin;
    Tmax=Fluid->LUT.Tmax;
    pmin=Fluid->LUT.pmin;
    pmax=Fluid->LUT.pmax;

	//Warn the user
    printf("Building Lookup Tables... Please wait...");
    for (i=0;i<nT;i++)
    {
        Tvec[i][k]=Tmin+i*(Tmax-Tmin)/(nT-1);
    }
    for (j=0;j<nP;j++)
    {
        pvec[j][k]=pmin+j*(pmax-pmin)/(nP-1);
    }
    for (i=0;i<nT;i++)
    {
        for (j=0;j<nP;j++)
        {
        	test1=Tvec[i][k]>Tc;
        	test2=(Fluid->Type==FLUIDTYPE_REFRIGERANT_PSEUDOPURE && (pvec[j][k]>p_dp(Tvec[i][k]) || pvec[j][k]<p_bp(Tvec[i][k])));
        	test3=(Fluid->Type==FLUIDTYPE_REFRIGERANT_PURE && fabs((pvec[j][k]-Fluid->funcs.psat(Tvec[i][k]))/pvec[j][k]-1)>0.001);

            if (test1 || test2 || test3)
            {					
                rhomat[i][j][k]=Props('D','T',Tvec[i][k],'P',pvec[j][k],Ref);
                hmat[i][j][k]=Props('H','T',Tvec[i][k],'D',rhomat[i][j][k],Ref);
                smat[i][j][k]=Props('S','T',Tvec[i][k],'D',rhomat[i][j][k],Ref);
                umat[i][j][k]=Props('U','T',Tvec[i][k],'D',rhomat[i][j][k],Ref);
                cpmat[i][j][k]=Props('C','T',Tvec[i][k],'D',rhomat[i][j][k],Ref);
                cvmat[i][j][k]=Props('O','T',Tvec[i][k],'D',rhomat[i][j][k],Ref);
                viscmat[i][j][k]=Props('V','T',Tvec[i][k],'D',rhomat[i][j][k],Ref);
                kmat[i][j][k]=Props('L','T',Tvec[i][k],'D',rhomat[i][j][k],Ref);
            }
            else
            {
                rhomat[i][j][k]=_HUGE;
                hmat[i][j][k]=_HUGE;
                smat[i][j][k]=_HUGE;
                umat[i][j][k]=_HUGE;
                cpmat[i][j][k]=_HUGE;
                cvmat[i][j][k]=_HUGE;
                viscmat[i][j][k]=_HUGE;
                kmat[i][j][k]=_HUGE;
            }
        }
    }
    UseSinglePhaseLUT(OldUseLUT);
    printf("Tables Built\n");
    return k;
}

double LookupValue(char Prop, double T, double p, char *Ref, struct fluidParamsVals *Fluid)
{
    int iPlow, iPhigh, iTlow, iThigh,L,R,M,ILUT,success;
    double T1, T2, T3, P1, P2, P3, y1, y2, y3, a1, a2, a3;
    double (*mat)[nT][nP][nLUT];
    double Tmin,Tmax,Pmin,Pmax;

    // Either build and get the index of, or just build the LUT
    ILUT=BuildLookupTable(Ref,Fluid);
    
    Tmin=Tvec[0][ILUT];
    Tmax=Tvec[nT-1][ILUT];
    Pmin=pvec[0][ILUT];
    Pmax=pvec[nP-1][ILUT];
    
    if (T>Tmax || T<Tmin || p>Pmax ||p<Pmin)
    {
        success=sprintf(CP_errString,"Input to LookupValue() for %c is out of bounds [T:%g p:%g]",Prop,T,p);
        if (success<0) printf("error writing error string");
        return -1e6;
    }

    L=0;
    R=nP-1;
    M=(L+R)/2;
    // Use interval halving to find the indices which bracket the temperature of interest
    while (R-L>1)
    {
        if (T>=Tvec[M][ILUT])
        { L=M; M=(L+R)/2; }
        if (T<Tvec[M][ILUT])
        { R=M; M=(L+R)/2; }
    }
    iTlow=L; iThigh=R;

    L=0;
    R=nP-1;
    M=(L+R)/2;
    // Use interval halving to find the indices which bracket the pressure of interest
    while (R-L>1)
    {
        if (p>=pvec[M][ILUT])
        { L=M; M=(L+R)/2; }
        if (p<pvec[M][ILUT])
        { R=M; M=(L+R)/2; }
    }
    iPlow=L; iPhigh=R;

    /* Depending on which property is desired, 
    make the matrix "mat" a pointer to the 
    desired property matrix */
    if (Prop=='D')
        mat=&rhomat;
    else if (Prop=='C')
        mat=&cpmat;
    else if (Prop=='O')
        mat=&cvmat;
    else if (Prop=='H')
        mat=&hmat;
    else if (Prop=='S')
        mat=&smat;
    else if (Prop=='U')
        mat=&umat;
    else if (Prop=='V')
        mat=&viscmat;
    else if (Prop=='L')
        mat=&kmat;
    else
    {
    	//ERROR
    	printf("Invalid output type [%c]",Prop);
    	return -1;
    }
    
    //At Low Temperature Index
    y1=(*mat)[iTlow][iPlow][ILUT];
    y2=(*mat)[iTlow][iPhigh][ILUT];
    y3=(*mat)[iTlow][iPhigh+1][ILUT];
    P1=pvec[iPlow][ILUT];
    P2=pvec[iPhigh][ILUT];
    P3=pvec[iPhigh+1][ILUT];
    a1=QuadInterp(P1,P2,P3,y1,y2,y3,p);

    //At High Temperature Index
    y1=(*mat)[iThigh][iPlow][ILUT];
    y2=(*mat)[iThigh][iPhigh][ILUT];
    y3=(*mat)[iThigh][iPhigh+1][ILUT];
    a2=QuadInterp(P1,P2,P3,y1,y2,y3,p);

    //At High Temperature Index+1 (for QuadInterp() )
    y1=(*mat)[iThigh+1][iPlow][ILUT];
    y2=(*mat)[iThigh+1][iPhigh][ILUT];
    y3=(*mat)[iThigh+1][iPhigh+1][ILUT];
    a3=QuadInterp(P1,P2,P3,y1,y2,y3,p);

    //At Final Interpolation
    T1=Tvec[iTlow][ILUT];
    T2=Tvec[iThigh][ILUT];
    T3=Tvec[iThigh+1][ILUT];
    
    return QuadInterp(T1,T2,T3,a1,a2,a3,T);	
}

void Append2ErrorString(char *string)
{
	ErrorFlag=FAIL;
	strcat(CP_errString,"||");
	strcat(CP_errString,string);
}
void PrintError(void)
{
	fprintf(stderr,"CoolProp error returned: %s\n",CP_errString);
}
