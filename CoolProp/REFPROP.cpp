#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#endif

#include <string>
#include "CoolProp.h"

#if defined(__ISWINDOWS__)
#include <windows.h>
#include "REFPROP.h"
#endif

#include <stdlib.h>
#include "string.h"
#include <stdio.h>
#include <iostream>

// Some constants for REFPROP... defined by macros for ease of use 
#define refpropcharlength 255
#define filepathlength 255
#define lengthofreference 3
#define errormessagelength 255
#define ncmax 20		// Note: ncmax is the max number of components
#define numparams 72 
#define maxcoefs 50

#if defined(__ISWINDOWS__)
// For C calling conventions, replaced all "double &" with "double *", and "long &" with "long *"
typedef void (__stdcall *fp_ABFL1dllTYPE)(double *,double *,double *,long *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_ABFL2dllTYPE)(double *,double *,double *,long *,long *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_ACTVYdllTYPE)(double *,double *,double *,double *);
typedef void (__stdcall *fp_AGdllTYPE)(double *,double *,double *,double *,double *);
typedef void (__stdcall *fp_CCRITdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_CP0dllTYPE)(double *,double *,double *);
typedef void (__stdcall *fp_CRITPdllTYPE)(double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_CSATKdllTYPE)(long *,double *,long *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_CV2PKdllTYPE)(long *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_CVCPKdllTYPE)(long *,double *,double *,double *,double *);
typedef void (__stdcall *fp_CVCPdllTYPE)(double *,double *,double *,double *,double *);
typedef void (__stdcall *fp_DBDTdllTYPE)(double *,double *,double *);
typedef void (__stdcall *fp_DBFL1dllTYPE)(double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_DBFL2dllTYPE)(double *,double *,double *,long *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_DDDPdllTYPE)(double *,double *,double *,double *);
typedef void (__stdcall *fp_DDDTdllTYPE)(double *,double *,double *,double *);
typedef void (__stdcall *fp_DEFLSHdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_DHD1dllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *);
typedef void (__stdcall *fp_DHFLSHdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_DIELECdllTYPE)(double *,double *,double *,double *);
typedef void (__stdcall *fp_DOTFILLdllTYPE)(long *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_DPDD2dllTYPE)(double *,double *,double *,double *);
typedef void (__stdcall *fp_DPDDKdllTYPE)(long *,double *,double *,double *);
typedef void (__stdcall *fp_DPDDdllTYPE)(double *,double *,double *,double *);
typedef void (__stdcall *fp_DPDTKdllTYPE)(long *,double *,double *,double *);
typedef void (__stdcall *fp_DPDTdllTYPE)(double *,double *,double *,double *);
typedef void (__stdcall *fp_DPTSATKdllTYPE)(long *,double *,long *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_DSFLSHdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_ENTHALdllTYPE)(double *,double *,double *,double *);
typedef void (__stdcall *fp_ENTROdllTYPE)(double *,double *,double *,double *);
typedef void (__stdcall *fp_ESFLSHdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_FGCTYdllTYPE)(double *,double *,double *,double *);
typedef void (__stdcall *fp_FPVdllTYPE)(double *,double *,double *,double *,double *);
typedef void (__stdcall *fp_GERG04dllTYPE)(long *,long *,long *,char*,long );
typedef void (__stdcall *fp_GETFIJdllTYPE)(char*,double *,char*,char*,long ,long ,long );
typedef void (__stdcall *fp_GETKTVdllTYPE)(long *,long *,char*,double *,char*,char*,char*,char*,long ,long ,long ,long ,long );
typedef void (__stdcall *fp_GIBBSdllTYPE)(double *,double *,double *,double *,double *);
typedef void (__stdcall *fp_HSFLSHdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_INFOdllTYPE)(long *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
typedef void (__stdcall *fp_LIMITKdllTYPE)(char*,long *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long ,long );
typedef void (__stdcall *fp_LIMITSdllTYPE)(char*,double *,double *,double *,double *,double *,long );
typedef void (__stdcall *fp_LIMITXdllTYPE)(char*,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long ,long );
typedef void (__stdcall *fp_MELTPdllTYPE)(double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_MELTTdllTYPE)(double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_MLTH2OdllTYPE)(double *,double *,double *);
typedef void (__stdcall *fp_NAMEdllTYPE)(long *,char*,char*,char*,long ,long ,long );
typedef void (__stdcall *fp_PDFL1dllTYPE)(double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_PDFLSHdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_PEFLSHdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_PHFL1dllTYPE)(double *,double *,double *,long *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_PHFLSHdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_PQFLSHdllTYPE)(double *,double *,double *,long *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_PREOSdllTYPE)(long *);
typedef void (__stdcall *fp_PRESSdllTYPE)(double *,double *,double *,double *);
typedef void (__stdcall *fp_PSFL1dllTYPE)(double *,double *,double *,long *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_PSFLSHdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_PUREFLDdllTYPE)(long *);
typedef void (__stdcall *fp_QMASSdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_QMOLEdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_SATDdllTYPE)(double *,double *,long *,long *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_SATEdllTYPE)(double *,double *,long *,long *,long *,double *,double *,double *,long *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_SATHdllTYPE)(double *,double *,long *,long *,long *,double *,double *,double *,long *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_SATPdllTYPE)(double *,double *,long *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_SATSdllTYPE)(double *,double *,long *,long *,long *,double *,double *,double *,long *,double *,double *,double *,long *,double *,double *,double *,long *,char*,long );
// subroutine SATT (t,x,kph,p,rhol,rhov,xliq,xvap,ierr,herr)
typedef void (__stdcall *fp_SATTdllTYPE)(double *,double *,long *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_SETAGAdllTYPE)(long *,char*,long );
typedef void (__stdcall *fp_SETKTVdllTYPE)(long *,long *,char*,double *,char*,long *,char*,long ,long ,long );
typedef void (__stdcall *fp_SETMIXdllTYPE)(char*,char*,char*,long *,char*,double *,long *,char*,long ,long ,long ,long ,long );
typedef void (__stdcall *fp_SETMODdllTYPE)(long *,char*,char*,char*,long *,char*,long ,long ,long ,long );
typedef void (__stdcall *fp_SETREFdllTYPE)(char*,long *,double *,double *,double *,double *,double *,long *,char*,long ,long );
typedef void (__stdcall *fp_SETUPdllTYPE)(long *,char*,char*,char*,long *,char*,long ,long ,long ,long );
typedef void (__stdcall *fp_SPECGRdllTYPE)(double *,double *,double *,double *);
typedef void (__stdcall *fp_SUBLPdllTYPE)(double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_SUBLTdllTYPE)(double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_SURFTdllTYPE)(double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_SURTENdllTYPE)(double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_TDFLSHdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_TEFLSHdllTYPE)(double *,double *,double *,long *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_THERM0dllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
typedef void (__stdcall *fp_THERM2dllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
typedef void (__stdcall *fp_THERM3dllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
typedef void (__stdcall *fp_THERMdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
typedef void (__stdcall *fp_THFLSHdllTYPE)(double *,double *,double *,long *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_TPFLSHdllTYPE)(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_TPRHOdllTYPE)(double *,double *,double *,long *,long *,double *,long *,char*,long );
typedef void (__stdcall *fp_TQFLSHdllTYPE)(double *,double *,double *,long *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_TRNPRPdllTYPE)(double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_TSFLSHdllTYPE)(double *,double *,double *,long *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,long *,char*,long );
typedef void (__stdcall *fp_VIRBdllTYPE)(double *,double *,double *);
typedef void (__stdcall *fp_VIRCdllTYPE)(double *,double *,double *);
typedef void (__stdcall *fp_WMOLdllTYPE)(double *,double *);
typedef void (__stdcall *fp_XMASSdllTYPE)(double *,double *,double *);
typedef void (__stdcall *fp_XMOLEdllTYPE)(double *,double *,double *);

//Define explicit function pointers
fp_ABFL1dllTYPE ABFL1dll;
fp_ABFL2dllTYPE ABFL2dll;
fp_ACTVYdllTYPE ACTVYdll;
fp_AGdllTYPE AGdll;
fp_CCRITdllTYPE CCRITdll;
fp_CP0dllTYPE CP0dll;
fp_CRITPdllTYPE CRITPdll;
fp_CSATKdllTYPE CSATKdll;
fp_CV2PKdllTYPE CV2PKdll;
fp_CVCPKdllTYPE CVCPKdll;
fp_CVCPdllTYPE CVCPdll;
fp_DBDTdllTYPE DBDTdll;
fp_DBFL1dllTYPE DBFL1dll;
fp_DBFL2dllTYPE DBFL2dll;
fp_DDDPdllTYPE DDDPdll;
fp_DDDTdllTYPE DDDTdll;
fp_DEFLSHdllTYPE DEFLSHdll;
fp_DHD1dllTYPE DHD1dll;
fp_DHFLSHdllTYPE DHFLSHdll;
fp_DIELECdllTYPE DIELECdll;
fp_DOTFILLdllTYPE DOTFILLdll;
fp_DPDD2dllTYPE DPDD2dll;
fp_DPDDKdllTYPE DPDDKdll;
fp_DPDDdllTYPE DPDDdll;
fp_DPDTKdllTYPE DPDTKdll;
fp_DPDTdllTYPE DPDTdll;
fp_DPTSATKdllTYPE DPTSATKdll;
fp_DSFLSHdllTYPE DSFLSHdll;
fp_ENTHALdllTYPE ENTHALdll;
fp_ENTROdllTYPE ENTROdll;
fp_ESFLSHdllTYPE ESFLSHdll;
fp_FGCTYdllTYPE FGCTYdll;
fp_FPVdllTYPE FPVdll;
fp_GERG04dllTYPE GERG04dll;
fp_GETFIJdllTYPE GETFIJdll;
fp_GETKTVdllTYPE GETKTVdll;
fp_GIBBSdllTYPE GIBBSdll;
fp_HSFLSHdllTYPE HSFLSHdll;
fp_INFOdllTYPE INFOdll;
fp_LIMITKdllTYPE LIMITKdll;
fp_LIMITSdllTYPE LIMITSdll;
fp_LIMITXdllTYPE LIMITXdll;
fp_MELTPdllTYPE MELTPdll;
fp_MELTTdllTYPE MELTTdll;
fp_MLTH2OdllTYPE MLTH2Odll;
fp_NAMEdllTYPE NAMEdll;
fp_PDFL1dllTYPE PDFL1dll;
fp_PDFLSHdllTYPE PDFLSHdll;
fp_PEFLSHdllTYPE PEFLSHdll;
fp_PHFL1dllTYPE PHFL1dll;
fp_PHFLSHdllTYPE PHFLSHdll;
fp_PQFLSHdllTYPE PQFLSHdll;
fp_PREOSdllTYPE PREOSdll;
fp_PRESSdllTYPE PRESSdll;
fp_PSFL1dllTYPE PSFL1dll;
fp_PSFLSHdllTYPE PSFLSHdll;
fp_PUREFLDdllTYPE PUREFLDdll;
fp_QMASSdllTYPE QMASSdll;
fp_QMOLEdllTYPE QMOLEdll;
fp_SATDdllTYPE SATDdll;
fp_SATEdllTYPE SATEdll;
fp_SATHdllTYPE SATHdll;
fp_SATPdllTYPE SATPdll;
fp_SATSdllTYPE SATSdll;
fp_SATTdllTYPE SATTdll;
fp_SETAGAdllTYPE SETAGAdll;
fp_SETKTVdllTYPE SETKTVdll;
fp_SETMIXdllTYPE SETMIXdll;
fp_SETMODdllTYPE SETMODdll;
fp_SETREFdllTYPE SETREFdll;
fp_SETUPdllTYPE SETUPdll;
fp_SPECGRdllTYPE SPECGRdll;
fp_SUBLPdllTYPE SUBLPdll;
fp_SUBLTdllTYPE SUBLTdll;
fp_SURFTdllTYPE SURFTdll;
fp_SURTENdllTYPE SURTENdll;
fp_TDFLSHdllTYPE TDFLSHdll;
fp_TEFLSHdllTYPE TEFLSHdll;
fp_THERM0dllTYPE THERM0dll;
fp_THERM2dllTYPE THERM2dll;
fp_THERM3dllTYPE THERM3dll;
fp_THERMdllTYPE THERMdll;
fp_THFLSHdllTYPE THFLSHdll;
fp_TPFLSHdllTYPE TPFLSHdll;
fp_TPRHOdllTYPE TPRHOdll;
fp_TQFLSHdllTYPE TQFLSHdll;
fp_TRNPRPdllTYPE TRNPRPdll;
fp_TSFLSHdllTYPE TSFLSHdll;
fp_VIRBdllTYPE VIRBdll;
fp_VIRCdllTYPE VIRCdll;
fp_WMOLdllTYPE WMOLdll;
fp_XMASSdllTYPE XMASSdll;
fp_XMOLEdllTYPE XMOLEdll;


std::string LoadedREFPROPRef;
HINSTANCE RefpropdllInstance=NULL;
long i;
char hfmix[] = "hmx.bnc";
char hrf[] = "DEF";

double REFPROP(char Output, char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
	return REFPROP(std::string(1,Output),std::string(1,Name1),Prop1,std::string(1,Name2),Prop2,std::string(Ref));
}
double REFPROP(std::string Output, std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref)
{
	long ierr=0,iOutput,iName1,iName2;
	char hf[refpropcharlength*ncmax], herr[errormessagelength+1];
	
	double x[ncmax],xliq[ncmax],xvap[ncmax];
	
	double T,p=0,d,dl,dv,dl_,dv_,q,e,h,s,cv,cp,w,MW,hl,hv,sl,sv,ul,
		uv,pl,pv,hjt,eta,tcx,Q,Tcrit,pcrit,dcrit,rho,sigma;
	std::string sRef;
	std::string RefString;
	// First create a pointer to an instance of the library
	// Then have windows load the library.
	
	// If REFPROP is not loaded
	if (RefpropdllInstance==NULL)
	{
		// Load it
		#if defined(UNICODE)
			RefpropdllInstance = LoadLibrary((LPCWSTR)"refprop.dll");
		#else
			RefpropdllInstance = LoadLibrary((LPCSTR)"refprop.dll");
		#endif
		if (RefpropdllInstance==NULL)
		{
			printf("Could not load REFPROP, not in current location or found on system PATH.  Add location of REFPROP to the PATH environmental variable\n");
			return -_HUGE;
		}

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
	}
	
	
	// If the fluid name does not start with the string "REFPROP-"
	if (Ref.find("REFPROP-") == std::string::npos)
	{
		// Fail and give error
		std::cout << "Invalid REFPROP string: " << Ref << std::endl;
	}
	// Chop off the "REFPROP-"
	else 
	{
		// Keep everything after the "REFPROP-"
		sRef = Ref.substr(8,Ref.size()-8);
	}
	
	// If the name of the refrigerant doesn't match 
	// that of the currently loaded refrigerant
	if (LoadedREFPROPRef.compare(Ref))
	{
		if (!strncmp(sRef.c_str(),"MIX",3))
		{
			// Sample sRef is "MIX:R32[0.697615]&R125[0.302385]" -  this is R410A
			
			// Chop off the MIX by keeping everything after the ':'
			std::string components_joined = strsplit(sRef,':')[1];

			// Split the components_joined into the components
			std::vector<std::string> components_split = strsplit(components_joined,'&');

			// Flush out the refrigerant string for REFPROP
			RefString.clear();

			for (unsigned int j=0;j<components_split.size();j++)
			{	
				// Get component name and mole fraction (as strings)
				std::vector<std::string> comp_fraction = strsplit(components_split[j],'[');
				
				// Build the refrigerant string
				if (j == 0){
					RefString = comp_fraction[0]+".fld";
				}
				else{
					RefString += "|"+comp_fraction[0]+".fld";
				}
				// Convert the mole fraction (as string) to a number
				x[j] = strtod(comp_fraction[1].c_str(),NULL);

				// Update the number of components
				i = j+1;
			}
		}
		else if (!sRef.compare("Air") || !sRef.compare("R507A") || !sRef.compare("R404A") || !sRef.compare("R410A") || !sRef.compare("R407C") || !sRef.compare("SES36"))
		{
			i=1;
			RefString = std::string(sRef)+std::string(".ppf");
			x[0]=1.0;     //Pseudo-Pure fluid
		}
		/*else if (!strcmp(Ref,"R507A"))
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
		}*/
		else
		{
			i=1;
			RefString = std::string(sRef)+std::string(".fld");
			x[0]=1.0;     //Pure fluid
		}

		strcpy(hf,RefString.c_str());

		ierr=999;
		//...Call SETUP to initialize the program
		SETUPdll(&i, hf, hfmix, hrf, &ierr, herr,
			refpropcharlength*ncmax,refpropcharlength,
			lengthofreference,errormessagelength);
		if (ierr != 0) printf("REFPROP setup gives this error during SETUP: %s\n",herr);
		//Copy the name of the loaded refrigerant back into the temporary holder
		LoadedREFPROPRef = std::string(Ref);
	}
	
	strcpy(herr,"Ok");
	
	iOutput = get_param_index(Output);
	iName1 = get_param_index(Name1);
	iName2 = get_param_index(Name2);
	
	// Get the molar mass of the fluid
	WMOLdll(x,&MW);	
	if (iOutput == iTcrit)
	{
		// Critical temperature
		CRITPdll(x,&Tcrit,&pcrit,&dcrit,&ierr,herr,255);
		return Tcrit;
	}
	else if (iOutput==iMM)
	{
		// mole mass
		return MW;
	}
	else if (iOutput==iPcrit)
	{
		// Critical pressure
		CRITPdll(x,&Tcrit,&pcrit,&dcrit,&ierr,herr,255);
		return pcrit;
	}
	else if (iOutput ==iRhocrit)
	{
		// Critical density
		CRITPdll(x,&Tcrit,&pcrit,&dcrit,&ierr,herr,255);
		return dcrit*MW;
		
	}
	else if (iOutput == iTmin)
	{
		// Minimum temperature
		double tmin,tmax,Dmax,pmax;
		LIMITSdll("EOS",x,&tmin,&tmax,&Dmax,&pmax,255);
		return tmin;
	}
	else if (iOutput == iAccentric)
	{
		double wmm,Ttriple,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas;
		// Accentric factor
		if (i>1)
		{
			fprintf(stderr,"Error: Accentric factor only defined for pure fluids\n");
			return _HUGE;
		}
		INFOdll(&i,&wmm,&Ttriple,&tnbpt,&tc,&pc,&Dc,&Zc,&acf,&dip,&Rgas);
		return acf;
	}
	else if (iOutput ==iDipole)
	{
		double wmm,Ttriple,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas;
		// Dipole moment
		if (i>1)
		{
			fprintf(stderr,"Error: Dipole moment only defined for pure fluids\n");
			return _HUGE;
		}
		INFOdll(&i,&wmm,&Ttriple,&tnbpt,&tc,&pc,&Dc,&Zc,&acf,&dip,&Rgas);
		return dip;
	}
	else if (iOutput==iTtriple)
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
		INFOdll(&i,&wmm,&Ttriple,&tnbpt,&tc,&pc,&Dc,&Zc,&acf,&dip,&Rgas);
		return Ttriple;
	}
	else if (iOutput==iI)
	{
		if (iName1==iT){
			SURFTdll(&Prop1,&dl,x,&sigma,&i,herr,errormessagelength);
			return sigma;
		}
		else{
			std::cout<< "If surface tension is the output, temperature must be the first input" << std::endl;
			return _HUGE;
		}
	}
	
	else if ((iName1==iT && iName2 == iP) || (iName2==  iT && iName1== iP))
	{
		// T in K, P in kPa
		if (iName1 == iP){
			std::swap(Prop1,Prop2);
		}
		T = Prop1; p = Prop2;

		// Use flash routine to find properties
		TPFLSHdll(&T,&p,x,&d,&dl,&dv,xliq,xvap,&q,&e,&h,&s,&cv,&cp,&w,&ierr,herr,errormessagelength);
		if (iOutput==iH) return h/MW;
		else if (iOutput==iD) return d*MW;
		else if (iOutput==iS) return s/MW;
		else if (iOutput==iU) return e/MW;
		else if (iOutput==iC) return cp/MW;
		else if (iOutput==iO) return cv/MW;
		else if (iOutput==iP) return p;
		else if (iOutput==iA) return w;
		else if (iOutput==iV) 
		{
			TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
			return eta/1.0e6; //uPa-s to Pa-s
		} 
		else if (iOutput==iL)
		{
			TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
			return tcx/1000.0; //W/m-K to kW/m-K
		}
		else
			return _HUGE;
	}
	else if ((iName1==iT && iName2==iD) || (iName2==iT && iName1==iD))
	{
		// T in K, D in kg/m^3
		if (iName2 == iT){
			std::swap(Prop1,Prop2);
		}
		T = Prop1; rho = Prop2/MW;
		
		// This is the explicit formulation of the EOS
		TDFLSHdll(&T,&rho,x,&p,&dl,&dv,xliq,xvap,&q,&e,&h,&s,&cv,&cp,&w,&ierr,herr,errormessagelength);

		if (iOutput==iP)
		{
			return p;
		}
		if (iOutput==iH)
		{
			return h/MW;
		}
		else if (iOutput==iA)
		{
			return w;
		}
		else if (iOutput==iS)
		{
			return s/MW;
		}
		else if (iOutput==iU)
		{
			return (h-p/rho)/MW;
		}
		else if (iOutput==iC)
		{
			return cp/MW;
		}
		else if (iOutput==iO)
		{
			return cv/MW;
		}
		else if (iOutput==iV) 
		{
			TRNPRPdll(&T,&rho,x,&eta,&tcx,&ierr,herr,errormessagelength);
			return eta/1.0e6; //uPa-s to Pa-s
		} 
		else if (iOutput==iL)
		{
			TRNPRPdll(&T,&rho,x,&eta,&tcx,&ierr,herr,errormessagelength);
			return tcx/1000.0; //W/m-K to kW/m-K
		}
		else if (iOutput==iD)
		{
			return rho*MW;
		}
		else
			return _HUGE;
	}
	else if ((iName1==iT && iName2==iQ) || (iName2==iT && iName1==iQ))
	{

		long ic;
		if (iName2 == iT){
			std::swap(Prop1,Prop2);
		}
		T = Prop1; Q = Prop2;
		
		// Saturation Density
		ic=1;
		SATTdll(&T,x,&ic,&pl,&dl,&dv_,xliq,xvap,&ierr,herr,errormessagelength);
		ic=2;
		SATTdll(&T,x,&ic,&pv,&dl_,&dv,xliq,xvap,&ierr,herr,errormessagelength);
		if (iOutput==iD) 
		{
			return 1/(Q/dv+(1-Q)/dl)*MW;
		}
		else if (iOutput==iP) 
		{
			return (pv*Q+pl*(1-Q));
		}
		else if (iOutput==iA)
		{
			rho=1/(Q/dv+(1-Q)/dl);
			THERMdll(&T,&rho,x,&p,&e,&h,&s,&cv,&cp,&w,&hjt);
			return w;
		}
		else if (iOutput==iH) 
		{
			ENTHALdll(&T,&dl,xliq,&hl);
			ENTHALdll(&T,&dv,xvap,&hv);
			return (hv*Q+hl*(1-Q))/MW; // J/kg to kJ/kg
		}
		else if (iOutput==iS) 
		{
			ENTROdll(&T,&dl,xliq,&sl);
			ENTROdll(&T,&dv,xvap,&sv);
			return (sv*Q+sl*(1-Q))/MW; // J/kg-K to kJ/kg-K
		}
		else if (iOutput==iU) 
		{
			ENTHALdll(&T,&dl,xliq,&hl);
			ENTHALdll(&T,&dv,xvap,&hv);
			p=pv*Q+pl*(1-Q);
			ul=hl-p/dl;
			uv=hv-p/dv;
			return (uv*Q+ul*(1-Q))/MW; // J/kg to kJ/kg
		}
		else if (iOutput==iC) 
		{
			d=1/(Q/dv+(1-Q)/dl);
			CVCPdll(&T,&d,x,&cv,&cp);
			return cp/MW; // J/kg-K to kJ/kg-K
		}
		else if (iOutput==iO) 
		{
			d=1/(Q/dv+(1-Q)/dl);
			CVCPdll(&T,&d,x,&cv,&cp);
			return cv/MW; // J/kg-K to kJ/kg-K
		}
		else if (iOutput==iV) 
		{
			d=1/(Q/dv+(1-Q)/dl);
			TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
			return eta/1.0e6; //uPa-s to Pa-s
		}
		else if (iOutput==iL) 
		{
			d=1/(Q/dv+(1-Q)/dl);
			TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
			return tcx/1000.0; //W/m-K to kW/m-K
		}
		else
			return _HUGE;
	}
	else if ((iName1==iP && iName2==iQ) || (iName2==iP && iName1==iQ))
	{
		if (iName2 == iP){
			std::swap(Prop1,Prop2);
		}
		p = Prop1; Q = Prop2;

		double dummy;
		
		// Saturation Density for the liquid
		long kph = 1;
		SATPdll(&p,x,&kph,&T,&dl,&dummy,xliq,xvap,&ierr,herr,errormessagelength);
		// Saturation density for the vapor
		kph = 2;
		SATPdll(&p,x,&kph,&T,&dummy,&dv,xliq,xvap,&ierr,herr,errormessagelength);
		if (iOutput==iT)
		{
			return T;
		}
		else if (iOutput==iD) 
		{
			return 1/(Q/dv+(1-Q)/dl)*MW;
		}
		else if (iOutput==iP) 
		{
			PRESSdll(&T,&dl,xliq,&pl);
			PRESSdll(&T,&dv,xvap,&pv);
			return (pv*Q+pl*(1-Q));
		}
		else if (iOutput==iH) 
		{
			ENTHALdll(&T,&dl,xliq,&hl);
			ENTHALdll(&T,&dv,xvap,&hv);
			return (hv*Q+hl*(1-Q))/MW; // J/kg to kJ/kg
		}
		else if (iOutput==iS) 
		{
			ENTROdll(&T,&dl,xliq,&sl);
			ENTROdll(&T,&dv,xvap,&sv);
			return (sv*Q+sl*(1-Q))/MW; // J/kg-K to kJ/kg-K
		}
		else if (iOutput==iU) 
		{
			ENTHALdll(&T,&dl,xliq,&hl);
			ENTHALdll(&T,&dv,xvap,&hv);
			ul=hl-p/dl;
			uv=hv-p/dv;
			return (uv*Q+ul*(1-Q))/MW; // J/kg to kJ/kg
		}
		else if (iOutput==iC) 
		{
			d=1/(Q/dv+(1-Q)/dl);
			CVCPdll(&T,&d,x,&cv,&cp);
			return cp/MW; // J/kg-K to kJ/kg-K
		}
		else if (iOutput==iO) 
		{
			d=1/(Q/dv+(1-Q)/dl);
			CVCPdll(&T,&d,x,&cv,&cp);
			return cv/MW; // J/kg-K to kJ/kg-K
		}
		else if (iOutput==iA)
		{
			rho=1/(Q/dv+(1-Q)/dl);
			THERMdll(&T,&rho,x,&p,&e,&h,&s,&cv,&cp,&w,&hjt);
			return w;
		}
		else if (iOutput==iV) 
		{
			d=1/(Q/dv+(1-Q)/dl);
			TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
			return eta/1.0e6; //uPa-s to Pa-s
		}
		else if (iOutput==iL) 
		{
			d=1/(Q/dv+(1-Q)/dl);
			TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
			return tcx/1000.0; //W/m-K to kW/m-K
		}
		else
			return _HUGE;
	}
	else if ((iName1==iP && iName2==iH) || (iName2==iP && iName1==iH))
	{
		// p in kPa, h in kJ/kg
		if (iName2 == iP){
			std::swap(Prop1,Prop2);
		}
		p = Prop1; h = Prop2*MW;
		
		// Use flash routine to find properties
		PHFLSHdll(&p,&h,x,&T,&d,&dl,&dv,xliq,xvap,&q,&e,&s,&cv,&cp,&w,&ierr,herr,errormessagelength);
		if (iOutput==iH) return h/MW;
		else if (iOutput==iT) return T;
		else if (iOutput==iD) return d*MW;
		else if (iOutput==iS) return s/MW;
		else if (iOutput==iU) return e/MW;
		else if (iOutput==iC) return cp/MW;
		else if (iOutput==iO) return cv/MW;
		else if (iOutput==iP) return p;
		else if (iOutput==iA) return w;
		else if (iOutput==iV) 
		{
			TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
			return eta/1.0e6; //uPa-s to Pa-s
		} 
		else if (iOutput==iL)
		{
			TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
			return tcx/1000.0; //W/m-K to kW/m-K
		}
		else
			return _HUGE;
	}
	else if ((iName1==iP && iName2==iS) || (iName2==iP && iName1==iS))
	{
		// p in kPa, h in kJ/kg
		if (iName2 == iP){
			std::swap(Prop1,Prop2);
		}
		p = Prop1; s = Prop2*MW;
		
		// Use flash routine to find properties
		PSFLSHdll(&p,&s,x,&T,&d,&dl,&dv,xliq,xvap,&q,&e,&h,&cv,&cp,&w,&ierr,herr,errormessagelength);
		if (iOutput==iH) return h/MW;
		else if (iOutput==iT) return T;
		else if (iOutput==iD) return d*MW;
		else if (iOutput==iS) return s/MW;
		else if (iOutput==iU) return e/MW;
		else if (iOutput==iC) return cp/MW;
		else if (iOutput==iO) return cv/MW;
		else if (iOutput==iP) return p;
		else if (iOutput==iA) return w;
		else if (iOutput==iV) 
		{
			TRNPRPdll(&T,&d,x,&eta,&tcx,&ierr,herr,errormessagelength);
			return eta/1.0e6; //uPa-s to Pa-s
		} 
		else if (iOutput==iL)
		{
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
