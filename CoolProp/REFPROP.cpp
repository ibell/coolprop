#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#endif

#include "CoolProp.h"

#if defined(__ISWINDOWS__)
#include <windows.h>
#include "REFPROP.h"
#endif

#include <stdlib.h>
#include "string.h"
#include <stdio.h>

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

char LoadedREFPROPRef[2550];
HINSTANCE RefpropdllInstance=NULL;

double REFPROP(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
	int j;
	long i,ierr=0;
	char hf[refpropcharlength*ncmax], hrf[lengthofreference+1],
	herr[errormessagelength+1],hfmix[refpropcharlength+1];
	
	double x[ncmax],xliq[ncmax],xvap[ncmax];
	char RefString[255];
	double T,p=0,d,dl,dv,dl_,dv_,q,e,h,s,cv,cp,w,MW,hl,hv,sl,sv,ul,
		uv,pl,pv,hjt,eta,tcx,Q,Tcrit,pcrit,dcrit,rho,sigma;

	// First create a pointer to an instance of the library
	// Then have windows load the library.
	
	// If REFPROP is not loaded, try to load it
	if (RefpropdllInstance==NULL)
	{
		#if defined(UNICODE)
			RefpropdllInstance = LoadLibrary((LPCSTRW)"refprop.dll");
		#else
			RefpropdllInstance = LoadLibrary((LPCSTR)"refprop.dll");
		#endif
		if (RefpropdllInstance==NULL)
		{
			printf("Could not load REFPROP, not in current location or found on system PATH.  Add location of REFPROP to the PATH environmental variable\n");
			return -_HUGE;
		}
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
		
		// If the fluid name starts with the string "REFPROP-", chop off the "REFPROP-"
		if (!strncmp(Ref,"REFPROP-",8))
		{
		char *REFPROPRef=NULL,*RefCopy=NULL;
		double prop;
			
		// Allocate space for refrigerant name
			RefCopy=(char *)malloc(strlen(Ref)+1);
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

			// Allocate space for refrigerant name
			RefCopy=(char *)malloc(strlen(Ref)+1);
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
		else if (!strcmp(Ref,"Air") || !strcmp(Ref,"R507A") || !strcmp(Ref,"R404A") || !strcmp(Ref,"R410A") || !strcmp(Ref,"R407C"))
		{
			i=1;
			strcpy(RefString,"");
			strcat(RefString,Ref);
			strcat(RefString,".ppf");
			x[0]=1.0;     //Pseudo-Pure fluid
		}
		else if (!strcmp(Ref,"R507A"))
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
			ierr=999;
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
		else if (Output=='I')
		{
			if (Name1=='T'){
				SURFTdll(&Prop1,&dl,x,&sigma,&i,herr,errormessagelength);
				return sigma;
			}
			else{
				std::cout<< "If surface tension is the output, temperature must be the first input" << std::endl;
				return _HUGE;
			}
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
			else if (Output=='P') return p;
			else if (Output=='A') return w;
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
			else if (Output=='A')
			{
				return w;
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
			else if (Output=='D')
			{
				return rho*MW;
			}
			else
				return _HUGE;
		}
		else if (Name1=='T' && Name2=='Q')
		{
			T=Prop1;
			Q=Prop2;
			
			// Saturation Density
			i=1;
			SATTdll(&T,x,&i,&pl,&dl,&dv_,xliq,xvap,&ierr,herr,errormessagelength);
			i=2;
			SATTdll(&T,x,&i,&pv,&dl_,&dv,xliq,xvap,&ierr,herr,errormessagelength);
			if (Output=='D') 
			{
				return 1/(Q/dv+(1-Q)/dl)*MW;
			}
			else if (Output=='P') 
			{
				return (pv*Q+pl*(1-Q));
			}
			else if (Output=='A')
			{
				rho=1/(Q/dv+(1-Q)/dl);
				THERMdll(&T,&rho,x,&p,&e,&h,&s,&cv,&cp,&w,&hjt);
				return w;
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
				p=pv*Q+pl*(1-Q);
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
			else if (Output=='A')
			{
				rho=1/(Q/dv+(1-Q)/dl);
				THERMdll(&T,&rho,x,&p,&e,&h,&s,&cv,&cp,&w,&hjt);
				return w;
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
