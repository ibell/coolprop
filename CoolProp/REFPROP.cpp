#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#include <sys/stat.h>
#else
#include <sys/stat.h>
#endif

#include <string>
#include "CoolProp.h"

#if defined(__ISWINDOWS__)
#include <windows.h>
#elif defined(__ISLINUX__)
#include <dlfcn.h>
#elif defined(__ISAPPLE__)
#include <dlfcn.h>
#endif

#include "REFPROP_lib.h"
#include "REFPROP.h"
#include "CoolPropTools.h"

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

// Check windows
#if _WIN32 || _WIN64
   #if _WIN64
     #define ENV64BIT
  #else
    #define ENV32BIT
  #endif
#endif

// Check GCC
#if __GNUC__
  #if __x86_64__ || __ppc64__
    #define ENV64BIT
  #else
    #define ENV32BIT
  #endif
#endif

std::vector<double> x(ncmax,0), LoadedREFPROPx(ncmax,0);

std::string LoadedREFPROPRef;

#if defined(__ISWINDOWS__)
HINSTANCE RefpropdllInstance=NULL;
#elif defined(__ISLINUX__)
void *RefpropdllInstance=NULL;
#elif defined(__ISAPPLE__)
void *RefpropdllInstance=NULL;
#else
void *RefpropdllInstance=NULL;
#endif

// Define functions as pointers and initialise them to NULL
// Declare the functions for direct access
 RPVersion_POINTER RPVersion;
 SETPATHdll_POINTER SETPATHdll;
 ABFL1dll_POINTER ABFL1dll;
 ABFL2dll_POINTER ABFL2dll;
 ACTVYdll_POINTER ACTVYdll;
 AGdll_POINTER AGdll;
 CCRITdll_POINTER CCRITdll;
 CP0dll_POINTER CP0dll;
 CRITPdll_POINTER CRITPdll;
 CSATKdll_POINTER CSATKdll;
 CV2PKdll_POINTER CV2PKdll;
 CVCPKdll_POINTER CVCPKdll;
 CVCPdll_POINTER CVCPdll;
 DBDTdll_POINTER DBDTdll;
 DBFL1dll_POINTER DBFL1dll;
 DBFL2dll_POINTER DBFL2dll;
 DDDPdll_POINTER DDDPdll;
 DDDTdll_POINTER DDDTdll;
 DEFLSHdll_POINTER DEFLSHdll;
 DHD1dll_POINTER DHD1dll;
 DHFLSHdll_POINTER DHFLSHdll;
 DHFL1dll_POINTER DHFL1dll;
 DHFL2dll_POINTER DHFL2dll;
 DIELECdll_POINTER DIELECdll;
 DOTFILLdll_POINTER DOTFILLdll;
 DPDD2dll_POINTER DPDD2dll;
 DPDDKdll_POINTER DPDDKdll;
 DPDDdll_POINTER DPDDdll;
 DPDTKdll_POINTER DPDTKdll;
 DPDTdll_POINTER DPDTdll;
 DPTSATKdll_POINTER DPTSATKdll;
 DSFLSHdll_POINTER DSFLSHdll;
 DSFL1dll_POINTER DSFL1dll;
 DSFL2dll_POINTER DSFL2dll;
 ENTHALdll_POINTER ENTHALdll;
 ENTROdll_POINTER ENTROdll;
 ESFLSHdll_POINTER ESFLSHdll;
 FGCTYdll_POINTER FGCTYdll;
 FPVdll_POINTER FPVdll;
 GERG04dll_POINTER GERG04dll;
 GETFIJdll_POINTER GETFIJdll;
 GETKTVdll_POINTER GETKTVdll;
 GIBBSdll_POINTER GIBBSdll;
 HSFLSHdll_POINTER HSFLSHdll;
 INFOdll_POINTER INFOdll;
 LIMITKdll_POINTER LIMITKdll;
 LIMITSdll_POINTER LIMITSdll;
 LIMITXdll_POINTER LIMITXdll;
 MELTPdll_POINTER MELTPdll;
 MELTTdll_POINTER MELTTdll;
 MLTH2Odll_POINTER MLTH2Odll;
 NAMEdll_POINTER NAMEdll;
 PDFL1dll_POINTER PDFL1dll;
 PDFLSHdll_POINTER PDFLSHdll;
 PEFLSHdll_POINTER PEFLSHdll;
 PHFL1dll_POINTER PHFL1dll;
 PHFLSHdll_POINTER PHFLSHdll;
 PQFLSHdll_POINTER PQFLSHdll;
 PREOSdll_POINTER PREOSdll;
 PRESSdll_POINTER PRESSdll;
 PSFL1dll_POINTER PSFL1dll;
 PSFLSHdll_POINTER PSFLSHdll;
 PUREFLDdll_POINTER PUREFLDdll;
 QMASSdll_POINTER QMASSdll;
 QMOLEdll_POINTER QMOLEdll;
 RESIDUALdll_POINTER RESIDUALdll;
 SATDdll_POINTER SATDdll;
 SATEdll_POINTER SATEdll;
 SATHdll_POINTER SATHdll;
 SATPdll_POINTER SATPdll;
 SATSdll_POINTER SATSdll;
 SATTdll_POINTER SATTdll;
 SETAGAdll_POINTER SETAGAdll;
 SETKTVdll_POINTER SETKTVdll;
 SETMIXdll_POINTER SETMIXdll;
 SETMODdll_POINTER SETMODdll;
 SETREFdll_POINTER SETREFdll;
 SETUPdll_POINTER SETUPdll;
//  SPECGRdll_POINTER SPECGRdll; // not found in library
 SUBLPdll_POINTER SUBLPdll;
 SUBLTdll_POINTER SUBLTdll;
 SURFTdll_POINTER SURFTdll;
 SURTENdll_POINTER SURTENdll;
 TDFLSHdll_POINTER TDFLSHdll;
 TEFLSHdll_POINTER TEFLSHdll;
 THERM0dll_POINTER THERM0dll;
 THERM2dll_POINTER THERM2dll;
 THERM3dll_POINTER THERM3dll;
 THERMdll_POINTER THERMdll;
 THFLSHdll_POINTER THFLSHdll;
 TPFLSHdll_POINTER TPFLSHdll;
 TPFL2dll_POINTER TPFL2dll;
 TPRHOdll_POINTER TPRHOdll;
 TQFLSHdll_POINTER TQFLSHdll;
 TRNPRPdll_POINTER TRNPRPdll;
 TSFLSHdll_POINTER TSFLSHdll;
 VIRBdll_POINTER VIRBdll;
 VIRCdll_POINTER VIRCdll;
 WMOLdll_POINTER WMOLdll;
 XMASSdll_POINTER XMASSdll;
 XMOLEdll_POINTER XMOLEdll;

void *getFunctionPointer(char * name)
{
	#if defined(__ISWINDOWS__)
		return (void *) GetProcAddress(RefpropdllInstance,name);
	#elif defined(__ISLINUX__)
		return dlsym(RefpropdllInstance,name);
	#elif defined(__ISAPPLE__)
		return dlsym(RefpropdllInstance,name);
	#else
		throw NotImplementedError("This function should not be called.");
		return NULL;
    #endif
}

//Moved pointer handling to a function, helps to maintain
//an overview and structures OS dependent parts
double setFunctionPointers()
{
	if (RefpropdllInstance==NULL)
	{
		printf("REFPROP is not loaded, make sure you call this function after loading the library.\n");
		return -_HUGE;
	}
	// set the pointers, platform independent
	RPVersion = (RPVersion_POINTER) getFunctionPointer((char *)RPVersion_NAME);
	ABFL1dll = (ABFL1dll_POINTER) getFunctionPointer((char *)ABFL1dll_NAME);
	ABFL2dll = (ABFL2dll_POINTER) getFunctionPointer((char *)ABFL2dll_NAME);
	ACTVYdll = (ACTVYdll_POINTER) getFunctionPointer((char *)ACTVYdll_NAME);
	AGdll = (AGdll_POINTER) getFunctionPointer((char *)AGdll_NAME);
	CCRITdll = (CCRITdll_POINTER) getFunctionPointer((char *)CCRITdll_NAME);
	CP0dll = (CP0dll_POINTER) getFunctionPointer((char *)CP0dll_NAME);
	CRITPdll = (CRITPdll_POINTER) getFunctionPointer((char *)CRITPdll_NAME);
	CSATKdll = (CSATKdll_POINTER) getFunctionPointer((char *)CSATKdll_NAME);
	CV2PKdll = (CV2PKdll_POINTER) getFunctionPointer((char *)CV2PKdll_NAME);
	CVCPKdll = (CVCPKdll_POINTER) getFunctionPointer((char *)CVCPKdll_NAME);
	CVCPdll = (CVCPdll_POINTER) getFunctionPointer((char *)CVCPdll_NAME);
	DBDTdll = (DBDTdll_POINTER) getFunctionPointer((char *)DBDTdll_NAME);
	DBFL1dll = (DBFL1dll_POINTER) getFunctionPointer((char *)DBFL1dll_NAME);
	DBFL2dll = (DBFL2dll_POINTER) getFunctionPointer((char *)DBFL2dll_NAME);
	DDDPdll = (DDDPdll_POINTER) getFunctionPointer((char *)DDDPdll_NAME);
	DDDTdll = (DDDTdll_POINTER) getFunctionPointer((char *)DDDTdll_NAME);
	DEFLSHdll = (DEFLSHdll_POINTER) getFunctionPointer((char *)DEFLSHdll_NAME);
	DHD1dll = (DHD1dll_POINTER) getFunctionPointer((char *)DHD1dll_NAME);
	DHFLSHdll = (DHFLSHdll_POINTER) getFunctionPointer((char *)DHFLSHdll_NAME);
	DIELECdll = (DIELECdll_POINTER) getFunctionPointer((char *)DIELECdll_NAME);
	DOTFILLdll = (DOTFILLdll_POINTER) getFunctionPointer((char *)DOTFILLdll_NAME);
	DPDD2dll = (DPDD2dll_POINTER) getFunctionPointer((char *)DPDD2dll_NAME);
	DPDDKdll = (DPDDKdll_POINTER) getFunctionPointer((char *)DPDDKdll_NAME);
	DPDDdll = (DPDDdll_POINTER) getFunctionPointer((char *)DPDDdll_NAME);
	DPDTKdll = (DPDTKdll_POINTER) getFunctionPointer((char *)DPDTKdll_NAME);
	DPDTdll = (DPDTdll_POINTER) getFunctionPointer((char *)DPDTdll_NAME);
	DPTSATKdll = (DPTSATKdll_POINTER) getFunctionPointer((char *)DPTSATKdll_NAME);
	DSFLSHdll = (DSFLSHdll_POINTER) getFunctionPointer((char *)DSFLSHdll_NAME);
	ENTHALdll = (ENTHALdll_POINTER) getFunctionPointer((char *)ENTHALdll_NAME);
	ENTROdll = (ENTROdll_POINTER) getFunctionPointer((char *)ENTROdll_NAME);
	ESFLSHdll = (ESFLSHdll_POINTER) getFunctionPointer((char *)ESFLSHdll_NAME);
	FGCTYdll = (FGCTYdll_POINTER) getFunctionPointer((char *)FGCTYdll_NAME);
	FPVdll = (FPVdll_POINTER) getFunctionPointer((char *)FPVdll_NAME);
	GERG04dll = (GERG04dll_POINTER) getFunctionPointer((char *)GERG04dll_NAME);
	GETFIJdll = (GETFIJdll_POINTER) getFunctionPointer((char *)GETFIJdll_NAME);
	GETKTVdll = (GETKTVdll_POINTER) getFunctionPointer((char *)GETKTVdll_NAME);
	GIBBSdll = (GIBBSdll_POINTER) getFunctionPointer((char *)GIBBSdll_NAME);
	HSFLSHdll = (HSFLSHdll_POINTER) getFunctionPointer((char *)HSFLSHdll_NAME);
	INFOdll = (INFOdll_POINTER) getFunctionPointer((char *)INFOdll_NAME);
	LIMITKdll = (LIMITKdll_POINTER) getFunctionPointer((char *)LIMITKdll_NAME);
	LIMITSdll = (LIMITSdll_POINTER) getFunctionPointer((char *)LIMITSdll_NAME);
	LIMITXdll = (LIMITXdll_POINTER) getFunctionPointer((char *)LIMITXdll_NAME);
	MELTPdll = (MELTPdll_POINTER) getFunctionPointer((char *)MELTPdll_NAME);
	MELTTdll = (MELTTdll_POINTER) getFunctionPointer((char *)MELTTdll_NAME);
	MLTH2Odll = (MLTH2Odll_POINTER) getFunctionPointer((char *)MLTH2Odll_NAME);
	NAMEdll = (NAMEdll_POINTER) getFunctionPointer((char *)NAMEdll_NAME);
	PDFL1dll = (PDFL1dll_POINTER) getFunctionPointer((char *)PDFL1dll_NAME);
	PDFLSHdll = (PDFLSHdll_POINTER) getFunctionPointer((char *)PDFLSHdll_NAME);
	PEFLSHdll = (PEFLSHdll_POINTER) getFunctionPointer((char *)PEFLSHdll_NAME);
	PHFL1dll = (PHFL1dll_POINTER) getFunctionPointer((char *)PHFL1dll_NAME);
	PHFLSHdll = (PHFLSHdll_POINTER) getFunctionPointer((char *)PHFLSHdll_NAME);
	PQFLSHdll = (PQFLSHdll_POINTER) getFunctionPointer((char *)PQFLSHdll_NAME);
	PREOSdll = (PREOSdll_POINTER) getFunctionPointer((char *)PREOSdll_NAME);
	PRESSdll = (PRESSdll_POINTER) getFunctionPointer((char *)PRESSdll_NAME);
	PSFL1dll = (PSFL1dll_POINTER) getFunctionPointer((char *)PSFL1dll_NAME);
	PSFLSHdll = (PSFLSHdll_POINTER) getFunctionPointer((char *)PSFLSHdll_NAME);
	PUREFLDdll = (PUREFLDdll_POINTER) getFunctionPointer((char *)PUREFLDdll_NAME);
	RESIDUALdll = (RESIDUALdll_POINTER) getFunctionPointer((char *)RESIDUALdll_NAME);
	QMASSdll = (QMASSdll_POINTER) getFunctionPointer((char *)QMASSdll_NAME);
	QMOLEdll = (QMOLEdll_POINTER) getFunctionPointer((char *)QMOLEdll_NAME);
	SATDdll = (SATDdll_POINTER) getFunctionPointer((char *)SATDdll_NAME);
	SATEdll = (SATEdll_POINTER) getFunctionPointer((char *)SATEdll_NAME);
	SATHdll = (SATHdll_POINTER) getFunctionPointer((char *)SATHdll_NAME);
	SATPdll = (SATPdll_POINTER) getFunctionPointer((char *)SATPdll_NAME);
	SATSdll = (SATSdll_POINTER) getFunctionPointer((char *)SATSdll_NAME);
	SATTdll = (SATTdll_POINTER) getFunctionPointer((char *)SATTdll_NAME);
	SETAGAdll = (SETAGAdll_POINTER) getFunctionPointer((char *)SETAGAdll_NAME);
	SETKTVdll = (SETKTVdll_POINTER) getFunctionPointer((char *)SETKTVdll_NAME);
	SETMIXdll = (SETMIXdll_POINTER) getFunctionPointer((char *)SETMIXdll_NAME);
	SETMODdll = (SETMODdll_POINTER) getFunctionPointer((char *)SETMODdll_NAME);
	SETREFdll = (SETREFdll_POINTER) getFunctionPointer((char *)SETREFdll_NAME);
	SETUPdll = (SETUPdll_POINTER) getFunctionPointer((char *)SETUPdll_NAME);
//		SPECGRdll = (SPECGRdll_POINTER) getFunctionPointer((char *)SPECGRdll_NAME); // not in library
	SUBLPdll = (SUBLPdll_POINTER) getFunctionPointer((char *)SUBLPdll_NAME);
	SUBLTdll = (SUBLTdll_POINTER) getFunctionPointer((char *)SUBLTdll_NAME);
	SURFTdll = (SURFTdll_POINTER) getFunctionPointer((char *)SURFTdll_NAME);
	SURTENdll = (SURTENdll_POINTER) getFunctionPointer((char *)SURTENdll_NAME);
	TDFLSHdll = (TDFLSHdll_POINTER) getFunctionPointer((char *)TDFLSHdll_NAME);
	TEFLSHdll = (TEFLSHdll_POINTER) getFunctionPointer((char *)TEFLSHdll_NAME);
	THERM0dll = (THERM0dll_POINTER) getFunctionPointer((char *)THERM0dll_NAME);
	THERM2dll = (THERM2dll_POINTER) getFunctionPointer((char *)THERM2dll_NAME);
	THERM3dll = (THERM3dll_POINTER) getFunctionPointer((char *)THERM3dll_NAME);
	THERMdll = (THERMdll_POINTER) getFunctionPointer((char *)THERMdll_NAME);
	THFLSHdll = (THFLSHdll_POINTER) getFunctionPointer((char *)THFLSHdll_NAME);
	TPFLSHdll = (TPFLSHdll_POINTER) getFunctionPointer((char *)TPFLSHdll_NAME);
	TPRHOdll = (TPRHOdll_POINTER) getFunctionPointer((char *)TPRHOdll_NAME);
	TQFLSHdll = (TQFLSHdll_POINTER) getFunctionPointer((char *)TQFLSHdll_NAME);
	TRNPRPdll = (TRNPRPdll_POINTER) getFunctionPointer((char *)TRNPRPdll_NAME);
	TSFLSHdll = (TSFLSHdll_POINTER) getFunctionPointer((char *)TSFLSHdll_NAME);
	VIRBdll = (VIRBdll_POINTER) getFunctionPointer((char *)VIRBdll_NAME);
	VIRCdll = (VIRCdll_POINTER) getFunctionPointer((char *)VIRCdll_NAME);
	WMOLdll = (WMOLdll_POINTER) getFunctionPointer((char *)WMOLdll_NAME);
	XMASSdll = (XMASSdll_POINTER) getFunctionPointer((char *)XMASSdll_NAME);
	XMOLEdll = (XMOLEdll_POINTER) getFunctionPointer((char *)XMOLEdll_NAME);
	return COOLPROP_OK;
}



static long i;
static char hfmix[] = "HMX.BNC";
static char hrf[] = "DEF";

#if defined(__ISWINDOWS__)
char refpropPath[] = "";
#elif defined(__ISLINUX__)
char refpropPath[] = "/opt/refprop";
#elif defined(__ISAPPLE__)
char refpropPath[] = "/opt/refprop";
#else
char refpropPath[] = "";
#endif

std::string get_REFPROP_fluid_path()
{
	std::string rpPath (refpropPath);
	#if defined(__ISWINDOWS__)
		return rpPath;
	#elif defined(__ISLINUX__)
		return rpPath + std::string("/fluids/");
	#elif defined(__ISAPPLE__)
		return rpPath + std::string("/fluids/");
	#else
		throw NotImplementedError("This function should not be called.");
		return rpPath;
	#endif
}
bool load_REFPROP()
{
	// If REFPROP is not loaded
	if (RefpropdllInstance==NULL)
	{
		// Load it
		#if defined(__ISWINDOWS__)
			#if defined(ENV64BIT)
				// 64-bit code here.
				TCHAR refpropdllstring[100] = TEXT("refprp64.dll");
				RefpropdllInstance = LoadLibrary(refpropdllstring);
			#elif defined (ENV32BIT)
				// 32-bit code here.
				TCHAR refpropdllstring[100] = TEXT("refprop.dll");
				RefpropdllInstance = LoadLibrary(refpropdllstring);
			#else
				// INCREASE ROBUSTNESS. ALWAYS THROW AN ERROR ON THE ELSE.
				#error "Must define either ENV32BIT or ENV64BIT"
			#endif

			
		#elif defined(__ISLINUX__)
			RefpropdllInstance = dlopen ("librefprop.so", RTLD_LAZY);
		#elif defined(__ISAPPLE__)
			RefpropdllInstance = dlopen ("librefprop.dylib", RTLD_LAZY);
		#else
			throw NotImplementedError("We should not reach this point.");
			RefpropdllInstance = NULL;
		#endif

		if (RefpropdllInstance==NULL)
		{
			#if defined(__ISWINDOWS__)
//				int  dw            = ::GetLastError();
//				char lpBuffer[256] = _T("?");
//				if(dwLastError != 0) {
//				    ::FormatMessage(FORMAT_MESSAGE_FROM_SYSTEM,       // Has to be a system error
//				                     NULL,                            // No formatter
//				                     dw,                              // Get error message for this int
//				                     MAKELANGID(LANG_NEUTRAL,SUBLANG_DEFAULT),  // Use system language
//				                     lpBuffer,                        // Write output
//				                     STR_ELEMS(lpBuffer)-1,           // Length of output
//				                     NULL);
//				}
//				printf(lpBuffer);
//				printf("\n");
				              printf("Could not load refprop.dll \n\n");
				throw AttributeError("Could not load refprop.dll, make sure it is in your system search path. In case you run 64bit and you have a REFPROP license, try installing the 64bit DLL from NIST.");
			#elif defined(__ISLINUX__)
				fputs (dlerror(), stderr);
				              printf("Could not load librefprop.so \n\n");
				throw AttributeError("Could not load librefprop.so, make sure it is in your system search path.");
			#elif defined(__ISAPPLE__)
				fputs (dlerror(), stderr);
				              printf("Could not load librefprop.dylib \n\n");
				throw AttributeError("Could not load librefprop.dylib, make sure it is in your system search path.");
			#else
				throw NotImplementedError("Something is wrong with the platform definition, you should not end up here.");
			#endif
			return false;
		}

		#if defined(__ISWINDOWS__)
		
		// Get data associated with path using the windows libraries, 
		// and if you can (result == 0), the path exists
		#ifdef __MINGW32__
			struct stat buf;
			if ( stat( "c:\\Program Files\\REFPROP\\fluids", &buf) != 0){
				throw ValueError("REFPROP fluid files must be copied to c:\\Program Files\\REFPROP\\fluids");
			}
		#else
			struct _stat buf;
			if ( _stat( "c:\\Program Files\\REFPROP\\fluids", &buf) != 0){
				throw ValueError("REFPROP fluid files must be copied to c:\\Program Files\\REFPROP\\fluids");
			}
		#endif
		#endif

		if (setFunctionPointers()!=COOLPROP_OK)
		{
			              printf("There was an error setting the REFPROP function pointers, check types and names in header file.\n");
			throw AttributeError("There was an error setting the REFPROP function pointers, check types and names in header file.");
			return false;
		}
		return true;
	}
	return true;
}

bool set_REFPROP_fluid(std::string Ref, std::vector<double> &x)
{
	long ierr=0;
	char hf[refpropcharlength*ncmax], herr[errormessagelength+1];
	std::string sRef, components_joined;
	std::string RefString;
	std::string fdPath = get_REFPROP_fluid_path();

	// Check platform support
	if(!REFPROPFluidClass::refpropSupported()){
		throw NotImplementedError("You cannot use the REFPROPFluidClass.");
	}

	// Load REFPROP if it isn't loaded yet
	load_REFPROP();
	
	// If the name of the refrigerant doesn't match 
	// that of the currently loaded refrigerant
	if (LoadedREFPROPRef.compare(Ref))
	{
		// If the fluid name starts with the string "REFPROP-MIX:"
		if (Ref.find("REFPROP-MIX:") == 0)
		{
			// Keep everything after the "REFPROP-MIX:"
			components_joined = Ref.substr(12,Ref.size()-12);
		
			// Sample sRef is "R32[0.697615]&R125[0.302385]" -  this is R410A
			// Or you could do "R410A.mix" to use the full mixture model for this predefined mixture
			
			// Try to process predefined mixtures with .mix or .MIX in the file name
			if (components_joined.find(".mix") != std::string::npos || components_joined.find(".MIX") != std::string::npos)
			{
				char hf[255];
				char hfiles[10000];
				char herr[255];
				double xx[ncmax];
				strcpy(hf,components_joined.c_str());

				SETMIXdll(hf, hfmix, hrf, 
						  &i, hfiles, xx,
						  &ierr, herr,
						  255,
						  255,
						  3, // lengthofreference
						  10000,
						  255);
				// c-string needs to be 0-terminated
				for (unsigned int j = 0; j < 255*ncmax; j++)
				{
					if (hfiles[j] == 32) // empty char
					{
						hfiles[j] = 0;
						break;
					}
				}
                // Resize the vector of mole fractions
				x.resize(i);

				RefString = std::string(hfiles,strlen(hfiles)+1);
				for (int j = 0; j < i; j++)
				{
					x[j] = xx[j];
				}
			}
			else
			{
				// Split the components_joined into the components
				std::vector<std::string> components_split = strsplit(components_joined,'&');

				if (components_split.size() == 1)
				{
					throw ValueError(format("REFPROP mixture specified composition desired [%s], but only one component found",components_joined.c_str()).c_str());
				}

				// Flush out the refrigerant string for REFPROP
				RefString.clear();

				// Resize the vector of mole fractions
				x.resize(components_split.size());

				for (unsigned int j=0;j<components_split.size();j++)
				{	
					// Get component name and mole fraction (as strings)
					std::vector<std::string> comp_fraction = strsplit(components_split[j],'[');

					if (comp_fraction.size() != 2)
					{
						throw ValueError(format("Could not parse name[molefraction] [%s]",components_split[j].c_str()).c_str());
					}
					
					// Build the refrigerant string
					if (j == 0){
						RefString = fdPath + comp_fraction[0]+".fld";
					}
					else{
						RefString += "|" + fdPath + comp_fraction[0]+".fld";
					}
					// Convert the mole fraction (as string) to a number
					x[j] = strtod(comp_fraction[1].c_str(),NULL);

					// Update the number of components
					i = j+1;
				}
			}
		}
		// Name starts with REFPROP-
		else if (Ref.find("REFPROP-") == 0)
		{
			// Keep everything after the "REFPROP-"
			sRef = Ref.substr(8,Ref.size()-8);

			if (!sRef.compare("Air") || !sRef.compare("R507A") || !sRef.compare("R404A") || !sRef.compare("R410A") || !sRef.compare("R407C") || !sRef.compare("SES36"))
			{
				i=1;
				RefString = fdPath + std::string(sRef)+std::string(".ppf");
				x[0]=1.0;     //Pseudo-Pure fluid
			}
			else
			{
				i=1;
				RefString = fdPath + std::string(sRef)+std::string(".fld");
				x[0]=1.0;     //Pure fluid
			}
		}
		else
		{
			throw ValueError(format("REFPROP fluid string [%s] is invalid", Ref.c_str()));
		}

		ierr=999;
		// Set path to fluid files
//		// std::string rpPath (refpropfluidpath);
//		if (rpPath.length()>0)
//		{
//			printf("Setting REFPROP path to: %s\n",rpPath.c_str());
//			char refproppath[refpropcharlength+1];
//			strcpy(refproppath,rpPath.c_str());
//			SETPATHdll(refproppath);
//			free(refproppath);
//		}

		char* hfm = (char*) calloc(refpropcharlength+8, sizeof(char));
		strcpy(hfm,fdPath.c_str());
		strcat(hfm,hfmix);
		strcpy(hf,RefString.c_str());

		//...Call SETUP to initialize the program
		SETUPdll(&i, hf, hfm, hrf, &ierr, herr,
		  refpropcharlength*ncmax,refpropcharlength,
		  lengthofreference,errormessagelength);

		if (ierr > 0){
			//...Call SETUP with capital letters
			for(unsigned int j = 0; j < strlen(hrf); j++)
			{
				hrf[j] = toupper(hrf[j]);
			}
			for(unsigned int j = 0; j < strlen(hfm); j++)
			{
				hfm[j] = toupper(hfm[j]);
			}
			for(unsigned int j = 0; j < strlen(hf); j++)
			{
				hf[j] = toupper(hf[j]);
			}
			SETUPdll(&i, hf, hfm, hrf, &ierr, herr,
			  refpropcharlength*ncmax,refpropcharlength,
			  lengthofreference,errormessagelength);
		}

		if (ierr > 0){
			//...Call SETUP with lower case letters
			for(unsigned int j = 0; j < strlen(hrf); j++)
			{
				hrf[j] = tolower(hrf[j]);
			}
			for(unsigned int j = 0; j < strlen(hfm); j++)
			{
				hfm[j] = tolower(hfm[j]);
			}
			for(unsigned int j = 0; j < strlen(hf); j++)
			{
				hf[j] = tolower(hf[j]);
			}
			SETUPdll(&i, hf, hfm, hrf, &ierr, herr,
			  refpropcharlength*ncmax,refpropcharlength,
			  lengthofreference,errormessagelength);
		}

		free (hfm);

		if (ierr > 0){
			throw ValueError(format("REFPROP: %s",herr).c_str());
			return false;
		}
		else if (ierr < 0)
		{
			set_warning(herr);
		}
		//Copy the name of the loaded refrigerant back into the temporary holder
		LoadedREFPROPRef = std::string(Ref);
		
		unsigned int jmax;
		for (jmax = 0; jmax < ncmax; jmax++)
		{
			if (jmax == x.size())
			{
				break;
			}
			if (x[jmax] < 1e-13)
			{
				break;
			}
		}
		x.resize(jmax);
		LoadedREFPROPx = x;
		return true;
	}
	else
	{
		x = LoadedREFPROPx;
	}
	return true;
}
double REFPROP(char Output, char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
	return REFPROP(std::string(1,Output),std::string(1,Name1),Prop1,std::string(1,Name2),Prop2,std::string(Ref));
}
/*!
From REFPROP:
temperature                     K
pressure, fugacity              kPa
density                         mol/L
composition                     mole fraction
quality                         mole basis (moles vapor/total moles)
enthalpy, internal energy       J/mol
Gibbs, Helmholtz free energy    J/mol
entropy, heat capacity          J/(mol.K)
speed of sound                  m/s
Joule-Thomson coefficient       K/kPa
d(p)/d(rho)                     kPa.L/mol
d2(p)/d(rho)2                   kPa.(L/mol)^2
viscosity                       microPa.s (10^-6 Pa.s)
thermal conductivity            W/(m.K)
dipole moment                   debye
surface tension                 N/m
*/
double REFPROP(std::string Output, std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref)
{
    long iOutput = get_param_index(Output);
	long iName1 = get_param_index(Name1);
	long iName2 = get_param_index(Name2);

    // Return the input if the output is the same as the input
	if (iOutput == iName1){return Prop1;}
	if (iOutput == iName2){return Prop2;}

    // Convert input values to SI
	Prop1 = convert_from_unit_system_to_SI(iName1,Prop1,get_standard_unit_system());
	Prop2 = convert_from_unit_system_to_SI(iName2,Prop2,get_standard_unit_system());

    // Call REFPROPSI
    double output_val = REFPROPSI(iOutput, iName1, Prop1, iName2, Prop2, Ref);

    // Convert back to desired unit system
    return convert_from_SI_to_unit_system(iOutput,output_val,get_standard_unit_system());

}
double REFPROPSI(long iOutput, long iName1, double Prop1, long iName2, double Prop2, std::string Ref)
{
	if(!REFPROPFluidClass::refpropSupported()){
		printf("You cannot use REFPROP, returning.");
		return _HUGE;
	}

	long ierr=0;
	char herr[errormessagelength+1];
	double xliq[ncmax],xvap[ncmax], dummyv[ncmax], output_val;
	
	double TL, TV, dummy;
	double T,p=0,d,dl,dv,q,e,h,s,cv,cp,w,MW,hl,hv,sl,sv,ul,
		uv,pl,pv,hjt,eta,tcx,Q,Tcrit,pcrit,dcrit,sigma;

	// First create a pointer to an instance of the library
	load_REFPROP();
	
	set_REFPROP_fluid(Ref, x);
	
	strcpy(herr,"Ok");
	
	// Return the input if the output is the same as the input
	if (iOutput == iName1){return Prop1;}
	if (iOutput == iName2){return Prop2;}
	
	// Get the molar mass of the fluid
	WMOLdll(&(x[0]),&MW);
	if (iOutput == iTcrit)
	{
		// Critical temperature
		CRITPdll(&(x[0]),&Tcrit,&pcrit,&dcrit,&ierr,herr,255); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
		output_val = Tcrit;
	}
	else if (iOutput==iMM)
	{
		// mole mass
		output_val = MW;
	}
	else if (iOutput==iPcrit)
	{
		// Critical pressure
		CRITPdll(&(x[0]),&Tcrit,&pcrit,&dcrit,&ierr,herr,255); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
		output_val = pcrit*1000;
	}
	else if (iOutput ==iRhocrit)
	{
		// Critical density
		CRITPdll(&(x[0]),&Tcrit,&pcrit,&dcrit,&ierr,herr,255); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
		output_val = dcrit*MW;
		
	}
	else if (iOutput == iTmin)
	{
		// Minimum temperature
		double tmin,tmax,Dmax,pmax;
		LIMITSdll((char *)"EOS",&(x[0]),&tmin,&tmax,&Dmax,&pmax,255); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
		output_val = tmin;
	}
	else if (iOutput == iAccentric)
	{
		double wmm,Ttriple,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas;
		// Accentric factor
		if (i>1)
		{
			fprintf(stderr,"Error: Accentric factor only defined for pure fluids\n");
			output_val = _HUGE;
		}
		INFOdll(&i,&wmm,&Ttriple,&tnbpt,&tc,&pc,&Dc,&Zc,&acf,&dip,&Rgas); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
		output_val = acf;
	}
	else if (iOutput ==iDipole)
	{
		double wmm,Ttriple,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas;
		// Dipole moment
		if (i>1)
		{
			fprintf(stderr,"Error: Dipole moment only defined for pure fluids\n");
			output_val = _HUGE;
		}
		INFOdll(&i,&wmm,&Ttriple,&tnbpt,&tc,&pc,&Dc,&Zc,&acf,&dip,&Rgas); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
		output_val = dip;
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
		INFOdll(&i,&wmm,&Ttriple,&tnbpt,&tc,&pc,&Dc,&Zc,&acf,&dip,&Rgas); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
		output_val = Ttriple;
	}
	else if (iOutput==iI)
	{
		if (iName1==iT){
			SURFTdll(&Prop1,&dl,&(x[0]),&sigma,&i,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			output_val = sigma;
		}
		else{
			std::cout<< "If surface tension is the output, temperature must be the first input" << std::endl;
			output_val = _HUGE;
		}
	}
	
	// Inputs that do not involve Saturation calls with quality imposed directly
	else if (iName1!=iQ && iName2 != iQ)
	{
		if ((iName1==iT && iName2 == iP) || (iName2==  iT && iName1== iP))
		{
			// T in K, P in Pa
			if (iName1 == iP){ std::swap(Prop1,Prop2); }

			T = Prop1; p = Prop2/1000.0; // Want p in [kPa]

			// Use flash routine to find properties
			TPFLSHdll(&T,&p,&(x[0]),&d,&dl,&dv,xliq,xvap,&q,&e,&h,&s,&cv,&cp,&w,&ierr,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
		}
		else if ((iName1==iT && iName2==iD) || (iName2==iT && iName1==iD))
		{
			// T in K, D in kg/m^3
			if (iName2 == iT){
				std::swap(Prop1,Prop2);
			}
			T = Prop1; d = Prop2/MW;
			
			// This is the explicit formulation of the EOS
			TDFLSHdll(&T,&d,&(x[0]),&p,&dl,&dv,xliq,xvap,&q,&e,&h,&s,&cv,&cp,&w,&ierr,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
		}
		else if ((iName1==iP && iName2==iH) || (iName2==iP && iName1==iH))
		{
			// p in Pa, h in J/kg
			if (iName2 == iP){
				std::swap(Prop1,Prop2);
			}
			p = Prop1/1000.0; h = Prop2*MW/1000; // Want h in J/mol, want p in kPa
			
			// Use flash routine to find properties
			PHFLSHdll(&p,&h,&(x[0]),&T,&d,&dl,&dv,xliq,xvap,&q,&e,&s,&cv,&cp,&w,&ierr,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
		}
		else if ((iName1==iP && iName2==iS) || (iName2==iP && iName1==iS))
		{
			// p in Pa, h in J/kg
			if (iName2 == iP){
				std::swap(Prop1,Prop2);
			}
			p = Prop1/1000.0; s = Prop2*MW/1000.0;
			
			// Use flash routine to find properties
			PSFLSHdll(&p,&s,&(x[0]),&T,&d,&dl,&dv,xliq,xvap,&q,&e,&h,&cv,&cp,&w,&ierr,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
		}
		else if ((iName1==iH && iName2==iS) || (iName2==iH && iName1==iS))
		{
			// H in kJ/kg, s in kJ/kg/K
			if (iName2 == iH){
				std::swap(Prop1,Prop2);
			}
			h = Prop1*MW/1000.0; s = Prop2*MW/1000.0;
			
			// Use flash routine to find properties
			HSFLSHdll(&h,&s,&(x[0]),&T,&p,&d,&dl,&dv,xliq,xvap,&q,&e,&cv,&cp,&w,&ierr,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
		}
		else if ((iName1==iP && iName2==iD) || (iName2==iP && iName1==iD))
		{
			// p in Pa, rho in kg/m^3
			if (iName2 == iP){
				std::swap(Prop1,Prop2);
			}
			p = Prop1/1000.0; d = Prop2/MW;
			
			// Use flash routine to find properties
			// from REFPROP: subroutine PDFLSH (p,D,z,t,Dl,Dv,x,y,q,e,h,s,cv,cp,w,ierr,herr)
			PDFLSHdll(&p,&d,&(x[0]),&T,&dl,&dv,xliq,xvap,&q,&e,&h,&s,&cv,&cp,&w,&ierr,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
		}

		// Get the output parameter and convert it to SI units
		switch (iOutput)
		{
		case iT: output_val = T; break;
		case iP: output_val = p*1000; break;
		case iH: output_val = h/MW*1000; break;
		case iD: output_val = d*MW; break;
		case iS: output_val = s/MW*1000; break;
		case iU: output_val = e/MW*1000; break;
		case iC: output_val = cp/MW*1000; break;
		case iO: output_val = cv/MW*1000; break;
		case iA: output_val = w; break;
		case iQ:
			output_val = (1/d-1/dl)/(1/dv-1/dl); break;
		case iV:
			TRNPRPdll(&T,&d,&(x[0]),&eta,&tcx,&ierr,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			output_val = eta/1.0e6; //uPa-s to Pa-s
			break;
		case iL:
			TRNPRPdll(&T,&d,&(x[0]),&eta,&tcx,&ierr,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			output_val = tcx;
			break;
		default:
			throw ValueError(format("Output parameter [%d] is invalid for REFPROP fluid in REFPROP.cpp",iOutput).c_str());
		}
	}
	// Now the two-phase inputs, either T-Q or P-Q
	else if ((iName1==iT && iName2==iQ) || (iName2==iT && iName1==iQ) ||
			 (iName1==iP && iName2==iQ) || (iName2==iP && iName1==iQ)
			 )
	{
		if (iName1 == iT || iName2 == iT)
		{
			if (iName2 == iT){ std::swap(Prop1,Prop2); }

			T = Prop1; Q = Prop2;

			// Saturation Density
			long ic;
			ic=1;
			SATTdll(&T,&(x[0]),&ic,&pl,&dl,&dummy,xliq,dummyv,&ierr,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			ic=2;
			SATTdll(&T,&(x[0]),&ic,&pv,&dummy,&dv,dummyv,xvap,&ierr,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

			p = (pv*Q+pl*(1-Q))*1000; // [Pa]
		}
		else
		{
			if (iName2 == iP){ std::swap(Prop1,Prop2); }

			p = Prop1/1000; Q = Prop2; // p should be in kPa for REFPROP

			// Saturation Density for the liquid
			
			long ic = 1;
			SATPdll(&p,&(x[0]),&ic,&TL,&dl,&dummy,xliq,dummyv,&ierr,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			// Saturation density for the vapor
			ic = 2;
			SATPdll(&p,&(x[0]),&ic,&TV,&dummy,&dv,dummyv,xvap,&ierr,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

			T = (TV*Q+TL*(1-Q));
		}

		if (iOutput == iT)
		{
			output_val = (TV*Q+TL*(1-Q));
		}
		else if (iOutput==iD) 
		{
			output_val = 1/(Q/dv+(1-Q)/dl)*MW;
		}
		else if (iOutput==iP) 
		{
			output_val = p;
		}
		else if (iOutput==iA)
		{
			d=1/(Q/dv+(1-Q)/dl);
			THERMdll(&T,&d,&(x[0]),&p,&e,&h,&s,&cv,&cp,&w,&hjt); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			output_val = w;
		}
		else if (iOutput==iH) 
		{
			ENTHALdll(&T,&dl,xliq,&hl); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			ENTHALdll(&T,&dv,xvap,&hv); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			output_val = (hv*Q+hl*(1-Q))/MW*1000; // kJ/kg to J/kg
		}
		else if (iOutput==iS) 
		{
			ENTROdll(&T,&dl,xliq,&sl); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			ENTROdll(&T,&dv,xvap,&sv); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			output_val = (sv*Q+sl*(1-Q))/MW*1000; // kJ/kg-K to J/kg-K
		}
		else if (iOutput==iU) 
		{
			ENTHALdll(&T,&dl,xliq,&hl); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			ENTHALdll(&T,&dv,xvap,&hv); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			p=pv*Q+pl*(1-Q);
			ul=hl-p/dl;
			uv=hv-p/dv;
			output_val = (uv*Q+ul*(1-Q))/MW*1000; // kJ/kg to J/kg
		}
		else if (iOutput==iC) 
		{
			d = 1/(Q/dv+(1-Q)/dl);
			CVCPdll(&T,&d,&(x[0]),&cv,&cp); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			output_val = cp/MW*1000; // kJ/kg-K to J/kg-K
		}
		else if (iOutput==iO) 
		{
			d = 1/(Q/dv+(1-Q)/dl);
			CVCPdll(&T,&d,&(x[0]),&cv,&cp); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			output_val = cv/MW*1000; // kJ/kg-K to J/kg-K
		}
		else if (iOutput==iV) 
		{
			d = 1/(Q/dv+(1-Q)/dl);
			TRNPRPdll(&T,&d,&(x[0]),&eta,&tcx,&ierr,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			output_val = eta/1.0e6; //uPa-s to Pa-s
		}
		else if (iOutput==iL) 
		{
			d = 1/(Q/dv+(1-Q)/dl);
			TRNPRPdll(&T,&d,&(x[0]),&eta,&tcx,&ierr,herr,errormessagelength); if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}
			output_val = tcx;
		}
		else
		{
			output_val = _HUGE;
		}
	}
	else
	{
		output_val = _HUGE;
	}

	return output_val;
}

REFPROPFluidClass::REFPROPFluidClass(std::string FluidName, std::vector<double> xmol)
{

	long ierr,ic;
	char herr[errormessagelength+1];
	std::vector<double> xliq = xmol, xvap = xmol;
	double Tcrit,dcrit,pcrit,MW,Ttriple,tnbpt,acf,Zcrit,dip,Rgas, dummy1, dummy2;
	
	// Check platform support
	if(!REFPROPFluidClass::refpropSupported()){
	    throw NotImplementedError("You cannot use the REFPROPFluidClass.");
	  }

	// Load REFPROP if not already loaded
	load_REFPROP();

	// Set the fluid
	set_REFPROP_fluid(FluidName, xmol);

	// Copy the molar fractions
	this->xmol = xmol;

	// Molar mass
	WMOLdll(&(xmol[0]),&MW);
	params.molemass = MW;

	// Other parameters
	INFOdll(&i,&MW,&Ttriple,&tnbpt,&Tcrit,&pcrit,&dcrit,&Zcrit,&acf,&dip,&Rgas);
	crit.T = Tcrit;
	crit.rho = dcrit*MW;
	crit.p = PressureUnit(pcrit,UNIT_KPA);

	params.accentricfactor = acf;
	params.R_u = Rgas; // J/(mol*K)
	params.Ttriple = Ttriple;
	limits.Tmin = Ttriple;
	
	ic = xmol.size();
	xliq.resize(ic);
	xvap.resize(ic);
	SATTdll(&(Ttriple),&(xmol[0]),&(ic),&(params.ptriple),&dummy1,&dummy2,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);

	name.assign(FluidName);
	//aliases.push_back(FluidName);
	REFPROPname.assign(FluidName);

	// Set the reducing values from the pointer
	reduce = *preduce;
}

bool REFPROPFluidClass::supported = true; // initialise with true
bool REFPROPFluidClass::refpropSupported () {
	/*
	 * Here we build the bridge from macro definitions
	 * into the actual code. This is also going to be
	 * the central place to handle error messages on
	 * unsupported platforms.
	 */

	// Abort check if Refprop has been loaded.
	if (RefpropdllInstance!=NULL) return true;

	// Store result of previous check.
	if (REFPROPFluidClass::supported) {
		// Either Refprop is supported or it is the first check.
		std::string rpv(RPVersion_NAME);
		if (rpv.compare("NOTAVAILABLE")!=0) {
			// Function names were defined in "REFPROP_lib.h",
			// This platform theoretically supports Refprop.
			if (load_REFPROP()) {
				return true;
			}
			else {
				printf("Good news: It is possible to use REFPROP on your system! However, the library \n");
				printf("could not be loaded. Please make sure that REFPROP is available on your system.\n\n");
				printf("Neither found in current location nor found in system PATH.\n");
				printf("If you already obtained a copy of REFPROP from http://www.nist.gov/srd/, \n");
				printf("add location of REFPROP to the PATH environment variable or your library path.\n\n");
				printf("In case you do not use Windows, have a look at https://github.com/jowr/librefprop.so \n");
				printf("to find instructions on how to compile your own version of the REFPROP library.\n\n");
				REFPROPFluidClass::supported = false;
				return false;
			}
		} else {
			// No definition of function names, we do not expect
			// the Refprop library to be available.
			REFPROPFluidClass::supported = false;
			return false;
		}
	} else {
		return false;
	}
	return false;
}



double REFPROPFluidClass::dphir_dDelta(double tau, double delta)
{
	double p,T,rho,rhobar;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	PRESSdll(&T, &rhobar, &(xmol[0]), &p);
	return 1/delta*(p/(rho*R()*T)-1);
}
double REFPROPFluidClass::dphir_dTau(double tau, double delta)
{
	double T,rho,rhobar,pr,er,hr,sr,cvr,cpr,Ar,Gr;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	RESIDUALdll(&T,&rhobar,&(xmol[0]),&pr,&er,&hr,&sr,&cvr,&cpr,&Ar,&Gr);
	return er/(params.R_u*T*tau);
}
double REFPROPFluidClass::d2phi0_dTau2(double tau, double delta)
{
	double T,rho,rhobar,p0,e0,h0,s0,cv0,cp0,w0,A0,G0;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	THERM0dll(&T,&rhobar,&(xmol[0]),&p0,&e0,&h0,&s0,&cv0,&cp0,&w0,&A0,&G0);
	// cv0 is in J/mol-K, so is R_u
	return -cv0/(params.R_u*tau*tau);
}
double REFPROPFluidClass::d2phir_dTau2(double tau, double delta)
{
	double T,rho,rhobar,pr,er,hr,sr,cvr,cpr,Ar,Gr;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	RESIDUALdll(&T,&rhobar,&(xmol[0]),&pr,&er,&hr,&sr,&cvr,&cpr,&Ar,&Gr);
	// cvr is in J/mol-K, so is R_u
	return -cvr/(params.R_u*tau*tau);
}

double REFPROPFluidClass::d2phir_dDelta_dTau(double tau, double delta)
{
	double T,rho,rhobar,dpdT_constrho;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	DPDTdll(&T, &rhobar, &(xmol[0]), &dpdT_constrho);
	return -1/(delta*tau)*(1/(rho*R())*dpdT_constrho-1-delta*this->dphir_dDelta(tau,delta));
}
double REFPROPFluidClass::d2phir_dDelta2(double tau, double delta)
{
	double T,rho,rhobar,dpdrhobar_constT;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	DPDDdll(&T, &rhobar, &(xmol[0]), &dpdrhobar_constT);
	return 1/(delta*delta)*(1/(params.R_u*T)*dpdrhobar_constT-1-2*delta*this->dphir_dDelta(tau,delta));
}
double REFPROPFluidClass::phir(double tau, double delta)
{
	double T,rho,rhobar,pr,er,hr,sr,cvr,cpr,Ar,Gr;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	RESIDUALdll(&T,&rhobar,&(xmol[0]),&pr,&er,&hr,&sr,&cvr,&cpr,&Ar,&Gr);
	return tau*this->dphir_dTau(tau,delta)-sr/params.R_u;
}
double REFPROPFluidClass::dphi0_dTau(double tau, double delta)
{
	double T,rho,rhobar,p0,e0,h0,s0,cv0,cp0,w0,A0,G0;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	THERM0dll(&T,&rhobar,&(xmol[0]),&p0,&e0,&h0,&s0,&cv0,&cp0,&w0,&A0,&G0);
	return e0/(R()*T*tau)/params.molemass;
}
double REFPROPFluidClass::phi0(double tau, double delta)
{
	double T,rho,rhobar,p0,e0,h0,s0,cv0,cp0,w0,A0,G0;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	THERM0dll(&T,&rhobar,&(xmol[0]),&p0,&e0,&h0,&s0,&cv0,&cp0,&w0,&A0,&G0);
	return (h0-T*s0)/params.R_u/T-1;
}



double REFPROPFluidClass::surface_tension_T(double T)
{
	double dl = 1, sigma;
	long ierr = 0;
	char herr[errormessagelength+1];
	SURFTdll(&T,&dl,&(xmol[0]),&sigma,&ierr,herr,errormessagelength);
	return sigma;
}
double REFPROPFluidClass::viscosity_Trho(double T, double rho)
{
	long ierr = 0;
	char herr[errormessagelength+1];
	double eta,tcx,rhobar = rho/params.molemass;
	
	TRNPRPdll(&T,&rhobar,&(xmol[0]),&eta,&tcx,&ierr,herr,errormessagelength);
	return eta/1e6; //[Pa-s]
}
double REFPROPFluidClass::conductivity_Trho(double T, double rho)
{
	long ierr = 0;
	char herr[errormessagelength+1];
	double eta,tcx,rhobar = rho/params.molemass;
	
	TRNPRPdll(&T,&rhobar,&(xmol[0]),&eta,&tcx,&ierr,herr,errormessagelength);
	return tcx; //[W/m/K]
}

void REFPROPFluidClass::saturation_T(double T, bool UseLUT, double &psatLout, double &psatVout, double &rhosatLout, double &rhosatVout)
{
	long ic,ierr;
	char herr[errormessagelength+1];
	std::vector<double> xliq = xmol, xvap = xmol;
	double dummy;

	ic=1;
	SATTdll(&T,&(xmol[0]),&ic,&psatLout,&rhosatLout,&dummy,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);
	ic=2;
	SATTdll(&T,&(xmol[0]),&ic,&psatVout,&dummy,&rhosatVout,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);

	// Unit conversions
	rhosatLout *= params.molemass;
	rhosatVout *= params.molemass;
	psatLout *= 1000; // 1000 to go from kPa to Pa
	psatVout *= 1000; // 1000 to go from kPa to Pa
}
void REFPROPFluidClass::saturation_p(double p, bool UseLUT, double &TsatLout, double &TsatVout, double &rhosatLout, double &rhosatVout)
{
	long ic,ierr;
	char herr[errormessagelength+1];
	p /= 1000; // 1000 to go from Pa to kPa

	std::vector<double> xliq = xmol,xvap = xmol;
	double dummy;
	ic=1;
	SATPdll(&p,&(xmol[0]),&ic,&TsatLout,&rhosatLout,&dummy,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);
	ic=2;
	SATPdll(&p,&(xmol[0]),&ic,&TsatVout,&dummy,&rhosatVout,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);
	rhosatLout *= params.molemass;
	rhosatVout *= params.molemass;
}
void REFPROPFluidClass::temperature_ph(double p, double h, double &Tout, double &rho, double &rhoL, double &rhoV, double &TsatL, double &TsatV, double T0, double rho0)
{
	double q,e,s,cv,cp,w,hbar = h*params.molemass, dummy1, dummy2;
	long ierr;
	char herr[errormessagelength+1];
	std::vector<double> xliq = xmol, xvap = xmol;

	p /= 1000; // 1000 to go from Pa to kPa
    hbar /= 1000; // 1000 to go from J/kmol to J/mol
	
	PHFLSHdll(&p,&hbar,&(xmol[0]),&Tout,&rho,&rhoL,&rhoV,&(xliq[0]),&(xvap[0]),&q,&e,&s,&cv,&cp,&w,&ierr,herr,errormessagelength);
	rho *= params.molemass;
	rhoL *= params.molemass;
	rhoV *= params.molemass;
	
	// If a single-phase solution, yield impossible saturation densities so that the phase will be properly calculated
	if (double_equal(rhoL, rhoV)){
		rhoL = -2;
		rhoV = -1;
	}
	else{
		// *1000 to go from kPa to Pa
		this->saturation_p(p*1000,false,TsatL,TsatV,dummy1,dummy2);
	}
}
void REFPROPFluidClass::temperature_ps(double p, double s, double &Tout, double &rho, double &rhoL, double &rhoV, double &TsatL, double &TsatV)
{
    double q,e,h,cv,cp,w,sbar = s*params.molemass, dummy1,dummy2;
	long ierr;
    std::vector<double> xliq = xmol, xvap = xmol;
	char herr[errormessagelength+1];

	p /= 1000; // 1000 to go from Pa to kPa
    sbar /= 1000; // 1000 to go from J/kmol/K to J/mol/K
	
	PSFLSHdll(&p,&sbar,&(xmol[0]),&Tout,&rho,&rhoL,&rhoV,&(xliq[0]),&(xvap[0]),&q,&e,&h,&cv,&cp,&w,&ierr,herr,errormessagelength);
	rho *= params.molemass;
	rhoL *= params.molemass;
	rhoV *= params.molemass;
	// If a single-phase solution, yield impossible saturation densities so that the phase will be properly calculated
	if (double_equal(rhoL, rhoV))
	{
		rhoL = -2;
		rhoV = -1;
	}
	else{
		// *1000 to go from kPa to Pa
		this->saturation_p(p*1000,false,TsatL,TsatV,dummy1,dummy2);
	}
}
void REFPROPFluidClass::temperature_hs(double h, double s, double &Tout, double &rho, double &rhoL, double &rhoV, double &TsatL, double &TsatV)
{
    double q,e,cv,cp,w,p,sbar = s*params.molemass, hbar = h*params.molemass, dummy1, dummy2;
	long ierr;
    std::vector<double> xliq = xmol, xvap = xmol;
	char herr[errormessagelength+1];

	hbar /= 1000; // 1000 to go from J/kmol to J/mol
	sbar /= 1000; // 1000 to go from J/kmol/K to J/mol/K
	
	HSFLSHdll(&hbar,&sbar,&(xmol[0]),&Tout,&p,&rho,&rhoL,&rhoV,&(xliq[0]),&(xvap[0]),&q,&e,&cv,&cp,&w,&ierr,herr,errormessagelength);
	rho *= params.molemass;
	rhoL *= params.molemass;
	rhoV *= params.molemass;

    // If a single-phase solution, yield impossible saturation densities so that the phase will be properly calculated
	if (double_equal(rhoL, rhoV))
	{
		rhoL = -2;
		rhoV = -1;
	}
	else{
		// *1000 to go from kPa to Pa
		this->saturation_p(p*1000,false,TsatL,TsatV,dummy1,dummy2);
	}
}
double REFPROPFluidClass::density_Tp(double T, double p)
{
	return this->density_Tp(T,p,0); // guess value is neglected
}
double REFPROPFluidClass::density_Tp(double T, double p, double rho_guess)
{
	long ierr;
	char herr[errormessagelength+1];
	p /= 1000; // 1000 to go from Pa to kPa
	std::vector<double> xliq = xmol,xvap = xmol;
	double q,e,s,cv,cp,w,h,rho,rhoL,rhoV;
	TPFLSHdll(&T,&p,&(xmol[0]),&rho,&rhoL,&rhoV,&(xliq[0]),&(xvap[0]),&q,&e,&h,&s,&cv,&cp,&w,&ierr,herr,errormessagelength);

	return rho*params.molemass;
}
double REFPROPFluidClass::psat(double T)
{
	long ierr,ic;
	char herr[errormessagelength+1];
	std::vector<double> xliq = xmol,xvap = xmol;
	double dummy1,dummy2,psatval;

	ic=1;
	SATTdll(&T,&(xmol[0]),&ic,&psatval,&dummy1,&dummy2,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);

	return psatval*1000; // 1000 to go from kPa to Pa
}
double REFPROPFluidClass::rhosatV(double T)
{
	long ierr,ic;
	char herr[errormessagelength+1];
	std::vector<double> xliq = xmol,xvap = xmol;
	double rhoV,dummy1,psatval;

	ic=2;
	SATTdll(&T,&(xmol[0]),&ic,&psatval,&dummy1,&rhoV,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);
	return rhoV*params.molemass;
}
double REFPROPFluidClass::rhosatL(double T)
{
	long ierr,ic;
	char herr[errormessagelength+1];
	std::vector<double> xliq = xmol,xvap = xmol;
	double rhoL,dummy1,psatval;

	ic=1;
	SATTdll(&T,&(xmol[0]),&ic,&psatval,&rhoL, &dummy1,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);
	return rhoL*params.molemass;
}




/*
!
!
!
!             CATCH TESTS
!
!
!
*/



#ifndef DISABLE_CATCH
#include "Catch/catch.hpp"
TEST_CASE("REFPROP Fluid Class Helmholtz derivatives check", "[helmholtz],[fast]")
{
	std::vector<double> x(1,1);
	REFPROPFluidClass fluid = REFPROPFluidClass("REFPROP-Water",x);
	double eps = sqrt(DBL_EPSILON);

	SECTION("dDelta")
	{
		double ANA = fluid.dphir_dDelta(0.5, 0.5);
		double NUM = (fluid.phir(0.5, 0.5+eps) - fluid.phir(0.5,0.5-eps))/(2*eps);
		REQUIRE(abs(NUM-ANA) < 1e-6);
	}
	SECTION("dTau")
	{
		double ANA = fluid.dphir_dTau(0.5, 0.5);
		double NUM = (fluid.phir(0.5+eps, 0.5) - fluid.phir(0.5-eps,0.5))/(2*eps);
		REQUIRE(abs(NUM-ANA) < 1e-6);
	}
	SECTION("dDelta2")
	{
		double ANA = fluid.d2phir_dDelta2(0.5, 0.5);
		double NUM = (fluid.dphir_dDelta(0.5, 0.5+eps) - fluid.dphir_dDelta(0.5,0.5-eps))/(2*eps);
		REQUIRE(abs(NUM-ANA) < 1e-6);
	}
	SECTION("dTau2")
	{
		double ANA = fluid.d2phir_dTau2(0.5, 0.5);
		double NUM = (fluid.dphir_dTau(0.5+eps, 0.5) - fluid.dphir_dTau(0.5-eps,0.5))/(2*eps);
		REQUIRE(abs(NUM-ANA) < 1e-6);
	}
	SECTION("dDeltadTau")
	{
		double ANA = fluid.d2phir_dDelta_dTau(0.5, 0.5);
		double NUM = (fluid.dphir_dTau(0.5, 0.5+eps) - fluid.dphir_dTau(0.5,0.5-eps))/(2*eps);
		REQUIRE(abs(NUM-ANA) < 1e-6);
	}
}

TEST_CASE("REFPROP Fluid Class check saturation consistency", "")
{
	std::vector<double> x(1,1);
	REFPROPFluidClass fluid = REFPROPFluidClass("REFPROP-Water",x);

	SECTION("sat")
	{
		double T = 313;
		double p, T2, psatV, TsatV,rhoL,rhoV;
		fluid.saturation_T(T, false, p, psatV, rhoL, rhoV);
		fluid.saturation_p(p, false, T2, TsatV, rhoL, rhoV);
		REQUIRE(fabs(T2-T) < 1e-5);
	}
}

TEST_CASE("Check fluid names", "[fast]")
{
	if (REFPROPFluidClass::refpropSupported()) {
		std::vector<double> x1(1,1), x2(2,0.5);
		SECTION("REFPROP-R134")
		{
			REQUIRE_THROWS(set_REFPROP_fluid("REFPROP-R134",x1));
		}
		SECTION("REFPROP-R134a")
		{
			REQUIRE_NOTHROW(set_REFPROP_fluid("REFPROP-R134a",x1));
		}
		SECTION("REFPROP-MIX:R410.m")
		{
			REQUIRE_THROWS(set_REFPROP_fluid("REFPROP-MIX:R410.m",x2));
		}
		SECTION("REFPROP-MIX:R410")
		{
			REQUIRE_THROWS(set_REFPROP_fluid("REFPROP-MIX:R410",x2));
		}
		SECTION("REFPROP-MIX:R32[0.5]R125[0.5]")
		{
			REQUIRE_THROWS(set_REFPROP_fluid("REFPROP-MIX:R32[0.5]R125[0.5]",x2));
		}
		SECTION("REFPROP-MIX:R32[0,5]R125[0,5]")
		{
			REQUIRE_THROWS(set_REFPROP_fluid("REFPROP-MIX:R32[0,5]R125[0,5]",x2));
		}
		SECTION("REFPROP-MIX:R32[A]&R125[B]")
		{
			REQUIRE_THROWS(set_REFPROP_fluid("REFPROP-MIX:R32[A]R125[B]",x2));
		}
		SECTION("REFPROP-MIX:R32[0.697614699375863]&R125[0.302385300624138]")
		{
			REQUIRE_NOTHROW(set_REFPROP_fluid("REFPROP-MIX:R32[0.697614699375863]&R125[0.302385300624138]",x2));
		}
		SECTION("REFPROP-MIXLR410A.mix")
		{
			REQUIRE_THROWS(set_REFPROP_fluid("REFPROP-MIXLR410A.mix",x2));
		}
	}
}

TEST_CASE("Fluid class for bad fluid", "[fast]")
{
	SECTION("REFPROP-R134")
	{
		if (REFPROPFluidClass::refpropSupported()) {
			REQUIRE_THROWS(REFPROPFluidClass("REFPROP-R134",std::vector<double>(1,1)));
		}
	}
	SECTION("REFPROP-R134a")
	{
		if (REFPROPFluidClass::refpropSupported()) {
			REQUIRE_NOTHROW(REFPROPFluidClass("REFPROP-R134a",std::vector<double>(1,1)));
		}
	}
}

#endif
