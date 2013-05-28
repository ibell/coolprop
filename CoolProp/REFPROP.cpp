#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
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
	return OK;
}



long i;
char hfmix[] = "HMX.BNC";
char hrf[] = "DEF";

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
			#if defined(UNICODE)
				RefpropdllInstance = LoadLibrary((LPCWSTR)"refprop.dll");
			#else
				RefpropdllInstance = LoadLibrary((LPCSTR)"refprop.dll");
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
				printf("Could not load refprop.dll\n\n");
			#elif defined(__ISLINUX__)
				fputs (dlerror(), stderr);
				printf("Could not load librefprop.so\n\n");
			#else
				throw NotImplementedError("You should not be here.");
			#endif
			return false;
		}

		if (setFunctionPointers()!=OK)
		{
			printf("There was an error setting the REFPROP function pointers, check types and names in header file.\n");
			return false;
		}
		return true;
	}
	return true;
}

bool set_REFPROP_fluid(std::string Ref, double *x)
{
	long ierr=0;
	char hf[refpropcharlength*ncmax], herr[errormessagelength+1];
	std::string sRef;
	std::string RefString;
	std::string fdPath = get_REFPROP_fluid_path();

	// Check platform support
	if(!REFPROPFluidClass::refpropSupported()){
		throw NotImplementedError("You cannot use the REFPROPFluidClass.");
	}

	// Load REFPROP if it isn't loaded yet
	load_REFPROP();

	// If the fluid name does not start with the string "REFPROP-"
	if (Ref.find("REFPROP-") == std::string::npos)
	{
		// Fail and give error
		std::cout << "Invalid REFPROP string: " << Ref.c_str() << std::endl;
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
//		std::string refpropFluidPath(refpropfluidpath);
//		std::string hfmixStr = refpropFluidPath + std::string(hfmix);
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
		else if (!sRef.compare("Air") || !sRef.compare("R507A") || !sRef.compare("R404A") || !sRef.compare("R410A") || !sRef.compare("R407C") || !sRef.compare("SES36"))
		{
			i=1;
			RefString = fdPath + std::string(sRef)+std::string(".ppf");
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
			RefString = fdPath + std::string(sRef)+std::string(".fld");
			x[0]=1.0;     //Pure fluid
		}
		
		strcpy(hf,RefString.c_str());

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

		//...Call SETUP to initialize the program
		char* hfm = (char*) calloc(refpropcharlength+8, sizeof(char));
		strcpy(hfm,fdPath.c_str());
		strcat(hfm,hfmix);

		SETUPdll(&i, hf, hfm, hrf, &ierr, herr,
				refpropcharlength*ncmax,refpropcharlength,
				lengthofreference,errormessagelength);
		free (hfm);
		

		if (ierr != 0) printf("REFPROP setup gives this error during SETUP: %s\n",herr);
		//Copy the name of the loaded refrigerant back into the temporary holder
		LoadedREFPROPRef = std::string(Ref);
		return true;
	}
	return true;
}
double REFPROP(char Output, char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
	return REFPROP(std::string(1,Output),std::string(1,Name1),Prop1,std::string(1,Name2),Prop2,std::string(Ref));
}
double REFPROP(std::string Output, std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref)
{
	if(!REFPROPFluidClass::refpropSupported()){
		printf("You cannot use REFPROP, returning.");
		return _HUGE;
	}

	long ierr=0,iOutput,iName1,iName2;
	char herr[errormessagelength+1];
	double x[ncmax],xliq[ncmax],xvap[ncmax];
	
	double T,p=0,d,dl,dv,dl_,dv_,q,e,h,s,cv,cp,w,MW,hl,hv,sl,sv,ul,
		uv,pl,pv,hjt,eta,tcx,Q,Tcrit,pcrit,dcrit,rho,sigma;

	// First create a pointer to an instance of the library
	load_REFPROP();
	
	set_REFPROP_fluid(Ref, x);
	
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
		LIMITSdll((char *)"EOS",x,&tmin,&tmax,&Dmax,&pmax,255);
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
		else if (iOutput == iQ)
		{
			return q;
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
		
		if (iOutput == iQ){return Q;}

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
		if (iOutput == iQ){return Q;}

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
	else if ((iName1==iH && iName2==iS) || (iName2==iH && iName1==iS))
	{
		// H in kJ/kg, s in kJ/kg/K
		if (iName2 == iH){
			std::swap(Prop1,Prop2);
		}
		h = Prop1*MW; s = Prop2*MW;
		
		// Use flash routine to find properties
		HSFLSHdll(&h,&s,x,&T,&p,&d,&dl,&dv,xliq,xvap,&q,&e,&cv,&cp,&w,&ierr,herr,errormessagelength);
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

REFPROPFluidClass::REFPROPFluidClass(std::string FluidName, std::vector<double> xmol)
{

	long ierr,ic;
	char herr[errormessagelength+1];
	std::vector<double> xliq = std::vector<double>(1,1),xvap = std::vector<double>(1,1);
	double Tcrit,dcrit,pcrit,MW,Ttriple,tnbpt,acf,Zcrit,dip,Rgas, dummy1, dummy2;
	
	// Check platform support
	if(!REFPROPFluidClass::refpropSupported()){
	    throw NotImplementedError("You cannot use the REFPROPFluidClass.");
	  }

	// Copy the molar fractions
	this->xmol = xmol;

	// Load REFPROP if not already loaded
	load_REFPROP();

	// Set the fluid
	set_REFPROP_fluid(FluidName, &(xmol[0]));

	// Molar mass
	WMOLdll(&(xmol[0]),&MW);
	params.molemass = MW;

	// Other parameters
	INFOdll(&i,&MW,&Ttriple,&tnbpt,&Tcrit,&pcrit,&dcrit,&Zcrit,&acf,&dip,&Rgas);
	crit.T = Tcrit;
	crit.rho = dcrit*MW;
	crit.p = pcrit;

	params.accentricfactor = acf;
	params.R_u = Rgas;
	params.Ttriple = Ttriple;
	limits.Tmin = Ttriple;
	
	ic = 1;
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
	double p,T,rho,R,rhobar;

	R = params.R_u/params.molemass;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	PRESSdll(&T,&rhobar,&(xmol[0]),&p);
	return 1/delta*(p/(rho*R*T)-1);
}
double REFPROPFluidClass::dphir_dTau(double tau, double delta)
{
	double T,rho,R,rhobar,pr,er,hr,sr,cvr,cpr,Ar,Gr;

	R = params.R_u/params.molemass;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	RESIDUALdll(&T,&rhobar,&(xmol[0]),&pr,&er,&hr,&sr,&cvr,&cpr,&Ar,&Gr);
	return er/(R*T*tau)/params.molemass;
}
double REFPROPFluidClass::d2phi0_dTau2(double tau, double delta)
{
	double T,rho,R,rhobar,p0,e0,h0,s0,cv0,cp0,w0,A0,G0;

	R = params.R_u/params.molemass;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	THERM0dll(&T,&rhobar,&(xmol[0]),&p0,&e0,&h0,&s0,&cv0,&cp0,&w0,&A0,&G0);
	return -cv0/(R*tau*tau)/params.molemass;
}
double REFPROPFluidClass::d2phir_dTau2(double tau, double delta)
{
	double T,rho,R,rhobar,pr,er,hr,sr,cvr,cpr,Ar,Gr;

	R = params.R_u/params.molemass;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	RESIDUALdll(&T,&rhobar,&(xmol[0]),&pr,&er,&hr,&sr,&cvr,&cpr,&Ar,&Gr);
	return -cvr/(R*tau*tau)/params.molemass;
}

double REFPROPFluidClass::d2phir_dDelta_dTau(double tau, double delta)
{
	double T,rho,R,rhobar,dpdT_constrho;

	R = params.R_u/params.molemass;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	DPDTdll(&T, &rhobar, &(xmol[0]), &dpdT_constrho);
	return -1/(delta*tau)*(1/(rho*R)*dpdT_constrho-1-delta*this->dphir_dDelta(tau,delta));
}
double REFPROPFluidClass::d2phir_dDelta2(double tau, double delta)
{
	double T,rho,R,rhobar,dpdrhobar_constT;

	R = params.R_u/params.molemass;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	DPDDdll(&T, &rhobar, &(xmol[0]), &dpdrhobar_constT);
	return 1/(delta*delta)*(1/(params.R_u*T)*dpdrhobar_constT-1-2*delta*this->dphir_dDelta(tau,delta));
}
double REFPROPFluidClass::phir(double tau, double delta)
{
	double T,rho,R,rhobar,pr,er,hr,sr,cvr,cpr,Ar,Gr;

	R = params.R_u/params.molemass;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	RESIDUALdll(&T,&rhobar,&(xmol[0]),&pr,&er,&hr,&sr,&cvr,&cpr,&Ar,&Gr);

	return tau*this->dphir_dTau(tau,delta)-sr/R/params.molemass;
}
double REFPROPFluidClass::dphi0_dTau(double tau, double delta)
{
	double T,rho,R,rhobar,p0,e0,h0,s0,cv0,cp0,w0,A0,G0;

	R = params.R_u/params.molemass;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	THERM0dll(&T,&rhobar,&(xmol[0]),&p0,&e0,&h0,&s0,&cv0,&cp0,&w0,&A0,&G0);
	return e0/(R*T*tau)/params.molemass;
}
double REFPROPFluidClass::phi0(double tau, double delta)
{
	double T,rho,R,rhobar,p0,e0,h0,s0,cv0,cp0,w0,A0,G0;

	R = params.R_u/params.molemass;
	rho = delta*reduce.rho;
	rhobar = rho/params.molemass;
	T = reduce.T/tau;
	
	THERM0dll(&T,&rhobar,&(xmol[0]),&p0,&e0,&h0,&s0,&cv0,&cp0,&w0,&A0,&G0);
	return (h0-T*s0)/params.molemass/R/T-1;
}

double REFPROPFluidClass::viscosity_Trho(double T, double rho)
{
	long ierr = 0;
	char herr[errormessagelength+1];
	double eta,tcx,rhobar = rho/params.molemass;;
	
	TRNPRPdll(&T,&rhobar,&(xmol[0]),&eta,&tcx,&ierr,herr,errormessagelength);
	return eta/1e6;
}
double REFPROPFluidClass::conductivity_Trho(double T, double rho)
{
	long ierr = 0;
	char herr[errormessagelength+1];
	double eta,tcx,rhobar = rho/params.molemass;
	
	TRNPRPdll(&T,&rhobar,&(xmol[0]),&eta,&tcx,&ierr,herr,errormessagelength);
	return tcx/1e3;
}

void REFPROPFluidClass::saturation_T(double T, bool UseLUT, double *psatLout, double *psatVout, double *rhosatLout, double *rhosatVout)
{
	long ic,ierr;
	char herr[errormessagelength+1];
	std::vector<double> xliq = std::vector<double>(1,1),xvap = std::vector<double>(1,1);
	double dummy;
	ic=1;
	SATTdll(&T,&(xmol[0]),&ic,psatLout,rhosatLout,&dummy,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);
	ic=2;
	SATTdll(&T,&(xmol[0]),&ic,psatVout,&dummy,rhosatVout,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);
	*rhosatLout *= params.molemass;
	*rhosatVout *= params.molemass;
}
void REFPROPFluidClass::saturation_p(double p, bool UseLUT, double *TsatLout, double *TsatVout, double *rhosatLout, double *rhosatVout)
{
	long ic,ierr;
	char herr[errormessagelength+1];
	std::vector<double> xliq = std::vector<double>(1,1),xvap = std::vector<double>(1,1);
	double dummy;
	ic=1;
	SATPdll(&p,&(xmol[0]),&ic,TsatLout,rhosatLout,&dummy,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);
	ic=2;
	SATPdll(&p,&(xmol[0]),&ic,TsatVout,&dummy,rhosatVout,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);
	*rhosatLout *= params.molemass;
	*rhosatVout *= params.molemass;
}
void REFPROPFluidClass::temperature_ph(double p, double h, double *Tout, double *rhoout, double *rhoLout, double *rhoVout, double *TsatLout, double *TsatVout, double T0, double rho0)
{
	long ierr;
	char herr[errormessagelength+1];
	std::vector<double> xliq = std::vector<double>(1,1),xvap = std::vector<double>(1,1);
	double q,e,s,cv,cp,w,hbar = h*params.molemass, dummy1, dummy2;
	PHFLSHdll(&p,&hbar,&(xmol[0]),Tout,rhoout,rhoLout,rhoVout,&(xliq[0]),&(xvap[0]),&q,&e,&s,&cv,&cp,&w,&ierr,herr,errormessagelength);
	*rhoout *= params.molemass;
	*rhoLout *= params.molemass;
	*rhoVout *= params.molemass;
	this->saturation_p(p,false,TsatLout,TsatVout,&dummy1,&dummy2);
}
void REFPROPFluidClass::temperature_ps(double p, double s, double *Tout, double *rhoout, double *rhoLout, double *rhoVout, double *TsatLout, double *TsatVout)
{
	long ierr;
	char herr[errormessagelength+1];
	std::vector<double> xliq = std::vector<double>(1,1),xvap = std::vector<double>(1,1);
	double q,e,h,cv,cp,w,sbar = s*params.molemass;
	PSFLSHdll(&p,&sbar,&(xmol[0]),Tout,rhoout,rhoLout,rhoVout,&(xliq[0]),&(xvap[0]),&q,&e,&h,&cv,&cp,&w,&ierr,herr,errormessagelength);
	*rhoout *= params.molemass;
	*rhoLout *= params.molemass;
	*rhoVout *= params.molemass;
}
void REFPROPFluidClass::temperature_hs(double h, double s, double *Tout, double *rhoout, double *rhoLout, double *rhoVout, double *TsatLout, double *TsatVout)
{
	long ierr;
	char herr[errormessagelength+1];
	std::vector<double> xliq = std::vector<double>(1,1),xvap = std::vector<double>(1,1);
	double q,e,cv,cp,w,p,sbar = s*params.molemass, hbar = h*params.molemass, dummy1, dummy2;
	HSFLSHdll(&hbar,&sbar,&(xmol[0]),Tout,&p,rhoout,rhoLout,rhoVout,&(xliq[0]),&(xvap[0]),&q,&e,&cv,&cp,&w,&ierr,herr,errormessagelength);
	*rhoout *= params.molemass;
	*rhoLout *= params.molemass;
	*rhoVout *= params.molemass;
	this->saturation_p(p,false,TsatLout,TsatVout,&dummy1,&dummy2);
}
double REFPROPFluidClass::density_Tp(double T, double p)
{
	return this->density_Tp(T,p,0);
}
double REFPROPFluidClass::density_Tp(double T, double p, double rho_guess)
{
	long ierr;
	char herr[errormessagelength+1];
	std::vector<double> xliq = std::vector<double>(1,1),xvap = std::vector<double>(1,1);
	double q,e,s,cv,cp,w,h,rho,rhoL,rhoV;
	TPFLSHdll(&T,&p,&(xmol[0]),&rho,&rhoL,&rhoV,&(xliq[0]),&(xvap[0]),&q,&e,&h,&s,&cv,&cp,&w,&ierr,herr,errormessagelength);

	return rho*params.molemass;
}
double REFPROPFluidClass::psat(double T)
{
	long ierr,ic;
	char herr[errormessagelength+1];
	std::vector<double> xliq = std::vector<double>(1,1),xvap = std::vector<double>(1,1);
	double dummy1,dummy2,psatval;

	ic=1;
	SATTdll(&T,&(xmol[0]),&ic,&psatval,&dummy1,&dummy2,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);

	return psatval;
}
double REFPROPFluidClass::rhosatV(double T)
{
	long ierr,ic;
	char herr[errormessagelength+1];
	std::vector<double> xliq = std::vector<double>(1,1),xvap = std::vector<double>(1,1);
	double rhoV,dummy1,psatval;

	ic=2;
	SATTdll(&T,&(xmol[0]),&ic,&psatval,&dummy1,&rhoV,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);
	return rhoV*params.molemass;
}
double REFPROPFluidClass::rhosatL(double T)
{
	long ierr,ic;
	char herr[errormessagelength+1];
	std::vector<double> xliq = std::vector<double>(1,1),xvap = std::vector<double>(1,1);
	double rhoL,dummy1,psatval;

	ic=1;
	SATTdll(&T,&(xmol[0]),&ic,&psatval,&rhoL, &dummy1,&(xliq[0]),&(xvap[0]),&ierr,herr,errormessagelength);
	return rhoL*params.molemass;
}
