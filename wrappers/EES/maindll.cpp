//============================================================================================//
//                                                                                            //
//                                  EES - CoolProp interface                                  //
//                                  -------------------------                                 //
//                                                                                            //
//  This dll is an interface between EES and CoolProp.                                        //
//  In EES, external functions need to be implemented in dynamic librairies.  The first       //
//  argument sent by EES is a 256 characters char variable. The second argument is pointer    //
//  structure conaining "double" values.  The third argument is a linked list for the         //
//  input data                                                                                //
//                                                                                            //
//  The arguments are defined as follows :                                                    //
//  - The string variable contains the the definition of the fluids and of their              //
//    concentrations with the input strings concatenated to the fluid name joined by |        //
//    (e.g. "R134a|T|P|D" or "REFPROP-R134a|O|T|P" or                                         //
//    "REFPROP-MIX:R32[0.697615]&R125[0.302385]|V|P|H" (R410A))                               //
//  - mode, which is -1 if to return a default form of the call as string, normal mode        //
//    otherwise                                                                               //
//  - The last value is a linked list of the input values                                     //
//																							  //
//  The file needs to be built in coolprop_ees.dlf, which is the standard extension           //
//  for EES external functions.  If CoolProp has been built to the static library             //
//  CoolPropStaticLibrary.lib, you can build (with visual studio) CoolProp_EES.dlf with:      //
//                                                                                            //
//     link /DEBUG /DLL main.obj CoolPropStaticLibrary.lib /OUT:COOLPROP_EES.dlf              //
//																							  //
//  Only one unit system is used (modified SI - see help). Future versions might              //
//  include a detection of EES current unit system and its definition in the dll              //
//																							  //
//  Ian Bell                                                                                  //
//  Thermodynamics Laboratory                                                                 //
//  University of Liège                                                                       //
//                                                                                            //
//  January 2013                                                                              //
//============================================================================================//

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h> 
#include <string.h>
#include <vector>
       
#include "CoolPropTools.h"
#include <algorithm>
#include <string>

#include <windows.h>

HINSTANCE CoolPropdllInstance;

typedef double (__stdcall *fp_PropsSdllTYPE)(char*, char*, double, char*, double, char*);
fp_PropsSdllTYPE PropsSdll;

typedef void (__stdcall *fp_set_debug_leveldllTYPE)(int);
fp_set_debug_leveldllTYPE set_debug_leveldll;

typedef long (__stdcall *fp_get_global_param_stringdllTYPE)(char*, char*);
fp_get_global_param_stringdllTYPE get_global_param_stringdll;

typedef long (__stdcall *fp_redirect_stdoutdllTYPE)(char*);
fp_redirect_stdoutdllTYPE redirect_stdoutdll;

static const bool EES_DEBUG = false;

// Structure for handling ees calling syntax
struct EesParamRec {
  double value;
  struct EesParamRec *next;
};

// Tell C++ to use the "C" style calling conventions rather than the C++ mangled names
extern "C"  
{
	__declspec (dllexport)  double COOLPROP_EES(char fluid[256], int mode, struct EesParamRec *input_rec)
	{        
		double In1 = _HUGE, In2 = _HUGE, out;  // Two inputs, one output
		int NInputs;           // Ninputs is the number of inputs
		char NInputs_string[3], err_str[1000];
		std::string fluid_string = fluid;
        
		std::string ErrorMsg, Outstr, In1str, In2str, Fluidstr;
		std::vector<std::string> fluid_split;

		if (mode==-1) {	
			strcpy(fluid,"D = coolprop('D','P',101.325,'Q',0,'R134a')"); 
			return 0;
		}

		// Split the string that is passed in at the | delimiter that was used to join it
		fluid_split = strsplit(fluid_string,'|');
		if (fluid_split.size() != 4) 
		{
			sprintf(err_str,"fluid[%s] length[%d] not 4 elements long",fluid_string.c_str(),fluid_split.size()); 
			strcpy(fluid,err_str);
            if (EES_DEBUG)
            {
                FILE *fp;
                fp = fopen("log.txt","a+");
                fprintf(fp,"%s %s %g %s %g %s\n",Outstr.c_str(),In1str.c_str(),In1,In2str.c_str(),In2,Fluidstr.c_str());
                fprintf(fp,"%s\n",err_str);
                fclose(fp);
            }
			return 0;
		}
		else{
			Fluidstr = upper(fluid_split[0]);
			Outstr = upper(fluid_split[1]);
			In1str = upper(fluid_split[2]);
			In2str = upper(fluid_split[3]);
		}
		
		// Check the number of inputs
		NInputs = 0;
		EesParamRec * aninput_rec = input_rec;
		while (aninput_rec != 0)
		{
			aninput_rec = aninput_rec->next;
			NInputs++;
		};
		if (NInputs != 2) {
			sprintf(NInputs_string,"%d", NInputs);
			strcpy(fluid,NInputs_string);
			return 0;
		}

		// Get the inputs from the pointer structure sent by EES:
		In1 = input_rec->value;
		input_rec = input_rec->next;
		In2 = input_rec->value;
		
		//This block can be used to debug the code by writing output or intermediate values to a text file
		if (EES_DEBUG)
        {
            FILE *fp;
            fp = fopen("log.txt","a+");
            fprintf(fp,"%s %s %17.16g %s %17.16g %s\n",Outstr.c_str(),In1str.c_str(),In1,In2str.c_str(),In2,Fluidstr.c_str());
            fclose(fp);
        }

		// 32-bit code here.
		TCHAR coolpropdllstring[100] = TEXT("CoolProp.dll");
		CoolPropdllInstance = LoadLibrary(coolpropdllstring);

        if (CoolPropdllInstance == NULL)
        {
            strcpy(fluid, "CoolProp.dll could not be loaded");
            FILE *fp;
            fp = fopen("log.txt","a+");
            fprintf(fp,"CoolProp.dll could not be loaded\n");
            fclose(fp);
            return 0.0;
        }
        PropsSdll = (fp_PropsSdllTYPE) GetProcAddress(CoolPropdllInstance,"PropsS");
        set_debug_leveldll = (fp_set_debug_leveldllTYPE) GetProcAddress(CoolPropdllInstance,"set_debug_level");
        get_global_param_stringdll = (fp_get_global_param_stringdllTYPE) GetProcAddress(CoolPropdllInstance,"get_global_param_string");
        redirect_stdoutdll = (fp_redirect_stdoutdllTYPE) GetProcAddress(CoolPropdllInstance,"redirect_stdout");

        if (PropsSdll == NULL)
        {
            strcpy(fluid, "PropsSdll could not be loaded");
            FILE *fp;
            fp = fopen("log.txt","a+");
            fprintf(fp,"PropsSdll could not be loaded\n");
            fclose(fp);
            return 0.0;
        }
        
        if (set_debug_leveldll == NULL)
        {
            strcpy(fluid, "set_debug_leveldll could not be loaded");
            FILE *fp;
            fp = fopen("log.txt","a+");
            fprintf(fp,"set_debug_leveldll could not be loaded\n");
            fclose(fp);
            return 0.0;
        }
        
        if (get_global_param_stringdll == NULL)
        {
            strcpy(fluid, "get_global_param_stringdll could not be loaded");
            FILE *fp;
            fp = fopen("log.txt","a+");
            fprintf(fp,"get_global_param_stringdll could not be loaded\n");
            fclose(fp);
            return 0.0;
        }
        
        if (redirect_stdoutdll == NULL)
        {
            strcpy(fluid, "redirect_stdoutdll could not be loaded");
            FILE *fp;
            fp = fopen("log.txt","a+");
            fprintf(fp,"redirect_stdoutdll could not be loaded\n");
            fclose(fp);
            return 0.0;
        }
        
		if (EES_DEBUG)
        {
            // This redirects all standard output to file
            redirect_stdoutdll("log_stdout_dll.txt");
            set_debug_leveldll(10); // Maximum debugging
        }
        
        out = PropsSdll((char*)Outstr.c_str(),(char*)In1str.c_str(),In1,(char*)In2str.c_str(),In2,(char*)Fluidstr.c_str());

		if (fabs(out)>1e90)
		{
            char err_chars[10000];
            get_global_param_stringdll("errstring",err_chars);
            std::string err_str = err_chars;
            // There was an error
            if (EES_DEBUG)
            {
                FILE *fp;
                fp = fopen("log.txt","a+");
                fprintf(fp,"Error: %s \n",err_str.c_str());
                fclose(fp);
            }
			strcpy(fluid,err_str.c_str());
			return 0.0;
		}
		else
        {
            // Check if there was a warning
            char warn_chars[10000];
            get_global_param_stringdll("warnstring",warn_chars);
            std::string warn_string = warn_chars;
            if (!warn_string.empty())
            {
                if (EES_DEBUG)
                {
                    FILE *fp;
                    fp = fopen("log.txt","a+");
                    fprintf(fp,"Warning: %s \n",warn_string.c_str());
                    fclose(fp);
                }
                // There was a warning, write it back
                strcpy(fluid, warn_string.c_str());
            }
			return out;
        }
	}
};


