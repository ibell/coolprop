
#include "../CoolProp.h"
#include "CoolProp_dll.h"

__declspec(dllexport) void  __stdcall Help_dll(void)
{
	Help();
	return; 
}

__declspec(dllexport) double  __stdcall Props_dll(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
{
	return Props(Output,Name1,Prop1,Name2,Prop2,Ref);
}

__declspec(dllexport) int  __stdcall errCode_dll(char * Ref)
{
	return errCode(Ref);
}

__declspec(dllexport) double  __stdcall pcrit_dll(char *Ref)
{
	return pcrit(Ref);
}

__declspec(dllexport) double  __stdcall Tcrit_dll(char *Ref)
{
	return Tcrit(Ref);
}

__declspec(dllexport) double  __stdcall T_hp_dll(char *Ref, double h, double p, double T_guess)
{
	return T_hp(Ref,h,p,T_guess);
}

__declspec(dllexport) double  __stdcall h_sp_dll(char *Ref, double s, double p, double T_guess)
{
	return h_sp(Ref,s,p,T_guess);
}

__declspec(dllexport) double  __stdcall Tsat_dll(char *Ref, double p, double Q, double T_guess)
{
	return Tsat(Ref,p,Q,T_guess);
}