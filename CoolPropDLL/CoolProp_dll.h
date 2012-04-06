__declspec(dllexport) double  __stdcall Props_dll(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref);
__declspec(dllexport) double  __stdcall pcrit_dll(char *Ref);
__declspec(dllexport) double  __stdcall Tcrit_dll(char *Ref);
__declspec(dllexport) double  __stdcall T_hp_dll(char *Ref, double h, double p, double T_guess);
__declspec(dllexport) double  __stdcall h_sp_dll(char *Ref, double s, double p, double T_guess);
__declspec(dllexport) double  __stdcall Tsat_dll(char *Ref, double p, double Q, double T_guess);
__declspec(dllexport) double  __stdcall Props_dll(char Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref);
__declspec(dllexport) void  __stdcall Help_dll();
