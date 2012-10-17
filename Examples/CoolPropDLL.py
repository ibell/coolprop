## REFPROP DLL is compiled using the __stdcall calling convention, therefore, use windll
from ctypes import *
cp=windll.LoadLibrary("CoolProp.dll")

Output = c_char_p("D")
Name1 = c_char('T')
Prop1 = c_double(400.0)
Name2 = c_char('P')
Prop2 = c_double(1000.0)
Ref = c_char_p("R134a")

T = c_double(47.0)
cp.F2K.restype = c_double
val= cp.F2K(T)
print '47F in K is',val,'K'

## c++ prototype: double Props(char *Output,char Name1, double Prop1, char Name2, double Prop2, char * Ref)
cp.Props.restype = c_double
Return_value = cp.Props(Output, Name1, Prop1, Name2, Prop2, Ref)
print Return_value