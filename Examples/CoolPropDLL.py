import os
os.environ['PATH']+=';..\\lib' #Add the path to the CoolProp DLL

#REFPROP DLL is compiled using the __cdecl calling convention, therefore, use cdll
from ctypes import *
cp=cdll.LoadLibrary("CoolProp.dll")

Output = c_char_p("D")
Name1 = c_char('T')
Prop1 = c_double(400.0)
Name2 = c_char('P')
Prop2 = c_double(1000.0)
Ref = c_char_p("REFPROP-R134a")

T = c_double(47.0)
cp.F2K.restype = c_double
val= cp.F2K(T)
print '47F in K is',val,'K'

cp.Props.restype = c_double
Return_value = cp.Props(Output, Name1, Prop1, Name2, Prop2, Ref)
print Return_value