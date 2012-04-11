import os
os.environ['PATH']+=';..\\CoolPropDLL' #Add the path to the CoolPropDLL

from ctypes import *
cp=windll.LoadLibrary("CoolProp_dll.dll")

Output = c_char('D')
Name1 = c_char('T')
Prop1 = c_double(400)
Name2 = c_char('P')
Prop2 = c_double(1000)
Ref = c_char_p("REFPROP-R134a")

cp.Props_dll.restype = c_double

Return_value = cp.Props_dll(Output, Name1, Prop1, Name2, Prop2, Ref)
print Return_value