function [ out ] = CoolProp(Output,Variable1,Value1,Variable2,Value2,Ref)
%COOLPROP Wrapper function for CoolProp DLL
    loadlibrary('CoolPropDLL','CoolProp_dll.h')
    out=calllib('CoolPropDLL','Props_dll',int8(Output),int8(Variable1),Value1,int8(Variable2),Value2,Ref);
    unloadlibrary CoolPropDLL

end

