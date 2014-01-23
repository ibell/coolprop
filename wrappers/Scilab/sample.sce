link('CoolProp_x64.dll',['PropsSIScilab','plus_one'],'c');
//link('show');
[two] = call("plus_one",1.3,1,"d","out",[1,1],2,"d")
disp(two)

funcprot(0)
function [out]=Props(Output,Input1,Value1,Input2,Value2,Name)
  out = call("PropsSIScilab",Output,1,"c",Input1,2,"c",Value1,3,"d",Input2,4,"c",Value2,5,"d",Name,6,"c","out",[1,1],7,"d")
endfunction

[rho] = Props("D","T",298.15,"P",101325.0,"Air")
disp(rho)
ulink()
