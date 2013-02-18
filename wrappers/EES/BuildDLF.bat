REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

REM ******* compile all the sources from CoolProp***************
cl /c /EHsc /I../../CoolProp /DEBUG ../../CoolProp/*.cpp
cl /c /EHsc /I../../CoolProp /DEBUG ../../CoolProp/purefluids/*.cpp
cl /c /EHsc /I../../CoolProp /DEBUG ../../CoolProp/pseudopurefluids/*.cpp
cl /c /EHsc /I../../CoolProp /DEBUG main.cpp
link /DLL main.obj ../../CoolProp/*.obj ../../CoolProp/purefluids/*.obj ../../CoolProp/pseudopurefluids/*.obj /OUT:COOLPROP_EES.dlf /WX

REM ~ erase COOLPROP_EES.exp
erase COOLPROP_EES.lib
erase *.obj