REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

REM ******* compile all the sources from CoolProp***************
cl /c /EHsc /I../../CoolProp ../../CoolProp/*.cpp
cl /c /EHsc /I../../CoolProp ../../CoolProp/purefluids/*.cpp
cl /c /EHsc /I../../CoolProp ../../CoolProp/pseudopurefluids/*.cpp
cl /c /EHsc /I../../CoolProp main.cpp
link /DLL main.obj *.obj /OUT:COOLPROP_EES.dlf

erase COOLPROP_EES.exp
erase COOLPROP_EES.lib
erase *.obj