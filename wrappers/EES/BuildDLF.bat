REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"

cl /c /DEBUG /EHsc /I../../CoolProp main.cpp
link /DEBUG /DLL main.obj CoolPropStaticLibrary.lib /OUT:COOLPROP_EES.dlf

copy CoolProp.lib c:\EES32\Userlib\CoolProp_EES
copy CoolProp_EES.dlf c:\EES32\Userlib\CoolProp_EES

erase COOLPROP_EES.exp
erase COOLPROP_EES.ilk
erase COOLPROP_EES.lib
erase COOLPROP_EES.pdb
erase main.obj