REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

REM ******* compile all the sources from CoolProp***************
cl /c /MP4 /EHsc /I../../CoolProp ../../CoolProp/*.cpp
cl /c /EHsc /I../../CoolProp main.cpp
link /DLL main.obj *.obj /OUT:COOLPROP_EES.dlf

erase COOLPROP_EES.exp
erase COOLPROP_EES.lib
erase *.obj