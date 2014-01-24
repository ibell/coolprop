REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

REM ******* compile all the sources from CoolProp***************
cl /c /EHsc /I../../CoolProp maindll.cpp ../../CoolProp/CoolPropTools.cpp
link /DLL maindll.obj CoolPropTools.obj /OUT:COOLPROP_EES.dlf

erase *.obj