REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"

cd ../../CoolProp

REM ******* compile all the sources ***************
cl /c /I. /EHsc /DCOOLPROP_LIB *.cpp
cl /c /I. /EHsc /DCOOLPROP_LIB purefluids/*.cpp
cl /c /I. /EHsc /DCOOLPROP_LIB pseudopurefluids/*.cpp

link /DLL CoolProp.obj *.obj /OUT:../lib/CoolProp.dll
erase *.obj