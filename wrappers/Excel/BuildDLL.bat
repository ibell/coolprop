REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

cd ../../CoolProp

REM ******* compile all the sources ***************
cl /c /I. /EHsc /DCOOLPROP_LIB *.cpp
cl /c /I. /EHsc /DCOOLPROP_LIB purefluids/*.cpp
cl /c /I. /EHsc /DCOOLPROP_LIB pseudopurefluids/*.cpp

link /DLL CoolProp.obj *.obj /OUT:../lib/CoolProp.dll
dumpbin /EXPORTS ../lib/CoolProp.dll > ../wrappers/Excel/exports.txt
erase *.obj
