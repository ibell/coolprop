REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

REM ******* compile all the sources ***************
cl /c /I../../CoolProp /EHsc /MT /DCOOLPROP_LIB /DCONVENTION=__cdecl ../../CoolProp/*.cpp
cl /c /I../../CoolProp /EHsc /MT /DCOOLPROP_LIB /DCONVENTION=__cdecl ../../CoolProp/purefluids/*.cpp
cl /c /I../../CoolProp /EHsc /MT /DCOOLPROP_LIB /DCONVENTION=__cdecl ../../CoolProp/pseudopurefluids/*.cpp

link /DLL CoolProp.obj *.obj /OUT:CoolProp.dll
dumpbin /EXPORTS CoolProp.dll > exports.txt
erase *.obj
erase *.lib
erase *.exp
