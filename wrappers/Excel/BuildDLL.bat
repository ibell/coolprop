REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

REM ******* compile all the sources ***************
cl /c /MP3 /I../../CoolProp /EHsc /D "CONVENTION=__stdcall"  /D "EXPORT_CODE=extern \"C\" __declspec(dllexport)" ../../CoolProp/*.cpp

link /DLL *.obj /OUT:CoolProp.dll
dumpbin /EXPORTS CoolProp.dll > exports.txt
erase *.obj
erase *.exp
erase *.lib