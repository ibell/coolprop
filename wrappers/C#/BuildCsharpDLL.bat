REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"

erase *_wrap.cpp

swig.exe -csharp -dllimport "CoolProp" -c++ -outcurrentdir ../../CoolProp/CoolProp.i
cl /c /I../../CoolProp /EHsc CoolProp_wrap.cxx

REM ******* compile all the sources ***************
cl /c /I../../CoolProp /EHsc ../../CoolProp/*.cpp
link /DLL CoolProp_wrap.obj *.obj /OUT:CoolProp.dll
dumpbin /EXPORTS CoolProp.dll > exports.txt
erase *.obj
erase CoolProp_wrap.cxx
erase CoolProp.lib
erase CoolProp.exp
move *.cs VSCsharp
move CoolProp.dll VSCsharp

rem **** Make a zip file using 7-zip ***
7z a -r VSCsharp.zip VSCsharp/*.*