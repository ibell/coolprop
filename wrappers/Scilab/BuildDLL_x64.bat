REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat" amd64
REM ******* compile all the sources ***************

cl /MP3 /O2 /Oi /GL /I "..\..\CoolProp" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_USRDLL" /D "LABVIEW_RT_EXPORTS" /D "COOLPROP_LIB" /D "CONVENTION=__cdecl" /D "_WINDLL" /FD /EHsc /MT /Gy /W3 /c /Zi /TP ..\..\CoolProp\*.cpp

link /OUT:"CoolProp_x64.dll" /INCREMENTAL:NO /NOLOGO /DLL /MANIFEST:NO /DEBUG /SUBSYSTEM:WINDOWS /OPT:REF /OPT:ICF /LTCG /DYNAMICBASE /NXCOMPAT /ERRORREPORT:PROMPT  *.obj

dumpbin /EXPORTS CoolProp_x64.dll > exports.txt
erase *.obj
erase *.pdb
erase *.idb
erase *.lib
erase *.exp
