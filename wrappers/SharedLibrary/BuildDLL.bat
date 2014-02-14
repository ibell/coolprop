
@echo off

call:defineEnv x86
cl /c /MP3 /I../../CoolProp /EHsc /DCOOLPROP_LIB ../../CoolProp/*.cpp
link /DLL *.obj /OUT:CoolProp.dll
dumpbin /EXPORTS CoolProp.dll > exports.txt
erase *.obj
erase *.exp
erase *.lib

REM call:defineEnv amd64
REM cl /c /MP /I../../CoolProp /EHsc /DCOOLPROP_LIB ../../CoolProp/*.cpp

REM link /DLL CoolProp.obj *.obj /OUT:CoolProp_x64.dll
REM dumpbin /EXPORTS CoolProp_x64.dll > exports_x64.txt
REM erase *.obj
REM erase *.exp
REM erase *.lib

goto:eof

rem ******** define some general functions ************
:defineEnv    - set the variables, accepts one argument
set stdpaths="C:\Program Files (x86)\Microsoft Visual Studio ","C:\Program Files\Microsoft Visual Studio "
rem this order assures that the latest version is used...
set versions="8.0","9.0","10.0","11.0","12.0"
set relPaths="\VC\vcvarsall.bat"
set filename=""
for %%i in (%stdpaths%) do (
  for %%j in (%versions%) do (
    for %%k in (%relPaths%) do (
      call:loadScript "%%~i%%~j%%~k" %~1
    )
  )
)
goto:eof

:loadScript
rem echo "%~1" "%~2"
if exist "%~1" (
  echo Calling "%~1" %~2 
  call "%~1" %~2 
) else (
  echo Could not find "%~1"
)
goto:eof