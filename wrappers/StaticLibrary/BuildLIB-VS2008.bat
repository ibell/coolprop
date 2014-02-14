REM ******** set the variables ************
REM call both to ensure that one works
call "C:\Program Files\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"

REM ******* compile all the sources from CoolProp ***************
cl /c /MP3 /I../../CoolProp /EHsc /MDd ../../CoolProp/*.cpp
lib *.obj /OUT:VS2008/CoolPropLib_MDd.lib
erase *.obj

REM ******* compile all the sources from CoolProp ***************
cl /c /MP3 /I../../CoolProp /EHsc /MD ../../CoolProp/*.cpp
lib *.obj /OUT:VS2008/CoolPropLib_MD.lib
erase *.obj

REM ******* compile all the sources from CoolProp ***************
cl /c /MP3 /I../../CoolProp /EHsc /MT ../../CoolProp/*.cpp
lib *.obj /OUT:VS2008/CoolPropLib_MT.lib
erase *.obj

REM ******* compile all the sources from CoolProp ***************
cl /c /MP3 /I../../CoolProp /EHsc /MTd ../../CoolProp/*.cpp
lib *.obj /OUT:VS2008/CoolPropLib_MTd.lib
erase *.obj