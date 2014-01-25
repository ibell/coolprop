REM ******** set the variables ************
REM call both to ensure that one works
call "C:\Program Files\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"
call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"

REM ******* compile all the sources from CoolProp ***************
cl /c /Ox /fp:fast /I../../../CoolProp /EHsc Example.cpp
cl /c /Ox /MP3 /fp:fast /I../../../CoolProp /EHsc ../../../CoolProp/*.cpp

link *.obj /OUT:Example.exe
erase *.obj

call Example > Output.txt
erase Example.exe