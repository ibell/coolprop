REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

cd ../../CoolProp

REM ******* compile all the sources ***************
cl /c /I. /EHsc *.cpp
cl /c /I. /EHsc purefluids/*.cpp
cl /c /I. /EHsc pseudopurefluids/*.cpp

lib CoolProp.obj *.obj /OUT:../wrappers/EES/CoolPropStaticLibrary.lib
erase *.obj
