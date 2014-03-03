copy ..\..\..\wrappers\Octave\example.m
copy ..\..\..\wrappers\Octave\CoolProp_wrap.cpp

call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"
REM : %%~nf gives just the file name, no path or extension
REM : %%f gives the full path and extension
for %%f in (..\..\..\CoolProp\*.cpp) do mkoctfile -v -c -I..\..\CoolProp -o %%~nf.o %%f
cl /c /I../../../CoolProp /EHsc CoolProp_wrap.cxx
mkoctfile -v -c -I..\..\..\CoolProp -o CoolProp_wrap.o CoolProp_wrap.cpp
cl /c /I../../../CoolProp /EHsc CoolProp_wrap.cxx
mkoctfile -v -o CoolProp.oct *.o
cl /c /I../../../CoolProp /EHsc CoolProp_wrap.cxx
erase *.o
erase CoolProp.lib
erase CoolProp.exp
erase CoolProp_wrap.cpp
octave Example.m > Output.txt
cl /c /I../../../CoolProp /EHsc CoolProp_wrap.cxx