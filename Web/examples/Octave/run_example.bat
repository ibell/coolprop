copy ..\..\..\wrappers\Octave\3.6.4\CoolProp.oct .
copy ..\..\..\wrappers\Octave\Example.m .

REM ~ c:\swigwin-3.0.0\swig -octave -c++ -outcurrentdir -o CoolProp_wrap.cpp ../../../CoolProp/CoolProp.i
    REM ~ if %errorlevel% neq 0 exit /b %errorlevel%

REM ~ REM : %%~nf gives just the file name, no path or extension
REM ~ REM : %%f gives the full path and extension
REM ~ for %%f in (..\..\..\CoolProp\*.cpp) do mkoctfile -v -c -I..\..\CoolProp -o %%~nf.o %%f
   REM ~ if %errorlevel% neq 0 exit /b %errorlevel%
REM ~ mkoctfile -v -c -I..\..\..\CoolProp -o CoolProp_wrap.o CoolProp_wrap.cpp
   REM ~ if %errorlevel% neq 0 exit /b %errorlevel%
REM ~ mkoctfile -v -o CoolProp.oct *.o
   REM ~ if %errorlevel% neq 0 exit /b %errorlevel%
REM ~ erase *.o
REM ~ erase CoolProp.lib
REM ~ erase CoolProp.exp
REM ~ erase CoolProp_wrap.cpp

c:\MinGWOctave3.6.4\Octave3.6.4_gcc4.6.2\bin\octave --path %~dp0 %~dp0Example.m > Output.txt
    if %errorlevel% neq 0 exit /b %errorlevel%
erase CoolProp.oct
erase Example.m