
mkdir 3.8.1

swig -octave -c++ -outcurrentdir -o CoolProp_wrap.cpp ../../CoolProp/CoolProp.i

rem Change this for your installation
set root_path=C:\octave-3.8.1-5-portable\octave-3.8.1\bin

REM : %%~nf gives just the file name, no path or extension
REM : %%f gives the full path and extension
%root_path%\mkoctfile -v -c -I..\..\CoolProp -IC:\octave-3.8.1-5-portable\octave-3.8.1\include\octave-3.8.1 -IC:\octave-3.8.1-5-portable\octave-3.8.1\include\octave-3.8.1\octave -IC:\octave-3.8.1-5-portable\octave-3.8.1\include -o CoolProp_wrap.o CoolProp_wrap.cpp
if %errorlevel% neq 0 exit /b %errorlevel%
for %%f in (..\..\CoolProp\*.cpp) do %root_path%\mkoctfile -v -c -I..\..\CoolProp -o %%~nf.o %%f
if %errorlevel% neq 0 exit /b %errorlevel%
%root_path%\mkoctfile -v -LC:\octave-3.8.1-5-portable\octave-3.8.1\lib\octave\3.8.1 -o 3.8.1\CoolProp.oct *.o
if %errorlevel% neq 0 exit /b %errorlevel%
erase *.o