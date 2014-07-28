
mkdir 3.6.4

swig -octave -c++ -outcurrentdir -o CoolProp_wrap.cpp ../../CoolProp/CoolProp.i

rem Change this for your installation
set root_path=C:\MinGWOctave3.6.4\Octave3.6.4_gcc4.6.2\bin

REM : %%~nf gives just the file name, no path or extension
REM : %%f gives the full path and extension
%root_path%\mkoctfile -v -c -I..\..\CoolProp -o CoolProp_wrap.o CoolProp_wrap.cpp
if %errorlevel% neq 0 exit /b %errorlevel%
for %%f in (..\..\CoolProp\*.cpp) do %root_path%\mkoctfile -v -c -I..\..\CoolProp -o %%~nf.o %%f
if %errorlevel% neq 0 exit /b %errorlevel%
%root_path%\mkoctfile -v -o 3.6.4\CoolProp.oct *.o
if %errorlevel% neq 0 exit /b %errorlevel%
erase *.o