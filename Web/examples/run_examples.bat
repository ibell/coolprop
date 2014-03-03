cd C++
call run_example.bat
if %errorlevel% neq 0 exit /b %errorlevel%
cd ..\CSharp
call run_example.bat
if %errorlevel% neq 0 exit /b %errorlevel%
cd ..\Java
call run_example.bat
if %errorlevel% neq 0 exit /b %errorlevel%
cd ..\MATLAB
call run_example.bat
if %errorlevel% neq 0 exit /b %errorlevel%
cd ..\Octave
call run_example.bat
if %errorlevel% neq 0 exit /b %errorlevel%
cd ..\Python
call run_example.bat
if %errorlevel% neq 0 exit /b %errorlevel%