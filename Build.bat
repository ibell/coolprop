@echo off

python setup.py build_ext --inplace

rmdir /S /Q build

 rem --compiler=mingw32
REM ~ set INPUT=
REM ~ set /P INPUT=Any key to exit.... %=%