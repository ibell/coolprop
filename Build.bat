@echo off

python setup.py build_ext --inplace
python setup.py build_ext --inplace --debug

rmdir /S /Q build

 rem --compiler=mingw32
REM ~ set INPUT=
REM ~ set /P INPUT=Any key to exit.... %=%