@echo off

rmdir /S /Q build
python setup.py build --quiet
REM ~ python setup.py build --debug
python setup.py bdist_msi
python setup.py install 
python testMix.py

REM ~ rmdir /S /Q build


 rem --compiler=mingw32
REM ~ set INPUT=
REM ~ set /P INPUT=Any key to exit.... %=%