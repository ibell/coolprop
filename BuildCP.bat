@echo off

rmdir /S /Q build
REM ~ python setup.py build --quiet
python setup.py bdist_msi
python setup.py install
REM ~ python testMix.py