@echo off

rmdir /S /Q build
python setup.py build
python setup.py bdist_msi
python setup.py install
python setup.py clean
rmdir /S /Q build

