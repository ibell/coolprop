
Pre-requisites:

* A compiler - Visual Studio C++ 2008 Express works well and is required for Python 2.7 since Python 2.7.x built with it.  Visual studio 2010 is not backwards compatible for Python
* Python - you can get by with one of the Python 2.7.x installs from http://www.python.org/download/

1. Run the CoolProp build script in the root of the CoolProp source by calling "python setup.py build".  This will put object files in build_lib folder
2. Run the build script in this folder by calling "python build.py" or by double-clicking on build.py
3. Copy the file REFPROP_wrapper.lib to %DYMOLADIR%\\BIN\\LIB\ (%DYMOLADIR% is DYMOLA's program directory)
4. Copy the file REFPROP_WRAPPER.H to %DYMOLADIR%\\SOURCE\\