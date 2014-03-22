Building
========

Requirements
------------
SWIG (http://www.swig.org/)
A compiler (the Visual Studio 2010 C++ Express version should be fine)

To Build
--------

**Windows**:

Run the script build_x64.bat - adjust the paths if necessary to the include folders for your java installation

If on 32-bit windows, run the build_win32.bat file

Each script will put the DLL in the corresponding folder (win32 for 32-bit, x64 for 64-bit)

----

**Linux**: 

**dependencies**: ``swig``, ``g++``, ``JDK(openjdk or oracle)``

Use your package manager to install the dependencies, then run ``./build_linux.sh``. It will build and install a system-wide shared library.

This script will first use the ``$JAVA_HOME`` variable as JDK path. If not set, will search ``/usr/lib/jvm/java-7(6)-openjdk-*`` directory to find out whether ``openjdk`` 7 or 6 is installed.

If you are using ``oracle JDK`` or your ``openjdk`` is not in ``/usr/lib/jvm/java-7(6)-openjdk-*``, you should set ``$JAVA_HOME``.

**Note**: If you are using ``oracle JDK``, please copy the file ``oracle_jdk_path/include/linux/jni_md.h`` to ``oracle_jdk_path/include/``.

Example::

    export JAVA_HOME="/usr/local/jdk1.7.0_45"
    ./build_linux.sh

Running
=======
At the console, run::

    javac *.java
    java runme
    
which should output::

    702.820647167934
    
Hiccups
=======
If the bin folder of the installation for java is not on the path, you may need to add it.


