
@echo on
copy ..\..\..\wrappers\Java\Example.java Example.java

REM ******** set the variables ************
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
call "C:\Program Files\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"

REM ******* compile all the sources ***************
swig -java -c++ -outcurrentdir ../../../CoolProp/CoolProp.i
cl /c /I../../../CoolProp /I"C:\Program Files\Java\jdk1.7.0_40\include" /I"C:\Program Files\Java\jdk1.7.0_40\include\win32" /EHsc *.cxx
cl /c /I../../../CoolProp /I"C:\Program Files\Java\jdk1.7.0_40\include" /I"C:\Program Files\Java\jdk1.7.0_40\include\win32" /EHsc ../../../CoolProp/*.cpp
link /DLL *.obj /OUT:CoolProp.dll
erase *.obj
erase *.exp
erase *.lib

javac *.java
java Example > Output.txt

erase *.java
erase *.class
erase CoolProp.dll
erase CoolProp_wrap.cxx
copy ..\..\..\wrappers\Java\Example.java Example.java