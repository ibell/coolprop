call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat"
octave _OctaveBuilder.m
move ..\..\CoolProp\CoolProp.oct 3.6.1\CoolProp.oct

call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"
c:\octave3_6_2\bin\octave _OctaveBuilder.m
move ..\..\CoolProp\CoolProp.oct 3.6.2\CoolProp.oct