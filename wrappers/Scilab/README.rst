
Using Scilab and CoolProp
=========================

SWIG and Scilab should theoretically work together, but they don't.  So the easiest thing is just to use a dynamically linked library.  

Instructions
============
Download the right DLL for your platform (CoolProp.dll for 32-bit Scilab on windows, CoolProp_x64.dll for 64-bit Scilab on Windows)
    
To use the PropsSI function, do
    
link('CoolProp_x64.dll',['Props','K2F'],'c')


More info
=========
http://mailinglists.scilab.org/Linking-Delphi-DLL-with-Scilab-td3967008.html