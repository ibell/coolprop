copy ..\..\..\wrappers\Python\Example.py
cl /c /I../../../CoolProp /EHsc CoolProp_wrap.cxx
python Example.py > Output.txt
cl /c /I../../../CoolProp /EHsc CoolProp_wrap.cxx