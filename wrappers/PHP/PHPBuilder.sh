swig -php -outcurrentdir -c++ -o CoolProp_wrap.cxx ../../CoolProp/CoolProp.i

g++ `php-config --includes` -c -fpic -Wall -I../../CoolProp CoolProp_wrap.cxx

# ******* compile all the sources from CoolProp***************
g++ -c -fpic -Wall -I../../CoolProp ../../CoolProp/*.cpp
g++ -shared -o CoolProp.so *.o
sudo cp CoolProp.so `php-config --extension-dir`
rm *.o 
