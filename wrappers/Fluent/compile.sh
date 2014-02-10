#!/bin/bash

# Warning: this script deletes directory libudf before creating a new one

FLUENT_COMPILER=g++

HOME=`pwd`

# Copy CoolProp sources to a local folder called coolprop
cp ../CoolProp coolprop

cd coolprop
g++ -c CoolProp/*.cpp -ICoolProp -fPIC -D_GNU_SOURCE -ansi -O -Wall -DPTR_RESTRICT= 
cd $HOME

echo exit > exit.jou
fluent 3ddp -g -env -i exit.jou > FluentEnvironment.dat 2>&1
rm exit.jou
FLUENT_INC=`cat FluentEnvironment.dat | grep FLUENT_INC | sed 's/FLUENT[_A-Z]*=//'`
FLUENT_ARCH=`cat FluentEnvironment.dat | grep FLUENT_ARCH | sed 's/FLUENT[_A-Z]*=//'`
FLUENT_PROD_DIR=`cat FluentEnvironment.dat | grep FLUENT_PROD_DIR | sed 's/FLUENT[_A-Z]*=//'`
PATH_MAKEFILE1=$FLUENT_PROD_DIR"/src/makefile.udf"
PATH_MAKEFILE2=$FLUENT_PROD_DIR"/src/makefile.udf2"
PATH_LIBRARY="libudf/"$FLUENT_ARCH"/3ddp"
rm -rf libudf/
mkdir libudf
cp $PATH_MAKEFILE2 libudf/Makefile
mkdir libudf/src
cp $PATH_MAKEFILE1 libudf/src/makefile	
mkdir libudf/$FLUENT_ARCH
mkdir $PATH_LIBRARY
mkdir ${PATH_LIBRARY}"_host"
mkdir ${PATH_LIBRARY}"_node"
echo CSOURCES=UDF.c > ${PATH_LIBRARY}/user.udf 
echo HSOURCES= >> ${PATH_LIBRARY}/user.udf 
echo FLUENT_INC=$FLUENT_INC >> ${PATH_LIBRARY}/user.udf
echo GPU_SUPPORT=off >> ${PATH_LIBRARY}/user.udf
cp ${PATH_LIBRARY}/user.udf ${PATH_LIBRARY}_host/user.udf
cp ${PATH_LIBRARY}/user.udf ${PATH_LIBRARY}_node/user.udf
cp UDF.c libudf/src/UDF.c
cd libudf && make "FLUENT_ARCH=$FLUENT_ARCH" "SOURCES=UDF.c" "FLUENT_INC=$FLUENT_INC" "CC=$FLUENT_COMPILER" "CFLAGS_LNAMD64=-D_lnamd64 -D_GNU_SOURCE -fpic -shared -ansi -O -Wall -DPTR_RESTRICT= "
cd $HOME
cp coolprop/*.o ${PATH_LIBRARY}
cp coolprop/*.o ${PATH_LIBRARY}"_host"
cp coolprop/*.o ${PATH_LIBRARY}"_node"
cd ${PATH_LIBRARY} && g++ -shared -lm -ldl *.o  -o libudf.so
cd $HOME
cd ${PATH_LIBRARY}"_host" && g++ -shared -lm -ldl *.o  -o libudf.so
cd $HOME
cd ${PATH_LIBRARY}"_node" && g++ -shared -lm -ldl *.o  -o libudf.so

rm -rf coolprop
