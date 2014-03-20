#!/bin/bash
#
#cd ../../CoolProp
NAME=CoolProp
SRCD="../../${NAME}"
JDIR=/usr/lib/jvm/java-7-openjdk-i386/include
#
INCL="-I${JDIR} -I${SRCD}"
#
swig -c++ -java -outcurrentdir ${SRCD}/${NAME}.i
g++ -fpic  ${INCL} -c ${NAME}_wrap.cxx
g++ -fpic  ${INCL} -c ${SRCD}/*.cpp
g++ -shared *.o -o lib${NAME}.so
javac *.java
#rm *.o
java -Djava.library.path=. Example