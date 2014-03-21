#!/bin/bash
#
#cd ../../CoolProp
NAME=CoolProp
SRCD="../../${NAME}"
JDK_PATH=''
find_JDK(){
    if [ -n "$JAVA_HOME" ]; then
        echo "Find \$JAVA_HOME variable, use this JDK path: $JAVA_HOME"
        JDK_PATH="$JAVA_HOME"
        return 0
    else
        for jvmdir in /usr/lib/jvm/java-7-openjdk-*
        do
            if [ -d "${jvmdir}" -a "${jvmdir}" != "/usr/lib/jvm/java-7-openjdk-common" ]
            then
                OPENJDKS=$jvmdir
                echo "Find openjdk7 path, use this JDK path: $OPENJDKS"
                JDK_PATH="$OPENJDKS"
                return 0
            fi
        done
        for jvmdir in /usr/lib/jvm/java-6-openjdk-*
        do
            if [ -d "${jvmdir}" -a "${jvmdir}" != "/usr/lib/jvm/java-6-openjdk-common" ]
            then
                OPENJDKS="${OPENJDKS} ${jvmdir}"
                echo "Find openjdk6 path, use this JDK path: $OPENJDKS"
                JDK_PATH="$OPENJDKS"
                return 0
            fi
        done
    fi
    echo "Not find \$JAVA_HOME variable and openjdk. Please install a JDK and set \$JAVA_HOME!"
    exit 1
}

find_JDK
JDIR="$JDK_PATH/include"
JAVA="$JDK_PATH/bin/java"
JAVAC="$JDK_PATH/bin/javac"

INCL="-I${JDIR} -I${SRCD}"
#
swig -c++ -java -outcurrentdir ${SRCD}/${NAME}.i
g++ -fpic  ${INCL} -c ${NAME}_wrap.cxx
g++ -fpic  ${INCL} -c ${SRCD}/*.cpp
g++ -shared *.o -o lib${NAME}.so
"$JAVAC" *.java
##rm *.o
"$JAVA" -Djava.library.path=. Example
