#!/bin/bash
#
#cd ../../CoolProp
NAME=CoolProp
SRCD="../../${NAME}"
JDK_PATH=''
find_JDK(){
    if [ -n "$JAVA_HOME" ]; then
        echo -e '\033[1;92m' "Find \$JAVA_HOME variable, use this JDK path: $JAVA_HOME" '\033[0m'
        JDK_PATH="$JAVA_HOME"
        return 0
    else
        for jvmdir in /usr/lib/jvm/java-7-openjdk-*
        do
            if [ -d "${jvmdir}" -a "${jvmdir}" != "/usr/lib/jvm/java-7-openjdk-common" ]
            then
                OPENJDKS=$jvmdir
                echo -e '\033[1;92m' "Find openjdk7 path, use this JDK path: $OPENJDKS" '\033[0m'
                JDK_PATH="$OPENJDKS"
                return 0
            fi
        done
        for jvmdir in /usr/lib/jvm/java-6-openjdk-*
        do
            if [ -d "${jvmdir}" -a "${jvmdir}" != "/usr/lib/jvm/java-6-openjdk-common" ]
            then
                OPENJDKS="${OPENJDKS} ${jvmdir}"
                echo -e '\033[1;92m' "Find openjdk6 path, use this JDK path: $OPENJDKS" '\033[0m'
                JDK_PATH="$OPENJDKS"
                return 0
            fi
        done
    fi
    echo -e '\033[1;91m' "Not find \$JAVA_HOME variable and openjdk. Please install a JDK and set \$JAVA_HOME!" '\033[0m'
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
#rm *.o
"$JAVA" -Djava.library.path=. Example

# a simple way of installing a system-wide shared library.
if [ $? -eq 0 ]; then
    echo ''
    echo -e '\033[1;92m' "Build successful. Will install a system-wide shared library. Need sudo privilege." '\033[0m'
    echo ''
    cd ../SharedLibrary
    make
    sudo make install
else
    echo -e '\033[1;91m' "Build failed. Please check the dependencies and run again." '\033[0m'
fi
