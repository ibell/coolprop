. ${MATLAB}/bin/mexopts.sh

CFLAGS='-D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread'
CXXFLAGS='-D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread'
CLIBS="-ldl ${CLIBS}"
CXXLIBS="-ldl ${CXXLIBS}"