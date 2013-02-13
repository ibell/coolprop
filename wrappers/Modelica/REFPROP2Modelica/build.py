
import subprocess,shutil,os,sys,glob
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from Cython.Distutils.extension import Extension as CyExtension
from distutils.sysconfig import get_python_inc
from distutils.ccompiler import new_compiler 
from distutils.dep_util import newer_group

#The path to the root of the source - two levels up
CProot = os.path.join('..','..')

# Include folders for build
include_dirs = [os.path.join(CProot,'CoolProp'),
                os.path.join(CProot,'CoolProp','purefluids'),
                os.path.join(CProot,'CoolProp','pseudopurefluids'),
                get_python_inc(False)]
                
#Instantiate the compiler
CC = new_compiler(verbose=True)

#Compile the wrapper file to an object file
objs = CC.compile(['REFPROP_wrapper.cpp'],'.',include_dirs=include_dirs)

# Collect all the compiled object files that are in the build_lib folder 
# relative to the root of the CoolProp source
for root, dirs, files in os.walk(os.path.join(CProot,'build_lib')):
    for file in files:
        objs.append(os.path.join(root,file))

#Create a static library
CC.create_static_lib(objs, 'REFPROP_wrapper','.',target_lang='c++')