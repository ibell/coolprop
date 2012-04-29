from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import sys,os
from glob import glob
if len(sys.argv)==1:
    sys.argv+=['build_ext','--inplace']

CPRoot=os.path.join('..','..')
purefluids      =glob(os.path.join(CPRoot,'src','purefluids','*.cpp'))
pseudopurefluids=glob(os.path.join(CPRoot,'src','pseudopurefluids','*.cpp'))
others          =glob(os.path.join(CPRoot,'src','*.cpp'))
Sources=["CoolPropCython.pyx"]+purefluids+pseudopurefluids+others

CPinclude=[os.path.join(CPRoot,'src'),os.path.join(CPRoot,'src','purefluids'),
    os.path.join(CPRoot,'src','pseudopurefluids')]

ext_modules = [Extension("CoolPropCy",Sources,language="c++",include_dirs=CPinclude)]

#Get the numpy include folder
import numpy
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

setup(
    name = 'Hello world app',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
    include_dirs = [numpy_include]+CPinclude,
)

import subprocess
subprocess.call("python -c 'from main import integrate'",shell=True)