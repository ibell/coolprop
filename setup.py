
from distutils.core import setup, Extension
import subprocess,shutil,os,sys

# Obtain the numpy include directory.  This logic works across numpy versions.
## import numpy
## try:
##     numpy_include = numpy.get_include()
## except AttributeError:
##     numpy_include = numpy.get_numpy_include()
numpy_include=[]

version='1.4.0'

#Unpack the __init__.py file and add the version number to the __init__ file
shutil.copy2('__init__.py.template','__init__.py')
lines=open('__init__.py','r').readlines()
lines=['__version__=\''+version+'\'\n']+lines
fp=open('__init__.py','w')
for line in lines:
    fp.write(line)
fp.close



import glob
#This will automagically find all the fluid sources as long as they are in the right folders
#Pure fluids should all be in the src/purefluids folder relative to setup.py
#Pseudo-Pure fluids should all be in the src/pseudopurefluids folder relative to setup.py

purefluids=glob.glob(os.path.join('src','purefluids','*.cpp'))
pseudopurefluids=glob.glob(os.path.join('src','pseudopurefluids','*.cpp'))
others=[]
FluidSources=purefluids+pseudopurefluids+others

CoolProp_module = Extension('_CoolProp',
                           sources=['CoolProp.i', 'CoolProp.cpp','CoolPropTools.cpp','REFPROP.cpp']+FluidSources,
                           #swig_opts=['-builtin']
                           swig_opts=['-builtin','-c++'],
                           #swig_opts=['-c++'],
                           include_dirs = [numpy_include],
                           )

FloodProp_module = Extension('_FloodProp',
                           sources=['FloodProp.i', 'FloodProp.cpp','CoolProp.cpp','CoolPropTools.cpp','REFPROP.cpp']+FluidSources,
                           swig_opts=['-c++']
                           #swig_opts=['-builtin']
                           )
                           
HumidAirProp_module = Extension('_HumidAirProp',
                           sources=['HumidAirProp.i','HumAir.cpp','CoolProp.cpp','CoolPropTools.cpp','Ice.cpp','REFPROP.cpp']+FluidSources,
                           #swig_opts=['-builtin']
                           swig_opts=['-c++']
                           )                           

setup (name = 'CoolProp',
       version = version,
       author      = "Ian Bell",
       author_email='ian.h.bell@gmail.com',
       url='http://coolprop.sourceforge.net',
       description = """Properties of R134a, R744, R410A, R290, R717, R32, R404A and brines""",
       packages = ['CoolProp','CoolProp.Plots'],
       ext_package = 'CoolProp',
       ext_modules = [CoolProp_module,FloodProp_module,HumidAirProp_module],
       package_dir = {'CoolProp':'.'}
       )

## Clean up the intermediate files that SWIG generates
if 'clean' in sys.argv:
    FileList=['__init__.py','CoolProp.py','CoolProp_wrap.cpp','FloodProp.py','FloodProp_wrap.cpp','HumidAirProp.py','HumidAirProp_wrap.cpp']
    for file in FileList:
        try:
            os.remove(file)
        except:
            print "Sorry, couldn't remove the file "+file
