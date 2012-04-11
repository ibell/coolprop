
from distutils.core import setup, Extension
import subprocess,shutil,os,sys

try:
    os.remove('CoolProp.pyc')
except:
    pass
shutil.copy2('__init__.py.template','__init__.py')

# If you add a fluid, update this list of fluids
FluidSources = ['R134a.cpp','R744.cpp','R290.cpp','R410A.cpp',
               'Brine.cpp','R32.cpp','R717.cpp','R404A.cpp','R407C.cpp',
               'R507A.cpp','Argon.cpp','Nitrogen.cpp','Water.cpp','Air.cpp']
                           
CoolProp_module = Extension('_CoolProp',
                           sources=['CoolProp.i', 'CoolProp.cpp','CoolPropTools.cpp','REFPROP.cpp']+FluidSources,
                           #swig_opts=['-builtin']
                           swig_opts=['-c++']
                           )

FloodProp_module = Extension('_FloodProp',
                           sources=['FloodProp.i', 'FloodProp.cpp','CoolProp.cpp','CoolPropTools.cpp','REFPROP.cpp']+FluidSources,
                           swig_opts=['-c++']
                           )
                           
HumidAirProp_module = Extension('_HumidAirProp',
                           sources=['HumidAirProp.i','HumAir.cpp','CoolProp.cpp','CoolPropTools.cpp','Ice.cpp','REFPROP.cpp']+FluidSources,
                           #swig_opts=['-builtin']
                           swig_opts=['-c++']
                           )                           

setup (name = 'CoolProp',
       version = '1.4.0',
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
