
from distutils.core import setup, Extension
import subprocess,shutil,os,sys

try:
    os.remove('CoolProp.pyc')
except:
    pass
shutil.copy2('__init__.py.template','__init__.py')

# If you add a fluid, update this list of fluids
FluidSources = ['R134a.c','R744.c','R290.c','R410A.c',
               'Brine.c','R32.c','R717.c','R404A.c','R407C.c',
               'R507A.c','Argon.c','Nitrogen.c','Water.c','Air.c']
                           
CoolProp_module = Extension('_CoolProp',
                           sources=['CoolProp.i', 'CoolProp.c']+FluidSources,
                           swig_opts=['-builtin']
                           )

FloodProp_module = Extension('_FloodProp',
                           sources=['FloodProp.i', 'FloodProp.c',
                           'R134a.c','R744.c','R290.c','R410A.c',
                           'Brine.c','R32.c','R717.c','R404A.c','Nitrogen.c','Argon.c'],
                           )
                           
HumidAirProp_module = Extension('_HumidAirProp',
                           sources=['HumidAirProp.i','HumAir.c','CoolProp.c','Ice.cpp','SolverFunctions.c']+FluidSources,
                           swig_opts=['-builtin']
                           )                           

setup (name = 'CoolProp',
       version = '1.3.2',
       author      = "Ian Bell",
       author_email='ian.h.bell@gmail.com',
       url='http://coolprop.sourceforge.net',
       description = """Properties of R134a, R744, R410A, R290, R717, R32, R404A and brines""",
       packages = ['CoolProp','CoolProp.Plots'],
       ext_package = 'CoolProp',
       ext_modules = [CoolProp_module,FloodProp_module,HumidAirProp_module], #PUT ME BACK!!!
       package_dir = {'CoolProp':'.'}
       )

## Clean up the intermediate files that SWIG generates
if 'clean' in sys.argv:
    FileList=['__init__.py','CoolProp.py','CoolProp_wrap.c','FloodProp.py','FloodProp_wrap.c','HumidAirProp.py','HumidAirProp_wrap.c']
    for file in FileList:
        try:
            os.remove(file)
        except:
            print "Sorry, couldn't remove the file "+file