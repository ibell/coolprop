
from distutils.core import setup, Extension
import subprocess,shutil,os,sys

try:
    os.remove('CoolProp.pyc')
except:
    pass
shutil.copy2('__init__.py.template','__init__.py')

CoolProp_module = Extension('_CoolProp',
                           sources=['CoolProp.i', 'CoolProp.c',
                           'R134a.c','R744.c','R290.c','R410A.c',
                           'Brine.c','R32.c','R717.c','R404A.c','R407C.c',
                           'R507A.c','Argon.c','Nitrogen.c','Water.c'],
                           )

FloodProp_module = Extension('_FloodProp',
                           sources=['FloodProp.i', 'FloodProp.c',
                           'R134a.c','R744.c','R290.c','R410A.c',
                           'Brine.c','R32.c','R717.c','R404A.c','Nitrogen.c','Argon.c'],
                           )
                           
HumidAirProp_module = Extension('_HumidAirProp',
                           sources=['HumidAirProp.i','HumAir.c'],
                           )                           

setup (name = 'CoolProp',
       version = '1.2.2',
       author      = "Ian Bell",
       author_email='ian.h.bell@gmail.com',
       url='http://coolprop.sourceforge.net',
       description = """Properties of R134a, R744, R410A, R290, R717, R32, R404A and brines""",
       packages = ['CoolProp','CoolProp.Plots'],
       ext_package = 'CoolProp',
       ext_modules = [CoolProp_module,FloodProp_module,HumidAirProp_module],
       package_dir = {'CoolProp':'.'}
       )

if 'install' in sys.argv:
    os.remove('__init__.py')
    os.remove('CoolProp.py')
    os.remove('CoolProp_wrap.c')
    os.remove('FloodProp.py')
    os.remove('FloodProp_wrap.c')
    os.remove('HumidAirProp.py')
    os.remove('HumidAirProp_wrap.c')