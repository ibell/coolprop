
from distutils.core import setup, Extension
import subprocess,shutil,os,sys
from Cython.Distutils import build_ext
from Cython.Distutils.extension import Extension as CyExtension

if len(sys.argv)==1:
    sys.argv+=['install','install']
    
badfiles = [os.path.join('CoolProp','__init__.pyc'),os.path.join('CoolProp','__init__.py')]
for file in badfiles:
    try:
        os.remove(file)
    except:
        pass
    
def availableFluids():
    try:
        FL = subprocess.check_output(['python','-c','import CoolProp; print CoolProp.CoolProp.FluidsList()']).rstrip()
        line = '__fluids__=' + str(FL.split(',')) +'\n' 
    except:
        line=''
    print line
    return line

version='1.5.0'

#This will automagically find all the fluid sources as long as they are in the right folders
#Pure fluids should all be in the src/purefluids folder relative to setup.py
#Pseudo-Pure fluids should all be in the src/pseudopurefluids folder relative to setup.py
import glob
purefluids=glob.glob(os.path.join('CoolProp','purefluids','*.cpp'))
pseudopurefluids=glob.glob(os.path.join('CoolProp','pseudopurefluids','*.cpp'))
others=glob.glob(os.path.join('CoolProp','*.cpp'))
Sources=purefluids+pseudopurefluids+others

#Build a list of all the available fluids - added to __init__.py
FluidsList=availableFluids()

#Unpack the __init__.py file template and add some things to the __init__ file
lines=open('__init__.py.template','r').readlines()
    
for i in range(len(lines)-1,-1,-1):
    if lines[i].strip().startswith('#') or len(lines[i].strip())==0: 
        lines.pop(i)
lines=[FluidsList]+['__version__=\''+version+'\'\n']+lines
fp=open(os.path.join('CoolProp','__init__.py'),'w')
for line in lines:
    fp.write(line)
fp.close()

## ==== start of SWIG code =======

swig_opts=['-c++','-python','-builtin']

"""
In this block of code, all the files that require SWIG are rebuilt on an as needed basis.  
If the target file X_wrap.cpp doesn't exist, or is older than the X.i file, SWIG
is called to rebuilt the file

swig_sources[(source,target)]

"""
swig_sources=[(os.path.join('CoolProp','CoolProp.i'),os.path.join('CoolProp','CoolProp_wrap.cpp')),
              (os.path.join('CoolProp','FloodProp.i'),os.path.join('CoolProp','FloodProp_wrap.cpp')),
              (os.path.join('CoolProp','HumidAirProp.i'),os.path.join('CoolProp','HumidAirProp_wrap.cpp'))]

for source,target in swig_sources:
    
    def rebuild_swig(source,target):
        swig_call=['swig']+swig_opts+['-o',target,source]
        print 'Swigging '+source+' to '+target+' ....'
        print swig_call
        subprocess.call(swig_call)
        
    #if the target doesn't exist or the wrapped C++ code is newer
    if not os.path.exists(target) or os.path.getmtime(source)>os.path.getmtime(target):
        rebuild_swig(source,target)
    else:
        print 'No SWIG required for '+source+' --> '+target+' (up-to-date)'

#from distutils.sysconfig import get_python_inc
#from distutils.ccompiler import new_compiler
#include_dirs = ['CoolProp',os.path.join('CoolProp','purefluids'),os.path.join('CoolProp','pseudopurefluids'),get_python_inc(False)]
#def StaticLibBuilder(sources):
#    CC = new_compiler()
#    objs = CC.compile(sources,'build_lib',[('COOLPROP_LIB',None)],include_dirs=include_dirs)
#    print objs
#    CC.create_static_lib(objs, 'CoolProp','lib')    
#StaticLibBuilder(Sources)            

## ==== end of SWIG code =======

CoolProp_module = Extension('CoolProp._CoolProp',
                           sources=[os.path.join('CoolProp','CoolProp_wrap.cpp')]+Sources,
                           #swig_opts=['-builtin']
                           swig_opts=['-builtin','-c++'],
                           #swig_opts=['-c++'],
                           include_dirs = ['CoolProp',os.path.join('CoolProp','purefluids'),os.path.join('CoolProp','pseudopurefluids')],
                           )

FloodProp_module = Extension('CoolProp._FloodProp',
                           sources=[os.path.join('CoolProp','FloodProp_wrap.cpp')]+Sources,
                           swig_opts=['-c++'],
                           include_dirs = ['CoolProp',os.path.join('CoolProp','purefluids'),os.path.join('CoolProp','pseudopurefluids')],
                           )

HASources = [
     os.path.join('CoolProp','HumidAirProp_wrap.cpp'),
     os.path.join('CoolProp','pseudopurefluids','Air.cpp'),
     os.path.join('CoolProp','purefluids','Water.cpp'),
     os.path.join('CoolProp','purefluids','R134a.cpp'),
     os.path.join('CoolProp','REFPROP.cpp'),
     os.path.join('CoolProp','Brine.cpp'),
     os.path.join('CoolProp','CoolProp.cpp'),
     os.path.join('CoolProp','CoolPropTools.cpp'),
     os.path.join('CoolProp','FluidClass.cpp'),
     os.path.join('CoolProp','Helmholtz.cpp'),
     os.path.join('CoolProp','HumAir.cpp'),
     os.path.join('CoolProp','Ice.cpp'),
     os.path.join('CoolProp','PengRobinson.cpp'),
     os.path.join('CoolProp','Solvers.cpp'),
     ]
HumidAirProp_module = Extension('CoolProp._HumidAirProp',
                           sources=HASources,
                           swig_opts=['-c++'],define_macros=[('ONLY_AIR_WATER',None)],
                           include_dirs = ['CoolProp',os.path.join('CoolProp','purefluids'),os.path.join('CoolProp','pseudopurefluids')],
                           )                     

#State_module = CyExtension('CoolProp.State',[os.path.join('CoolProp','State.pyx')],language='c++')                

State_module = CyExtension('CoolProp.State',[os.path.join('CoolProp','State.pyx')],language='c++',libraries=['CoolProp'],
                        library_dirs=['lib'],
                        include_dirs = ['CoolProp',os.path.join('CoolProp','purefluids'),os.path.join('CoolProp','pseudopurefluids')])
#                        
setup (name = 'CoolProp',
       version = version, #look above for the definition of version variable - don't modify it here
       author = "Ian Bell",
       author_email='ian.h.bell@gmail.com',
       url='http://coolprop.sourceforge.net',
       description = """Properties of pure fluids, pseudo-pure fluids and brines""",
       packages = ['CoolProp','CoolProp.Plots','CoolProp.tests'],
       ext_modules = [State_module,CoolProp_module,FloodProp_module,HumidAirProp_module],
       package_dir = {'CoolProp':'CoolProp',},
       cmdclass={'build_ext': build_ext}
       )

badfiles = [os.path.join('CoolProp','__init__.pyc'),os.path.join('CoolProp','__init__.py')]
for file in badfiles:
    try:
        os.remove(file)
    except:
        pass