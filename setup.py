
from distutils.core import setup, Extension
import subprocess,shutil,os,sys,glob
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from Cython.Distutils.extension import Extension as CyExtension
from distutils.sysconfig import get_python_inc
from distutils.ccompiler import new_compiler 
from distutils.dep_util import newer_group

import Cython
#This will generate HTML to show where there are still pythonic bits hiding out
Cython.Compiler.Options.annotate = True

## If the file is run directly without any parameters, build and install
if len(sys.argv)==1:
    #sys.argv+=['build_ext','--inplace']
    sys.argv+=['install']
    
if '--DLL' in sys.argv:
    packDLL = True
    sys.argv.pop(sys.argv.index('--DLL'))
else:
    packDLL = False
    
if '--no-static-lib' in sys.argv or 'sdist' in sys.argv:
    useStaticLib = False
    if '--no-static-lib' in sys.argv:
        sys.argv.pop(sys.argv.index('--no-static-lib'))
else:
    useStaticLib=True
    

if '--update-fluids' in sys.argv:
    FL = subprocess.Popen(['python','-c','import CoolProp; print CoolProp.CoolProp.FluidsList()'], stdout=subprocess.PIPE).communicate()[0].rstrip()
    fp = open('fluids.txt','w')
    fp.write(FL)
    fp.close()
    sys.argv.pop(sys.argv.index('--update-fluids'))
    
badfiles = [os.path.join('CoolProp','__init__.pyc'),os.path.join('CoolProp','__init__.py')]
for file in badfiles:
    try:
        os.remove(file)
    except:
        pass
    
def availableFluids():
    #First see if there is a file called fluids.txt - a comma-separated list of fluids
    if os.path.exists('fluids.txt'):
        FL = open('fluids.txt','r').read()
        line = '__fluids__=' + str(FL.split(',')) +'\n'
        print line
        return line
        
    try:
        ## Can't load the CoolProp module in this file because then you can't install the files :(
        ## Run a subprocess to do it instead
        FL = subprocess.Popen(['python','-c','import CoolProp; print CoolProp.CoolProp.FluidsList()'], stdout=subprocess.PIPE).communicate()[0].rstrip()
        ## Take the comma,separated string and turn it into a Python-friendly list for writing to __init__.py file
        line = '__fluids__=' + str(FL.split(',')) +'\n' 
    except:
        line=''
    print line
    return line

version='2.0.4'

#########################
## __init__.py builder ##
#########################

#Build a list of all the available fluids - added to __init__.py
FluidsList=availableFluids()

#Unpack the __init__.py file template and add some things to the __init__ file
lines=open('__init__.py.template','r').readlines()

import subprocess
try:
    subprocess.Popen(['svn','update'], stdout=subprocess.PIPE).communicate()
    SVNInfo = subprocess.Popen(['svn','info'], stdout=subprocess.PIPE).communicate()[0].split('\n')
    for line in SVNInfo:
        if line.startswith('Revision'):
            svnstring='__svnrevision__ = '+line.strip().split(':')[1].strip()+'\n'
            break
    svnstring=''
except:
    svnstring=''
    
for i in range(len(lines)-1,-1,-1):
    if lines[i].strip().startswith('#') or len(lines[i].strip())==0: 
        lines.pop(i)
lines=[FluidsList]+['__version__=\''+version+'\'\n']+[svnstring]+lines
fp=open(os.path.join('CoolProp','__init__.py'),'w')
for line in lines:
    fp.write(line)
fp.close()

######################################
######################################
######################################
##### SWIG!! SWIG!!!! SWIG!!! ########
######################################
##       start of SWIG code         ## 
######################################
swig_opts=['-c++','-python']

"""
In this block of code, all the files that require SWIG are rebuilt on an as needed basis.  
If the target file X_wrap.cpp doesn't exist, or is older than the X.i file, SWIG
is called to rebuilt the file.  Does not check the c++ headers.  Delete the _wrap.cpp to force build

swig_sources[(source,target)]

"""
swig_sources=[(os.path.join('CoolProp','FloodProp.i'),os.path.join('CoolProp','FloodProp_wrap.cpp')),
              (os.path.join('CoolProp','HumidAirProp.i'),os.path.join('CoolProp','HumidAirProp_wrap.cpp'))]

for source,target in swig_sources:
    header = source.rsplit('.',1)[0]+'.h'
    def rebuild_swig(source,target):
        swig_call=['swig']+swig_opts+['-o',target,source]
        print 'Swigging '+source+' to '+target+' ....'
        subprocess.call(swig_call)
        
    #if the target doesn't exist or the wrapped C++ code is newer
    if (not os.path.exists(target) 
        or os.path.getmtime(source)>os.path.getmtime(target) 
        or os.path.getmtime(header)>os.path.getmtime(target)
        ):
        
        rebuild_swig(source,target)
    else:
        print 'No SWIG required for '+source+' --> '+target+' (up-to-date)'

################################
##  CoolProp Static Library ####
################################
## Build a static library of CoolProp

#This will automagically find all the fluid sources as long as they are in the right folders
#Pure fluids should all be in the CoolProp/purefluids folder relative to setup.py
#Pseudo-Pure fluids should all be in the CoolProp/pseudopurefluids folder relative to setup.py
purefluids=glob.glob(os.path.join('CoolProp','purefluids','*.cpp'))
pseudopurefluids=glob.glob(os.path.join('CoolProp','pseudopurefluids','*.cpp'))
CPCore=['CoolProp.cpp','CoolPropTools.cpp','FloodProp.cpp','FluidClass.cpp',
        'Helmholtz.cpp','PengRobinson.cpp','REFPROP.cpp','Solvers.cpp',
        'Brine.cpp']
others=[os.path.join('CoolProp',f) for f in CPCore]
Sources=purefluids+pseudopurefluids+others

### Include folders for build
include_dirs = ['CoolProp',os.path.join('CoolProp','purefluids'),os.path.join('CoolProp','pseudopurefluids'),get_python_inc(False)]

def StaticLibBuilder(sources,LibName='CoolProp',build_path='build_lib',lib_path='lib',force=False,DLL=True,StaticLib=True):
    CC = new_compiler(verbose=True)
    #Default to not build
    buildCPLib=False
    # The full path to the library to be built
    CPLibPath=CC.library_filename(os.path.join(lib_path,LibName),lib_type='static')    
    
    if sys.platform.startswith('linux'):
        extra_compile_args=['-fPIC']
        MACROS = None
    elif sys.platform.startswith('darwin'):
        extra_compile_args=['']
        MACROS = None
    else:
        extra_compile_args=['/EHsc']
        MACROS = [('COOLPROP_LIB',None)]
            
    if not os.path.exists(build_path) or not os.path.exists(lib_path):
        if not os.path.exists(build_path): os.mkdir(build_path)
        if not os.path.exists(lib_path): os.mkdir(lib_path)
        #Force rebuild of all
        buildCPLib=True
    elif newer_group(sources,CPLibPath):
        buildCPLib=True
    
    if buildCPLib or DLL:
        objs=CC.compile(sources,build_path,MACROS,include_dirs=include_dirs,extra_postargs=extra_compile_args)
        CC.create_static_lib(objs, LibName,lib_path,target_lang='c++')
        print 'Built the static library in',CPLibPath
        if DLL:
            CC.link_shared_lib(objs, LibName,lib_path,target_lang='c++')
            print 'Built the shared library in',os.path.join(lib_path,LibName)
    else:
        print 'No build of CoolProp static library required.'
  
if useStaticLib or packDLL:
    StaticLibBuilder(Sources,DLL=packDLL)
    
if packDLL:
    ZIPfilePath = 'dist/CoolPropDLL-'+version+'.zip'
    from zipfile import ZipFile
    with ZipFile(ZIPfilePath,'w') as z:
        z.write(os.path.join('lib','CoolProp.dll'),arcname='CoolProp.dll')
        z.write(os.path.join('CoolProp','CoolProp.h'),arcname='CoolProp.h')
        z.write(os.path.join('Examples','CoolPropDLL.py'),arcname='CoolPropDLL.py')
        z.write('DLLREADME.txt',arcname='README.txt')
    print 'DLL file has been packed into file', ZIPfilePath

print 'UseStaticLib is',useStaticLib
##Now come in and build the modules themselves
    
if useStaticLib==True:
    FloodProp_module = Extension('CoolProp._FloodProp',
                           sources=[os.path.join('CoolProp','FloodProp_wrap.cpp')],
                           include_dirs = include_dirs,language='c++',
                           libraries=['CoolProp'],library_dirs=['lib']
                           )
else:
    FloodProp_module = Extension('CoolProp._FloodProp',
                           sources=[os.path.join('CoolProp','FloodProp_wrap.cpp')]+Sources,
                           include_dirs = include_dirs,language='c++'
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
     os.path.join('CoolProp','HumidAirProp.cpp'),
     os.path.join('CoolProp','Ice.cpp'),
     os.path.join('CoolProp','PengRobinson.cpp'),
     os.path.join('CoolProp','Solvers.cpp'),
     ]
     
HumidAirProp_module = Extension('CoolProp._HumidAirProp',
                        sources=HASources,
                        define_macros=[('ONLY_AIR_WATER',None)],
                        include_dirs = include_dirs,
                        language='c++'
                        )                            

if useStaticLib==True:
    State_module = CyExtension('CoolProp.State',
                        [os.path.join('CoolProp','State.pyx')],
                        include_dirs = include_dirs,language='c++',
                        libraries=['CoolProp'],
                        library_dirs=['lib']
                        )
else:
    State_module = CyExtension('CoolProp.State',
                        [os.path.join('CoolProp','State.pyx')]+Sources,
                        include_dirs = include_dirs,
                        language='c++'
                        )

if useStaticLib==True:
    CoolProp2_module = CyExtension('CoolProp.CoolProp',
                        [os.path.join('CoolProp','CoolProp.pyx')],
                        include_dirs = include_dirs,
                        language='c++',
                        libraries=['CoolProp'],
                        library_dirs=['lib'],
                        cython_c_in_temp = True
                        )
else:
    CoolProp2_module = CyExtension('CoolProp.CoolProp',
                        [os.path.join('CoolProp','CoolProp.pyx')]+Sources,
                        include_dirs = include_dirs,
                        language='c++'
                        )
                        
setup (name = 'CoolProp',
       version = version, #look above for the definition of version variable - don't modify it here
       author = "Ian Bell",
       author_email='ian.h.bell@gmail.com',
       url='http://coolprop.sourceforge.net',
       description = """Open-source thermodynamic and transport properties database""",
       packages = ['CoolProp','CoolProp.Plots','CoolProp.tests'],
       ext_modules = [FloodProp_module,HumidAirProp_module,State_module,CoolProp2_module],
       package_dir = {'CoolProp':'CoolProp',},
       package_data = {'CoolProp':['State.pxd']},
       cmdclass={'build_ext': build_ext},
       
       classifiers = [
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules"
        ],

       )

badfiles = [os.path.join('CoolProp','__init__.pyc'),os.path.join('CoolProp','__init__.py')]
for file in badfiles:
    try:
        os.remove(file)
    except:
        pass