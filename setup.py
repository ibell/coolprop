
#Check for cython >= 0.17 due to the use of STL containers
try:
    import Cython
except ImportError:
    raise ImportError("Cython not found")
major,minor = Cython.__version__.split('.')[0:2]
#Convert major to integer
major = int(major)
iEnd = 0
while minor[iEnd].isdigit():
    iEnd+=1
    if iEnd == len(minor):
        break

minor = int(minor[0:iEnd])
if not(major > 0 or minor >= 17):
    raise ImportError('Cython version >= 0.17 required due to the use of STL wrappers.  Please update your version of cython')
    
from distutils.core import setup, Extension
import subprocess,shutil,os,sys,glob
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from Cython.Distutils.extension import Extension as CyExtension
from distutils.sysconfig import get_python_inc
from distutils.ccompiler import new_compiler 
from distutils.dep_util import newer_group

#This will generate HTML to show where there are still pythonic bits hiding out
Cython.Compiler.Options.annotate = True

version = '2.5'

if __name__=='__main__':

    def svnrev_to_file():
        """
        If a svn repo, use subversion to update the file in revision
        """
        try:
            subprocess.call(['svn','update'], shell = True)
            p = subprocess.Popen(['svn','info'], 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()
            for line in str(stdout).split('\n'):
                if line.startswith('Revision'):
                    rev = line.split(':')[1].strip()
                    svnstring = 'long svnrevision = '+rev+';'
                    #Check if it is different than the current version
                    f = open('CoolProp/svnrevision.h','r')
                    current_svn = f.read()
                    f.close()
                    if not current_svn.strip() == svnstring.strip():                
                        f = open('CoolProp/svnrevision.h','w')
                        f.write(svnstring)
                        f.close()
                    break
        except subprocess.CalledProcessError:
            pass
        
    svnrev_to_file()

    def version_to_file():
        string_for_file = 'char version [] ="{v:s}";'.format(v = version)
        f = open('CoolProp/version.h','r')
        current_version = f.read()
        f.close()
        if not current_version.strip() == string_for_file.strip():
            f = open('CoolProp/version.h','w')
            f.write(string_for_file)
            f.close()
        
    version_to_file()
    
    #This will generate HTML to show where there are still pythonic bits hiding out
    Cython.Compiler.Options.annotate = True

    ## If the file is run directly without any parameters, build and install
    if len(sys.argv)==1:
        #sys.argv+=['build_ext','--inplace']
        #sys.argv+=['build','--compiler=mingw32','install']
        sys.argv+=['install']
        
    if '--DLL' in sys.argv:
        packDLL = True
        sys.argv.pop(sys.argv.index('--DLL'))
    else:
        packDLL = False
        
    if '--no-static-lib' in sys.argv or 'sdist' in sys.argv:
        useStaticLib = False
        if '--no-static-lib' in sys.argv:
            sys.argv.remove('--no-static-lib')
    else:
        useStaticLib=True
        
    badfiles = [os.path.join('CoolProp','__init__.pyc'),os.path.join('CoolProp','__init__.py')]
    for file in badfiles:
        try:
            os.remove(file)
        except:
            pass

    #########################
    ## __init__.py builder ##
    #########################

    #Unpack the __init__.py file template and add some things to the __init__ file
    lines=open('__init__.py.template','r').readlines()

    f = open('CoolProp/svnrevision.h','r')
    rev = f.read().strip().split('=')[1].strip(';').strip()
    f.close()
    svnstring = '__svnrevision__ ='+rev+'\n'
        
    for i in range(len(lines)-1,-1,-1):
        if lines[i].strip().startswith('#') or len(lines[i].strip())==0: 
            lines.pop(i)
    lines=['__version__=\''+version+'\'\n']+[svnstring]+lines
    fp=open(os.path.join('CoolProp','__init__.py'),'w')
    for line in lines:
        fp.write(line)
    fp.close()

    ################################
    ##  CoolProp Static Library ####
    ################################
    ## Build a static library of CoolProp

    #This will automagically find all the fluid sources as long as they are in the right folders
    #Pure fluids should all be in the CoolProp/purefluids folder relative to setup.py
    #Pseudo-Pure fluids should all be in the CoolProp/pseudopurefluids folder relative to setup.py
    purefluids=glob.glob(os.path.join('CoolProp','purefluids','*.cpp'))
    pseudopurefluids=glob.glob(os.path.join('CoolProp','pseudopurefluids','*.cpp'))
    others = glob.glob(os.path.join('CoolProp','*.cpp'))
    
    #Remove the _wrap files from the build
    for f in glob.glob(os.path.join('CoolProp','*_wrap.cpp')):
        others.remove(f)
    
    Sources=purefluids + pseudopurefluids + others
         
    ### Include folders for build
    include_dirs = ['CoolProp',
                    os.path.join('CoolProp','purefluids'),
                    os.path.join('CoolProp','pseudopurefluids'),
                    get_python_inc(False)]

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
            if DLL:
                MACROS = [('COOLPROP_LIB',None)]
            else:
                MACROS = None
                
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
            print('Built the static library in',CPLibPath)
            if DLL:
                CC.link_shared_lib(objs, LibName,lib_path,target_lang='c++')
                print('Built the shared library in',os.path.join(lib_path,LibName))
        else:
            print('No build of CoolProp static library required.')
      
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
        print('DLL file has been packed into file', ZIPfilePath)
        quit()

    print('UseStaticLib is',useStaticLib)
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

    def touch(fname):
        open(fname, 'a').close()
        os.utime(fname, None)

    #If library has been updated but the cython sources haven't been,
    #force cython to build by touching the cython sources
    cython_sources = [os.path.join('CoolProp','State.pyx'),
                      os.path.join('CoolProp','CoolProp.pyx')]
                      
    if useStaticLib:
        CC = new_compiler()
        CPLibPath=CC.library_filename(os.path.join('lib','CoolProp'),lib_type='static')
        if not newer_group(cython_sources, CPLibPath):
            for source in cython_sources:
                print('touching',source)
                touch(source)
        else:
            print('no touching of cython sources needed')

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
                            language='c++',
                            cython_c_in_temp = True
                            )
                            
    setup (name = 'CoolProp',
           version = version, #look above for the definition of version variable - don't modify it here
           author = "Ian Bell",
           author_email='ian.h.bell@gmail.com',
           url='http://coolprop.sourceforge.net',
           description = """Open-source thermodynamic and transport properties database""",
           packages = ['CoolProp','CoolProp.Plots','CoolProp.tests','CoolProp.GUI'],
           ext_modules = [CoolProp2_module],
           package_dir = {'CoolProp':'CoolProp',},
           package_data = {'CoolProp':['State.pxd','CoolProp.pxd']},
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

    badfiles = [os.path.join('CoolProp','__init__.pyc'),
                os.path.join('CoolProp','__init__.py'),
                os.path.join('CoolProp','HumidAirProp_wrap.cpp'),
                os.path.join('CoolProp','FloodProp_wrap.cpp'),
                os.path.join('CoolProp','State.cpp')
                ]
    for file in badfiles:
        try:
            os.remove(file)
        except:
            pass
            
    touch('setup.py')