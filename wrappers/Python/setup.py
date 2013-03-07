from __future__ import print_function
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

version = open(os.path.join('..','..','version.txt'),'r').read().strip()

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
                    f = open('../../CoolProp/svnrevision.h','r')
                    current_svn = f.read()
                    f.close()
                    if not current_svn.strip() == svnstring.strip():                
                        f = open('../../CoolProp/svnrevision.h','w')
                        f.write(svnstring)
                        f.close()
                    break
        except subprocess.CalledProcessError:
            pass
        
    svnrev_to_file()

    def version_to_file():
        string_for_file = 'char version [] ="{v:s}";'.format(v = version)
        f = open('../../CoolProp/version.h','r')
        current_version = f.read()
        f.close()
        if not current_version.strip() == string_for_file.strip():
            f = open('../../CoolProp/version.h','w')
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

    f = open('../../CoolProp/svnrevision.h','r')
    rev = f.read().strip().split('=')[1].strip(';').strip()
    f.close()
    svnstring = '__svnrevision__ ='+rev+'\n'
    
    for i in range(len(lines)-1,-1,-1):
        if lines[i].strip().startswith('#') or len(lines[i].strip())==0: 
            lines.pop(i)
    lines=['from __future__ import absolute_import\n']+['__version__=\''+version+'\'\n']+[svnstring]+lines
    fp=open(os.path.join('CoolProp','__init__.py'),'w')
    for line in lines:
        fp.write(line)
    fp.close()

    Sources = glob.glob(os.path.join('..','..','CoolProp','*.cpp'))
    
    ### Include folders for build
    include_dirs = [os.path.join('..','..','CoolProp'),get_python_inc(False)]

    def touch(fname):
        open(fname, 'a').close()
        os.utime(fname, None)

    #If library has been updated but the cython sources haven't been,
    #force cython to build by touching the cython sources
    cython_sources = [os.path.join('CoolProp','CoolProp.pyx')]
                      
##     if useStaticLib:
##         CC = new_compiler()
##         CPLibPath=CC.library_filename(os.path.join('lib','CoolProp'),lib_type='static')
##         if not newer_group(cython_sources, CPLibPath):
##             for source in cython_sources:
##                 print('touching',source)
##                 touch(source)
##         else:
##             print('no touching of cython sources needed')

    CoolProp_module = CyExtension('CoolProp.CoolProp',
                        [os.path.join('CoolProp','CoolProp.pyx')]+Sources,
                        include_dirs = include_dirs,
                        language='c++',
                        cython_c_in_temp = True
                        )
        
    param_constants_module = CyExtension('CoolProp.param_constants',
                            [os.path.join('CoolProp','param_constants.pyx')],
                            include_dirs = include_dirs,
                            language='c++',
                            cython_c_in_temp = True
                            )
                            
    phase_constants_module = CyExtension('CoolProp.phase_constants',
                            [os.path.join('CoolProp','phase_constants.pyx')],
                            include_dirs = include_dirs,
                            language='c++',
                            cython_c_in_temp = True
                            )
                            
    #Collect all the header files in the main folder into an include folder
    try:
        os.mkdir(os.path.join('CoolProp','include'))
    except:
        pass
        
    for header in glob.glob(os.path.join('..','..','CoolProp','*.h')):
        pth,fName = os.path.split(header)
        shutil.copy2(header,os.path.join('CoolProp','include',fName))
    
    setup (name = 'CoolProp',
           version = version, #look above for the definition of version variable - don't modify it here
           author = "Ian Bell",
           author_email='ian.h.bell@gmail.com',
           url='http://coolprop.sourceforge.net',
           description = """Open-source thermodynamic and transport properties database""",
           packages = ['CoolProp','CoolProp.Plots','CoolProp.tests','CoolProp.GUI'],
           ext_modules = [CoolProp_module,param_constants_module,phase_constants_module],
           package_dir = {'CoolProp':'CoolProp',},
           package_data = {'CoolProp':['State.pxd','CoolProp.pxd','param_constants.pxd','include/*.h']},
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
    
    #Clean up the 
    shutil.rmtree(os.path.join('CoolProp','include'), ignore_errors = True)
            
    touch('setup.py')
    