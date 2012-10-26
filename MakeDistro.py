import subprocess,os,shutil

#These should be paths to python executables that you want want use to build versions of CoolProp
PYTHONVERSIONS=['_python\py27\python.exe',
                '_python\py27_x64\python.exe',
                ]

if not os.path.exists('_deps'):
    os.mkdir('_deps')
    
def remove_coolprop_lib():
    #Force the remove of the CoolProp library if one exists
    if os.path.exists(os.path.join('lib','CoolProp.lib')):
        os.remove(os.path.join('lib','CoolProp.lib'))
        
def InstallPrereqs():
    """ Get the requirements for CoolProp """
    #Collect the source for Cython and put in _deps/cython
    subprocess.call(['easy_install','-U','--editable','--build-directory','_deps','Cython'], shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    #Collect the source for setuptools and put in _deps/distribute
    subprocess.call(['easy_install','-U','--editable','--build-directory','_deps','setuptools'], shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        
    for python_install in PYTHONVERSIONS:
        for cwd in ['_deps/cython','_deps/distribute']:
            print subprocess.check_output([python_install, 'setup.py', 'install'], cwd = cwd)
    
def PYPI():
    subprocess.call(['python','setup.py','sdist','upload'])
    
def Source():
    python_install = PYTHONVERSIONS[0]
    print subprocess.check_output([python_install,'setup.py','sdist','--dist-dir=dist_temp/Python'],shell=True,cwd='.')

def DLL():
    """ Build a DLL using __stdcall calling convention """
    subprocess.call(['python','setup.py','--DLL'])
    #Collect the zip file and p
    try:
        os.makedirs(os.path.join('dist_temp','Excel and DLL'))
    except os.error:
        pass
        
    shutil.copy2(os.path.join('lib','CoolProp.dll'),os.path.join('dist_temp','Excel and DLL','CoolProp.dll'))
    shutil.copy2(os.path.join('CoolProp','CoolProp.h'),os.path.join('dist_temp','Excel and DLL','CoolProp.h'))
    shutil.copy2('DLLREADME.txt',os.path.join('dist_temp','Excel and DLL','README.txt'))
    shutil.copy2(os.path.join('Examples','CoolPropDLL.py'),os.path.join('dist_temp','Excel and DLL','CoolPropDLL.py'))
    shutil.copy2(os.path.join('wrappers','Excel','CoolProp.xlam'),os.path.join('dist_temp','Excel and DLL','CoolProp.xlam'))
    shutil.copy2(os.path.join('wrappers','Excel','TestExcel.xlsx'),os.path.join('dist_temp','Excel and DLL','TestExcel.xlsx'))
    shutil.copy2(os.path.join('wrappers','Excel','ExcelInstructions.txt'),os.path.join('dist_temp','Excel and DLL','ExcelInstructions.txt'))
    
def Octave():
    try:
        os.makedirs(os.path.join('dist_temp','Octave'))
    except os.error: pass
        
    subprocess.check_output(['OctaveBuilder.bat'],shell=True,cwd=os.path.join('wrappers','Octave'))
    shutil.copy2(os.path.join('wrappers','Octave','CoolProp.oct'),os.path.join('dist_temp','Octave','CoolProp.oct'))
    shutil.copy2(os.path.join('wrappers','Octave','README.txt'),os.path.join('dist_temp','Octave','README.txt'))
    
def MATLAB():
    try:
        os.makedirs(os.path.join('dist_temp','MATLAB'))
    except os.error: pass
        
    process = subprocess.Popen(['C:\\MATLAB_32bit\\bin\\matlab','-r','MATLABBuilder'],shell=True,cwd=os.path.join('wrappers','MATLAB'))
    process.wait()
    process = subprocess.Popen(['matlab','-nojvm','-nodesktop','-nosplash','-r','MATLABBuilder'],shell=True,cwd=os.path.join('wrappers','MATLAB'))
    process.wait()
    shutil.copy2(os.path.join('wrappers','MATLAB','Props.mexw64'),os.path.join('dist_temp','MATLAB','Props.mexw64'))
    shutil.copy2(os.path.join('wrappers','MATLAB','HAProps.mexw64'),os.path.join('dist_temp','MATLAB','HAProps.mexw64'))
    shutil.copy2(os.path.join('wrappers','MATLAB','Props.mexw32'),os.path.join('dist_temp','MATLAB','Props.mexw32'))
    shutil.copy2(os.path.join('wrappers','MATLAB','HAProps.mexw32'),os.path.join('dist_temp','MATLAB','HAProps.mexw32'))
    shutil.copy2(os.path.join('wrappers','MATLAB','README.txt'),os.path.join('dist_temp','MATLAB','README.txt'))
    shutil.copy2(os.path.join('wrappers','MATLAB','MATLAB_sample.m'),os.path.join('dist_temp','MATLAB','MATLAB_sample.m'))
    
def PythonInstallers():
    
    for python_install in PYTHONVERSIONS:
        remove_coolprop_lib()
        print subprocess.check_output([python_install,'setup.py','bdist','--format=wininst','--dist-dir=dist_temp/Python'],shell=True,cwd='.')
    #For good measure, clean up after ourselves
    remove_coolprop_lib()
    
def UploadSourceForge():
    #Rename folder to version number
    from setup import version
    try:
        shutil.copytree('dist_temp',version)
    except WindowsError: pass
    
    call_str = ['pscp','README.txt','ibell,coolprop@frs.sf.net:/home/pfs/project/c/co/coolprop/CoolProp/']
    print 'Calling: '+' '.join(call_str)
    print subprocess.check_output(call_str,shell=True)
    
    call_str = ['pscp','-r','-v',version,'ibell,coolprop@frs.sf.net:/home/pfs/project/c/co/coolprop/CoolProp/']
    print 'Calling: '+' '.join(call_str)
    print subprocess.check_output(call_str,shell=True)
    
    os.remove('dist_temp')
    os.remove(version)
    
if __name__=='__main__':
##     InstallPrereqs()
##     PythonInstallers()
##     DLL()
##     Source()
##     Octave()
    MATLAB()
##     PYPI()
##     UploadSourceForge()