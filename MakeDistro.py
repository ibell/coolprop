import subprocess,os,shutil

#These should be paths to python executables that you want want use to build versions of CoolProp
PYTHONVERSIONS=['python.exe', #This is python 2.7 on my computer
                'c:\\python\\py27_x64\\python.exe',
                'c:\\python\\py32\\python.exe',
                'c:\\python\\py32_x64\\python.exe',
                'c:\\python\\py33\\python.exe',
                'c:\\python\\py33_x64\\python.exe',
                ]

if not os.path.exists('_deps'):
    os.mkdir('_deps')
        
def InstallPrereqs():
    """ Get the requirements for CoolProp """
    #Collect the source for Cython and put in _deps/cython-master
    import urllib,zipfile
    print 'getting cython sources'
    urllib.urlretrieve('https://github.com/cython/cython/archive/master.zip', filename = 'master.zip')
    with zipfile.ZipFile('master.zip', 'r') as myzip:
        myzip.extractall(path='_deps')
    os.remove('master.zip')
    for python_install in PYTHONVERSIONS:
        for cwd in ['_deps/cython-master']:
            print subprocess.check_output([python_install, 'setup.py', 'install'], cwd = cwd)
            
    
def PYPI():
    subprocess.call(['python','setup.py','sdist','upload'],cwd=os.path.join('wrappers','Python'))
    
def Source():
    print subprocess.check_output(['python','setup.py','sdist','--dist-dir=../../dist_temp/Python'],shell=True,cwd=os.path.join('wrappers','Python'))

def DLL_and_Excel():
    """ Build a DLL using __stdcall calling convention """
    subprocess.check_output(['BuildDLL'],shell=True,cwd=os.path.join('wrappers','Excel'))
    subprocess.check_output(['BuildDLLx64'],shell=True,cwd=os.path.join('wrappers','Excel'))
    #Collect the zip file and p
    try:
        os.makedirs(os.path.join('dist_temp','Excel and DLL'))
    except os.error:
        pass
    
    shutil.copy2(os.path.join('CoolProp','CoolProp.h'),os.path.join('dist_temp','Excel and DLL','CoolProp.h'))    
    shutil.copy2(os.path.join('wrappers','Excel','CoolProp.dll'),os.path.join('dist_temp','Excel and DLL','CoolProp.dll'))
    shutil.copy2(os.path.join('wrappers','Excel','CoolProp_x64.dll'),os.path.join('dist_temp','Excel and DLL','CoolProp_x64.dll'))
    shutil.copy2(os.path.join('wrappers','Excel','CoolProp.xlam'),os.path.join('dist_temp','Excel and DLL','CoolProp.xlam'))
    shutil.copy2(os.path.join('wrappers','Excel','CoolProp.xla'),os.path.join('dist_temp','Excel and DLL','CoolProp.xla'))
    shutil.copy2(os.path.join('wrappers','Excel','TestExcel.xlsx'),os.path.join('dist_temp','Excel and DLL','TestExcel.xlsx'))
    shutil.copy2(os.path.join('wrappers','Excel','README.rst'),os.path.join('dist_temp','Excel and DLL','README.rst'))
    
def Octave():
    try:
        os.makedirs(os.path.join('dist_temp','Octave'))
        os.makedirs(os.path.join('dist_temp','Octave','3.6.1'))
        os.makedirs(os.path.join('dist_temp','Octave','3.6.2'))
    except os.error: pass
        
    subprocess.check_output(['OctaveBuilder.bat'],shell=True,cwd=os.path.join('wrappers','Octave'))
    shutil.copy2(os.path.join('wrappers','Octave','3.6.1','CoolProp.oct'),os.path.join('dist_temp','Octave','3.6.1','CoolProp.oct'))
    shutil.copy2(os.path.join('wrappers','Octave','3.6.2','CoolProp.oct'),os.path.join('dist_temp','Octave','3.6.2','CoolProp.oct'))
    shutil.copy2(os.path.join('wrappers','Octave','sample_code.m'),os.path.join('dist_temp','Octave','sample_code.m'))
    shutil.copy2(os.path.join('wrappers','Octave','README.rst'),os.path.join('dist_temp','Octave','README.rst'))
    
def Csharp():
    try:
        os.makedirs(os.path.join('dist_temp','C#'))
    except os.error: pass
        
    subprocess.check_output(['BuildCsharpDLL.bat'],shell=True,cwd=os.path.join('wrappers','C#'))
    shutil.copy2(os.path.join('wrappers','C#','readme.txt'),os.path.join('dist_temp','C#','readme.txt'))
    shutil.copy2(os.path.join('wrappers','C#','VSCsharp.zip'),os.path.join('dist_temp','C#','VSCsharp.zip'))
    
def MATLAB():
    try:
        os.makedirs(os.path.join('dist_temp','MATLAB'))
    except os.error: pass
        
    process = subprocess.Popen(['C:\\MATLAB_32bit\\bin\\matlab','-wait','-nodesktop','-nojvm','-r','MATLABBuilder'],shell=True,cwd=os.path.join('wrappers','MATLAB'))
    process.wait()
    process = subprocess.Popen(['matlab','-nojvm','-nodesktop','-nosplash','-wait','-r','MATLABBuilder'],shell=True,cwd=os.path.join('wrappers','MATLAB'))
    process.wait()
    shutil.copy2(os.path.join('wrappers','MATLAB','Props.mexw64'),os.path.join('dist_temp','MATLAB','Props.mexw64'))
    shutil.copy2(os.path.join('wrappers','MATLAB','HAProps.mexw64'),os.path.join('dist_temp','MATLAB','HAProps.mexw64'))
    shutil.copy2(os.path.join('wrappers','MATLAB','Props.mexw32'),os.path.join('dist_temp','MATLAB','Props.mexw32'))
    shutil.copy2(os.path.join('wrappers','MATLAB','HAProps.mexw32'),os.path.join('dist_temp','MATLAB','HAProps.mexw32'))
    shutil.copy2(os.path.join('wrappers','MATLAB','README.rst'),os.path.join('dist_temp','MATLAB','README.rst'))
    shutil.copy2(os.path.join('wrappers','MATLAB','MATLAB_sample.m'),os.path.join('dist_temp','MATLAB','MATLAB_sample.m'))
    
def Labview():
    import CoolProp
    version = CoolProp.__version__
    try:
        os.makedirs(os.path.join('dist_temp','Labview'))
    except os.error: pass
    
    shutil.copy2(os.path.join('wrappers','Labview','CoolProp.dll'),os.path.join('dist_temp','Labview','CoolProp.dll'))
    shutil.copy2(os.path.join('wrappers','Labview','CoolProp.vi'),os.path.join('dist_temp','Labview','CoolProp.vi'))
    shutil.copy2(os.path.join('wrappers','Labview','README.rst'),os.path.join('dist_temp','Labview','README.rst'))

def EES():
    import CoolProp
    version = CoolProp.__version__
    try:
        os.makedirs(os.path.join('dist_temp','EES'))
    except os.error: pass
        
    process = subprocess.Popen(['BuildLIB.bat'],shell=True,cwd=os.path.join('wrappers','EES'))
    process.wait()
    process = subprocess.Popen(['BuildDLF.bat'],shell=True,cwd=os.path.join('wrappers','EES'))
    process.wait()
    process = subprocess.Popen(['BuildMSI.bat'],shell=True,cwd=os.path.join('wrappers','EES'))
    process.wait()
    
    shutil.copy2(os.path.join('wrappers','EES','Debug','CoolProp_EES_installer.msi'),os.path.join('dist_temp','EES','CoolProp_EES_installer.msi'))
    shutil.copy2(os.path.join('wrappers','EES','CoolProp.htm'),os.path.join('dist_temp','EES','CoolProp.htm'))
    shutil.copy2(os.path.join('wrappers','EES','README.rst'),os.path.join('dist_temp','EES','README.rst'))
    
def Python():
    
    for python_install in PYTHONVERSIONS:
        print subprocess.check_output([python_install,'setup.py','bdist','--format=wininst','--dist-dir=../../dist_temp/Python'],shell=True,cwd=os.path.join('wrappers','Python'))
    
def Modelica():
    try:
        os.makedirs(os.path.join('dist_temp','Modelica'))
    except os.error: pass
        
    process = subprocess.Popen(['BuildLIB-VS2008.bat'],shell=True,cwd=os.path.join('wrappers','Modelica')); process.wait()
    process = subprocess.Popen(['BuildLIB-VS2010.bat'],shell=True,cwd=os.path.join('wrappers','Modelica')); process.wait()
        
    shutil.copy2(os.path.join('wrappers','Modelica','README.rst'),os.path.join('dist_temp','Modelica','README.rst'))
    shutil.copy2(os.path.join('wrappers','Modelica','src_modelica','CoolProp2Modelica.mo'),os.path.join('dist_temp','Modelica','CoolProp2Modelica.mo'))
    shutil.copytree(os.path.join('wrappers','Modelica','bin','VS2008'),os.path.join('dist_temp','Modelica','VS2008'))
    shutil.copytree(os.path.join('wrappers','Modelica','bin','VS2010'),os.path.join('dist_temp','Modelica','VS2010'))
    
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
    
def BuildDocs():
    #Open Doxyfile, and update the version number in the file
    lines = open('Doxyfile','r').readlines()
    import CoolProp
    for i in range(len(lines)):
        if lines[i].startswith('PROJECT_NUMBER'):
            line = lines[i].split('=')[0]+' = '+CoolProp.__version__+'\n'
            lines[i]=line
            break
    open('Doxyfile','w').write(''.join(lines))
    print subprocess.check_output(['doxygen','Doxyfile'],shell=True)
    print subprocess.check_output(['BuildCPDocs.bat'],shell=True,cwd='Web')
    
def UploadDocs():
    call_str = ['pscp','-r','-v','Web/_build/html/*.*','ibell@web.sourceforge.net:/home/groups/coolprop/htdocs']
    print 'Calling: '+' '.join(call_str)
    print subprocess.check_output(call_str,shell=True)
    
if __name__=='__main__':
    
##     InstallPrereqs()  #This is optional if you think any of the pre-reqs have been updated

##     DLL_and_Excel()
    Source()
##     Python()
##     Csharp()
##     Octave()
##     MATLAB()
##     EES()
##     Labview()
##     Modelica()
##     PYPI()
##     UploadSourceForge()

##     BuildDocs()
##     UploadDocs()