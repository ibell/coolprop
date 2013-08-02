"""
A script for making nightly builds.

Uses whatever the default python distro is on the building computer

To Use:
 - SSH key should be uploaded to SF servers and key managed by pageant on building computer
 - Add a shortcut to the .ppk private key file to the Startup folder in windows so that pageant will be started and key loaded at logon - you will be asked for password when you log in
 - Windows Task should be setup in the task scheduler that will run python.exe with argument nightly_build.py in the source folder
 
Change ibell to your sourceforge name below
"""

def make_temp_folder():
    while True:
        try:
            os.mkdir('dist_tmp')
            break
        except WindowsError:
            pass
            
def delete_temp_folder():
    if os.path.exists('dist_tmp'):
        shutil.rmtree('dist_tmp')
            
import subprocess, os, shutil, glob

delete_temp_folder()
make_temp_folder()

call_str = ['python','setup.py','bdist_msi','--dist-dir=../../dist_tmp']
print 'Calling: '+' '.join(call_str)
print subprocess.check_output(call_str,cwd=os.path.join('wrappers','Python'))

fpath = glob.glob('dist_tmp/*.msi')[0]
fName = fpath.split(os.sep,1)[1]
call_str = ['pscp',fpath,'ibell,coolprop@frs.sf.net:/home/pfs/project/c/co/coolprop/CoolProp/Nightly/'+fName]
print 'Calling: '+' '.join(call_str)
print subprocess.check_output(call_str)

delete_temp_folder()