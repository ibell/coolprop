
import subprocess, os
import glob

exports = ['-s','EXPORTED_FUNCTIONS=\"[\'_main\',\'_F2K\',\'_HAProps\',\'_Props1SI\',\'_PropsSI\',\'_get_global_param_string\']\"','-s','DISABLE_EXCEPTION_CATCHING=0']
optimization = '-O2'

def compile_sources():
    for f in glob.glob(os.path.join('..','..','CoolProp','*.cpp')):
        
        # Don't compile CoolPropDLL in the compile pass to avoid duplicate symbols
        if f.find('CoolPropDLL.cpp') > -1: 
            continue 
        
        call = [r'em++',optimization,f,'-I../../CoolProp','-c','-DEXTERNC']+ exports
        print 'Calling:',' '.join(call)
        subprocess.check_output(' '.join(call), shell = True)

def link():
    call = [r'em++',optimization,'-o','coolprop.js','../../CoolProp/CoolPropDLL.cpp']+glob.glob('*.o')+['-I../../CoolProp','-DEXTERNC']  +  exports
    print 'Calling:',' '.join(call)
    subprocess.check_output(' '.join(call), shell = True)

def closure_compiler():
    call = ['java','-Xmx1024m','-jar','compiler.jar','--js','coolprop.js','--js_output_file','coolprop2.js','--compilation_level','ADVANCED_OPTIMIZATIONS','--language_in','ECMASCRIPT5']
    print 'Using the closure compiler, this will take a while...   (from https://developers.google.com/closure/compiler/)'
    print 'Calling:',' '.join(call)
    subprocess.check_output(' '.join(call), shell = True)

def cleanup():
    for file in glob.glob('*.o'):
        print 'removing',file
        os.remove(file)

def run():
    os.startfile('index.html')

if __name__=='__main__':
    compile_sources()
    link()
    cleanup()
#     run()
#    closure_compiler()