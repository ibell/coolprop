from CoolProp.Plots.Plots import hs
import CoolProp
from CoolProp.CoolProp import Props
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import json

from CAS_data_generator import get_CAS

env_json = json.loads(open('environmental.json','r').read())

def get_environmental_data(fluid):
    if fluid in env_json:
        return env_json[fluid]
    else:
        print fluid,'not found in env_json'
        return {}
    
def hsatVmax(fluid):
    Tmin = Props(fluid, 'Tmin')
    Tmax = Props(fluid, 'Tcrit')
    
    def OBJECTIVE(T):
        return -Props('H','T',T,'Q',1,fluid)
    
    T = scipy.optimize.minimize_scalar(OBJECTIVE,bounds = (Tmin,Tmax),method = 'Bounded').x
    
    h = Props('H','T',T,'Q',1,fluid)
    s = Props('S','T',T,'Q',1,fluid)
    rho = Props('D','T',T,'Q',1,fluid)
    
    return h,T,s,rho
    
def fit_hs(fluid):
    T = np.linspace(Props(fluid, 'Tmin'), Props(fluid, 'Tcrit'))
    sL = Props('S', 'T', T, 'Q', 0, fluid)
    hL = Props('H', 'T', T, 'Q', 0, fluid)
    a = np.polyfit(sL,hL,4)
    n = range(4,-1,-1)
    
    d = dict(a_hs_satL = list(a),
             n_hs_satL = list(n)
             )
    return d
    
################## GENERATE THE JSON PRECURSOR DICTIONARY ######################
################## GENERATE THE JSON PRECURSOR DICTIONARY ######################
################## GENERATE THE JSON PRECURSOR DICTIONARY ######################

code = {}
for fluid in CoolProp.__fluids__:
    code[fluid] = {}
    code[fluid]['CAS'] = get_CAS(fluid)
    code[fluid]['hsatVmax'],code[fluid]['T_hsatVmax'],code[fluid]['s_hsatVmax'],code[fluid]['rho_hsatVmax'] = hsatVmax(fluid)
    code[fluid].update(get_environmental_data(fluid))
    if fluid in ['ParaHydrogen','MethylLinoleate','MethylLinolenate','R407C','R507A']:    
        continue
    code[fluid].update(fit_hs(fluid))
    
    
####################### WRITE THE FILE #################################
####################### WRITE THE FILE #################################
####################### WRITE THE FILE #################################

# Dump all to a string
s = json.dumps(code, sort_keys=True, indent=2, separators=(',', ': '))

f = open('JSON_code.h','w')

# Header
print >> f, 'const char JSON_code[] = "' + s.split('\n')[0].replace('"','\\"') + '"'

# Modify the string to replace " with \"
for line in s[1::].split('\n'):
    print >> f, '"' + line.replace('"','\\"') + '"'

print >> f, ";\n"
f.close()