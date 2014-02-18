from CoolProp.Plots.Plots import hs
import CoolProp
from CoolProp.CoolProp import Props
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import json

from build_DTU_JSON import RP2CAS

env_json = json.loads(open('DTU_environmental.json','r').read())

def get_environmental_data(fluid):
    if fluid in env_json:
        return env_json[fluid]
    else:
        print fluid,'not found in env_json, filling with empty values'
        return dict(GWP100 = -1, 
                    GWP20 = -1, 
                    GWP500 = -1, 
                    ODP = -1, 
                    HH = -1,
                    FH = -1,
                    PH = -1,
                    ASHRAE34 = "UNKNOWN"
                    )
    
def hsatVmax(fluid):
    Tmin = Props(fluid, 'Tmin')
    Tmax = Props(fluid, 'Tcrit')
    
    def OBJECTIVE(T):
        return -Props('H','T',T,'Q',1,fluid)
    
    T = scipy.optimize.minimize_scalar(OBJECTIVE,bounds = (Tmin,Tmax),method = 'Bounded').x
    
    h = Props('H','T',T,'Q',1,fluid)
    s = Props('S','T',T,'Q',1,fluid)
    rho = Props('D','T',T,'Q',1,fluid)
    
    return h*1000,T,s*1000,rho
    
def fit_hs(fluid):
    T = np.linspace(Props(fluid, 'Tmin'), Props(fluid, 'Tcrit')-2)
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
    RPName = CoolProp.CoolProp.get_REFPROPname(fluid)
    try:
        CAS = RP2CAS[RPName.upper()]
    except KeyError:
        NOT_IN_REFPROP_CAS = {'R1234ze(Z)':'29118-25-0',
                              'ParaDeuterium':'7782-39-0p',
                              'OrthoDeuterium':'7782-39-0o',
                              'R407F':'R407F.ppf',
                              }
        CAS = NOT_IN_REFPROP_CAS[fluid]
    
    code[CAS] = {}
    if CAS.upper().endswith('.PPF'):
        code[CAS]['CAS'] = 'N/A'
    else:
        code[CAS]['CAS'] = CAS
    print fluid, RPName, code[CAS]['CAS']
        
    code[CAS]['hsatVmax'],code[CAS]['T_hsatVmax'],code[CAS]['s_hsatVmax'],code[CAS]['rho_hsatVmax'] = hsatVmax(fluid)
    code[CAS].update(get_environmental_data(CAS))
    
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