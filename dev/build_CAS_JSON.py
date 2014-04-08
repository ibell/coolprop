import CoolProp, json
from build_DTU_JSON import RP2CAS

code = {}
for fluid in CoolProp.__fluids__:
    RPName = CoolProp.CoolProp.get_REFPROPname(fluid)
    try:
        CAS = RP2CAS[RPName.upper()]
    except KeyError:
        NOT_IN_REFPROP_CAS = {'R1234ze(Z)':'29118-25-0',
                              'ParaDeuterium':'7782-39-0p',
                              'OrthoDeuterium':'7782-39-0o',
                              'R407F': 'R407F.ppf',
                              'AceticAcid':'64-19-7'
                              }
        CAS = NOT_IN_REFPROP_CAS[fluid]
    
    code[fluid] = CAS
    
####################### WRITE THE FILE #################################
####################### WRITE THE FILE #################################
####################### WRITE THE FILE #################################

# Dump all to a string
s = json.dumps(code, sort_keys=True, indent=2, separators=(',', ': '))

f = open('JSON_CAS.h','w')

# Header
print >> f, 'const char JSON_cas[] = "' + s.split('\n')[0].replace('"','\\"') + '"'

# Modify the string to replace " with \"
for line in s[1::].split('\n'):
    print >> f, '"' + line.replace('"','\\"') + '"'

print >> f, ";\n"