from CoolProp.Plots.Plots import hs
import CoolProp
from CoolProp.CoolProp import Props
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import json

from build_DTU_JSON import RP2CAS

# CAS ; Name ; Tmin [K] ; Tmax [K] ; pmax [Pa]
limits_data = """7732-18-5;Water;273.16;1273;1000000000
811-97-2;R134a;169.85;455;70000000
7440-59-7;Helium;2.1768;2000;1000000000
7782-44-7;Oxygen;54.361;300;80000000
1333-74-0;Hydrogen;13.957;1000;2000000000
1333-74-0p;ParaHydrogen;13.8033;1000;2000000000
1333-74-0o;OrthoHydrogen;14.008;1000;2000000000
7440-37-1;Argon;83.806;700;1000000000
124-38-9;CarbonDioxide;216.592;1100;800000000
7727-37-9;Nitrogen;63.151;1100;2200000000
74-98-6;n-Propane;85.525;650;1000000000
7664-41-7;Ammonia;195.495;700;1000000000
754-12-1;R1234yf;220;410;30000000
29118-24-9;R1234ze(E);240;410;15000000
75-10-5;R32;136.34;435;70000000
75-45-6;R22;115.73;550;60000000
SES36.ppf;SES36;278.674;725;500000000
74-85-1;Ethylene;103.989;450;300000000
2551-62-4;SulfurHexafluoride;223.555;650;150000000
64-17-5;Ethanol;159.1;650;280000000
115-10-6;DimethylEther;131.66;550;50000000
616-38-6;DimethylCarbonate;277.06;600;60000000
420-46-2;R143a;161.34;450;150000000
75-46-7;R23;118.02;425;120000000
112-40-3;n-Dodecane;263.6;800;200000000
115-07-1;Propylene;87.953;575;1000000000
287-92-3;Cyclopentane;179.7;550;250000000
690-39-1;R236FA;179.6;400;70000000
431-63-0;R236EA;240;412;6000000
431-89-0;R227EA;146.35;475;60000000
406-58-6;R365MFC;239;500;35000000
353-36-6;R161;130;450;5000000
421-14-7;HFE143m;240;4220;7200000
71-43-2;Benzene;278.674;725;500000000
1120-21-4;n-Undecane;247.541;700;500000000
354-33-6;R125;172.52;500;60000000
75-19-4;CycloPropane;273;473;28000000
7440-01-9;Neon;24.56;723;700000000
2837-89-0;R124;100;470;40000000
74-99-7;Propyne;323;474;31800000
7782-41-4;Fluorine;53.4811;300;20000000
67-56-1;Methanol;175.61;580;500000000
115-25-3;RC318;233.35;623;60000000
75-43-4;R21;200;473;138000000
76-14-2;R114;273.15;507;21000000
75-72-9;R13;98.15;450;50000000
75-73-0;R14;120;623;51000000
75-71-8;R12;116.099;525;200000000
76-13-1;R113;236.93;525;200000000
29118-25-0;R1234ze(Z);273;430;6000000
64-19-7;AceticAcid;289.8;500;30000000
460-73-1;R245fa;171.05;440;200000000
593-53-3;R41;129.82;425;70000000
630-08-0;CarbonMonoxide;68.16;500;100000000
463-58-1;CarbonylSulfide;134.3;650;50000000
124-18-5;n-Decane;243.5;675;800000000
7783-06-4;HydrogenSulfide;187.7;760;170000000
78-78-4;Isopentane;112.65;500;1000000000
463-82-1;Neopentane;256.6;550;200000000
107-83-5;Isohexane;119.6;500;1000000000
7439-90-9;Krypton;115.77;750;200000000
111-84-2;n-Nonane;219.7;600;800000000
108-88-3;Toluene;178;700;500000000
7440-63-3;Xenon;161.4;750;700000000
76-16-4;R116;173.1;425;50000000
67-64-1;Acetone;178.5;550;700000000
10024-97-2;NitrousOxide;182.33;525;50000000
7446-09-5;SulfurDioxide;197.7;525;35000000
1717-00-6;R141b;169.68;500;400000000
75-68-3;R142b;142.72;470;60000000
76-19-7;R218;125.45;440;20000000
74-82-8;Methane;90.6941;625;1000000000
74-84-0;Ethane;90.368;675;900000000
106-97-8;n-Butane;134.895;575;12000000
75-28-5;IsoButane;113.73;575;35000000
109-66-0;n-Pentane;143.47;573;69000000
110-54-3;n-Hexane;177.83;548;92000000
142-82-5;n-Heptane;182.55;600;100000000
111-65-9;n-Octane;216.37;548;96000000
75-37-6;R152A;154.56;471;58000000
306-83-2;R123;166;523;76000000
75-69-4;R11;162.68;595;100000000
107-51-7;MDM;187.2;673;30000000
141-62-8;MD2M;205.2;673;30000000
141-63-9;MD3M;192;673;30000000
540-97-6;D6;270.2;673;30000000
107-46-0;MM;273;673;30000000
107-52-8;MD4M;300;673;30000000
556-67-2;D4;290.25;673;30000000
541-02-6;D5;273;673;30000000
106-98-9;1-Butene;87.8;525;50000000
115-11-7;IsoButene;132.4;525;50000000
590-18-1;cis-2-Butene;134.3;525;50000000
624-64-6;trans-2-Butene;167.6;525;50000000
112-39-0;MethylPalmitate;302.71;700;50000000
112-61-8;MethylStearate;311.84;700;50000000
112-62-9;MethylOleate;253.47;700;50000000
112-63-0;MethylLinoleate;238.1;700;50000000
301-00-8;MethylLinolenate;218.65;700;50000000
95-47-6;o-Xylene;247.985;700;70000000
108-38-3;m-Xylene;225.3;700;200000000
106-42-3;p-Xylene;286.4;700;200000000
100-41-4;EthylBenzene;178.2;700;60000000
7782-39-0;Deuterium;18.724;600;2000000000
7782-39-0p;ParaDeuterium;18.724;600;2000000000
7782-39-0o;OrthoDeuterium;18.724;600;2000000000
AIR.PPF;Air;60;2000;2000000000
R404A.PPF;R404A;200;450;50000000
R410A.PPF;R410A;200;450;50000000
R407C.PPF;R407C;200;450;50000000
R507A.PPF;R507A;200;450;50000000
102687-65-0;R1233zdE;195.15;550;100000000
110-82-7;Cyclohexane;279.47;700;250000000
R407F.ppf;R407F;200;450;50000000"""

limits = {}

for line in limits_data.split('\n'):
    CAS,Name,Tmin,Tmax,pmax = line.split(';')
    el = dict(Name = Name, Tmin = float(Tmin), Tmax = float(Tmax), pmax = float(pmax))
    limits[CAS] = el

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
                              'AceticAcid':'64-19-7'
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
    
    code[CAS]['Tmin'] = limits[CAS]['Tmin']
    code[CAS]['Tmax'] = limits[CAS]['Tmax']
    code[CAS]['pmax'] = limits[CAS]['pmax']
    
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