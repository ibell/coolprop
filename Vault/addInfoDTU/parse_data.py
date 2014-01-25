
fluids = """:'1BUTENE.FLD','ACETONE.FLD','AIR.PPF','AMMONIA.FLD','ARGON.FLD',
:'BENZENE.FLD','BUTANE.FLD','C1CC6.FLD','C2BUTENE.FLD','C3CC6.FLD',
:'C4F10.FLD','C5F12.FLD','C12.FLD','CF3I.FLD','CO.FLD','CO2.FLD',
:'COS.FLD','CYCLOHEX.FLD','CYCLOPEN.FLD','CYCLOPRO.FLD','D2.FLD',
:'D2O.FLD','D4.FLD','D5.FLD','D6.FLD','DECANE.FLD','DMC.FLD',
:'DME.FLD','ETHANE.FLD','ETHANOL.FLD','ETHYLENE.FLD','FLUORINE.FLD'
:,'H2S.FLD','HELIUM.FLD','HEPTANE.FLD','HEXANE.FLD','HYDROGEN.FLD',
:'IBUTENE.FLD','IHEXANE.FLD','IPENTANE.FLD','ISOBUTAN.FLD',
:'KRYPTON.FLD','MD2M.FLD','MD3M.FLD','MDM.FLD','METHANE.FLD',
:'METHANOL.FLD','MLINOLEA.FLD','MLINOLEN.FLD','MM.FLD',
:'MOLEATE.FLD','MPALMITA.FLD','MSTEARAT.FLD','N2O.FLD','NEON.FLD',
:'NEOPENTN.FLD','NF3.FLD','NITROGEN.FLD','NONANE.FLD',
:'OCTANE.FLD','ORTHOHYD.FLD','OXYGEN.FLD','PARAHYD.FLD',
:'PENTANE.FLD','PROPANE.FLD','PROPYLEN.FLD','PROPYNE.FLD',
:'R32.FLD','R41.FLD','R115.FLD','R116.FLD','R124.FLD','R125.FLD',
:'R141B.FLD','R142B.FLD','R143A.FLD','R161.FLD','R218.FLD',
:'R227EA.FLD','R236EA.FLD','R236FA.FLD','R245CA.FLD','R245FA.FLD',
:'R365MFC.FLD','R507A.PPF','R1234YF.FLD','R1234ZE.FLD','SF6.FLD',
:'SO2.FLD','T2BUTENE.FLD','TOLUENE.FLD','WATER.FLD','XENON.FLD',
:'R11.FLD','R12.FLD','R13.FLD','R14.FLD','R21.FLD','R22.FLD',
:'R23.FLD','R113.FLD','R114.FLD','R123.FLD','R134A.FLD','R152A.FLD',
:'R404A.PPF','R407C.PPF','R410A.PPF','RC318.FLD'"""

HH = """:'0','2','0','3','0',
:'2','1','2','0','NA',
:'1','NA','2','1','1','1',
:'3','1','2','2','NA',
:'NA','NA','NA','NA','2','2',
:'1','1','2','2','4'
:,'4','0','1','2','0',
:'1','2','1','1',
:'0','1','1','1','0',
:'2','NA','NA','2',
:'2','1','0','1','0',
:'1','1','0','2',
:'2','NA','0','NA',
:'2','1','1','1',
:'1','2','1','1','1','1',
:'1','1','1','NA','2',
:'1','NA','1','NA','2',
:'NA','1','1','1','1',
:'3','0','2','0','0',  
:'1','1','NA','NA','NA','1',
:'1','1','1','2','1','1',
:'1','1','1','1'"""

FH = """:'4','3','0','1','0',
:'3','4','3','4','NA',
:'0','NA','2','0','4','0',
:'4','3','3','2','NA',
:'NA','NA','NA','NA','2','3',
:'4','4','3','4','3'
:,'4','0','3','3','4',
:'4','3','4','4',
:'0','2','2','2','4',
:'3','NA','NA','3',
:'1','0','1','0','0',
:'4','0','0','3',
:'3','NA','0','NA',
:'4','4','4','4',
:'4','3','0','0','1','1',
:'1','1','1','NA','1',
:'0','NA','0','NA','0',
:'NA','1','2','2','0',
:'0','4','3','0','0',  
:'1','1','NA','NA','NA','1',
:'1','0','0','1','1','4',
:'1','1','1','0' """

PH = """:'0','0','0','0','0',
:'0','0','0','1','NA',
:'0','NA','0','0','3','0',
:'1','0','1','0','NA',
:'NA','NA','NA','NA','0','0',
:'2','0','0','2','0'
:,'0','0','0','0','0',
:'2','1','0','0',
:'0','0','0','0','0',
:'0','NA','NA','1',
:'0','0','0','0','0',
:'0','3','0','0',
:'0','NA','0','NA',
:'0','0','1','1',
:'1','2','1','1','0','0',
:'0','0','0','NA','1',
:'0','NA','1','NA','1',
:'NA','0','0','0','0',
:'0','1','0','0','3',  
:'0','0','NA','NA','NA','0',
:'0','0','0','0','0','1',
:'0','0','0','2'"""

pp_fluids = fluids.replace(':','').replace('\n','').replace('.FLD','').replace('.PPF','').replace("'",'').split(",")
pp_HH = HH.replace(':','').replace('\n','').replace("'",'').split(",")
pp_FH = FH.replace(':','').replace('\n','').replace("'",'').split(",")
pp_PH = PH.replace(':','').replace('\n','').replace("'",'').split(",")

from fluid_lookup import *

HH_dict = {RP2CAS[k]:v for k,v in zip(pp_fluids,pp_HH)}
FH_dict = {RP2CAS[k]:v for k,v in zip(pp_fluids,pp_FH)}
PH_dict = {RP2CAS[k]:v for k,v in zip(pp_fluids,pp_PH)}

def get_env_data(fluid):
    a = dict(
        HH = HH_dict[fluid],
        FH = FH_dict[fluid],
        PH = PH_dict[fluid],
        ODP = ODP_dict[fluid],
        GWP20 = GWP20_dict[fluid],
        GWP100 = GWP100_dict[fluid],
        GWP500 = GWP500_dict[fluid]
        )
        
    for k,v in a.iteritems():
        try:
            a[k] = int(v)
        except ValueError:
            try:
                a[k] = float(v)
            except ValueError: 
                a[k] = None
                
    for term in ['GWP100','GWP20','GWP500','ODP']:
        try:
            a[term] = float(a[term])
        except TypeError:
            a[term] = -1
        
    if fluid in ASHRAE34_dict:
        a['ASHRAE34'] = ASHRAE34_dict[fluid]
    else:
        a['ASHRAE34'] = 'UNKNOWN'
    
    return a
        
import json

code = {}
for fluid in pp_fluids:
    fluid = RP2CAS[fluid]
    if name_dict[fluid]:
        code[fluid] = get_env_data(fluid)
    else:
        continue

if __name__=='__main__':
    print "DO NOT USE THIS SCRIPT! Use the new one in dev/build_JSON_code.py  instead"
    f = open('environmental.json','w')
    print>>f, json.dumps(code, sort_keys=True, indent=2, separators=(',', ': '))
    f.close()