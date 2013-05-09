header = """// This file contains the environmental data provided by the Danish Techological University (DTU)
//  
//

void syaml_build()
{
    std::map<std::string, double> fluid_map;
    
"""

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

HH_dict = {k:v for k,v in zip(pp_fluids,pp_HH)}
FH_dict = {k:v for k,v in zip(pp_fluids,pp_FH)}
PH_dict = {k:v for k,v in zip(pp_fluids,pp_PH)}
    
f = open('../../CoolProp/EnvironmentalData.h','w')

f.write(header)

template = """    fluid_map.clear();
    fluid_map["HH"] = {HH:s};
    fluid_map["FH"] = {FH:s};
    fluid_map["PH"] = {PH:s};
    fluid_map["GWP20"] = {GWP20:s};
    fluid_map["GWP100"] = {GWP100:s};
    fluid_map["GWP500"] = {GWP500:s};
    fluid_map["ODP"] = {ODP:s};
    ASHRAE34_map()["{fluid:s}"] = "{ASHRAE34:s}";
    syaml_environmental_map()["{fluid:s}"] = fluid_map;

"""

from fluid_lookup import *

for fluid in pp_fluids:
    if name_dict[fluid]:
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
                a[k] = str(float(v))
            except ValueError:
                a[k] = '_HUGE'
        
        if a['ODP'] == '_HUGE': a['ODP'] = '0'
            
        if fluid in ASHRAE34_dict:
            a['ASHRAE34'] = ASHRAE34_dict[fluid]
        else:
            a['ASHRAE34'] = 'UNKNOWN'
        
        chunk = template.format(fluid = name_dict[fluid],
                                HH = a['HH'],
                                FH = a['FH'],
                                PH = a['PH'],
                                ODP = a['ODP'],
                                GWP20 = a['GWP20'],
                                GWP100 = a['GWP100'],
                                GWP500 = a['GWP500'],
                                ASHRAE34 = a['ASHRAE34'],
                                )
        f.write(chunk)
                            
f.write('}')
f.close()