import matplotlib.pyplot as plt
import numpy as np
from CoolProp.CoolProp import Props
import textwrap
'''
Functional forms are based on the ancillary equations from R134a
'''

# Key is the REFPROP name(here case sensitive), then tuple of CoolProp name followed by aliases
RefLookup={'acetone':('Acetone',),
            'CO':('CarbonMonoxide','CO'),
            'COS':('CarbonylSulfide','COS'),
            'decane':('Decane',),
            'H2S':('HydrogenSulfide','H2S'),
            'ipentane':('Isopentane',),
            'neopentn':('Neopentane',),
            'ihexane':('Isohexane',),
            'krypton':('Krypton',),
            'N2O':('NitrousOxide','N2O'),
            'nonane':('Nonane',),
            'SO2':('SulfurDioxide','SO2'),
            'toluene':('Toluene',),
            'xenon':('Xenon','Xe'),
            'R116':('R116',),
            'R141b':('R141b',),
            'R142b':('R142b',),
            'R218':('R218',),
            'R245fa':('R245fa',),
            'R41':('R41',),
            }
Tt={}
Tbp={}
MM={}
omega={}
Tmax={}
pmax={}
Tc={}
pc={}
rhoc={}
c0={}
c1={}
c2={}
a1={}
a2={}
n={}
v={}
u={}

def fluidconstants():
    lines=open('crit.txt','r').readlines()
    for line in lines:
        if line.startswith('#'):
            continue
        vals=line.split(',')
        Tc[vals[0]]=float(vals[1])
        pc[vals[0]]=float(vals[2])
        rhoc[vals[0]]=float(vals[3])
        
    lines=open('constants.txt','r').readlines()
    for line in lines:
        if line.startswith('#'):
            continue
        vals=line.split(',')
        Tt[vals[0]]=float(vals[1])
        Tbp[vals[0]]=float(vals[2])
        MM[vals[0]]=float(vals[3])
        omega[vals[0]]=float(vals[4])
        Tmax[vals[0]]=float(vals[5])
        pmax[vals[0]]=float(vals[6])
        
    lines=open('idealgas.txt','r').readlines()
    for line in lines:
        if line.startswith('#'):
            continue
        vals=line.split(',')
        if len(vals)==4:
            c0[vals[0]]=float(vals[1])
            c1[vals[0]]=0.0
            c2[vals[0]]=0.0
            a1[vals[0]]=float(vals[2])
            a2[vals[0]]=float(vals[3])
        elif len(vals)==6:
            c0[vals[0]]=float(vals[1])
            c1[vals[0]]=float(vals[2])
            c2[vals[0]]=float(vals[3])
            a1[vals[0]]=float(vals[4])
            a2[vals[0]]=float(vals[5])
        else:
            print 'error',len(vals)
            
    lines = open('nonpolar_transposed.txt','r').readlines()
    for line in lines:
        if line.startswith('#'):
            continue
        vals=line.split('\t')
        n[vals[0]]=vals[1::]
        
    lines = open('VU.txt','r').readlines()
    for line in lines:
        if line.startswith('#'):
            continue
        vals=line.rstrip().split(',')
        v[vals[0]]=vals[1::2]
        u[vals[0]]=vals[2::2]
        
from scipy.optimize import curve_fit
f1=plt.figure()
f2=plt.figure()
f3=plt.figure()

def Ancillary(Ref,plot=False):
    print Ref
    Tcrit=Props(Fluid,'Tcrit')
    Ttrip=Props(Fluid,'Ttriple')
    pcrit = Props(Fluid,'pcrit')
    rhoc = Props('D','T',Tcrit,'P',pcrit,Fluid)
    Tv=np.linspace(Ttrip,0.995*Tcrit,500)
    rhoL = np.zeros_like(Tv)
    rhoV = np.zeros_like(Tv)
    p = np.zeros_like(Tv)
    for i,T in enumerate(Tv):
        rhoL[i] = Props('D','T',T,'Q',0,Fluid)
        rhoV[i] = Props('D','T',T,'Q',1,Fluid)
        p[i] = Props('P','T',T,'Q',1,Fluid)
    
    def fun_rhoL(x,a,b,c,d,e):
        return a+b*np.power(x,c)+d*np.power(x,e)
        
    popt,pcov=curve_fit(fun_rhoL,1-Tv/Tcrit,rhoL,[500,800,0.33,500,0.66])
    
    Refbare = RefLookup[Fluid.split('-')[1]][0]
    rhoL_code=('double {0:s}Class::rhosatL(double T) \n'.format(Refbare)
    +'{\n    double THETA = 1-T/crit.T;\n'
    +'    return {0:g}{1:+g}*pow(THETA,{2:g}){3:+g}*pow(THETA,{4:g});\n'.format(*popt)
    +'}')
    rhoL_fit = fun_rhoL(1-Tv/Tcrit,*popt)
    print 'Fluid {0:s} max error {1:g} % in rhoL'.format(Fluid.split('-')[1],np.max(np.abs(rhoL_fit/rhoL*100-100)))
    if plot==True:
        f1.gca().plot(1-Tv/Tcrit,rhoL,'-',1-Tv/Tcrit,rhoL_fit,'--')
    
    def fun_rhoV(x,a,b,c,d,e,f,g):
        return a+b*np.power(x,c)+d*np.power(x,e)+f*np.power(x,g)
        
    popt,pcov=curve_fit(fun_rhoV,1-Tv/Tcrit,np.log(rhoV/rhoc),[1,0,0.33,0,0.66,1,2.0])
    rhoV_code=('double {0:s}Class::rhosatV(double T) \n'.format(Refbare)
    +'{\n    double THETA = 1-T/crit.T;\n'
    +'    double RHS = {0:g}{1:+g}*pow(THETA,{2:g}){3:+g}*pow(THETA,{4:g}){5:+g}*pow(THETA,{6:g});\n'.format(*popt)
    +'    return exp(RHS)*crit.rho;'
    +'\n}')
    
    rhoV_fit = np.exp(fun_rhoV(1-Tv/Tcrit,*popt))*rhoc
    print 'Fluid {0:s} max error {1:g} % in rhoV'.format(Fluid.split('-')[1],np.max(np.abs(rhoV_fit/rhoV*100-100)))
    if plot==True:
        f2.gca().plot(1-Tv/Tcrit,np.log(rhoV/rhoc),1-Tv/Tcrit,fun_rhoV(1-Tv/Tcrit,*popt))
    
    def fun_psat(x,a,b,c,d,e):
        return a+b*np.power(x,c)+d*np.power(x,e)
        
    popt,pcov=curve_fit(fun_psat,1-Tv/Tcrit,Tv/Tcrit*np.log(p/pcrit),[1,1,0.33,1,0.66])
    psat_code=('double {0:s}Class::psat(double T) \n'.format(Refbare)
    +'{\n    double THETA = 1-T/crit.T;\n'
    +'    double RHS = {0:g}{1:+g}*pow(THETA,{2:g}){3:+g}*pow(THETA,{4:g});\n'.format(*popt)
    +'    return exp(crit.T/T*RHS)*crit.p;'
    +'\n}')
    p_fit = pcrit*np.exp(fun_psat(1-Tv/Tcrit,*popt)*Tcrit/Tv)
    print 'Fluid {0:s} max error {1:g} % in psat'.format(Fluid.split('-')[1],np.max(np.abs(p_fit/p*100-100)))
    if plot==True:
        f3.gca().plot(1-Tv/Tcrit,-Tv/Tcrit*np.log(p/pcrit))
        f3.gca().set_xlabel('1-T/Tc [-]')
        f3.gca().set_ylabel('-T/Tc*log(p/pc) [-]')
    
    return rhoL_code, rhoV_code, psat_code
    
def buildEOS(Ref,CPPName,polar=True):
    #Ref is the REFPROP name
    Fluid = Ref.split('-')[1]; ## The CoolProp name

    if polar==True:
        a='a_polar'
        d='d_polar'
        t='t_polar'
        c='c_polar'
        
    else:
        a='a_nonpolar'
        d='d_nonpolar'
        t='t_nonpolar'
        c='c_nonpolar'
    
    key = RefLookup[Fluid][0]
    
    nstring='const double n[]={0.0,'+','.join(map(str,n[key])).rstrip()+'};'
    v0string='const double v0[]={0.0,'+','.join(map(str,v[key])).rstrip()+'};'
    u0string='const double u0[]={0.0,'+','.join(map(str,u[key])).rstrip()+'};'
    
    if c1[key]==0:
        c1string='0.0'
    else:
        c1string='-{c1}*pow(crit.T,{c2})/({c2}*({c2}+1))'.format(c1=c1[key],c2=c2[key])
    
    print str(c1[key]),str(c1[key])=='0.0'
    if str(c1[key])=='0.0':
        phi0_power_def = ''
    else:
        phi0_power_def = 'phi_BC * phi0_power_ = new phi0_power({c1string},-{c2});\n\tphi0list.push_back(phi0_power_);'.format(c1string=c1string,c2=c2[key]);
        
    if abs(float(u[key][0]))<1e-12:
        phi0_PE_def = ''
    else:
        phi0_PE_def = 'phi_BC * phi0_Planck_Einstein_ = new phi0_Planck_Einstein(v0_v,u0_v,1,v0_v.size()-1);\n\tphi0list.push_back(phi0_Planck_Einstein_);';
        
    rhoL_code,rhoV_code,psat_code = Ancillary(Ref)
    CPPString=textwrap.dedent(
"""

{Fluid:s}Class::{Fluid:s}Class()
{{
    {nstring}
    {u0string}
    {v0string}
    
    // Critical parameters
    crit.rho = {rhoc};
    crit.p = {pc};
    crit.T = {Tc};
    crit.v = 1.0/crit.rho;
    
    // Other fluid parameters
    params.molemass = {MM};
    params.Ttriple = {Tt};
    params.accentricfactor = {omega};
    params.R_u = 8.314472;

    // Limits of EOS
    limits.Tmin = {Tmin};
    limits.Tmax = {Tmax};
    limits.pmax = {pmax};
    limits.rhomax = {rhomax_molar}*params.molemass;    
    
    std::vector<double> n_v(n,n+sizeof(n)/sizeof(double));
    std::vector<double> d_v({d},{d}+sizeof({d})/sizeof(double));
    std::vector<double> t_v({t},{t}+sizeof({t})/sizeof(double));
    std::vector<double> l_v({c},{c}+sizeof({c})/sizeof(double));
    std::vector<double> u0_v(u0,u0+sizeof(u0)/sizeof(double));
    std::vector<double> v0_v(v0,v0+sizeof(v0)/sizeof(double));
    
    for (unsigned int i=0;i<u0_v.size();i++) {{ u0_v[i]/=crit.T; }}

    phi_BC * phir_ = new phir_power(n_v,d_v,t_v,l_v,1,12);
    phirlist.push_back(phir_);

    phi_BC * phi0_lead_ = new phi0_lead({a1},{a2});
    phi0list.push_back(phi0_lead_);
    
    phi_BC * phi0_logtau_ = new phi0_logtau({c0}-1);
    phi0list.push_back(phi0_logtau_);
    
    {phi0_power_def}
    {phi0_PE_def}
    
    EOSReference.assign("Lemmon, E.W., and R. Span, \\"Short Fundamental Equations of State for 20 Industrial Fluids,\\", J. Chem. Eng. Data, 51:785-850, 2006.");
    TransportReference.assign(\"Using ECS\");

    name.assign(\"{Fluid:s}\");
    {aliases}
    
}}
double {Fluid}Class::viscosity_Trho(double tau, double delta)
{{
    //throw NotImplementedError();
    return _HUGE;
}}
double {Fluid}Class::conductivity_Trho(double tau, double delta)
{{
    //throw NotImplementedError();
    return _HUGE;
}}
{rhoL_code}
{rhoV_code}
{psat_code}

""".format(Fluid=key,d=d,t=t,c=c,rhoc=rhoc[key]*MM[key],pc=pc[key]*1000 ,Tc=Tc[key],MM=MM[key],Tt=Tt[key],omega=omega[key],Tmin=Tt[key],Tmax=Tmax[key],pmax=pmax[key]*1000,rhomax_molar=1e6,u0string=u0string, v0string=v0string,a1=a1[key],a2=a2[key],c0=c0[key],c1string=c1string,c2=c2[key],aliases='//Aliases not implemented yet',nstring=nstring,rhoL_code=rhoL_code,rhoV_code=rhoV_code,psat_code=psat_code,phi0_power_def=phi0_power_def,phi0_PE_def=phi0_PE_def))
    fp = open(CPPName,'a')
    fp.write(CPPString)
    fp.close()
    
def AddHeader(Ref,HeadName):
    #Ref is the REFPROP name
    Fluid = Ref.split('-')[1]; ## The CoolProp name
    
    print RefLookup[Fluid][0]
    HeadString=textwrap.dedent("""
    class {Fluid}Class : public Fluid {{

    public:
        {Fluid}Class();
        ~{Fluid}Class(){{}};
        virtual double conductivity_Trho(double, double);
        virtual double viscosity_Trho(double, double);
        double psat(double);
        double rhosatL(double);
        double rhosatV(double);
    }};
    """.format(Fluid=RefLookup[Fluid][0])
    )
    
    fp = open(HeadName,'a')
    fp.write(HeadString)
    fp.close()
    
CPPheader = textwrap.dedent("""

// **** WARNING ******
// **** WARNING ******
// **** WARNING ******

// Do NOT modify this file.  It is created by a script in the industrialfluidsbuilder folder within the source

#include \"CoolProp.h\"
#include "IndustrialFluids.h"
#include <vector>
#include "CPExceptions.h"

static const double d_nonpolar[] =
{
0,
1.0, //[1]
1.0, //[2]
1.0, //[3]
2.0, //[4]
3.0, //[5]
7.0, //[6]
2.0, //[7]
5.0, //[8]
1.0, //[9]
4.0, //[10]
3.0, //[11]
4.0, //[12]
};

static const double t_nonpolar[] =
{
0,
0.25,  //[1]
1.125, //[2]
1.5,   //[3]
1.375, //[4]
0.25,  //[5]
0.875, //[6]
0.625, //[7]
1.75,  //[8]
3.625, //[9]
3.625, //[10]
14.5,  //[11]
12.0,  //[12]
};

static const double c_nonpolar[] =
{
0,
0.0, //[1]
0.0, //[2]
0.0, //[3]
0.0, //[4]
0.0, //[5]
0.0, //[6]
1.0, //[7]
1.0, //[8]
2.0, //[9]
2.0, //[10]
3.0, //[11]
3.0, //[12]
};

static const double d_polar[] =
{
0,
1.0, //[1]
1.0, //[2]
1.0, //[3]
3.0, //[4]
7.0, //[5]
1.0, //[6]
2.0, //[7]
5.0, //[8]
1.0, //[9]
1.0, //[10]
4.0, //[11]
2.0, //[12]
};

static const double t_polar[] =
{
0,
0.25,  //[1]
1.25,  //[2]
1.5,   //[3]
0.25,  //[4]
0.875, //[5]
2.375, //[6]
2.0,   //[7]
2.125, //[8]
3.5,   //[9]
6.5,   //[10]
4.75,  //[11]
12.5,  //[12]
};

static const double c_polar[] =
{
0,
0.0, //[1]
0.0, //[2]
0.0, //[3]
0.0, //[4]
0.0, //[5]
1.0, //[6]
1.0, //[7]
1.0, //[8]
2.0, //[9]
2.0, //[10]
2.0, //[11]
3.0, //[12]
};

""")

##print header
    
fluidconstants()

nonpolarfluids =['REFPROP-CO','REFPROP-COS','REFPROP-decane','REFPROP-H2S','REFPROP-ipentane','REFPROP-neopentn','REFPROP-ihexane','REFPROP-krypton','REFPROP-nonane','REFPROP-toluene','REFPROP-xenon','REFPROP-R116']

polarfluids =['REFPROP-acetone','REFPROP-N2O','REFPROP-SO2','REFPROP-R141b','REFPROP-R142b','REFPROP-R218','REFPROP-R245fa','REFPROP-R41']

import os
CPPName = os.path.join('..','CoolProp','purefluids','IndustrialFluids.cpp');
fp = open(CPPName,'w')
fp.write(CPPheader)
fp.close()

HName = os.path.join('..','CoolProp','purefluids','IndustrialFluids.h');
fp = open(HName,'w')
fp.write('#ifndef INDUSTRIALFLUIDS_H\n#define INDUSTRIALFLUIDS_H\n\n')
fp.close()


for Fluid in nonpolarfluids:
    buildEOS(Fluid,CPPName,polar=False)
    
for Fluid in nonpolarfluids:
    AddHeader(Fluid,HName)

for Fluid in polarfluids:
    buildEOS(Fluid,CPPName,polar=True)
    
for Fluid in polarfluids:
    AddHeader(Fluid,HName)
    
fp = open(HName,'a')
fp.write('#endif\n')
fp.close()

for Ref in nonpolarfluids+polarfluids:
    
    Fluid = Ref.split('-')[1]; ## The CoolProp name
    FluidString = RefLookup[Fluid][0]
    print 'FluidsList.push_back(new {Ref}Class());'.format(Ref=FluidString)
#plt.show()