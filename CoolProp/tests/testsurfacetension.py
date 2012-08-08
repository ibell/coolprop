
import CoolProp

errdict={}
wL=[0.1,0.3,0.5,0.7,0.9,0.95]
for Fluid in CoolProp.__fluids__:
    if Fluid=='Neopentane':
        #REFPROP doesn't have values for neopentane for some reason
        continue
    REFPROPname = CoolProp.CoolProp.get_REFPROPname(Fluid)
    Tc = CoolProp.CoolProp.Props(Fluid,'Tcrit')
    Tt = CoolProp.CoolProp.Props(Fluid,'Ttriple')
    print Fluid
    try:
        err = []
        tL = []
        for w in wL:
            T=w*Tc+(1-w)*Tt
            t=1-T/Tc
            IR=CoolProp.CoolProp.Props('I','T',T,'Q',0,'REFPROP-'+REFPROPname)
            IC=CoolProp.CoolProp.Props('I','T',T,'Q',0,Fluid)
            tL.append(t)
            err.append( (IR/IC-1)*100.0 )
            
        errdict[Fluid]=(err,tL)
        if Fluid=='Isopentane':
            print IR
            print IC
            print err
        
    except ValueError:
        print 'Error with',Fluid
    
import matplotlib.pyplot as plt
fig = plt.figure()
print '{:20s}{:s}'.format('Fluid','Max abs error (up to 90% of critical)')
for k,(v,tL) in errdict.iteritems():
    print '{:20s}{:f} %'.format(k,max([abs(_) for _ in v]))
    if k=='R717':
        plt.plot(tL,v,'-o')
    
plt.gca().set_xlabel(r'$t=1-T/Tc$')
plt.gca().set_ylabel('Error [%]')
plt.show()