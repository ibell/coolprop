import CoolProp
import CoolProp.CoolProp as CP

print 'Fluid | EOS | CP0 | VISCOSITY | CONDUCTIVITY | ECS_LENNARD_JONES | ECS_FITS | SURFACE_TENSION'

for fluid in CoolProp.__fluids__:
    
    print '{f:20s}'.format(f=fluid),'|',
    for key in ['EOS','CP0','VISCOSITY','CONDUCTIVITY','ECS_LENNARD_JONES','ECS_FITS','SURFACE_TENSION']:
        
        k = CP.get_BibTeXKey(fluid,key)
        if k and not k.startswith('__'):
            print 0,'|',
        elif k and k.startswith('__'):
            print 'X','|',
        else:
            print ' ','|',
    print ''
    
print '----------------------------------------------------'
for fluid in CoolProp.__fluids__:
    
    print '{f:20s}'.format(f=fluid),' & ',
    for i, key in enumerate(['EOS','CP0','VISCOSITY','CONDUCTIVITY','ECS_LENNARD_JONES','ECS_FITS','SURFACE_TENSION']):
        
        k = CP.get_BibTeXKey(fluid,key)
        if k and not k.startswith('__'):
            print '\cite{{{k:s}}}'.format(k=k),' & ',
        else:
            print ' ',' & ',
    print '\\\\'