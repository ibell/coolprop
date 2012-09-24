import pylab
import numpy as np

class FluidClass(object):
    pass

def is_float(a):
    try:
        a = float(a)
        return True
    except ValueError:
        return False
     
def parse_Melinder_nonaqueous():
    
    lines = open('Melinder coefficients.csv').readlines()
    
    #find #Fluid# dividers
    dividers = [i for i in range(len(lines)) if lines[i].find('#Fluid#') > -1]
    dividers.append(len(lines))
    
    for i in range(len(dividers)-1):
        Fluid = FluidClass()
        L = dividers[i]
        R = dividers[i+1]
        
        l = lines[L:R]
        l.pop(0) # #Fluid# remove
        Fluid.name, Fluid.longname = l.pop(0).split(',')[0].split('::')
        Fluid.description = l.pop(0).split(',')[0]
        l.pop(0) #empty line
        l.pop(0) #header
        l.pop(0) #units
        
        Fluid.t,Fluid.rho,Fluid.cp,Fluid.k,Fluid.mu = [],[],[],[],[]
        for line in l:
            ls = line.split(',')
            if is_float(ls[0]):
                tf,t,rho,cp,k,mu =  [float(thing) for thing in line.split(',')]
            else:
                if not is_float(ls[1]):
                    break
                t,rho,cp,k,mu =  [float(thing) for thing in line.split(',')[1::]]
            Fluid.t.append(t)
            Fluid.rho.append(rho)
            Fluid.cp.append(cp)
            Fluid.k.append(k)
            Fluid.mu.append(mu)
        pylab.semilogy(Fluid.t, Fluid.mu)
    
    pylab.show()
                

parse_Melinder_nonaqueous()

#for i in range(len(lines)):
#    if lines[i].strip() == '#Fluid#':
#        print lines[i],lines[i+1],lines[i+2]