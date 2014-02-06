from __future__ import division

import string
import random
import numpy as np
from scipy.odr import *

import matplotlib.pyplot as plt

LIBRARY = [i/6.0 for i in range(1,151)]+[0.35+i/2000 for i in range(1,100)]+[0.05+0.001*i for i in range(1,100)]+[i+0.5 for i in range(10)]
    
class Sample(object):
    def __init__(self,v):
        self.v = v
    
class GeneticAncillaryFitter(object):
    def __init__(self,
               num_samples = 500, # Have this many chromos in the sample group
               num_selected = 30, # Have this many chromos in the selected group
               mutation_factor = 2, # Randomly mutate 1/n of the chromosomes
               num_powers = 8, # How many powers in the fit
               Ref = 'AceticAcid',
               value = 'rhoV',
               addTr = True,
               values = None
                ):
        self.num_samples = num_samples
        self.num_selected = num_selected
        self.mutation_factor = mutation_factor
        self.num_powers = num_powers
        self.addTr = addTr
        self.value = value
        self.Ref = Ref 
        
        #Thermodynamics
        from CoolProp.CoolProp import Props
        
        if values is None:
            self.Tc = Props(Ref,'Tcrit')
            self.pc = Props(Ref,'pcrit')
            self.rhoc = Props(Ref,'rhocrit')
            self.Tmin = Props(Ref,'Tmin')
            self.T = np.append(np.linspace(self.Tmin+1e-14, self.Tc-1,150), np.logspace(np.log10(self.Tc-1), np.log10(self.Tc-0.0001),10))
            self.p = [Props('P','T',T,'Q',0,Ref) for T in self.T]
            self.rhoL = [Props('D','T',T,'Q',0,Ref) for T in self.T]
            self.rhoV = [Props('D','T',T,'Q',1,Ref) for T in self.T]
        else:
            self.Tc = values['Tcrit']
            self.pc = values['pcrit']
            self.rhoc = values['rhocrit']
            self.Tmin = values['Tmin']
            self.T = values['T']
            self.p = values['p']
            self.rhoL = values['rhoL']
            self.rhoV = values['rhoV']
        
        self.logppc = (np.log(self.p)-np.log(self.pc))
        self.rhoLrhoc = np.array(self.rhoL)/self.rhoc
        self.rhoVrhoc = np.array(self.rhoV)/self.rhoc
        self.logrhoLrhoc = np.log(self.rhoL)-np.log(self.rhoc)
        self.logrhoVrhoc = np.log(self.rhoV)-np.log(self.rhoc)
        
        self.x = 1.0-self.T/self.Tc
        
        if self.value == 'p':
            self.LHS = self.logppc.copy()
        elif self.value == 'rhoL':
            self.LHS = self.logrhoLrhoc.copy()
        elif self.value == 'rhoV':
            self.LHS = self.logrhoVrhoc.copy()
        elif self.value == 'rhoLnoexp':
            self.LHS = (self.rhoLrhoc-1).copy()
        else:
            raise ValueError
            
        if self.addTr:
            self.LHS *= self.T/self.Tc

    def generate_random_chromosomes(self,):
        ''' 
        Create a list of random chromosomes to seed our alogrithm.
        '''
        chromos = []
        while len(chromos) < self.num_samples:
            chromos.append(Sample(sorted(random.sample(LIBRARY,self.num_powers))))
        return chromos

    def fitness(self, chromo):
        ''' 
        Fitness of a chromo is the sum of the squares of the error of the correlation
        '''
            
        output = np.zeros_like(self.T)
        
        def f_RHS(B,x):
            
            # see http://stackoverflow.com/questions/21309816/how-can-i-efficiently-use-numpy-to-carry-out-iterated-summation
            return np.sum(B * x.reshape(-1, 1)**chromo.v, axis = 1)
        
        linear = Model(f_RHS)
        mydata = Data(self.x, self.LHS)
        myodr = ODR(mydata, linear, beta0=[0.0]*self.num_powers)#, maxit = 100, sstol = 1e-10, partol = 1e-10)
        myoutput = myodr.run()
        
        chromo.beta = myoutput.beta
    
        if self.addTr:
            RHS = f_RHS(myoutput.beta, self.x)*self.Tc/self.T
        else:
            RHS = f_RHS(myoutput.beta, self.x)
            
        if self.value == 'p':
            fit_value = np.exp(RHS)*self.pc
            EOS_value = self.p
        elif self.value in ['rhoL','rhoV']:
            fit_value = np.exp(RHS)*self.rhoc
            if self.value == 'rhoL':
                EOS_value = self.rhoL
            elif self.value == 'rhoV':
                EOS_value = self.rhoV
        else:
            fit_value = self.rhoc*(1+RHS)
            EOS_value = self.rhoL
                
        max_abserror = np.max(np.abs((fit_value/EOS_value)-1)*100)
        
#        import matplotlib.pyplot as plt
#        plt.plot(x,fit_value,x,EOS_value)
#        plt.show()        
        
        chromo.fitness = max_abserror
        chromo.sum_square = myoutput.sum_square
        chromo.max_abserror = max_abserror
        
        print '.',
        
        return chromo.fitness

    def tourny_select_chromo(self, samples):
        ''' 
        Randomly select two chromosomes from the samples, then return the one
        with the best fitness.
        '''
        a = random.choice(samples)
        b = random.choice(samples)
        if a.fitness < b.fitness:
            return a
        else:
            return b

    def breed(self, a, b):
        ''' 
        Breed two chromosomes by splicing them in a random spot and combining
        them together to form two new chromos.
        '''
        splice_pos = random.randrange(len(a.v))
        new_a = a.v[:splice_pos] + b.v[splice_pos:]
        new_b = b.v[:splice_pos] + a.v[splice_pos:]
        return Sample(sorted(new_a)), Sample(sorted(new_b))

    def mutate(self, chromo):
        ''' 
        Mutate a chromosome by changing one of the parameters, but only if it improves the fitness
        '''
        v = chromo.v
        if hasattr(chromo,'fitness'):
            old_fitness = chromo.fitness
        else:
            old_fitness = self.fitness(chromo)
        
        for i in range(10):
            pos = random.randrange(len(chromo.v))
            chromo.v[pos] = random.choice(LIBRARY)
            new_fitness = self.fitness(chromo)
            if new_fitness < old_fitness:
                return chromo
            else:
                return Sample(sorted(v))

    def run(self):
        # Create a random sample of chromos
        samples = self.generate_random_chromosomes()
        
        # Calculate the fitness for the initial chromosomes
        for chromo in samples:
            self.fitness(chromo)
        print '#'
    
        decorated = [(sample.fitness,sample) for sample in samples]
        decorated.sort()
        samples = [s for sv,s in decorated]
        values = [sv for sv,s in decorated]
        plt.plot(values[0:len(values)//2])
        plt.close()

        # Main loop: each generation select a subset of the sample and breed from
        # them.
        generation = -1
        while generation < 0 or samples[0].fitness > 0.02 or generation < 3 and generation < 15:
            generation += 1
                
            # Generate the selected group from sample- take the top 10% of samples
            # and tourny select to generate the rest of selected.
            ten_percent = int(len(samples)*.1)
            selected = samples[:ten_percent]
            while len(selected) < self.num_selected:
                selected.append(self.tourny_select_chromo(samples))
            
            # Generate the solution group by breeding random chromos from selected
            solution = []
            while len(solution) < self.num_samples:
                solution.extend(self.breed(random.choice(selected),
                                         random.choice(selected)))
            
            # Apply a mutation to a subset of the solution set
            mutate_indices = random.sample(range(len(solution)), len(solution)//self.mutation_factor)
            for i in mutate_indices:
                solution[i] = self.mutate(solution[i])
            
            for chromo in solution:
                self.fitness(chromo)
            print '#'

            decorated = [(sample.fitness,sample) for sample in solution]
            decorated.sort()
            samples = [s for sv,s in decorated]
            
            print '------------------  Top 10 values  ---------------'
            for sample in samples[0:10]:
                print sample.v, sample.fitness, sample.max_abserror
            
            print '// Max error is ',samples[0].max_abserror,'% between',np.min(self.T),'and',np.max(self.T),'K'
            print str(samples[0].v)[1::], str(list(samples[0].beta))[1::]
                
            # Print useful stats about this generation
            (min, median, max) =  [samples[0].fitness, samples[len(samples)//2].fitness, samples[-1].fitness]
            print("{0} best value: {1}. fitness: best {2}, median {3}, worst {4}".format(generation, samples[0].v, min, median, max))
    
        return samples[0]
  
def main():
    gaf = GeneticAncillaryFitter()
    r = gaf.run()
    
    print r.v
    print list(r.beta)

if __name__ == "__main__":
    
    values = dict(Tcrit = 590.70, rhocrit = 351, pcrit = 5817526.2731115920, Tmin = 289.8)
    
    T,p,rhoL,rhoV = zip(*[[float(el) for el in line.strip().split(',')] for line in open('../aceticancillary.csv','r').readlines()])
    values['T'] = np.array(T)
    values['p'] = np.array(p)
    values['rhoL'] = np.array(rhoL)
    values['rhoV'] = np.array(rhoV)
    gaf = GeneticAncillaryFitter(value = 'rhoLnoexp', addTr = False, values = values)
    gaf.run()
    main() 