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
               num_samples = 100, # Have this many chromos in the sample group
               num_selected = 20, # Have this many chromos in the selected group
               mutation_factor = 2, # Randomly mutate 1/n of the chromosomes
               num_powers = 5, # How many powers in the fit
               Ref = 'R407C',
               value = 'rhoV',
               addTr = True,
               values = None,
               Tlims = None
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
            if Tlims is None:
                self.T = np.append(np.linspace(self.Tmin+1e-14, self.Tc-1,150), np.logspace(np.log10(self.Tc-1), np.log10(self.Tc)-1e-15,40))
            else:
                self.T = np.linspace(Tlims[0],Tlims[1])
            self.p = Props('P','T',self.T,'Q',0,Ref)
            self.rhoL = Props('D','T',self.T,'Q',0,Ref)
            self.rhoV = Props('D','T',self.T,'Q',1,Ref)
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
            self.EOS_value = self.p
        elif self.value == 'rhoL':
            self.LHS = self.logrhoLrhoc.copy()
            self.EOS_value = self.rhoL
        elif self.value == 'rhoV':
            self.LHS = self.logrhoVrhoc.copy()
            self.EOS_value = self.rhoV
        elif self.value == 'rhoLnoexp':
            self.LHS = (self.rhoLrhoc-1).copy()
            self.EOS_value = self.rhoL
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
        
        # theta^t where the i-th row is equal to theta^t[i]
        # We need these terms later on to build the A and b matrices
        theta_t = (self.x.reshape(-1, 1)**chromo.v).T
        
        # TODO: more numpy broadcasting should speed this up even more
        # Another few times faster ought to be possible
        I = len(chromo.v)
        A = np.zeros((I,I))
        b = np.zeros((I,1))
        for i in range(I):
            for j in range(I):
                A[i,j] = np.sum(theta_t[i]*theta_t[j])
            b[i] = np.sum(theta_t[i]*self.LHS)
        
        # If you end up with a singular matrix, quit this run
        try:
            n = np.linalg.solve(A, b).T
        except np.linalg.linalg.LinAlgError as E:
            chromo.fitness = 1e99
            return
        
        chromo.beta = n
        RHS = np.sum(n * self.x.reshape(-1, 1)**chromo.v, axis = 1)
            
        if self.addTr:
            RHS *= self.Tc/self.T
            
        if self.value == 'p':
            fit_value = np.exp(RHS)*self.pc
        elif self.value in ['rhoL','rhoV']:
            fit_value = np.exp(RHS)*self.rhoc
        else:
            fit_value = self.rhoc*(1+RHS)
        
        max_abserror = np.max(np.abs((fit_value/self.EOS_value)-1)*100)
        
        chromo.fitness = max_abserror
        chromo.fit_value = fit_value
        chromo.max_abserror = max_abserror
        
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
        while generation < 0 or samples[0].fitness > 0.02 or (generation < 3 and generation < 15):
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
            print str(samples[0].v), samples[0].beta.tolist()
                
            # Print useful stats about this generation
            (min, median, max) =  [samples[0].fitness, samples[len(samples)//2].fitness, samples[-1].fitness]
            print("{0} best value: {1}. fitness: best {2}, median {3}, worst {4}".format(generation, samples[0].v, min, median, max))
            
            # If the string representations of all the chromosomes are the same, stop
            if len(set([str(s.v) for s in samples[0:5]])) == 1:
                break
                
        self.fitness(samples[0])
    
        return samples[0]

# def test_water():
#     
#     gaf = GeneticAncillaryFitter(Ref = 'Water', value = 'p', addTr = False, num_powers = 6)
#     gaf.run()
#     
#     gaf = GeneticAncillaryFitter(Ref = 'Water', value = 'rhoL', addTr = False, num_powers = 6)
#     gaf.run()
#     
#     gaf = GeneticAncillaryFitter(Ref = 'Water', value = 'p', addTr = False, num_powers = 6)
#     gaf.run()

if __name__ == "__main__":
    
    gaf = GeneticAncillaryFitter(Ref = 'AceticAcid', value = 'rhoLnoexp', addTr = False, num_powers = 5, Tlims = (290,590))
    ch = gaf.run()
    plt.plot(gaf.T, ch.fit_value)
    plt.plot(gaf.T, gaf.EOS_value)
    plt.show()
    
#     import nose
#     nose.runmodule()