from __future__ import division

import string
import random
import numpy as np
from scipy.odr import *

from summer import sum_function

LIBRARY = [i/6.0 for i in range(1,151)]+[0.35+i/2000 for i in range(1,100)]+[0.05+0.001*i for i in range(1,100)]

class Sample(object):
    def __init__(self,v):
        self.v = v
    
class GeneticHelloWorld(object):
    def __init__(self,
               num_samples = 60, # Have 60 chromos in the sample group
               num_selected = 20, # Have 10 chromos in the selected group
               mutation_factor = 3, # Randomly mutate 1/n of the chromosomes
               num_powers = 6, # How many powers in the fit
               Ref = 'REFPROP-R161',
               value = 'rhoL',
               addTr = False
                ):
        self.num_samples = num_samples
        self.num_selected = num_selected
        self.mutation_factor = mutation_factor
        self.num_powers = num_powers
        self.addTr = addTr
        self.value = value
        self.Ref = Ref 
        
        #THermodynamics
        from CoolProp.CoolProp import Props
        self.Tc = Props(Ref,'Tcrit')
        self.pc = Props(Ref,'pcrit')
        self.rhoc = Props(Ref,'rhocrit')
        self.Tmin = Props(Ref,'Tmin')
        
        self.T = np.linspace(self.Tmin+1e-6, self.Tc-0.000001,10000)
        self.p = [Props('P','T',T,'Q',0,Ref) for T in self.T]
        self.rhoL = [Props('D','T',T,'Q',0,Ref) for T in self.T]
        self.rhoV = [Props('D','T',T,'Q',1,Ref) for T in self.T]
        self.logppc = (np.log(self.p)-np.log(self.pc))
        self.rhoLrhoc = np.array(self.rhoL)/self.rhoc
        self.rhoVrhoc = np.array(self.rhoV)/self.rhoc
        self.logrhoLrhoc = np.log(self.rhoL)-np.log(self.rhoc)
        self.logrhoVrhoc = np.log(self.rhoV)-np.log(self.rhoc)

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
        x = 1.0-self.T/self.Tc
        
        if self.value == 'p':
            LHS = self.logppc
        elif self.value == 'rhoL':
            LHS = self.logrhoLrhoc
        elif self.value == 'rhoV':
            LHS = self.logrhoVrhoc
            
        if self.addTr:
            LHS *= self.T/self.Tc
        
        def f_coeffs(n,x):
        # The outer optimizer that will find the coefficients
        
            def f_RHS(B, x, n):
                return sum_function(B,x,np.array(chromo.v))
            
            linear = Model(f_RHS, extra_args = (n,))
            mydata = Data(x, LHS)
            myodr = ODR(mydata, linear, beta0=[0.0]*self.num_powers)
            myoutput = myodr.run()
            
            chromo.beta = myoutput.beta
        
            if self.addTr:
                RHS = f_RHS(myoutput.beta,x,n)*self.Tc/self.T
            else:
                RHS = f_RHS(myoutput.beta,x,n)
                
            if self.value == 'p':
                fit_value = np.exp(RHS)*self.pc
                EOS_value = self.p
            elif self.value in ['rhoL','rhoV']:
                fit_value = np.exp(RHS)*self.rhoc
                if self.value == 'rhoL':
                    EOS_value = self.rhoL
                else:
                    EOS_value = self.rhoV
                    
            max_abserror = np.max(np.abs((fit_value/EOS_value)-1)*100)
            
            #print self.Ref, 'rhoL SSE =', myoutput.sum_square, myoutput.beta, max_abserror,'%'
            print '.',
            return (myoutput.sum_square, max_abserror)
        
        chromo.fitness,chromo.max_abserror = f_coeffs(chromo,x)
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
        Breed two chromosomes by splicng them in a random spot and combining
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

        # Main loop: each generation select a subset of the sample and breed from
        # them.
        generation = -1
        while generation < 0 or samples[0].fitness > 1e-5 or generation < 3 and generation < 35:
            generation += 1
                
            # Generate the selected group from sample- take the top 5% of samples
            # and tourny select to generate the rest of selected.
            ten_percent = int(len(samples)*.05)
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
            
            print samples[0].v, list(samples[0].beta)
                
            # Print useful stats about this generation
            (min, median, max) =  [samples[0].fitness, samples[len(samples)//2].fitness, samples[-1].fitness]
            print("{0} best value: {1}. fitness: best {2}, median {3}, worst {4}".format(generation, samples[0].v, min, median, max))
    
        return samples[0]
  
def main():
    ghw = GeneticHelloWorld()
    r = ghw.run()
    
    print r.v
    print list(r.beta)

if __name__ == "__main__":
  main() 