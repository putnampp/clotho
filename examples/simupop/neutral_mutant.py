#   Copyright 2015 Patrick Putnam
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#!/usr/bin/env python
#

import simuOpt, os, sys, json

# use 'mutant module' for representing long sequences 
simuOpt.setOptions(optimized=True, alleleType='mutant')

import simuPOP as sim
import random
from itertools import izip
from itertools import islice

from pprint import pprint
from bisect import bisect
from optparse import OptionParser

import my_timer

# equal fitness function
def eq_fit():
    return 1

full_t= my_timer.timer()
tt = my_timer.timer()

def record_stop(pop):
    tt.stop()
    pop.dvars().perform.append(tt.elapsed_long())
    return True

def record_start(pop):
    global tt
    tt.start()
    return True

def record_full_start(pop):
    global full_t
    full_t.start()
    return True

def record_full_stop(pop):
    full_t.stop()
    pop.dvars().rt=full_t.elapsed_long()

def block_counter( pop ):
    tot = 0
    for ind in pop.individuals():
        for m in ind.mutants(0):
            tot += 1
        for m in ind.mutants(1):
            tot += 1

    pop.dvars().mem_usage.append(tot)
    return True
        

def myParentChooser(pop, subPop):
    psize=pop.subPopSize( subPop )
    while True:
        rank1=sim.getRNG().randInt(psize)
        rank2=sim.getRNG().randInt(psize)

        while rank2 == rank1:
            rank2 = sim.getRNG().randInt(psize)

        yield pop.individual(rank1, subPop), pop.individual( rank2, subPop)

def harmonic_number( s ):
    tot = .0
    for i in range(1, s - 1):
        tot += 1.0 / i
    return tot

def logConfiguration( config, output ):
    tmp_dict = {}
    tmp_dict['configuration'] = { x:y for x,y in config.__dict__.iteritems() }
    with open( output + ".config.json", 'wt') as lfile:
        res = json.dump( tmp_dict, lfile, sort_keys=True, indent=1)
    return True

class logPerformance(sim.PyOperator):
    def __init__(self, step=1, output="", config_opts=[], final=False, *args, **kwargs):
        self.step=step
        self.called=step
        self.block=0
        self.output=output
        self.config_opts=config_opts
        self.final=final
        sim.PyOperator.__init__(self, func=self.log, *args, **kwargs)

    def log(self, pop):
        self.called -= 1
        if self.called != 0:
            return True
        self.block += 1

        tmp_dict={}
        tmp_dict['perform'] = pop.dvars().perform
        tmp_dict['mem_usage'] = pop.dvars().mem_usage
        tmp_dict['runtime'] = pop.dvars().rt

        with open( self.output + "." + str(int(self.block * self.step)) + ".perform.json", 'wt') as lfile:
            res = json.dump( tmp_dict, lfile, sort_keys=True, indent=1)

        self.called = self.step
        return True


class logPopulation(sim.PyOperator):
    def __init__(self, step=1, output="", *args, **kwargs):
        self.step=step
        self.called=step
        self.block=0
        self.output=output
        sim.PyOperator.__init__(self, func=self.log, *args, **kwargs)

    def log(self, pop):
        self.called -= 1
        if self.called != 0:
            return True
        self.block += 1

        sim.stat( pop, alleleFreq=sim.ALL_AVAIL, vars=['alleleNum'] )

        tmp_dict={}
        tmp_dict['summary'] = self.population_summary(pop,0)
        tmp_dict['loci'] = self.allele_summary(pop, 0)

#        if self.final:
#            tmp_dict['perform'] = pop.dvars().perform
#            tmp_dict['configuration'] = { x:y for x,y in self.config_opts.__dict__.iteritems() }
#            tmp_dict['runtime'] = pop.dvars().rt
#            tmp_dict['mem_usage'] = pop.dvars().mem_usage
        
        with open( self.output + "." + str(int(self.block * self.step)) + ".pop.json", 'wt') as lfile:
            res = json.dump( tmp_dict, lfile, sort_keys=True, indent=1)

        self.called = self.step
        return True

    def allele_summary( self, pop, subPop ):
        tmp = [0] * (pop.popSize() * pop.ploidy() + 1)
        mut_loci={}

        tmp[0] = pop.popSize() * pop.ploidy()
        for x,y in pop.dvars().alleleNum.iteritems():
            if 1 in y:
                n = int(y[1])
                tmp[n] += 1
                tmp[0] -= n
                mut_loci[x] = y

        dist = { x: y for x,y in enumerate(tmp) if y != 0 }

        return { 'allele_distribution': dist, 'alleles_per_locus': mut_loci }

    def population_summary(self, pop, subPop, genotypes=False, site_distribution=True):
        if not (genotypes or site_distribution):
            return True

        genos = {}
        sites = {}

        for idx, x in enumerate(pop.individuals(subPop)):
            for y in range( pop.ploidy() ):
                count=sum(x.genotype(y))    # assumes alleles are 0 or 1; summation results in count of non-wildtype (mutated) loci
                if genotypes:
                    if idx in genos:
                        genos[idx][y] = count
                    else:
                        genos[idx] = {y: count}

                if site_distribution:
                    if count in sites:
                        sites[count] = sites[count] + 1
                    else:
                        sites[count] = 1

        tmp = {}
        if genotypes:
            tmp['genotypes'] = genos

        if site_distribution:
            tmp['sequence_distribution'] = sites
        return tmp


###################################################################################
##  Main program
###################################################################################

parser = OptionParser()
parser.add_option( "-G", "--generations", dest="generations", help="Generations", default=100, type="int")
parser.add_option( "-P", "--pop_size", dest="pop_size", help="Population size", default=10000, type="int")
parser.add_option( "-m", "--mu", dest="mu", help="Mutation Rate per chromosome", default=0.001, type="float")
parser.add_option( "-r", "--rho", dest="rho", help="Recombination Rate per chromosome", default=0.001, type="float")
parser.add_option( "-l", "--log_period", dest="log_period", help="Logging step size", default=100, type="int")
parser.add_option( "-p", "--prefix", dest="log_prefix", help="Log path prefix", default='data/test', type="string")

(options, args) = parser.parse_args()

log_period=options.log_period   # log period
nGen=options.generations        # generations
popSize=options.pop_size        # population size
mu=options.mu                   # mutation per chromosome
mu_base=1e-8                    # mutation per base
rho=options.rho                 # recombination rate per chromosome
rho_base=1e-8                   # recombination rate per base

chromLen= int(mu / mu_base)  # (mutation/chromosome) / (mutation/base) = (base/chromosome)

theta=int(2 * popSize * mu)  # mutation events per generation in a diploid population

exp_seg_sites = 4 * popSize * mu * harmonic_number( popSize )   # Tajima's D

# padding the number of expected segregation sites by 50%
pop=sim.Population( size=popSize, ploidy=2, loci=chromLen, infoFields=['fitness'] )

# Generate Population


logConfiguration( options, options.log_prefix )

pop.evolve(
    initOps=[   sim.PyExec( 'perform=[]'),
                sim.PyExec( 'mem_usage=[]' ),
                sim.PyExec( 'rt=0' ),
#                sim.InitGenotype( genotype=0 ), # increases start up cost
                sim.PyOperator(func=record_full_start)
            ],
    preOps=[    sim.PyOperator( func=record_start ),
                sim.PySelector(func=eq_fit) ], 
    matingScheme=sim.RandomSelection(ops=sim.Recombinator( rates=rho_base )),
    postOps=[ sim.SNPMutator( u=mu_base, v=0 ),     # apply mutation after recombination
                sim.PyOperator(func=record_stop )#,
#                logPopulation( step=log_period, output=options.log_prefix )
                , sim.PyOperator( func=block_counter )
                , logPerformance( step=log_period, output=options.log_prefix )
             ],
    finalOps=[  sim.PyOperator(func=record_full_stop )
                , logPopulation( step=1, output=options.log_prefix + ".final")
                , logPerformance( step=1, output=options.log_prefix + ".final")
            ],
    gen=nGen
)
