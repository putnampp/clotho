#!/usr/bin/env python
#

import simuOpt, os, sys, time

# use 'mutant module' for representing long sequences 
simuOpt.setOptions(optimized=True, alleleType='long')

import simuPOP as sim
import random

# equal fitness function
def eq_fit():
    return 1

# modified version of model provided in User Guide (pg. 86-87)
def infSitesMutate(pop, param):
    (startPos, endPos, rate) = param

    # determine the number of mutation events occurring in this generation
    nEvents = sim.getRNG().randPoisson( rate )

#    print "Generating %d mutation events" % nEvents
    main_pop_idx = pop.subPopByName( "pop" )
    sites_ind = pop.individual( 0, pop.subPopByName( "metadata" ) )

    max_size = pop.subPopSize( main_pop_idx ) * pop.ploidy()
    while nEvents > 0:
        nEvents -= 1
        # select an individual at random
        next_avail_idx = 0
        try:
            next_avail_idx = sites_ind.genotype(0).index(0)
        except ValueError:
            print('Insufficient loci allocated %d; %d; %d' % (pop.totNumLoci(), pop.ancestralGens(), nEvents) )
            raise

        offset = sim.getRNG().randInt( max_size )   # randInt range: [0, max_size)

        indIdx = offset / pop.ploidy()
        pIdx= offset % pop.ploidy()
        
        loc = sim.getRNG().randInt( endPos ) + 1    # generate the location of a mutation event; range: [1, endPos]
        while loc in sites_ind.genotype(0):         # generate a new location if mutation at location already mutated in population
            loc = sim.getRNG().randInt( endPos ) + 1

        pop.individual( indIdx, main_pop_idx ).setAllele( 1, next_avail_idx, ploidy=pIdx )

        sites_ind.setAllele( loc, next_avail_idx, ploidy=0) # use metadata population individual as a 'free list' tracker
#        sites_ind.setAllele( 1, next_avail_idx, ploidy=1)   # allele frequency tracker

#    print sites_ind.genotype(0)
    return True

class updateAlleleFreq(sim.PyOperator):
    def __init__(self,  *args, **kwargs):
        sim.PyOperator.__init__(self, func=self.update, *args, **kwargs)

    def update(self, pop):
        site_ind = pop.individual(0, pop.subPopByName("metadata"))
        for x in range(pop.totNumLoci()):
#            print "(%d, %d) -> %d; (%d, %d) -> %d" % (x, 0, pop.vars(0)['alleleNum'][x][0], x, 1, pop.vars(0)['alleleNum'][x][1])
            if pop.vars(0)['alleleNum'][x][1] == 0 or pop.vars(0)['alleleNum'][x][0] == 0:  # [x][1] == 0 => lost; [x][0] == 0 => fixed
                site_ind.setAllele(0, x, ploidy=0)
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

nGen=100                # generations
popSize=10000           # population size
mu=0.001                # mutation rate per chromosome
mu_base=1e-8            # mutation rate per base
rho=0.001               # recombination rate

chromLen= int(mu / mu_base)  # (rate/chromosome) / (rate/base) = (base/chromosome)

theta=int(4 * popSize * mu)  # mutation events per generation in a diploid population

exp_seg_sites = theta * harmonic_number( popSize )   # Tajima's D

print "Expected number of segregation sites: %d" % exp_seg_sites
print "Chromosome Length: %d" % chromLen

pop=sim.Population( size=[popSize, 1], ploidy=2, loci=int(1.5 * exp_seg_sites), subPopNames=["pop", "metadata"], infoFields=['fitness'] ) # pad the number of expected segregation sites by 50%

pop.evolve(
    initOps=[sim.InitGenotype( genotype=0 )],
    preOps=[ sim.PyOperator( func=infSitesMutate, param=(1, chromLen, theta)),  sim.PySelector(func=eq_fit) ], 
    matingScheme=sim.HeteroMating([sim.HomoMating( sim.PyParentsChooser( myParentChooser ), sim.OffspringGenerator(ops=sim.MendelianGenoTransmitter(subPops=0)), subPopSize=popSize, subPops=0),
                                    sim.CloneMating( subPopSize=1, subPops=1)] ),
    postOps=[ sim.Stat( alleleFreq=sim.ALL_AVAIL, subPops=0, vars=['alleleNum_sp'] ),
#              sim.PyOperator( func=updateAlleleFreq, param=sim.ALL_AVAIL )
#             sim.PyEval(r"'%d: (2,0) => %d; (2,1) => %d\n' % (gen, subPop[0]['alleleNum'][2][0], subPop[0]['alleleNum'][2][1])")
                updateAlleleFreq()
             ],
    gen=nGen
)

sim.Stat(pop, alleleFreq=sim.ALL_AVAIL, subPops=0 )

nSites = 0
for x in range( pop.totNumLoci() ):
    print "(%d, %d) -> %d; (%d, %d) -> %d" % (x, 0, pop.vars(0)['alleleNum'][x][0], x, 1, pop.vars(0)['alleleNum'][x][1])
    if pop.vars(0)['alleleNum'][x][1] > 0:
        nSites += 1

print "Total Segregation Sites: %d" % nSites
