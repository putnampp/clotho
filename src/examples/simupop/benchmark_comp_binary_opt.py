#!/usr/bin/env python
#

import simuOpt, os, sys, time

# use 'binary module' for representing long sequences 
simuOpt.setOptions(optimized=True, alleleType='binary')

import simuPOP as sim
import random
from itertools import izip

from pprint import pprint
from bisect import bisect

# equal fitness function
def eq_fit():
    return 1

# modified version of model provided in User Guide (pg. 86-87)
def infSitesMutate(pop, param):
    (startPos, endPos, rate) = param

    # determine the number of mutation events occurring in this generation
    nEvents = sim.getRNG().randPoisson( rate )

#    print "Generating %d mutation events" % nEvents
    max_size = pop.popSize()
    if len( pop.dvars().allele_id ) == 0:
        pop.dvars().allele_id = [0] * pop.totNumLoci()

    lociset=set( pop.dvars().allele_id )
    while nEvents > 0:
        nEvents -= 1
        # select an individual at random
        next_avail_idx = 0
        try:
            next_avail_idx = pop.dvars().allele_id.index(0)
        except ValueError:
            print('Insufficient loci allocated %d; %d; %d' % (pop.totNumLoci(), pop.ancestralGens(), nEvents) )
            raise

        offset = sim.getRNG().randInt( max_size )   # randInt range: [0, max_size)

        indIdx = offset / pop.ploidy()
        pIdx= offset % pop.ploidy()
        
        loc = sim.getRNG().randInt( endPos ) + 1    # generate the location of a mutation event; range: [1, endPos]
        while loc in lociset:                       # generate a new location if mutation at location already mutated in population
            loc = sim.getRNG().randInt( endPos ) + 1

        pop.individual( indIdx ).setAllele( 1, next_avail_idx, ploidy=pIdx )

#        sites_ind.setAllele( loc, next_avail_idx, ploidy=0) # use metadata population individual as a 'free list' tracker
        pop.dvars().allele_id[ next_avail_idx ] = loc
        lociset.add(loc)
        if( next_avail_idx >= pop.dvars().max_locus ):
            pop.dvars().max_locus = next_avail_idx + 1

    return True

class updateAlleleFreq(sim.PyOperator):
    def __init__(self,  *args, **kwargs):
        sim.PyOperator.__init__(self, func=self.update, *args, **kwargs)

    def update(self, pop):
#        site_ind = pop.individual(0, pop.subPopByName("metadata"))

        for x, y in pop.dvars().alleleNum.iteritems():
            if x >= pop.dvars().max_locus:
                break

            if (1 not in y) or (y[1] == 0): # lost
                pop.dvars().allele_id[x] = 0
            elif (y[0] == 0): # fixed
                pop.dvars().allele_id[x] = 0
                pop.dvars().nFixed += 1

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

        tmp_dict={}
        tmp_dict['summary'] = self.population_summary(pop,0)
        tmp_dict['loci'] = self.allele_summary(pop, 0)
        tmp_dict['loci']['fixed'] = pop.dvars().nFixed
        
        lfile=open( self.output + "." + str(int(self.block * self.step)) + ".json", 'w')
        pprint( tmp_dict, lfile )
        lfile.close()

        self.called = self.step
        return True

    def allele_summary( self, pop, subPop ):
        dist = {}
        for x in range(pop.totNumLoci()):
            n = pop.dvars().alleleNum[x][1]
            if n in dist:
                dist[n] = dist[n] + 1
            else:
                dist[n] = 1
        return { 'allele_distribution': dist, 'alleles_per_locus': pop.dvars().alleleNum }

    def population_summary(self, pop, subPop, genotypes=False, site_distribution=True):
        if not (genotypes or site_distribution):
            return True

        genos = {}
        sites = {}

#        idx = 0
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
#            idx += 1

        tmp = {}
        if genotypes:
            tmp['genotypes'] = genos

        if site_distribution:
            tmp['sequence_distribution'] = sites
        return tmp

class RegionRecombinator( sim.PyOperator ):
    def __init__( self, rate, param, *args, **kwargs ):
        self.rate = rate
        (self.startPos, self.endPos) = param
        self.chromLen = self.endPos - self.startPos + 1

        self.copier = sim.GenoTransmitter()
        self.recombiner=self.recParent

        sim.PyOperator.__init__(self, func=self.transmitGenotype, *args, **kwargs)

    def transmitGenotype( self, pop, off, dad, mom ):
        self.recombiner(pop, off, dad, 0)
        self.recombiner(pop, off, mom, 1)

        return True

    def recParent( self, pop, off, parent, ploidy ):
        nEvents = sim.getRNG().randPoisson( self.rate )
        srcIdx = sim.getRNG().randInt( 2 )

        if nEvents == 0:
            self.copier.copyChromosomes(parent, srcIdx , off, ploidy)
        else:
            seen = set()
            while len(seen) < nEvents:
                seen.add( sim.getRNG().randInt( self.chromLen ) )
            
            breaks = list(seen)
            parent_ploidy = [ ((srcIdx + i) % 2) for i in range( nEvents + 1) ]

            lociIdx=0
            for (x,y,z) in izip(parent.genotype(0), parent.genotype(1), pop.dvars().allele_id):
                if x == y or parent_ploidy[ bisect(breaks,z) ] == 0:
                    off.setAllele( x, lociIdx, ploidy=ploidy)
                else:
                    off.setAllele( y, lociIdx, ploidy=ploidy)

                lociIdx += 1
                if lociIdx >= pop.dvars().max_locus:
                    break

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


###################################################################################
##  Main program
###################################################################################

log_period=100
nGen=100             # generations
popSize=10000           # population size
mu=0.001                # mutation per chromosome
mu_base=1e-8            # mutation per base
rho=30.0                # recombination rate per chromosome

chromLen= int(mu / mu_base)  # (mutation/chromosome) / (mutation/base) = (base/chromosome)

theta=int(2 * popSize * mu)  # mutation events per generation in a diploid population

exp_seg_sites = 4 * popSize * mu * harmonic_number( popSize )   # Tajima's D

print "Expected number of segregation sites: %d" % exp_seg_sites
print "Chromosome Length: %d" % chromLen

 # pad the number of expected segregation sites by 50%
pop=sim.Population( size=popSize, ploidy=2, loci=int(1.5 * exp_seg_sites), infoFields=['fitness'] )

# Generate Population
pop.evolve(
    initOps=[sim.InitGenotype( genotype=0 ), sim.PyExec('nFixed=0'), sim.PyExec('max_locus=0'), sim.PyExec('allele_id=[]')],
#    preOps=[ sim.PyOperator( func=infSitesMutate, param=(1, chromLen, theta)),  sim.PySelector(func=eq_fit) ], 
    preOps=[ sim.PySelector(func=eq_fit) ], 
    matingScheme=sim.HeteroMating([
                                    #sim.HomoMating( sim.PyParentsChooser( myParentChooser ), sim.OffspringGenerator(ops=sim.MendelianGenoTransmitter(subPops=0)), subPopSize=popSize, subPops=0),
                                    sim.HomoMating( sim.PyParentsChooser( myParentChooser ), sim.OffspringGenerator(ops=RegionRecombinator(rate=rho, param=(1, chromLen))), subPopSize=popSize)] ),
    postOps=[ sim.PyOperator( func=infSitesMutate, param=(1, chromLen, theta)),     # apply mutation after recombination
                sim.Stat( alleleFreq=sim.ALL_AVAIL, vars=['alleleNum'] ),
                updateAlleleFreq(),
                logPopulation( step=log_period, output="data/test" )
             ],
    gen=nGen
)

# Print data
sim.Stat(pop, alleleFreq=sim.ALL_AVAIL, subPops=0 )

nSites = 0
for x in range( pop.totNumLoci() ):
    if pop.dvars().alleleNum[x][1] > 0:
        nSites += 1

print "Total Segregation Sites: %d" % nSites
