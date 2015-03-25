#!/usr/bin/env python
#

import simuOpt, os, sys, time

# use 'binary module' for representing long sequences 
simuOpt.setOptions(optimized=True, alleleType='binary')

import simuPOP as sim
import random
from itertools import izip
from itertools import islice

from pprint import pprint
from bisect import bisect
from optparse import OptionParser

# equal fitness function
def eq_fit():
    return 1

def record_time(pop, param=False):
    t = time.clock()
    if param:
        pop.dvars().perform[-1] = t - pop.dvars().perform[-1]
    else:
        pop.dvars().perform.append(t)
    return True
    

# modified version of model provided in User Guide (pg. 86-87)
def infSitesMutate(pop, param):
    (startPos, endPos, rate) = param

    # determine the number of mutation events occurring in this generation
    nEvents = sim.getRNG().randPoisson( rate )

#    print "Generating %d mutation events" % nEvents
    max_size = pop.popSize()
    if len( pop.dvars().variable_sites ) == 0:
        pop.dvars().variable_sites = [0] * pop.totNumLoci()

    lociset=set( pop.dvars().variable_sites )
    while nEvents > 0:
        nEvents -= 1
        # select an individual at random
        next_avail_idx = 0
        try:
            next_avail_idx = pop.dvars().variable_sites.index(0)
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

        # variable_sites serves a dual purpose
        # 1) it is a list of the genetic positions which have mutant alleles within the population
        # 2) those indices which are '0' in the list amount to 'free' loci within the population which can
        #    be used for new alleles
        pop.dvars().variable_sites[ next_avail_idx ] = loc
        lociset.add(loc)

    return True

class updateAlleleFreq(sim.PyOperator):
    def __init__(self,  *args, **kwargs):
        sim.PyOperator.__init__(self, func=self.update, *args, **kwargs)

    def update(self, pop):
        flist = []
        m = 0
        for x, y in pop.dvars().alleleNum.iteritems():
            if (1 not in y) or (y[1] == 0): # lost
                pop.dvars().variable_sites[x] = 0
            elif (y[0] == 0):               # fixed
                pop.dvars().variable_sites[x] = 0
                pop.dvars().nFixed += 1
                flist.append(x)
                print "Fixed %d" % x
            else:   # variable
                m = x   # m = last variable locus index

        pop.dvars().max_locus = m + 1 # variable locus index range = [0, max_locus) => contains all variable loci within population

        # loci with fixed alleles have to be reset
        # before loci can be re-used for next allele
        if len(flist) > 0:
            for ind in pop.individuals():
                for p in range(pop.ploidy()):
                    for x in flist:
                        ind.setAllele(0, x, ploidy=p)

        return True

class logPopulation(sim.PyOperator):
    def __init__(self, step=1, output="", final=False, *args, **kwargs):
        self.step=step
        self.called=step
        self.block=0
        self.output=output
        self.final=final
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

        if self.final:
            tmp_dict['perform'] = pop.dvars().perform
        
        lfile=open( self.output + "." + str(int(self.block * self.step)) + ".json", 'w')
        pprint( tmp_dict, lfile )
        lfile.close()

        self.called = self.step
        return True

    def allele_summary( self, pop, subPop ):
        tmp = [0] * pop.popSize() * pop.ploidy()
        for x in range(pop.totNumLoci()):
            tmp[int(pop.dvars().alleleNum[x][1])] += 1

        dist = { x: y for x,y in enumerate(tmp) if y != 0 }
        return { 'allele_distribution': dist, 'alleles_per_locus': pop.dvars().alleleNum }

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

        self.copier.copyChromosomes(parent, srcIdx , off, ploidy)
        if nEvents != 0:
            seen = set()
            while len(seen) < nEvents:
                seen.add( sim.getRNG().randInt( self.chromLen ) )
            
            breaks = list(seen)
            parent_ploidy = [ ((srcIdx + i) % 2) for i in range( nEvents + 1) ]

            # This loop iterates over all possible loci which have a variable allele
            # The list of loci padded to be larger than the expected number of loci (segregation sites)
            # In effect, we iterate over a longer list than is actually necessary
            # Therefore, we slice off the portion of the list after the last locus to have a non-zero allele
            # Although, the benefit of this 'early breaking' may be limited, if the number of segregations
            # sites approaches the padded size
            # It may be more efficient, in general, to simply assume that the padding will always be 'small'
            #  (ie. we have chosen the padding size such that it only increases the amount of work by a small
            #   amount)
            # Therefore iterating over the padded portion, although unnecessary, is more efficient because
            # it removes any additional conditional checking performed to determine 'early break' case
            for lociIdx,(x,y,z) in enumerate(islice(izip(parent.genotype((srcIdx % 2)), parent.genotype((srcIdx + 1) % 2), pop.dvars().variable_sites), 0, pop.dvars().max_locus)):

                # skip locus if the locus is free (z == 0) or parent is homozygous at locus (x == y)
                if z == 0 or x == y:
                    continue

                # if the allele falls within a recombination region not from the source chromosome
                if parent_ploidy[ bisect(breaks,z) ] != srcIdx:
                    off.setAllele( y, lociIdx, ploidy=ploidy)

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

parser = OptionParser()
parser.add_option( "-G", "--generations", dest="generations", help="Generations", default=100)
parser.add_option( "-P", "--pop_size", dest="pop_size", help="Population size", default=10000)
parser.add_option( "-m", "--mu", dest="mu", help="Mutation Rate per chromosome", default=0.001)
parser.add_option( "-r", "--rho", dest="rho", help="Recombination Rate per chromosome", default=0.001)
parser.add_option( "-l", "--log_period", dest="log_period", help="Logging step size", default=100)
parser.add_option( "-p", "--prefix", dest="log_prefix", help="Log path prefix", default='data/test')

(options, args) = parser.parse_args()

log_period=options.log_period   # log period
nGen=options.generations        # generations
popSize=options.pop_size        # population size
mu=options.mu                   # mutation per chromosome
mu_base=1e-8                    # mutation per base
rho=options.rho                 # recombination rate per chromosome

chromLen= int(mu / mu_base)  # (mutation/chromosome) / (mutation/base) = (base/chromosome)

theta=int(2 * popSize * mu)  # mutation events per generation in a diploid population

exp_seg_sites = 4 * popSize * mu * harmonic_number( popSize )   # Tajima's D

print "Expected number of segregation sites: %d" % exp_seg_sites
print "Chromosome Length: %d" % chromLen

# padding the number of expected segregation sites by 50%
pop=sim.Population( size=popSize, ploidy=2, loci=int(1.5 * exp_seg_sites), infoFields=['fitness'] )

# Generate Population

# variable_sites serves a dual purpose
# 1) it is a list of the genetic positions which have mutant alleles within the population
# 2) those indices which are '0' in the list amount to 'free' loci within the population which can
#    be used for new alleles
pop.evolve(
    initOps=[sim.InitGenotype( genotype=0 ), sim.PyExec('nFixed=0'), sim.PyExec('max_locus=0'), sim.PyExec('variable_sites=[]'), sim.PyExec( 'perform=[]')],
    preOps=[ sim.PyOperator( func=record_time ), sim.PySelector(func=eq_fit) ], 
    matingScheme=sim.HeteroMating([
                                    sim.HomoMating( sim.PyParentsChooser( myParentChooser ), sim.OffspringGenerator(ops=RegionRecombinator(rate=rho, param=(1, chromLen))), subPopSize=popSize)] ),
    postOps=[ sim.PyOperator( func=infSitesMutate, param=(1, chromLen, theta)),     # apply mutation after recombination
                sim.Stat( alleleFreq=sim.ALL_AVAIL, vars=['alleleNum'] ),
                updateAlleleFreq(),
                sim.PyOperator(func=record_time, param=True ),
                logPopulation( step=log_period, output="data/test" )
             ],
    finalOps=[ logPopulation( step=1, output="data/test.final", final=True) ],
    gen=nGen
)

# Print data
sim.Stat(pop, alleleFreq=sim.ALL_AVAIL, subPops=0 )

nSites = 0
for x in range( pop.totNumLoci() ):
    if pop.dvars().alleleNum[x][1] > 0:
        nSites += 1

print "Total Segregation Sites: %d" % nSites