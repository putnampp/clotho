#!/usr/bin/env python
#

# quick and dirty MS formatted file parser
# calculates Allele Number (Count) distribution ['allele_num']
# calculates Sequence Weight (alleles/sample) distribution ['sequence_dist']

import sys,re

from optparse import OptionParser
from pprint import pprint

parser = OptionParser()
parser.add_option("-i", "--input", dest="file_path", help="Input MS formatted file" )

(options, args) = parser.parse_args()

print options

f=open( options.file_path, 'r')

ss_regex=re.compile('segsites: (\d+)')
one_regex=re.compile('1')

nReps = 0
nSeqs = 0
distributions = {}

allele_dist = []
seq_dist = {}

nSites = 0 
for line in f:
    line=line.rstrip('\n')
    if line.startswith('//'):
        nReps+=1
    elif len(line) == 0:
        distributions[nReps] = {}
        distributions[nReps]['allele_num'] = { x:y for x,y in enumerate(allele_dist) if y != 0 }
        distributions[nReps]['seq_dist'] = seq_dist
        distributions[nReps]['sequence_count'] = nSeqs
        nSeqs=0
        seq_dist = {}
        allele_dist = []
    elif ss_regex.match( line ):
        nSites=int(ss_regex.findall(line)[0])
        allele_dist = [0] * nSites
    elif len(line) == nSites:
        nSeqs+=1
        idx = line.find("1")
        n=0
        while idx != -1:
            n += 1
            allele_dist[ idx ] += 1
            idx =line.find("1",idx + 1)
        if n in seq_dist:
            seq_dist[n] += 1
        else:
            seq_dist[n] = 1
f.close()

print nReps

pprint(distributions)
