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
parser.add_option("-g", "--gnu", dest="gnu_path", help="Path for CSV format (as GNU Plot tables)")
parser.add_option("-j", "--json", dest="json_path", help="Path JSON format")

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

# parse header to first empty line
for line in f:
    line=line.rstrip('\n')
    if len(line) == 0:
        break

in_sample=False
nSites = 0 
for line in f:
    line=line.rstrip('\n')
    if line.startswith('//'):
        nReps+=1
        in_sample=True
    elif len(line) == 0:
        distributions[nReps] = {}
        distributions[nReps]['allele_num'] = { x:y for x,y in enumerate(allele_dist) if y != 0 }
        distributions[nReps]['seq_dist'] = seq_dist
        distributions[nReps]['sequence_count'] = nSeqs

        tmp = [0] * nSeqs
        for y in allele_dist:
            tmp[y] += 1

        distributions[nReps]['allele_dist'] = { x:y for x,y in enumerate(tmp) if y != 0}

        nSeqs=0
        seq_dist = {}
        allele_dist = []
        in_sample=False
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

if in_sample:
    distributions[nReps] = {}
    distributions[nReps]['allele_num'] = { x:y for x,y in enumerate(allele_dist) if y != 0 }
    distributions[nReps]['seq_dist'] = seq_dist
    distributions[nReps]['sequence_count'] = nSeqs

    tmp = [0] * nSeqs
    for y in allele_dist:
        tmp[y] += 1

    distributions[nReps]['allele_dist'] = { x:y for x,y in enumerate(tmp) if y != 0}

print nReps

if len( options.json_path) != 0:
    w=open( options.json_path + ".json", "w")
    pprint(distributions,w)
    w.close()

if len( options.gnu_path ) != 0:
    for x,y in distributions.iteritems():
        # x = repeat
        # y = {...}
        if 'allele_num' in y:
            path = options.gnu_path + ".allele_num." + str(x) + ".csv"
            w=open( path, 'w')
            w.write( "#Allele,Count\n")
            for a,b in sorted(y['allele_num'].iteritems()):
                w.write( "%d,%d\n" % (a, b))
            w.close()

        if 'seq_dist' in y:
            path = options.gnu_path + ".seq_dist." + str(x) + ".csv"
            w=open( path, 'w')
            w.write( "#Allele Count,Sequence Count\n")
            for a,b in sorted(y['seq_dist'].iteritems()):
                w.write( "%d,%d\n" % (a, b))
            w.close()

        if 'allele_dist' in y:
            path = options.gnu_path + ".allele_dist." + str(x) + ".csv"
            w=open( path, 'w')
            w.write( "#Sequence Count,Alleles\n")
            for a,b in sorted(y['allele_dist'].iteritems()):
                w.write( "%d,%d\n" % (a, b))
            w.close()
