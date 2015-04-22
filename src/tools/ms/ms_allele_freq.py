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

# quick and dirty MS formatted file parser
# calculates Allele Number (Count) distribution ['allele_num']
# calculates Sequence Weight (alleles/sample) distribution ['sequence_dist']

import sys,re,json

import numpy as np
import itertools

from optparse import OptionParser
from pprint import pprint

def pd( s0, s1 ):
    s0_it = len(s0)
    s1_it = len(s1)

    i = 0
    j = 0
    n=0
    while True:
        if i == s0_it:
            n += (s1_it - j)
            break

        if j == s1_it:
            n += (s0_it - i)
            break

        if (s0[i] < s1[j]):
            n+=1
            i += 1
        elif( s1[j] < s0[i]):
            n+=1
            j += 1
        else:
            i += 1
            j += 1

    return n

def computePairwiseDiff( seq_map, count ):
    if count == 0:
        return {}

    klist = []
    for x,y in seq_map.iteritems():
        if x == 'unique':
            continue
        klist.append(x)
        print x
        if 'dup' in y:
            for i in range(y['dup']):
                klist.append(x)
    np.random.shuffle(klist)

    keys = []
    if count > 0:
        for i in np.random.randint( len(klist), size=count):
            keys.append(klist[i])
    else:
        keys=klist

    diffs = {}

    s=len(keys)
    n=((s - 1) * s) / 2
    tot=0
    for i,x in enumerate(keys):
        if not x in diffs:
            diffs[x] = {}

        for y in itertools.islice(keys, i, s):
            if x == y:
                continue
            if y in diffs[x]:
                tot += diffs[x][y]
            elif y in diffs:
                if x in diffs[y]:
                    tot += diffs[y][x]
                else:
                    d = pd(seq_map[x]['pos'], seq_map[y]['pos'])
                    tot += d
                    diffs[y][x] = d
            else:
                d = pd( seq_map[x]['pos'], seq_map[y]['pos'] )
                tot += d
                diffs[x][y]=d

    return { "mean": (tot / n) }

def updateLogs( distributions, params ):
    (nReps, allele_dist, seq_dist, nSeqs, seq_map, seq_uniq, pairwise) = params

    distributions[nReps] = {}
    distributions[nReps]['allele_num'] = { x:y for x,y in enumerate(allele_dist) if y != 0 }
    distributions[nReps]['seq_dist'] = seq_dist
    distributions[nReps]['sequence_count'] = nSeqs
    distributions[nReps]['seq_map'] = seq_map
    distributions[nReps]['seq_map']['unique']=seq_uniq

    tmp = [0] * nSeqs
    for y in allele_dist:
        tmp[y] += 1

    distributions[nReps]['allele_dist'] = { x:y for x,y in enumerate(tmp) if y != 0}

    if not 'pairwise' in distributions[nReps]:
        distributions[nReps]['pairwise'] = {}

    m = 0
    for i in range(20):
        d = computePairwiseDiff( seq_map, pairwise )
        distributions[nReps]['pairwise'][i] = d
        m += d['mean']
    distributions[nReps]['pairwise']['average'] = (m / 20)

    return True

#################################################
## MAIN
#################################################

parser = OptionParser()
parser.add_option("-i", "--input", dest="file_path", help="Input MS formatted file", type="string", default="" )
parser.add_option("-g", "--gnu", dest="gnu_path", help="Path for CSV format (as GNU Plot tables)", type="string", default="")
parser.add_option("-j", "--json", dest="json_path", help="Path JSON format", type="string", default="")
parser.add_option("-p", "--pairwise", dest="pairwise", help="Random number of samples to compute pairwise difference", type="int", default=-1)

(options, args) = parser.parse_args()

if len(options.file_path) == 0:
    sys.exit(1)

print options

f=open( options.file_path, 'r')

ss_regex=re.compile('segsites: (\d+)')
one_regex=re.compile('1')

nReps = 0
nSeqs = 0
distributions = {}

allele_dist = []
seq_dist = {}

seq_map = {}

# parse header to first empty line
for line in f:
    line=line.rstrip('\n')
    if len(line) == 0:
        break

in_sample=False
nSites = 0
seq_uniq = 0
for line in f:
    line=line.rstrip('\n')
    if line.startswith('//'):
        nReps+=1
        in_sample=True
    elif len(line) == 0:
        updateLogs( distributions, (nReps, allele_dist, seq_dist, nSeqs, seq_map, seq_uniq, options.pairwise))

        # reset objects for next log
        nSeqs=0
        seq_dist = {}
        allele_dist = []
        seq_map = {}
        seq_uniq = 0
        in_sample=False
    elif ss_regex.match( line ):
        nSites=int(ss_regex.findall(line)[0])
        allele_dist = [0] * nSites
    elif len(line) == nSites:
        nSeqs+=1
        idx = line.find("1")
        n=0
        pos=[]
        while idx != -1:
            pos.append(idx)
            n += 1
            allele_dist[ idx ] += 1
            idx =line.find("1",idx + 1)

        dupped=False
        for x,y in seq_map.iteritems():
            if( cmp(y['pos'], pos) == 0 ):
                if 'dup' in y:
                    y['dup'] += 1
                else:
                    y['dup'] = 1
                dupped=True
                break

        if not dupped:
            seq_map[nSeqs] = { 'pos': pos }
            seq_uniq += 1

        if n in seq_dist:
            seq_dist[n] += 1
        else:
            seq_dist[n] = 1
f.close()

if in_sample:
    updateLogs( distributions, (nReps, allele_dist, seq_dist, nSeqs, seq_map, seq_uniq, options.pairwise))

print nReps

if len( options.json_path) != 0:
    with open(options.json_path + ".json", "w") as fp:
        json.dump( distributions, fp, indent=2, sort_keys=True )

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

        if 'seq_map' in y:
            path = options.gnu_path + ".seq_map." + str(x) + ".csv"
            w=open(path, 'w')
            w.write( "#Sequence Index,Allele List\n")
            for a, b in sorted( y['seq_map'].iteritems() ):
                w.write("%d" % a )
                for x in b['pos']:
                    w.write( ",%d", x)
                w.write("\n")
            w.close()
